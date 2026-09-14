#!/usr/bin/env python3
"""Validate the compiled paper, its six result tables, and AI disclosure PDF.

Requires PyMuPDF (``import fitz``). This is a deterministic text/geometry audit,
not a substitute for visual inspection or the authors' factual review.
"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import sys
import traceback


BASE = Path(__file__).resolve().parent
MARGIN_PT = 25 / 25.4 * 72
GLYPH_TOLERANCE_PT = 2.0
TABLES = {
    "q1_temperature": ("tab:q1-temperature", "1", ("30分钟", "温度")),
    "q1_moisture": ("tab:q1-moisture", "2", ("30分钟", "干基含水率")),
    "q2_temperature": ("tab:q2-temperature", "3", ("3小时", "温度")),
    "q2_moisture": ("tab:q2-moisture", "4", ("3小时", "干基含水率")),
    "q3_moisture": ("tab:q3-moisture", "5", ("固定半径", "干基含水率")),
    "q4_moisture": ("tab:q4-moisture", "6", ("收缩过程", "干基含水率")),
}
IDENTITY_PATTERNS = {
    "known_account": re.compile(r"(?i)(?<![a-z0-9])(?:Hel10o|libai)(?![a-z0-9])"),
    "windows_absolute_path": re.compile(r"(?<![A-Za-z0-9])[A-Za-z]:[\\/]"),
    "unix_user_path": re.compile(r"/(?:home|Users)/[^\s/]+/"),
    "windows_unc_path": re.compile(r"\\\\[A-Za-z0-9_.-]+\\[A-Za-z0-9_$.-]+"),
}
CELL = re.compile(r"^(?:[+-]?\d+(?:\.\d+)?|--|—|–|−)$")


def compact(text):
    return re.sub(r"\s+", "", text)


def input_record(path):
    result = {"path": str(path), "exists": path.is_file()}
    if path.is_file():
        result.update(bytes=path.stat().st_size,
                      sha256=hashlib.sha256(path.read_bytes()).hexdigest())
    return result


def group_words(words, y_tolerance=3.0):
    """Spatial rows; words within a row are always read from left to right."""
    rows = []
    for word in sorted(words, key=lambda w: ((w[1] + w[3]) / 2, w[0])):
        center = (word[1] + word[3]) / 2
        match = next((r for r in reversed(rows) if abs(r["y"] - center) <= y_tolerance), None)
        if match is None:
            rows.append({"y": center, "words": [word]})
        else:
            match["words"].append(word)
    for row in rows:
        row["words"].sort(key=lambda w: w[0])
        row["text"] = " ".join(w[4] for w in row["words"])
        row["bbox"] = [min(w[0] for w in row["words"]), min(w[1] for w in row["words"]),
                       max(w[2] for w in row["words"]), max(w[3] for w in row["words"])]
    return sorted(rows, key=lambda row: row["y"])


def page_spans(page):
    for block in page.get_text("dict")["blocks"]:
        for line in block.get("lines", []):
            yield from line["spans"]


def footer_word(word, page):
    # With generous bottom margins the footer can itself lie above the official
    # 25 mm line. Search the bottom 100 pt and select the lowest centered row.
    return (word[4].isdigit() and word[1] >= page.rect.height - 100
            and abs((word[0] + word[2]) / 2 - page.rect.width / 2) < 45)


def horizontal_rules(page):
    result = []
    for drawing in page.get_drawings():
        rect = drawing["rect"]
        if rect.width > 100 and rect.height < 1.5:
            rule = [float(rect.x0), float((rect.y0 + rect.y1) / 2), float(rect.x1)]
            if not any(max(abs(a - b) for a, b in zip(rule, old)) < 0.8 for old in result):
                result.append(rule)
    return sorted(result, key=lambda rule: rule[1])


def table_segment(page, caption_bottom=None, previous_bounds=None):
    """Locate a booktabs table by its vector rules, never by expected values."""
    rules = horizontal_rules(page)
    if caption_bottom is not None:
        candidates = [r for r in rules if caption_bottom < r[1] < caption_bottom + 65]
    else:
        candidates = [r for r in rules if MARGIN_PT - 3 < r[1] < MARGIN_PT + 150
                      and abs(r[0] - previous_bounds[0]) < 2
                      and abs(r[2] - previous_bounds[1]) < 2]
    if not candidates:
        raise ValueError("Caption/continuation has no matching upper horizontal table rule")
    top = candidates[0]
    full = [r for r in rules if r[1] >= top[1] - 0.2
            and abs(r[0] - top[0]) < 1.5 and abs(r[2] - top[2]) < 1.5]
    if len(full) < 3:
        raise ValueError("Cannot identify top, header separator, and bottom table rules")
    # The first equal-width line following the top separates the header from data.
    # Partial cmidrules are deliberately excluded by the x-coordinate condition.
    mid, bottom = full[1], full[2]
    words = [w for w in page.get_text("words") if top[0] - 0.5 <= w[0]
             and w[2] <= top[2] + 0.5 and mid[1] < (w[1] + w[3]) / 2 < bottom[1]]
    rows = group_words(words)
    extracted = []
    for row in rows:
        tokens = [compact(w[4]) for w in row["words"]]
        if not all(CELL.fullmatch(token) for token in tokens):
            raise ValueError(f"Non-cell text inside the bounded table body: {row['text']!r}")
        extracted.append(["--" if token in ("—", "–", "−") else token for token in tokens])
    return {"clip": [top[0], mid[1], top[2], bottom[1]], "rows": extracted,
            "row_bboxes": [row["bbox"] for row in rows],
            "bounds": [top[0], top[2]], "bottom_y": bottom[1]}


class Audit:
    def __init__(self, args):
        self.args = args
        self.report = {"schema_version": 1, "created_at_utc": datetime.now(timezone.utc).isoformat(),
                       "validation_scope": "body_preview" if args.body_only else "full_document",
                       "full_document_validated": False, "submission_ready": False,
                       "inputs": {name: input_record(getattr(args, name))
                                  for name in ("pdf", "aux", "ai_pdf", "expected")},
                       "checks": [], "table_results": {}, "documents": {},
                       "limitations": [
                           "Text and geometry audit does not replace rendered-page visual inspection.",
                           "Unknown author/school/team identifiers require human anonymity review.",
                           "AI disclosure truthfulness and complete human review cannot be certified by this script.",
                           "Table checks compare rendered strings with the supplied expected manifest; numerical model validity is outside scope."]}

    def check(self, key, passed, message, evidence=None):
        item = {"id": key, "status": "pass" if passed else "fail", "message": message}
        if evidence is not None:
            item["evidence"] = evidence
        self.report["checks"].append(item)

    def info(self, key, message, evidence=None):
        self.report["checks"].append({"id": key, "status": "info", "message": message,
                                      "evidence": evidence})

    def inspect_document(self, doc, name):
        meta = doc.metadata
        self.report["documents"][name] = {"page_count": len(doc), "metadata": meta}
        path = self.args.pdf if name == "paper" else self.args.ai_pdf
        self.check(f"{name}.size", path.stat().st_size < 20_000_000,
                   "PDF must be smaller than 20 MB (conservative decimal-byte threshold).",
                   {"bytes": path.stat().st_size, "limit_bytes_exclusive": 20_000_000})
        a4_bad, footers_bad, blank, bounds_bad, identity = [], [], [], [], []
        for pno, page in enumerate(doc, 1):
            if (abs(page.rect.width - 595.276) > 1 or abs(page.rect.height - 841.890) > 1
                    or page.rotation != 0):
                a4_bad.append({"page": pno, "size_pt": list(page.rect), "rotation": page.rotation})
            words = page.get_text("words")
            candidates = [w for w in words if footer_word(w, page)]
            lowest_y = max((w[1] for w in candidates), default=-1)
            feet = [w for w in candidates if abs(w[1] - lowest_y) < 1.5]
            if len(feet) != 1 or feet[0][4] != str(pno):
                footers_bad.append({"page": pno, "found": [w[4] for w in feet]})
            nonfooter = [w for w in words if w not in feet]
            if not nonfooter and not page.get_images() and not page.get_drawings():
                blank.append(pno)
            for span in page_spans(page):
                text, bbox = span["text"], span["bbox"]
                if not text.strip():
                    continue
                # Only a correctly positioned page-number span is exempt, not the bottom strip.
                if (text == str(pno) and any(abs(bbox[0] - w[0]) < 1.5
                                             and abs(bbox[1] - w[1]) < 1.5 for w in feet)):
                    continue
                excess = max(MARGIN_PT - bbox[0], MARGIN_PT - bbox[1],
                             bbox[2] - (page.rect.width - MARGIN_PT),
                             bbox[3] - (page.rect.height - MARGIN_PT), 0)
                if excess > GLYPH_TOLERANCE_PT:
                    bounds_bad.append({"page": pno, "text": text[:160], "bbox": list(bbox),
                                       "outside_25mm_box_pt": round(excess, 3)})
            text = page.get_text()
            for pattern_name, pattern in IDENTITY_PATTERNS.items():
                for match in pattern.finditer(text):
                    identity.append({"page": pno, "type": pattern_name,
                                     "context": text[max(0, match.start()-35):match.end()+55]})
        meta_text = json.dumps(meta, ensure_ascii=False)
        for pattern_name, pattern in IDENTITY_PATTERNS.items():
            if pattern.search(meta_text):
                identity.append({"page": "metadata", "type": pattern_name})
        if meta.get("author", "").strip():
            identity.append({"page": "metadata", "type": "nonempty_author", "value": meta["author"]})
        self.check(f"{name}.a4", not a4_bad, "All pages are unrotated portrait A4.", a4_bad)
        self.check(f"{name}.footer_sequence", not footers_bad,
                   "Every physical page has exactly one centered Arabic footer matching its 1-based index.", footers_bad)
        self.check(f"{name}.blank_pages", not blank, "No structurally empty pages.", blank)
        self.check(f"{name}.known_anonymity", not identity,
                   "No known account identifiers, absolute machine/user paths, or nonempty PDF author metadata.", identity)
        self.check(f"{name}.text_margins", not bounds_bad,
                   "Text lies inside the 25 mm margin box, allowing 2 pt for font bounding-box overhang; centered footer only is exempt.",
                   {"violation_count": len(bounds_bad), "violations": bounds_bad,
                    "margin_pt": MARGIN_PT, "glyph_bbox_tolerance_pt": GLYPH_TOLERANCE_PT})

    def parse_aux(self):
        text = self.args.aux.read_text(encoding="utf-8")
        records = re.findall(r"\\newlabel\{([^}]+)\}\{\{([^{}]*)\}\{([^{}]*)\}", text)
        labels = {key: {"number": number, "page": int(page) if page.isdigit() else page}
                  for key, number, page in records}
        duplicates = [key for key, count in Counter(key for key, _, _ in records).items() if count > 1]
        self.check("aux.unique_labels", not duplicates, "AUX labels have no duplicate definitions.", duplicates)
        cited = {key.strip() for group in re.findall(r"\\citation\{([^}]+)\}", text)
                 for key in group.split(",") if key.strip()}
        bibliography = re.findall(r"\\bibcite\{([^}]+)\}", text)
        undefined = sorted(cited - set(bibliography))
        repeated = [key for key, count in Counter(bibliography).items() if count > 1]
        self.check("aux.citation_keys", not undefined and not repeated,
                   "Every cited AUX key has exactly one bibliography definition.",
                   {"cited": sorted(cited), "undefined": undefined, "duplicate_bibcite": repeated})
        self.report["aux_labels"] = {key: value for key, value in labels.items()
                                     if key in ("body:start", "appendix:start", "document:end", "ai:declaration")
                                     or key.startswith("tab:q")}
        return labels

    def inspect_structure(self, doc, labels):
        body = labels.get("body:start", {}).get("page")
        app = labels.get("appendix:start", {}).get("page")
        end = labels.get("document:end", {}).get("page")
        ai = labels.get("ai:declaration", {}).get("page")
        self.check("paper.abstract_one_page", body == 2,
                   "AUX body:start must be physical page 2, leaving exactly one abstract page.", {"body_start": body})
        first = compact(doc[0].get_text()) if len(doc) else ""
        self.check("paper.abstract_content", "摘要" in first and "关键词" in first,
                   "The first page contains the abstract and keywords.")
        self.check("paper.no_submission_cover_pages", all(x not in first for x in ("承诺书", "编号专用页")),
                   "The electronic-paper PDF does not begin with a commitment or numbering cover.")
        if self.args.body_only:
            self.info("paper.appendix_not_validated", "Explicit body-only preview: appendix presence, completeness, and final ordering are not validated.")
            body_end = len(doc)
            self.check("paper.preview_boundary", app is None,
                       "A body-only preview must not silently contain an unvalidated appendix.", {"appendix_start": app})
        else:
            self.check("paper.appendix_boundary", isinstance(app, int) and isinstance(body, int)
                       and body < app <= len(doc),
                       "Full validation requires an appendix:start label after the body and within the PDF.",
                       {"body_start": body, "appendix_start": app, "pdf_pages": len(doc)})
            body_end = app - 1 if isinstance(app, int) else len(doc)
        self.check("paper.end_boundary", end == len(doc),
                   "document:end agrees with the actual final physical PDF page.", {"label_page": end, "pdf_pages": len(doc)})
        body_pages = body_end - body + 1 if isinstance(body, int) else None
        self.report["body_page_count"] = body_pages
        self.check("paper.body_page_limit", isinstance(body_pages, int) and 1 <= body_pages <= 30,
                   "Body, references, and AI declaration together occupy at most 30 pages, excluding abstract and appendix.",
                   {"body_start": body, "body_end": body_end, "body_pages": body_pages, "maximum": 30})
        self.check("paper.ai_declaration_location", isinstance(ai, int) and isinstance(body, int)
                   and body <= ai <= body_end,
                   "AI declaration is included in the counted body pages.", {"ai_declaration_page": ai})
        selected = list(doc)[:max(0, min(body_end, len(doc)))]
        unresolved, toc, refs = [], [], []
        fonts = Counter()
        for pno, page in enumerate(selected, 1):
            text = page.get_text()
            if "??" in text or "[?]" in text or "\ufffd" in text:
                unresolved.append(pno)
            for row in group_words(page.get_text("words")):
                line = compact(row["text"])
                if line.casefold() in ("目录", "contents", "tableofcontents"):
                    toc.append(pno)
                if line == "参考文献":
                    refs.append({"page": pno, "y": row["bbox"][1]})
            if isinstance(body, int) and pno >= body:
                for span in page_spans(page):
                    count = len(re.findall(r"[\u4e00-\u9fff]", span["text"]))
                    if count:
                        fonts[round(span["size"], 2)] += count
        self.check("paper.references_resolved", not unresolved,
                   "Abstract/body contain no ??, [?], or Unicode replacement characters; literal appendix source is excluded.", unresolved)
        self.check("paper.no_contents_page", not toc, "Abstract/body contain no table-of-contents heading.", toc)
        ai_text = compact(doc[ai-1].get_text()) if isinstance(ai, int) and 1 <= ai <= len(doc) else ""
        ref_after = len(refs) == 1 and isinstance(ai, int) and refs[0]["page"] >= ai
        if ref_after and refs[0]["page"] == ai:
            ref_after = ai_text.find("AI工具使用声明") >= 0 and ai_text.find("AI工具使用声明") < ai_text.find("参考文献")
        self.check("paper.ai_before_references", ref_after,
                   "Exactly one references heading follows the AI declaration, with appendix afterwards.", refs)
        self.info("paper.body_font_distribution", "CJK text-size distribution, weighted by character count; equations, captions, and tables may use other sizes.",
                  {"cjk_character_counts_by_pt": dict(fonts.most_common()),
                   "dominant_cjk_size_pt": fonts.most_common(1)[0][0] if fonts else None})

    def inspect_tables(self, doc, labels, expected):
        self.check("tables.manifest_keys", set(expected.get("tables", {})) == set(TABLES),
                   "The expected manifest defines exactly the six requested result tables.",
                   {"found": list(expected.get("tables", {})), "required": list(TABLES)})
        for key, (label, number, snippets) in TABLES.items():
            result = {"aux_label": label, "expected_table_number": number, "segments": []}
            self.report["table_results"][key] = result
            try:
                table = expected["tables"][key]
                anchor = labels[label]
                if anchor["number"] != number:
                    raise ValueError(f"AUX number {anchor['number']!r} does not match {number!r}")
                pno = anchor["page"]
                if not isinstance(pno, int) or not 1 <= pno <= len(doc):
                    raise ValueError("AUX table page is absent or outside this PDF")
                captions = [row for row in group_words(doc[pno-1].get_text("words"))
                            if re.match(rf"^表\s*{number}(?!\d)", row["text"])
                            and all(s in compact(row["text"]) for s in snippets)]
                if len(captions) != 1:
                    raise ValueError(f"Expected one spatial caption matching table {number}; found {len(captions)}")
                result["caption"] = {"page": pno, "text": captions[0]["text"], "bbox": captions[0]["bbox"]}
                rows = []
                bounds = None
                for segment_index in range(3):
                    page = doc[pno-1]
                    segment = table_segment(page, captions[0]["bbox"][3] if segment_index == 0 else None, bounds)
                    result["segments"].append({"page": pno, "clip": segment["clip"],
                                               "row_bboxes": segment["row_bboxes"], "rows": segment["rows"]})
                    rows.extend(segment["rows"])
                    bounds = segment["bounds"]
                    if len(rows) >= len(table["rows"]):
                        break
                    # Only page-bottom segments may continue. Missing rows in a small
                    # ordinary table are a failure, not permission to search later text.
                    if segment["bottom_y"] < page.rect.height - MARGIN_PT - 100 or pno >= len(doc):
                        break
                    pno += 1
                result["expected_row_count"] = len(table["rows"])
                result["actual_row_count"] = len(rows)
                result["column_count"] = table["column_count"]
                mismatches = []
                if len(rows) != len(table["rows"]):
                    mismatches.append({"type": "row_count", "expected": len(table["rows"]), "actual": len(rows)})
                for r, (wanted, actual) in enumerate(zip(table["rows"], rows), 1):
                    if len(actual) != table["column_count"]:
                        mismatches.append({"type": "column_count", "row": r,
                                           "expected": table["column_count"], "actual": len(actual)})
                    for c in range(max(len(wanted), len(actual))):
                        want = wanted[c] if c < len(wanted) else None
                        got = actual[c] if c < len(actual) else None
                        if want != got:
                            mismatches.append({"type": "cell", "row": r, "column": c+1, "expected": want, "actual": got})
                result["mismatches"] = mismatches
                result["matched_cells"] = sum(1 for wanted, actual in zip(table["rows"], rows)
                                               for want, got in zip(wanted, actual) if want == got)
                result["status"] = "fail" if mismatches else "pass"
                self.check(f"tables.{key}", not mismatches,
                           "Caption- and vector-rule-bounded PDF cells match all expected strings, including times and domain-outside dashes.",
                           {"pages": [s["page"] for s in result["segments"]],
                            "matched_cells": result["matched_cells"], "mismatch_count": len(mismatches)})
            except Exception as exc:
                result.update(status="fail", error=f"{type(exc).__name__}: {exc}")
                self.check(f"tables.{key}", False, "Table could not be located or checked; no global-token fallback was used.", result["error"])

    def run(self):
        try:
            import fitz
            self.report["runtime"] = {"python": sys.version, "executable": sys.executable,
                                      "pymupdf": fitz.VersionBind}
        except ImportError as exc:
            self.check("runtime.pymupdf", False, "PyMuPDF is required; install it or use the documented project validation interpreter.", str(exc))
            return
        for name in ("pdf", "aux", "ai_pdf", "expected"):
            self.check(f"input.{name}", getattr(self.args, name).is_file(), f"Required input {name} exists.")
        if not all(getattr(self.args, name).is_file() for name in ("pdf", "aux", "expected")):
            self.info("paper.not_run", "Paper validation was not run because a required PDF/AUX/expected input is missing.")
            return
        with fitz.open(self.args.pdf) as doc:
            self.inspect_document(doc, "paper")
            labels = self.parse_aux()
            self.inspect_structure(doc, labels)
            expected = json.loads(self.args.expected.read_text(encoding="utf-8-sig"))
            self.inspect_tables(doc, labels, expected)
        if self.args.ai_pdf.is_file():
            with fitz.open(self.args.ai_pdf) as ai_doc:
                self.inspect_document(ai_doc, "ai_details")
                text = compact("\n".join(page.get_text() for page in ai_doc))
                required = {"title": "AI工具使用详情" in text,
                            "tool_and_version": "工具" in text and ("版本" in text or "模型" in text),
                            "purpose_or_stage": any(term in text for term in ("目的", "环节", "分工")),
                            "input_and_interaction": "输入" in text and "交互" in text,
                            "adoption_and_review": "采纳" in text and any(term in text for term in ("核验", "验证", "审核"))}
                self.check("ai_details.required_topics", all(required.values()),
                           "AI-details PDF contains the expected disclosure topic labels; this does not certify the statements' truth or completeness.", required)
                self.info("ai_details.human_review_external", "Actual authors must verify AI outputs and complete truthful records; automated topic detection does not establish human review.")

    def save(self):
        changed = [name for name, original in self.report["inputs"].items()
                   if input_record(getattr(self.args, name)) != original]
        self.check("inputs.stable_during_validation", not changed,
                   "No input file changed while this audit was running.", changed)
        fails = [item["id"] for item in self.report["checks"] if item["status"] == "fail"]
        self.report["failed_check_ids"] = fails
        self.report["passed"] = not fails
        self.report["full_document_validated"] = not fails and not self.args.body_only
        self.report["status"] = "failed" if fails else ("preview_passed" if self.args.body_only else "automated_checks_passed")
        self.report["check_counts"] = dict(Counter(item["status"] for item in self.report["checks"]))
        self.args.output.parent.mkdir(parents=True, exist_ok=True)
        self.args.output.write_text(json.dumps(self.report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        lines = [f"PDF validation: {self.report['status']}",
                 f"Scope: {self.report['validation_scope']}",
                 f"Paper: {self.args.pdf}", f"AI details: {self.args.ai_pdf}",
                 f"Checks: {self.report['check_counts']}",
                 f"Counted body pages: {self.report.get('body_page_count', 'not checked')}",
                 f"Full document validated: {self.report['full_document_validated']}",
                 "Submission-ready certification: not provided (visual and human factual review remain external)."]
        for check in self.report["checks"]:
            if check["status"] == "fail":
                lines.append(f"FAIL {check['id']}: {check['message']}")
                if check.get("evidence") is not None:
                    lines.append(json.dumps(check["evidence"], ensure_ascii=False))
        for key, result in self.report["table_results"].items():
            lines.append(f"TABLE {key}: {result['status']}; matched cells={result.get('matched_cells', 0)}")
        lines.extend(["", "Limitations:"] + self.report["limitations"])
        summary = "\n".join(lines) + "\n"
        self.args.output.with_suffix(".txt").write_text(summary, encoding="utf-8")
        print(summary)
        print(f"JSON: {self.args.output}")
        return 1 if fails else 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdf", type=Path, default=BASE / "build/main.pdf")
    parser.add_argument("--aux", type=Path, default=BASE / "build/main.aux")
    parser.add_argument("--ai-pdf", type=Path, default=BASE / "build/ai_details.pdf")
    parser.add_argument("--expected", type=Path, default=BASE / "evidence/table_expected.json")
    parser.add_argument("--output", type=Path, default=BASE / "evidence/pdf_validation.json")
    parser.add_argument("--body-only", action="store_true", help="Explicit preview audit: appendix is unavailable and full validation cannot pass.")
    args = parser.parse_args()
    for key in ("pdf", "aux", "ai_pdf", "expected", "output"):
        setattr(args, key, getattr(args, key).resolve())
    audit = Audit(args)
    try:
        audit.run()
    except Exception as exc:
        audit.check("validator.exception", False, "Validation aborted; unchecked items must not be considered passed.",
                    {"error": f"{type(exc).__name__}: {exc}", "traceback": traceback.format_exc()})
    return audit.save()


if __name__ == "__main__":
    raise SystemExit(main())
