#!/usr/bin/env python3

import argparse
import json
import math
import shutil
import sys
from decimal import Decimal, ROUND_HALF_UP
from pathlib import Path

import pandas as pd


MISSING_VALUES = {"", "na", "nan", "none", "null"}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run GSEApy Enrichr over-representation analysis and summarize ranked terms."
    )
    parser.add_argument("--gene-list", required=True, help="Input gene list or table.")
    parser.add_argument("--gene-column", default="1", help="One-based column number or column name.")
    parser.add_argument(
        "--delimiter",
        choices=["line", "tab", "comma"],
        default="line",
        help="Input gene-list delimiter.",
    )
    parser.add_argument(
        "--has-header",
        choices=["true", "false"],
        default="false",
        help="Whether the gene-list table has a header row.",
    )
    parser.add_argument("--gene-sets", required=True, help="Enrichr library name or GMT file path.")
    parser.add_argument("--organism", default="human", help="Organism label passed to GSEApy.")
    parser.add_argument("--background", default="", help="Optional background gene list.")
    parser.add_argument(
        "--background-column",
        default="1",
        help="One-based column number or column name for the background file.",
    )
    parser.add_argument(
        "--background-delimiter",
        choices=["line", "tab", "comma"],
        default="line",
        help="Background gene-list delimiter.",
    )
    parser.add_argument(
        "--background-has-header",
        choices=["true", "false"],
        default="false",
        help="Whether the background table has a header row.",
    )
    parser.add_argument("--cutoff", type=float, default=1.0, help="GSEApy result cutoff.")
    parser.add_argument(
        "--sort-by",
        choices=["input_order", "p_value", "adjusted_p_value", "combined_score"],
        default="p_value",
        help="Ranking rule used for output and top-N summaries.",
    )
    parser.add_argument("--top-n", type=int, default=20, help="Number of ranked terms to report.")
    parser.add_argument(
        "--term-contains",
        default="",
        help="Optional case-insensitive substring to count among the top ranked terms.",
    )
    parser.add_argument("--output-results", required=True, help="Full standardized result table.")
    parser.add_argument("--output-top", required=True, help="Top-N standardized result table.")
    parser.add_argument("--output-summary", required=True, help="Term-substring summary table.")
    parser.add_argument("--metadata", required=True, help="Run metadata in JSON format.")
    return parser.parse_args()


def sep_from_name(delimiter):
    if delimiter == "tab":
        return "\t"
    if delimiter == "comma":
        return ","
    raise ValueError(f"Unsupported table delimiter for pandas parsing: {delimiter}")


def clean_value(value):
    if value is None:
        return None
    if isinstance(value, float) and math.isnan(value):
        return None
    text = str(value).strip()
    if text.lower() in MISSING_VALUES:
        return None
    return text


def select_column(frame, column):
    if column.isdigit():
        idx = int(column) - 1
        if idx < 0 or idx >= len(frame.columns):
            raise ValueError(f"Column {column} is outside the input table width")
        return frame.iloc[:, idx]
    if column not in frame.columns:
        raise ValueError(f"Column {column!r} was not found in the input table")
    return frame[column]


def read_gene_values(path, column, delimiter, has_header):
    path = Path(path)
    if delimiter == "line":
        values = path.read_text().splitlines()
        if has_header == "true" and values:
            values = values[1:]
    else:
        header = 0 if has_header == "true" else None
        frame = pd.read_csv(path, sep=sep_from_name(delimiter), header=header, dtype=str)
        values = select_column(frame, column).tolist()

    cleaned = []
    seen = set()
    for value in values:
        text = clean_value(value)
        if text is None:
            continue
        if text not in seen:
            cleaned.append(text)
            seen.add(text)
    if not cleaned:
        raise ValueError(f"No genes were read from {path}")
    return cleaned


def normalize_column_name(name):
    return (
        str(name)
        .strip()
        .lower()
        .replace("-", "_")
        .replace(" ", "_")
        .replace("/", "_")
        .replace(".", "_")
    )


def normalize_results(frame):
    column_map = {}
    for column in frame.columns:
        normalized = normalize_column_name(column)
        if normalized in {"gene_set", "geneset"}:
            column_map[column] = "gene_set"
        elif normalized == "term":
            column_map[column] = "term"
        elif normalized == "overlap":
            column_map[column] = "overlap"
        elif normalized in {"p_value", "pvalue"}:
            column_map[column] = "p_value"
        elif normalized in {"adjusted_p_value", "adjusted_pvalue", "fdr", "padj"}:
            column_map[column] = "adjusted_p_value"
        elif normalized == "odds_ratio":
            column_map[column] = "odds_ratio"
        elif normalized == "combined_score":
            column_map[column] = "combined_score"
        elif normalized == "genes":
            column_map[column] = "genes"
    normalized = frame.rename(columns=column_map).copy()

    required = {"term", "p_value"}
    missing = sorted(required - set(normalized.columns))
    if missing:
        raise ValueError("GSEApy result table is missing required columns: " + ", ".join(missing))

    for column in ["p_value", "adjusted_p_value", "odds_ratio", "combined_score"]:
        if column in normalized.columns:
            normalized[column] = pd.to_numeric(normalized[column], errors="coerce")

    ordered = [
        "gene_set",
        "term",
        "overlap",
        "p_value",
        "adjusted_p_value",
        "odds_ratio",
        "combined_score",
        "genes",
    ]
    for column in ordered:
        if column not in normalized.columns:
            normalized[column] = ""
    return normalized[ordered]


def sort_results(frame, sort_by):
    sorted_frame = frame.copy()
    if sort_by == "input_order":
        pass
    elif sort_by == "p_value":
        sorted_frame = sorted_frame.sort_values(
            ["p_value", "adjusted_p_value", "combined_score"],
            ascending=[True, True, False],
            na_position="last",
            kind="mergesort",
        )
    elif sort_by == "adjusted_p_value":
        sorted_frame = sorted_frame.sort_values(
            ["adjusted_p_value", "p_value", "combined_score"],
            ascending=[True, True, False],
            na_position="last",
            kind="mergesort",
        )
    elif sort_by == "combined_score":
        sorted_frame = sorted_frame.sort_values(
            ["combined_score", "p_value", "adjusted_p_value"],
            ascending=[False, True, True],
            na_position="last",
            kind="mergesort",
        )
    sorted_frame = sorted_frame.reset_index(drop=True)
    sorted_frame.insert(0, "rank", range(1, len(sorted_frame) + 1))
    return sorted_frame


def decimal_ratio(numerator, denominator):
    if denominator == 0:
        return Decimal("NaN")
    return Decimal(numerator) / Decimal(denominator)


def write_term_summary(top, args):
    keyword = args.term_contains.strip()
    denominator = len(top)
    if keyword:
        matches = top["term"].astype(str).str.contains(keyword, case=False, na=False)
        matching_terms = top.loc[matches, "term"].astype(str).tolist()
    else:
        matching_terms = []
    fraction = decimal_ratio(len(matching_terms), denominator)
    rounded = "" if fraction.is_nan() else str(fraction.quantize(Decimal("0.1"), rounding=ROUND_HALF_UP))
    summary = pd.DataFrame(
        [
            {
                "term_contains": keyword,
                "top_n_requested": args.top_n,
                "terms_considered": denominator,
                "matching_terms": len(matching_terms),
                "fraction": "" if fraction.is_nan() else str(fraction),
                "rounded_fraction_1_decimal": rounded,
                "matching_term_names": ";".join(matching_terms),
            }
        ]
    )
    summary.to_csv(args.output_summary, sep="\t", index=False)


def run_enrichr(args, genes, background):
    import gseapy as gp

    outdir = Path(args.metadata).parent / "gseapy_enrichr_report"
    outdir.mkdir(parents=True, exist_ok=True)
    gene_sets = args.gene_sets
    if Path(gene_sets).is_file():
        staged_gmt = outdir / "gene_sets.gmt"
        shutil.copyfile(gene_sets, staged_gmt)
        gene_sets = str(staged_gmt)
    result = gp.enrichr(
        gene_list=genes,
        gene_sets=gene_sets,
        organism=args.organism,
        background=background,
        outdir=str(outdir),
        cutoff=args.cutoff,
        no_plot=True,
        verbose=False,
    )
    if result.res2d is None or result.res2d.empty:
        raise ValueError("GSEApy returned no enrichment results")
    return result.res2d, gp.__version__


def main():
    args = parse_args()
    genes = read_gene_values(args.gene_list, args.gene_column, args.delimiter, args.has_header)
    background = None
    if args.background:
        background = read_gene_values(
            args.background,
            args.background_column,
            args.background_delimiter,
            args.background_has_header,
        )

    raw_results, gseapy_version = run_enrichr(args, genes, background)
    standardized = normalize_results(raw_results)
    ranked = sort_results(standardized, args.sort_by)

    top_n = max(args.top_n, 0)
    top = ranked.head(top_n).copy()

    ranked.to_csv(args.output_results, sep="\t", index=False)
    top.to_csv(args.output_top, sep="\t", index=False)
    write_term_summary(top, args)

    metadata = {
        "gseapy_version": gseapy_version,
        "gene_count": len(genes),
        "background_count": None if background is None else len(background),
        "gene_sets": args.gene_sets,
        "organism": args.organism,
        "cutoff": args.cutoff,
        "sort_by": args.sort_by,
        "top_n": top_n,
        "term_contains": args.term_contains,
        "result_terms": int(len(ranked)),
    }
    Path(args.metadata).write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)
