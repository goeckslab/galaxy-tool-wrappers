#!/usr/bin/env python
"""Feature-wise correlation tests across matched tabular matrices."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import numpy as np
from scipy import stats


MISSING = {"", "na", "nan", "none", "null", "."}


def one_based(value: str) -> int:
    out = int(value)
    if out < 1:
        raise argparse.ArgumentTypeError("Column numbers are one-based and must be positive")
    return out - 1


def parse_bool(value: str) -> bool:
    return value.lower() in {"true", "t", "yes", "y", "1"}


def delimiter(value: str) -> str:
    if value == "tab":
        return "\t"
    if value == "comma":
        return ","
    raise argparse.ArgumentTypeError("delimiter must be tab or comma")


def clean_float(value: str) -> float:
    value = value.strip()
    if value.lower() in MISSING:
        return math.nan
    return float(value)


def read_first_row(path: Path, sep: str) -> list[str]:
    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter=sep)
        try:
            return next(reader)
        except StopIteration as exc:
            raise ValueError(f"{path} is empty") from exc


def feature_names(path: Path, sep: str, start_col: int, *, has_header: bool) -> list[str]:
    first_row = read_first_row(path, sep)
    if start_col >= len(first_row):
        raise ValueError("Feature start column is beyond the table width")
    if has_header:
        names = [value.strip() for value in first_row[start_col:]]
        empty_positions = [str(start_col + i + 1) for i, name in enumerate(names) if not name]
        if empty_positions:
            sample = ", ".join(empty_positions[:10])
            raise ValueError(f"Empty feature names are ambiguous: columns {sample}")
    else:
        names = [f"column_{i + 1}" for i in range(start_col, len(first_row))]
    seen: set[str] = set()
    duplicates: set[str] = set()
    for name in names:
        if name in seen:
            duplicates.add(name)
        seen.add(name)
    if duplicates:
        sample = ", ".join(sorted(duplicates)[:10])
        raise ValueError(f"Duplicate feature names are ambiguous: {sample}")
    return names


def read_matrix(
    path: Path,
    sep: str,
    id_col: int,
    selected_feature_indices: list[int],
    *,
    has_header: bool,
) -> tuple[list[str], np.ndarray]:
    rows: list[np.ndarray] = []
    ids: list[str] = []
    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter=sep)
        if has_header:
            next(reader, None)
        for row_number, row in enumerate(reader, start=2 if has_header else 1):
            if not row or all(cell.strip() == "" for cell in row):
                continue
            if id_col >= len(row):
                raise ValueError(f"Row {row_number} has no identifier column {id_col + 1}")
            obs_id = row[id_col].strip()
            if not obs_id:
                raise ValueError(f"Row {row_number} has an empty observation identifier")
            values = []
            for col in selected_feature_indices:
                if col >= len(row):
                    values.append(math.nan)
                else:
                    values.append(clean_float(row[col]))
            ids.append(obs_id)
            rows.append(np.asarray(values, dtype=float))
    if not rows:
        raise ValueError(f"{path} has no data rows")
    if len(set(ids)) != len(ids):
        raise ValueError(f"{path} contains duplicate observation identifiers")
    return ids, np.vstack(rows)


def align_rows(
    ids_a: list[str],
    matrix_a: np.ndarray,
    ids_b: list[str],
    matrix_b: np.ndarray,
) -> tuple[list[str], np.ndarray, np.ndarray]:
    b_pos = {obs_id: i for i, obs_id in enumerate(ids_b)}
    keep_a: list[int] = []
    keep_b: list[int] = []
    aligned_ids: list[str] = []
    for i, obs_id in enumerate(ids_a):
        j = b_pos.get(obs_id)
        if j is None:
            continue
        keep_a.append(i)
        keep_b.append(j)
        aligned_ids.append(obs_id)
    if not keep_a:
        raise ValueError("No observation identifiers overlap between the two matrices")
    return aligned_ids, matrix_a[keep_a, :], matrix_b[keep_b, :]


def benjamini_hochberg(pvalues: list[float]) -> list[float]:
    adjusted = [math.nan] * len(pvalues)
    finite = [(i, p) for i, p in enumerate(pvalues) if math.isfinite(p)]
    m = len(finite)
    if m == 0:
        return adjusted
    ranked = sorted(finite, key=lambda item: item[1])
    running = 1.0
    for rank_from_end, (idx, pvalue) in enumerate(reversed(ranked), start=1):
        rank = m - rank_from_end + 1
        running = min(running, pvalue * m / rank)
        adjusted[idx] = min(running, 1.0)
    return adjusted


def format_float(value: float) -> str:
    if not math.isfinite(value):
        return "nan"
    return f"{value:.17g}"


def run(args: argparse.Namespace) -> None:
    sep_a = delimiter(args.delimiter_a)
    sep_b = delimiter(args.delimiter_b)
    path_a = Path(args.matrix_a)
    path_b = Path(args.matrix_b)
    names_a = feature_names(path_a, sep_a, args.feature_start_a, has_header=args.has_header_a)
    names_b = feature_names(path_b, sep_b, args.feature_start_b, has_header=args.has_header_b)
    if args.has_header_a != args.has_header_b:
        raise ValueError("Both matrices must either have header rows or no header rows")
    if args.has_header_a:
        a_feature_pos = {name: i for i, name in enumerate(names_a)}
        b_feature_pos = {name: i for i, name in enumerate(names_b)}
        common_features = [name for name in names_a if name in b_feature_pos]
        if not common_features:
            raise ValueError("No feature names overlap between the two matrices")
        indices_a = [args.feature_start_a + a_feature_pos[name] for name in common_features]
        indices_b = [args.feature_start_b + b_feature_pos[name] for name in common_features]
    else:
        if len(names_a) != len(names_b):
            raise ValueError("No-header inputs must contain the same number of selected feature columns")
        common_features = names_a
        indices_a = list(range(args.feature_start_a, args.feature_start_a + len(names_a)))
        indices_b = list(range(args.feature_start_b, args.feature_start_b + len(names_b)))
    ids_a, matrix_a = read_matrix(path_a, sep_a, args.id_column_a, indices_a, has_header=args.has_header_a)
    ids_b, matrix_b = read_matrix(path_b, sep_b, args.id_column_b, indices_b, has_header=args.has_header_b)
    aligned_ids, matrix_a, matrix_b = align_rows(ids_a, matrix_a, ids_b, matrix_b)

    if args.transform_a == "negate":
        matrix_a = -matrix_a
    if args.transform_b == "negate":
        matrix_b = -matrix_b

    statistics: list[float] = []
    pvalues: list[float] = []
    pair_counts: list[int] = []
    for feature_index in range(len(common_features)):
        x = matrix_a[:, feature_index]
        y = matrix_b[:, feature_index]
        mask = np.isfinite(x) & np.isfinite(y)
        n = int(mask.sum())
        pair_counts.append(n)
        if n < args.min_pairs:
            statistics.append(math.nan)
            pvalues.append(math.nan)
            continue
        if args.method == "spearman":
            result = stats.spearmanr(x[mask], y[mask], alternative=args.alternative)
        else:
            result = stats.pearsonr(x[mask], y[mask], alternative=args.alternative)
        statistics.append(float(result.statistic))
        pvalues.append(float(result.pvalue))

    adjusted = benjamini_hochberg(pvalues)
    finite_tests = sum(math.isfinite(p) for p in pvalues)
    significant_count = sum(math.isfinite(q) and q <= args.alpha for q in adjusted)
    percent = 100.0 * significant_count / finite_tests if finite_tests else math.nan

    with Path(args.per_feature_output).open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["feature", "n_pairs", "statistic", "p_value", "p_adjust_bh", "significant"])
        for feature, n, stat, pvalue, qvalue in zip(common_features, pair_counts, statistics, pvalues, adjusted, strict=True):
            writer.writerow([
                feature,
                n,
                format_float(stat),
                format_float(pvalue),
                format_float(qvalue),
                "true" if math.isfinite(qvalue) and qvalue <= args.alpha else "false",
            ])

    summary_rows = [
        ("method", args.method),
        ("alternative", args.alternative),
        ("alpha", format_float(args.alpha)),
        ("matched_observations", str(len(aligned_ids))),
        ("matched_features", str(len(common_features))),
        ("finite_tests", str(finite_tests)),
        ("significant_count", str(significant_count)),
        ("significant_percent", format_float(percent)),
    ]
    with Path(args.summary_output).open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["field", "value"])
        writer.writerows(summary_rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix-a", required=True)
    parser.add_argument("--matrix-b", required=True)
    parser.add_argument("--delimiter-a", choices=["tab", "comma"], default="tab")
    parser.add_argument("--delimiter-b", choices=["tab", "comma"], default="tab")
    parser.add_argument("--has-header-a", type=parse_bool, default=True)
    parser.add_argument("--has-header-b", type=parse_bool, default=True)
    parser.add_argument("--id-column-a", type=one_based, default=0)
    parser.add_argument("--id-column-b", type=one_based, default=0)
    parser.add_argument("--feature-start-a", type=one_based, default=1)
    parser.add_argument("--feature-start-b", type=one_based, default=1)
    parser.add_argument("--transform-a", choices=["none", "negate"], default="none")
    parser.add_argument("--transform-b", choices=["none", "negate"], default="none")
    parser.add_argument("--method", choices=["spearman", "pearson"], default="spearman")
    parser.add_argument("--alternative", choices=["two-sided", "less", "greater"], default="two-sided")
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--min-pairs", type=int, default=3)
    parser.add_argument("--per-feature-output", required=True)
    parser.add_argument("--summary-output", required=True)
    return parser


def main() -> int:
    parser = build_parser()
    run(parser.parse_args())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
