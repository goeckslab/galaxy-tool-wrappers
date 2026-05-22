#!/usr/bin/env python
"""Run selected nonparametric rank tests from tabular inputs."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

from scipy import stats


def column_index(value: str) -> int:
    idx = int(value)
    if idx < 1:
        raise argparse.ArgumentTypeError("Column numbers are one-based and must be positive")
    return idx - 1


def parse_bool(value: str) -> bool:
    return value.lower() in {"true", "t", "yes", "y", "1"}


def clean_float(value: str, *, row_number: int, column_number: int) -> float | None:
    value = value.strip()
    if value == "" or value.lower() in {"na", "nan", "none", "null"}:
        return None
    try:
        out = float(value)
    except ValueError as exc:
        raise ValueError(
            f"Non-numeric value at row {row_number}, column {column_number}: {value!r}"
        ) from exc
    if math.isnan(out):
        return None
    return out


def iter_rows(path: Path, has_header: bool):
    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        if has_header:
            next(reader, None)
        for row_number, row in enumerate(reader, start=2 if has_header else 1):
            if not row or all(cell.strip() == "" for cell in row):
                continue
            yield row_number, row


def values_from_column(path: Path, value_col: int, has_header: bool) -> list[float]:
    values: list[float] = []
    for row_number, row in iter_rows(path, has_header):
        if value_col >= len(row):
            raise ValueError(f"Row {row_number} has no column {value_col + 1}")
        value = clean_float(row[value_col], row_number=row_number, column_number=value_col + 1)
        if value is not None:
            values.append(value)
    return values


def paired_values_from_columns(
    path_a: Path,
    path_b: Path,
    value_col_a: int,
    value_col_b: int,
    has_header_a: bool,
    has_header_b: bool,
) -> tuple[list[float], list[float]]:
    rows_a = list(iter_rows(path_a, has_header_a))
    rows_b = list(iter_rows(path_b, has_header_b))
    if len(rows_a) != len(rows_b):
        raise ValueError("Wilcoxon signed-rank test requires the same number of paired rows")

    values_a: list[float] = []
    values_b: list[float] = []
    for (row_number_a, row_a), (row_number_b, row_b) in zip(rows_a, rows_b):
        if value_col_a >= len(row_a):
            raise ValueError(f"Group A row {row_number_a} has no column {value_col_a + 1}")
        if value_col_b >= len(row_b):
            raise ValueError(f"Group B row {row_number_b} has no column {value_col_b + 1}")

        value_a = clean_float(row_a[value_col_a], row_number=row_number_a, column_number=value_col_a + 1)
        value_b = clean_float(row_b[value_col_b], row_number=row_number_b, column_number=value_col_b + 1)
        if value_a is None or value_b is None:
            raise ValueError(
                "Wilcoxon signed-rank test requires complete paired observations; "
                f"missing value at Group A row {row_number_a} or Group B row {row_number_b}"
            )
        values_a.append(value_a)
        values_b.append(value_b)

    return values_a, values_b


def values_from_grouped_table(
    path: Path,
    value_col: int,
    group_col: int,
    has_header: bool,
    group_a: str,
    group_b: str,
) -> tuple[list[float], list[float], str, str]:
    groups: dict[str, list[float]] = {}
    order: list[str] = []
    for row_number, row in iter_rows(path, has_header):
        max_col = max(value_col, group_col)
        if max_col >= len(row):
            raise ValueError(f"Row {row_number} has no column {max_col + 1}")
        label = row[group_col].strip()
        if label not in groups:
            groups[label] = []
            order.append(label)
        value = clean_float(row[value_col], row_number=row_number, column_number=value_col + 1)
        if value is not None:
            groups[label].append(value)

    if (group_a == "") != (group_b == ""):
        raise ValueError("Provide both group labels, or leave both blank to auto-detect groups")

    if group_a == "" and group_b == "":
        if len(order) != 2:
            raise ValueError(
                "Exactly two groups are required when group labels are not supplied; "
                f"found {len(order)} groups"
            )
        group_a, group_b = order

    missing = [label for label in (group_a, group_b) if label not in groups]
    if missing:
        raise ValueError(f"Requested group label(s) not found: {', '.join(missing)}")
    return groups[group_a], groups[group_b], group_a, group_b


def format_float(value: float) -> str:
    return f"{value:.17g}"


def normalize_wilcoxon_method(method: str) -> str:
    if method == "asymptotic":
        return "approx"
    return method


def run_test(args):
    if args.input_mode == "grouped_table":
        if args.test == "wilcoxon_signed_rank":
            raise ValueError(
                "Wilcoxon signed-rank test requires paired observations; "
                "use two_tables input mode"
            )
        values_a, values_b, label_a, label_b = values_from_grouped_table(
            Path(args.input),
            args.value_column,
            args.group_column,
            args.has_header,
            args.group_a,
            args.group_b,
        )
    else:
        if args.test == "wilcoxon_signed_rank":
            values_a, values_b = paired_values_from_columns(
                Path(args.input_a),
                Path(args.input_b),
                args.value_column_a,
                args.value_column_b,
                args.has_header_a,
                args.has_header_b,
            )
        else:
            values_a = values_from_column(Path(args.input_a), args.value_column_a, args.has_header_a)
            values_b = values_from_column(Path(args.input_b), args.value_column_b, args.has_header_b)
        label_a = args.label_a
        label_b = args.label_b

    if not values_a or not values_b:
        raise ValueError("Both groups must contain at least one numeric value")

    if args.test == "mann_whitney_u":
        result = stats.mannwhitneyu(
            values_a,
            values_b,
            alternative=args.alternative,
            method=args.method,
            use_continuity=args.use_continuity,
        )
        statistic_name = "u_statistic"
        method = args.method
    else:
        if len(values_a) != len(values_b):
            raise ValueError("Wilcoxon signed-rank test requires paired groups with equal length")
        method = normalize_wilcoxon_method(args.wilcoxon_method)
        result = stats.wilcoxon(
            values_a,
            values_b,
            alternative=args.alternative,
            zero_method=args.zero_method,
            correction=args.correction,
            method=method,
        )
        statistic_name = "w_statistic"

    rows = [
        ("test", args.test),
        ("group_a", label_a),
        ("group_b", label_b),
        ("n_a", str(len(values_a))),
        ("n_b", str(len(values_b))),
        ("alternative", args.alternative),
        ("method", method),
        (statistic_name, format_float(float(result.statistic))),
        ("p_value", format_float(float(result.pvalue))),
    ]
    with Path(args.output).open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["field", "value"])
        writer.writerows(rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("--test", choices=["mann_whitney_u", "wilcoxon_signed_rank"], required=True)
    parser.add_argument("--input-mode", choices=["grouped_table", "two_tables"], required=True)
    parser.add_argument("--input")
    parser.add_argument("--value-column", type=column_index)
    parser.add_argument("--group-column", type=column_index)
    parser.add_argument("--has-header", type=parse_bool, default=False)
    parser.add_argument("--group-a", default="")
    parser.add_argument("--group-b", default="")
    parser.add_argument("--input-a")
    parser.add_argument("--input-b")
    parser.add_argument("--value-column-a", type=column_index)
    parser.add_argument("--value-column-b", type=column_index)
    parser.add_argument("--has-header-a", type=parse_bool, default=False)
    parser.add_argument("--has-header-b", type=parse_bool, default=False)
    parser.add_argument("--label-a", default="group_a")
    parser.add_argument("--label-b", default="group_b")
    parser.add_argument("--alternative", choices=["two-sided", "less", "greater"], default="two-sided")
    parser.add_argument("--method", choices=["auto", "asymptotic", "exact"], default="auto")
    parser.add_argument("--use-continuity", type=parse_bool, default=True)
    parser.add_argument("--wilcoxon-method", choices=["auto", "approx", "asymptotic", "exact"], default="auto")
    parser.add_argument("--zero-method", choices=["wilcox", "pratt", "zsplit"], default="wilcox")
    parser.add_argument("--correction", type=parse_bool, default=False)
    parser.add_argument("--output", required=True)
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    run_test(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
