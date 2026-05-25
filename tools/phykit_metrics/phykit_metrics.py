#!/usr/bin/env python

import argparse
import csv
import os
import re
import subprocess
import sys
import tempfile
import zipfile
from dataclasses import dataclass
from pathlib import PurePosixPath


TREE_METRICS = {
    "evolutionary_rate": "evolutionary_rate",
    "total_tree_length": "total_tree_length",
    "long_branch_score": "long_branch_score",
    "treeness": "treeness",
}
ALIGNMENT_METRICS = {
    "relative_composition_variability": "relative_composition_variability",
    "relative_composition_variability_taxon": "relative_composition_variability_taxon",
}
PER_TAXON_METRICS = {"relative_composition_variability_taxon"}


@dataclass
class Member:
    path: str
    name: str
    group: str
    selected: bool = True
    archive_name: str = ""

    @property
    def file_id(self):
        base = os.path.basename(self.name)
        if "." in base and base.split(".", 1)[0]:
            return base.split(".", 1)[0]
        stem, _ = os.path.splitext(base)
        return stem or base


def clean_error(text):
    return " | ".join(line.strip() for line in text.splitlines() if line.strip())


def parse_numeric_token(tokens):
    for token in reversed(tokens):
        try:
            float(token)
        except ValueError:
            continue
        return token
    return None


def parse_scalar_output(stdout):
    lines = [line.strip() for line in stdout.splitlines() if line.strip()]
    if not lines:
        raise ValueError("PhyKIT produced no output")
    value = parse_numeric_token(lines[-1].split())
    if value is None:
        raise ValueError("PhyKIT output did not contain a numeric value")
    return value


def parse_taxon_output(stdout):
    rows = []
    for line in stdout.splitlines():
        line = line.strip()
        if not line:
            continue
        fields = line.split()
        value = parse_numeric_token(fields[1:])
        if len(fields) < 2 or value is None:
            raise ValueError(f"Could not parse per-taxon PhyKIT output line: {line}")
        rows.append((fields[0], value))
    if not rows:
        raise ValueError("PhyKIT produced no per-taxon rows")
    return rows


def run_phykit(metric, path, verbose):
    command = TREE_METRICS.get(metric) or ALIGNMENT_METRICS.get(metric)
    if command is None:
        raise ValueError(f"Unsupported metric: {metric}")
    cmd = ["phykit", command, path]
    if metric == "long_branch_score" and verbose:
        cmd.append("--verbose")
    return subprocess.run(cmd, capture_output=True, text=True)


def per_taxon_metric(metric, verbose):
    return metric in PER_TAXON_METRICS or (metric == "long_branch_score" and verbose)


def parse_members(args):
    paths = args.member or []
    names = args.member_name or []
    groups = args.member_group or []
    if len(paths) != len(names):
        raise ValueError("Internal error: collection member paths and names differ")
    if args.mode == "grouped_collections" and len(paths) != len(groups):
        raise ValueError("Internal error: grouped collection member paths and labels differ")

    members = []
    for i, path in enumerate(paths):
        name = names[i] or os.path.basename(path)
        if args.mode == "collection":
            group = args.collection_label
        else:
            group = groups[i]
        members.append(Member(path=path, name=name, group=group))
    return members


def safe_zip_member_name(name):
    path = PurePosixPath(name)
    if path.is_absolute() or any(part in {"", ".", ".."} for part in path.parts):
        raise ValueError(f"Unsafe ZIP member path: {name}")
    return path


def archive_members(args):
    archives = args.archive or []
    names = args.archive_name or []
    groups = args.archive_group or []
    if names and len(archives) != len(names):
        raise ValueError("Internal error: archive paths and names differ")
    if args.mode == "grouped_zip_archives" and len(archives) != len(groups):
        raise ValueError("Internal error: grouped archive paths and labels differ")

    members = []
    for i, archive in enumerate(archives):
        archive_name = names[i] if names else os.path.basename(archive)
        group = args.archive_label if args.mode == "zip_archive" else groups[i]
        extract_root = tempfile.mkdtemp(prefix="phykit_metrics_zip_")
        with zipfile.ZipFile(archive) as archive_handle:
            for info in archive_handle.infolist():
                if info.is_dir():
                    continue
                member_path = safe_zip_member_name(info.filename)
                output_path = os.path.join(extract_root, *member_path.parts)
                os.makedirs(os.path.dirname(output_path), exist_ok=True)
                with archive_handle.open(info) as source, open(output_path, "wb") as target:
                    target.write(source.read())
                members.append(
                    Member(
                        path=output_path,
                        name=info.filename,
                        group=group,
                        archive_name=archive_name,
                    )
                )
    return members


def apply_filters(members, suffix, regex):
    pattern = re.compile(regex) if regex else None
    for member in members:
        if suffix and not member.name.endswith(suffix):
            member.selected = False
        if pattern is not None and not pattern.search(member.name):
            member.selected = False
    return members


def write_single(args):
    result = run_phykit(args.metric, args.single_input, args.verbose)
    with open(args.output, "w") as output:
        output.write(result.stdout)
    if result.returncode != 0:
        if result.stderr:
            sys.stderr.write(result.stderr)
        return result.returncode
    return 0


def write_batch(args):
    if args.mode in {"zip_archive", "grouped_zip_archives"}:
        if not args.filter_suffix and not args.filter_regex:
            sys.stderr.write(
                "WARNING: ZIP mode is running without a suffix or regular-expression "
                "filter, so every regular archive member will be processed.\n"
            )
        members = archive_members(args)
    else:
        members = parse_members(args)
    members = apply_filters(members, args.filter_suffix, args.filter_regex)
    per_taxon = per_taxon_metric(args.metric, args.verbose)
    metric_header = ["group", "file_id", "file_name", "metric", "value", "status"]
    if per_taxon:
        metric_header = ["group", "file_id", "file_name", "metric", "taxon", "value", "status"]
    group_counts = {}

    for member in members:
        group_counts.setdefault(
            member.group,
            {"total": 0, "selected": 0, "processed": 0, "failed": 0, "skipped": 0},
        )
        group_counts[member.group]["total"] += 1
        if member.selected:
            group_counts[member.group]["selected"] += 1
        else:
            group_counts[member.group]["skipped"] += 1

    with open(args.output, "w", newline="") as output_handle, open(
        args.failed_output, "w", newline=""
    ) as failed_handle:
        output_writer = csv.writer(output_handle, delimiter="\t", lineterminator="\n")
        failed_writer = csv.writer(failed_handle, delimiter="\t", lineterminator="\n")
        output_writer.writerow(metric_header)
        failed_writer.writerow(["group", "file_id", "file_name", "metric", "exit_code", "error"])

        for member in members:
            if not member.selected:
                continue
            result = run_phykit(args.metric, member.path, args.verbose)
            if result.returncode != 0:
                group_counts[member.group]["failed"] += 1
                error = clean_error(result.stderr) or clean_error(result.stdout)
                failed_writer.writerow(
                    [member.group, member.file_id, member.name, args.metric, result.returncode, error]
                )
                if per_taxon:
                    output_writer.writerow(
                        [member.group, member.file_id, member.name, args.metric, "", "", "failed"]
                    )
                else:
                    output_writer.writerow(
                        [member.group, member.file_id, member.name, args.metric, "", "failed"]
                    )
                continue

            try:
                if per_taxon:
                    rows = parse_taxon_output(result.stdout)
                    for taxon, value in rows:
                        output_writer.writerow(
                            [member.group, member.file_id, member.name, args.metric, taxon, value, "ok"]
                        )
                else:
                    output_writer.writerow(
                        [
                            member.group,
                            member.file_id,
                            member.name,
                            args.metric,
                            parse_scalar_output(result.stdout),
                            "ok",
                        ]
                    )
            except ValueError as error:
                group_counts[member.group]["failed"] += 1
                failed_writer.writerow(
                    [member.group, member.file_id, member.name, args.metric, "0", str(error)]
                )
                if per_taxon:
                    output_writer.writerow(
                        [member.group, member.file_id, member.name, args.metric, "", "", "failed"]
                    )
                else:
                    output_writer.writerow(
                        [member.group, member.file_id, member.name, args.metric, "", "failed"]
                    )
                continue
            group_counts[member.group]["processed"] += 1

    with open(args.summary_output, "w", newline="") as summary_handle:
        writer = csv.writer(summary_handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            ["group", "total_files", "selected_files", "processed_files", "failed_files", "skipped_files"]
        )
        for group in sorted(group_counts):
            counts = group_counts[group]
            writer.writerow(
                [
                    group,
                    counts["total"],
                    counts["selected"],
                    counts["processed"],
                    counts["failed"],
                    counts["skipped"],
                ]
            )
    return 0


def build_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--metric", required=True)
    parser.add_argument(
        "--mode",
        choices=[
            "single",
            "collection",
            "grouped_collections",
            "zip_archive",
            "grouped_zip_archives",
        ],
        required=True,
    )
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--output", required=True)
    parser.add_argument("--failed-output", required=True)
    parser.add_argument("--summary-output", required=True)
    parser.add_argument("--single-input")
    parser.add_argument("--member", action="append")
    parser.add_argument("--member-name", action="append")
    parser.add_argument("--member-group", action="append")
    parser.add_argument("--collection-label", default="collection")
    parser.add_argument("--archive", action="append")
    parser.add_argument("--archive-name", action="append")
    parser.add_argument("--archive-group", action="append")
    parser.add_argument("--archive-label", default="archive")
    parser.add_argument("--filter-suffix", default="")
    parser.add_argument("--filter-regex", default="")
    return parser


def main():
    args = build_parser().parse_args()
    if args.mode == "single":
        if not args.single_input:
            raise ValueError("Single-file mode requires --single-input")
        return write_single(args)
    return write_batch(args)


if __name__ == "__main__":
    sys.exit(main())
