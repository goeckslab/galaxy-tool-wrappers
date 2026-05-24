#!/usr/bin/env python3
"""KEGG over-representation analysis for tabular gene lists."""

from __future__ import annotations

import argparse
import csv
import math
import sys
import urllib.request
from collections import defaultdict
from pathlib import Path


MISSING = {"", "na", "nan", "none", "null"}


def clean_gene(value: str) -> str:
    value = value.strip()
    if ":" in value:
        value = value.split(":", 1)[1]
    return value


def read_gene_list(path: Path, column: int, has_header: bool) -> list[str]:
    genes: list[str] = []
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        if has_header:
            next(reader, None)
        for row_number, row in enumerate(reader, start=2 if has_header else 1):
            if not row or all(cell.strip() == "" for cell in row):
                continue
            if len(row) < column:
                raise ValueError(
                    f"Gene-list row {row_number} has {len(row)} columns, "
                    f"but selected gene column is {column}"
                )
            value = clean_gene(row[column - 1])
            if value.lower() in MISSING:
                continue
            genes.append(value)
    if not genes:
        raise ValueError(f"No genes were read from {path}")
    return genes


def request_text(url: str) -> str:
    with urllib.request.urlopen(url, timeout=120) as response:
        return response.read().decode("utf-8")


def read_kegg_rest_mapping(organism: str) -> tuple[dict[str, set[str]], dict[str, str]]:
    organism = organism.strip()
    if not organism:
        raise ValueError("KEGG organism code is required for KEGG REST mode")
    link_text = request_text(f"https://rest.kegg.jp/link/pathway/{organism}")
    list_text = request_text(f"https://rest.kegg.jp/list/pathway/{organism}")

    pathway_to_genes: dict[str, set[str]] = defaultdict(set)
    for line in link_text.splitlines():
        if not line.strip():
            continue
        gene_raw, pathway_raw = line.rstrip("\n").split("\t")[:2]
        pathway = pathway_raw.replace("path:", "")
        pathway_to_genes[pathway].add(clean_gene(gene_raw))

    pathway_names: dict[str, str] = {}
    for line in list_text.splitlines():
        if not line.strip():
            continue
        pathway_raw, name = line.rstrip("\n").split("\t", 1)
        pathway_names[pathway_raw.replace("path:", "")] = name

    return dict(pathway_to_genes), pathway_names


def read_mapping_table(
    path: Path,
    gene_column: int,
    pathway_column: int,
    name_column: int | None,
    has_header: bool,
) -> tuple[dict[str, set[str]], dict[str, str]]:
    pathway_to_genes: dict[str, set[str]] = defaultdict(set)
    pathway_names: dict[str, str] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        if has_header:
            next(reader, None)
        for row_number, row in enumerate(reader, start=2 if has_header else 1):
            if len(row) < max(gene_column, pathway_column, name_column or 1):
                raise ValueError(f"Mapping row {row_number} has too few columns")
            gene = clean_gene(row[gene_column - 1])
            pathway = row[pathway_column - 1].strip().replace("path:", "")
            if gene.lower() in MISSING or pathway.lower() in MISSING:
                continue
            pathway_to_genes[pathway].add(gene)
            if name_column is not None:
                name = row[name_column - 1].strip()
                if name:
                    pathway_names[pathway] = name
    if not pathway_to_genes:
        raise ValueError("No pathway-to-gene mappings were read")
    return dict(pathway_to_genes), pathway_names


def hypergeom_sf(k: int, n: int, K: int, N: int) -> float:
    """P(X >= k) for X ~ Hypergeometric(N, K, n)."""
    if k <= 0:
        return 1.0
    lo = max(k, 0)
    hi = min(K, n)
    denom = math.comb(N, n)
    total = 0
    for x in range(lo, hi + 1):
        if n - x <= N - K:
            total += math.comb(K, x) * math.comb(N - K, n - x)
    return min(1.0, total / denom)


def bh_adjust(p_values: list[float]) -> list[float]:
    m = len(p_values)
    if m == 0:
        return []
    indexed = sorted(enumerate(p_values), key=lambda item: item[1])
    adjusted = [1.0] * m
    running = 1.0
    for rank, (idx, p_value) in reversed(list(enumerate(indexed, start=1))):
        running = min(running, p_value * m / rank)
        adjusted[idx] = min(1.0, running)
    return adjusted


def write_outputs(
    output: Path,
    contributing: Path,
    rows: list[dict[str, object]],
) -> None:
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "pathway_id",
                "pathway_name",
                "foreground_hits",
                "pathway_genes_in_universe",
                "foreground_genes_in_universe",
                "universe_genes",
                "p_value",
                "adjusted_p_value",
                "foreground_hit_genes",
            ]
        )
        for row in rows:
            writer.writerow(
                [
                    row["pathway_id"],
                    row["pathway_name"],
                    row["k"],
                    row["K"],
                    row["n"],
                    row["N"],
                    f"{row['p_value']:.12g}",
                    f"{row['adjusted_p_value']:.12g}",
                    ",".join(row["hit_genes"]),
                ]
            )

    with contributing.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["pathway_id", "pathway_name", "gene_id"])
        for row in rows:
            for gene in row["hit_genes"]:
                writer.writerow([row["pathway_id"], row["pathway_name"], gene])


def run(args: argparse.Namespace) -> None:
    foreground = set(read_gene_list(Path(args.foreground), args.foreground_column, args.foreground_has_header == "true"))
    if args.background:
        background = set(read_gene_list(Path(args.background), args.background_column, args.background_has_header == "true"))
    else:
        background = None

    if args.mapping_source == "kegg_rest":
        pathway_to_genes, pathway_names = read_kegg_rest_mapping(args.organism)
    else:
        pathway_to_genes, pathway_names = read_mapping_table(
            Path(args.mapping_table),
            args.mapping_gene_column,
            args.mapping_pathway_column,
            args.mapping_name_column if args.mapping_name_column > 0 else None,
            args.mapping_has_header == "true",
        )

    all_mapped_genes = set().union(*pathway_to_genes.values())
    universe = all_mapped_genes if background is None else background & all_mapped_genes
    foreground_in_universe = foreground & universe
    N = len(universe)
    n = len(foreground_in_universe)
    if N == 0:
        raise ValueError("The universe has zero genes after intersecting with pathway mappings")
    if n == 0:
        raise ValueError("The foreground has zero genes after intersecting with the universe and pathway mappings")

    rows: list[dict[str, object]] = []
    p_values: list[float] = []
    for pathway_id, raw_genes in sorted(pathway_to_genes.items()):
        pathway_genes = raw_genes & universe
        K = len(pathway_genes)
        hit_genes = sorted(pathway_genes & foreground_in_universe)
        k = len(hit_genes)
        if K == 0:
            continue
        p_value = hypergeom_sf(k, n, K, N)
        p_values.append(p_value)
        rows.append(
            {
                "pathway_id": pathway_id,
                "pathway_name": pathway_names.get(pathway_id, pathway_id),
                "k": k,
                "K": K,
                "n": n,
                "N": N,
                "p_value": p_value,
                "adjusted_p_value": 1.0,
                "hit_genes": hit_genes,
            }
        )

    for row, adjusted in zip(rows, bh_adjust(p_values), strict=True):
        row["adjusted_p_value"] = adjusted

    rows.sort(key=lambda row: (row["adjusted_p_value"], row["p_value"], row["pathway_id"]))
    write_outputs(Path(args.output), Path(args.contributing_genes), rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("--foreground", required=True)
    parser.add_argument("--foreground-column", type=int, required=True)
    parser.add_argument("--foreground-has-header", choices=["true", "false"], required=True)
    parser.add_argument("--background", default="")
    parser.add_argument("--background-column", type=int, default=1)
    parser.add_argument("--background-has-header", choices=["true", "false"], default="true")
    parser.add_argument("--mapping-source", choices=["kegg_rest", "mapping_table"], required=True)
    parser.add_argument("--organism", default="")
    parser.add_argument("--mapping-table", default="")
    parser.add_argument("--mapping-gene-column", type=int, default=1)
    parser.add_argument("--mapping-pathway-column", type=int, default=2)
    parser.add_argument("--mapping-name-column", type=int, default=3)
    parser.add_argument("--mapping-has-header", choices=["true", "false"], default="true")
    parser.add_argument("--output", required=True)
    parser.add_argument("--contributing-genes", required=True)
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    try:
        run(args)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
