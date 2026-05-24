# Galaxy Tools

Galaxy wrappers and wrapper prototypes for biomedical data analysis.

This repository is intended to live under `goeckslab` and follow the broad
pattern used by community Galaxy wrapper repositories: one tool directory per
wrapper, Tool Shed metadata in each directory, Planemo lint/test support, and
small committed test data where possible.

## Tool Index

| Tool | Purpose | Dependency strategy |
| --- | --- | --- |
| `featurewise_correlation` | Compute one Spearman or Pearson correlation test per matched feature across two matrices, with multiple-testing correction. | Explicit public BioContainer runtime plus `scipy=1.17.1` requirement. |
| `gseapy_enrichr` | Run GSEApy Enrichr-style over-representation analysis from a Galaxy gene list, named Enrichr libraries, or uploaded GMT files. | `gseapy=1.2.1` package requirement. |
| `kegg_ora` | Run KEGG-style pathway over-representation analysis from foreground/background gene lists and a gene-to-pathway mapping. | `python=3.11` package requirement, resolving to the Python BioContainer. |
| `nonparametric_rank_tests` | Run independent Mann-Whitney U tests and paired Wilcoxon signed-rank tests on tabular data. | Explicit public BioContainer runtime plus `scipy=1.17.1` requirement. |
| `phykit_metrics` | Expose selected PhyKIT tree and alignment metrics not covered by the current public Tool Shed PhyKIT wrappers. | `phykit=2.1.93` package requirement. |
| `rds_to_tabular` | Convert RDS/RData objects containing rectangular R data into Galaxy tabular datasets. | `bioconductor-deseq2=1.42.0` package requirement for a BioContainer-backed R environment. |

## Development

Use the project virtual environment from the repository root:

```bash
.venv/bin/planemo --version
```

Run Planemo lint for every wrapper:

```bash
for tool in tools/*/*.xml; do
    .venv/bin/planemo lint "$tool"
done
```

Run Tool Shed lint for every tool directory:

```bash
for repo in tools/*; do
    .venv/bin/planemo shed_lint "$repo"
done
```

Run a focused wrapper test:

```bash
.venv/bin/planemo test --conda_auto_install tools/kegg_ora/kegg_ora.xml
```

For usegalaxy.org-oriented checks, prefer testing the same container resolution
path Galaxy will use. For example, `kegg_ora` should resolve `python=3.11` to
`quay.io/biocontainers/python:3.11`:

```bash
TMPDIR="$HOME/.tmp/planemo-galaxy-tools" \
    .venv/bin/planemo test --docker tools/kegg_ora/kegg_ora.xml
```

## Current Status

- `featurewise_correlation`: Planemo lint and fixture tests pass locally; the
  wrapper uses an explicit public GSEApy BioContainer runtime that includes
  SciPy, avoiding Conda fallback.
- `gseapy_enrichr`: Planemo lint, Tool Shed lint, and the committed GMT fixture
  test pass locally; dependency resolution can use the existing GSEApy
  BioContainer.
- `kegg_ora`: Planemo lint, Tool Shed lint, direct fixture checks, and Docker
  Planemo tests pass locally; malformed non-empty gene-list rows fail clearly.
- `nonparametric_rank_tests`: Planemo lint, Tool Shed lint, and fixture tests
  pass locally; grouped-table Wilcoxon input is rejected because it has no pair
  identifier, and the wrapper uses an explicit public GSEApy BioContainer
  runtime that includes SciPy, avoiding Conda fallback.
- `phykit_metrics`: Planemo lint, Tool Shed lint, and fixture tests pass
  locally with Galaxy-managed dependencies.
- `rds_to_tabular`: Planemo lint, Tool Shed lint, and fixture tests pass locally;
  list extraction and RData object-name handling have targeted regression
  coverage.
