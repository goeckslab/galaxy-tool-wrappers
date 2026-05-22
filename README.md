# Galaxy Tools

Galaxy wrappers and wrapper prototypes for biomedical data analysis.

This repository is intended to live under `goeckslab` and follow the same broad
pattern as community Galaxy wrapper repositories: one tool directory per wrapper,
Planemo lint/test support, and small test data where possible.

## Current Priorities

1. `rds_to_tabular`
   - Converts `.rds` objects into Galaxy tabular datasets using R `readRDS()`.
   - Supports common rectangular R objects such as data frames, matrices,
     Bioconductor DataFrame-like objects, and lists containing a rectangular
     table.

2. `phykit_metrics`
   - Exposes PhyKIT operations missing from the current ToolShed `padge/phykit`
     wrapper, especially long branch score, total tree length, and standalone
     relative composition variability.

3. `nonparametric_rank_tests`
   - Modern rank-test wrapper using SciPy; currently supports Mann-Whitney U
     and Wilcoxon signed-rank tests.
   - Replaces the need to install the legacy `bebatut/compute_wilcoxon_test`
     wrapper on usegalaxy.org.

4. `gseapy_enrichr`
   - Runs GSEApy Enrichr over-representation analysis from a Galaxy gene list.
   - Supports named Enrichr libraries or uploaded GMT files, standardized
     result tables, ranked top terms, and term-substring summaries.

## Development

Run lint checks with Planemo from the repository root:

```bash
planemo lint tools/rds_to_tabular/rds_to_tabular.xml
planemo lint tools/phykit_metrics/phykit_metrics.xml
planemo lint tools/nonparametric_rank_tests/nonparametric_rank_tests.xml
planemo lint tools/gseapy_enrichr/gseapy_enrichr.xml
```

Run committed fixture tests:

```bash
planemo test --no_dependency_resolution tools/rds_to_tabular/rds_to_tabular.xml
planemo test --conda_auto_install tools/phykit_metrics/phykit_metrics.xml
planemo test --conda_auto_install tools/nonparametric_rank_tests/nonparametric_rank_tests.xml
planemo test --conda_auto_install tools/gseapy_enrichr/gseapy_enrichr.xml
```

Current local status:

- `rds_to_tabular`: Planemo lint passed; committed fixture test passed.
- `phykit_metrics`: Planemo lint passed; committed fixture tests passed with
  Galaxy-managed Conda dependencies; `planemo shed_lint` passed.
- `nonparametric_rank_tests`: Planemo lint passed; ToolShed `shed_lint` passed;
  committed fixture tests passed with Galaxy-managed Conda dependencies.
- `gseapy_enrichr`: Planemo lint passed; ToolShed `shed_lint` passed;
  committed GMT fixture test passed with Galaxy-managed Conda dependencies.
