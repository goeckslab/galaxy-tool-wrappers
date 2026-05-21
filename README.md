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

2. `phykit_extended`
   - Exposes PhyKIT operations missing from the current ToolShed `padge/phykit`
     wrapper, especially long branch score, total tree length, and standalone
     relative composition variability.

3. `mann_whitney_wilcoxon`
   - Modern Mann-Whitney U and Wilcoxon signed-rank wrapper using SciPy.
   - Replaces the need to install the legacy `bebatut/compute_wilcoxon_test`
     wrapper on usegalaxy.org.

## Development

Run lint checks with Planemo from the repository root:

```bash
planemo lint tools/rds_to_tabular/rds_to_tabular.xml
planemo lint tools/phykit_extended/phykit_extended.xml
planemo lint tools/mann_whitney_wilcoxon/mann_whitney_wilcoxon.xml
```

Run committed fixture tests:

```bash
planemo test --no_dependency_resolution tools/rds_to_tabular/rds_to_tabular.xml
planemo test --no_dependency_resolution tools/phykit_extended/phykit_extended.xml
planemo test --conda_auto_install tools/mann_whitney_wilcoxon/mann_whitney_wilcoxon.xml
```

Current local status:

- `rds_to_tabular`: Planemo lint passed; committed fixture test passed.
- `phykit_extended`: Planemo lint passed; committed fixture tests passed;
  `planemo shed_lint` passed.
- `mann_whitney_wilcoxon`: Planemo lint passed; ToolShed `shed_lint` passed;
  committed fixture tests passed with Galaxy-managed Conda dependencies.
