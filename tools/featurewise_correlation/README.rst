Feature-wise Correlation Tests
==============================

Galaxy wrapper for per-feature correlations between matched tabular matrices.

The tool aligns observations by an identifier column, aligns features by column
name for matrices with headers, computes one Pearson or Spearman correlation
test per feature, drops missing values within each feature pair, applies
Benjamini-Hochberg correction, and reports both a per-feature table and a
summary table.

This is useful when two assays are measured across the same samples or models
and the scientific question is one test per matched feature.
