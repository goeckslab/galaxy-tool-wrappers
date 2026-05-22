Row-wise association
====================

Galaxy wrapper for feature-wise association tests across two tabular matrices.

The tool aligns observations by an identifier column, aligns features by column
name, computes one correlation test per feature, drops missing values within
each feature pair, applies Benjamini-Hochberg correction, and reports both a
per-feature table and a summary table.

This is useful when two assays are measured across the same samples or models
and the scientific question is one test per matched feature.
