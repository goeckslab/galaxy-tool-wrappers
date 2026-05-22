GSEApy Enrichr
==============

Galaxy wrapper for GSEApy's Enrichr over-representation analysis.

The wrapper accepts a gene list and either a named Enrichr library or an
uploaded GMT file. Outputs include a standardized enrichment table, a ranked
top-terms table, a simple term-substring summary, and JSON metadata.

Use a named Enrichr library when the Galaxy job environment can reach the
Enrichr service. Use an uploaded GMT file when reproducibility or offline
execution is required.
