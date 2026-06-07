KEGG ORA
========

Run KEGG-style pathway over-representation analysis from tabular gene lists and
a gene-to-pathway mapping table, or from KEGG REST organism mappings when
network access is available.

The wrapper supports reproducible supplied or pinned mapping tables, explicit
background universes, and selectable Benjamini-Hochberg adjustment families:
all pathways in the universe or only pathways with foreground hits. A summary
output records the mapping source, mapping hash, universe size, foreground size,
pathway counts, and adjustment scope used for the run.
