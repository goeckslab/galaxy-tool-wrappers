# PhyKIT Metrics

This wrapper exposes selected PhyKIT operations through Galaxy for either one
file or many files in Galaxy collections.

Initial operations:

- evolutionary rate
- total tree length
- long branch score
- treeness
- relative composition variability

The existing ToolShed `padge/phykit` wrapper should still be installed for its
covered operations, especially parsimony informative sites, patristic
distances, and treeness over RCV. This wrapper complements those tools by
exposing a small set of direct metric commands with focused tests.

Input modes:

- single file: preserves the direct PhyKIT output for one tree or alignment
- one collection: emits one combined metric table for every selected collection
  member
- grouped collections: emits one long table with group labels, file names,
  metric values, and run status for downstream statistics
- ZIP archive: convenience mode for one archive-heavy input
- grouped ZIP archives: convenience mode for group-labelled archive-heavy
  inputs

Collection modes also write a failed-files table and a per-group summary table
with total, selected, processed, failed, and skipped file counts. Use Galaxy
archive or unzip tools to extract archive members before running this wrapper
when you need to inspect or reuse extracted files. ZIP modes extract local
Galaxy ZIP datasets inside the job, reject unsafe archive paths, and warn when
no suffix or regular-expression filter is supplied.
