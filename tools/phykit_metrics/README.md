# PhyKIT Metrics

This wrapper exposes PhyKIT operations that are present in the PhyKIT command
line program but are not exposed by the current ToolShed `padge/phykit` wrapper.

Initial operations:

- total tree length
- long branch score
- treeness
- relative composition variability

The existing ToolShed `padge/phykit` wrapper should still be installed for its
covered operations, especially parsimony informative sites, evolutionary rate,
patristic distances, and treeness over RCV. This wrapper complements those
tools by exposing metrics that are absent from the published wrapper XMLs.
