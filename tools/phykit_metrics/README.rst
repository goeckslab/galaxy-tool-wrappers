PhyKIT metrics
==============

Expose selected PhyKIT tree and alignment metrics in Galaxy, including
evolutionary rate, total tree length, long branch score, treeness, and relative
composition variability. This wrapper complements existing ToolShed PhyKIT
wrappers with focused direct metric interfaces.

The wrapper supports single-file mode, collection modes, and ZIP convenience
modes. Single-file mode preserves direct PhyKIT output. Collection,
grouped-collection, ZIP, and grouped-ZIP modes emit a combined metric table plus
failed-file and per-group summary tables so many tree or alignment files can be
compared downstream without silently dropping failed PhyKIT runs.

Collection mode remains the preferred Galaxy-native route when users need to
inspect or reuse extracted files. ZIP modes extract local Galaxy ZIP datasets
inside the job, reject unsafe archive paths, and warn when no suffix or
regular-expression filter is supplied.
