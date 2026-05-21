# RDS To Tabular

Converts R serialized objects to Galaxy tabular datasets.

This wrapper is for `.rds` files, where the correct R API is `readRDS()`. It is
not the same as older `.RData` extraction wrappers that call `load()`.

Supported rectangular outputs:

- `data.frame`
- matrix
- tibble-like objects coercible with `as.data.frame`
- S4 `DataFrame`-like objects coercible with `as.data.frame`
- lists containing a rectangular object when `first rectangular table` is chosen

Rows names can be preserved as the first column when needed.
