# Convert basic sightings spreadsheet into all `dodo` formats

Converts a sightings spreadsheet into all the data formats used in
`dodo` models.

## Usage

``` r
convert_dodo(
  x,
  init.time,
  test.time = as.numeric(format(Sys.Date(), "%Y")),
  threshold = 0.9,
  unique = TRUE,
  aggregate = "pci_prod",
  time,
  certainty,
  certainty_lower,
  certainty_upper
)
```

## Arguments

- x:

  sightings spreadsheet to be converted. Must be a `data.frame` with all
  of the requisite columns included. See Table S1 in Buchanan et
  al. (2025) for an example.

- init.time:

  start of the observation period.

- test.time:

  end of the observation period, typically the present day (defaults to
  the current year).

- threshold:

  cutoff certainty value for which sightings to consider certain (for
  binary certain/uncertain models). Defaults to `0.9`.

- unique:

  whether to deduplicate the dataset (to mitigate nonidependence
  issues). Defaults to `TRUE`, in which case the dataset is grouped by
  every column except `time`, with the most recent row kept for each
  group.

- aggregate:

  what function to use when aggregating records into `iucn` format.
  Default is `"pci_prod"`, which returns \\1 - \prod(1 - x)\\ where
  \\x\\ is the vector of certainty values for each time step. Other
  sensible functions include `"min`, `"mean"`, or `"max"`.

- time:

  name of the column with time values.

- certainty:

  name of the column with certainty values.

- certainty_lower:

  name of the column with certainty lower bound values.

- certainty_upper:

  name of the column with certainty upper bound values.

## Value

a `list` object with 11 elements: the original parameters `init.time`,
`test.time`, `threshold`, and `unique`, and the 7 converted datasets.

## Data Formats

**`ccon`**: *continuous certain sightings*. The simplest format, a
vector of all the (certain) sighting times of the taxon of interest. See
e.g.
[`SO93F1`](https://louis-backstrom.github.io/dodo/reference/SO93F1.md).

**`cbin`**: *binary certain sightings*. A vector of 0/1 (certain)
sighting states at all discrete (integer) time intervals between
`init.time` and `test.time` for the taxon of interest. Similar to
`ccon`, but note that multiple records in an interval collapse to 1. See
e.g.
[`CB14B1`](https://louis-backstrom.github.io/dodo/reference/CB14B1.md).

**`cdis`**: *discrete certain sightings*. A vector of (certain) sighting
counts at all discrete (integer) time intervals between `init.time` and
`test.time` for the taxon of interest. Similar to `cbin`, but note that
multiple records in an interval *do not* collapse to 1. See e.g.
[`BU95F1`](https://louis-backstrom.github.io/dodo/reference/BU95F1.md).

**`ucon`**: *continuous uncertain sightings*. A `data.frame` of all the
sighting times of the taxon of interest, with certainty scores between 0
and 1 associated with each record. See e.g.
[`JR14F1`](https://louis-backstrom.github.io/dodo/reference/JR14F1.md).

**`ubin`**: *binary uncertain sightings*. A `data.frame` with two
columns, the first being the vector of 0/1 certain sighting states at
all discrete (integer) time intervals between `init.time` and
`test.time` for the taxon of interest, and the second being the
equivalent vector for uncertain sightings. See e.g.
[`KO20B2`](https://louis-backstrom.github.io/dodo/reference/KO20B2.md).

**`umcd`**: *multi-class discrete uncertain sightings*. A `data.frame`
of sighting counts at all discrete (integer) time intervals between
`init.time` and `test.time` for the taxon of interest. Sightings are
split up into classes based on their `certainty` score. See e.g.
[`TH13B1`](https://louis-backstrom.github.io/dodo/reference/TH13B1.md).

**`iucn`**: *IUCN-format sightings*. A `data.frame` of sighting
certainty (\\p(ci)\\) bounds (lower and upper) at all discrete (integer)
time intervals between `init.time` and `test.time` for the taxon of
interest. If there are multiple sightings in a single time interval,
they are assumed to be independent of one another and the overall
\\p(ci)\\ score for that period is defined as \\1 - \prod 1 - p(ci)\\.
See e.g.
[`TH17I1`](https://louis-backstrom.github.io/dodo/reference/TH17I1.md).

## Examples

``` r
# Convert the raw Slender-billed Curlew data
convert_dodo(
  x = curlew_raw, init.time = 1817, test.time = 2022, threshold = 0.9,
  unique = TRUE, aggregate = "pci_prod", time = "year", certainty = "p_ci",
  certainty_lower = "p_ci_min", certainty_upper = "p_ci_max"
)
#> Error in if (class(x) != "data.frame") {    stop("Sightings spreadsheet must be a data.frame object!")}: the condition has length > 1
```
