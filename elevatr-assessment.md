# elevatr: assessment and a GDAL-native upgrade pathway

Assessment of `USEPA/elevatr` at commit HEAD of the default branch (version 0.99.1,
DESCRIPTION dated after the 2025-09-09 NEWS entry). Source read in full:
`R/get_elev_raster.R`, `R/internal.R`, `R/get_elev_point.R`, `R/zzz.R`. The open
issue tracker and open PRs were reviewed and are cross-referenced in section 9.

## 1. Summary

elevatr's raster path is built around "enumerate slippy tiles, download every tile
with an HTTP client, mosaic the files, then reproject". Every one of those four
steps is something GDAL already does natively, lazily, and in one pass. The result
is a package that:

- issues roughly two HTTP round trips per tile, the first of them in a serial loop
  before any data is fetched;
- downloads the full tile grid at a user-chosen zoom regardless of the resolution
  the user actually wants;
- writes the entire mosaic to disk twice (an intermediate GeoTIFF, then the warped
  output) before returning;
- exposes zoom level `z` as the only resolution control, with no way to request a
  target extent, dimension and CRS directly;
- carries `future`, `furrr`, `purrr`, `progressr`, `slippymath`, `httr` and `raster`
  to reimplement machinery that `gdalwarp` contains.

The core recommendation is to stop treating the AWS terrain tiles as a set of files
to download and start treating them as a single virtual raster that GDAL can warp.
This is not a rewrite of the package's public interface; it is a replacement of the
private fetch layer, plus new optional arguments on the existing entry point.

## 2. What the current pipeline costs

Concrete cost for a real request, computed from the tile arithmetic in
`get_tilexy()` (assuming 256 x 256 tiles at 4 bytes per pixel):

| AOI | tiles | HTTP requests | pixels | payload |
|---|---|---|---|---|
| Tasmania, z=11 | 600 | 1200 | 39.3 M | 150 MB |
| Tasmania, z=12 | 2350 | 4700 | 154.0 M | 588 MB |
| 0.1 deg catchment, z=13 | 12 | 24 | 0.8 M | 3 MB |
| Svalbard 78N, 4 deg, z=11 | 1320 | 2640 | 86.5 M | 330 MB |

Correction after checking the Tilezen format documentation: the `geotiff`
endpoint that elevatr uses serves 512 x 512 tiles ("optimized to reduce transfer
costs in 512x512 tile sizes but with internal 256x256 image pyramiding", Web
Mercator). Raw pixel counts are therefore 4x the table above - Tasmania z=12 is
about 616 M pixels uncompressed - and the effective ground resolution at a given
`z` is one power of two finer than the 256-px arithmetic suggests. This makes the
over-download problem worse than the table shows, and it also means
`estimate_raster_size()`'s degree-based arithmetic (calibrated to 256-px tiles) is
off by an extra factor of 4 on top of the latitude error described below.

The request count doubles because of this loop in `get_aws_terrain()`:

```r
for(i in urls){
  if(httr::http_error(i)) {
    message("An AWS URL is invalid.")
    return(NULL)
  }
}
```

That is one extra serial round trip per tile before a single byte of elevation data
is fetched. At 60 ms latency the Tasmania z=12 case spends over two minutes on
liveness checks alone, and the check is not even useful: it aborts the whole request
if any single tile 404s, which is exactly what happens over ocean where the tile set
is sparse. The subsequent download loop already errors on a non-tif response, so the
pre-flight loop adds latency and a failure mode without adding safety.

Peak disk for the z=12 case is also worth noting: 2350 tile files, plus the
intermediate mosaic GeoTIFF, plus the final warped GeoTIFF. Over 1.2 GB on disk for
a 588 MB payload.

### The resolution mismatch

This is the point that matters most for usability. Suppose the user wants a 512 x 512
grid over Tasmania. There is no way to ask for that. They must choose `z`, and the
choice determines the download:

| z | ground res at 42 S | tiles | payload |
|---|---|---|---|
| 7 | 909 m | 6 | 1.5 MB |
| 8 | 454 m | 16 | 4 MB |
| 10 | 114 m | 156 | 39 MB |
| 12 | 28 m | 2350 | 588 MB |

A 512 x 512 grid over that extent is about 660 m per pixel, so z=8 is the correct
source level. A user who picks z=12 because they want "good data" downloads 588 MB
and then throws away 99.7 percent of it in the resample. There is nothing in the API
or documentation that connects the zoom argument to the grid they are going to end
up with, so this happens routinely.

### The download size estimator is wrong at high latitude

`estimate_raster_size()` computes pixel counts by dividing the bounding box in
*degrees* by a constant degrees-per-pixel derived from the equatorial ground
resolution. Web Mercator tile resolution varies with latitude, so the estimate is
increasingly wrong away from the equator. The Svalbard row above: the function
reports 65 MB, the actual payload is 330 MB, a 5x under-estimate. It is also more
work than necessary, since `get_tilexy()` has already computed the exact tile grid;
`nrow(tiles) * 256L * 256L * 4L` is exact and free.

## 3. The two-pass warp, and why the VRT fixes it properly

`merge_rasters()` runs `gdalwarp` twice:

```r
sf::gdal_utils("warp", source = files, destination = destfile,
               options = c("-r", method))
# Using two steps now as gdal with one step introduced NA's along seams
# Slower but more accurate!
sf::gdal_utils("warp", source = destfile, destination = destfile2,
               options = c("-r", method, "-t_srs", sf::st_crs(target_prj)$wkt))
```

The comment identifies a real bug. When `gdalwarp` is given N separate source files
it warps each one independently into the destination; at a tile edge the resampling
kernel has no neighbouring pixels, so edge cells are computed from a truncated
kernel or dropped to nodata. That is the seam artefact.

I reproduced this offline. Four adjacent 256 x 256 Float32 tiles on a common
EPSG:3857 grid, carrying a smooth analytic surface, reprojected to EPSG:4326 two
ways and compared against the analytic truth evaluated on the destination grid:

```
per-tile warp        max|err| = 0.021 m   RMSE = 0.0007 m
VRT then one warp    max|err| = 0.000 m   RMSE = 0.0000 m

cells where the two differ by more than 1 cm: 412
all 412 lie on row 255 - the internal tile seam
```

The error is exactly the seam, and building a VRT over the tiles first removes it
exactly, because the warper then sees one contiguous dataset and the kernel has its
neighbours. On a smooth synthetic surface the magnitude is small; on real terrain
with steep gradients, and with the tiles' nodata handling in play, this is the NA
seam the maintainer observed.

So the fix that was reached for (a second full warp) is the expensive workaround,
and the cheap correct fix (`gdalbuildvrt`, an XML file of a few KB, zero pixel I/O)
was available. Replacing pass one with a VRT is a self-contained change to
`merge_rasters()` that removes a full-resolution disk write and read per call and is
strictly more accurate. **If nothing else in this document is adopted, this one
change is worth making.**

## 4. The architectural fix: the tile set is already one GDAL dataset

GDAL's WMS driver has a TMS minidriver that addresses a slippy tile pyramid as a
single raster with a full overview stack. I verified this against GDAL 3.10.3:

```xml
<GDAL_WMS>
  <Service name="TMS">
    <ServerUrl>https://s3.amazonaws.com/elevation-tiles-prod/geotiff/${z}/${x}/${y}.tif</ServerUrl>
  </Service>
  <DataWindow>
    <UpperLeftX>-20037508.34</UpperLeftX>
    <UpperLeftY>20037508.34</UpperLeftY>
    <LowerRightX>20037508.34</LowerRightX>
    <LowerRightY>-20037508.34</LowerRightY>
    <TileLevel>14</TileLevel>
    <TileCountX>1</TileCountX>
    <TileCountY>1</TileCountY>
    <YOrigin>top</YOrigin>
  </DataWindow>
  <Projection>EPSG:3857</Projection>
  <BlockSizeX>256</BlockSizeX>
  <BlockSizeY>256</BlockSizeY>
  <BandsCount>1</BandsCount>
  <DataType>Int16</DataType>
  <ZeroBlockHttpCodes>403,404</ZeroBlockHttpCodes>
  <ZeroBlockOnServerException>true</ZeroBlockOnServerException>
  <Cache><Path>~/.cache/elevatr/aws</Path></Cache>
</GDAL_WMS>
```

Opening that file gives:

```
driver           : WMS
size             : 4194304 x 4194304
crs              : EPSG:3857
dtype            : int16, block 256 x 256
overview factors : 2 4 8 16 32 64 128 256 512 1024 2048 4096 8192 16384
```

Fourteen overviews, one per zoom level from 14 down to 0. Everything follows from
this:

- **`z` becomes an implementation detail.** `gdalwarp` selects the overview level
  from the requested output resolution (`-ovr AUTO` is the default). The user asks
  for a grid; GDAL picks the zoom. The 588 MB accident becomes structurally
  impossible.
- **Fetching is lazy and block-level.** Only the 256 x 256 blocks the warper touches
  are requested, at the chosen overview. No tile enumeration, no tempfile per tile,
  no mosaic step.
- **Caching is free.** The `<Cache>` element gives a persistent on-disk block cache
  across calls and across sessions. Repeated work in the same area stops hitting the
  network at all.
- **Missing tiles stop being fatal.** `ZeroBlockHttpCodes` turns a 404 over ocean
  into a zero block instead of aborting the request.
- **Parallelism comes from GDAL.** `-multi -wo NUM_THREADS=ALL_CPUS` plus
  `GDAL_HTTP_MULTIPLEX=YES` and `GDAL_HTTP_VERSION=2` replaces the whole
  `future`/`furrr`/`progressr` layer, and does not clobber the user's `future::plan()`.

Two caveats on the XML. I could not reach `s3.amazonaws.com` from this sandbox
(egress proxy returns 403), so this is validated structurally against GDAL, not
against live tiles; `DataType` in particular should be confirmed against the
`geotiff` endpoint rather than taken from the block above. And the AWS terrain tile
set is sparse over ocean, which is what `ZeroBlockOnServerException` is there for.

The equivalent for a VRT-of-`/vsicurl/` sources is a reasonable fallback if the WMS
driver turns out to be awkward, but it needs one header read per tile and gives up
the overview pyramid, so it is strictly worse.

One real trade-off to flag: the current download loop harvests the
`x-amz-meta-x-imagery-sources` response header from each tile and attaches the
composited source datasets as a `sources` attribute on the result. A GDAL-side
fetch does not surface per-tile HTTP headers, and open issue #109 asks for exactly
this attribution information to be made *more* accessible, not less. The fix is
cheap - the header is still fetchable with one HEAD request per intersecting tile,
no pixel download needed, so a small `get_elev_sources()` helper (which is roughly
what #109 asks for anyway) preserves and improves the feature - but it should be
planned, not discovered after the header path is deleted.

## 5. Proposed interface

The GDAL warp arguments are the natural API here, because they are what users
actually want to specify: an extent, a shape, and a CRS.

```r
get_elev_raster(
  locations,            # unchanged: used to derive extent when extent is NULL
  extent    = NULL,     # c(xmin, xmax, ymin, ymax) or an sf bbox
  dimension = NULL,     # c(ncol, nrow)
  res       = NULL,     # target resolution in crs units; alternative to dimension
  crs       = NULL,     # target crs; defaults to crs(locations)
  z         = NULL,     # retained, now a resolution hint mapped to an overview
  resample  = "bilinear",
  src       = "aws",
  ...
)
```

Any two of `extent`, `dimension`, `res` determine the grid; `z` is kept for
back-compatibility and simply maps to a source overview. `clip = "bbox"` becomes a
no-op because the extent *is* the request, and `clip = "locations"` reduces to a
`-cutline` on the warp rather than a post-hoc `crop` and `mask` on a fully
materialised raster.

Sketch of the fetch layer, using only `sf::gdal_utils()`, which the package already
depends on:

```r
aws_tms_xml <- function(cache_dir = tools::R_user_dir("elevatr", "cache")) {
  f <- file.path(tempdir(), "elevatr_aws_tms.xml")
  writeLines(sprintf(TMS_TEMPLATE, normalizePath(cache_dir, mustWork = FALSE)), f)
  f
}

warp_dem <- function(extent, dimension, crs, resample = "bilinear",
                     src_xml = aws_tms_xml(), dst = tempfile(fileext = ".tif")) {
  opts <- c(
    "-t_srs", sf::st_crs(crs)$wkt,
    "-te", format(extent[c(1, 3, 2, 4)], scientific = FALSE),
    "-ts", format(dimension),
    "-r",  resample,
    "-of", "GTiff",
    "-co", "TILED=YES", "-co", "COMPRESS=DEFLATE",
    "-multi", "-wo", "NUM_THREADS=ALL_CPUS",
    "-ovr", "AUTO"
  )
  sf::gdal_utils("warp", source = src_xml, destination = dst, options = opts)
  terra::rast(dst)
}
```

That is the whole raster path. No tile enumeration, no `httr`, no `furrr`, no
mosaic, no clip, no second warp.

For `get_elev_point()` with `src = "aws"`, the same virtual dataset removes the
"download a DEM then extract" detour entirely: open the TMS XML with
`terra::rast()` and `terra::extract()` the points. GDAL reads only the blocks
containing points, which for scattered points is a small fraction of what the
current code downloads. `gdallocationinfo` semantics via `vapour::vapour_raster_read`
or `gdalraster` would be equivalent if a lighter path is preferred.

## 6. Suggested staging

Each stage is independently shippable and the first two are invisible to users.

**Stage 1 - VRT instead of the first warp pass.** Change `merge_rasters()` to build
a VRT over the downloaded tiles and warp once. Removes one full-resolution write and
read per call, and fixes the seam problem properly rather than papering over it.
Contained, low risk, immediately measurable.

**Stage 2 - remove the pre-flight loop and fix the estimator.** Delete the
`for(i in urls) httr::http_error(i)` loop; the download path already handles bad
responses. Replace `estimate_raster_size()` with the exact tile-count arithmetic.
Halves the request count and makes the size warning correct at high latitude. Also
remove the dead `dir <- tempdir()` in `get_aws_terrain()`, guard the `readline()`
size prompt with `interactive()` (in a script it currently returns `""` and the
function returns `NULL` silently), and restore rather than reset the user's
`future::plan()` on exit.

**Stage 3 - TMS source behind the existing interface.** Add the WMS/TMS XML as an
internal source and route `src = "aws"` through a single warp. Keep `z` working by
mapping it to an overview level, so no user code changes. This is where the
`future`/`furrr`/`purrr`/`progressr`/`slippymath`/`httr` dependencies can come out.

**Stage 4 - expose the grid.** Add `extent`, `dimension`, `res`, `crs` and
`resample`. Deprecate `clip` in favour of `extent` plus an optional cutline. Switch
the return value to `SpatRaster` and drop `raster`, which the startup message has
been promising since 0.4.5.

**Stage 5 - unify point extraction.** Route `get_elev_point(src = "aws")` through
`terra::extract()` on the virtual dataset.

## 7. Smaller items noted while reading

- `zzz.R` announces "elevatr v0.99.0" from a 0.99.1 package; the string is
  hard-coded rather than read from `DESCRIPTION`.
- `TODO.md` still lists "change prj to crs everywhere". Worth doing in the same
  release as the interface change in stage 4, with `prj` accepted and deprecated.
- `neg_to_na` does `raster_elev[raster_elev < 0] <- NA` on a `SpatRaster`, which
  materialises the whole grid in memory. `terra::classify()` or a warp-time
  `-srcnodata` would stream.
- `get_opentopo()` builds its URL by string concatenation with an unencoded API key.
  `httr2::url_modify_query()` or `utils::URLencode()` would be safer, and the key
  currently ends up in the URL that `httr::http_error()` logs on failure.
- `set_opentopo_key()` writes to `~/.Renviron` via `cat(..., append = TRUE)` with no
  newline check, so it can concatenate onto an existing unterminated final line.
- The `locations` argument does a lot of work: `data.frame`, `sf`, `sfc`, `sfg`,
  `SpatRaster`, `SpatVector`, all normalised to an sf POINT set even when only the
  bounding box is wanted. For the raster path only `sf::st_bbox()` and a CRS are
  needed, so `loc_check()` could be bypassed entirely there.

## 8. Prior art worth pointing the maintainer at

Several R packages already implement the "extent, dimension, CRS, one warp" pattern
against remote tiled sources and are worth reading before committing to a design:
`ceramic` (slippy tile fetch with a persistent cache), `whatarelief` (elevation by
extent/dimension/CRS with no zoom argument in the interface), `vapour` (thin GDAL
warp and read bindings), and `gdalraster` (direct GDAL VSI and warp access, on
CRAN). `terra::project()` with a target template grid covers much of the same ground
for the final step.

## 9. Cross-check against the issue tracker

The repository has 32 open issues and 4 open PRs. The recent ones line up with this
assessment closely enough to be worth citing when proposing the changes.

**Issues the GDAL-native path resolves structurally:**

- **#105, #106, #107 (Bengtsson, Feb 2026)** - three issues from the author of the
  `future` package: use `parallelly::availableCores()`, cap worker counts properly,
  and above all undo changes to `future::plan()` (elevatr currently sets the plan
  and resets it to `sequential`, clobbering whatever the user had; #107 points at
  the exact lines and supplies the `with(plan(...), local = TRUE)` idiom). Stage 2
  fixes this in place; stage 3 makes all three issues moot by removing the
  `future`/`furrr` layer entirely - parallelism moves into the warper, where it
  does not interact with the user's session state.
- **#104 (courtiol, Jan 2026)** - `clip = "bbox"` broken (degenerate LinearRing from
  a two-point bbox polygon, behaviour shifted under newer sf/terra). Symptomatic of
  the clip-after-materialise design: when the extent is part of the warp request
  (stage 4), this entire class of geometry-juggling bug disappears.
- **#102 (courtiol, Sep 2024)** - `z = 15` fails with "An AWS URL is invalid"
  followed by `attempt to set an attribute on NULL`. That is the pre-flight
  `httr::http_error()` loop identified in section 2 returning `NULL` and the caller
  not checking. Stage 2 removes the loop; stage 3 removes user-facing `z`
  validation as a failure mode altogether (an out-of-range zoom cannot be
  requested when the overview level is derived from the target grid).
- **#101 (courtiol, Sep 2024)** - elevatr unfriendly outside the US: the default
  point source (epqs) is US-only and fails slowly and cryptically for the rest of
  the world. The maintainer agreed but was wary of switching defaults. Stage 5 is
  relevant here: with point extraction reading blocks lazily from the global tile
  set, `src = "aws"` becomes cheap enough to be a defensible worldwide default,
  with epqs retained for US users who want the authoritative value.
- **#109 (rdjmacklin, Jun 2026)** - make sourcing/attribution accessible. See the
  trade-off note in section 4; a `get_elev_sources()` helper serves this and
  survives the architecture change.

**Direct support for the stage 4 interface:**

- **PR #108 (rhijmans, Mar 2026, open)** - the author of `terra` has an open PR
  titled "consider SpatExtent": accept a `SpatExtent` as `locations`, plus a fix
  for `crs` being ignored for `SpatVector` input. That is the extent-first
  interface of section 5 arriving from the outside. Merging it (or superseding it
  with the fuller `extent`/`dimension`/`res`/`crs` signature) would be pushing on
  an open door.
- **PR #99 (elipousson, Aug 2024, open)** and **#93** - `get_elev_profile()`,
  elevation along a line. Trivially served by the virtual-dataset path:
  densify the line, `terra::extract()` on the TMS dataset, no DEM download.

**One caution the tracker adds:**

- **#103 (psychemedia, Feb 2025)** - elevatr fails in webR/ShinyLive because
  `httr`/`curl` are unavailable there; maintainer is open to a PR using
  `download.file()`. Worth knowing that the GDAL-native path does not help this
  environment either: sf/terra in webR are built without GDAL network drivers, so
  neither `/vsicurl/` nor the WMS driver works in wasm. If webR support matters, a
  minimal `download.file()` fallback for the tile fetch would need to survive as a
  side path. This is the one open issue that argues for keeping a plain-HTTP code
  path around.

## 10. What remains unverified

- No live requests to `s3.amazonaws.com` or `portal.opentopography.org`: the
  sandbox egress policy denies both hosts (CONNECT 403), and that denial is not
  something to route around. Payload figures are arithmetic from tile counts, not
  measured transfers, and the `DataType` in the TMS XML still needs confirming
  against a real tile from the `geotiff` endpoint (the Tilezen docs confirm
  512 x 512 Web Mercator but do not state the pixel type on the page checked).
- The TMS XML was validated by opening it in GDAL 3.10.3 (driver, size, CRS, block
  size, full overview stack) but not by reading live pixels; sparse-ocean
  `ZeroBlockHttpCodes` behaviour is likewise untested against the real endpoint.
- The seam reproduction used synthetic tiles via rasterio's `reproject` (the same
  underlying warper as `gdalwarp`), not a byte-for-byte rerun of elevatr's
  `sf::gdal_utils()` call.
- Issue review covered the 12 most recent open issues, all open PRs, and the
  individual threads cited above; the 20 older open issues were not read
  (the paginated list is blocked from this environment).
