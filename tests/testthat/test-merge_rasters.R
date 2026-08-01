context("merge_rasters")
library(terra)
library(elevatr)

# Offline tests: build a set of adjacent tiles on a common grid carrying a
# smooth analytic surface, and check that merge_rasters() mosaics and
# reprojects them without introducing seam artifacts.  No network access is
# required.

# Analytic surface so expected values are known everywhere
surf <- function(x, y){
  1000 * sin(x / 20000) * cos(y / 17000) + 0.002 * x - 0.001 * y
}

# Write a 2x2 block of adjacent 256x256 tiles (EPSG:3857, 100 m pixels) and a
# single pre-merged raster covering the same extent, and return their paths
make_tiles <- function(dir){
  n <- 256
  res <- 100
  crs_merc <- "EPSG:3857"
  x0 <- 0
  y1 <- 2 * n * res  # top of the block

  tile_files <- character(0)
  for(ti in 0:1){
    for(tj in 0:1){
      xmin <- x0 + ti * n * res
      ymax <- y1 - tj * n * res
      r <- terra::rast(nrows = n, ncols = n,
                       xmin = xmin, xmax = xmin + n * res,
                       ymin = ymax - n * res, ymax = ymax,
                       crs = crs_merc)
      xy <- terra::xyFromCell(r, seq_len(terra::ncell(r)))
      terra::values(r) <- surf(xy[, 1], xy[, 2])
      f <- file.path(dir, sprintf("tile_%s_%s.tif", ti, tj))
      terra::writeRaster(r, f, overwrite = TRUE)
      tile_files <- c(tile_files, f)
    }
  }

  big <- terra::rast(nrows = 2 * n, ncols = 2 * n,
                     xmin = x0, xmax = x0 + 2 * n * res,
                     ymin = 0, ymax = y1,
                     crs = crs_merc)
  xy <- terra::xyFromCell(big, seq_len(terra::ncell(big)))
  terra::values(big) <- surf(xy[, 1], xy[, 2])
  big_file <- file.path(dir, "premerged.tif")
  terra::writeRaster(big, big_file, overwrite = TRUE)

  list(tiles = tile_files, premerged = big_file)
}

tmp <- file.path(tempdir(), "merge_rasters_test")
dir.create(tmp, showWarnings = FALSE)
tt <- make_tiles(tmp)

test_that("merge_rasters mosaics and reprojects without seam artifacts", {
  got <- suppressMessages(
    merge_rasters(tt$tiles, target_prj = 4326, tmp_dir = tmp)
  )

  # reference: identical single warp of the pre-merged raster
  ref_file <- file.path(tmp, "reference.tif")
  sf::gdal_utils(util = "warp",
                 source = tt$premerged,
                 destination = ref_file,
                 options = c("-r", "bilinear",
                             "-t_srs", sf::st_crs(4326)$wkt))
  ref <- terra::rast(ref_file)

  expect_s4_class(got, "SpatRaster")
  expect_equal(as.character(terra::crs(got, describe = TRUE)$code), "4326")

  # same grid as the reference warp
  expect_equal(dim(got), dim(ref))
  expect_equal(as.vector(terra::ext(got)), as.vector(terra::ext(ref)),
               tolerance = 1e-8)

  gv <- terra::values(got)
  rv <- terra::values(ref)

  # no NA seams introduced by the mosaic (the reference has none)
  expect_equal(sum(is.na(gv)), sum(is.na(rv)))

  # pixel-identical to warping one contiguous raster: this is the seam test.
  # Warping the tiles as separate sources leaves errors on the seam rows and
  # columns; the VRT mosaic must not.
  expect_equal(max(abs(gv - rv), na.rm = TRUE), 0, tolerance = 1e-6)
})

test_that("merge_rasters works with a single input file", {
  got <- suppressMessages(
    merge_rasters(tt$tiles[1], target_prj = 4326, tmp_dir = tmp)
  )
  expect_s4_class(got, "SpatRaster")
  expect_equal(as.character(terra::crs(got, describe = TRUE)$code), "4326")
})

test_that("merge_rasters returnRaster = FALSE returns a readable path", {
  got <- suppressMessages(
    merge_rasters(tt$tiles, target_prj = 4326, returnRaster = FALSE,
                  tmp_dir = tmp)
  )
  expect_type(got, "character")
  expect_true(file.exists(got))
  expect_s4_class(terra::rast(got), "SpatRaster")
})

test_that("merge_rasters takes the target projection from the tiles when NULL", {
  got <- suppressMessages(
    merge_rasters(tt$tiles, target_prj = NULL, tmp_dir = tmp)
  )
  expect_equal(as.character(terra::crs(got, describe = TRUE)$code), "3857")
  # in-CRS mosaic must reproduce the premerged raster exactly
  ref <- terra::rast(tt$premerged)
  expect_equal(dim(got), dim(ref))
  expect_equal(max(abs(terra::values(got) - terra::values(ref)), na.rm = TRUE),
               0, tolerance = 1e-6)
})
