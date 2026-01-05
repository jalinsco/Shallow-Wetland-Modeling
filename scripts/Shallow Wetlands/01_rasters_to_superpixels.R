#############################################X
#----------- Prepare Superpixels ------------X
#############################################X
# ----- RGB rasters to l*A*B superpixels ----X
#############################################X



# Load packages ----------------------

# required R packages
required_pkgs <- c("terra", 
                   "farver",
                   "supercells",
                   "future",
                   "future.apply")

# load 
missing <- required_pkgs[!vapply(required_pkgs, 
                                 requireNamespace, 
                                 logical(1), 
                                 quietly = TRUE)]

# check if missing
if (length(missing) > 0) {
  stop("Missing packages: ", paste(missing, collapse = ", "))
}


# Load raster list  ----------------------

# code designed for 4-band Planet-Scope surface reflectance rasters
# here, rasters cropped to small tiles
# tile tifs not included in this repo 

# load list of tiled rasters to process
# files <- list.files(paste0(w, "Planet Spring Tiles/"), full.names=TRUE)


# Convert RGB to l*A*B  ----------------------

# function to covert
lab_rasters = function(x){
  
  # load
  r <- rast(x)
  
  # process 
  r <- clamp(r, lower=0, upper=10000, values=TRUE)
  NAflag(r) <- 0
  r = r/10000 # scale
  names(r) <- c("blue", "green", "red", "nir")
  
  # transform to l*A*B colorspace
  r <- terra::subset(r, c(3, 2, 1))
  vals <- values(r)
  vals <- vals*255
  new_vals <- convert_colour(vals, from = "rgb", to = "lab")
  values(r) <- new_vals
  names(r) <- c("l", "a", "b")
  
  # save to new file
  out_name = ('path/to/LAB.tif')
  writeRaster(r, filename=out_name) #overwrite=TRUE
  print(paste0(tn, " is finished"))
}

# apply function sequentially over the list
lapply(files, lab_rasters)


# Run SLIC algorithm ----------------------


# function to run
slic_subtiles <- function(x) {
  
  # create extent objects
  sub_ext <- ext(as.numeric(x[3:6]))
  sub_ext_sf <- st_bbox(sub_ext) %>% st_as_sfc() %>% st_set_crs(32614)
  
  # load raster
  r <- rast('path/to/LAB.tif')
  
  msg <- tryCatch({
    
    # segment
    result <- supercells(r, 
                         step = 10, 
                         compactness = 2, 
                         clean = TRUE, 
                         chunks = FALSE, 
                         verbose = 1)
    
    # clean up
    rm(r)
    gc()
    
    # load surface water (dswe) polygons
    # (must update as needed)
    polys <- st_read('path/to/dswe_polygons.gpkg')
    
    # filter dswe supercells
    mat <- st_intersects(result, polys, sparse = TRUE)
    idx <- apply(mat, 1, any)
    result <- result[idx,]
    
    # clean up
    rm(mat, idx, polys)
    gc()
    
    # save
    # (must update as needed)
    st_write(result, 'path/to/superpixels.gpkg')
    
    # clean up environment
    rm(result, e, idx)
    gc()
    
    # report success
    return(paste0(t, " ", s, " for ", yr, " is finished."))
    
    # catch & report errors  
  }, error = function(e) {
    
    return(paste0(t, " ", s, " for ", yr, " failed. Ruh-roh."))
    
  })
  
  return(msg)
  
}

# run in parallel
plan(multisession, workers = 2)
results <- future_apply(sub_df, 1, slic_subtiles, future.seed = NULL)

# code to process a full study area is available on request!
