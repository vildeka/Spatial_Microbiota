Load Spatial data
================
5/13/24

### Load libraries

``` r
##################
# LOAD LIBRARIES #
##################
library(tidyverse)
library(tidyseurat)
library(Seurat)
library(hdf5r)

library(xml2) # loads the image 
library(sp)
library(terra)
library(cowplot)

source("../bin/help_functions.R")
source("../bin/spatial_visualization.R")

#################
# COLOR PALETTS #
#################
pal <- rep(c(RColorBrewer::brewer.pal(9,"Set1"),
         RColorBrewer::brewer.pal(9,"Pastel1"),
         RColorBrewer::brewer.pal(8,"Accent"),
         RColorBrewer::brewer.pal(8,"Set2"),
         RColorBrewer::brewer.pal(8,"Pastel2") ,
         scales::hue_pal()(8)),99)
```

### Load Visium data

``` r
#########
# PATHS #
#########
input_dir <- "../data/spatial_data"
result_dir <- "../results/00_load_st_data/"
if( isFALSE(dir.exists(result_dir)) ) { dir.create(result_dir,recursive = TRUE) }

#############
# LODA DATA #
#############
meta <- read_csv("../data/ST-samples_metadata.csv")
h5_files <- list.dirs(path = input_dir,
                      full.names = T, recursive = T) %>%
            grep("P\\d\\d\\d$", ., value = TRUE) %>%
            set_names(., str_extract(., "P\\d\\d\\d$"))
sample_id <- c("P020","P045", "P050","P057","P008", "P026", "P031", "P044", "P080", "P105", "P001", "P004", "P014", "P018", "P087", "P108", "P118","P021", "P024", "P067", "P081", "P117" ) %>% set_names()

h5_files <- h5_files[names(h5_files) %in% sample_id]
# the new version of seurat v5 does not allow to read in any other than low res img
# to circumvent this rename the hires to low res and move the lowres file
image <- map(h5_files, 
             ~Read10X_Image(
               filter.matrix = T,
               image.dir = paste0(.x, "/spatial"),
               image.name = "tissue_hires_image.png"))

# Read in h5 files and create Seurat Object
seuratObj_list <- pmap(list(h5_files, image, names(h5_files)),
                       ~Load10X_Spatial(
                         filename = "filtered_feature_bc_matrix.h5",
                         filter.matrix = T,
                         assay = "RNA",
                         data.dir = ..1,
                         image =  ..2,
                         slice = ..3)) 

sample_id <- names(h5_files) %>% set_names()
```

### Tidy up the seurat object

``` r
####################################
# RENAME SAMPLES, SPOTS AND IMAGES #
####################################
seuratObj_list <- seuratObj_list %>%
  imap(., ~AddMetaData(object = .x, 
                       metadata = rep(.y, length(Idents(.x))), 
                       col.name = "orig.ident")) %>%
  map(.,  ~SetIdent(., value = .@meta.data$orig.ident)) %>%
  imap(., ~RenameCells(.x, 
                      new.names = paste0(.y,"_", gsub("-.*","",colnames(.x[["RNA"]])))) ) %>%
  imap(., ~ {.x@images <- set_names(.@images,.y); .x})

##################
# MERGE SAMPLES #
#################
# Merge datasets into one single seurat object

DATA  <- merge(seuratObj_list[[1]], y = seuratObj_list[2:length(seuratObj_list)])
# DATA <- JoinLayers(DATA) # Seurat version 5.0
```

``` r
##########################
# ADD MANUAL ANNOTATION #
#########################
### Load morphology annotation
a2 <- map(sample_id, ~xml2::as_list( read_xml( paste0( input_dir,"/",.x,"/",.x,".svg")) ) )
dd <- map(a2, ~as.numeric( strsplit( attributes(.x$svg)$viewBox , " ")[[1]] ))
dd2 <- map(sample_id, ~dim(DATA@images[[.x]]@image))

# sample_id <- "P097"
# a2 <- a2[["P097"]]
# dd <- dd[["P097"]]
# dd2 <- dd2[["P097"]]

# add image coordinates to the seurat object
get_img_coord <- function(DATA, sample_id){
  img_coord <- DATA@images[[sample_id]]@coordinates
  DATA@images[[sample_id]]@coordinates <<- img_coord 
}

# add manual spatial annotation
get_sp_annot <- function(a2, dd, dd2, sample_id){
  scale.factor <- DATA@images[[sample_id]]@scale.factors$hires
  img_coord <- DATA@images[[sample_id]]@coordinates
  
  id <- a2$svg %>%
    map_chr(., ~attr(.x,"id")) %>% 
    set_names(seq_along(.), .)
  
  annot_coord <- id %>%
    map(., ~get_shape(a2$svg[[.x]]) ) %>%
    map(., ~as_tibble(.x, .name_repair="unique"))  %>%
    # convert the view box from the svg to the pixel dimensions of the raster img:
    map(., ~mutate(.x, x = .x[[1]]*dd2[2]/dd[3],
                       y = .x[[2]]*dd2[1]/dd[4] )) %>%
    map(., ~rowid_to_column(., var = "path_idx")) %>%
    {. ->>  temp} %>%
    bind_rows(., .id="name") %>%
    group_by(., name) %>%
    mutate(elem_idx = cur_group_id()) %>% # group_indices(., name)
    ungroup() %>%
    mutate(colour = ifelse(grepl("fov|zoom|full_image",.$name),
                           "transparent", "black")) %>%
    list(.) %>%
    set_names(., sample_id[1])
  
  DATA@tools <<- append(DATA@tools, annot_coord )
  
  img_coord <- temp %>%
    list_modify("fov" = NULL, "full_image" = NULL, "zoom" = NULL) %>%
    compact() %>%
    imap(., ~mutate(img_coord, !!.y := sp::point.in.polygon(
         point.x = img_coord$imagecol*scale.factor,
         point.y = img_coord$imagerow*scale.factor,
         pol.x = .x$x,
         pol.y = .x$y )) ) %>%
    map(., ~select(.x, last_col())) %>%
    cbind(img_coord, .)
  
  DATA@images[[sample_id]]@coordinates <<- img_coord 
  
  sp_annot <- img_coord %>%
    rownames_to_column(., var = "barcodes") %>%
    mutate(across(7:ncol(.), ~ifelse(. == 0, NA, .)) ) %>%
    pivot_longer(., cols = 7:ncol(.), names_to ="sp_annot", values_to = "count") %>%
    filter(!(is.na(count))) %>%
    group_by(barcodes) %>%
    mutate(dupp = row_number()) %>%
    ungroup() %>%
    filter(., .$dupp == 1) %>%
    select(., barcodes, sp_annot)
  
  #DATA <<- left_join( DATA, sp_annot, by=c(".cell"="barcodes", "sp_annot"="sp_annot")) 
  
  return(list(coord=annot_coord, annot=sp_annot))
}
annot <- list(a2, dd, dd2, names(a2)) %>%
  pmap(., ~get_sp_annot(..1, ..2, ..3, ..4)) %>%
  set_names(., names(a2))

sp <- map(annot, 2) %>% bind_rows()

DATA <-  DATA %>%
  left_join(., sp, by=c(".cell"="barcodes")) %>%
  mutate(sp_annot2 = .$sp_annot) %>%
  mutate(sp_annot = ifelse(grepl("epi", .$sp_annot2), "epi", 
                           ifelse(grepl("SubMuc", .$sp_annot2), "SubMuc", .$sp_annot2 ))) %>%
  select(-sp_annot2)

DATA
```

    # A Seurat-tibble abstraction: 31,541 × 5
    # [90mFeatures=36601 | Cells=31541 | Active assay=RNA | Assays=RNA[0m
       .cell                 orig.ident nCount_RNA nFeature_RNA sp_annot
       <chr>                 <chr>           <dbl>        <int> <chr>   
     1 P001_AAACAAGTATCTCCCA P001             3996         1900 SubMuc  
     2 P001_AAACAGGGTCTATATT P001            27612         4795 epi     
     3 P001_AAACCCGAACGAAATC P001            27933         6549 epi     
     4 P001_AAACCGTTCGTCCAGG P001             4066         1761 SubMuc  
     5 P001_AAACCTCATGAAGTTG P001            19936         4449 epi     
     6 P001_AAACGAGACGGTTGAT P001             2906         1579 SubMuc  
     7 P001_AAACTGCTGGCTCCAA P001             4523         2203 SubMuc  
     8 P001_AAACTTGCAAACGTAT P001            18680         3955 epi     
     9 P001_AAAGGCTCTCGCGCCG P001             2868         1396 SubMuc  
    10 P001_AAAGTAGCATTGCTCA P001             4265         2086 SubMuc  
    # ℹ 31,531 more rows

``` r
# adds a matrix with zeros outside the polygon boundary to DATA@misc$alpha
# when plotting it is added as a 4th dim to the image matrix  
# to specify alpha for each pixel in the image
mask_background.fun <- function(df, im, sampleid, manual = TRUE){
  # get polygon coordinates for individual polygons:
  df <- df %>%
    select( any_of(c("name","path_idx", "x", "y"))) %>% 
    filter(grepl("epi|SubMuc", .$name)) %>% 
    split(~name) 
  
  # Function to create SpatialPolygons from a dataframe
  createSpatialPolygons <- function(df) {
    # Convert dataframe to SpatialPointsDataFrame
    #coordinates(df) <- c("y","x")
    coordinates(df) <- c("x", "y")
    
    # Order the points based on 'path_idx'
    df <- df[order(df$path_idx, decreasing = F), ]
    #df <- df %>% arrange(desc(path_idx))
    
    # Create a SpatialPolygons object
    poly <- SpatialPolygons(list(Polygons(list(Polygon(df)), ID = df$name[1])))
    
    return(poly)
  }
  
  # Apply the function to each element of the list
  polygons_list <- map(df, ~createSpatialPolygons(.x))
  
  # Combine the SpatialPolygons into a single object
  combined_polygon <- do.call("rbind", polygons_list)
  
  # Convert it to a polygon (SpatVector)
  polygon <- vect(combined_polygon)
  
  # get image matrix and convert to rgb before rasterizing
  img <- matrix(
      rgb(im[,,1],im[,,2],im[,,3], im[4,,]* 1), nrow=dim(im)[1])
  raster <- rast(img)
  
  rasterized <- rasterize(polygon, raster)
  rasterized <- flip(rasterized, direction="vertical")
  
  # Plot the polygon
  # plot(polygon)

  # Crop the raster using the polygon
  cropped_raster <- mask(raster, rasterized, updatevalue=0)
  
  # Plot the cropped raster to verify
  # plot(cropped_raster)

  cropped_raster <- ifel(cropped_raster > 0, 1, cropped_raster)
  
  image_array <- as.array(cropped_raster)[,,1] # keep first dim only 
  #NA_idx <- which(is.na(image_array[,,1]))
  
  # minpulate image matrix
  #im <- replace(im, list = which(image_array != 0)), values = 1)
  DATA@misc$alpha <<- append(DATA@misc$alpha, set_names(list(image_array), sampleid) )
}

walk(sample_id, ~mask_background.fun(DATA@tools[[.x]], DATA@images[[.x]]@image, .x) )
```

``` r
df_ <- DATA@images[["P031"]]@coordinates %>%
    mutate(scale_fact = DATA@images[["P031"]]@scale.factors$hires) %>%
    mutate(imagecol = .$imagecol * .$scale_fact) %>%
    mutate(imagerow = .$imagerow * .$scale_fact) %>%
  rename( y = "imagerow", x = "imagecol") %>%
  rename(flightlineID="tissue") %>% select(-row, -col, -scale_fact)

ggplot(df_, aes(x=x, y=y)) +
  #geom_point(size=3) + 
  theme_classic() +
    ggstar::geom_star(
      starshape = "hexagon",
      fill = "red",
      size = 2,
      #stroke = 0,
      alpha = 1
    )
  geom_encircle()
  stat_contour_fill(aes(z = tissue))
  geom_density_2d(data = df, aes(x=x, y=y), contour_var = "count", contour = T)
  
m = list(x=df$x, y=df$y, z=matrix(1:12, 4, 3))
  
  resx <- ( m$x[length(m$x)] - m$x[1] ) / (length(m$x)-1)
  resy <- ( m$y[length(m$y)] - m$y[1] ) / (length(m$y)-1)
    xmn <- min(m$x) - 0.5 * resx
    xmx <- max(m$x) + 0.5 * resx
    ymn <- min(m$y) - 0.5 * resy
    ymx <- max(m$y) + 0.5 * resy
  z <- t(m$z)
  z <- z[nrow(z):1, ]

las_df <- df_
  
las_df1 <- las_df[which(las_df$flightlineID == 1), ]
#las_df2 <- las_df[which(las_df$flightlineID == 2L), ]
las_vect1 <- vect(las_df1, geom = c('x', 'y'), crs = 'EPSG:32755')
#las_vect2 <- vect(las_df2, geom = c('X', 'Y'), crs = 'EPSG:32755')
las_rast <- rast(xmin=0, nrow = length(unique(las_df$x)), ncol = length(unique(las_df$y)), crs='EPSG:32755')
set.ext(las_rast, c(min(las_df$x), max(las_df$x), min(las_df$y), max(las_df$y)))
pts1_rast <- rasterize(las_vect1, las_rast, fun = length)
pts2_rast <- rasterize(las_vect2, las_rast, fun = length)
pts1_pts2_rast <- c(pts1_rast, pts2_rast)
names(pts1_pts2_rast) <- c('lyr.1', 'lyr.2') # have to attend to this as both lyr.1 after `c(`
plot(pts1_pts2_rast$lyr.1, col = 'red')

plot(pts1_rast)
bo <- boundaries(pts1_rast, inner=T, classes=T, directions=8, falseval=0,)
bo <- patches(bo, directions=8, zeroAsNA=T, allowGaps=TRUE)
plot(bo)

img <- rast("/Users/vilkal/work/Brolidens_work/Projects/Spatial_Microbiota/data/Spatial_data/P008/spatial/detected_tissue_image.jpg")
plot(img)
bo <- plotRGB(img, 3, 1, 1)
bo <- rast(bo)
bo <- boundaries(bo, inner=T, classes=T, directions=8, falseval=0,)
plot(bo)

contour(pts1_rast, filled=TRUE, nlevels=1)
raster <- contour(img, filled=TRUE, nlevels=1)
raster <- ifel(raster > 200, NA, raster)
plot(raster)

r <- contour(img, filled=TRUE, nlevels=3)
contour(img, add=TRUE)

plot(img)
lines(v)
v <- as.contour(img) %>% st_as_sf()
plot(v)
r <- rast(v)

a <- aggregate(v)

a <- as.polygons(v)
bo <- boundaries(a, inner=F, classes=T)

plot(img)
plot(r)
lines(v)

template <- disagg(rast(raster), 10)
rr <- resample(raster, template)
rr <- floor(rr/100) * 100
v <- as.polygons(rr)
plot(v, 1, col=terrain.colors(7))

cropped_raster <- mask(img, a, updatevalue=0)
```

``` r
https://stackoverflow.com/questions/74153263/calculate-distance-of-points-to-polygon-boundary-using-terra-package-in-r
# can identify the boundary of the tissue and use it to crop
# it gives a more pixelated border than when using a hand drawn polygon as in the function
# above "mask_background.fun()"
find_boundary.fun <- function(im){
  img <- matrix(
      rgb(im[,,1],im[,,2],im[,,3], im[4,,]* 1), nrow=dim(im)[1])
  img <- rast(img)
  #plot(img)
  
  i <- im[,,1]+im[,,2]+im[,,3]
  #i <- replace(im[,,1], list = which(im[,,1] > .7), values = 0)
  
  #raster <- rast(im[,,1])
  raster <- rast(i)
  # plot(raster)
  raster <- ifel(raster > 2.9, NA, raster)
  contour(raster, filled=TRUE, nlevels=5)
  r <- sieve(raster, directions=8, threshold=500)
  # plot(r)
  
  rc2 <- classify(r, matrix(c(0, NA,
                              2, 1,
                              3, 1), ncol=2, byrow=TRUE), 
                  include.lowest=TRUE, brackets=F)
  plot(rc2)
  
  # Extract only the boundary
  bo <- boundaries(rc2, inner=F, classes=T)
  # plot(bo)
  
  # Convert SpatRaster to SpatialPolygonsDataFrame
  boundary_polygon <- as.polygons(bo) 
  # plot(boundary_polygon)
  
  # Crop the raster using the polygon
  cropped_raster <- mask(img, boundary_polygon, updatevalue=0)
  plot(cropped_raster)
  
  return(boundary_polygon)
}

crop_with_poly.fun <-(im, polygon) {
  # get image matrix and convert to rgb before rasterizing
  img <- matrix(
      rgb(im[,,1],im[,,2],im[,,3], im[4,,]* 1), nrow=dim(im)[1])
  raster <- rast(img)
  
  rasterized <- rasterize(polygon, raster)
  rasterized <- flip(rasterized, direction="vertical")
  
  # Plot the polygon
  # plot(polygon)

  # Crop the raster using the polygon
  cropped_raster <- mask(raster, rasterized, updatevalue=0)
  
  # Plot the cropped raster to verify
  # plot(cropped_raster)

  cropped_raster <- ifel(cropped_raster > 0, 1, cropped_raster)
  
  image_array <- as.array(cropped_raster)[,,1] # keep first dim only 
  #NA_idx <- which(is.na(image_array[,,1]))
  
  # minpulate image matrix
  #im <- replace(im, list = which(image_array != 0)), values = 1)
  DATA@misc$alpha <<- append(DATA@misc$alpha, set_names(list(image_array), sampleid) )
}

#boundary_polygon <- find_boundary.fun(DATA@images[["P057"]]@image)
boundary_polygon <- map(sample_id, ~find_boundary.fun(DATA@images[[.x]]@image))
walk(sample_id, ~crop_with_poly.fun(DATA@images[[.x]]@image, .x) )

im <- DATA@images[["P081"]]@image
```

``` r
# dev.new(height=12.5, width=12.5, noRStudioGD = TRUE)
plot_spatial.fun(DATA, 
      sampleid = sample_id,
      sp_annot = T,
      lab = F,
      save_space = T,
      alpha = 0, # hides the spots
      geneid = "nFeature_RNA",
      zoom = NULL,
      ncol = 5,
      img_alpha = 1) + theme(legend.position = "none") # theme_nothing() #
```

<img src="../Figures/00/00a_plot_img.png" data-fig-align="center" />

### Add meta data

``` r
meta <- meta %>%
  select(orig.ident="ID", groups=Luminal_gr_v3, Tissue_gr_v3) %>%
  filter(orig.ident %in% sample_id)

DATA <-  DATA %>%
  left_join(., meta, by="orig.ident") %>%
  select(groups, Tissue_gr_v3, sp_annot, everything())
```

``` r
# dev.new(height=4.5, width=5, noRStudioGD = TRUE)
(p <- plot_spatial.fun(DATA, 
      sampleid = c("P020", "P067", "P044", "P026"), #sample_id,
      sp_annot = T, 
      lab = T,
      save_space = F,
      alpha = 1,
      geneid = "nFeature_RNA",
      zoom = "zoom",
      ncol = 2,
      annot_line = .2,
      point_size = .3,
      img_alpha = 1) +  theme(legend.position = "none")) 

ggsave("./fig_1.pdf", p)

# dev.new(height=12.5, width=12.5, noRStudioGD = TRUE)
p <- DATA %>%
  plot_st_meta.fun(.,  
          assay="RNA",
          feat = "sp_annot",
          zoom = "full_image",
          ncol = 5,
          annot_line = .1,
          annot_col = "black",
          img_alpha = 0,
          point_size = 0.7
        )
```

### Cluster resolutions on tissue

``` r
# dev.new(width=6.6929133858, height=14, noRStudioGD = TRUE)
plots <- DATA %>%
  arrange(groups) %>%
  mutate(id = orig.ident, gr = groups) %>%
  nest(., data = -c("gr", "id")) %>%
  mutate( "sp_annot" = pmap(., 
    ~plot_spatial.fun(..3, sampleid=..1, geneid="sp_annot", img_alpha = .5, colors = c("#FBB4AE", "#B3CDE3", "red"),
                      point_size = 0.8, annot_line = .1, zoom="zoom", save_space = F) + theme(legend.position="none"))) %>%
  mutate( "nFeature_RNA" = pmap(., 
    ~plot_spatial.fun(..3, sampleid=..1, geneid="nFeature_RNA",
                      point_size = 0.8, annot_line = .1, zoom="zoom", save_space = F) + theme(legend.position="none")))

ncols <- table(plots$gr)
combined <- split(plots, ~gr) %>%
  map(., ~plot_grid( plotlist = c(.x$sp_annot, .x$nFeature_RNA), nrow = 2, byrow = T ))
#ggsave("./fig_L1&L3.pdf", combined, width = 28, height = 5)
#ggsave("./fig_.pdf", combined, width = 5, height = 5)
```

<img src="../Figures/00/00b_L1.png" data-fig-align="center" />

<img src="../Figures/00/00c_L2.png" data-fig-align="center" />

<img src="../Figures/00/00d_L3.png" data-fig-align="center" />

<img src="../Figures/00/00e_L4.png" data-fig-align="center" />

## Plot spots to be removed

``` r
# dev.new(width=2.5*5, height=2.5*5, noRStudioGD = TRUE)

DATA <-  DATA %>%
  mutate(filt = case_when(is.na(.$sp_annot)  ~ 'filt',
                          TRUE ~ 'keep')) 
DATA %>%
  plot_st_meta.fun(., 
          lvls = sample_id,
          assay="RNA",
          feat = "filt",
          zoom = "zoom",
          ncol = 5,
          annot_line = .1,
          img_alpha = 0,
          point_size = .8
        )
```

<img src="../Figures/00/00f_plot_spots_to_remove.png"
data-fig-align="center" />

### Remove NA spots

``` r
dim(DATA)
```

    [1] 36601 31541

``` r
# Filter spots outside manual annotation
DATA <- DATA[, !(is.na(DATA$sp_annot))]
DATA$filt <- NULL

dim(DATA)
```

    [1] 36601 30794

``` r
DATA
```

    # A Seurat-tibble abstraction: 30,794 × 7
    # [90mFeatures=36601 | Cells=30794 | Active assay=RNA | Assays=RNA[0m
       .cell         groups Tissue_gr_v3 sp_annot orig.ident nCount_RNA nFeature_RNA
       <chr>         <chr>  <chr>        <chr>    <chr>           <dbl>        <int>
     1 P001_AAACAAG… L3     T3           SubMuc   P001             3996         1900
     2 P001_AAACAGG… L3     T3           epi      P001            27612         4795
     3 P001_AAACCCG… L3     T3           epi      P001            27933         6549
     4 P001_AAACCGT… L3     T3           SubMuc   P001             4066         1761
     5 P001_AAACCTC… L3     T3           epi      P001            19936         4449
     6 P001_AAACGAG… L3     T3           SubMuc   P001             2906         1579
     7 P001_AAACTGC… L3     T3           SubMuc   P001             4523         2203
     8 P001_AAACTTG… L3     T3           epi      P001            18680         3955
     9 P001_AAAGGCT… L3     T3           SubMuc   P001             2868         1396
    10 P001_AAAGTAG… L3     T3           SubMuc   P001             4265         2086
    # ℹ 30,784 more rows

``` r
# dev.new(width=12.5*5, height=12.5*5, noRStudioGD = TRUE)
DATA %>%
  plot_st_meta.fun(.,  
          assay="RNA",
          feat = "sp_annot",
          zoom = "zoom",
          ncol = 5,
          annot_line = .1,
          annot_col = "black",
          img_alpha = 0,
          point_size = 0.8
        )
```

<img src="../Figures/00/00g_final-sp_annot.png"
data-fig-align="center" />

## Save seurat object

``` r
##################################
# SAVE INTERMEDIATE SEURAT OJECT #
##################################
saveRDS(DATA, paste0(result_dir,"seuratObj_merged.RDS"))
# DATA <- readRDS(paste0(result_dir,"seuratObj_merged.RDS")) # seuratObj_merged_seurat_v5.RDS
```

### Session info

``` r
sessionInfo()
```

    R version 4.3.3 (2024-02-29)
    Platform: x86_64-apple-darwin20 (64-bit)
    Running under: macOS Sonoma 14.4.1

    Matrix products: default
    BLAS:   /Library/Frameworks/R.framework/Versions/4.3-x86_64/Resources/lib/libRblas.0.dylib 
    LAPACK: /Library/Frameworks/R.framework/Versions/4.3-x86_64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0

    locale:
    [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

    time zone: Europe/Stockholm
    tzcode source: internal

    attached base packages:
    [1] stats     graphics  grDevices utils     datasets  methods   base     

    other attached packages:
     [1] RColorBrewer_1.1-3 cowplot_1.1.3      terra_1.7-71       xml2_1.3.6        
     [5] hdf5r_1.3.9        Seurat_4.3.0       tidyseurat_0.8.0   SeuratObject_5.0.1
     [9] sp_2.1-3           ttservice_0.4.0    lubridate_1.9.3    forcats_1.0.0     
    [13] stringr_1.5.1      dplyr_1.1.4        purrr_1.0.2        readr_2.1.5       
    [17] tidyr_1.3.1        tibble_3.2.1       ggplot2_3.5.1      tidyverse_2.0.0   

    loaded via a namespace (and not attached):
      [1] rstudioapi_0.16.0      jsonlite_1.8.8         magrittr_2.0.3        
      [4] spatstat.utils_3.0-4   farver_2.1.1           rmarkdown_2.26        
      [7] fs_1.6.4               vctrs_0.6.5            ROCR_1.0-11           
     [10] spatstat.explore_3.2-6 htmltools_0.5.8.1      sctransform_0.4.1     
     [13] parallelly_1.37.0      KernSmooth_2.23-22     htmlwidgets_1.6.4     
     [16] ica_1.0-3              plyr_1.8.9             plotly_4.10.4         
     [19] zoo_1.8-12             igraph_2.0.2           mime_0.12             
     [22] lifecycle_1.0.4        pkgconfig_2.0.3        Matrix_1.6-5          
     [25] R6_2.5.1               fastmap_1.1.1          fitdistrplus_1.1-11   
     [28] future_1.33.1          shiny_1.8.0            digest_0.6.35         
     [31] colorspace_2.1-0       patchwork_1.2.0        tensor_1.5            
     [34] irlba_2.3.5.1          labeling_0.4.3         progressr_0.14.0      
     [37] fansi_1.0.6            spatstat.sparse_3.0-3  timechange_0.3.0      
     [40] httr_1.4.7             polyclip_1.10-6        abind_1.4-5           
     [43] compiler_4.3.3         bit64_4.0.5            withr_3.0.0           
     [46] MASS_7.3-60.0.1        tools_4.3.3            lmtest_0.9-40         
     [49] httpuv_1.6.14          future.apply_1.11.1    goftest_1.2-3         
     [52] glue_1.7.0             nlme_3.1-164           promises_1.2.1        
     [55] grid_4.3.3             Rtsne_0.17             cluster_2.1.6         
     [58] reshape2_1.4.4         generics_0.1.3         gtable_0.3.5          
     [61] spatstat.data_3.0-4    tzdb_0.4.0             data.table_1.15.0     
     [64] hms_1.1.3              utf8_1.2.4             spatstat.geom_3.2-8   
     [67] RcppAnnoy_0.0.22       ggrepel_0.9.5          RANN_2.6.1            
     [70] pillar_1.9.0           spam_2.10-0            vroom_1.6.5           
     [73] later_1.3.2            splines_4.3.3          lattice_0.22-6        
     [76] survival_3.5-8         bit_4.0.5              deldir_2.0-2          
     [79] tidyselect_1.2.1       miniUI_0.1.1.1         pbapply_1.7-2         
     [82] knitr_1.46             gridExtra_2.3          scattermore_1.2       
     [85] xfun_0.43              matrixStats_1.2.0      stringi_1.8.3         
     [88] lazyeval_0.2.2         yaml_2.3.8             evaluate_0.23         
     [91] codetools_0.2-19       cli_3.6.2              uwot_0.1.16           
     [94] xtable_1.8-4           reticulate_1.35.0      munsell_0.5.1         
     [97] Rcpp_1.0.12            globals_0.16.2         spatstat.random_3.2-2 
    [100] png_0.1-8              parallel_4.3.3         ellipsis_0.3.2        
    [103] dotCall64_1.1-1        listenv_0.9.1          viridisLite_0.4.2     
    [106] scales_1.3.0           ggridges_0.5.6         crayon_1.5.2          
    [109] leiden_0.4.3.1         rlang_1.1.3           
