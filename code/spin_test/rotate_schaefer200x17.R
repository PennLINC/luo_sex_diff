library(freesurfer)
library(freesurferformats)

# Set project root depending on whether code is running on PARCC or locally
if (dir.exists("/ceph/projects/sattertt/pennlinc-parcc/network_replication")) {
  root <- "/ceph/projects/sattertt/pennlinc-parcc/network_replication"
} else if (dir.exists("~/parcc/network_replication")) {
  root <- path.expand("~/parcc/network_replication")
} else {
  stop("Could not find network_replication project.")
}

rotate_dir <- file.path(root, "software", "rotate_parcellation")
out_dir <- file.path(root, "covariate_analyses", "sex_diff", "software", "rotate_parcellation")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# rotate.parcellation.R and perm.sphere.p.R are downloaded from https://github.com/frantisekvasa/rotate_parcellation
source(file.path(root, "software", "rotate.parcellation.R"))
source(file.path(root, "software", "perm.sphere.p.R"))

# Extract parcel centroid coordinates on the fsaverage sphere from Schaefer .annot files
# load freesurfer fsaverage surface
lh <- read.fs.surface(file.path(rotate_dir, "lh.sphere"))
rh <- read.fs.surface(file.path(rotate_dir, "rh.sphere"))

extract_centroid_coords <- function(schaefer_atlas, network_parcellation) {
  lh.schaefer.annot <- read_annotation(file.path(rotate_dir, "annot", sprintf("lh.Schaefer2018_%sParcels_%sNetworks_order.annot", schaefer_atlas, network_parcellation)))
  rm <- which(lh.schaefer.annot$colortable$label == "Background+FreeSurfer_Defined_Medial_Wall")
  lhcoords <- t(sapply(seq_along(lh.schaefer.annot$colortable$label)[-rm], function(x) colMeans(lh$vertices[which(lh.schaefer.annot$label == lh.schaefer.annot$colortable$code[x]), ])))
  
  rh.schaefer.annot <- read_annotation(file.path(rotate_dir, "annot", sprintf("rh.Schaefer2018_%sParcels_%sNetworks_order.annot", schaefer_atlas, network_parcellation)))
  rm <- which(rh.schaefer.annot$colortable$label == "Background+FreeSurfer_Defined_Medial_Wall")
  rhcoords <- t(sapply(seq_along(rh.schaefer.annot$colortable$label)[-rm], function(x) colMeans(rh$vertices[which(rh.schaefer.annot$label == rh.schaefer.annot$colortable$code[x]), ])))
  
  schaefer.fsaverage.coords <- rbind(lhcoords, rhcoords)
  write.csv(schaefer.fsaverage.coords, file.path(out_dir, sprintf("schaefer%sx%s_sphericalcoords.csv", schaefer_atlas, network_parcellation)), row.names = FALSE)
}

extract_centroid_coords(200, 17) # schaefer200x17

# label the coordinates of schaefer200x17.coords parcel centroids on the freesurfer sphere
schaefer200x17.coords <- read.csv(file.path(out_dir, "schaefer200x17_sphericalcoords.csv"))
lh.schaefer.annot <- read_annotation(file.path(rotate_dir, "annot", "lh.Schaefer2018_200Parcels_17Networks_order.annot"))
rh.schaefer.annot <- read_annotation(file.path(rotate_dir, "annot", "rh.Schaefer2018_200Parcels_17Networks_order.annot"))
label <- data.frame(c(lh.schaefer.annot$colortable[2:101, 1], rh.schaefer.annot$colortable[2:101, 1]))
schaefer200x17.coords <- cbind(label, schaefer200x17.coords)
names(schaefer200x17.coords)[1] <- "label"
write.csv(schaefer200x17.coords, file.path(out_dir, "schaefer200x17_sphericalcoords.csv"), row.names = FALSE)

# rotate the Schaefer200x17 parcellation 10,000 times on the FreeSurfer sphere
set.seed(10)
perm.id.full <- rotate.parcellation(coord.l = as.matrix(schaefer200x17.coords[1:100, 2:4]), coord.r = as.matrix(schaefer200x17.coords[101:200, 2:4]), nrot = 10000)
saveRDS(perm.id.full, file.path(out_dir, "schaefer200x17.coords_sphericalrotations_N10k_seed10.rds"))
