# GAIA-BD GFI Forest Inventory Outlier Filtering
#
# This script identifies multivariate outliers in GFI plot-level forest inventory data
# using bivariate bagplot analysis of basal area and tree density.
# Diagnostic bagplots with whiskers are generated to visualize the distribution
# of observations and detected outliers.
#
# Plot-level outlier PlotIDs are then used to remove corresponding records
# from the tree-level forest inventory data.
#
# The workflow is shown using the North America dataset as an example.
# The same procedure was applied separately by continent using region-specific
# bagplot factors.

rm(list = ls())

library(data.table)
library(aplpack)
library(car)

# ------------------------------------------------------------
# Helper functions for drawing bagplot whiskers
# ------------------------------------------------------------

cut.z.pg <- function(zx, zy, p1x, p1y, p2x, p2y) {
  
  a2 <- (p2y - p1y) / (p2x - p1x)
  a1 <- zy / zx
  
  sx <- (p1y - a2 * p1x) / (a1 - a2)
  sy <- a1 * sx
  sxy <- cbind(sx, sy)
  
  h <- any(is.nan(sxy)) || any(is.na(sxy)) || any(Inf == abs(sxy))
  
  if (h) {
    if (!exists("verbose")) verbose <- FALSE
    
    h <- 0 == zx
    sx <- ifelse(h, zx, sx)
    sy <- ifelse(h, p1y - a2 * p1x, sy)
    
    a1 <- ifelse(abs(a1) == Inf, sign(a1) * 123456789 * 1E10, a1)
    a2 <- ifelse(abs(a2) == Inf, sign(a2) * 123456789 * 1E10, a2)
    
    h <- 0 == (a1 - a2) & sign(zx) == sign(p1x)
    sx <- ifelse(h, p1x, sx)
    sy <- ifelse(h, p1y, sy)
    
    h <- 0 == (a1 - a2) & sign(zx) != sign(p1x)
    sx <- ifelse(h, p2x, sx)
    sy <- ifelse(h, p2y, sy)
    
    h <- p1x == p2x & zx != p1x & p1x != 0
    sx <- ifelse(h, p1x, sx)
    sy <- ifelse(h, zy * p1x / zx, sy)
    
    h <- p1x == p2x & zx != p1x & p1x == 0
    sx <- ifelse(h, p1x, sx)
    sy <- ifelse(h, 0, sy)
    
    h <- p1x == p2x & zx == p1x & p1x != 0
    sx <- ifelse(h, zx, sx)
    sy <- ifelse(h, zy, sy)
    
    h <- p1x == p2x & zx == p1x & p1x == 0 & sign(zy) == sign(p1y)
    sx <- ifelse(h, p1x, sx)
    sy <- ifelse(h, p1y, sy)
    
    h <- p1x == p2x & zx == p1x & p1x == 0 & sign(zy) != sign(p1y)
    sx <- ifelse(h, p1x, sx)
    sy <- ifelse(h, p2y, sy)
    
    h <- zx == p1x & zy == p1y
    sx <- ifelse(h, p1x, sx)
    sy <- ifelse(h, p1y, sy)
    
    h <- zx == p2x & zy == p2y
    sx <- ifelse(h, p2x, sx)
    sy <- ifelse(h, p2y, sy)
    
    h <- zx == 0 & zy == 0
    sx <- ifelse(h, 0, sx)
    sy <- ifelse(h, 0, sy)
    
    sxy <- cbind(sx, sy)
  }
  
  return(sxy)
}

find.cut.z.pg <- function(z, pg, center = c(0, 0)) {
  
  win <- function(dx, dy) atan2(y = dy, x = dx)
  
  if (!is.matrix(z)) z <- rbind(z)
  if (nrow(pg) == 1) return(matrix(center, nrow(z), 2, TRUE))
  
  n.pg <- nrow(pg)
  
  z <- cbind(z[, 1] - center[1], z[, 2] - center[2])
  pg <- cbind(pg[, 1] - center[1], pg[, 2] - center[2])
  
  apg <- win(pg[, 1], pg[, 2])
  apg[is.nan(apg)] <- 0
  
  a <- order(apg)
  apg <- apg[a]
  pg <- pg[a, ]
  
  az <- win(z[, 1], z[, 2])
  
  segm.no <- apply((outer(apg, az, "<")), 2, sum)
  segm.no <- ifelse(segm.no == 0, n.pg, segm.no)
  next.no <- 1 + (segm.no %% length(apg))
  
  cuts <- cut.z.pg(
    z[, 1], z[, 2],
    pg[segm.no, 1], pg[segm.no, 2],
    pg[next.no, 1], pg[next.no, 2]
  )
  
  cuts <- cbind(cuts[, 1] + center[1], cuts[, 2] + center[2])
  
  return(cuts)
}

draw_whiskers <- function(bag_trans, color = "gold3") {
  
  xy <- bag_trans$xy
  pxy.outer <- bag_trans$pxy.outer
  hull.bag <- bag_trans$hull.bag
  center <- bag_trans$center
  
  if (length(pxy.outer) > 0) {
    
    if ((n <- length(xy[, 1])) < 15) {
      segments(
        xy[, 1], xy[, 2],
        rep(center[1], n), rep(center[2], n),
        col = color
      )
    } else {
      pkt.cut <- find.cut.z.pg(pxy.outer, hull.bag, center = center)
      segments(
        pxy.outer[, 1], pxy.outer[, 2],
        pkt.cut[, 1], pkt.cut[, 2],
        col = color
      )
    }
  }
}

# ------------------------------------------------------------
# Load plot-level forest inventory data
# ------------------------------------------------------------

setwd("F:/GAIA_datapaper/final_code/Outliers")

tree <- fread("Plot_level_FIA.csv")

# Region-specific bagplot factor
# Africa = 4.4; North America = 4.5; South America = 5.0;
# Oceania = 3.2; Europe = 4.0; Asia = 5.5
bagplot_factor <- 4.5

# ------------------------------------------------------------
# Identify multivariate outliers
# ------------------------------------------------------------

# Retain records with non-missing basal area and tree density
id_not_na <- which(!is.na(tree$B) & !is.na(tree$N))

# Create bivariate input using basal area and tree density
a <- with(tree, as.matrix(cbind(B[id_not_na], N[id_not_na])))

if (is.matrix(a)) {
  a <- list(
    x = a[, 1],
    y = a[, 2]
  )
}

# Estimate Box-Cox transformation parameters
p1 <- powerTransform(cbind(y, x) ~ 1, a)

# Identify multivariate outliers using bivariate bagplot analysis
bag_trans <- with(
  a,
  compute.bagplot(
    x^coef(p1)["x"],
    y^coef(p1)["y"],
    approx.limit = 10000,
    factor = bagplot_factor,
    precision = 1,
    dkmethod = 2,
    debug.plots = TRUE
  )
)

# ------------------------------------------------------------
# Plot bagplot diagnostics
# ------------------------------------------------------------

plot(
  bag_trans,
  xlim = range(a$x^coef(p1)["x"]),
  ylim = range(a$y^coef(p1)["y"]),
  xlab = "Transformed basal area",
  ylab = "Transformed tree density",
  cex.lab = 1.5,
  main = "GFI multivariate outlier detection"
)

draw_whiskers(bag_trans, color = "gold3")

outlier_points <- bag_trans$pxy.outlier

if (!is.null(outlier_points) && nrow(outlier_points) > 0) {
  points(
    outlier_points[, 1],
    outlier_points[, 2],
    pch = 19,
    col = "red",
    cex = 0.2
  )
}

# ------------------------------------------------------------
# Extract and remove outlier records
# ------------------------------------------------------------

k <- which(
  bag_trans$xydata[, 1] %in% bag_trans$pxy.outlier[, 1] &
    bag_trans$xydata[, 2] %in% bag_trans$pxy.outlier[, 2]
)

id_outliers <- unlist(
  sapply(
    k,
    function(i) which(tree$B == a$x[i] & tree$N == a$y[i])
  )
)

tree_outliers <- tree[id_outliers, ]
tree_clean <- tree[-id_outliers, ]

# ------------------------------------------------------------
# Export outlier and cleaned datasets
# ------------------------------------------------------------

fwrite(
  tree_outliers,
  "F:/GAIA_datapaper/final_code/Outliers/Outliers_Plot_level_FIA3.csv"
)

fwrite(
  tree_clean,
  "F:/GAIA_datapaper/final_code/Outliers/Cleaned_Outliers_Plot_level_FIA3.csv"
)


# ------------------------------------------------------------
# Apply plot-level outlier removal to tree-level records
# ------------------------------------------------------------

# Load North America tree-level forest inventory data
tree_level <- fread(
  "F:/GAIA_datapaper/final_code/Outliers/Tree_level_FIA.csv"
)

# Remove tree-level records belonging to plot-level outlier PlotIDs
tree_level_clean <- tree_level[!(PlotID %in% tree_outliers$PlotID)]

# Extract tree-level records belonging to plot-level outlier PlotIDs
tree_level_outliers <- tree_level[PlotID %in% tree_outliers$PlotID]

# Export cleaned tree-level data
fwrite(
  tree_level_clean,
  "F:/GAIA_datapaper/final_code/Outliers/Cleaned_Outliers_tree_level_FIA3.csv"
)

# Export removed tree-level outlier records
fwrite(
  tree_level_outliers,
  "F:/GAIA_datapaper/final_code/Outliers/Outliers_treelevel_FIA3.csv"
)
