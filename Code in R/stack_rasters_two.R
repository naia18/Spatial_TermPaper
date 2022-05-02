# ========================================================================
# stack_rasters_two.R    -   Naia Ormaza Zulueta   -  Apr 2022
# In this file:
# - I do the same as in stack_rasters.R, but considering only LMICs
# - The raster for the Water Issue is given by the access to water 
# sanitation facilities
# ========================================================================

# Clear the environment
rm(list=ls()) 
# load libraries
packages <- c("terra", "raster", "tidyverse", "rasterVis", "ncdf4", 
              "lattice", "foreign", "rworldmap")
lapply(packages, require, character=TRUE)


# Open raster data
airQuality <- rast("/Users/naiacasina/Documents/IDEA SECOND/Sem 3/ENVS/Codes and Data/Rasters/airQuality.tif")
water <- rast("/Users/naiacasina/Documents/IDEA SECOND/Sem 3/ENVS/Codes and Data/Water/Percent_sanitationAccess/IHME_LMIC_WASH_2000_2017_S_IMP_PERCENT_MEAN_2017_Y2020M06D02.TIF")
biodiv <- rast("/Users/naiacasina/Documents/IDEA SECOND/Sem 3/ENVS/Codes and Data/Rasters/biodiv.tif")
pop <- rast("/Users/naiacasina/Documents/IDEA SECOND/Sem 3/ENVS/Codes and Data/Rasters/pop.tif")
heat <- rast("/Users/naiacasina/Documents/IDEA SECOND/Sem 3/ENVS/Codes and Data/Rasters/heat.tif")
flood <- rast("/Users/naiacasina/Documents/IDEA SECOND/Sem 3/ENVS/Codes and Data/Rasters/flood.tif")

newres <- 0.1
water <- aggregate(pop,newres/res(pop),fun=sum)
water <- resample(water, pop)

s <- c(airQuality, water, biodiv, heat, flood, pop)

df <- as.data.frame(s, xy=TRUE, na.rm=FALSE)
df$floodMapGL_rp10y[is.na(df$floodMapGL_rp10y)] <- 0
df <- na.omit(df)

# Absolute threshold
df$air <- 0
df$water <- 0
df$biodiv <- 0
df$safecl <- 0
for (i in 1:dim(df)[1]) {
  if(df$`V5GL02.HybridPM25c_0p10.Global.202001-202012`[i]>35){df$air[i] <- 1}
  if(df$IHME_LMIC_WASH_2000_2017_S_IMP_PERCENT_MEAN_2017_Y2020M06D02[i]<70){df$water[i] <- 1}
  if(df$lbii[i]<0.80){df$biodiv[i] <- 1}
  h <- 0
  f <- 0
  if(df$heat[i]>55){h <- 1}
  if(df$floodMapGL_rp10y[i]>0){f <- 1}
  df$safecl[i] <- (h+f)
}

# Save Data Frame
save(df,file="/Users/naiacasina/Documents/IDEA SECOND/Sem 3/ENVS/Codes and Data/R/Data/stack_layers_two.Rdata")
df_two <- df
save(df_two, file= "/Users/naiacasina/Documents/IDEA SECOND/Sem 3/ENVS/Codes and Data/R/Data/stack_layers_two.Rdata")
# Human Rights Violations
df$HRV <- rowSums(df[,c(9:12)])
df$logpop <- log(df$gpw_v4_population_count_rev11_2020_2pt5_min)
#df$logpop[which(df$logpop==-Inf)] <- 0
#df$logpop[which(df$logpop<0)] <- 0

# --- Dfs into rasters -------
# Separate
df_rasterx <- df[,c('x','y','HRV')]
df_rastery <- df[,c('x','y','gpw_v4_population_count_rev11_2020_2pt5_min')]
# Ens
df_rast <- df[,c('x','y','HRV','gpw_v4_population_count_rev11_2020_2pt5_min')]

# ------ Create rasters -------
# Separate
rasterx <- rast(df_rasterx, type="xyz")
rastery <- rast(df_rastery, type="xyz")
# Ens
raster_HRV <- rast(df_rast, type="xyz")
# write rasters
writeRaster(raster_HRV, filename="/Users/naiacasina/Documents/IDEA SECOND/Sem 3/ENVS/Results/two_HRV_2.tif", overwrite=TRUE)
writeRaster(rasterx, filename="/Users/naiacasina/Documents/IDEA SECOND/Sem 3/ENVS/Results/two_HRV_Layer1.tif", overwrite=TRUE)
writeRaster(rastery, filename="/Users/naiacasina/Documents/IDEA SECOND/Sem 3/ENVS/Results/two_Population_Layer2.tif", overwrite=TRUE)

rasterx <- rast("/Users/naiacasina/Documents/IDEA SECOND/Sem 3/ENVS/Results/two_HRV_Layer1.tif")
rastery <- rast("/Users/naiacasina/Documents/IDEA SECOND/Sem 3/ENVS/Results/two_Population_Layer2.tif")

# --------------------- BIVARIATE PLOTTING ---------------------

library(data.table)
library(tidyverse)
library(raster)
library(classInt)
library(patchwork)

# Clip to Globe
clipExt <- extent(-180, 180, -55, 70)
# Clip to Europe
clipExt_eur <- extent(-15, 40, 30, 68)
# Clip to North America and Central America
clipExt_NA <- extent(-160,-40, 5, 65)
# Clip to SE Asia
clipExt_SEA <- extent(45,140,-10,45)


# Define the number of breaks
nBreaks <- 4

# Create the colour matrix
col.matrixQ <- colmat(nbreaks = nBreaks, breakstyle = "quantile",
                      xlab = "HRV", ylab = "Population", 
                      bottomright = "#F7900A", upperright = "#993A65",
                      bottomleft = "#44B360", upperleft = "#3A88B5",
                      saveLeg = FALSE, plotLeg = TRUE)

# create the bivariate raster
rasterx_br <- brick(rasterx)
rastery_br <- brick (rastery)

bivmapQ <- bivariate.map(rasterx = rasterx_br, rastery =  rastery_br,
                         export.colour.matrix = FALSE,
                         colourmatrix = col.matrixQ)

# Convert to dataframe for plotting with ggplot
bivMapDFQ <- setDT(as.data.frame(bivmapQ, xy = TRUE))
colnames(bivMapDFQ)[3] <- "BivValue"
bivMapDFQ <- melt(bivMapDFQ, id.vars = c("x", "y"),
                  measure.vars = "BivValue",
                  value.name = "bivVal",
                  variable.name = "Variable")

# Make the map using ggplot
map_q <- ggplot(bivMapDFQ, aes(x = x, y = y)) +
  geom_raster(aes(fill = bivVal)) +
  scale_y_continuous(breaks = seq(-55, 75, by = 20), 
                     labels = paste0(seq(-55, 75, 20), "°")) +
  scale_x_continuous(breaks = seq(-180,180,30), 
                     labels = paste0(seq(-180,180,30), "°")) +
  scale_fill_gradientn(colours = col.matrixQ, na.value = "transparent") + 
  theme_bw() +
  theme(text = element_text(size = 10, colour = "black")) +
  borders(colour = "black", size = 0.2) +
  coord_quickmap(expand = FALSE, xlim = clipExt[1:2], ylim = clipExt[3:4]) +
  theme(legend.position = "none",
        plot.background = element_blank(),
        strip.text = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(angle = 90, hjust = 0.5),
        axis.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 12, colour = "black")) +
  labs(x = "Longitude", y = "Latitude")


map_q_eur <- ggplot(bivMapDFQ, aes(x = x, y = y)) +
  geom_raster(aes(fill = bivVal)) +
  scale_y_continuous(breaks = seq(30,68, by = 10), 
                     labels = paste0(seq(30,68, 10), "°")) +
  scale_x_continuous(breaks = seq(-15, 40, 5), 
                     labels = paste0(seq(-15, 40, 5), "°")) +
  scale_fill_gradientn(colours = col.matrixQ, na.value = "transparent") + 
  theme_bw() +
  theme(text = element_text(size = 10, colour = "black")) +
  borders(colour = "black", size = 0.2) +
  coord_quickmap(expand = FALSE, xlim = clipExt_eur[1:2], ylim = clipExt_eur[3:4]) +
  theme(legend.position = "none",
        plot.background = element_blank(),
        strip.text = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(angle = 90, hjust = 0.5),
        axis.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 12, colour = "black")) +
  labs(x = "Longitude", y = "Latitude")

map_q_na <- ggplot(bivMapDFQ, aes(x = x, y = y)) +
  geom_raster(aes(fill = bivVal)) +
  scale_y_continuous(breaks = seq(5, 65, by = 20), 
                     labels = paste0(seq(5, 65, 20), "°")) +
  scale_x_continuous(breaks = seq(-160,-40,20), 
                     labels = paste0(seq(-160,-40,20), "°")) +
  scale_fill_gradientn(colours = col.matrixQ, na.value = "transparent") + 
  theme_bw() +
  theme(text = element_text(size = 10, colour = "black")) +
  borders(colour = "black", size = 0.2) +
  coord_quickmap(expand = FALSE, xlim = clipExt_NA[1:2], ylim = clipExt_NA[3:4]) +
  theme(legend.position = "none",
        plot.background = element_blank(),
        strip.text = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(angle = 90, hjust = 0.5),
        axis.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 12, colour = "black")) +
  labs(x = "Longitude", y = "Latitude")


map_q_sea <- ggplot(bivMapDFQ, aes(x = x, y = y)) +
  geom_raster(aes(fill = bivVal)) +
  scale_y_continuous(breaks = seq(-10, 45, by = 20), 
                     labels = paste0(seq(-10, 45, 20), "°")) +
  scale_x_continuous(breaks = seq(45,140,20), 
                     labels = paste0(seq(45,140,20), "°")) +
  scale_fill_gradientn(colours = col.matrixQ, na.value = "transparent") + 
  theme_bw() +
  theme(text = element_text(size = 10, colour = "black")) +
  borders(colour = "black", size = 0.2) +
  coord_quickmap(expand = FALSE, xlim = clipExt_SEA[1:2], ylim = clipExt_SEA[3:4]) +
  theme(legend.position = "none",
        plot.background = element_blank(),
        strip.text = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(angle = 90, hjust = 0.5),
        axis.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 12, colour = "black")) +
  labs(x = "Longitude", y = "Latitude")

# Using different breaks algorithm -----------------------------------
# Create the colour matrix
col.matrixF <- colmat(nbreaks = nBreaks, breakstyle = "fisher",
                      xlab = "Temperature", ylab = "Precipiation", 
                      bottomright = "#F7900A", upperright = "#993A65",
                      bottomleft = "#44B360", upperleft = "#3A88B5",
                      saveLeg = FALSE, plotLeg = TRUE)

# create the bivariate raster
bivmapF <- bivariate.map(rasterx = rasterx_br, rastery = rastery_br,
                         export.colour.matrix = FALSE,
                         colourmatrix = col.matrixF)
bivMapDFF <- setDT(as.data.frame(bivmapF, xy = TRUE))
colnames(bivMapDFF)[3] <- "BivValue"
bivMapDFF <- melt(bivMapDFF, id.vars = c("x", "y"),
                  measure.vars = "BivValue",
                  value.name = "bivVal",
                  variable.name = "Variable")

# Make the map using ggplot
map_F <- ggplot(bivMapDFF, aes(x = x, y = y)) +
  geom_raster(aes(fill = bivVal)) +
  scale_y_continuous(breaks = seq(-20, 70, by = 20), 
                     labels = paste0(seq(-20, 70, 20), "°")) +
  scale_x_continuous(breaks = seq(50,175,25), 
                     labels = paste0(seq(50,175,25), "°")) +
  scale_fill_gradientn(colours = col.matrixF, na.value = "transparent") + 
  theme_bw() +
  theme(text = element_text(size = 10, colour = "black")) +
  borders(colour = "black", size = 0.2) +
  coord_quickmap(expand = FALSE, xlim = clipExt[1:2], ylim = clipExt[3:4]) +
  theme(legend.position = "left",
        plot.background = element_blank(),
        strip.text = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(angle = 90, hjust = 0.5),
        axis.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 12, colour = "black")) +
  labs(x = "Longitude", y = "Latitude")

fig_global <- {map_q + ggtitle("Human Rights Violations and Population with Quantile breaks")} + 
  inset_element(BivLegend + theme(plot.background = element_rect(fill = "white",
                                                                 colour = NA)), 
                left = 0.1, bottom = 0.2, right = 0.3, top = 0.5,
                align_to = "full") +
  plot_annotation(caption = "Capt")
fig_global

fig_eur <- {map_q_eur + ggtitle("HRV Europe")} + 
  inset_element(BivLegend + theme(plot.background = element_rect(fill = "white",
                                                                 colour = NA)), 
                left = 0.15, bottom = 0.77, right = 0.3, top = 0.92,
                align_to = "full") +
  plot_annotation(caption = "Capt")
fig_eur

fig_na <- {map_q_na + ggtitle("HRV North and Central America")} + 
  inset_element(BivLegend + theme(plot.background = element_rect(fill = "white",
                                                                 colour = NA)), 
                left = 0.73, bottom = 0.65, right = 0.98, top = 0.90,
                align_to = "full") +
  plot_annotation(caption = "Capt")
fig_na

fig_sea <- {map_q_sea + ggtitle("HRV South and Southeast Asia")} + 
  inset_element(BivLegend + theme(plot.background = element_rect(fill = "white",
                                                                 colour = NA)), 
                left = 0.1, bottom = 0.2, right = 0.3, top = 0.5,
                align_to = "full") +
  plot_annotation(caption = "Capt")
fig_sea

# Save
ggsave(plot = fig_global,
       filename = "two_BivariatePlot_Global_ggsave.png",
       device = "png", path = "/Users/naiacasina/Documents/IDEA SECOND/Sem 3/ENVS/Results/",
       width = 6, height = 7, units = "in",
       dpi = 320)

ggsave(plot = fig_eur,
       filename = "two_BivariatePlot_Europe.png",
       device = "png", path = "/Users/naiacasina/Documents/IDEA SECOND/Sem 3/ENVS/Results/",
       width = 6, height = 7, units = "in",
       dpi = 320)

ggsave(plot = fig_na,
       filename = "two_BivariatePlot_NorthCentralAmerica.png",
       device = "png", path = "/Users/naiacasina/Documents/IDEA SECOND/Sem 3/ENVS/Results/",
       width = 6, height = 7, units = "in",
       dpi = 320)

ggsave(plot = fig_sea,
       filename = "two_BivariatePlot_SEAsia.png",
       device = "png", path = "/Users/naiacasina/Documents/IDEA SECOND/Sem 3/ENVS/Results/",
       width = 6, height = 7, units = "in",
       dpi = 320)

# ------------------------ BIVARIATE PLOTTING FUNCTIONS ----------------------

# The function that produces the colour matrix
colmat <- function(nbreaks = 3, breakstyle = "quantile",
                   upperleft = "#0096EB", upperright = "#820050", 
                   bottomleft = "#BEBEBE", bottomright = "#FFE60F",
                   xlab = "x label", ylab = "y label", plotLeg = TRUE,
                   saveLeg = FALSE) {
  # TODO - replace any tidyr, dplyr etc. functions with data.table #
  library(tidyverse)
  require(ggplot2)
  require(classInt)
  if (breakstyle == "sd") {
    warning("SD breaks style cannot be used.\nWill not always return the correct number of breaks.\nSee classInt::classIntervals() for details.\nResetting to quantile",
            call. = FALSE, immediate. = FALSE)
    breakstyle <- "quantile"}
  
  my.data <- seq(0, 1, .01)
  my.class <- classInt::classIntervals(my.data,
                                       n = nbreaks,
                                       style = breakstyle,
  )
  my.pal.1 <- classInt::findColours(my.class, c(upperleft, bottomleft))
  my.pal.2 <- classInt::findColours(my.class, c(upperright, bottomright))
  col.matrix <- matrix(nrow = 101, ncol = 101, NA)
  for (i in 1:101) {
    my.col <- c(paste(my.pal.1[i]), paste(my.pal.2[i]))
    col.matrix[102 - i, ] <- classInt::findColours(my.class, my.col)
  }
  ## need to convert this to data.table at some stage.
  col.matrix.plot <- col.matrix %>%
    as.data.frame(.) %>% 
    mutate("Y" = row_number()) %>%
    mutate_at(.tbl = ., .vars = vars(starts_with("V")), .funs = list(as.character)) %>% 
    pivot_longer(data = ., cols = -Y, names_to = "X", values_to = "HEXCode") %>% 
    mutate("X" = as.integer(sub("V", "", .$X))) %>%
    distinct(as.factor(HEXCode), .keep_all = TRUE) %>%
    mutate(Y = rev(.$Y)) %>% 
    dplyr::select(-c(4)) %>%
    mutate("Y" = rep(seq(from = 1, to = nbreaks, by = 1), each = nbreaks),
           "X" = rep(seq(from = 1, to = nbreaks, by = 1), times = nbreaks)) %>%
    mutate("UID" = row_number())
  # Use plotLeg if you want a preview of the legend
  if (plotLeg) {
    p <- ggplot(col.matrix.plot, aes(X, Y, fill = HEXCode)) +
      geom_tile() +
      scale_fill_identity() +
      coord_equal(expand = FALSE) +
      theme_void() +
      theme(aspect.ratio = 1,
            axis.title = element_text(size = 12, colour = "black",hjust = 0.5, 
                                      vjust = 1),
            axis.title.y = element_text(angle = 90, hjust = 0.5)) +
      xlab(bquote(.(xlab) ~  symbol("\256"))) +
      ylab(bquote(.(ylab) ~  symbol("\256")))
    print(p)
    assign(
      x = "BivLegend",
      value = p,
      pos = .GlobalEnv
    )
  }
  # Use saveLeg if you want to save a copy of the legend
  if (saveLeg) {
    ggsave(filename = "bivLegend.pdf", plot = p, device = "pdf",
           path = "./", width = 4, height = 4, units = "in",
           dpi = 300)
  }
  seqs <- seq(0, 100, (100 / nbreaks))
  seqs[1] <- 1
  col.matrix <- col.matrix[c(seqs), c(seqs)]
  attr(col.matrix, "breakstyle") <- breakstyle
  attr(col.matrix, "nbreaks") <- nbreaks
  return(col.matrix)
}



# Function to assign colour-codes to raster data
# As before, by default assign tercile breaks
bivariate.map <- function(rasterx = rasterx_br, rastery = rastery_br, colourmatrix = col.matrix,
                          export.colour.matrix = TRUE,
                          outname = paste0("colMatrix_rasValues", names(rasterx))) {
  # TO DO - replace raster with terra #
  require(raster)
  require(classInt)
  # export.colour.matrix will export a data.frame of rastervalues and RGB codes 
  # to the global environment outname defines the name of the data.frame
  quanx <- getValues(rasterx)
  tempx <- data.frame(quanx, quantile = rep(NA, length(quanx)))
  brks <- with(tempx, classIntervals(quanx,
                                     n = attr(colourmatrix, "nbreaks"),
                                     style = attr(colourmatrix, "breakstyle"))$brks)
  ## Add (very) small amount of noise to all but the first break
  ## https://stackoverflow.com/a/19846365/1710632
  brks[-1] <- brks[-1] + seq_along(brks[-1]) * .Machine$double.eps
  r1 <- within(tempx, quantile <- cut(quanx,
                                      breaks = brks,
                                      labels = 2:length(brks),
                                      include.lowest = TRUE))
  quantr <- data.frame(r1[, 2])
  quany <- getValues(rastery)
  tempy <- data.frame(quany, quantile = rep(NA, length(quany)))
  brksy <- with(tempy, classIntervals(quany,
                                      n = attr(colourmatrix, "nbreaks"),
                                      style = attr(colourmatrix, "breakstyle"))$brks)
  brksy[-1] <- brksy[-1] + seq_along(brksy[-1]) * .Machine$double.eps
  r2 <- within(tempy, quantile <- cut(quany,
                                      breaks = brksy,
                                      labels = 2:length(brksy),
                                      include.lowest = TRUE
  ))
  quantr2 <- data.frame(r2[, 2])
  as.numeric.factor <- function(x) {
    as.numeric(levels(x))[x]
  }
  col.matrix2 <- colourmatrix
  cn <- unique(colourmatrix)
  for (i in 1:length(col.matrix2)) {
    ifelse(is.na(col.matrix2[i]),
           col.matrix2[i] <- 1, col.matrix2[i] <- which(
             col.matrix2[i] == cn
           )[1]
    )
  }
  # Export the colour.matrix to data.frame() in the global env
  # Can then save with write.table() and use in ArcMap/QGIS
  # Need to save the output raster as integer data-type
  if (export.colour.matrix) {
    # create a dataframe of colours corresponding to raster values
    exportCols <- as.data.frame(cbind(
      as.vector(col.matrix2), as.vector(colourmatrix),
      t(col2rgb(as.vector(colourmatrix)))
    ))
    # rename columns of data.frame()
    colnames(exportCols)[1:2] <- c("rasValue", "HEX")
    # Export to the global environment
    assign(
      x = outname,
      value = exportCols,
      pos = .GlobalEnv
    )
  }
  cols <- numeric(length(quantr[, 1]))
  for (i in 1:length(quantr[, 1])) {
    a <- as.numeric.factor(quantr[i, 1])
    b <- as.numeric.factor(quantr2[i, 1])
    cols[i] <- as.numeric(col.matrix2[b, a])
  }
  r <- rasterx
  r[1:length(r)] <- cols
  return(r)
}
