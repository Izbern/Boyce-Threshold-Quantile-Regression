# Read environment variables
library(raster)
library(virtualspecies)
library(rasterVis)
raster_brick <- brick("./EnvVar/princip.tif")
raster_stack <- stack(raster_brick)
layer_names <- c("var1", "var2", "var3")
names(raster_stack) <- layer_names
band1 <- raster_brick[[1]]
band2 <- raster_brick[[2]]
band3 <- raster_brick[[3]]          
histogram(band1, main = "Band 1 Histogram")
histogram(band2, main = "Band 2 Histogram")
histogram(band3, main = "Band 3 Histogram")          
sd_vals <- cellStats(raster_brick, sd)

#generate virtual species
num_species<- 3
output_folder <- "response_curves"
dir.create(output_folder, showWarnings = FALSE)  # Create the folder if it doesn't exist
species_list <- list()
for (i in num_species:num_species) {
  my.parameters <- formatFunctions(
    x = raster_stack,
    rescale = TRUE,
    var1 = c(fun = "dnorm", mean = random_var1_mean, sd = random_var1_sd),
    var2 = c(fun = "linearFun", a = random_var2_mean, b = random_var2_sd),
    var3 = c(fun = "logisticFun", alpha = random_var3_mean, beta = random_var3_sd)
  )
  # Generate the virtual species using the current parameters
  sp <- generateSpFromFun(raster_stack, my.parameters, species.type = "additive",)
  # Store the virtual species in the list
  species_list[[i]] <- sp
  plot_file <- file.path(output_folder, paste("ResponsePlot_Species", i, ".png", sep = ""))
  png(file = plot_file)
  response_plot <- plotResponse(sp)
  dev.off()
}


 # Set the prevalence of virtual species
prevalence_list <- list(0.1,0.25,0.5,0.75)
# Convert to presence-absence map
pa_folder <- "presence_absence_maps"
SuToProba_folder <- "SuToProb_Plot"
dir.create(pa_folder, showWarnings = FALSE)  
dir.create(SuToProba_folder, showWarnings = FALSE)  
rds_main_folder <- "RDS"
dir.create(rds_main_folder, showWarnings = FALSE)  
num_repetitions <- 1
i <- 3
prevalence <- 0.5
# Loop for each virtual species
if (TRUE) {
  # Create a subfolder for this species
  species_pa_folder <- file.path(pa_folder, paste("Species", i))
  species_SuToProba_folder <- file.path(SuToProba_folder, paste("Species", i))
  dir.create(species_pa_folder, showWarnings = FALSE)
  dir.create(species_SuToProba_folder, showWarnings = FALSE)
  rds_folder <- file.path(rds_main_folder, paste("Species", i))
  dir.create(rds_folder, showWarnings = FALSE)
  # Loop for each repetition
  for (j in 1:num_repetitions) {
    # Generate presence-absence map using convertToPA
    pa_map <- convertToPA(species_list[[i]], PA.method = "probability", 
                          prob.method = "logistic", species.prevalence = prevalence)
    # Save the presence-absence map as a PNG file
    plot_file <- file.path(species_pa_folder, paste("PA_Species", i, "_Rep", j, ".png", sep = ""))
    png(file = plot_file)
    plot(pa_map)
    dev.off()
    plot_file <- file.path(species_SuToProba_folder, paste("SuToProba_Species", i, "_Rep", j, ".png", sep = ""))
    png(file = plot_file)
    plotSuitabilityToProba(pa_map)
    dev.off()
    rds_file <- file.path(rds_folder, paste("PA_Species", i, "_Rep", j, ".rds", sep = ""))
    saveRDS(pa_map, file = rds_file)
  }
}


# Generate validation dataset
library(dplyr)
# Create an output folder to save the sampled points
sampled_points_folder <- "validation_points/Species3"
dir.create(sampled_points_folder, showWarnings = FALSE) 
env_raster <- raster("./EnvVar/princip.tif")
raster_brick <- brick("./EnvVar/princip.tif")
presence_num_list <- list(1000,2000,4000,8000)
# Loop for each virtual species
for (i in 1:4) { 
  rds_folder <- file.path("RDS", paste("Species", i))
  rds_file <- file.path(rds_folder, paste("PA_Species", i, "_Rep1", ".rds", sep = ""))
  pa_map <- readRDS(rds_file)
  for (j in 1:10){
    for (k in 1:4){
      set.seed(as.numeric(Sys.time()))
      # Sample presence and absence points using sampleOccurrences
      sampled <- sampleOccurrences(pa_map, 
                                   n = presence_num_list[[k]], # The number of points to sample
                                   type = "presence only",
                                   plot = FALSE,
                                   correct.by.suitability = TRUE)
      sampled_points_df <- sampled$sample.points[,c(1,2)]
      sampled_points_df$Species <- 'species1'
      sampled_points <- SpatialPointsDataFrame(coords = sampled_points_df[,c(1,2)], 
                                               proj4string = crs(raster_brick), data = sampled_points_df)
      shp_file <- file.path(sampled_points_folder, paste("Sp_Species", i, "_testRep", j, "_Presence", k, ".shp", sep = ""))
      shapefile(sampled_points,  shp_file, overwrite = TRUE)
    } 
  }
}


#Generate a training dataset
library(dplyr)
sampled_points_folder <- "TrainData/Species4"
dir.create(sampled_points_folder, showWarnings = FALSE)  # Create the folder if it doesn't exist
# Number of repetitions for sampling
env_raster <- raster("./EnvVar/princip.tif")
raster_brick <- brick("./EnvVar/princip.tif")
train_num <- 100
for (i in 4:4) {  # Assuming there are 6 species
  # Loop for each repetition
  rds_folder <- file.path("RDS", paste("Species", i))
  for (j in 1:10) {
      # Read the RDS file for this species and repetition
      rds_file <- file.path(rds_folder, paste("PA_Species", i, "_Rep1", ".rds", sep = ""))
      pa_map <- readRDS(rds_file)
      set.seed(as.numeric(Sys.time()))
      # Sample presence and absence points using sampleOccurrences
      sampled <- sampleOccurrences(
        pa_map,
        n = train_num,
        # The number of points to sample
        type = "presence only",
        plot = TRUE,
        correct.by.suitability = TRUE
      )
      sampled_points_df <- sampled$sample.points[, c(1, 2)]
      sampled_points_df$Species <- paste("Species", i)
      sampled_points <-
        SpatialPointsDataFrame(
          coords = sampled_points_df[, c(1, 2)],
          proj4string = crs(raster_brick),
          data = sampled_points_df
        )
      shp_file <-
        file.path(
          sampled_points_folder,
          paste("Sp_Species", i, "_Rep", j,".shp", sep = "")
        )
      shapefile(sampled_points,  shp_file, overwrite = TRUE)
  }
}


#Generate a testdataset
env_raster <- raster("./EnvVar/princip.tif")
total_pixels <- length(env_raster) 
for(sp_index in 3:3){
  for(rep_index in 1:1){
      dir.create(paste0("Test_dataset/Species",sp_index), showWarnings = FALSE) 
      testdataset_file <- file.path(paste0("Test_dataset/Species",sp_index),paste0("Species",sp_index,"_Rep",rep_index,".csv"))
      rdsfile <- file.path("RDS",paste0("Species ",sp_index),
                           paste0("PA_Species",sp_index,"_Rep1",".rds"))
      virtualsp <- readRDS(rdsfile)
      test_points <- sampleOccurrences(virtualsp,
                                       n = 50000,
                                       type = "presence-absence",
                                       sample.prevalence = 0.5,
                                       plot = FALSE)
      test_df <- test_points$sample.points[,c(1,2,3)]
      write.csv(test_df,testdataset_file,row.names = FALSE) 
  }
}

