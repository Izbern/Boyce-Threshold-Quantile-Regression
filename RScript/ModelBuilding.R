library(terra)
library(biomod2)
library(tidyterra)
library(dismo)
library(dplyr)
library(modEvA)
library(chngpt)
library(irr)
library(readxl)
library(tidyr)
library(virtualspecies)
library(openxlsx)
library(PresenceAbsence)

# Function: Plot P/E curve
PE_Plot <- function(plot_folder, PETable_folder, p_df, proj_result, sp_i, sample_i, run_i, presence_num_i){
  PETable <- as.data.frame(matrix(ncol = 12, nrow = 100))
  names(PETable) <- c('GLM','GLM.ratio','GAM','GAM.ratio','GBM','GBM.ratio','MARS','MARS.ratio','MAXENT','MAXENT.ratio','RF','RF.ratio')
  for (i in 1:6){
    plot_file <- file.path(plot_folder, paste0(colnames(PETable)[2*i-1],"_Species",sp_i,"_Sample",sample_i,"_run",run_i,"_presence",presence_num_i,".png"))
    png(file = plot_file)
    BoyceTable<-Boyce(obs = p_df[,c(-3)], pred = proj_result[[i]], n.bins = NA,
                      bin.width = "default", res = 100, method = "spearman", rm.dup.classes = FALSE,
                      rm.dup.points = TRUE, plot = TRUE, plot.lines = TRUE, plot.values = TRUE,
                      plot.digits = 3, na.rm = TRUE)[[1]]
    dev.off()
    x <- nrow(BoyceTable)
    PETable[1:x, c(2*i-1)]<-BoyceTable[, c(2)]
    PETable[1:x, c(2*i)]<-BoyceTable[, c(7)]
  }
  write.csv(PETable, file =  file.path(PETable_folder, 
                                       paste0("PETable_Species",sp_i,
                                              "_Sample",sample_i,
                                              "_run",run_i,
                                              "_presence",presence_num_i, ".csv")), row.names = FALSE)
  return(PETable)
}

# Building species distribution models
env_raster <- raster("./EnvVar/princip.tif")
Env_variable <- rast("./EnvVar/princip.tif")
trainData_folder <- "TrainData/Species4"
for(sp_index in 4:4){
  for(sample_index in 1:10){
    trainData_name <-
      file.path(trainData_folder,
                paste0("Sp_Species", sp_index,
                       "_Rep", sample_index, ".shp"))
    trainData_points <- shapefile(trainData_name)
    train_presence_df <- trainData_points@data[, c(1, 2)]
    absence_number <- 
    train_pre_spvector <-
      vect(
        train_presence_df,
        geom = c("x", "y"),
        crs = crs(Env_variable, proj = TRUE)
      )
    SpeciesPA <- BIOMOD_FormatingData(
      resp.var = train_pre_spvector,
      expl.var = Env_variable,
      resp.name = paste0("Species", sp_index),
      PA.nb.rep = 1,
      PA.nb.absences = 200,
      # PA.strategy = 'disk',
      # PA.dist.min = 5000,
      # PA.strategy = "random",
      PA.strategy = 'sre',
      PA.sre.quant = 0.005
    )
    myBiomodOptions <-
      BIOMOD_ModelingOptions(
        # RF = list(
        #   do.classif = TRUE,
        #   ntrees = 1000,
        #   nodesize = 1,
        #   mtry = 2)
      )
    myBiomodModelOut <- BIOMOD_Modeling(
      bm.format = SpeciesPA,
      modeling.id = paste0("Sample", sample_index),
      models =  c("RF"),
      bm.options = myBiomodOptions,
      CV.strategy = "random",
      CV.perc = 0.2,
      do.full.models = FALSE,
    )
    myBiomodProjectionOut <- BIOMOD_Projection(
      bm.mod = myBiomodModelOut,
      proj.name = paste0("Sample", sample_index),
      new.env = Env_variable,
      models.chosen = 'all',
      binary.meth = NULL,
      compress = 'xz',
      on_0_1000 = FALSE,
      output.format = '.tif',
    )
  }
}
plot(myBiomodProjectionOut)


# Plot P/E curve
sp_i <- 4
pe_i <- 1
for(sp_index in sp_i:sp_i){
  model_folder <- paste0("Species", sp_index)
  PE_sp_folder <- paste0("PE_Data/Species", sp_index)
  dir.create(PE_sp_folder, showWarnings = FALSE) 
  PE_Plot_File <-  file.path(PE_sp_folder,paste0("PE_Plot",pe_i))
  dir.create(PE_Plot_File, showWarnings = FALSE) 
  PE_Table_File <-  file.path(PE_sp_folder,paste0("PE_Table",pe_i))
  dir.create(PE_Table_File, showWarnings = FALSE) 
  Validation_Data_File <- paste0("validation_points/Species", sp_index)
  for(sample_index in 1:10){
    proj_file <- file.path(model_folder, paste0("proj_Sample", sample_index),
                           paste0("Species",sp_index,".Sample",sample_index,".projection.out"))
    myBiomodProjectionOut <- get(load(proj_file))
    var_name <- paste0("Species", sp_index, ".Sample", sample_index,".projection.out")
    eval(parse(text = paste0("rm(", var_name, ")")))
    Proj_Raster <- terra::unwrap(myBiomodProjectionOut@proj.out@val)
    for(test_index in 1:10){
      for(presence_num_index in 1:4){
        remaining_points <- shapefile(file.path(Validation_Data_File,
                                                paste0("Sp_Species", sp_index,
                                                       "_testRep", test_index,"_Presence",presence_num_index,".shp")))
        remaining_df<-remaining_points@data[,c(1,2)]
        validation_p_df <- remaining_df
        model_name <- myBiomodProjectionOut@models.projected
        PETable <- as.data.frame(matrix(ncol = length(model_name)*2, nrow = 100))
        for(i in 1:length(model_name)){
          model_n <- tail(strsplit(model_name[i], "_")[[1]], n=1)
          plot_file <- file.path(PE_Plot_File, 
                                 paste0(model_n,"_Species",sp_index,"_Sample",sample_index,"_run",test_index,"_presence",presence_num_index,".png"))
          png(file = plot_file)
          BoyceTable<-Boyce(obs = validation_p_df, pred = Proj_Raster[[i]], n.bins = NA,
                            bin.width = "default", res = 100, method = "spearman", rm.dup.classes = FALSE,
                            rm.dup.points = TRUE, plot = TRUE, plot.lines = TRUE, plot.values = TRUE,
                            plot.digits = 3, na.rm = TRUE)[[1]]
          dev.off()
          x <- nrow(BoyceTable)
          PETable[1:x, c(2*i-1)]<-BoyceTable[, c(4)]
          PETable[1:x, c(2*i)]<-BoyceTable[, c(7)]
          names(PETable)[2*i-1]<-model_n
          names(PETable)[2*i]<-paste0(model_n,".ratio")
          # final_table <- cbind(origin_table, PETable)
        }
        write.csv(PETable, file =  file.path(PE_Table_File , 
                                             paste0("PETable_Species",sp_index,
                                                    "_Sample",sample_index,
                                                    "_run",test_index,
                                                    "_presence",presence_num_index, ".csv")), row.names = FALSE)
      } 
    }
  }
}




# Threshold value was calculated by the method except BTQR method
sp_i<-4
threshold_i<-1
for(sp_index in sp_i:sp_i){
  threshold_folder <- paste0("Threshold/Species", sp_index)
  dir.create(threshold_folder, showWarnings = FALSE) 
  threshold_folder <- paste0("Threshold/Species", sp_index, "/threshold", threshold_i)
  dir.create(threshold_folder, showWarnings = FALSE)  # Create the folder if it doesn't exist
  Validation_Data_File <- paste0("validation_points/Species", sp_index)
  Env_variable <- rast("./EnvVar/princip.tif")
  model_folder <- file.path(paste0("Species", sp_index))
  for(sample_index in 1:10){
    proj_file <- file.path(model_folder, paste0("proj_Sample", sample_index),
                           paste0("Species",sp_index,".Sample",sample_index,".projection.out"))
    myBiomodProjectionOut <- get(load(proj_file))
    var_name <- paste0("Species", sp_index, ".Sample", sample_index,".projection.out")
    eval(parse(text = paste0("rm(", var_name, ")")))
    Proj_Raster <- terra::unwrap(myBiomodProjectionOut@proj.out@val)
    for(test_index in 1:10){
      for(presence_num_index in 1:4){
        remaining_points <- shapefile(file.path(Validation_Data_File,
                                                paste0("Sp_Species", sp_index,
                                                       "_testRep", test_index,"_Presence",presence_num_index,".shp")))
        remaining_df<-remaining_points@data[,c(1,2)]
        absence_number <- 10000
        validation_p_df <- remaining_df
        validation_ab_spvector <- spatSample(Proj_Raster, size = absence_number, as.points = TRUE, values = TRUE,method = 'random')
        validation_ab <- terra::extract(Proj_Raster, as(validation_ab_spvector, "SpatVector"))[,-1]
        validation_p_spvector <- vect(validation_p_df, geom = c("x", "y"), crs = crs(Env_variable, proj = TRUE))
        validation_p <- terra::extract(Proj_Raster, as(validation_p_spvector, "SpatVector"))[,-1]
        validation_ab <- data.frame(observed = 0, validation_ab)
        validation_p <- data.frame(observed = 1, validation_p)
        names(validation_p)[2] <- "value"
        names(validation_ab)[2] <- "value"
        validation_data <- rbind(validation_p,validation_ab)
        validation_data <- data.frame(plotID = seq_len(nrow(validation_data)), validation_data)
        thresholdTable<-optimal.thresholds(validation_data,threshold = 800,opt.methods = c(2,3,4,6,8,9,10),req.sens=0.9)
        thresholdTable[8,]  <- optimal.thresholds(validation_data,threshold = 800,opt.methods = c(10),req.sens=0.75)
        thresholdTable[9,]  <- optimal.thresholds(validation_data,threshold = 800,opt.methods = c(10),req.sens=0.5)
        thresholdTable$Method[7:9] <- c("Se0.9", "Se0.75", "Se0.5")
        thresholdTable  <- t(thresholdTable)
        thresholdTable <- data.frame(thresholdTable)
        colnames(thresholdTable) <- thresholdTable[1, ]
        thresholdTable <- thresholdTable[-1, ]
        thresholdTable$BTQR <- 0
        origin_threshold <- read.xlsx(file.path(threshold_folder, 
                                                paste0("threshold_Species",sp_index,
                                                       "_Sample",sample_index,
                                                       "_run",test_index,
                                                       "_presence",presence_num_index, ".xlsx")))[-6,]
        names(thresholdTable)[c(1, 2, 4)] <- c("Sens.Spec","MaxSens.Spec","PredPrev.Obs")
        completed_threshold <- rbind(origin_threshold, thresholdTable)
        openxlsx::write.xlsx (completed_threshold, file = file.path(threshold_folder, 
                                                               paste0("threshold_Species",sp_index,
                                                                      "_Sample",sample_index,
                                                                      "_run",test_index,
                                                                      "_presence",presence_num_index, ".xlsx")), rowNames = FALSE)
      } 
    }
  }
}

#calculate Kappa
sp_i <- 4
kappa_i <- 1
threshold_i <- 1
for(sp_index in sp_i:sp_i){
  spkappafolder <- paste0("Kappa/Species", sp_index)
  dir.create(spkappafolder,showWarnings = FALSE)
  kappafolder <- file.path(spkappafolder,paste0("Kappa", kappa_i))
  dir.create(kappafolder,showWarnings = FALSE)
  thresholdfolder <- paste0("Threshold/Species", sp_index, "/threshold", threshold_i)
  test_df_file <- file.path(paste0("Test_dataset/Species",sp_index),paste0("Species",sp_index,"_Rep1", ".csv"))
  test_df <- read.csv(test_df_file)
  test_xy <- test_df[,c(1,2)]
  for(sample_index in 1:10){
    for(run_index in 1:10){
      for(presence_num_index in 1:4){
        threshold_file <- file.path(thresholdfolder,
                                    paste0("threshold_Species",sp_index, 
                                           "_Sample",sample_index,
                                           "_run",run_index,
                                           "_presence", presence_num_index, ".xlsx")) 
        kappa_file <- file.path(kappafolder,
                                paste0("kappa_Species",sp_index,
                                       "_Sample",sample_index,
                                       "_run",run_index,
                                       "_presence", presence_num_index, ".csv")) 
        threshold_table <- read_xlsx(threshold_file)
        # kappa_table <- threshold_table
        proj_file <- file.path(paste0("Species", sp_index), paste0("proj_Sample", sample_index),
                               paste0("Species",sp_index,".Sample",sample_index,".projection.out"))
        myBiomodProjectionOut <- get(load(proj_file))
        var_name <- paste0("Species", sp_index, ".Sample", sample_index, ".projection.out")
        eval(parse(text = paste0("rm(", var_name, ")")))
        Proj_Raster <- terra::unwrap(myBiomodProjectionOut@proj.out@val)
        test_points_value <- terra::extract(Proj_Raster, test_xy, xy=TRUE)
        kappa_table <- read.csv(kappa_file)
        threshold_list <- as.numeric(threshold_table[6,])
        for(i in 1:length(threshold_list)){
          threshold_value <- threshold_list[i]
          Method_name <- names(threshold_table)[i]
          if (is.na(threshold_value)){
            kappa_table$Kappa[kappa_table$Model=='RF' & kappa_table$Method==Method_name] <- NA
          } else {
            test_points_value$reclassified_value <- ifelse(test_points_value[,2] < threshold_list[i], 0, 1)
            kappa_table$Kappa[kappa_table$Model=='RF' & kappa_table$Method==Method_name] <- kappa2(data.frame(test_points_value$reclassified_value, test_df$Real))$value
            kappa_table$Threshold[kappa_table$Model=='RF' & kappa_table$Method==Method_name] <- threshold_value
          }
        }
        write.csv(kappa_table,kappa_file,row.names = FALSE)
      }
    }
  }
}