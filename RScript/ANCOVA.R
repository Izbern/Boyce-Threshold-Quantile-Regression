library(dplyr)
library(ggbreak)
#read data
kappa_table <- data.frame()
for (j in 1:4){
  for (i in 1:4){
    df_name <- paste0("Kappa/Species", j, "/allKappa1/allKappa_Species", j, "_presence", i,".csv")
    df <- read.csv(df_name)
    kappa_table <- rbind(kappa_table, df)
  }
}
for (j in 1:4){
  df_name <- paste0("../Pandas_threshold/Kappa/allKappa/allKappa_pandas_presence", j, ".csv")
  df <- read.csv(df_name)[,c(-1)]
  kappa_table <- rbind(kappa_table, df)
}

kappa_table <- kappa_table %>%
  mutate(Method = case_when(
    Method == "MaxSens.Spec" ~ "MSS",
    Method == "Sens.Spec" ~ "ESS",
    Method == "PredPrev.Obs" ~ "EqualPrev",
    Method == "MinROCdist" ~ "MinROC",
    TRUE ~ Method  
  ))

#F test
methods <- levels(as.factor(kappa_table$Method))
models <- levels(as.factor(kappa_table$Model))
Sp <- levels(as.factor(kappa_table$Species))
method<- "BTQR"
Ftable <- data.frame()
for (species in Sp){
  for (method in methods){
    for (model in models){
      threshold_table <- kappa_table %>% filter(Model == model, Method == method, Species == species) 
      fit <- aov(Threshold~Presence, threshold_table)
      summary <- summary(fit)[[1]]
      Fvalue <- summary$`F value`[1]
      Pr <- summary$`Pr(>F)`[1]
      Ftable_row <- data.frame(Species = species,
                               Method = method,
                               Model = model,
                               Fvalue = Fvalue,
                               Pr = Pr)
      Ftable <- rbind(Ftable,Ftable_row)
    }
  }
}



  


