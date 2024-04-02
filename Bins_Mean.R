# This script computes the averages of the binary maps.
library(raster)
library(dplyr)
library(foreach)

spp <- read.table("./spp.csv", header = T, sep = ',') 
dim(spp)
head(spp, 10)

table(spp$sp)

especies <- unique(spp$sp)
especies
# for (especie in especies) {
foreach(especie = especies, .packages = c("raster", "dplyr")) %dopar% {
  ev = read.csv(paste0(
    "./outputs/",
    especie, "/", especie, "_Evaluation.all.csv"))[["TSS"]]
  
  # Presente #
  Cur = stack(paste0(
    "./outputs/",
    especie,
    ".Presente/",
    especie,
    "_Current.bin.tif"
  ))
  
  Cur_mean = mean(Cur)
  
  Cur_Wmean = weighted.mean(Cur, ev)
  
  writeRaster(Cur_mean, paste0(
    "./Bins_Mean/",
    especie,
    "_Current.mean.bin.tif"
  ), formato = "GTiff", overwrite = T)
  
  writeRaster(Cur_Wmean, paste0(
    "./Bins_Mean/",
    especie,
    "_Current.Wmean.bin.tif"
  ), formato = "GTiff")
  
  # Futuro 45 2030 #
  Fut.45.2030 = stack(
    paste0(
      "./outputs/",
      especie,
      ".RCP4.5.2030/",
      especie,
      "_Future.45.2030.bin.tif"
    )
  )
  
  Fut.45.2030_mean = mean(Fut.45.2030)
  
  Fut.45.2030_Wmean = weighted.mean(Fut.45.2030, ev)
  
  writeRaster(Fut.45.2030_mean, paste0(
    "./Bins_Mean/",
    especie,
    "_Future.45.2030.mean.bin.tif"
  ), formato = "GTiff", overwrite = T)
  
  writeRaster(Fut.45.2030_Wmean, paste0(
    "./Bins_Mean/",
    especie,
    "_Future.45.2030.Wmean.bin.tif"
  ), formato = "GTiff")
  
  # Futuro 85 2030 #
  Fut.85.2030 = stack(
    paste0(
      "./outputs/",
      especie,
      ".RCP8.5.2030/",
      especie,
      "_Future.85.2030.bin.tif"
    )
  )
  
  Fut.85.2030_mean = mean(Fut.85.2030)
  
  Fut.85.2030_Wmean = weighted.mean(Fut.85.2030, ev)
  
  writeRaster(Fut.85.2030_mean, paste0(
    "./Bins_Mean/",
    especie,
    "_Future.85.2030.mean.bin.tif"
  ), formato = "GTiff", overwrite = T)
  
  writeRaster(Fut.85.2030_Wmean, paste0(
    "./Bins_Mean/",
    especie,
    "_Future.85.2030.Wmean.bin.tif"
  ), formato = "GTiff")
  
  # Futuro 45 2050 #
  Fut.45.2050 = stack(
    paste0(
      "./outputs/",
      especie,
      ".RCP4.5.2050/",
      especie,
      "_Future.45.2050.bin.tif"
    )
  )
  
  Fut.45.2050_mean = mean(Fut.45.2050)
  
  Fut.45.2050_Wmean = weighted.mean(Fut.45.2050, ev)
  
  writeRaster(Fut.45.2050_mean, paste0(
    "./Bins_Mean/",
    especie,
    "_Future.45.2050.mean.bin.tif"
  ), formato = "GTiff", overwrite = T)
  
  writeRaster(Fut.45.2050_Wmean, paste0(
    "./Bins_Mean/",
    especie,
    "_Future.45.2050.Wmean.bin.tif"
  ), formato = "GTiff")
  
  # Futuro 85 2050 #
  Fut.85.2050 = stack(
    paste0(
      "./outputs/",
      especie,
      ".RCP8.5.2050/",
      especie,
      "_Future.85.2050.bin.tif"
    )
  )
  
  Fut.85.2050_mean = mean(Fut.85.2050)
  
  Fut.85.2050_Wmean = weighted.mean(Fut.85.2050, ev)
  
  writeRaster(Fut.85.2050_mean, paste0(
    "./Bins_Mean/",
    especie,
    "_Future.85.2050.mean.bin.tif"
  ), formato = "GTiff", overwrite = T)
  
  writeRaster(Fut.85.2050_Wmean, paste0(
    "./Bins_Mean/",
    especie,
    "_Future.85.2050.Wmean.bin.tif"
  ), formato = "GTiff")
  
  # Futuro 45 2070 #
  Fut.45.2070 = stack(
    paste0(
      "./outputs/",
      especie,
      ".RCP4.5.2070/",
      especie,
      "_Future.45.2070.bin.tif"
    )
  )
  
  Fut.45.2070_mean = mean(Fut.45.2070)
  
  Fut.45.2070_Wmean = weighted.mean(Fut.45.2070, ev)
  
  writeRaster(Fut.45.2070_mean, paste0(
    "./Bins_Mean/",
    especie,
    "_Future.45.2070.mean.bin.tif"
  ), formato = "GTiff", overwrite = T)
  
  writeRaster(Fut.45.2070_Wmean, paste0(
    "./Bins_Mean/",
    especie,
    "_Future.45.2070.Wmean.bin.tif"
  ), formato = "GTiff")
  
  # Futuro 85 2070 #
  Fut.85.2070 = stack(
    paste0(
      "./outputs/",
      especie,
      ".RCP8.5.2070/",
      especie,
      "_Future.85.2070.bin.tif"
    )
  )
  
  Fut.85.2070_mean = mean(Fut.85.2070)
  
  Fut.85.2070_Wmean = weighted.mean(Fut.85.2070, ev)
  
  writeRaster(Fut.85.2070_mean, paste0(
    "./Bins_Mean/",
    especie,
    "_Future.85.2070.mean.bin.tif"
  ), formato = "GTiff", overwrite = T)
  
  writeRaster(Fut.85.2070_Wmean, paste0(
    "./Bins_Mean/",
    especie,
    "_Future.85.2070.Wmean.bin.tif"
  ), formato = "GTiff")
  
  # Futuro 45 2090 #
  Fut.45.2090 = stack(
    paste0(
      "./outputs/",
      especie,
      ".RCP4.5.2090/",
      especie,
      "_Future.45.2090.bin.tif"
    )
  )
  
  Fut.45.2090_mean = mean(Fut.45.2090)
  
  Fut.45.2090_Wmean = weighted.mean(Fut.45.2090, ev)
  
  writeRaster(Fut.45.2090_mean, paste0(
    "./Bins_Mean/",
    especie,
    "_Future.45.2090.mean.bin.tif"
  ), formato = "GTiff", overwrite = T)
  
  writeRaster(Fut.45.2090_Wmean, paste0(
    "./Bins_Mean/",
    especie,
    "_Future.45.2090.Wmean.bin.tif"
  ), formato = "GTiff")
  
  # Futuro 85 2090 #
  Fut.85.2090 = stack(
    paste0(
      "./outputs/",
      especie,
      ".RCP8.5.2090/",
      especie,
      "_Future.85.2090.bin.tif"
    )
  )
  
  Fut.85.2090_mean = mean(Fut.85.2090)
  
  Fut.85.2090_Wmean = weighted.mean(Fut.85.2090, ev)
  
  writeRaster(Fut.85.2090_mean, paste0(
    "./Bins_Mean/",
    especie,
    "_Future.85.2090.mean.bin.tif"
  ), formato = "GTiff", overwrite = T)
  
  writeRaster(Fut.85.2090_Wmean, paste0(
    "./Bins_Mean/",
    especie,
    "_Future.85.2090.Wmean.bin.tif"
  ), formato = "GTiff")
  # remove all files #
  remove(list = ls()[c(grep(
    "Cur", as.factor(ls())))])
  gc(reset = TRUE, full = T)
  remove(list = ls()[c(grep(
    "Fut", as.factor(ls())))])
  gc(reset = TRUE, full = T)
}
