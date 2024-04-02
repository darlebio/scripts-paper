
###                 Acknowledgments                ##########################
### Dr. Matheus de Souza Lima-Ribeiro's team of Universidade Federal de Jataí.
### Dr. João Carlos Pires de Oliveira 

## Install and Loading packages ####
# install.packages("rgdal")
# install.packages("sp")
# install.packages("raster")
# install.packages("maps")
# install.packages("mnormt") #mnormt / psych
# install.packages("psych")
# #install.packages("permut")#permut / vegan ... permut desatualizado
# install.packages("vegan")
# install.packages("dismo")
# install.packages("kernlab")
# install.packages("rJava")
# install.packages("randomForest")
# install.packages("earth")
# install.packages("mgcv")
# install.packages("nnet")
# install.packages("beepr")
# install.packages("doParallel")
# install.packages("biomod2")
# install.packages("sdmvspecies")
# install.packages("filesstrings")
# install.packages("dplyr")
# install.packages("lhs")
# install.packages("ggplot2")
library(sp)
library(rgdal)
library(tcltk2)
library(raster)
library(maps)
library(psych)
library(vegan)
library(mnormt)
library(dismo)
library(rJava)
library(randomForest)
library(beepr)
library(doParallel)
library(biomod2)
library(sdmvspecies)
library(filesstrings)
library(dplyr)
library(lhs)
library(gam)
library(ggplot2)
# Criating output dir #
if (dir.exists("outputs") == F) {
  dir.create("outputs")
}

# Environmental Niche Ovelap ###

envOver = function (
    model.1,
    model.2,
    env,
    tolerance = 0.001,
    max.reps = 9999,
    cor.method = "spearman",
    chunk.size = 1e+5,
    recal.model.1 = NA,
    recal.model.2 = NA,
    verbose = FALSE){
  
  model.1 <- model.1
  model.2 <- model.2
  
  continue <- FALSE
  n.reps <- 1
  while (continue == FALSE & n.reps < max.reps) {
    gens <- chunk.size
    pred1 <- NA
    pred2 <- NA
    pred1.recal <- NA
    pred2.recal <- NA
    this.lhs <- randomLHS(chunk.size, length(names(env)))
    if (inherits(env,
                 c("raster", "RasterStack",
                   "RasterBrick", "RasterLayer"))) {
      mins <- minValue(env)
      maxes <- maxValue(env)
    }
    predict.table <- t(t(this.lhs) * (maxes - mins) + mins)
    colnames(predict.table) <- names(env)
    if (inherits(model.1, "DistModel")) {
      pred1 <- as.numeric(predict(model.1, x = data.frame(predict.table), 
                                  type = "response"))
    }
    else {
      if (inherits(model.1, "ranger")) {
        pred1 <-
          as.numeric(predict(
            model.1,
            data = data.frame(predict.table)
          )$predictions[, 2, drop = TRUE])
      }
      else {
        pred1 <-
          as.numeric(predict(
            model.1,
            newdata = data.frame(predict.table),
            type = "response"
          ))
      }
    }
    if (inherits(model.2, "DistModel")) {
      pred2 <- as.numeric(predict(model.2, x = data.frame(predict.table),
                                  type = "response"))
    }
    else {
      if (inherits(model.2, "ranger")) {
        pred2 <-
          as.numeric(predict(
            model.2,
            data = data.frame(predict.table)
          )$predictions[, 2, drop = TRUE])
      }
      else {
        pred2 <-
          as.numeric(predict(
            model.2,
            newdata = data.frame(predict.table), 
            type = "response"
          ))
      }
    }
    pred1[pred1 < 0] <- 0
    pred2[pred2 < 0] <- 0
    
    if (sd(pred1) == 0 | sd(pred2) == 0) {
      n.reps <- n.reps + 1
      next
    }
    
    this.d <- 1 - sum(abs(pred1 / sum(pred1) - pred2 / (sum(pred2)))) /
      2
    this.i <-
      1 - sum((sqrt(pred1 / sum(pred1)) - sqrt(pred2 / sum(pred2))) ^ 2) / 2
    this.cor <- cor(pred1, pred2, method = cor.method)
    if (!is.nan(this.d) & !is.nan(this.i)) {
      continue <- TRUE
    }
    else {
      n.reps <- n.reps + 1
    }
  }
  
  delta <- 1
  while (delta > tolerance) {
    this.lhs <- randomLHS(chunk.size, length(names(env)))
    predict.table <- t(t(this.lhs) * (maxes - mins) +
                         mins)
    colnames(predict.table) <- names(env)
    if (inherits(model.1, "DistModel")) {
      pred1 <- as.numeric(predict(
        model.1,
        x = data.frame(predict.table), type = "response"
      ))
    }
    else {
      if (inherits(model.1, "ranger")) {
        pred1 <-
          as.numeric(predict(
            model.1,
            data = data.frame(predict.table)
          )$predictions[, 2,
                        drop = TRUE])
      }
      else {
        pred1 <-
          as.numeric(predict(
            model.1,
            newdata = data.frame(predict.table), 
            type = "response"
          ))
      }
    }
    if (inherits(model.2, "DistModel")) {
      pred2 <- as.numeric(predict(
        model.2,
        x = data.frame(predict.table), type = "response"
      ))
    }
    else {
      if (inherits(model.2, "ranger")) {
        pred2 <-
          as.numeric(predict(
            model.2,
            data = data.frame(predict.table)
          )$predictions[, 2,
                        drop = TRUE])
      }
      else {
        pred2 <-
          as.numeric(predict(
            model.2,
            newdata = data.frame(predict.table),
            type = "response"
          ))
      }
    }
    pred1[pred1 < 0] <- 0
    pred2[pred2 < 0] <- 0
    if (sd(pred1) == 0 | sd(pred2) == 0) {
      next
    }
    else {
      gens <- c(gens, max(gens) + chunk.size)
    }
    
    n <- length(this.d)
    old.d <- this.d[n]
    new.d <- 1 - sum(abs(pred1 / sum(pred1) - pred2 / (sum(pred2)))) /
      2
    old.i <- this.i[n]
    new.i <-
      1 - sum((sqrt(pred1 / sum(pred1)) - sqrt(pred2 / sum(pred2))) ^ 2) / 2
    this.d <- c(this.d, old.d * (n / (n + 1)) + new.d *
                  1 / (n + 1))
    this.i <- c(this.i, old.i * (n / (n + 1)) + new.i *
                  1 / (n + 1))
    if (sd(pred1) == 0 | sd(pred2) == 0) {
      this.cor <- c(this.cor, NA)
    }
    else {
      old.cor <- this.cor[n]
      n <- n - length(which(is.na(this.cor)))
      new.cor <- cor(pred1, pred2, method = cor.method)
      this.cor <- c(this.cor, old.cor * (n / (n + 1)) +
                      new.cor * 1 / (n + 1))
    }
    delta <-
      max(c(abs(mean(this.d) - mean(this.d[-length(this.d)])),
            abs(mean(this.i) - mean(this.i[-length(this.i)])),
            abs(mean(this.cor) - mean(this.cor[-length(this.cor)]))),
          na.rm = TRUE)
  }
  
  output <- NA
  if (all(is.na(recal.model.1)) & all(is.na(recal.model.2))) {
    output <- list(
      env.D = mean(this.d),
      env.I = mean(this.i),
      env.cor = mean(this.cor)
    )
  }
  
  return(output)
}

# Function to predict future predictions ###
prefut = function(rast, rast1 = NULL, model, GCM =  NULL) {
  if(missing(GCM)){
    GCM = 1
  }else{GCM = GCM}
  pre = foreach::foreach(
    gcm = 1:length(GCM),
    .combine = stack,
    .packages = c("raster","biomod2", 'sp', "kernlab","dismo",
                  "stats","randomForest", "nnet", "earth", "adehabitatHS", 
                  "gam" )
  ) %dopar% {
    predict.enfa <-function (object.enfa,baseline.climate,new.climate,nf = 2,...) {
      m.baseline <- apply(slot(baseline.climate, "data"), 2, mean)
      sd.baseline <-apply(slot(baseline.climate, "data"), 2, sd)
      Zli <-object.enfa$li[, 1:nf]
      f1 <-function(x) {rep(x, object.enfa$pr)}
      Sli <- apply(Zli, 2, f1)
      m <- apply(Sli, 2, mean)
      cov <-t(as.matrix(Sli)) %*% as.matrix(Sli) / nrow(Sli)
      if (!missing("new.climate")) {
        new.climate.scale <- sweep(slot(new.climate, "data"), 2, m.baseline)
        new.climate.scale <-
          as.matrix(new.climate.scale) %*% diag(1 / sd.baseline)
        Zli <-
          new.climate.scale %*% as.matrix(object.enfa$co)}
      maha <- mahalanobis(Zli, center = m, cov = cov)
      map <-rasterize(data.frame(new.climate@coords),rast[[1]],maha) * -1
      return(invisible(map))
    }
    if (!"madifa" %in% class(model)) {
      raster::predict(rast[gcm][[1]], model)
    } else if ("madifa" %in% class(model) && "list" %in% class(rast1)) {
      climaPres.spdf = na.omit(data.frame(xyFromCell(rast, 1:ncell(rast)),
                                          raster::values(rast)))
      climaFut.spdf = na.omit(data.frame(xyFromCell(rast1[gcm][[1]], 
                                                    1:ncell(rast1[gcm][[1]])),
                                         raster::values(rast1[gcm][[1]])))
      suppressWarnings(gridded(climaPres.spdf) <-~ x + y)
      suppressWarnings(gridded(climaFut.spdf) <- ~ x + y)
      predict.enfa( object.enfa = model,baseline.climate =climaPres.spdf,
                    new.climate = climaFut.spdf)}else{
                      climaPres.spdf = na.omit(data.frame(xyFromCell(rast, 1:ncell(rast)),
                                                          raster::values(rast)))
                      climaFut.spdf = na.omit(data.frame(xyFromCell(rast1, 
                                                                    1:ncell(rast1)),
                                                         raster::values(rast1)))
                      suppressWarnings(gridded(climaPres.spdf) <-~ x + y)
                      suppressWarnings(gridded(climaFut.spdf) <- ~ x + y)
                      predict.enfa( object.enfa = model,baseline.climate =climaPres.spdf,
                                    new.climate = climaFut.spdf)}
  }
  pre = mean(pre)
  return(pre)
  rm(new.climate.scale,maha,map,climaPres.spdf,climaFut.spdf)
  gc(reset = T, full = T)
}

# Function to scale maps -- DON'T CHANGE###
rescx = function(x){
  
  if(!require(dplyr)) install.packages("raster")
  library(dplyr)
  if(!require(dplyr)) install.packages("sp")
  library(dplyr)
  if(!require(dplyr)) install.packages("vegan")
  library(vegan)
  # if(!require(dplyr)) install.packages("reshape")
  # library(reshape)
  
  nb = nlayers(x)
  df = na.exclude(as.data.frame(rasterToPoints(x)))
  nomes = colnames(df[,-c(1,2)])
  ID = rep(nomes, each = nrow(df))
  co = df[,1:2]
  df1 = data.frame(ID = 1:nrow(df), df[,-c(1,2)])
  
  df1_lon = reshape(df1, direction="long",
                    idvar = "ID",
                    varying = colnames(df1[-1]),
                    v.names = c("layer"))
  
  df1_lon$layer = vegan::decostand(df1_lon$layer, MARGIN = 2, method = "range")
  
  df2 = reshape(df1_lon, direction = "wide", idvar = "ID",
                timevar = "time", v.names = "layer", sep= "_")
  df2 = df2[,-1]
  colnames(df2) = colnames(df1[,-1])
  df3 = data.frame(co, df2)
  # res = rasterFromXYZ(cbind(co, df2))
  # res = raster::rasterize(co, x[[1]], df2)
  suppressWarnings(sp::gridded(df3) <-~ x + y)
  res = stack(df3)
  return(res)
  
}
# Function to Evaluate All Models --- DON'T CHANGE  ####

eval.All.Model <- function(rast, dismoTestPrepared, tr ) {{
  p = raster::extract(rast, dismoTestPrepared[dismoTestPrepared[, 3] == 1, 1:2 ])
  a = raster::extract(rast, dismoTestPrepared[dismoTestPrepared[, 3] == 0, 1:2 ])}
  p <- stats::na.omit(p)
  a <- stats::na.omit(a)
  np <- length(p)
  na <- length(a)
  if (na == 0 | np == 0) {
    stop('cannot evaluate a model without absence and presence data that are not NA')
  }
  if (missing(tr)) {
    if (length(p) > 1000) {
      tr <- as.vector(quantile(p, 0:1000/1000))
    } else {
      tr <- p}
    if (length(a) > 1000) {
      tr <- c(tr, as.vector(quantile(a, 0:1000/1000)))
    } else {
      tr <- c(tr, a)}
    tr <- sort(unique( round(tr, 8)))
    tr <- c( tr - 0.0001, tr[length(tr)] + c(0, 0.0001))
  } else {
    tr <- sort(as.vector(tr))}
  N <- na + np
  xc <- new('ModelEvaluation')
  xc@presence = p
  xc@absence = a
  R <- sum(rank(c(p, a))[1:np]) - (np*(np+1)/2)
  xc@auc <- R / (as.numeric(na) * as.numeric(np))
  cr <- try( cor.test(c(p,a), c(rep(1, length(p)), rep(0, length(a))) ), 
             silent=TRUE )
  if (class(cr) != 'try-error') {
    xc@cor <- cr$estimate
    xc@pcor <- cr$p.value}
  res <- matrix(ncol=4, nrow=length(tr))
  colnames(res) <- c('tp', 'fp', 'fn', 'tn')
  xc@t <- tr
  for (i in 1:length(tr)) {
    res[i,1] <- length(p[p>=tr[i]])  # a  true positives
    res[i,2] <- length(a[a>=tr[i]])  # b  false positives
    res[i,3] <- length(p[p<tr[i]])    # c  false negatives
    res[i,4] <- length(a[a<tr[i]])}    # d  true negatives
  xc@confusion = res
  a = res[,1]
  b = res[,2]
  c = res[,3]
  d = res[,4]
  # after Fielding and Bell	
  xc@np <- as.integer(np)
  xc@na <- as.integer(na)
  xc@prevalence = (a + c) / N
  xc@ODP = (b + d) / N
  xc@CCR = (a + d) / N
  xc@TPR = a / (a + c)
  xc@TNR = d / (b + d)
  xc@FPR = b / (b + d)
  xc@FNR = c/(a + c)
  xc@PPP = a/(a + b)
  xc@NPP = d/(c + d)
  xc@MCR = (b + c)/N
  xc@OR = (a*d)/(c*b)
  prA = (a+d)/N
  prY = (a+b)/N * (a+c)/N
  prN = (c+d)/N * (b+d)/N
  prE = prY + prN
  xc@kappa = (prA - prE) / (1-prE)
  return(xc)}

quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 

### Whenever necessary:
# Parallel processing #
detectCores()
getDoParWorkers()
cl <-
  parallel::makeCluster(5, outfile = paste0("./outputs/", "log_models.log"))
registerDoParallel(cl)
getDoParWorkers()

## CHANGE raster TEMPORARY FILE DIRECTORY
## define the name of a temp directory where raster tmp files will be stored
raster_tmp_dir <- "raster_tmp"
if (dir.exists("raster_tmp") == F) {
  ## create the directory (only when starting modeling):
  dir.create(raster_tmp_dir,
             showWarnings = F,
             recursive = T)
}
#
### set raster options
rasterOptions(tmpdir = raster_tmp_dir)


# Loading variaveis ####

bio.crop<-
  list.files(
    "./Vars/Mundo/Presente/PCA/",  pattern = ".grd$",
    full.names = TRUE
  )
bio.crop
bio.crop <- raster::stack(bio.crop)
names(bio.crop)<- c("PCA1", "PCA2", "PCA3","PCA4", "PCA5","PCA6")
names(bio.crop)

#--------------------------------------------------#
#      LOADING FUTURE VARIABLES                ####
#------------------------------------------------#
### 2030 ####
###GCM 1: CanESM5

bio70_CA_45_2030 <-
  list.files(
    "./vars/Mundo/Future/2030/RCP45/CanESM5/ssp245/PCA/",
    pattern = ".grd$",
    full.names = TRUE
  )
bio70_CA_45_2030
bio70_CA_45_2030 <- raster::stack(bio70_CA_45_2030)
bio70_CA_45_2030
names(bio70_CA_45_2030) <-
  c("PCA1", "PCA2", "PCA3","PCA4", "PCA5","PCA6")
names(bio70_CA_45_2030)


###GCM 2: CNRM-CM6-1

bio70_CN_45_2030 <-
  list.files(
    "./vars/Mundo/Future/2030/RCP45/CNRM-CM6-1/ssp245/PCA/",
    pattern = ".grd$",
    full.names = TRUE
  )
bio70_CN_45_2030
bio70_CN_45_2030 <- raster::stack(bio70_CN_45_2030)
bio70_CN_45_2030
names(bio70_CN_45_2030) <-
  c("PCA1", "PCA2", "PCA3","PCA4", "PCA5","PCA6")
names(bio70_CN_45_2030)

###GCM 3: MIROC-ES2L

bio70_MI_45_2030 <-
  list.files(
    "./vars/Mundo/Future/2030/RCP45/MIROC-ES2L/ssp245/PCA/",
    pattern = ".grd$",
    full.names = TRUE
  )
bio70_MI_45_2030
bio70_MI_45_2030 <- raster::stack(bio70_MI_45_2030)
bio70_MI_45_2030
names(bio70_MI_45_2030) <-
  c("PCA1", "PCA2", "PCA3","PCA4", "PCA5","PCA6")
names(bio70_MI_45_2030)

# Loading raster layer with future projections 85 ###
#------------------------------------------------#

###GCM 1: CanESM5

bio70_CA_85_2030 <-
  list.files(
    "./vars/Mundo/Future/2030/RCP85/CanESM5/ssp585/PCA/",
    pattern = ".grd$",
    full.names = TRUE
  )
bio70_CA_85_2030
bio70_CA_85_2030 <- raster::stack(bio70_CA_85_2030)
bio70_CA_85_2030
names(bio70_CA_85_2030) <-
  c("PCA1", "PCA2", "PCA3","PCA4", "PCA5","PCA6")
names(bio70_CA_85_2030)


###GCM 2: CNRM-CM6-1

bio70_CN_85_2030 <-
  list.files(
    "./vars/Mundo/Future/2030/RCP85/CNRM-CM6-1/ssp585/PCA/",
    pattern = ".grd$",
    full.names = TRUE
  )
bio70_CN_85_2030
bio70_CN_85_2030 <- raster::stack(bio70_CN_85_2030)
bio70_CN_85_2030
names(bio70_CN_85_2030) <-
  c("PCA1", "PCA2", "PCA3","PCA4", "PCA5","PCA6")
names(bio70_CN_85_2030)

###GCM 3: MIROC-ES2L

bio70_MI_85_2030 <-
  list.files(
    "./vars/Mundo/Future/2030/RCP85/MIROC-ES2L/ssp585/PCA/",
    pattern = ".grd$",
    full.names = TRUE
  )
bio70_MI_85_2030
bio70_MI_85_2030 <- raster::stack(bio70_MI_85_2030)
bio70_MI_85_2030
names(bio70_MI_85_2030) <-
  c("PCA1", "PCA2", "PCA3","PCA4", "PCA5","PCA6")
names(bio70_MI_85_2030)

scenario.list.45_2030 <- list(
  bio70_CA_45_2030,
  bio70_CN_45_2030,
  bio70_MI_45_2030
)

scenario.list.85_2030 <- list(
  bio70_CA_85_2030,
  bio70_CN_85_2030,
  bio70_MI_85_2030
)

rm( bio70_CA_45_2030,
    bio70_CN_45_2030,
    bio70_MI_45_2030,
    bio70_CA_85_2030,
    bio70_CN_85_2030,
    bio70_MI_85_2030)

### 2050 ####
###GCM 1: CanESM5

bio70_CA_45_2050 <-
  list.files(
    "./vars/Mundo/Future/2050/RCP45/CanESM5/ssp245/PCA/",
    pattern = ".grd$",
    full.names = TRUE
  )
bio70_CA_45_2050
bio70_CA_45_2050 <- raster::stack(bio70_CA_45_2050)
bio70_CA_45_2050
names(bio70_CA_45_2050) <-
  c("PCA1", "PCA2", "PCA3","PCA4", "PCA5","PCA6")
names(bio70_CA_45_2050)


###GCM 2: CNRM-CM6-1

bio70_CN_45_2050 <-
  list.files(
    "./vars/Mundo/Future/2050/RCP45/CNRM-CM6-1/ssp245/PCA/",
    pattern = ".grd$",
    full.names = TRUE
  )
bio70_CN_45_2050
bio70_CN_45_2050 <- raster::stack(bio70_CN_45_2050)
bio70_CN_45_2050
names(bio70_CN_45_2050) <-
  c("PCA1", "PCA2", "PCA3","PCA4", "PCA5","PCA6")
names(bio70_CN_45_2050)

###GCM 3: MIROC-ES2L

bio70_MI_45_2050 <-
  list.files(
    "./vars/Mundo/Future/2050/RCP45/MIROC-ES2L/ssp245/PCA/",
    pattern = ".grd$",
    full.names = TRUE
  )
bio70_MI_45_2050
bio70_MI_45_2050 <- raster::stack(bio70_MI_45_2050)
bio70_MI_45_2050
names(bio70_MI_45_2050) <-
  c("PCA1", "PCA2", "PCA3","PCA4", "PCA5","PCA6")
names(bio70_MI_45_2050)

# Loading raster layer with future projections 85 ###
#------------------------------------------------#

###GCM 1: CanESM5

bio70_CA_85_2050 <-
  list.files(
    "./vars/Mundo/Future/2050/RCP85/CanESM5/ssp585/PCA/",
    pattern = ".grd$",
    full.names = TRUE
  )
bio70_CA_85_2050
bio70_CA_85_2050 <- raster::stack(bio70_CA_85_2050)
bio70_CA_85_2050
names(bio70_CA_85_2050) <-
  c("PCA1", "PCA2", "PCA3","PCA4", "PCA5","PCA6")
names(bio70_CA_85_2050)


###GCM 2: CNRM-CM6-1

bio70_CN_85_2050 <-
  list.files(
    "./vars/Mundo/Future/2050/RCP85/CNRM-CM6-1/ssp585/PCA/",
    pattern = ".grd$",
    full.names = TRUE
  )
bio70_CN_85_2050
bio70_CN_85_2050 <- raster::stack(bio70_CN_85_2050)
bio70_CN_85_2050
names(bio70_CN_85_2050) <-
  c("PCA1", "PCA2", "PCA3","PCA4", "PCA5","PCA6")
names(bio70_CN_85_2050)

###GCM 3: MIROC-ES2L

bio70_MI_85_2050 <-
  list.files(
    "./vars/Mundo/Future/2050/RCP85/MIROC-ES2L/ssp585/PCA/",
    pattern = ".grd$",
    full.names = TRUE
  )
bio70_MI_85_2050
bio70_MI_85_2050 <- raster::stack(bio70_MI_85_2050)
bio70_MI_85_2050
names(bio70_MI_85_2050) <-
  c("PCA1", "PCA2", "PCA3","PCA4", "PCA5","PCA6")
names(bio70_MI_85_2050)

scenario.list.45_2050 <- list(
  bio70_CA_45_2050,
  bio70_CN_45_2050,
  bio70_MI_45_2050
)

scenario.list.85_2050 <- list(
  bio70_CA_85_2050,
  bio70_CN_85_2050,
  bio70_MI_85_2050
)

rm( bio70_CA_45_2050,
    bio70_CN_45_2050,
    bio70_MI_45_2050,
    bio70_CA_85_2050,
    bio70_CN_85_2050,
    bio70_MI_85_2050)

### 2070 ####
###GCM 1: CanESM5

bio70_CA_45_2070 <-
  list.files(
    "./vars/Mundo/Future/2070/RCP45/CanESM5/ssp245/PCA/",
    pattern = ".grd$",
    full.names = TRUE
  )
bio70_CA_45_2070
bio70_CA_45_2070 <- raster::stack(bio70_CA_45_2070)
bio70_CA_45_2070
names(bio70_CA_45_2070) <-
  c("PCA1", "PCA2", "PCA3","PCA4", "PCA5","PCA6")
names(bio70_CA_45_2070)


###GCM 2: CNRM-CM6-1

bio70_CN_45_2070 <-
  list.files(
    "./vars/Mundo/Future/2070/RCP45/CNRM-CM6-1/ssp245/PCA/",
    pattern = ".grd$",
    full.names = TRUE
  )
bio70_CN_45_2070
bio70_CN_45_2070 <- raster::stack(bio70_CN_45_2070)
bio70_CN_45_2070
names(bio70_CN_45_2070) <-
  c("PCA1", "PCA2", "PCA3","PCA4", "PCA5","PCA6")
names(bio70_CN_45_2070)

###GCM 3: MIROC-ES2L

bio70_MI_45_2070 <-
  list.files(
    "./vars/Mundo/Future/2070/RCP45/MIROC-ES2L/ssp245/PCA/",
    pattern = ".grd$",
    full.names = TRUE
  )
bio70_MI_45_2070
bio70_MI_45_2070 <- raster::stack(bio70_MI_45_2070)
bio70_MI_45_2070
names(bio70_MI_45_2070) <-
  c("PCA1", "PCA2", "PCA3","PCA4", "PCA5","PCA6")
names(bio70_MI_45_2070)

# Loading raster layer with future projections 85 ###
#------------------------------------------------#

###GCM 1: CanESM5

bio70_CA_85_2070 <-
  list.files(
    "./vars/Mundo/Future/2070/RCP85/CanESM5/ssp585/PCA/",
    pattern = ".grd$",
    full.names = TRUE
  )
bio70_CA_85_2070
bio70_CA_85_2070 <- raster::stack(bio70_CA_85_2070)
bio70_CA_85_2070
names(bio70_CA_85_2070) <-
  c("PCA1", "PCA2", "PCA3","PCA4", "PCA5","PCA6")
names(bio70_CA_85_2070)


###GCM 2: CNRM-CM6-1

bio70_CN_85_2070 <-
  list.files(
    "./vars/Mundo/Future/2070/RCP85/CNRM-CM6-1/ssp585/PCA/",
    pattern = ".grd$",
    full.names = TRUE
  )
bio70_CN_85_2070
bio70_CN_85_2070 <- raster::stack(bio70_CN_85_2070)
bio70_CN_85_2070
names(bio70_CN_85_2070) <-
  c("PCA1", "PCA2", "PCA3","PCA4", "PCA5","PCA6")
names(bio70_CN_85_2070)

###GCM 3: MIROC-ES2L

bio70_MI_85_2070 <-
  list.files(
    "./vars/Mundo/Future/2070/RCP85/MIROC-ES2L/ssp585/PCA/",
    pattern = ".grd$",
    full.names = TRUE
  )
bio70_MI_85_2070
bio70_MI_85_2070 <- raster::stack(bio70_MI_85_2070)
bio70_MI_85_2070
names(bio70_MI_85_2070) <-
  c("PCA1", "PCA2", "PCA3","PCA4", "PCA5","PCA6")
names(bio70_MI_85_2070)

scenario.list.45_2070 <- list(
  bio70_CA_45_2070,
  bio70_CN_45_2070,
  bio70_MI_45_2070
)

scenario.list.85_2070 <- list(
  bio70_CA_85_2070,
  bio70_CN_85_2070,
  bio70_MI_85_2070
)

rm( bio70_CA_45_2070,
    bio70_CN_45_2070,
    bio70_MI_45_2070,
    bio70_CA_85_2070,
    bio70_CN_85_2070,
    bio70_MI_85_2070)

### 2090 ####
###GCM 1: CanESM5

bio70_CA_45_2090 <-
  list.files(
    "./vars/Mundo/Future/2090/RCP45/CanESM5/ssp245/PCA/",
    pattern = ".grd$",
    full.names = TRUE
  )
bio70_CA_45_2090
bio70_CA_45_2090 <- raster::stack(bio70_CA_45_2090)
bio70_CA_45_2090
names(bio70_CA_45_2090) <-
  c("PCA1", "PCA2", "PCA3","PCA4", "PCA5","PCA6")
names(bio70_CA_45_2090)


###GCM 2: CNRM-CM6-1

bio70_CN_45_2090 <-
  list.files(
    "./vars/Mundo/Future/2090/RCP45/CNRM-CM6-1/ssp245/PCA/",
    pattern = ".grd$",
    full.names = TRUE
  )
bio70_CN_45_2090
bio70_CN_45_2090 <- raster::stack(bio70_CN_45_2090)
bio70_CN_45_2090
names(bio70_CN_45_2090) <-
  c("PCA1", "PCA2", "PCA3","PCA4", "PCA5","PCA6")
names(bio70_CN_45_2090)

###GCM 3: MIROC-ES2L

bio70_MI_45_2090 <-
  list.files(
    "./vars/Mundo/Future/2090/RCP45/MIROC-ES2L/ssp245/PCA/",
    pattern = ".grd$",
    full.names = TRUE
  )
bio70_MI_45_2090
bio70_MI_45_2090 <- raster::stack(bio70_MI_45_2090)
bio70_MI_45_2090
names(bio70_MI_45_2090) <-
  c("PCA1", "PCA2", "PCA3","PCA4", "PCA5","PCA6")
names(bio70_MI_45_2090)

# Loading raster layer with future projections 85 ###
#------------------------------------------------#

###GCM 1: CanESM5

bio70_CA_85_2090 <-
  list.files(
    "./vars/Mundo/Future/2090/RCP85/CanESM5/ssp585/PCA/",
    pattern = ".grd$",
    full.names = TRUE
  )
bio70_CA_85_2090
bio70_CA_85_2090 <- raster::stack(bio70_CA_85_2090)
bio70_CA_85_2090
names(bio70_CA_85_2090) <-
  c("PCA1", "PCA2", "PCA3","PCA4", "PCA5","PCA6")
names(bio70_CA_85_2090)


###GCM 2: CNRM-CM6-1

bio70_CN_85_2090 <-
  list.files(
    "./vars/Mundo/Future/2090/RCP85/CNRM-CM6-1/ssp585/PCA/",
    pattern = ".grd$",
    full.names = TRUE
  )
bio70_CN_85_2090
bio70_CN_85_2090 <- raster::stack(bio70_CN_85_2090)
bio70_CN_85_2090
names(bio70_CN_85_2090) <-
  c("PCA1", "PCA2", "PCA3","PCA4", "PCA5","PCA6")
names(bio70_CN_85_2090)

###GCM 3: MIROC-ES2L

bio70_MI_85_2090 <-
  list.files(
    "./vars/Mundo/Future/2090/RCP85/MIROC-ES2L/ssp585/PCA/",
    pattern = ".grd$",
    full.names = TRUE
  )
bio70_MI_85_2090
bio70_MI_85_2090 <- raster::stack(bio70_MI_85_2090)
bio70_MI_85_2090
names(bio70_MI_85_2090) <-
  c("PCA1", "PCA2", "PCA3","PCA4", "PCA5","PCA6")
names(bio70_MI_85_2090)

scenario.list.45_2090 <- list(
  bio70_CA_45_2090,
  bio70_CN_45_2090,
  bio70_MI_45_2090
)

scenario.list.85_2090 <- list(
  bio70_CA_85_2090,
  bio70_CN_85_2090,
  bio70_MI_85_2090
)

rm( bio70_CA_45_2090,
    bio70_CN_45_2090,
    bio70_MI_45_2090,
    bio70_CA_85_2090,
    bio70_CN_85_2090,
    bio70_MI_85_2090)

##### Occurrences ####
# Select you species data matrix # 
spp <- read.table("./spp.csv", header = T, sep = ',')  

# plot all your occurence points #
# plot(bio.crop[[1]]) 
# points(spp[,c("lon","lat")], pch = 19)


dim(spp)
head(spp, 10)

table(spp$sp)

especies <- unique(spp$sp)
especies


# Pseudo-absence Set

PAs <- 1
RUNs = 10

# 2030 #
GCM45s_2030 <- c("bio70_CA_45_2030", "bio70_CN_45_2030", "bio70_MI_45_2030")
GCM85s_2030 <- c("bio70_CA_85_2030", "bio70_CN_85_2030", "bio70_MI_85_2030")
# 2050 #
GCM45s_2050 <- c("bio70_CA_45_2050", "bio70_CN_45_2050", "bio70_MI_45_2050")
GCM85s_2050 <- c("bio70_CA_85_2050", "bio70_CN_85_2050", "bio70_MI_85_2050")
# 2070 #
GCM45s_2070 <- c("bio70_CA_45_2070", "bio70_CN_45_2070", "bio70_MI_45_2070")
GCM85s_2070 <- c("bio70_CA_85_2070", "bio70_CN_85_2070", "bio70_MI_85_2070")
# 2090 #
GCM45s_2090 <- c("bio70_CA_45_2090", "bio70_CN_45_2090", "bio70_MI_45_2090")
GCM85s_2090 <- c("bio70_CA_85_2090", "bio70_CN_85_2090", "bio70_MI_85_2090")

# especie = especies[1]
# Species Loop ####
# For sequential loop (One species) ###
# for (especie in especies) {

## For species in parallel ###
foreach(especie = especies, # For parallel looping (Multiple Species)
        .packages = c("raster", "biomod2",'sp',"sdmvspecies", "filesstrings",
                      "rgdal","maps","mnormt","kernlab","dismo","doParallel",
                      "stats","rJava","randomForest","lhs","psych", "earth",
                      "gam", "adehabitatHS"
        ),
        .verbose = F,
        .errorhandling = "stop") %do% {
          
          if (dir.exists(paste0("./temp_output/",especie,"/")) == F) {
            dir.create(paste0("./temp_output/"))
            dir.create(paste0("./temp_output/",especie,"/"))
          }
          
          if (dir.exists(paste0("./Models/",especie,"/")) == F) {
            dir.create(paste0("./Models/"))
            dir.create(paste0("./Models/",especie,"/"))
          }
          
          ini1 = Sys.time()
          print(paste0("Starting", " ", especie, " " ,"modeling"))
          
          cl1 <-
            parallel::makeCluster(3)
          registerDoParallel(cl1)
          
          PA = 1
          # Creating empty objects to store results ###
          bioclim.e<-
            # maha.e <-
            enfa.e <-
            glm.e <- 
            gam.e <- 
            maxent.e <-
            rf.e <- NULL
          
          
          occs <- spp[spp$sp == especie, c("lon", "lat")]
          
          # Data checking and prepation
          ocor.val <- raster::extract(bio.crop, occs, cellnumbers = T)
          
          sum(is.na(ocor.val[, 1]))
          
          # plot(predictors[[1]], colNA = "red")
          
          # points(ocor.all[, -1], pch = 19, cex = 0.5)
          
          ocor.val <- cbind(occs, ocor.val)
          
          ocor.val <- na.omit(ocor.val)
          
          id <- duplicated(ocor.val[, "cells"]) # Checking dumplicate points
          
          sum(id == T)
          
          ocor <-
            ocor.val[id == F, c("lon", "lat")] # Removing duplicate points
          
          #------------------------------------------#
          #           SELECT PAs                  ###
          #----------------------------------------#
          
          try({
            coord1 = ocor
            sp::coordinates(coord1) <- ~ lon + lat
            raster::crs(coord1) <- raster::crs(bio.crop)
            
            dist.mean <- mean(sp::spDists(
              x = coord1,
              longlat = T,
              segments = FALSE
            ))
            dist.min = 5
            dist.min <-  min(sp::spDists(
              x = coord1,
              longlat = T,
              segments = F
            ))
            dist.min = 5
            
            # write.table(
            #   c(dist.min, dist.mean),
            #   paste0('./outputs/',
            #          especie, "_","dist", ".csv"),
            #   row.names = F,
            #   sep = ";"
            # )
          })
          
          PA.number <- nrow(ocor)
          PA.number #número de pontos de ocorrência espacialmente únicos
          
          # write.csv(data.frame(sp = especie, ocor), row.names = F, 
          #           paste0(especie, ".csv"))
          
          diretorio = paste0("Occurrence.", especie)
          
          #### Loop PAs ###
          # for (PA in seq(PAs)) {
          
          p = rep(1, times = nrow(ocor))
          
          invisible(capture.output(sel.PA <- biomod2::BIOMOD_FormatingData(
            resp.var = p,
            expl.var = raster::stack(bio.crop),
            resp.xy = ocor,
            resp.name = diretorio,
            PA.nb.rep = 1,
            PA.nb.absences = PA.number,
            PA.strategy = "disk",
            PA.dist.min = dist.min * 1000,
            PA.dist.max = dist.mean * 1000,
            na.rm = TRUE
          )))
          
          li.p <- grep("pa", rownames(sel.PA@coord))
          
          
          pa <- sel.PA@coord[li.p,]
          
          invisible(capture.output(sel.back <- biomod2::BIOMOD_FormatingData(
            resp.var = p,
            expl.var = raster::stack(bio.crop),
            resp.xy = ocor,
            resp.name = diretorio,
            PA.nb.rep = 1,
            PA.nb.absences = 10000,
            PA.strategy = "disk",
            PA.dist.min = dist.min * 1000,
            PA.dist.max = dist.mean * 1000,
            na.rm = TRUE
          )))
          
          li.b <- grep("pa", rownames(sel.back@coord))
          
          
          back <- sel.back@coord[li.b,]
          
          
          rm(sel.back, sel.PA)
          # set.seed(0)
          areaToral <- nrow(rasterToPoints(bio.crop))
          
          ### Loop RUN ####
          # for (RUN in seq(RUNs)) {
          
          # List Models #
          bioclim.models = list()
          # maha.models = list()
          enfa.models = list()
          glm.models = list()
          gam.models = list()
          maxent.models = list()
          rf.models = list()
          
          # Train/Test Data #
          
          training.pa_list = list()
          training.b_list = list()
          test.pa_list = list()
          
          for (RUN in 1:RUNs) {
            
            print(paste0(especie," ","PASet", PA," ","RUN", RUN))
            
            # Separating test/ training data
            
            # prepare data 70/30
            id.training.pa <-
              sample(1:nrow(ocor), round(0.7 * nrow(ocor), 0)) 
            # prepare data 70/30
            id.training.b <-
              sample(1:nrow(back), round(0.7 * nrow(back), 0)) 
            
            training.b <-
              na.omit(dismo::prepareData(bio.crop, p = ocor[id.training.pa, ], 
                                         b = back[id.training.b,], xy = T))
            # head(training.b)
            # tail(training.b)
            training.pa <-
              na.omit(dismo::prepareData(bio.crop, p = ocor[id.training.pa, ], 
                                         b = pa[id.training.pa, ], xy = T))
            test.pa <-
              na.omit(dismo::prepareData(bio.crop, p = ocor[-id.training.pa, ], 
                                         b = pa[-id.training.pa, ], xy = T))
            training.b_list[[RUN]] = training.b
            training.pa_list[[RUN]] = training.pa
            test.pa_list[[RUN]] = test.pa
            ## Modeling ####
            
            ### Bioclim  ####
            print(paste0(especie, " ","Modeling_Bioclim"," ",  
                         "PA",  PA," ","RUN", RUN))
            
            bioclim_model <-
              bioclim(x = training.b[training.b[, "pb"] == 1, -c(1:3)], )
            bioclim.models[[RUN]] = bioclim_model
            # Save Model Object #
            saveRDS(bioclim_model,
                    paste0("./Models/",
                           especie, "/", especie,
                           "_bioclim_",
                           "PA_",PA,
                           "_",
                           "RUN_",RUN,
                           ".rds"))
            
            ### ENFA  ####
            print(paste0(especie, " ","Modeling_enfa"," ",
                         "PA",  PA," ","RUN", RUN))
            
            climaPres <- raster::stack(bio.crop)
            # names(climaPres)# <- paste0("PC",1:6)
            
            climaPres.values <- raster::values(climaPres)
            climaPres.spdf <- na.omit(data.frame(xyFromCell(climaPres, 1:ncell(climaPres)), 
                                                 climaPres.values))
            
            suppressWarnings(gridded(climaPres.spdf) <- ~x+y)
            climaPres <- raster::stack(climaPres.spdf)
            climaPres.values <- raster::values(climaPres)
            media.climaPres <- apply(slot(climaPres.spdf, "data"), 2, mean)
            sd.climaPres <- apply(slot(climaPres.spdf, "data"), 2, sd)
            climaPres.scale<- sweep(slot(climaPres.spdf, "data"),2, media.climaPres)
            climaPres.scale<- as.matrix(climaPres.scale) %*% diag(1/sd.climaPres)
            
            #adjustment of the ENFA model
            
            pr.cell <- raster::extract(climaPres, 
                                       training.pa[training.pa[,"pb"]==1,1:2], 
                                       cellnumber=T)
            pr <- data.frame(pr= rep(0, ncell(climaPres)), climaPres.values)
            pr[pr.cell[,"cells"], 1] <- 1
            pr <- na.omit(pr)
            pr <- pr[,1]
            enfa_model <- adehabitatHS::madifa(ade4::dudi.pca(climaPres.scale, 
                                                              center=F, scale=F, scannf=F), 
                                               pr, scannf=F)
            
            
            enfa.models[[RUN]] = enfa_model
            # Save Model Object #
            saveRDS(enfa_model,
                    paste0("./Models/",
                           especie, "/", especie,
                           "_enfa_",
                           "PA_",PA,
                           "_",
                           "RUN_",RUN,
                           ".rds"))
            
            ### GLM  ####
            print(paste0(especie, " ","Modeling_GLM"," ", 
                         "PA",  PA," ","RUN", RUN))
            
            pb <- training.pa$pb
            
            glm_model <- glm(
              formula = pb ~ .,
              data = training.pa[, -c(1:3)],
              type = 'quadratic',
              interaction.level = 0,
              myFormula = NULL,
              test = 'AIC',
              family = binomial(link = 'logit'),
              control = glm.control(
                epsilon = 1e-08,
                maxit = 50,
                trace = FALSE
              )
            )
            glm.models[[RUN]] = glm_model
            # Save Model Object #
            saveRDS(glm_model, 
                    paste0("./Models/", 
                           especie, "/", especie,
                           "_glm_",
                           "PA_",PA,
                           "_",
                           "RUN_",RUN,
                           ".rds"))
            
            ### GAM ####
            print(paste0(especie, " ","Modeling_GAM"," ", 
                         "PA",  PA," ","RUN", RUN))
            
            mf = formula(paste0("pb~", 
                                paste0("s(",
                                       colnames(training.pa[,-c(1:3)]),")",
                                       collapse = "+")))
            gam_model <-
              suppressWarnings(gam::gam(mf,
                                        data = training.pa[,-c(1, 2)],
                                        family=binomial(link = "logit")
              ))
            gam.models[[RUN]] = gam_model
            # Save Model Object #
            saveRDS(gam_model, 
                    paste0("./Models/",
                           especie, "/", especie,
                           "_gam_",
                           "PA_",PA,
                           "_",
                           "RUN_",RUN,
                           ".rds"))
            
            ### Maxent #### -> ku.enm ####
            print(paste0(especie, " ","Modeling_MAXENT"," ", 
                         "PA",  PA," ","RUN", RUN))
            
            maxent_model <- quiet(maxent(
              x = training.b[, -c(1:3)],
              p = training.b[, 'pb'],
              memory_allocated = 512,
              background_data_dir = 'default',
              maximumbackground = 'default',
              maximumiterations = 200,
              visible = FALSE,
              linear = TRUE,
              quadratic = TRUE,
              product = TRUE,
              threshold = TRUE,
              hinge = TRUE,
              lq2lqptthreshold = 80,
              l2lqthreshold = 10,
              hingethreshold = 15,
              beta_threshold = -1,
              beta_categorical = -1,
              beta_lqp = -1,
              beta_hinge = -1,
              betamultiplier = 1,
              defaultprevalence = 0.5
            ))
            maxent.models[[RUN]] = maxent_model
            # Save Model Object #
            saveRDS(maxent_model, 
                    paste0("./Models/",
                           especie, "/", especie,
                           "_maxent_",
                           "PA_",PA,
                           "_",
                           "RUN_",RUN,
                           ".rds"))
            
            ### Random Forest ####
            print(paste0(especie, " ","Modeling_Random Forest"," ", 
                         "PA",  PA," ","RUN", RUN))
            
            suppressWarnings(rf_model <-
                               randomForest::randomForest(
                                 pb ~ .,
                                 data = training.pa[,-c(1:2)],
                                 ntree = 500,
                                 nodesize = 5, importance = T
                               ))
            rf.models[[RUN]] = rf_model
            # Save Model Object #
            saveRDS(rf_model, 
                    paste0("./Models/",
                           especie, "/", especie,
                           "_rf_",
                           "PA_",PA,
                           "_",
                           "RUN_",RUN,
                           ".rds"))
            
          }
          saveRDS(training.b_list, 
                  paste0("./temp_output/", especie,"/",
                         especie,"_training.b_list.rds"))
          
          saveRDS(training.pa_list, 
                  paste0("./temp_output/", especie,"/",
                         especie,"_training.pa_list.rds"))
          
          saveRDS(test.pa_list, 
                  paste0("./temp_output/", especie,"/",
                         especie,"_test.pa_list.rds"))
          
          # Building Preojections ####
          
          for (RUN in seq(RUNs)) {
            print(paste0(especie, " ","Proj_Bioclim"," ",
                         "PA",  PA," ","RUN", RUN))
            ##### bioclim ###
            # Current #
            bioclim_Cur <-
              raster::predict(bio.crop, bioclim.models[[RUN]],
                              type = "response")
            
            saveRDS(
              bioclim_Cur,
              paste0("./temp_output/",especie,"/", "bioclim_Cur","PA_",PA,
                     "RUN_",RUN,".rds"))
            
            
            # Furure 2030 45 #
            bioclim_Fut.45_2030 <- prefut(rast = scenario.list.45_2030,
                                          model =   bioclim.models[[RUN]],
                                          GCM = GCM45s_2030)
            
            saveRDS(
              bioclim_Fut.45_2030,
              paste0("./temp_output/",especie,"/","bioclim_Fut.45_2030",
                     "PA_",PA,"RUN_",RUN,".rds"))
            
            
            # Furure 2030 85 #
            bioclim_Fut.85_2030 <- prefut(rast = scenario.list.85_2030,
                                          model =   bioclim.models[[RUN]],
                                          GCM = GCM85s_2030)
            
            saveRDS(
              bioclim_Fut.85_2030,
              paste0("./temp_output/",especie,"/","bioclim_Fut.85_2030",
                     "PA_",PA,"RUN_",RUN,".rds"))
            
            # Furure 2050 45 #
            bioclim_Fut.45_2050 <- prefut(rast = scenario.list.45_2050,
                                          model =   bioclim.models[[RUN]],
                                          GCM = GCM45s_2050)
            
            saveRDS(
              bioclim_Fut.45_2050,
              paste0("./temp_output/",especie,"/","bioclim_Fut.45_2050",
                     "PA_",PA,"RUN_",RUN,".rds"))
            
            
            # Furure 2050 85 #
            bioclim_Fut.85_2050 <- prefut(rast = scenario.list.85_2050,
                                          model =   bioclim.models[[RUN]],
                                          GCM = GCM85s_2050)
            
            saveRDS(
              bioclim_Fut.85_2050,
              paste0("./temp_output/",especie,"/","bioclim_Fut.85_2050",
                     "PA_",PA,"RUN_",RUN,".rds"))
            
            # Furure 2070 45 #
            bioclim_Fut.45_2070 <- prefut(rast = scenario.list.45_2070,
                                          model =   bioclim.models[[RUN]],
                                          GCM = GCM45s_2070)
            
            saveRDS(
              bioclim_Fut.45_2070,
              paste0("./temp_output/",especie,"/","bioclim_Fut.45_2070",
                     "PA_",PA,"RUN_",RUN,".rds"))
            
            
            # Furure 2070 85 #
            bioclim_Fut.85_2070 <- prefut(rast = scenario.list.85_2070,
                                          model =   bioclim.models[[RUN]],
                                          GCM = GCM85s_2070)
            
            saveRDS(
              bioclim_Fut.85_2070,
              paste0("./temp_output/",especie,"/","bioclim_Fut.85_2070",
                     "PA_",PA,"RUN_",RUN,".rds"))
            
            # Furure 2090 45 #
            bioclim_Fut.45_2090 <- prefut(rast = scenario.list.45_2090,
                                          model =   bioclim.models[[RUN]],
                                          GCM = GCM45s_2090)
            
            saveRDS(
              bioclim_Fut.45_2090,
              paste0("./temp_output/",especie,"/","bioclim_Fut.45_2090",
                     "PA_",PA,"RUN_",RUN,".rds"))
            
            
            # Furure 2090 85 #
            bioclim_Fut.85_2090 <- prefut(rast = scenario.list.85_2090,
                                          model =   bioclim.models[[RUN]],
                                          GCM = GCM85s_2090)
            
            saveRDS(
              bioclim_Fut.85_2090,
              paste0("./temp_output/",especie,"/","bioclim_Fut.85_2090",
                     "PA_",PA,"RUN_",RUN,".rds"))
            

            #### enfa ###
            print(paste0(especie," ","Proj_ENFA"," ",
                         "PA", PA," ","RUN",RUN))
            # Current #
            enfa_Cur <- prefut(rast = bio.crop,rast1 = bio.crop,
                               model =  enfa.models[[RUN]])
            
            saveRDS(
              enfa_Cur,
              paste0("./temp_output/",especie,"/", "enfa_Cur","PA_",
                     PA,"RUN_",RUN,".rds"))
            
            
            # Furure 2030 45 #
            enfa_Fut.45_2030 <- prefut(rast = bio.crop,
                                       rast1 = scenario.list.45_2030,
                                       model =   enfa.models[[RUN]],
                                       GCM = GCM45s_2030)
            
            saveRDS(
              enfa_Fut.45_2030,
              paste0("./temp_output/",especie,"/","enfa_Fut.45_2030",
                     "PA_",PA,"RUN_",RUN,".rds"))
            
            # Furure 2030 85 #
            enfa_Fut.85_2030 <- prefut(rast = bio.crop,
                                       rast1 = scenario.list.85_2030,
                                       model =   enfa.models[[RUN]],
                                       GCM = GCM85s_2030)
            
            saveRDS(
              enfa_Fut.85_2030,
              paste0("./temp_output/",especie,"/","enfa_Fut.85_2030",
                     "PA_",PA,"RUN_",RUN,".rds"))
            
            # Furure 2050 45 #
            enfa_Fut.45_2050 <- prefut(rast = bio.crop,
                                       rast1 = scenario.list.45_2050,
                                       model =   enfa.models[[RUN]],
                                       GCM = GCM45s_2050)
            
            saveRDS(
              enfa_Fut.45_2050,
              paste0("./temp_output/",especie,"/","enfa_Fut.45_2050",
                     "PA_",PA,"RUN_",RUN,".rds"))
            
            # Furure 2050 85 #
            enfa_Fut.85_2050 <- prefut(rast = bio.crop,
                                       rast1 = scenario.list.85_2050,
                                       model =   enfa.models[[RUN]],
                                       GCM = GCM85s_2050)
            
            saveRDS(
              enfa_Fut.85_2050,
              paste0("./temp_output/",especie,"/","enfa_Fut.85_2050",
                     "PA_",PA,"RUN_",RUN,".rds"))
            
            
            # Furure 2070 45 #
            enfa_Fut.45_2070 <- prefut(rast = bio.crop,
                                       rast1 = scenario.list.45_2070,
                                       model =   enfa.models[[RUN]],
                                       GCM = GCM45s_2070)
            
            saveRDS(
              enfa_Fut.45_2070,
              paste0("./temp_output/",especie,"/","enfa_Fut.45_2070",
                     "PA_",PA,"RUN_",RUN,".rds"))
            
            # Furure 2070 85 #
            enfa_Fut.85_2070 <- prefut(rast = bio.crop,
                                       rast1 = scenario.list.85_2070,
                                       model =   enfa.models[[RUN]],
                                       GCM = GCM85s_2070)
            
            saveRDS(
              enfa_Fut.85_2070,
              paste0("./temp_output/",especie,"/","enfa_Fut.85_2070",
                     "PA_",PA,"RUN_",RUN,".rds"))
            
            # Furure 2090 45 #
            enfa_Fut.45_2090 <- prefut(rast = bio.crop,
                                       rast1 = scenario.list.45_2090,
                                       model =   enfa.models[[RUN]],
                                       GCM = GCM45s_2090)
            
            saveRDS(
              enfa_Fut.45_2090,
              paste0("./temp_output/",especie,"/","enfa_Fut.45_2090",
                     "PA_",PA,"RUN_",RUN,".rds"))
            
            # Furure 2090 85 #
            enfa_Fut.85_2090 <- prefut(rast = bio.crop,
                                       rast1 = scenario.list.85_2090,
                                       model =   enfa.models[[RUN]],
                                       GCM = GCM85s_2090)
            
            saveRDS(
              enfa_Fut.85_2090,
              paste0("./temp_output/",especie,"/","enfa_Fut.85_2090",
                     "PA_",PA,"RUN_",RUN,".rds"))
            
            ##### glm ###
            print(paste0(especie," ","Proj_GLM"," ",
                         "PA", PA," ","RUN",RUN))
            # Current #
            glm_Cur <-
              raster::predict(bio.crop, glm.models[[RUN]],
                              type = "response")
            
            saveRDS(
              glm_Cur,
              paste0("./temp_output/",especie,"/", "glm_Cur","PA_",PA,
                     "RUN_",RUN,".rds"))
            
            
            # Furure 2030 45 #
            glm_Fut.45_2030 <- prefut(rast = scenario.list.45_2030,
                                      model =   glm.models[[RUN]],
                                      GCM = GCM45s_2030)
            
            saveRDS(
              glm_Fut.45_2030,
              paste0("./temp_output/",especie,"/","glm_Fut.45_2030",
                     "PA_",PA,"RUN_",RUN,".rds"))
            
            
            # Furure 2030 85 #
            glm_Fut.85_2030 <- prefut(rast = scenario.list.85_2030,
                                      model =   glm.models[[RUN]],
                                      GCM = GCM85s_2030)
            
            saveRDS(
              glm_Fut.85_2030,
              paste0("./temp_output/",especie,"/","glm_Fut.85_2030",
                     "PA_",PA,"RUN_",RUN,".rds"))
            
            # Furure 2050 45 #
            glm_Fut.45_2050 <- prefut(rast = scenario.list.45_2050,
                                      model =   glm.models[[RUN]],
                                      GCM = GCM45s_2050)
            
            saveRDS(
              glm_Fut.45_2050,
              paste0("./temp_output/",especie,"/","glm_Fut.45_2050",
                     "PA_",PA,"RUN_",RUN,".rds"))
            
            
            # Furure 2050 85 #
            glm_Fut.85_2050 <- prefut(rast = scenario.list.85_2050,
                                      model =   glm.models[[RUN]],
                                      GCM = GCM85s_2050)
            
            saveRDS(
              glm_Fut.85_2050,
              paste0("./temp_output/",especie,"/","glm_Fut.85_2050",
                     "PA_",PA,"RUN_",RUN,".rds"))
            
            # Furure 2070 45 #
            glm_Fut.45_2070 <- prefut(rast = scenario.list.45_2070,
                                      model =   glm.models[[RUN]],
                                      GCM = GCM45s_2070)
            
            saveRDS(
              glm_Fut.45_2070,
              paste0("./temp_output/",especie,"/","glm_Fut.45_2070",
                     "PA_",PA,"RUN_",RUN,".rds"))
            
            
            # Furure 2070 85 #
            glm_Fut.85_2070 <- prefut(rast = scenario.list.85_2070,
                                      model =   glm.models[[RUN]],
                                      GCM = GCM85s_2070)
            
            saveRDS(
              glm_Fut.85_2070,
              paste0("./temp_output/",especie,"/","glm_Fut.85_2070",
                     "PA_",PA,"RUN_",RUN,".rds"))
            
            # Furure 2090 45 #
            glm_Fut.45_2090 <- prefut(rast = scenario.list.45_2090,
                                      model =   glm.models[[RUN]],
                                      GCM = GCM45s_2090)
            
            saveRDS(
              glm_Fut.45_2090,
              paste0("./temp_output/",especie,"/","glm_Fut.45_2090",
                     "PA_",PA,"RUN_",RUN,".rds"))
            
            
            # Furure 2090 85 #
            glm_Fut.85_2090 <- prefut(rast = scenario.list.85_2090,
                                      model =   glm.models[[RUN]],
                                      GCM = GCM85s_2090)
            
            saveRDS(
              glm_Fut.85_2090,
              paste0("./temp_output/",especie,"/","glm_Fut.85_2090",
                     "PA_",PA,"RUN_",RUN,".rds"))
            
            ##### gam ###
            print(paste0(especie," ","Proj_GAM"," ",
                         "PA", PA," ","RUN",RUN))
            # Current #
            gam_Cur <-
              raster::predict(bio.crop, gam.models[[RUN]],
                              type = "response")
            
            saveRDS(
              gam_Cur,
              paste0("./temp_output/",especie,"/", "gam_Cur","PA_",PA,
                     "RUN_",RUN,".rds"))
            
            
            # Furure 2030 45 #
            gam_Fut.45_2030 <- prefut(rast = scenario.list.45_2030,
                                      model =   gam.models[[RUN]],
                                      GCM = GCM45s_2030)
            
            saveRDS(
              gam_Fut.45_2030,
              paste0("./temp_output/",especie,"/","gam_Fut.45_2030",
                     "PA_",PA,"RUN_",RUN,".rds"))
            
            
            # Furure 2030 85 #
            gam_Fut.85_2030 <- prefut(rast = scenario.list.85_2030,
                                      model =   gam.models[[RUN]],
                                      GCM = GCM85s_2030)
            
            saveRDS(
              gam_Fut.85_2030,
              paste0("./temp_output/",especie,"/","gam_Fut.85_2030",
                     "PA_",PA,"RUN_",RUN,".rds"))
            
            # Furure 2050 45 #
            gam_Fut.45_2050 <- prefut(rast = scenario.list.45_2050,
                                      model =   gam.models[[RUN]],
                                      GCM = GCM45s_2050)
            
            saveRDS(
              gam_Fut.45_2050,
              paste0("./temp_output/",especie,"/","gam_Fut.45_2050",
                     "PA_",PA,"RUN_",RUN,".rds"))
            
            
            # Furure 2050 85 #
            gam_Fut.85_2050 <- prefut(rast = scenario.list.85_2050,
                                      model =   gam.models[[RUN]],
                                      GCM = GCM85s_2050)
            
            saveRDS(
              gam_Fut.85_2050,
              paste0("./temp_output/",especie,"/","gam_Fut.85_2050",
                     "PA_",PA,"RUN_",RUN,".rds"))
            
            # Furure 2070 45 #
            gam_Fut.45_2070 <- prefut(rast = scenario.list.45_2070,
                                      model =   gam.models[[RUN]],
                                      GCM = GCM45s_2070)
            
            saveRDS(
              gam_Fut.45_2070,
              paste0("./temp_output/",especie,"/","gam_Fut.45_2070",
                     "PA_",PA,"RUN_",RUN,".rds"))
            
            
            # Furure 2070 85 #
            gam_Fut.85_2070 <- prefut(rast = scenario.list.85_2070,
                                      model =   gam.models[[RUN]],
                                      GCM = GCM85s_2070)
            
            saveRDS(
              gam_Fut.85_2070,
              paste0("./temp_output/",especie,"/","gam_Fut.85_2070",
                     "PA_",PA,"RUN_",RUN,".rds"))
            
            # Furure 2090 45 #
            gam_Fut.45_2090 <- prefut(rast = scenario.list.45_2090,
                                      model =   gam.models[[RUN]],
                                      GCM = GCM45s_2090)
            
            saveRDS(
              gam_Fut.45_2090,
              paste0("./temp_output/",especie,"/","gam_Fut.45_2090",
                     "PA_",PA,"RUN_",RUN,".rds"))
            
            
            # Furure 2090 85 #
            gam_Fut.85_2090 <- prefut(rast = scenario.list.85_2090,
                                      model =   gam.models[[RUN]],
                                      GCM = GCM85s_2090)
            
            saveRDS(
              gam_Fut.85_2090,
              paste0("./temp_output/",especie,"/","gam_Fut.85_2090",
                     "PA_",PA,"RUN_",RUN,".rds"))
            
            ##### maxent ###
            print(paste0(especie," ","Proj_Maxent"," ",
                         "PA", PA," ","RUN",RUN))
            # Current #
            maxent_Cur <-
              quiet(raster::predict(bio.crop, maxent.models[[RUN]],
                                    type = "response"))
            
            saveRDS(
              maxent_Cur,
              paste0("./temp_output/",especie,"/", "maxent_Cur","PA_",PA,
                     "RUN_",RUN,".rds"))
            
            
            # Furure 2030 45 #
            maxent_Fut.45_2030 <- quiet(prefut(rast = scenario.list.45_2030,
                                               model =   maxent.models[[RUN]],
                                               GCM = GCM45s_2030))
            
            saveRDS(
              maxent_Fut.45_2030,
              paste0("./temp_output/",especie,"/","maxent_Fut.45_2030",
                     "PA_",PA,"RUN_",RUN,".rds"))
            
            
            # Furure 2030 85 #
            maxent_Fut.85_2030 <- quiet(prefut(rast = scenario.list.85_2030,
                                               model =   maxent.models[[RUN]],
                                               GCM = GCM85s_2030))
            
            saveRDS(
              maxent_Fut.85_2030,
              paste0("./temp_output/",especie,"/","maxent_Fut.85_2030",
                     "PA_",PA,"RUN_",RUN,".rds"))
            
            # Furure 2050 45 #
            maxent_Fut.45_2050 <- quiet(prefut(rast = scenario.list.45_2050,
                                               model =   maxent.models[[RUN]],
                                               GCM = GCM45s_2050))
            
            saveRDS(
              maxent_Fut.45_2050,
              paste0("./temp_output/",especie,"/","maxent_Fut.45_2050",
                     "PA_",PA,"RUN_",RUN,".rds"))
            
            
            # Furure 2050 85 #
            maxent_Fut.85_2050 <- quiet(prefut(rast = scenario.list.85_2050,
                                               model =   maxent.models[[RUN]],
                                               GCM = GCM85s_2050))
            
            saveRDS(
              maxent_Fut.85_2050,
              paste0("./temp_output/",especie,"/","maxent_Fut.85_2050",
                     "PA_",PA,"RUN_",RUN,".rds"))
            
            # Furure 2070 45 #
            maxent_Fut.45_2070 <- quiet(prefut(rast = scenario.list.45_2070,
                                               model =   maxent.models[[RUN]],
                                               GCM = GCM45s_2070))
            
            saveRDS(
              maxent_Fut.45_2070,
              paste0("./temp_output/",especie,"/","maxent_Fut.45_2070",
                     "PA_",PA,"RUN_",RUN,".rds"))
            
            
            # Furure 2070 85 #
            maxent_Fut.85_2070 <- quiet(prefut(rast = scenario.list.85_2070,
                                               model =   maxent.models[[RUN]],
                                               GCM = GCM85s_2070))
            
            saveRDS(
              maxent_Fut.85_2070,
              paste0("./temp_output/",especie,"/","maxent_Fut.85_2070",
                     "PA_",PA,"RUN_",RUN,".rds"))
            
            # Furure 2090 45 #
            maxent_Fut.45_2090 <- quiet(prefut(rast = scenario.list.45_2090,
                                               model =   maxent.models[[RUN]],
                                               GCM = GCM45s_2090))
            
            saveRDS(
              maxent_Fut.45_2090,
              paste0("./temp_output/",especie,"/","maxent_Fut.45_2090",
                     "PA_",PA,"RUN_",RUN,".rds"))
            
            
            # Furure 2090 85 #
            maxent_Fut.85_2090 <- quiet(prefut(rast = scenario.list.85_2090,
                                               model =   maxent.models[[RUN]],
                                               GCM = GCM85s_2090))
            
            saveRDS(
              maxent_Fut.85_2090,
              paste0("./temp_output/",especie,"/","maxent_Fut.85_2090",
                     "PA_",PA,"RUN_",RUN,".rds"))
            
            ##### rf ###
            print(paste0(especie," ","Proj_Random Forest"," ",
                         "PA", PA," ","RUN",RUN))
            # Current #
            rf_Cur <-
              raster::predict(bio.crop, rf.models[[RUN]],
                              type = "response")
            
            saveRDS(
              rf_Cur,
              paste0("./temp_output/",especie,"/", "rf_Cur","PA_",PA,
                     "RUN_",RUN,".rds"))
            
            
            # Furure 2030 45 #
            rf_Fut.45_2030 <- prefut(rast = scenario.list.45_2030,
                                     model =   rf.models[[RUN]],
                                     GCM = GCM45s_2030)
            
            saveRDS(
              rf_Fut.45_2030,
              paste0("./temp_output/",especie,"/","rf_Fut.45_2030",
                     "PA_",PA,"RUN_",RUN,".rds"))
            
            
            # Furure 2030 85 #
            rf_Fut.85_2030 <- prefut(rast = scenario.list.85_2030,
                                     model =   rf.models[[RUN]],
                                     GCM = GCM85s_2030)
            
            saveRDS(
              rf_Fut.85_2030,
              paste0("./temp_output/",especie,"/","rf_Fut.85_2030",
                     "PA_",PA,"RUN_",RUN,".rds"))
            
            # Furure 2050 45 #
            rf_Fut.45_2050 <- prefut(rast = scenario.list.45_2050,
                                     model =   rf.models[[RUN]],
                                     GCM = GCM45s_2050)
            
            saveRDS(
              rf_Fut.45_2050,
              paste0("./temp_output/",especie,"/","rf_Fut.45_2050",
                     "PA_",PA,"RUN_",RUN,".rds"))
            
            
            # Furure 2050 85 #
            rf_Fut.85_2050 <- prefut(rast = scenario.list.85_2050,
                                     model =   rf.models[[RUN]],
                                     GCM = GCM85s_2050)
            
            saveRDS(
              rf_Fut.85_2050,
              paste0("./temp_output/",especie,"/","rf_Fut.85_2050",
                     "PA_",PA,"RUN_",RUN,".rds"))
            
            # Furure 2070 45 #
            rf_Fut.45_2070 <- prefut(rast = scenario.list.45_2070,
                                     model =   rf.models[[RUN]],
                                     GCM = GCM45s_2070)
            
            saveRDS(
              rf_Fut.45_2070,
              paste0("./temp_output/",especie,"/","rf_Fut.45_2070",
                     "PA_",PA,"RUN_",RUN,".rds"))
            
            
            # Furure 2070 85 #
            rf_Fut.85_2070 <- prefut(rast = scenario.list.85_2070,
                                     model =   rf.models[[RUN]],
                                     GCM = GCM85s_2070)
            
            saveRDS(
              rf_Fut.85_2070,
              paste0("./temp_output/",especie,"/","rf_Fut.85_2070",
                     "PA_",PA,"RUN_",RUN,".rds"))
            
            # Furure 2090 45 #
            rf_Fut.45_2090 <- prefut(rast = scenario.list.45_2090,
                                     model =   rf.models[[RUN]],
                                     GCM = GCM45s_2090)
            
            saveRDS(
              rf_Fut.45_2090,
              paste0("./temp_output/",especie,"/","rf_Fut.45_2090",
                     "PA_",PA,"RUN_",RUN,".rds"))
            
            
            # Furure 2090 85 #
            rf_Fut.85_2090 <- prefut(rast = scenario.list.85_2090,
                                     model =   rf.models[[RUN]],
                                     GCM = GCM85s_2090)
            
            saveRDS(
              rf_Fut.85_2090,
              paste0("./temp_output/",especie,"/","rf_Fut.85_2090",
                     "PA_",PA,"RUN_",RUN,".rds"))
            
          }
          
          # Scaling ####
          { # bioclim #
            bioclim_Cur_stack = stack()
            for(PA in 1:PAs){
              for (RUN in 1:RUNs) {
                bioclim_Cur_stack = stack(bioclim_Cur_stack, readRDS(list.files(
                  paste0("./temp_output/", especie, "/"),
                  paste0('bioclim_Cur',"PA_",PA,"RUN_",RUN,".rds"),
                  full.names = T
                )))
              }
            }
            
            names(bioclim_Cur_stack)<- paste0("bioclim.Cur.", 
                                              1:nlayers(bioclim_Cur_stack))
            
            bioclim_Fut.45_2030_stack = stack()
            for(PA in 1:PAs){
              for (RUN in 1:RUNs) {
                bioclim_Fut.45_2030_stack = stack(bioclim_Fut.45_2030_stack, 
                                                  readRDS(list.files(
                                                    paste0("./temp_output/", especie, "/"),
                                                    paste0('bioclim_Fut.45_2030',"PA_",PA,"RUN_",RUN,".rds"),
                                                    full.names = T
                                                  )))
              }
            }
            
            names(bioclim_Fut.45_2030_stack)<- paste0("bioclim_Fut_2030_45.", 
                                                      1:nlayers(bioclim_Fut.45_2030_stack))
            
            bioclim_Fut.85_2030_stack = stack()
            for(PA in 1:PAs){
              for (RUN in 1:RUNs) {
                bioclim_Fut.85_2030_stack = stack(bioclim_Fut.85_2030_stack, 
                                                  readRDS(list.files(
                                                    paste0("./temp_output/", especie, "/"),
                                                    paste0('bioclim_Fut.85_2030',"PA_",PA,"RUN_",RUN,".rds"),
                                                    full.names = T
                                                  )))
              }
            }
            
            names(bioclim_Fut.85_2030_stack)<- paste0("bioclim_Fut_2030_85.", 
                                                      1:nlayers(bioclim_Fut.85_2030_stack))
            
            bioclim_Fut.45_2050_stack = stack()
            for(PA in 1:PAs){
              for (RUN in 1:RUNs) {
                bioclim_Fut.45_2050_stack = stack(bioclim_Fut.45_2050_stack, 
                                                  readRDS(list.files(
                                                    paste0("./temp_output/", 
                                                           especie, "/"),
                                                    paste0('bioclim_Fut.45_2050',"PA_",
                                                           PA,"RUN_",RUN,".rds"),
                                                    full.names = T
                                                  )))
              }
            }
            
            names(bioclim_Fut.45_2050_stack)<- paste0("bioclim_Fut_2050_45.", 
                                                      1:nlayers(bioclim_Fut.45_2050_stack))
            
            bioclim_Fut.85_2050_stack = stack()
            for(PA in 1:PAs){
              for (RUN in 1:RUNs) {
                bioclim_Fut.85_2050_stack = stack(bioclim_Fut.85_2050_stack, 
                                                  readRDS(list.files(
                                                    paste0("./temp_output/", 
                                                           especie, "/"),
                                                    paste0('bioclim_Fut.85_2050',"PA_",
                                                           PA,"RUN_",RUN,".rds"),
                                                    full.names = T
                                                  )))
              }
            }
            
            names(bioclim_Fut.85_2050_stack)<- paste0("bioclim_Fut_2050_85.", 
                                                      1:nlayers(bioclim_Fut.85_2050_stack))
            
            bioclim_Fut.45_2070_stack = stack()
            for(PA in 1:PAs){
              for (RUN in 1:RUNs) {
                bioclim_Fut.45_2070_stack = stack(bioclim_Fut.45_2070_stack, 
                                                  readRDS(list.files(
                                                    paste0("./temp_output/", 
                                                           especie, "/"),
                                                    paste0('bioclim_Fut.45_2070',"PA_",
                                                           PA,"RUN_",RUN,".rds"),
                                                    full.names = T
                                                  )))
              }
            }
            
            names(bioclim_Fut.45_2070_stack)<- paste0("bioclim_Fut_2070_45.", 
                                                      1:nlayers(bioclim_Fut.45_2070_stack))
            
            bioclim_Fut.85_2070_stack = stack()
            for(PA in 1:PAs){
              for (RUN in 1:RUNs) {
                bioclim_Fut.85_2070_stack = stack(bioclim_Fut.85_2070_stack, 
                                                  readRDS(list.files(
                                                    paste0("./temp_output/", 
                                                           especie, "/"),
                                                    paste0('bioclim_Fut.85_2070',"PA_",
                                                           PA,"RUN_",RUN,".rds"),
                                                    full.names = T
                                                  )))
              }
            }
            
            names(bioclim_Fut.85_2070_stack)<- paste0("bioclim_Fut_2070_85.", 
                                                      1:nlayers(bioclim_Fut.85_2070_stack))
            
            bioclim_Fut.45_2090_stack = stack()
            for(PA in 1:PAs){
              for (RUN in 1:RUNs) {
                bioclim_Fut.45_2090_stack = stack(bioclim_Fut.45_2090_stack, 
                                                  readRDS(list.files(
                                                    paste0("./temp_output/", 
                                                           especie, "/"),
                                                    paste0('bioclim_Fut.45_2090',"PA_",
                                                           PA,"RUN_",RUN,".rds"),
                                                    full.names = T
                                                  )))
              }
            }
            
            names(bioclim_Fut.45_2090_stack)<- paste0("bioclim_Fut_2090_45.", 
                                                      1:nlayers(bioclim_Fut.45_2090_stack))
            
            bioclim_Fut.85_2090_stack = stack()
            for(PA in 1:PAs){
              for (RUN in 1:RUNs) {
                bioclim_Fut.85_2090_stack = stack(bioclim_Fut.85_2090_stack, 
                                                  readRDS(list.files(
                                                    paste0("./temp_output/", 
                                                           especie, "/"),
                                                    paste0('bioclim_Fut.85_2090',"PA_",
                                                           PA,"RUN_",RUN,".rds"),
                                                    full.names = T
                                                  )))
              }
            }
            
            names(bioclim_Fut.85_2090_stack)<- paste0("bioclim_Fut_2090_85.", 
                                                      1:nlayers(bioclim_Fut.85_2090_stack))
            
            bioclim.std = rescx(stack(bioclim_Cur_stack, 
                                      bioclim_Fut.45_2030_stack,
                                      bioclim_Fut.85_2030_stack,
                                      bioclim_Fut.45_2050_stack,
                                      bioclim_Fut.85_2050_stack,
                                      bioclim_Fut.45_2070_stack,
                                      bioclim_Fut.85_2070_stack,
                                      bioclim_Fut.45_2090_stack,
                                      bioclim_Fut.85_2090_stack))
            rm(bioclim_Cur_stack, 
               bioclim_Fut.45_2030_stack,
               bioclim_Fut.85_2030_stack,
               bioclim_Fut.45_2050_stack,
               bioclim_Fut.85_2050_stack,
               bioclim_Fut.45_2070_stack,
               bioclim_Fut.85_2070_stack,
               bioclim_Fut.45_2090_stack,
               bioclim_Fut.85_2090_stack)
            gc(reset = T, full = T)
            saveRDS(bioclim.std, paste0("./temp_output/", especie,
                                        "/", especie,"_bioclim.std.rds"))
            

            # enfa #
            enfa_Cur_stack = stack()
            for(PA in 1:PAs){
              for (RUN in 1:RUNs) {
                enfa_Cur_stack = stack(enfa_Cur_stack, readRDS(list.files(
                  paste0("./temp_output/", especie, "/"),
                  paste0('enfa_Cur',"PA_",PA,"RUN_",RUN,".rds"),
                  full.names = T
                )))
              }
            }
            
            names(enfa_Cur_stack)<- paste0("enfa.Cur.", 
                                           1:nlayers(enfa_Cur_stack))
            
            enfa_Fut.45_2030_stack = stack()
            for(PA in 1:PAs){
              for (RUN in 1:RUNs) {
                enfa_Fut.45_2030_stack = stack(enfa_Fut.45_2030_stack, 
                                               readRDS(list.files(
                                                 paste0("./temp_output/", especie, "/"),
                                                 paste0('enfa_Fut.45_2030',"PA_",PA,"RUN_",RUN,".rds"),
                                                 full.names = T
                                               )))
              }
            }
            
            names(enfa_Fut.45_2030_stack)<- paste0("enfa_Fut_2030_45.", 
                                                   1:nlayers(enfa_Fut.45_2030_stack))
            
            enfa_Fut.85_2030_stack = stack()
            for(PA in 1:PAs){
              for (RUN in 1:RUNs) {
                enfa_Fut.85_2030_stack = stack(enfa_Fut.85_2030_stack, 
                                               readRDS(list.files(
                                                 paste0("./temp_output/", 
                                                        especie, "/"),
                                                 paste0('enfa_Fut.85_2030',"PA_",
                                                        PA,"RUN_",RUN,".rds"),
                                                 full.names = T
                                               )))
              }
            }
            
            names(enfa_Fut.85_2030_stack)<- paste0("enfa_Fut_2030_85.", 
                                                   1:nlayers(enfa_Fut.85_2030_stack))
            
            enfa_Fut.45_2050_stack = stack()
            for(PA in 1:PAs){
              for (RUN in 1:RUNs) {
                enfa_Fut.45_2050_stack = stack(enfa_Fut.45_2050_stack, 
                                               readRDS(list.files(
                                                 paste0("./temp_output/", 
                                                        especie, "/"),
                                                 paste0('enfa_Fut.45_2050',"PA_",
                                                        PA,"RUN_",RUN,".rds"),
                                                 full.names = T
                                               )))
              }
            }
            
            names(enfa_Fut.45_2050_stack)<- paste0("enfa_Fut_2050_45.", 
                                                   1:nlayers(enfa_Fut.45_2050_stack))
            
            enfa_Fut.85_2050_stack = stack()
            for(PA in 1:PAs){
              for (RUN in 1:RUNs) {
                enfa_Fut.85_2050_stack = stack(enfa_Fut.85_2050_stack, 
                                               readRDS(list.files(
                                                 paste0("./temp_output/", 
                                                        especie, "/"),
                                                 paste0('enfa_Fut.85_2050',"PA_",
                                                        PA,"RUN_",RUN,".rds"),
                                                 full.names = T
                                               )))
              }
            }
            
            names(enfa_Fut.85_2050_stack)<- paste0("enfa_Fut_2050_85.", 
                                                   1:nlayers(enfa_Fut.85_2050_stack))
            
            enfa_Fut.45_2070_stack = stack()
            for(PA in 1:PAs){
              for (RUN in 1:RUNs) {
                enfa_Fut.45_2070_stack = stack(enfa_Fut.45_2070_stack, 
                                               readRDS(list.files(
                                                 paste0("./temp_output/", 
                                                        especie, "/"),
                                                 paste0('enfa_Fut.45_2070',"PA_",
                                                        PA,"RUN_",RUN,".rds"),
                                                 full.names = T
                                               )))
              }
            }
            
            names(enfa_Fut.45_2070_stack)<- paste0("enfa_Fut_2070_45.", 
                                                   1:nlayers(enfa_Fut.45_2070_stack))
            
            enfa_Fut.85_2070_stack = stack()
            for(PA in 1:PAs){
              for (RUN in 1:RUNs) {
                enfa_Fut.85_2070_stack = stack(enfa_Fut.85_2070_stack, 
                                               readRDS(list.files(
                                                 paste0("./temp_output/", 
                                                        especie, "/"),
                                                 paste0('enfa_Fut.85_2070',"PA_",
                                                        PA,"RUN_",RUN,".rds"),
                                                 full.names = T
                                               )))
              }
            }
            
            names(enfa_Fut.85_2070_stack)<- paste0("enfa_Fut_2070_85.", 
                                                   1:nlayers(enfa_Fut.85_2070_stack))
            
            enfa_Fut.45_2090_stack = stack()
            for(PA in 1:PAs){
              for (RUN in 1:RUNs) {
                enfa_Fut.45_2090_stack = stack(enfa_Fut.45_2090_stack, 
                                               readRDS(list.files(
                                                 paste0("./temp_output/", 
                                                        especie, "/"),
                                                 paste0('enfa_Fut.45_2090',"PA_",
                                                        PA,"RUN_",RUN,".rds"),
                                                 full.names = T
                                               )))
              }
            }
            
            names(enfa_Fut.45_2090_stack)<- paste0("enfa_Fut_2090_45.", 
                                                   1:nlayers(enfa_Fut.45_2090_stack))
            
            enfa_Fut.85_2090_stack = stack()
            for(PA in 1:PAs){
              for (RUN in 1:RUNs) {
                enfa_Fut.85_2090_stack = stack(enfa_Fut.85_2090_stack, 
                                               readRDS(list.files(
                                                 paste0("./temp_output/", 
                                                        especie, "/"),
                                                 paste0('enfa_Fut.85_2090',"PA_",
                                                        PA,"RUN_",RUN,".rds"),
                                                 full.names = T
                                               )))
              }
            }
            
            names(enfa_Fut.85_2090_stack)<- paste0("enfa_Fut_2090_85.", 
                                                   1:nlayers(enfa_Fut.85_2090_stack))
            
            enfa.std = rescx(stack(enfa_Cur_stack,
                                   enfa_Fut.45_2030_stack,
                                   enfa_Fut.85_2030_stack,
                                   enfa_Fut.45_2050_stack,
                                   enfa_Fut.85_2050_stack,
                                   enfa_Fut.45_2070_stack,
                                   enfa_Fut.85_2070_stack,
                                   enfa_Fut.45_2090_stack,
                                   enfa_Fut.85_2090_stack))
            rm(enfa_Cur_stack,
               enfa_Fut.45_2030_stack,
               enfa_Fut.85_2030_stack,
               enfa_Fut.45_2050_stack,
               enfa_Fut.85_2050_stack,
               enfa_Fut.45_2070_stack,
               enfa_Fut.85_2070_stack,
               enfa_Fut.45_2090_stack,
               enfa_Fut.85_2090_stack)
            gc(reset = T, full = T)
            saveRDS(enfa.std, paste0("./temp_output/", especie,
                                     "/", especie,"_enfa.std.rds"))
            
            # glm #
            glm_Cur_stack = stack()
            for(PA in 1:PAs){
              for (RUN in 1:RUNs) {
                glm_Cur_stack = stack(glm_Cur_stack, readRDS(list.files(
                  paste0("./temp_output/", especie, "/"),
                  paste0('glm_Cur',"PA_",PA,"RUN_",RUN,".rds"),
                  full.names = T
                )))
              }
            }
            
            names(glm_Cur_stack)<- paste0("glm.Cur.", 
                                          1:nlayers(glm_Cur_stack))
            
            glm_Fut.45_2030_stack = stack()
            for(PA in 1:PAs){
              for (RUN in 1:RUNs) {
                glm_Fut.45_2030_stack = stack(glm_Fut.45_2030_stack, 
                                              readRDS(list.files(
                                                paste0("./temp_output/", especie, "/"),
                                                paste0('glm_Fut.45_2030',"PA_",PA,"RUN_",RUN,".rds"),
                                                full.names = T
                                              )))
              }
            }
            
            names(glm_Fut.45_2030_stack)<- paste0("glm_Fut_2030_45.", 
                                                  1:nlayers(glm_Fut.45_2030_stack))
            
            glm_Fut.85_2030_stack = stack()
            for(PA in 1:PAs){
              for (RUN in 1:RUNs) {
                glm_Fut.85_2030_stack = stack(glm_Fut.85_2030_stack, 
                                              readRDS(list.files(
                                                paste0("./temp_output/", especie, "/"),
                                                paste0('glm_Fut.85_2030',"PA_",PA,"RUN_",RUN,".rds"),
                                                full.names = T
                                              )))
              }
            }
            
            names(glm_Fut.85_2030_stack)<- paste0("glm_Fut_2030_85.", 
                                                  1:nlayers(glm_Fut.85_2030_stack))
            
            glm_Fut.45_2050_stack = stack()
            for(PA in 1:PAs){
              for (RUN in 1:RUNs) {
                glm_Fut.45_2050_stack = stack(glm_Fut.45_2050_stack, 
                                              readRDS(list.files(
                                                paste0("./temp_output/", 
                                                       especie, "/"),
                                                paste0('glm_Fut.45_2050',"PA_",
                                                       PA,"RUN_",RUN,".rds"),
                                                full.names = T
                                              )))
              }
            }
            
            names(glm_Fut.45_2050_stack)<- paste0("glm_Fut_2050_45.", 
                                                  1:nlayers(glm_Fut.45_2050_stack))
            
            glm_Fut.85_2050_stack = stack()
            for(PA in 1:PAs){
              for (RUN in 1:RUNs) {
                glm_Fut.85_2050_stack = stack(glm_Fut.85_2050_stack, 
                                              readRDS(list.files(
                                                paste0("./temp_output/", 
                                                       especie, "/"),
                                                paste0('glm_Fut.85_2050',"PA_",
                                                       PA,"RUN_",RUN,".rds"),
                                                full.names = T
                                              )))
              }
            }
            
            names(glm_Fut.85_2050_stack)<- paste0("glm_Fut_2050_85.", 
                                                  1:nlayers(glm_Fut.85_2050_stack))
            
            glm_Fut.45_2070_stack = stack()
            for(PA in 1:PAs){
              for (RUN in 1:RUNs) {
                glm_Fut.45_2070_stack = stack(glm_Fut.45_2070_stack, 
                                              readRDS(list.files(
                                                paste0("./temp_output/", 
                                                       especie, "/"),
                                                paste0('glm_Fut.45_2070',"PA_",
                                                       PA,"RUN_",RUN,".rds"),
                                                full.names = T
                                              )))
              }
            }
            
            names(glm_Fut.45_2070_stack)<- paste0("glm_Fut_2070_45.", 
                                                  1:nlayers(glm_Fut.45_2070_stack))
            
            glm_Fut.85_2070_stack = stack()
            for(PA in 1:PAs){
              for (RUN in 1:RUNs) {
                glm_Fut.85_2070_stack = stack(glm_Fut.85_2070_stack, 
                                              readRDS(list.files(
                                                paste0("./temp_output/", 
                                                       especie, "/"),
                                                paste0('glm_Fut.85_2070',"PA_",
                                                       PA,"RUN_",RUN,".rds"),
                                                full.names = T
                                              )))
              }
            }
            
            names(glm_Fut.85_2070_stack)<- paste0("glm_Fut_2070_85.", 
                                                  1:nlayers(glm_Fut.85_2070_stack))
            
            glm_Fut.45_2090_stack = stack()
            for(PA in 1:PAs){
              for (RUN in 1:RUNs) {
                glm_Fut.45_2090_stack = stack(glm_Fut.45_2090_stack, 
                                              readRDS(list.files(
                                                paste0("./temp_output/", 
                                                       especie, "/"),
                                                paste0('glm_Fut.45_2090',"PA_",
                                                       PA,"RUN_",RUN,".rds"),
                                                full.names = T
                                              )))
              }
            }
            
            names(glm_Fut.45_2090_stack)<- paste0("glm_Fut_2090_45.", 
                                                  1:nlayers(glm_Fut.45_2090_stack))
            
            glm_Fut.85_2090_stack = stack()
            for(PA in 1:PAs){
              for (RUN in 1:RUNs) {
                glm_Fut.85_2090_stack = stack(glm_Fut.85_2090_stack, 
                                              readRDS(list.files(
                                                paste0("./temp_output/", 
                                                       especie, "/"),
                                                paste0('glm_Fut.85_2090',"PA_",
                                                       PA,"RUN_",RUN,".rds"),
                                                full.names = T
                                              )))
              }
            }
            
            names(glm_Fut.85_2090_stack)<- paste0("glm_Fut_2090_85.", 
                                                  1:nlayers(glm_Fut.85_2090_stack))
            
            glm.std = rescx(stack(glm_Cur_stack, 
                                  glm_Fut.45_2030_stack,
                                  glm_Fut.85_2030_stack,
                                  glm_Fut.45_2050_stack,
                                  glm_Fut.85_2050_stack,
                                  glm_Fut.45_2070_stack,
                                  glm_Fut.85_2070_stack,
                                  glm_Fut.45_2090_stack,
                                  glm_Fut.85_2090_stack))
            
            rm(glm_Cur_stack, 
               glm_Fut.45_2030_stack,
               glm_Fut.85_2030_stack,
               glm_Fut.45_2050_stack,
               glm_Fut.85_2050_stack,
               glm_Fut.45_2070_stack,
               glm_Fut.85_2070_stack,
               glm_Fut.45_2090_stack,
               glm_Fut.85_2090_stack)
            gc(reset = T, full = T)
            saveRDS(glm.std, paste0("./temp_output/", especie,
                                    "/", especie,"_glm.std.rds"))
            
            # gam #
            gam_Cur_stack = stack()
            for(PA in 1:PAs){
              for (RUN in 1:RUNs) {
                gam_Cur_stack = stack(gam_Cur_stack, readRDS(list.files(
                  paste0("./temp_output/", especie, "/"),
                  paste0('gam_Cur',"PA_",PA,"RUN_",RUN,".rds"),
                  full.names = T
                )))
              }
            }
            
            names(gam_Cur_stack)<- paste0("gam.Cur.", 
                                          1:nlayers(gam_Cur_stack))
            
            gam_Fut.45_2030_stack = stack()
            for(PA in 1:PAs){
              for (RUN in 1:RUNs) {
                gam_Fut.45_2030_stack = stack(gam_Fut.45_2030_stack, 
                                              readRDS(list.files(
                                                paste0("./temp_output/", especie, "/"),
                                                paste0('gam_Fut.45_2030',"PA_",PA,"RUN_",RUN,".rds"),
                                                full.names = T
                                              )))
              }
            }
            
            names(gam_Fut.45_2030_stack)<- paste0("gam_Fut_2030_45.", 
                                                  1:nlayers(gam_Fut.45_2030_stack))
            
            gam_Fut.85_2030_stack = stack()
            for(PA in 1:PAs){
              for (RUN in 1:RUNs) {
                gam_Fut.85_2030_stack = stack(gam_Fut.85_2030_stack, 
                                              readRDS(list.files(
                                                paste0("./temp_output/", especie, "/"),
                                                paste0('gam_Fut.85_2030',"PA_",PA,"RUN_",RUN,".rds"),
                                                full.names = T
                                              )))
              }
            }
            
            names(gam_Fut.85_2030_stack)<- paste0("gam_Fut_2030_85.", 
                                                  1:nlayers(gam_Fut.85_2030_stack))
            
            gam_Fut.45_2050_stack = stack()
            for(PA in 1:PAs){
              for (RUN in 1:RUNs) {
                gam_Fut.45_2050_stack = stack(gam_Fut.45_2050_stack, 
                                              readRDS(list.files(
                                                paste0("./temp_output/", 
                                                       especie, "/"),
                                                paste0('gam_Fut.45_2050',"PA_",
                                                       PA,"RUN_",RUN,".rds"),
                                                full.names = T
                                              )))
              }
            }
            
            names(gam_Fut.45_2050_stack)<- paste0("gam_Fut_2050_45.", 
                                                  1:nlayers(gam_Fut.45_2050_stack))
            
            gam_Fut.85_2050_stack = stack()
            for(PA in 1:PAs){
              for (RUN in 1:RUNs) {
                gam_Fut.85_2050_stack = stack(gam_Fut.85_2050_stack, 
                                              readRDS(list.files(
                                                paste0("./temp_output/", 
                                                       especie, "/"),
                                                paste0('gam_Fut.85_2050',"PA_",
                                                       PA,"RUN_",RUN,".rds"),
                                                full.names = T
                                              )))
              }
            }
            
            names(gam_Fut.85_2050_stack)<- paste0("gam_Fut_2050_85.", 
                                                  1:nlayers(gam_Fut.85_2050_stack))
            
            gam_Fut.45_2070_stack = stack()
            for(PA in 1:PAs){
              for (RUN in 1:RUNs) {
                gam_Fut.45_2070_stack = stack(gam_Fut.45_2070_stack, 
                                              readRDS(list.files(
                                                paste0("./temp_output/", 
                                                       especie, "/"),
                                                paste0('gam_Fut.45_2070',"PA_",
                                                       PA,"RUN_",RUN,".rds"),
                                                full.names = T
                                              )))
              }
            }
            
            names(gam_Fut.45_2070_stack)<- paste0("gam_Fut_2070_45.", 
                                                  1:nlayers(gam_Fut.45_2070_stack))
            
            gam_Fut.85_2070_stack = stack()
            for(PA in 1:PAs){
              for (RUN in 1:RUNs) {
                gam_Fut.85_2070_stack = stack(gam_Fut.85_2070_stack, 
                                              readRDS(list.files(
                                                paste0("./temp_output/", 
                                                       especie, "/"),
                                                paste0('gam_Fut.85_2070',"PA_",
                                                       PA,"RUN_",RUN,".rds"),
                                                full.names = T
                                              )))
              }
            }
            
            names(gam_Fut.85_2070_stack)<- paste0("gam_Fut_2070_85.", 
                                                  1:nlayers(gam_Fut.85_2070_stack))
            
            gam_Fut.45_2090_stack = stack()
            for(PA in 1:PAs){
              for (RUN in 1:RUNs) {
                gam_Fut.45_2090_stack = stack(gam_Fut.45_2090_stack, 
                                              readRDS(list.files(
                                                paste0("./temp_output/", 
                                                       especie, "/"),
                                                paste0('gam_Fut.45_2090',"PA_",
                                                       PA,"RUN_",RUN,".rds"),
                                                full.names = T
                                              )))
              }
            }
            
            names(gam_Fut.45_2090_stack)<- paste0("gam_Fut_2090_45.", 
                                                  1:nlayers(gam_Fut.45_2090_stack))
            
            gam_Fut.85_2090_stack = stack()
            for(PA in 1:PAs){
              for (RUN in 1:RUNs) {
                gam_Fut.85_2090_stack = stack(gam_Fut.85_2090_stack, 
                                              readRDS(list.files(
                                                paste0("./temp_output/", 
                                                       especie, "/"),
                                                paste0('gam_Fut.85_2090',"PA_",
                                                       PA,"RUN_",RUN,".rds"),
                                                full.names = T
                                              )))
              }
            }
            
            names(gam_Fut.85_2090_stack)<- paste0("gam_Fut_2090_85.", 
                                                  1:nlayers(gam_Fut.85_2090_stack))
            
            gam.std = rescx(stack(gam_Cur_stack, 
                                  gam_Fut.45_2030_stack,
                                  gam_Fut.85_2030_stack,
                                  gam_Fut.45_2050_stack,
                                  gam_Fut.85_2050_stack,
                                  gam_Fut.45_2070_stack,
                                  gam_Fut.85_2070_stack,
                                  gam_Fut.45_2090_stack,
                                  gam_Fut.85_2090_stack))
            
            rm(gam_Cur_stack, 
               gam_Fut.45_2030_stack,
               gam_Fut.85_2030_stack,
               gam_Fut.45_2050_stack,
               gam_Fut.85_2050_stack,
               gam_Fut.45_2070_stack,
               gam_Fut.85_2070_stack,
               gam_Fut.45_2090_stack,
               gam_Fut.85_2090_stack)
            gc(reset = T, full = T)
            saveRDS(gam.std, paste0("./temp_output/", especie,
                                    "/", especie,"_gam.std.rds"))
            # maxent #
            maxent_Cur_stack = stack()
            for(PA in 1:PAs){
              for (RUN in 1:RUNs) {
                maxent_Cur_stack = stack(maxent_Cur_stack, readRDS(list.files(
                  paste0("./temp_output/", especie, "/"),
                  paste0('maxent_Cur',"PA_",PA,"RUN_",RUN,".rds"),
                  full.names = T
                )))
              }
            }
            
            names(maxent_Cur_stack)<- paste0("maxent.Cur.", 
                                             1:nlayers(maxent_Cur_stack))
            
            maxent_Fut.45_2030_stack = stack()
            for(PA in 1:PAs){
              for (RUN in 1:RUNs) {
                maxent_Fut.45_2030_stack = stack(maxent_Fut.45_2030_stack, 
                                                 readRDS(list.files(
                                                   paste0("./temp_output/", especie, "/"),
                                                   paste0('maxent_Fut.45_2030',"PA_",PA,"RUN_",RUN,".rds"),
                                                   full.names = T
                                                 )))
              }
            }
            
            names(maxent_Fut.45_2030_stack)<- paste0("maxent_Fut_2030_45.", 
                                                     1:nlayers(maxent_Fut.45_2030_stack))
            
            maxent_Fut.85_2030_stack = stack()
            for(PA in 1:PAs){
              for (RUN in 1:RUNs) {
                maxent_Fut.85_2030_stack = stack(maxent_Fut.85_2030_stack, 
                                                 readRDS(list.files(
                                                   paste0("./temp_output/", especie, "/"),
                                                   paste0('maxent_Fut.85_2030',"PA_",PA,"RUN_",RUN,".rds"),
                                                   full.names = T
                                                 )))
              }
            }
            
            names(maxent_Fut.85_2030_stack)<- paste0("maxent_Fut_2030_85.", 
                                                     1:nlayers(maxent_Fut.85_2030_stack))
            
            maxent_Fut.45_2050_stack = stack()
            for(PA in 1:PAs){
              for (RUN in 1:RUNs) {
                maxent_Fut.45_2050_stack = stack(maxent_Fut.45_2050_stack, 
                                                 readRDS(list.files(
                                                   paste0("./temp_output/", 
                                                          especie, "/"),
                                                   paste0('maxent_Fut.45_2050',"PA_",
                                                          PA,"RUN_",RUN,".rds"),
                                                   full.names = T
                                                 )))
              }
            }
            
            names(maxent_Fut.45_2050_stack)<- paste0("maxent_Fut_2050_45.", 
                                                     1:nlayers(maxent_Fut.45_2050_stack))
            
            maxent_Fut.85_2050_stack = stack()
            for(PA in 1:PAs){
              for (RUN in 1:RUNs) {
                maxent_Fut.85_2050_stack = stack(maxent_Fut.85_2050_stack, 
                                                 readRDS(list.files(
                                                   paste0("./temp_output/", 
                                                          especie, "/"),
                                                   paste0('maxent_Fut.85_2050',"PA_",
                                                          PA,"RUN_",RUN,".rds"),
                                                   full.names = T
                                                 )))
              }
            }
            
            names(maxent_Fut.85_2050_stack)<- paste0("maxent_Fut_2050_85.", 
                                                     1:nlayers(maxent_Fut.85_2050_stack))
            
            maxent_Fut.45_2070_stack = stack()
            for(PA in 1:PAs){
              for (RUN in 1:RUNs) {
                maxent_Fut.45_2070_stack = stack(maxent_Fut.45_2070_stack, 
                                                 readRDS(list.files(
                                                   paste0("./temp_output/", 
                                                          especie, "/"),
                                                   paste0('maxent_Fut.45_2070',"PA_",
                                                          PA,"RUN_",RUN,".rds"),
                                                   full.names = T
                                                 )))
              }
            }
            
            names(maxent_Fut.45_2070_stack)<- paste0("maxent_Fut_2070_45.", 
                                                     1:nlayers(maxent_Fut.45_2070_stack))
            
            maxent_Fut.85_2070_stack = stack()
            for(PA in 1:PAs){
              for (RUN in 1:RUNs) {
                maxent_Fut.85_2070_stack = stack(maxent_Fut.85_2070_stack, 
                                                 readRDS(list.files(
                                                   paste0("./temp_output/", 
                                                          especie, "/"),
                                                   paste0('maxent_Fut.85_2070',"PA_",
                                                          PA,"RUN_",RUN,".rds"),
                                                   full.names = T
                                                 )))
              }
            }
            
            names(maxent_Fut.85_2070_stack)<- paste0("maxent_Fut_2070_85.", 
                                                     1:nlayers(maxent_Fut.85_2070_stack))
            
            maxent_Fut.45_2090_stack = stack()
            for(PA in 1:PAs){
              for (RUN in 1:RUNs) {
                maxent_Fut.45_2090_stack = stack(maxent_Fut.45_2090_stack, 
                                                 readRDS(list.files(
                                                   paste0("./temp_output/", 
                                                          especie, "/"),
                                                   paste0('maxent_Fut.45_2090',"PA_",
                                                          PA,"RUN_",RUN,".rds"),
                                                   full.names = T
                                                 )))
              }
            }
            
            names(maxent_Fut.45_2090_stack)<- paste0("maxent_Fut_2090_45.", 
                                                     1:nlayers(maxent_Fut.45_2090_stack))
            
            maxent_Fut.85_2090_stack = stack()
            for(PA in 1:PAs){
              for (RUN in 1:RUNs) {
                maxent_Fut.85_2090_stack = stack(maxent_Fut.85_2090_stack, 
                                                 readRDS(list.files(
                                                   paste0("./temp_output/", 
                                                          especie, "/"),
                                                   paste0('maxent_Fut.85_2090',"PA_",
                                                          PA,"RUN_",RUN,".rds"),
                                                   full.names = T
                                                 )))
              }
            }
            
            names(maxent_Fut.85_2090_stack)<- paste0("maxent_Fut_2090_85.", 
                                                     1:nlayers(maxent_Fut.85_2090_stack))
            
            maxent.std = rescx(stack(maxent_Cur_stack, 
                                     maxent_Fut.45_2030_stack,
                                     maxent_Fut.85_2030_stack,
                                     maxent_Fut.45_2050_stack,
                                     maxent_Fut.85_2050_stack,
                                     maxent_Fut.45_2070_stack,
                                     maxent_Fut.85_2070_stack,
                                     maxent_Fut.45_2090_stack,
                                     maxent_Fut.85_2090_stack))
            
            rm(maxent_Cur_stack, 
               maxent_Fut.45_2030_stack,
               maxent_Fut.85_2030_stack,
               maxent_Fut.45_2050_stack,
               maxent_Fut.85_2050_stack,
               maxent_Fut.45_2070_stack,
               maxent_Fut.85_2070_stack,
               maxent_Fut.45_2090_stack,
               maxent_Fut.85_2090_stack)
            gc(reset = T, full = T)
            saveRDS(maxent.std, paste0("./temp_output/", especie,
                                       "/", especie,"_maxent.std.rds"))
            
            # rf #
            rf_Cur_stack = stack()
            for(PA in 1:PAs){
              for (RUN in 1:RUNs) {
                rf_Cur_stack = stack(rf_Cur_stack, readRDS(list.files(
                  paste0("./temp_output/", especie, "/"),
                  paste0('rf_Cur',"PA_",PA,"RUN_",RUN,".rds"),
                  full.names = T
                )))
              }
            }
            
            names(rf_Cur_stack)<- paste0("rf.Cur.", 
                                         1:nlayers(rf_Cur_stack))
            
            rf_Fut.45_2030_stack = stack()
            for(PA in 1:PAs){
              for (RUN in 1:RUNs) {
                rf_Fut.45_2030_stack = stack(rf_Fut.45_2030_stack, 
                                             readRDS(list.files(
                                               paste0("./temp_output/", especie, "/"),
                                               paste0('rf_Fut.45_2030',"PA_",PA,"RUN_",RUN,".rds"),
                                               full.names = T
                                             )))
              }
            }
            
            names(rf_Fut.45_2030_stack)<- paste0("rf_Fut_2030_45.", 
                                                 1:nlayers(rf_Fut.45_2030_stack))
            
            rf_Fut.85_2030_stack = stack()
            for(PA in 1:PAs){
              for (RUN in 1:RUNs) {
                rf_Fut.85_2030_stack = stack(rf_Fut.85_2030_stack, 
                                             readRDS(list.files(
                                               paste0("./temp_output/", especie, "/"),
                                               paste0('rf_Fut.85_2030',"PA_",PA,"RUN_",RUN,".rds"),
                                               full.names = T
                                             )))
              }
            }
            
            names(rf_Fut.85_2030_stack)<- paste0("rf_Fut_2030_85.", 
                                                 1:nlayers(rf_Fut.85_2030_stack))
            
            rf_Fut.45_2050_stack = stack()
            for(PA in 1:PAs){
              for (RUN in 1:RUNs) {
                rf_Fut.45_2050_stack = stack(rf_Fut.45_2050_stack, 
                                             readRDS(list.files(
                                               paste0("./temp_output/", 
                                                      especie, "/"),
                                               paste0('rf_Fut.45_2050',"PA_",
                                                      PA,"RUN_",RUN,".rds"),
                                               full.names = T
                                             )))
              }
            }
            
            names(rf_Fut.45_2050_stack)<- paste0("rf_Fut_2050_45.", 
                                                 1:nlayers(rf_Fut.45_2050_stack))
            
            rf_Fut.85_2050_stack = stack()
            for(PA in 1:PAs){
              for (RUN in 1:RUNs) {
                rf_Fut.85_2050_stack = stack(rf_Fut.85_2050_stack, 
                                             readRDS(list.files(
                                               paste0("./temp_output/", 
                                                      especie, "/"),
                                               paste0('rf_Fut.85_2050',"PA_",
                                                      PA,"RUN_",RUN,".rds"),
                                               full.names = T
                                             )))
              }
            }
            
            names(rf_Fut.85_2050_stack)<- paste0("rf_Fut_2050_85.", 
                                                 1:nlayers(rf_Fut.85_2050_stack))
            
            rf_Fut.45_2070_stack = stack()
            for(PA in 1:PAs){
              for (RUN in 1:RUNs) {
                rf_Fut.45_2070_stack = stack(rf_Fut.45_2070_stack, 
                                             readRDS(list.files(
                                               paste0("./temp_output/", 
                                                      especie, "/"),
                                               paste0('rf_Fut.45_2070',"PA_",
                                                      PA,"RUN_",RUN,".rds"),
                                               full.names = T
                                             )))
              }
            }
            
            names(rf_Fut.45_2070_stack)<- paste0("rf_Fut_2070_45.", 
                                                 1:nlayers(rf_Fut.45_2070_stack))
            
            rf_Fut.85_2070_stack = stack()
            for(PA in 1:PAs){
              for (RUN in 1:RUNs) {
                rf_Fut.85_2070_stack = stack(rf_Fut.85_2070_stack, 
                                             readRDS(list.files(
                                               paste0("./temp_output/", 
                                                      especie, "/"),
                                               paste0('rf_Fut.85_2070',"PA_",
                                                      PA,"RUN_",RUN,".rds"),
                                               full.names = T
                                             )))
              }
            }
            
            names(rf_Fut.85_2070_stack)<- paste0("rf_Fut_2070_85.", 
                                                 1:nlayers(rf_Fut.85_2070_stack))
            
            rf_Fut.45_2090_stack = stack()
            for(PA in 1:PAs){
              for (RUN in 1:RUNs) {
                rf_Fut.45_2090_stack = stack(rf_Fut.45_2090_stack, 
                                             readRDS(list.files(
                                               paste0("./temp_output/", 
                                                      especie, "/"),
                                               paste0('rf_Fut.45_2090',"PA_",
                                                      PA,"RUN_",RUN,".rds"),
                                               full.names = T
                                             )))
              }
            }
            
            names(rf_Fut.45_2090_stack)<- paste0("rf_Fut_2090_45.", 
                                                 1:nlayers(rf_Fut.45_2090_stack))
            
            rf_Fut.85_2090_stack = stack()
            for(PA in 1:PAs){
              for (RUN in 1:RUNs) {
                rf_Fut.85_2090_stack = stack(rf_Fut.85_2090_stack, 
                                             readRDS(list.files(
                                               paste0("./temp_output/", 
                                                      especie, "/"),
                                               paste0('rf_Fut.85_2090',"PA_",
                                                      PA,"RUN_",RUN,".rds"),
                                               full.names = T
                                             )))
              }
            }
            
            names(rf_Fut.85_2090_stack)<- paste0("rf_Fut_2090_85.", 
                                                 1:nlayers(rf_Fut.85_2090_stack))
            
            rf.std = rescx(stack(rf_Cur_stack, 
                                 rf_Fut.45_2030_stack,
                                 rf_Fut.85_2030_stack,
                                 rf_Fut.45_2050_stack,
                                 rf_Fut.85_2050_stack,
                                 rf_Fut.45_2070_stack,
                                 rf_Fut.85_2070_stack,
                                 rf_Fut.45_2090_stack,
                                 rf_Fut.85_2090_stack))
            
            rm(rf_Cur_stack, 
               rf_Fut.45_2030_stack,
               rf_Fut.85_2030_stack,
               rf_Fut.45_2050_stack,
               rf_Fut.85_2050_stack,
               rf_Fut.45_2070_stack,
               rf_Fut.85_2070_stack,
               rf_Fut.45_2090_stack,
               rf_Fut.85_2090_stack)
            gc(reset = T, full = T)
            saveRDS(rf.std, paste0("./temp_output/", especie,
                                   "/", especie,"_rf.std.rds"))
          }
          # Evaluating ####
          
          for (RUN in seq(RUNs)) {
            
            ##### bioclim ###
            bioclim_eval <-
              eval.All.Model(bioclim.std[[RUN]], 
                             test.pa_list[[RUN]])
            
            bioclim_th.spec_sens <-
              dismo::threshold(bioclim_eval, "spec_sens")
            
            bioclim_eval.spec_sens <-
              eval.All.Model(bioclim.std[[RUN]], test.pa_list[[RUN]],
                             tr = bioclim_th.spec_sens)
            
            bioclim_e.spec_sens <- c(
              AUC.SpecSens = bioclim_eval.spec_sens@auc,
              TPR.SpecSens = bioclim_eval.spec_sens@TPR,
              TNR.SpecSens = bioclim_eval.spec_sens@TNR,
              thr.SpecSens = bioclim_th.spec_sens,
              TSS.SpecSens = (bioclim_eval.spec_sens@TPR+bioclim_eval.spec_sens@TNR)-1
           )
            
            bioclim.e = rbind(bioclim.e, c(bioclim_e.spec_sens, 
                                           PA = PA, RUN = RUN))
            rownames(bioclim.e) = rep(paste0("bioclim"), 
                                      nrow(bioclim.e))
            
            write.csv(bioclim.e,
                      paste0("./temp_output/",especie,"/",
                             "bioclim_eval.all.csv"),
                      row.names = T)
            
            #### enfa ###
            enfa_eval <-
              eval.All.Model(enfa.std[[RUN]],
                             test.pa_list[[RUN]])
            
            enfa_th.spec_sens <-
              dismo::threshold(enfa_eval, "spec_sens")
            
            enfa_eval.spec_sens <-
              eval.All.Model(enfa.std[[RUN]], test.pa_list[[RUN]],
                             tr = enfa_th.spec_sens)
            
            enfa_e.spec_sens <- c(
              AUC.SpecSens = enfa_eval.spec_sens@auc,
              TPR.SpecSens = enfa_eval.spec_sens@TPR,
              TNR.SpecSens = enfa_eval.spec_sens@TNR,
              thr.SpecSens = enfa_th.spec_sens,
              TSS.SpecSens = (enfa_eval.spec_sens@TPR+enfa_eval.spec_sens@TNR)-1
            )
            
            enfa.e = rbind(enfa.e, c(enfa_e.spec_sens,
                                     PA = PA, RUN = RUN))
            rownames(enfa.e) = rep(paste0("enfa"),
                                   nrow(enfa.e))
            
            write.csv(enfa.e,
                      paste0("./temp_output/",especie,"/",
                             "enfa_eval.all.csv"),
                      row.names = T)
            
            #### glm ### 
            glm_eval <-
              eval.All.Model(glm.std[[RUN]], 
                             test.pa_list[[RUN]])
            
            glm_th.spec_sens <-
              dismo::threshold(glm_eval, "spec_sens")
            
            glm_eval.spec_sens <-
              eval.All.Model(glm.std[[RUN]], test.pa_list[[RUN]],
                             tr = glm_th.spec_sens)
            
            glm_e.spec_sens <- c(
              AUC.SpecSens = glm_eval.spec_sens@auc,
              TPR.SpecSens = glm_eval.spec_sens@TPR,
              TNR.SpecSens = glm_eval.spec_sens@TNR,
              thr.SpecSens = glm_th.spec_sens,
             TSS.SpecSens = (glm_eval.spec_sens@TPR+glm_eval.spec_sens@TNR)-1
           )
            
            glm.e = rbind(glm.e, c(glm_e.spec_sens,
                                   PA = PA, RUN = RUN))
            rownames(glm.e) = rep(paste0("glm"), 
                                  nrow(glm.e))
            
            write.csv(glm.e,
                      paste0("./temp_output/",especie,"/",
                             "glm_eval.all.csv"),
                      row.names = T)
            
            #### gam ### 
            gam_eval <-
              eval.All.Model(gam.std[[RUN]], 
                             test.pa_list[[RUN]])
            
            gam_th.spec_sens <-
              dismo::threshold(gam_eval, "spec_sens")
            
            gam_eval.spec_sens <-
              eval.All.Model(gam.std[[RUN]], test.pa_list[[RUN]],
                             tr = gam_th.spec_sens)
            
            gam_e.spec_sens <- c(
              AUC.SpecSens = gam_eval.spec_sens@auc,
              TPR.SpecSens = gam_eval.spec_sens@TPR,
              TNR.SpecSens = gam_eval.spec_sens@TNR,
              thr.SpecSens = gam_th.spec_sens,
           TSS.SpecSens = (gam_eval.spec_sens@TPR+gam_eval.spec_sens@TNR)-1
           )
            
            gam.e = rbind(gam.e, c(gam_e.spec_sens,
                                   PA = PA, RUN = RUN))
            rownames(gam.e) = rep(paste0("gam"), 
                                  nrow(gam.e))
            
            write.csv(gam.e,
                      paste0("./temp_output/",especie,"/",
                             "gam_eval.all.csv"),
                      row.names = T)
            
            #### maxent ### 
            maxent_eval <-
              eval.All.Model(maxent.std[[RUN]], 
                             test.pa_list[[RUN]])
            
            maxent_th.spec_sens <-
              dismo::threshold(maxent_eval, "spec_sens")
            
            maxent_eval.spec_sens <-
              eval.All.Model(maxent.std[[RUN]], test.pa_list[[RUN]],
                             tr = maxent_th.spec_sens)
            
            maxent_e.spec_sens <- c(
              AUC.SpecSens = maxent_eval.spec_sens@auc,
              TPR.SpecSens = maxent_eval.spec_sens@TPR,
              TNR.SpecSens = maxent_eval.spec_sens@TNR,
              thr.SpecSens = maxent_th.spec_sens,
           TSS.SpecSens = (maxent_eval.spec_sens@TPR+maxent_eval.spec_sens@TNR)-1
          )
            
            maxent.e = rbind(maxent.e, c(maxent_e.spec_sens,
                                         PA = PA, RUN = RUN))
            rownames(maxent.e) = rep(paste0("maxent"), 
                                     nrow(maxent.e))
            
            write.csv(maxent.e,
                      paste0("./temp_output/",especie,"/",
                             "maxent_eval.all.csv"),
                      row.names = T)
            
            
            #### rf ###
            rf_eval <-
              eval.All.Model(rf.std[[RUN]], 
                             test.pa_list[[RUN]])
            
            rf_th.spec_sens <-
              dismo::threshold(rf_eval, "spec_sens")
            
            rf_eval.spec_sens <-
              eval.All.Model(rf.std[[RUN]], test.pa_list[[RUN]],
                             tr = rf_th.spec_sens)
            
            rf_e.spec_sens <- c(
              AUC.SpecSens = rf_eval.spec_sens@auc,
              TPR.SpecSens = rf_eval.spec_sens@TPR,
              TNR.SpecSens = rf_eval.spec_sens@TNR,
              thr.SpecSens = rf_th.spec_sens,
            TSS.SpecSens = (rf_eval.spec_sens@TPR+rf_eval.spec_sens@TNR)-1
          )
            
            rf.e = rbind(rf.e, c(rf_e.spec_sens, PA = PA, RUN = RUN))
            rownames(rf.e) = rep(paste0("rf"), 
                                 nrow(rf.e))
            
            write.csv(rf.e,
                      paste0("./temp_output/",especie,"/", "rf_eval.all.csv"),
                      row.names = T)
            
            
          }
          
          
          remove(list = ls()[c(grep(
            "bioclim_", as.factor(ls())))])
          gc(reset = TRUE, full = T)
         
          remove(list = ls()[c(grep(
            "enfa_", as.factor(ls())))])
          gc(reset = TRUE, full = T)
          
          
          remove(list = ls()[c(grep(
            "glm_", as.factor(ls())))])
          gc(reset = TRUE, full = T)
          
          
          remove(list = ls()[c(grep(
            "gam_", as.factor(ls())))])
          gc(reset = TRUE, full = T)
          
          
          remove(list = ls()[c(grep(
            "maxent_", as.factor(ls())))])
          gc(reset = TRUE, full = T)
          
          
          remove(list = ls()[c(grep(
            "rf_", as.factor(ls())))])
          gc(reset = TRUE, full = T)
          
          # }##### End RUN loop ###
          # print(paste0(especie," ","PA set"," ", PA))
          
          # } ##### End PAs loop ###
          
          # Writing Evaluation Resultes #
          {
            # bioclim #
            bioclim_e = read.table(paste0("./temp_output/",
                                          especie,"/", "bioclim_eval.all.csv"), 
                                   header = T, sep = ',')
            
         
            # enfa #
            enfa_e = read.table(paste0("./temp_output/",
                                       especie,"/", "enfa_eval.all.csv"),
                                header = T, sep = ',')
            
            # glm #
            glm_e = read.table(paste0("./temp_output/",especie,"/",
                                      "glm_eval.all.csv"), 
                               header = T, sep = ',')
            
            # gam #
            gam_e = read.table(paste0("./temp_output/",
                                      especie,"/", "gam_eval.all.csv"), 
                               header = T, sep = ',')
            
            # maxent #
            maxent_e = read.table(paste0("./temp_output/",
                                         especie,"/", "maxent_eval.all.csv"), 
                                  header = T, sep = ',')
            
            
            # rf #
            rf_e = read.table(paste0("./temp_output/",
                                     especie,"/","rf_eval.all.csv"), 
                              header = T, sep = ',')
            
            Evaluation.all <- rbind(bioclim_e,
                                    enfa_e,
                                    glm_e, 
                                    gam_e, 
                                    maxent_e, 
                                    rf_e)
            
            rm(bioclim_e,
               enfa_e,
               glm_e, 
               gam_e, 
               maxent_e, 
               rf_e)  
            
            Evaluation.all <- cbind(Evaluation.all, ID = 1:nrow(Evaluation.all))
            
            colnames(Evaluation.all) <-
              c("Model",
                "AUC",
                "TPR",
                "TNR",
                "threshold",
                "TSS",
                "PA",
                "RUN",
                "ID")
            
            ## Write Evaluation ##
            write.csv(Evaluation.all,
                      paste0("./outputs/", especie, "_", "Evaluation.all.csv"),
                      row.names = T)
            
            
            gc(reset = TRUE, full = TRUE)
          }
          ## Write Rasters ##
          {
            # Current ####
            
            
            # Current Ensemble ####
            
            Current.All = stack(subset(bioclim.std, grep("Cur", 
                                                         names(bioclim.std))),
                                subset(enfa.std, grep("Cur",
                                                      names(enfa.std))),
                                subset(glm.std, grep("Cur", 
                                                     names(glm.std))), 
                                subset(gam.std, grep("Cur", 
                                                     names(gam.std))), 
                                subset(maxent.std, grep("Cur", 
                                                        names(maxent.std))), 
                                subset(rf.std, grep("Cur", 
                                                    names(rf.std))))
            saveRDS(
              Current.All,
              paste0("./outputs/", especie, "_", "Current.All.rds"))
            
            Current.mean <- 
              weighted.mean(Current.All, Evaluation.all[, "TSS"])
            writeRaster(
              Current.mean,
              paste0("./outputs/", especie, "_", "Current.mean.tif"), 
              formtat = "GTiff", overwrite = T)
            
            # Binary Transformation #
            Current.bin <- biomod2::BinaryTransformation(Current.All, 
                                                         Evaluation.all[,5])
            names(Current.bin) = names(Current.All)
            
            gc(reset = T, full = T)
            writeRaster(
              Current.bin,
              paste0("./outputs/", especie, "_", "Current.bin.tif"), 
              formtat = "GTiff", overwrite = T)
            rm(Current.mean, 
               Current.All, 
               Current.bin)
            gc(reset = TRUE, full = TRUE)
          }
          # Future 2030 45 ####
          {
            Fut_2030_45.All = stack(
              subset(bioclim.std, grep("_Fut_2030_45",
                                       names(bioclim.std))),
              subset(enfa.std, grep("_Fut_2030_45",
                                    names(enfa.std))),
              subset(glm.std, grep("_Fut_2030_45",
                                   names(glm.std))),
              subset(gam.std, grep("_Fut_2030_45",
                                   names(gam.std))),
              subset(maxent.std, grep("_Fut_2030_45",
                                      names(maxent.std))),
              subset(rf.std, grep("_Fut_2030_45",
                                  names(rf.std)))
            )
            saveRDS(
              Fut_2030_45.All,
              paste0("./outputs/", especie, "_", "Future.45.2030.All.rds"))
            
            ## Future Ensemble 2030 45 ####
            
            Future.45.2030.mean <- weighted.mean(Fut_2030_45.All,
                                                 Evaluation.all[, "TSS"])
            
            writeRaster(
              Future.45.2030.mean,
              paste0("./outputs/", especie, "_", "Future.45.2030.mean.tif"), 
              formtat = "GTiff", overwrite = T)
            
            # Binary Transformation #
            Future.45.2030.bin <- biomod2::BinaryTransformation(Fut_2030_45.All,
                                                                Evaluation.all[,5])
            names(Future.45.2030.bin) = names(Fut_2030_45.All)
            gc(reset = T)
            writeRaster(
              Future.45.2030.bin,
              paste0("./outputs/", especie, "_", "Future.45.2030.bin.tif"), 
              formtat = "GTiff", overwrite = T)
            rm(Fut_2030_45.All, 
               Future.45.2030.mean, 
               Future.45.2030.bin)
            gc(reset = TRUE)
          }
          
          # Future 2030 85 ####
          {
            Fut_2030_85.All = stack(
              subset(bioclim.std, grep("_Fut_2030_85",
                                       names(bioclim.std))),
              subset(enfa.std, grep("_Fut_2030_85",
                                    names(enfa.std))),
              subset(glm.std, grep("_Fut_2030_85",
                                   names(glm.std))),
              subset(gam.std, grep("_Fut_2030_85",
                                   names(gam.std))),
              subset(maxent.std, grep("_Fut_2030_85",
                                      names(maxent.std))),
              subset(rf.std, grep("_Fut_2030_85",
                                  names(rf.std)))
            )
            saveRDS(
              Fut_2030_85.All,
              paste0("./outputs/", especie, "_", "Future.85.2030.All.rds"))
            ## Future Ensemble 2030 85 ####
            
            Future.85.2030.mean <- weighted.mean(Fut_2030_85.All,
                                                 Evaluation.all[, "TSS"])
            
            writeRaster(
              Future.85.2030.mean,
              paste0("./outputs/", especie, "_", "Future.85.2030.mean.tif"), 
              formtat = "GTiff", overwrite = T)
            
            # Binary Transformation #
            Future.85.2030.bin <- biomod2::BinaryTransformation(
              Fut_2030_85.All, Evaluation.all[,5])
            names(Future.85.2030.bin) = names(Fut_2030_85.All)
            gc(reset = T)
            writeRaster(
              Future.85.2030.bin,
              paste0("./outputs/", especie, "_", "Future.85.2030.bin.tif"), 
              formtat = "GTiff", overwrite = T)
            
            rm(Fut_2030_85.All, 
               Future.85.2030.mean, 
               Future.85.2030.bin)
            gc(reset = TRUE)
          }
          # Future 2050 45 ####
          {
            Fut_2050_45.All = stack(
              subset(bioclim.std, grep("_Fut_2050_45",
                                       names(bioclim.std))),
              subset(enfa.std, grep("_Fut_2050_45",
                                    names(enfa.std))),
              subset(glm.std, grep("_Fut_2050_45",
                                   names(glm.std))),
              subset(gam.std, grep("_Fut_2050_45",
                                   names(gam.std))),
              subset(maxent.std, grep("_Fut_2050_45",
                                      names(maxent.std))),
              subset(rf.std, grep("_Fut_2050_45",
                                  names(rf.std)))
            )
            saveRDS(
              Fut_2050_45.All,
              paste0("./outputs/", especie, "_", "Future.45.2050.All.rds"))
            
            ## Future Ensemble 2050 45 ####
            
            Future.45.2050.mean <- weighted.mean(Fut_2050_45.All,
                                                 Evaluation.all[, "TSS"])
            
            writeRaster(
              Future.45.2050.mean,
              paste0("./outputs/", especie, "_", "Future.45.2050.mean.tif"), 
              formtat = "GTiff", overwrite = T)
            
            # Binary Transformation #
            Future.45.2050.bin <- biomod2::BinaryTransformation(Fut_2050_45.All,
                                                                Evaluation.all[,5])
            names(Future.45.2050.bin) = names(Fut_2050_45.All)
            gc(reset = T)
            writeRaster(
              Future.45.2050.bin,
              paste0("./outputs/", especie, "_", "Future.45.2050.bin.tif"), 
              formtat = "GTiff", overwrite = T)
            rm(Fut_2050_45.All, 
               Future.45.2050.mean, 
               Future.45.2050.bin)
            gc(reset = TRUE)
          }
          
          # Future 2050 85 ####
          {
            Fut_2050_85.All = stack(
              subset(bioclim.std, grep("_Fut_2050_85",
                                       names(bioclim.std))),
              subset(enfa.std, grep("_Fut_2050_85",
                                    names(enfa.std))),
              subset(glm.std, grep("_Fut_2050_85",
                                   names(glm.std))),
              subset(gam.std, grep("_Fut_2050_85",
                                   names(gam.std))),
              subset(maxent.std, grep("_Fut_2050_85",
                                      names(maxent.std))),
              subset(rf.std, grep("_Fut_2050_85",
                                  names(rf.std)))
            )
            saveRDS(
              Fut_2050_85.All,
              paste0("./outputs/", especie, "_", "Future.85.2050.All.rds"))
            ## Future Ensemble 2050 85 ####
            
            Future.85.2050.mean <- weighted.mean(Fut_2050_85.All,
                                                 Evaluation.all[, "TSS"])
            
            writeRaster(
              Future.85.2050.mean,
              paste0("./outputs/", especie, "_", "Future.85.2050.mean.tif"), 
              formtat = "GTiff", overwrite = T)
            
            # Binary Transformation #
            Future.85.2050.bin <- biomod2::BinaryTransformation(
              Fut_2050_85.All, Evaluation.all[,5])
            names(Future.85.2050.bin) = names(Fut_2050_85.All)
            gc(reset = T)
            writeRaster(
              Future.85.2050.bin,
              paste0("./outputs/", especie, "_", "Future.85.2050.bin.tif"), 
              formtat = "GTiff", overwrite = T)
            
            rm(Fut_2050_85.All, 
               Future.85.2050.mean, 
               Future.85.2050.bin)
            gc(reset = TRUE)
          }
          
          # Future 2070 45 ####
          {
            Fut_2070_45.All = stack(
              subset(bioclim.std, grep("_Fut_2070_45",
                                       names(bioclim.std))),
              subset(enfa.std, grep("_Fut_2070_45",
                                    names(enfa.std))),
              subset(glm.std, grep("_Fut_2070_45",
                                   names(glm.std))),
              subset(gam.std, grep("_Fut_2070_45",
                                   names(gam.std))),
              subset(maxent.std, grep("_Fut_2070_45",
                                      names(maxent.std))),
              subset(rf.std, grep("_Fut_2070_45",
                                  names(rf.std)))
            )
            saveRDS(
              Fut_2070_45.All,
              paste0("./outputs/", especie, "_", "Future.45.2070.All.rds"))
            
            ## Future Ensemble 2070 45 ####
            
            Future.45.2070.mean <- weighted.mean(Fut_2070_45.All,
                                                 Evaluation.all[, "TSS"])
            
            writeRaster(
              Future.45.2070.mean,
              paste0("./outputs/", especie, "_", "Future.45.2070.mean.tif"), 
              formtat = "GTiff", overwrite = T)
            
            # Binary Transformation #
            Future.45.2070.bin <- biomod2::BinaryTransformation(Fut_2070_45.All,
                                                                Evaluation.all[,5])
            names(Future.45.2070.bin) = names(Fut_2070_45.All)
            gc(reset = T)
            writeRaster(
              Future.45.2070.bin,
              paste0("./outputs/", especie, "_", "Future.45.2070.bin.tif"), 
              formtat = "GTiff", overwrite = T)
            rm(Fut_2070_45.All, 
               Future.45.2070.mean, 
               Future.45.2070.bin)
            gc(reset = TRUE)
          }
          
          # Future 2070 85 ####
          {
            Fut_2070_85.All = stack(
              subset(bioclim.std, grep("_Fut_2070_85",
                                       names(bioclim.std))),
              subset(enfa.std, grep("_Fut_2070_85",
                                    names(enfa.std))),
              subset(glm.std, grep("_Fut_2070_85",
                                   names(glm.std))),
              subset(gam.std, grep("_Fut_2070_85",
                                   names(gam.std))),
              subset(maxent.std, grep("_Fut_2070_85",
                                      names(maxent.std))),
              subset(rf.std, grep("_Fut_2070_85",
                                  names(rf.std)))
            )
            saveRDS(
              Fut_2070_85.All,
              paste0("./outputs/", especie, "_", "Future.85.2070.All.rds"))
            ## Future Ensemble 2070 85 ####
            
            Future.85.2070.mean <- weighted.mean(Fut_2070_85.All,
                                                 Evaluation.all[, "TSS"])
            
            writeRaster(
              Future.85.2070.mean,
              paste0("./outputs/", especie, "_", "Future.85.2070.mean.tif"), 
              formtat = "GTiff", overwrite = T)
            
            # Binary Transformation #
            Future.85.2070.bin <- biomod2::BinaryTransformation(
              Fut_2070_85.All, Evaluation.all[,5])
            names(Future.85.2070.bin) = names(Fut_2070_85.All)
            gc(reset = T)
            writeRaster(
              Future.85.2070.bin,
              paste0("./outputs/", especie, "_", "Future.85.2070.bin.tif"), 
              formtat = "GTiff", overwrite = T)
            
            rm(Fut_2070_85.All, 
               Future.85.2070.mean, 
               Future.85.2070.bin)
            gc(reset = TRUE)
          }
          
          # Future 2090 45 ####
          {
            Fut_2090_45.All = stack(
              subset(bioclim.std, grep("_Fut_2090_45",
                                       names(bioclim.std))),
              subset(enfa.std, grep("_Fut_2090_45",
                                    names(enfa.std))),
              subset(glm.std, grep("_Fut_2090_45",
                                   names(glm.std))),
              subset(gam.std, grep("_Fut_2090_45",
                                   names(gam.std))),
              subset(maxent.std, grep("_Fut_2090_45",
                                      names(maxent.std))),
              subset(rf.std, grep("_Fut_2090_45",
                                  names(rf.std)))
            )
            saveRDS(
              Fut_2090_45.All,
              paste0("./outputs/", especie, "_", "Future.45.2090.All.rds"))
            
            ## Future Ensemble 2090 45 ####
            
            Future.45.2090.mean <- weighted.mean(Fut_2090_45.All,
                                                 Evaluation.all[, "TSS"])
            
            writeRaster(
              Future.45.2090.mean,
              paste0("./outputs/", especie, "_", "Future.45.2090.mean.tif"), 
              formtat = "GTiff", overwrite = T)
            
            # Binary Transformation #
            Future.45.2090.bin <- biomod2::BinaryTransformation(Fut_2090_45.All,
                                                                Evaluation.all[,5])
            names(Future.45.2090.bin) = names(Fut_2090_45.All)
            gc(reset = T)
            writeRaster(
              Future.45.2090.bin,
              paste0("./outputs/", especie, "_", "Future.45.2090.bin.tif"), 
              formtat = "GTiff", overwrite = T)
            rm(Fut_2090_45.All, 
               Future.45.2090.mean, 
               Future.45.2090.bin)
            gc(reset = TRUE)
          }
          
          # Future 2090 85 ####
          {
            Fut_2090_85.All = stack(
              subset(bioclim.std, grep("_Fut_2090_85",
                                       names(bioclim.std))),
              subset(enfa.std, grep("_Fut_2090_85",
                                    names(enfa.std))),
              subset(glm.std, grep("_Fut_2090_85",
                                   names(glm.std))),
              subset(gam.std, grep("_Fut_2090_85",
                                   names(gam.std))),
              subset(maxent.std, grep("_Fut_2090_85",
                                      names(maxent.std))),
              subset(rf.std, grep("_Fut_2090_85",
                                  names(rf.std)))
            )
            saveRDS(
              Fut_2090_85.All,
              paste0("./outputs/", especie, "_", "Future.85.2090.All.rds"))
            ## Future Ensemble 2090 85 ####
            
            Future.85.2090.mean <- weighted.mean(Fut_2090_85.All,
                                                 Evaluation.all[, "TSS"])
            
            writeRaster(
              Future.85.2090.mean,
              paste0("./outputs/", especie, "_", "Future.85.2090.mean.tif"), 
              formtat = "GTiff", overwrite = T)
            
            # Binary Transformation #
            Future.85.2090.bin <- biomod2::BinaryTransformation(
              Fut_2090_85.All, Evaluation.all[,5])
            names(Future.85.2090.bin) = names(Fut_2090_85.All)
            gc(reset = T)
            writeRaster(
              Future.85.2090.bin,
              paste0("./outputs/", especie, "_", "Future.85.2090.bin.tif"), 
              formtat = "GTiff", overwrite = T)
            
            rm(Fut_2090_85.All, 
               Future.85.2090.mean, 
               Future.85.2090.bin)
            gc(reset = TRUE)
          }
          
          
          remove(list = ls()[c(grep(
            c("Cur"), as.factor(ls())))])
          
          remove(list = ls()[c(grep(
            c("Fut"), as.factor(ls())))])
          
          remove(list = ls()[c(grep(
            "bioclim", as.factor(ls())))])
          gc(reset = TRUE, full = T)
          
          remove(list = ls()[c(grep(
            "enfa", as.factor(ls())))])
          gc(reset = TRUE, full = T)
          
          
          remove(list = ls()[c(grep(
            "glm", as.factor(ls())))])
          gc(reset = TRUE, full = T)
          
          
          remove(list = ls()[c(grep(
            "gam", as.factor(ls())))])
          gc(reset = TRUE, full = T)
          
          
          remove(list = ls()[c(grep(
            "maxent", as.factor(ls())))])
          gc(reset = TRUE, full = T)
          
          
          remove(list = ls()[c(grep(
            "rf", as.factor(ls())))])
          gc(reset = TRUE, full = T)
          
          # unlink(especie,recursive = T, force = T)
          #--------------------#
          # Move the files  ###
          #------------------#
          file.move((list.files("./outputs/", paste0(especie, "_", "Current"),
                                full.names = TRUE)), 
                    (paste0("./outputs/", especie, '.', 'Presente')), 
                    overwrite = TRUE)
          
          filesstrings::file.move((list.files(
            "./outputs/", paste0(especie, "_", "Future.45.2030"),
            full.names = TRUE
          )), (paste0("./outputs/", especie, '.', 'RCP4.5.2030')), 
          overwrite = TRUE)
          
          filesstrings::file.move((list.files(
            "./outputs/",  paste0(especie, "_", "Future.85.2030"),
            full.names = TRUE
          )),(paste0("./outputs/", especie, '.', 'RCP8.5.2030')), 
          overwrite = TRUE)
          
          filesstrings::file.move((list.files(
            "./outputs/", paste0(especie, "_", "Future.45.2050"),
            full.names = TRUE
          )), (paste0("./outputs/", especie, '.', 'RCP4.5.2050')), 
          overwrite = TRUE)
          
          filesstrings::file.move((list.files(
            "./outputs/",  paste0(especie, "_", "Future.85.2050"),
            full.names = TRUE
          )),(paste0("./outputs/", especie, '.', 'RCP8.5.2050')), 
          overwrite = TRUE)
          
          filesstrings::file.move((list.files(
            "./outputs/", paste0(especie, "_", "Future.45.2070"),
            full.names = TRUE
          )), (paste0("./outputs/", especie, '.', 'RCP4.5.2070')), 
          overwrite = TRUE)
          
          filesstrings::file.move((list.files(
            "./outputs/",  paste0(especie, "_", "Future.85.2070"),
            full.names = TRUE
          )),(paste0("./outputs/", especie, '.', 'RCP8.5.2070')), 
          overwrite = TRUE)
          
          filesstrings::file.move((list.files(
            "./outputs/", paste0(especie, "_", "Future.45.2090"),
            full.names = TRUE
          )), (paste0("./outputs/", especie, '.', 'RCP4.5.2090')), 
          overwrite = TRUE)
          
          filesstrings::file.move((list.files(
            "./outputs/",  paste0(especie, "_", "Future.85.2090"),
            full.names = TRUE
          )),(paste0("./outputs/", especie, '.', 'RCP8.5.2090')), 
          overwrite = TRUE)
          
          filesstrings::file.move((list.files(
            "./outputs/", paste0(especie, "_"),
            full.names = TRUE
          )), (paste0("./outputs/", especie)), 
          overwrite = TRUE)
          
          parallel::stopCluster(cl1)
          doParallel::stopImplicitCluster()
          # Time Compute ###
          sink("./outputs/tempo.txt", append = T)
          print(Sys.time() - ini1)
          sink()
          print(paste0("Finishing", " ",especie," " ,"modeling"))
        }## End species loop ####

# beep(sound = 2)
# beep(sound = 2)
# # base::quit(save = "yes")
