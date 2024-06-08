### script ensemble weighted average ###

# Mauricio Humberto Vancine - mauricio.vancine@gmail.com
# 05/06/2017

###----------------------------------------------------------------------------###

# 1. clear memory and load packages 
# clear workspace and increase memory
rm(list = ls())
gc()
memory.limit(size = 1.75e13) 

# packages
if(!require("pacman")) install.packages("pacman")
pacman::p_load(raster, rgdal, data.table)

# verify packages
search()

###----------------------------------------------------------------------------###

# import data
# directory
setwd("C:\\Users\\Thadeu Sobral\\Documents\\MEGA\\_disciplina_enm_lavras_2018\\scripts_r\\enm\\01_dados\\modelos\\results")

# enms
# list files
tif <- list.files(patt = ".tif$")
tif

enm <- raster(tif[3])
enm
plot(enm)

# evaluate
txt <- list.files(patt = ".txt$")
txt

eva <- lapply(txt, read.table)
eva
names(eva) <- txt
eva

###-----------------------------------------------------------------------------###

## weighted average ensemble 
# lists
# species
sp <- ("B.balansae")
sp

# periods
pe <- c("0k", "6k", "21k")
pe


# data.table
da <- data.table()
da

# raster
ens <- enm[[1]]
ens[] <- NA
names(ens) <- "ens"
ens

# ensemble
dir.create("ensemble_wei")

for(i in sp){
  tif.sp <- grep(i, tif, value = T)
  eva.sp <- eva[grep(i, names(eva))]
  
  tss <- do.call("rbind", eva.sp)$TSS
  id.tss <- which(tss > .5)
  tss.05 <- tss[tss > .5]
  
  for(j in pe){
    tif.pe <- grep(j, tif.sp, value = T)
    da <- rbind(da, stack(tif.pe[id.tss])[], use.names = F)}
  
  da.r <- data.table(decostand(da, "range", na.rm = T)) 
  
  da.r.pe <- data.table(pe = rep(pe, each = ncell(enm)), da.r)
  
  for(k in pe){
    da.pe <- da.r.pe[pe == k, -1]
    ens[] <- apply(da.pe, 1, function (x) sum(x * tss.05) / sum(tss.05))
    
    setwd("ensemble_wei")
    
    writeRaster(ens, paste0("ensemble_wei_aver_", i, "_", k, ".tif"), 
                format = "GTiff")
    
    setwd("..")
    
    print(paste0("Nice! The ensemble ", i, " for ", k, " it's done!"))}
  
  da <- data.table()
  ens[] <- NA}

###----------------------------------------------------------------------------###
