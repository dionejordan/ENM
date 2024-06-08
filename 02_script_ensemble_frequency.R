### script frequency ensemble ###

# Maur?cio Humberto Vancine - mauricio.vancine@gmail.com
##Thadeu Sobral de souza - thadeusobral@gmail.com

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
setwd("C:\\Users\\tammyarai\\Documents\\UNICAMP\\_disciplina_enm_lavras_2018\\projeto_humpback_enm\\bio_oracle\\results")

# enms
# list files
## usando $ para carregar apenas os formatos.tif
tif <- list.files(patt = ".tif$")
tif

enm <- raster(tif[[1]])
enm
plot(enm)

# evaluate
txt <- list.files(patt = ".txt$")
txt

eva <- lapply(txt, read.table)
eva
names(eva) <- txt
eva[[1]]

###----------------------------------------------------------------------------###

## frequency ensemble 
# lists
# species
sp <- ("L_sang")
sp

# gcms
gc <- c("marspec")
gc

# periods
pe <- c("0k", "6k", "21k")
pe

# algorithms
al <- c("Bioclim", "Gower", "Maha", "Maxent", "SVM")
al

# replicates
re <- 1:10
re

# ensembles
ens.re <- enm
ens.re[] <- 0
names(ens.re) <- "ens.re"
ens.re

ens.al <- enm
ens.al[] <- 0
names(ens.al) <- "ens.al"
ens.al

# for
for(i in sp){		
  tif.sp <- grep(i, tif, value = T)
  eva.sp <- eva[grep(i, names(eva))]
  
  for(j in gc){		
    tif.gc <- grep(j, tif.sp, value = T)
    eva.gc <- eva.sp[grep(j, names(eva.sp))]
       
    for(k in pe){		
      tif.pe <- grep(k, tif.gc, value = T)

      for(l in al){
	      tif.al <- grep(l, tif.pe, value = T)
        eva.al <- eva.gc[grep(l, names(eva.gc))]
              	           
	      for(m in re){		
          enm.al <- stack(tif.al)
          ens.re <- sum(ens.re, enm.al[[m]] >= eva.al[[1]][m, 1])}
              
          dir.create("ensemble_freq")
          setwd("ensemble_freq")
	        writeRaster(ens.re, paste0("ensemble_freq_", i, "_", j, "_", k, "_", l, ".tif"), 
	                    format = "GTiff")
	        
	        setwd("..")

	        ens.al <- sum(ens.al, ens.re)
		  	
	        ens.re[] <- 0}
      
      setwd("ensemble_freq")
      writeRaster(ens.al, paste0("ensemble_freq_", i, "_", j, "_", k, ".tif"), 
                  format = "GTiff")
      writeRaster(ens.al / (length(al) * length(re)), paste0("ensemble_freq_", i, 
                                                             "_", j, "_", k, "_bin.tif"), 
                  format = "GTiff")
          
      setwd("..")
      
      ens.al[] <- 0}}}

###----------------------------------------------------------------------------###



