  install.packages(("maxnet"), dep=T)
  rm(list=ls())
  library('rgdal')
  library('raster')
  library('dismo')
  library('gam')
  library('randomForest')
  library('kernlab')
  library('rJava')

  bio <- 'bio'
  
  #climatic layers
  
  setwd('/home/iwasarai/Documents/ENMs_Ampithoe')
  getwd()
  #abrindo variaveis climaticas
  ti <- list.files(patt = "tif")
  ti
  
  #presente
  pres <- grep ("0k", ti, value = T)
  pres <- stack(pres)
  names(pres) <- paste0("bio", c("05", "11", "12", "16", "17", "18"))
  plot(pres)
  
  
  
  
  
  ## extraindo coordenadas
  id.na <- na.omit(cbind(1:ncell(pres), values(pres)[,1]))
  coords <- xyFromCell(pres, id.na[,1])
  colnames(coords) <- c('Long','Lat')
  head(coords)
  
  
  #species presences, absences and fossil records
  pr.points <- read.table('dist_lsang.txt', h=T)
  pr.points
  
  getwd() 
  dir.create("results") 
  setwd("results") 
  getwd() 
  

  #colocar quais algoritmos eu vou usar, o Thadeu colocou um algortimo de 
  #envelope, dois de aprendizado de maquina e dois de alguma outra coisa que 
  #eu esqueci o nome
 
  x<-factor(pr.points[,1])
  levels(x)
  x
  
  
   for(i in 1:length(levels(x))){
  
  #criar um objeto NULL (vazio), pra que no final do script eu saber quantas 
    #replicas eu vou salvar no final nessa pasta
    
  eval.Bioclim <- NULL
  eval.Gower <- NULL
  eval.Maha <- NULL
  eval.Maxent <- NULL
  eval.SVM <- NULL
  
  eval.names <- NULL
  
  #selecting presences and absences for each species
  #numero de ocorrencias que eu tenho
  #se eu tiver, coloco 1, se nao tiver coloca 0
  	id.specie <- levels(x)[i]
  	pr.specie <- pr.points[which(x == id.specie), 2:3]
  	
  	#aleatorizar o mesmo numero de pontos que o numero de especies nesse id.background
  	id.background <- sample(nrow(coords), nrow(pr.specie))
  	#informar quais sao os pontos de ocorrencia e quais sao os potos do modelo
  	bc.specie <- coords[id.background,]
  	
  #pra ver quantas replicas eu vou fazer, dando loop variando de 1 ate 10
  for(r in 1:10){	
  ##building ENMs
  #data		
    #alteatorizar 75% dos meus pontos de ocorrencia pra que crie um objeto de treino
    # o que nao entrar pra treino entra pra teste
  	pr.sample.train <- sample(nrow(pr.specie), round(0.75 * nrow(pr.specie)))
  	bc.sample.train <- sample(nrow(bc.specie), round(0.75 * nrow(bc.specie)))
  	test <- na.omit(prepareData(x= pres, p= pr.specie[-pr.sample.train,], b= bc.specie[-bc.sample.train,]))
  	train <- na.omit(prepareData(x= pres, p= pr.specie[pr.sample.train,], b= bc.specie[bc.sample.train,]))
  
  ##Bioclim	
  	Bioclim <- bioclim(train[which(train[,1]==1), -1])	
  	
  	writeRaster(predict(pres, Bioclim), paste(bio, '_Bioclim_0k_', id.specie, r, ".tif", sep=""), format= "GTiff")	
  	
  	eBioclim <- evaluate(p= test[test[,1]==1, -1], a= test[test[,1]==0, -1], model= Bioclim)
  	idBioclim <- which(eBioclim@t == as.numeric(threshold(eBioclim, 'spec_sens')))
  	#usei o threshold de maxima especificidade e sensibilidade ("spec_sens")
  	
  	eval.Bioclim.sp <- c(eBioclim@t[idBioclim], eBioclim@auc, (eBioclim@TPR[idBioclim]+eBioclim@TNR[idBioclim]-1))
  	eval.Bioclim <- rbind(eval.Bioclim, eval.Bioclim.sp)
  	
  	?threshold
  
  ##Gower	
  	Gower <- domain(train[which(train[,1]==1), -1])	
  	
  
    writeRaster(predict(pres, Gower), paste(bio, '_Gower_0k_', id.specie, r, ".tif", sep=""), format= "GTiff")	
    
  	
  	eGower <- evaluate(p= test[test[,1]==1, -1], a= test[test[,1]==0, -1], model= Gower)
  	idGower <- which(eGower@t == as.numeric(threshold(eGower, 'spec_sens')))
  	eval.Gower.sp <- c(eGower@t[idGower], eGower@auc, (eGower@TPR[idGower]+eGower@TNR[idGower]-1))
  	eval.Gower <- rbind(eval.Gower, eval.Gower.sp)
  
  
  ##Maha	
  	Maha <- mahal(train[which(train[,1]==1), -1])	
  	
    writeRaster(predict(pres, Maha), paste(bio, '_Maha_0k_', id.specie, r, ".tif", sep=""), format= "GTiff")	
   
  	eMaha <- evaluate(p= test[test[,1]==1, -1], a= test[test[,1]==0, -1], model= Maha)
  	idMaha <- which(eMaha@t == as.numeric(threshold(eMaha, 'spec_sens')))
  	eval.Maha.sp <- c(eMaha@t[idMaha], eMaha@auc, (eMaha@TPR[idMaha]+eMaha@TNR[idMaha]-1))
  	eval.Maha <- rbind(eval.Maha, eval.Maha.sp)
  	
  
  ##Maxent	
  	Maxent <- maxent(train[,-1], train[,1])	
  
    writeRaster(predict(pres, Maxent), paste(bio, '_Maxent_0k_', id.specie, r, ".tif", sep=""), format= "GTiff")  
    
  	eMaxent <- evaluate(p= test[test[,1]==1, -1], a= test[test[,1]==0, -1], model= Maxent)
  	idMaxent <- which(eMaxent@t == as.numeric(threshold(eMaxent, 'spec_sens')))
  	eval.Maxent.sp <- c(eMaxent@t[idMaxent], eMaxent@auc, (eMaxent@TPR[idMaxent]+eMaxent@TNR[idMaxent]-1))
  	eval.Maxent <- rbind(eval.Maxent, eval.Maxent.sp)
  	
  
  	eval.names <- c(eval.names, paste(id.specie, r, sep=""))		
  }#ends for'r'
  
  dimnames(eval.Bioclim) <- list(eval.names, c('thrs','AUC','TSS'))
  dimnames(eval.Gower) <- list(eval.names, c('thrs','AUC','TSS'))
  dimnames(eval.Maha) <- list(eval.names, c('thrs','AUC','TSS'))
  dimnames(eval.Maxent) <- list(eval.names, c('thrs','AUC','TSS'))

  
  write.table(eval.Bioclim, paste("zEval_", bio, '_Bioclim_', id.specie, ".txt", sep=""))
  write.table(eval.Gower, paste("zEval_", bio, '_Gower_', id.specie, ".txt", sep=""))
  write.table(eval.Maha, paste("zEval_", bio, '_Maha_', id.specie, ".txt", sep=""))
  write.table(eval.Maxent, paste("zEval_", bio, '_Maxent_', id.specie, ".txt", sep=""))
 
  }#ends for'i'
  
