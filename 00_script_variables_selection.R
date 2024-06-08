### script variable selection ### 

# Thadeu Sobral de Souza - thadeusobral@gmail.com 
# Mauricio Humberto Vancine - mauricio.vancine@gmail.com

###-----------------------------------------------------------------------------------------###

# 1. clear the memory and load the packages
# clear workspace and increase memory
rm(list = ls())
memory.limit(size = 1.75e13) 

# install packages
#install.packages(c("raster", "rgdal", "corrplot", "RStoolbox", "vegan", "psych"), dep = T)
#install.packages("RStoolbox") 
# load packages
library(raster) # sig 
library(rgdal) # sig
library(corrplot) # correlation
library(RStoolbox)# raster pca
library(vegan) # multivariates
library(psych) # factorial analysis
library(terra) # multivariates

# check loaded packets
search()

###-----------------------------------------------------------------------------------------###

# 2. import data
# directory
setwd("~/Documents/mestrado/WD/enm/variáveis")
getwd()

# list name of archives
tif <- list.files(patt = ".tif$")
tif

# select name
#pres=presente #0k=atual
pres <- grep("5m", tif, value = T)
pres
#hol=holoceno #6k=5mil anos
hol <- grep("6k", tif, value = T)
hol
#lgm=last glaciation...21k=21mil anos
lgm <- grep("21k", tif, value = T)
lgm

# load variables ans rename
pres.s <- stack(pres)
pres.s

#so renomeia aqui no r pra ficar mais curto
names(pres.s)
names(pres.s) <-c(paste0("pres_", "bio0", 1:9), paste0("pres_", "bio", 10:18))
names(pres.s)

plot(pres.s)
plot(pres.s[[1]])

# bioclim
# http://www.worldclim.org/bioclim

hol.s <- stack(hol)
hol.s
names(hol.s) <- c(paste0("hol_", "bio0", 1:9), paste0("hol_", "bio", 10:18))
plot(hol.s)

lgm.s <- stack(lgm)
lgm.s
names(lgm.s) <- c(paste0("lgm_", "bio0", 1:9), paste0("lgm", "bio", 10:18))
plot(lgm.s)


###-----------------------------------------------------------------------------------------###

# 3. adjust raster to mask
# 3.1. limit of interest
# extension
sa <- extent(c(-180, 180, -90, 90)) # xmin, xmax, ymin, ymax; x=long, y=lat
plot(sa, xlab = "long", ylab = "lat", type = "n")
plot(pres.s[[1]], add = T)
plot(sa, cex.lab = 1.3, col = "red", add = T)

# cut
pres.sa <- crop(pres.s, sa) #s eh o plot atual, sa eh o plot que eu quero, onde sa eh a extensao que eu estipulei antes
pres.sa
pres.s

plot(pres.s[[1]])

par(mfrow = c(1, 2))
plot(pres.s$pres_bio01)
plot(pres.sa$pres_bio01)
dev.off()

hol.sa <- crop(hol.s, sa)
hol.sa
plot(hol.sa)

lgm.sa <- crop(lgm.s, sa)
lgm.sa
plot(lgm.sa)


# 3.2 shapefile
# import shapefile
br <- shapefile("BR_ESTADOS_IBGE.shp")
br
plot(br, col = "gray", axes = T)

# download shapefile
#getData vai fazer o download e salvar na sua pasta com arquivo .rda
br.s <- getData("GADM", country = 'BRA', level = 0)
br.s
plot(br.s, col = "green", axes = T)

# adjust to limit - usar um shape pra cortar o outro
ma <- mask(pres.s, br)
ma
plot(ma)
plot(ma[[1]])
plot(br, add = T)

# adjust to extention
cr <- crop(pres.s, br)
cr
plot(cr)
plot(cr[[1]])
plot(br, add = T)

# adjust to limit and extention
#SEMPRE pegar a mascara, ajustar a extensao e colocar a mascara novamente
ma <- mask(pres.s, br)
cr <- crop(ma, br)
plot(cr[[1]])
plot(br, add = T)

cr <- crop(pres.s, br)
ma <- mask(cr, br)
plot(ma[[1]])
plot(br, add = T)

# together, um comando pra tres coisas que fizemos antes separado
pres.br <- mask(crop(pres.s, br), br)
plot(pres.br)

hol.br <- mask(crop(hol.s, br), br)
plot(hol.br)

lgm.br <- mask(crop(lgm.s, br), br)
plot(lgm.br)


###-----------------------------------------------------------------------------------------###

# 4. extract values of cells
pres.sa.v <- values(pres.sa)
pres.sa.r<-pres.sa[] #mesma funcao do comando anterior
head(pres.sa.v)
#head eh o cabecalho (o numero dentro dos parenteses eh o numero de linhas que eu quero que seja mostrado)
head(pres.sa.v, 10)
head(pres.sa.v, 50)

# number of lines
nrow(pres.sa.v)

# dimension
dim(pres.sa.v)

# omit NAs #NA sao as celulas que nao tem valor (vazias)
pres.sa.v.na <- na.omit(pres.sa.v)
head(pres.sa.v.na, 50)
dim(pres.sa.v.na)

# there's NAs?
any(is.na(pres.sa.v.na))

###-----------------------------------------------------------------------------------------###

# 5.correlation

# create directory
dir.create("analise_selecao_variaveis") 

setwd("analise_selecao_variaveis") 

dir.create("correlacao") 

setwd("correlacao") 

getwd() 

# correlation
corr <- cor(pres.sa.v.na)
corr
round(corr, 2) #arredondar os valores com duas casas decimais 
abs(round(corr, 2)) #tirar os valores negativos, deixar tudo positivo
ifelse(corr >= 0.7, "sim", "nao") #0.7 eh o valor que ele usou de correlacao #se o valor de correlacao for maior que 0.7, ele coloca sim, se form menor, vai colocar nao
ifelse(corr >= 0.7, 1, 0)

# export pro excel
write.table(abs(round(corr, 2)), "cor_pres.xls", row.names = T, sep = "\t")
write.table(ifelse(corr >= 0.7, "sim", "nao"), "cor_pres_afirmacao.xls", 
      row.names = T, sep = "\t")

# plot of correlation
corrplot(corr, type = "lower", diag = F, tl.srt = 45, mar = c(3, 0.5, 2, 1),
	  title = "Correlacao entre variaveis Bioclimaticas")

# blue
corrplot(abs(corr), type = "lower", diag = F, tl.srt = 45, mar = c(3, 0.5, 2, 1),
	  title = "Correlacao entre variaveis Bioclimaticas")

# red
corrplot(-1 * (abs(corr)), type = "lower", diag = F, tl.srt = 45, mar = c(3, 0.5, 2, 1),
	  title = "Correlacao entre variaveis Bioclimaticas")

# exporta figure
tiff("cor_ma.tif", w = 18, he = 18, units = "cm", res = 300, comp = "lzw")

corrplot(corr, type = "lower", diag = F, tl.srt = 45, mar = c(3, 0.5, 2, 1),
	  title = "Correlacao entre variaveis Bioclimaticas")

dev.off()

###-----------------------------------------------------------------------------------------###

# 6. pca

# directory
setwd("C:\\Users\\tammyarai\\Documents\\UNICAMP\\_disciplina_enm_lavras_2018\\scripts_r\\enm\\01_dados") #voltar pro diretorio inicial
getwd() 
dir.create("pca") 
setwd("pca") 
getwd() 

# 6.1. pca to choose variables
# pca 
pca <- prcomp(pres.sa.v.na, scale = T)
pca

# contribution of each axis (eigenvalues - autovalores)
summary(pca)

# screeplot
par(mar = c(3, 5, 5, 2))
screeplot(pca, main = "Contribuicao de cada PC", ylab = "Autovalores", cex.lab = 1.3)
abline(h = 1, col = "red", lty = 2)

tiff("screeplot.tif", wi = 18, he = 18, un = "cm", res = 300, comp = "lzw")
par(mar = c(3, 5, 5, 2))
screeplot(pca, main = "Contribuicao de cada PC", ylab = "Autovalores", cex.lab = 1.3)
abline(h = 1, col = "red", lty = 2)
dev.off()

# values of each axis (eigenvectors - autovetores - escores)
pca$x
dim(pca$x)

# relation of the variables with each axis (loadings - cargas)
pca$rotation[, 1:5]
round(pca$rotation[, 1:5], 2)
abs(round(pca$rotation[, 1:5], 2))

# exportar tabela com a contribuicao
write.table(abs(round(pca$rotation[, 1:5], 2)), "contr_pca.xls", row.names = T, sep = "\t")

# plot
biplot(pca, var.axes = T, xlabs = rep("o", nrow(pca$x)), ylabs = paste0("bio", 1:19), cex = .8,
	 expand = 1.2, xlab = "PC1 (43.52%)", ylab = "PC2 (23.51%)", main = "PCA Bioclimaticas", 
	 xlim = c(-.03, .04))

# 6.2. pca as new variables
# pca of raster
pca.sa <- rasterPCA(pres.sa, spca = T) 
pca.sa

# contribution of each axis 
summary(pca.sa$model)
summary(pca)

# screeplot
par(mar = c(3, 5, 5, 2))
screeplot(pca.sa$model, main = "Contribuicao de cada PC", ylab = "Autovalores", cex.lab = 1.3)
abline(h = 1, col = "red", lty = 2)

# comparation
par(mfrow = c(1, 2))
screeplot(pca.sa$model, main = "Contribuicao de cada PC", ylab = "Autovalores", cex.lab = 1.3)
abline(h = 1, col = "red", lty = 2)

screeplot(pca, main = "Contribuicao de cada PC", ylab = "Autovalores", cex.lab = 1.3)
abline(h = 1, col = "red", lty = 2)

tiff("screeplot_raster.tif", wi = 18, he = 18, un = "cm", res = 300, comp = "lzw")
par(mar = c(3, 5, 5, 2))
screeplot(pca.sa$model, main = "Contribuicao de cada PC", ylab = "Autovalores", cex.lab = 1.3)
abline(h = 1, col = "red", lty = 2)
dev.off()

# plot pf pcas as variables
plot(pca.sa$map)

plot(pca.sa$map[[1:5]])


# export
# export one by one
writeRaster(pca.sa$map[[1]], "pc1_e .tif", format = "GTiff")
writeRaster(pca.sa$map[[2]], "pc2_e .tif", format = "GTiff")
writeRaster(pca.sa$map[[3]], "pc3_e .tif", format = "GTiff")
writeRaster(pca.sa$map[[4]], "pc4_e .tif", format = "GTiff")
writeRaster(pca.sa$map[[5]], "pc5_e .tif", format = "GTiff")

# export usinf for
for(i in 1:5){
 writeRaster(pca.sa$map[[i]], paste0("pc", i, ".tif"), format = "GTiff", 
             overwrite = T)}

###-----------------------------------------------------------------------------------------###

# 7. factorial analysis

# directory
setwd("C:\\Users\\tammyarai\\Documents\\UNICAMP\\_disciplina_enm_lavras_2018\\scripts_r\\enm\\01_dados") 
getwd() 
dir.create("fatorial") 
setwd("fatorial") 
getwd() 

# preliminaries analysis
# kmo e bartlett
KMO(cor(pres.sa.v.na)) # > 0.5
cortest.bartlett(cor(pres.sa.v.na), n = nrow(pres.sa.v.na)) # p < 0.05

# number os axis
# screeplot
fa <- fa.parallel(pres.sa.v.na, fa = "fa", fm = "ml")
fa

# exportar screeplot
tiff("screeplot_fatorial.tif", wi = 18, he = 18, un = "cm", res = 300, comp = "lzw")
fa.parallel(pres.sa.v.na, fm = "ml", fa = "fa") 
dev.off()

# fatorial
fa.sa <- fa(pres.sa.v.na, nfactors = 5, rotate = "varimax", fm = "ml")
sa.loadings <- loadings(fa.sa)
sa.loadings

abs(round(sa.loadings, 2))

# exportar tabela dos resultados
write.table(abs(round(sa.loadings, 2)), "as_loadings.xls", row.names = T, sep = "\t")

# bios
# bio02, bio04, bio10, bio16, bio17

# significado das bios
# BIO01 = Temperatura media anual
# BIO02 = Variacao da media diurna (media por mes (temp max - temp min))
# BIO03 = Isotermalidade (BIO02/BIO07) (* 100)
# BIO04 = Sazonalidade da temperatura (desvio padrao deviation *100)
# BIO05 = Temperatura maxima do mes mais quente
# BIO06 = Temperatura minima do mes mais frio
# BIO07 = Variacao da temperatura anual (BIO5-BIO6)
# BIO08 = Temperatura media do trimestre mais chuvoso
# BIO09 = Temperatura media do trimestre mais seco
# BIO10 = Temperatura media do trimestre mais quente
# BIO11 = Temperatura media do trimestre mais frio
# BIO12 = Precipitacao anual
# BIO13 = Precipitacao do mes mais chuvoso
# BIO14 = Precipitacao do mes mais seco
# BIO15 = Sazonalidade da precipitacao (coeficiente de variacao)
# BIO16 = Precipitacao do trimestre mais chuvoso
# BIO17 = Precipitacao do trimestre mais seco
# BIO18 = Precipitacao do trimestre mais quente
# BIO19 = Precipitacao do trimestre mais frio

###-----------------------------------------------------------------------------------------###

# 8. export variables
# bios
# bio02, bio04, bio10, bio16, bio17

pres.sa
names(pres.sa)

lista <- c(05, 11, 12, 13, 16, 17, 18)
lista

pres.sa[[06]]



# present
for(i in lista){
 writeRaster(pres.sa[[i]], ifelse(i < 10, paste0("CCSM_0k_am_bio0", i, " .tif"), 
		 paste0("CCSM_0k_am_bio", i, " .tif")), format = "GTiff")}

# holoceno
for(i in lista){
 writeRaster(hol.sa[[i]], ifelse(i < 10, paste0("CCSM_6k_am_bio0", i, " .tif"), 
		 paste0("CCSM_6k_am_bio", i, " .tif")), format = "GTiff")}

# lgm
for(i in lista){
 writeRaster(lgm.sa[[i]], ifelse(i < 10, paste0("CCSM_21k_am_bio0", i, " .tif"), 
		 paste0("CCSM_21k_am_bio", i, " .tif")), format = "GTiff")}

###-----------------------------------------------------------------------------------------###


