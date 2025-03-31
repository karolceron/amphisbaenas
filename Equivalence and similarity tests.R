#################################
# Extracting Data from Rasters Using Coordinates
#################################
# loading the packages
library(raster)
library(maps)

# read WorldClim files
files <- list.files('G:/Meu Drive/Modelagem/script/wc5', pattern='.bil$', full.names=TRUE)
layers1 <-stack(files) 


# Cutting by the extent
e<-extent(-80, -50, -35, -5)
layers1 <- crop(layers1,e)

#plot
plot(layers1[[1]])

#loading the table with lat and long data (camura.txt or bolivica.txt)
species <- read.table("camura.txt", h=T, sep=",")
species

# Let us now organize coordinates into a spatial object:
coords<-data.frame(lon=species[,2], lat=species[,1])
coordinates(coords)<-c("lon","lat")

map()
points(coords, pch=16)

#extracting data
val<-extract(x=layers1, y=coords)
write.csv(val, file = "camura_bio.csv")

#######################################
########## Equivalence and Similarity tests
######################################


# loading the packages
library(ecospat)
library(ade4)
library(usdm)

camura<-read.csv("camura_bio.csv", header = T,sep = ",",stringsAsFactors=FALSE)
bolivica<-read.csv("bolivica_bio.csv", header = T, sep = ",",stringsAsFactors=FALSE) 

camura$X<-1
bolivica$X<-1

# environmental pca 
combined<-data.frame(rbind(camura[,2:20],bolivica[,2:20]))
pca.env<-dudi.pca(combined, scannf = FALSE,nf=2)

#same thing, different version (pca ambiental)
pca.env<-dudi.pca(rbind(camura,bolivica)[,2:20], scannf=FALSE, nf=2)
###################################################################

# Plot	Variables	Contribution
ecospat.plot.contrib(contrib=pca.env$co, eigen=pca.env$eig)

#	predict	the	scores	on	the	axes
scores.globclim<-pca.env$li #	PCA	scores	for	the	whole	study	area

dim(camura)
camura[,20]<-1
scores.sp.camura<-suprow(pca.env,camura[which(camura[,1]==1),2:20])$li

dim(bolivica)
bolivica[,20]<-1
scores.sp.bolivica<-suprow(pca.env,bolivica[which(bolivica[,1]==1),2:20])$li


#PCA	scores
scores.clim.camura<-suprow(pca.env,camura[,2:20])$li 

scores.clim.bolivica<-suprow(pca.env,bolivica[,2:20])$li 	


###	Calculate	the	Occurrence	Densities	Grid	with	ecospat.grid.clim.dyn()

#	For	camura
grid.clim.camura<-ecospat.grid.clim.dyn(glob=scores.globclim,
                                       glob1=scores.clim.camura, sp=scores.sp.camura, R=100)

#	For bolivica
grid.clim.bolivica<-ecospat.grid.clim.dyn(glob=scores.globclim,
                                       glob1=scores.clim.bolivica, sp=scores.sp.bolivica, R=100)

#	Calculate	Niche	Overlap
ecospat.niche.overlap (grid.clim.camura, grid.clim.bolivica, cor=TRUE)

#	Perform	the	Niche	Equivalency	Test
eq.test<-ecospat.niche.equivalency.test(grid.clim.camura,
                                        grid.clim.bolivica, rep=100, overlap.alternative = "higher", ncores = 1)
eq.test


#	Niche	Similarity	Test
sim.test<-ecospat.niche.similarity.test(grid.clim.camura,grid.clim.bolivica,
                                        rep=100, overlap.alternative = "higher")
sim.test

sim.test<-ecospat.niche.similarity.test(grid.clim.bolivica,grid.clim.camura,
                                        rep=100, overlap.alternative = "higher")
sim.test

#	Plot	Equivalency	test
ecospat.plot.overlap.test(eq.test, "D", "Equivalency")

#	Plot	Similarity	test
ecospat.plot.overlap.test(sim.test, "D", "Similarity")

#	Delimiting	niche	categories	and	quantifying	niche	dynamics	in	analogue	climates
ecospat.niche.dyn.index (grid.clim.camura, grid.clim.bolivica,
                         intersection=0.1)

#Visualizing	niche	categories,	niche	dynamics	and	climate	analogy	between	ranges	

ecospat.plot.niche.dyn(grid.clim.camura, grid.clim.bolivica, quant=0.8,
                       interest=1, title= "Niche Overlap", name.axis1="PC1",
                       name.axis2="PC2",  colZ1="green3", colZ2= "red3", col.unf =
                         "green", col.exp = "red", col.stab = "blue")

ecospat.plot.niche.dyn(grid.clim.bolivica, grid.clim.camura, quant=0.8,
                       interest=1, title= "Niche Overlap", name.axis1="PC1",
                       name.axis2="PC2",  colZ1="green3", colZ2= "red3")




