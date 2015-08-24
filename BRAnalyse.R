###############################################################################
###############################################################################
#Script for the analyses of the Microcyclus ulei populations in Mato Grosso 
###############################################################################
###############################################################################

library(maptools)
library(maps)
library(mapdata)
library(RColorBrewer)
library(adegenet)
library(raster)
library(RColorBrewer)
library(rgdal)
library(combinat)
library(pegas)
library(gdistance)
library(vegan)
library(Geneland)

setwd("~/work/Rfichiers/Githuber/BRAMU_data")


###############################################################################
#Loading the genetic data and preparing the datafile for other softwares
###############################################################################

#first of all, we load the genetic dataset
BRA<-read.table("BRAdata.txt",header=T,sep="\t")
#here is the structure of the datafile, for explanation of each columns, 
#see ReadMe.txt file in DRYAD repository
head(BRA)
#a summary of the different variables
summary(BRA)
#number of individuals in each sampled populations
table(BRA$pop_ID)
sum(table(BRA$pop_ID)) #474 individuals

#two clone-corrected datasets are build on the complete dataset. The first one 
#is a 'conservative' clone correction since we only keep on MLG per site
BRAcccons<-BRA[BRA$cc_cons==1,]
sum(table(BRAcccons$pop_ID)) #only 259 individuals left
#the second clone corrected dataset is less conservative, we only removed 
#over-represented multicopies MLG (see Mat & Meth for details)
BRAcc<-BRA[BRA$cc==1,]
sum(table(BRAcc$pop_ID)) #only 338 individuals left

#### STRUCTURE file format
#a function for exporting the file to a "STRUCTURE" file format, this function 
#is only working with the appropriate input datafile (ie a file format similar 
#to "BRAdata.txt")
structexport<-function(fileinput) {
  temp<-fileinput[,c(1:2,14:27)]
  temp[temp=="0"]<-"-9"
  levels(temp$pop_ID)<-c(9,4,5,6,8,15,1,2,13,14,3,7,12,10,11)
  write.table(temp,file="output.str",row.names=F)
}

#only a few edition to the output file are needed before running STRUCTURE 
#on it (remove "", and the two first column headers)
structexport(BRA)
#change the name of the 'output.str' file before running the function again, 
#or else the file will be overwrite
structexport(BRAcccons)
structexport(BRAcc)

#### TESS file format
#a function for exporting the file to a "TESS" file format, this function is 
#only working with the appropriate input datafile (ie a file format similar 
#to "BRAdata.txt")
tessexport<-function(fileinput) {
  temp<-fileinput[,c(1:2,4:5,14:27)]
  #we add some noise to geographic coordinates, so each individual has 
  #different coordinates
  temp[,3]<-jitter(temp[,3])
  temp[,4]<-jitter(temp[,4])
  #missing data should be formatted with negative number
  temp[temp=="0"]<-"-9"
  levels(temp$pop_ID)<-c(9,4,5,6,8,15,1,2,13,14,3,7,12,10,11)
  write.table(temp,file="output.tess",row.names=F)
}

tessexport(BRA)
#remove "" in the output datafile, then change the name of the 'output.tess' 
#file before running the function again, or else the file will be overwrite
tessexport(BRAcccons)
tessexport(BRAcc)

#### GENELAND file format
#a function for exporting the file to "GENELAND" file format, this function is 
#only working with the appropriate input datafile (ie a file format similar to 
#"BRAdata.txt"), it is only necessary if you want to use the graphical user 
#interface of Geneland ('Geneland.GUI()')
genelandexport<-function(fileinput) {
  tempgeo<-fileinput[,c(4,5)]
  tempgeo<-SpatialPoints(tempgeo,proj4string=CRS("+proj=longlat +datum=WGS84"))
  tempgeo<-spTransform(tempgeo, CRS("+proj=utm +zone=21 +datum=WGS84"))
  tempgen<-fileinput[,c(14:27)]
  #missing data should be formatted with "NA"
  tempgen[tempgen=="0"]<-"NA"
  write.table(tempgeo,file="outputgeo.geneland",row.names=F,col.names=FALSE,
              quote=FALSE,sep=" ")
  write.table(tempgen,file="outputgen.geneland",row.names=F,col.names=FALSE,
              quote=FALSE,sep=" ")
}

genelandexport(BRA)
#Change the name of the 'outputgeo.geneland' and 'outputgen.geneland' files 
#before running the function again, or else the file will be overwrite
genelandexport(BRAcccons)
genelandexport(BRAcc)


###############################################################################
#Identifying the best K for STRUCTURE run
###############################################################################

#Analyzes were performed using STRUCTURE2.3.4 software, with a model allowing 
#admixture and correlation of allele frequencies. Each run consisted of a 
#burn-in period of 100.000 iterations followed by 500.000 simulations. Fifteen 
#repetitions of each run were performed for K ranging from 1 to 15. before 
#importing the file, replace white space in the column header names with 
#underscore, replace "?1" by "alpha", and remove double white spaces or it
#will provoc importation problem or failure
resstr<-read.table(file="BRAstr.out", header=T,sep=" ",blank.lines.skip=T)
resccstr<-read.table(file="BRAccstr.out", header=T,sep=" ",blank.lines.skip=T)
resccconsstr<-read.table(file="BRAccconsstr.out", header=T,sep=" ",
                         blank.lines.skip=T)

#new version of the file with 15 repetitions of each run
resstr<-read.table(file="BRAoutput.out", header=T,sep=" ",blank.lines.skip=T)
resccstr<-read.table(file="BRAccoutput.out", header=T,sep=" ",
                     blank.lines.skip=T)
resccconsstr<-read.table(file="BRAccconsoutput.out", header=T,sep=" ",
                         blank.lines.skip=T)


#a function which compute delta K values, nb_K is the number of different K 
#considered, and nb_rep is the number of repetition of each K
chooseK<-function(str_out,nb_K,nb_rep) {
  datatable<-data.frame("K"=c(rep(1:nb_K,each=nb_rep)),"Ln(Pd)"=str_out[,4])
  Lprim<-c(rep("NA",nb_rep))
  for (i in ((nb_rep+1):(nb_K*nb_rep))) {
    Lprim<-c(Lprim,str_out[i,4]-str_out[i-nb_rep,4])
  }
  datatable<-data.frame(datatable,as.numeric(Lprim))
  Lsecond<-c(rep("NA",nb_rep))
  for (i in (((2*nb_rep)+1):(nb_K*nb_rep))) {
    Lsecond<-c(Lsecond,abs(datatable[i,3]-datatable[i-nb_rep,3]))
  }
  Lsecond<-c(Lsecond,rep("NA",nb_rep))
  datatable<-data.frame(datatable,as.numeric(Lsecond))
  reztable<-data.frame("K"=c(1:nb_K))
  meanL<-c()
  sdL<-c()
  for (i in (1:nb_K)) {
    meanL<-c(meanL,mean(datatable[datatable$K==i,2]))
    sdL<-c(sdL,sd(datatable[datatable$K==i,2]))
  }
  reztable<-data.frame(reztable,meanL,sdL)
  meanLprime<-c()
  sdLprime<-c()
  for (i in (1:nb_K)) {
    meanLprime<-c(meanLprime,mean(as.numeric(datatable[datatable$K==i,3])))
    sdLprime<-c(sdLprime,sd(datatable[datatable$K==i,3]))
  }
  reztable<-data.frame(reztable,meanLprime,sdLprime)
  meanLsecond<-c()
  sdLsecond<-c()
  for (i in (1:nb_K)) {
    meanLsecond<-c(meanLsecond,mean(as.numeric(datatable[datatable$K==i,4])))
    sdLsecond<-c(sdLsecond,sd(datatable[datatable$K==i,4]))
  }
  reztable<-data.frame(reztable,meanLsecond,sdLsecond)
  deltaK<-c()
  for (i in (1:nb_K)) {
    deltaK<-c(deltaK,reztable[reztable$K==i,6]/reztable[reztable$K==i,3])
  }
  reztable<-data.frame(reztable,deltaK)
  return(reztable)
}

deltastr<-chooseK(resstr,15,15)
deltaccstr<-chooseK(resccstr,15,15)
deltaccconsstr<-chooseK(resccconsstr,15,15)

#a function to plot variation of Delta K and Ln(P(X|K)) with K. 'datadeltak' 
#is the output file of 'chooseK' function, and nb_K is the number of different
#K considered
plotdeltaK<-function(datadeltaK,nb_K,titre){
  op<-par(pty="s")
  plot(datadeltaK[1:(nb_K-2),8],type="b",pch=24,cex=2.5,lwd=4,lty=1,
       col="transparent",bg="white",bty="n",ann=F)
  par(new=TRUE)
  plot(datadeltaK[1:(nb_K-2),8],type="b",pch=24,bty="n",xaxt="n",yaxt="n",
       ann=F,cex=2.5,lwd=4,lty=1)
  axis(side=1,at=seq(1,13,1),lwd=3,font.axis=2)
  axis(side=2,lwd=3,font.axis=2)
  title(ylab="Delta K",font.lab=2,cex.lab=1.5)
  par(new=TRUE)
  plot(datadeltaK[1:(nb_K-2),2],type="b",pch=22,cex=2.5,lwd=4,lty=2,
       col="grey50",bg="white",bty="n",xaxt="n",yaxt="n",ann=F)
  axis(side=4,lwd=3,font.axis=2,col="grey50")
  mtext("Ln(P(X|K))", side=4, line=4,font=2,cex=1,col="grey50")
  title(main=titre,xlab="K",font.lab=2,cex.lab=1.5,cex.main=2)
  par(op)
}

plotdeltaK(deltastr,15,"Complete dataset (n=474)")
plotdeltaK(deltaccstr,15,"Clone Corrected dataset (n=338)")
plotdeltaK(deltaccconsstr,15,"Conservative Clone Corrected dataset (n=259)")

#you can obtain the same figure as in the manuscript by exporting the plot to 
#pdf format, with a width of 12 inches and an height of 11 inches

#You can also obtain a combined plot for the three different dataset
op<-par(mfrow=c(1,3))
plotdeltaK(deltastr,15,"Complete dataset (n=474)")
plotdeltaK(deltaccstr,15,"Clone Corrected dataset (n=338)")
plotdeltaK(deltaccconsstr,15,"Conservative Clone Corrected dataset (n=259)")
par(op) 
#then export with a width of 32 inches and an height of 10 inches


###############################################################################
#Loading geographical data / informations / rasters / shapefiles
###############################################################################

#before importing the files, don't forget to set the correct working directory

#administrative limit and roads
municipioSH<-readShapePoly("mucipio_MT.shp",
                           proj4string=CRS("+proj=longlat +datum=WGS84"))
voies<-readShapeLines("voies.shp",
                      proj4string=CRS("+proj=longlat +datum=WGS84"))

#shapefile of the Hevea brasiliensis plantations identified by analysis of 
#satellite images
parcPoly<-readShapePoly("parc_MT.shp",
                        proj4string=CRS("+proj=longlat +datum=WGS84"))
#in order to have the area of the different plantations in square meters and 
#the distance in meters, we turn the coordinates of the object into planar 
#coordinates format
parcPoly.utm <- spTransform(parcPoly, CRS("+proj=utm +zone=21 +datum=WGS84"))
#exploring the structure of the "polygons" object
str(parcPoly[1,])
parcPoly.utm[1,]@polygons
parcPoints<-readShapePoints("parc_point.shp",
                            proj4string=CRS("+proj=longlat +datum=WGS84"))
#in order to have the distance in meters, we turn the coordinates of the 
#object into planar coordinates format
parcPoints.utm <- spTransform(parcPoints,
                              CRS("+proj=utm +zone=21 +datum=WGS84"))
sampPoints<-readShapePoints("prelevement.shp",
                            proj4string=CRS("+proj=longlat +datum=WGS84"))
#because there was no infection in 'Campo Verde' location, we remove this 
#sampling point
sampPoints<-sampPoints[-c(2),]
levels(sampPoints@data$name_id)[levels(sampPoints@data$name_id)=="SOR1"]<-"SOR"
#in order to have the distance in meters, we turn the coordinates of the 
#object into planar coordinates format
sampPoints.utm <- spTransform(sampPoints,
                              CRS("+proj=utm +zone=21 +datum=WGS84"))

#loading the altitudinal data, no need to clip the information, because it has 
#been already clipped
alt<-raster("alt_Clip1.img")

#loading the raster of the land cover, no need to clip the information, 
#because it has been already clipped
veget<-raster("couvert_veget1.img")
#simplification of the 'veget' raster to keep only 'veget' categories
#where Hevea brasiliensis can occur (40, 41 and 42, see 
#http://postel.mediasfrance.org/fr/)
temp<-c(13,39,NA,43,210,NA,40,42,40)
temp <- matrix(temp, ncol=3, byrow=TRUE)
veget<-reclassify(veget,temp)
#there is no H. brasiliensis under 14°S approximatly, 
temp<-c(extent(veget)@xmin,-14.3)
temp<-matrix(temp,ncol=2,byrow=T)
#determining the first cell of the grid we want to modify
extract(veget,temp,cellnumbers=TRUE)[1] 
#We turn classification that don't fit our criteria in missing data
veget[extract(veget,temp,cellnumbers=TRUE)[1]:ncell(veget)]<-NA

#extraction of information from the different raster using the coordinates of 
#sampling points
bioclim1<-raster("bio_1_Clip1.img")
bioclim2<-raster("bio_2_Clip1.img")
bioclim3<-raster("bio_3_Clip1.img")
bioclim4<-raster("bio_4_Clip1.img")
bioclim5<-raster("bio_5_Clip1.img")
bioclim6<-raster("bio_6_Clip1.img")
bioclim7<-raster("bio_7_Clip1.img")
bioclim8<-raster("bio_8_Clip1.img")
bioclim9<-raster("bio_9_Clip1.img")
bioclim10<-raster("bio_10_Clip1.img")
bioclim11<-raster("bio_11_Clip1.img")
bioclim12<-raster("bio_12_Clip1.img")
bioclim13<-raster("bio_13_Clip1.img")
bioclim14<-raster("bio_14_Clip1.img")
bioclim15<-raster("bio_15_Clip1.img")
bioclim16<-raster("bio_16_Clip1.img")
bioclim17<-raster("bio_17_Clip1.img")
bioclim18<-raster("bio_18_Clip1.img")
bioclim19<-raster("bio_19_Clip1.img")
site_table<-data.frame("site_ID"=sampPoints$name_id,
                       coordinates(sampPoints)[,1:2],
                       "Altitude"=extract(alt,
                                          coordinates(sampPoints)[,1:2]), 
                       "BioClim1"=extract(bioclim1,
                                          coordinates(sampPoints)[,1:2])/10,
                       "BioClim2"=extract(bioclim2,
                                          coordinates(sampPoints)[,1:2])/10, 
                       "BioClim3"=extract(bioclim3,
                                          coordinates(sampPoints)[,1:2]), 
                       "BioClim4"=extract(bioclim4,
                                          coordinates(sampPoints)[,1:2])/100, 
                       "BioClim5"=extract(bioclim5,
                                          coordinates(sampPoints)[,1:2])/10, 
                       "BioClim6"=extract(bioclim6,
                                          coordinates(sampPoints)[,1:2])/10, 
                       "BioClim7"=extract(bioclim7,
                                          coordinates(sampPoints)[,1:2])/10, 
                       "BioClim8"=extract(bioclim8,
                                          coordinates(sampPoints)[,1:2])/10, 
                       "BioClim9"=extract(bioclim9,
                                          coordinates(sampPoints)[,1:2])/10, 
                       "BioClim10"=extract(bioclim10,
                                           coordinates(sampPoints)[,1:2])/10, 
                       "BioClim11"=extract(bioclim11,
                                           coordinates(sampPoints)[,1:2])/10, 
                       "BioClim12"=extract(bioclim12,
                                           coordinates(sampPoints)[,1:2]), 
                       "BioClim13"=extract(bioclim13,
                                           coordinates(sampPoints)[,1:2]), 
                       "BioClim14"=extract(bioclim14,
                                           coordinates(sampPoints)[,1:2]), 
                       "BioClim15"=extract(bioclim15,
                                           coordinates(sampPoints)[,1:2]), 
                       "BioClim16"=extract(bioclim16,
                                           coordinates(sampPoints)[,1:2]), 
                       "BioClim17"=extract(bioclim17,
                                           coordinates(sampPoints)[,1:2]), 
                       "BioClim18"=extract(bioclim18,
                                           coordinates(sampPoints)[,1:2]), 
                       "BioClim19"=extract(bioclim19,
                                           coordinates(sampPoints)[,1:2]), 
                       stringsAsFactors = FALSE)
colnames(site_table)<-c("site_ID","Longitude","Latitude","Altitude",
                  "Annual Mean Temperature",
                  "Mean Diurnal Range (Mean of monthly (max temp - min temp))",
                  "Isothermality",
                  "Temperature Seasonality (standard deviation *100)", 
                  "Max Temperature of Warmest Month",
                  "Min Temperature of Coldest Month", 
                  "Temperature Annual Range",
                  "Mean Temperature of Wettest Quarter", 
                  "Mean Temperature of Driest Quarter", 
                  "Mean Temperature of Warmest Quarter", 
                  "Mean Temperature of Coldest Quarter",
                  "Annual Precipitation", 
                  "Precipitation of Wettest Month",
                  "Precipitation of Driest Month",
                  "Precipitation Seasonality (Coefficient of Variation)", 
                  "Precipitation of Wettest Quarter",
                  "Precipitation of Driest Quarter", 
                  "Precipitation of Warmest Quarter",
                  "Precipitation of Coldest Quarter")
site_table$site_ID<-as.character(site_table$site_ID)
site_table[site_table$site_ID=="SOR1",1]<-c("SOR")
#we ordered this table to match the organization of the genetic table
site_table<-site_table[order(site_table$site_ID),]
site_table


###############################################################################
#some plotting examples of the files of geographical data sets
###############################################################################

plot(municipioSH)
plot(voies,col="red",add=TRUE)
plot(parcPoly,col="red")
#usefull for illustration only, there may be conflict with the 
#zoom function of ape. Use detach(package:xx) to solve this.
raster::zoom(parcPoly,col="red")
plot(parcPoly,col="yellow",lwd=0.1)
plot(sampPoints,add=TRUE, col="black",bg="red",pch=21,cex=1.5)
text(coordinates(parcPoly),labels=round(plant_area), cex=0.2)
plot(alt)

#an example of a map that can be produced with the data loaded
op<-par(pty="s")
image(alt,col=brewer.pal(9,"Greys"))
image(veget,col="darkgreen",add=TRUE)
plot(municipioSH,add=TRUE)
plot(voies,add=TRUE,col="green")
#plot(parcPoly,add=TRUE, col="red") #not very usefull because area of 
#plantation are too small for the considered scale
plot(parcPoints,add=TRUE, col="black",bg="blue",pch=21)
plot(sampPoints,add=TRUE, col="black",bg="red",pch=21)
par(op)
#you can export the file with a large size of file (like 50 by 50 inches)


#submap of South America
op<-par(pty="s")
image(alt,col=brewer.pal(9,"Greys"))
image(veget,col="darkgreen",add=TRUE)
#plot(municipioSH,add=TRUE)
plot(parcPoints,add=TRUE, col="black",bg="blue",pch=21)
plot(sampPoints,add=TRUE, col="black",bg="red",pch=21)
polygon(x=c(-53.39,-53.13,-53.13,-53.39,-53.39),
        y=c(-13.29,-13.29,-13.03,-13.03,-13.29),col="white",
        density=35,border="white",lwd=1.5)
box()
par(op)
#the following code is adapted from: 
#http://www.stat.auckland.ac.nz/~paul/RGraphics/examples-map.R
#plotting an inset map of south america
maplocs <- map("worldHires", ylim=c(-60,12),xlim=c(-89,-35),
               col="lightblue1",fill=TRUE,add=TRUE,plot=FALSE)
xrange <- range(maplocs$x, na.rm=TRUE)
yrange <- range(maplocs$y, na.rm=TRUE)
aspect <- abs(diff(yrange))/abs(diff(xrange))
# customised to 6.5 by 4.5 figure size
par(fig=c(0.93 - 0.15, 0.929, 0.72, 0.72 + 0.15*aspect), mar=rep(0, 4), 
    ew=TRUE, pty="s")
plot.new()
plot.window(xlim=xrange, ylim=yrange)
polygon(x=c(-88.5,-34.5,-34.5,-88.5,-88.5),y=c(-58,-58,15,15,-58),
        col="lightblue1")
map("worldHires", ylim=c(-60,12),xlim=c(-89,-35), col="white",fill=T,add=TRUE)
polygon(x=c(-59.5,-52.5,-52.5,-59.5,-59.5),y=c(-17.5,-17.5,-10.5,-10.5,-17.5),
        col="grey")
#plotting an inset map of a zoom on a part of Mato Grosso plantations
par(fig=c(0.93 - 0.15, 0.9298, 0.4, 0.4 + 0.15), mar=rep(0, 4),new=TRUE,
    pty="s")
plot.new()
plot.window(xlim=c(-53.39,-53.13), ylim=c(-13.29,-13.03))
polygon(x=c(-53.39,-53.13,-53.13,-53.39,-53.39),
        y=c(-13.29,-13.29,-13.03,-13.03,-13.29),col="white")
plot(parcPoly,col="red",xlim=c(-53.39,-53.13),ylim=c(-13.29,-13.03),
     lwd="0.2",add=TRUE)
#export in pdf format with a size of 9.5*9.5 inches in order to obtain the 
#same figure

map("world",regions=c("brazil","argentina","uruguay","paraguay","chile",
                      "bolivia","peru","ecuador","colombia","venezuela",
                      "guyana","suriname","french guiana","panama",
                      "costa rica","nicaragua"))


###############################################################################
#Distance matix, Minimum spanning tree network or least cost pathway analysis
###############################################################################

#first we built the combined table of coordinates of plantations and sampling 
#sites
xy.utm<-rbind(coordinates(parcPoints.utm),coordinates(sampPoints.utm)[,c(1,2)])
row.names(xy.utm)[577:591]<-c("BBU","DEN","DAQ1","DAQ2","DAQ3","NHO","PEM1",
                              "PEM2","PGA1","PGA2","RND","ROS","SRC1","SRC2",
                              "SOR")
colnames(xy.utm)<-c("Longitude", "Latitude")
matdist.utm<-vegdist(xy.utm,method="euclidean")
mateucl.utm<-vegdist(xy.utm,method="euclidean")

#Euclidean distances beetween sites
mateucl.utm<-as.matrix(mateucl.utm)
dist_eucl<-mateucl.utm[577:591,577:591]
dist_eucl<-dist_eucl/1000 #we turn the distance in kilometers
#we ordered this table to match the organization of the genetic table
dist_eucl<-dist_eucl[order(colnames(dist_eucl)),order(colnames(dist_eucl))]

#Distances from the Minimum spanning tree network between Hevea brasiliensis 
#plantation. We use this function to find the minimum spanning tree
MinSpaTre.utm<-spantree(matdist.utm) 
#it is then possible to compute the distances between the nodes of the tree, 
#following a tree path
coph.utm<-cophenetic(MinSpaTre.utm)
#we just want the distance between the sampling points
matcoph.utm<-as.matrix(coph.utm)
dist_MinSpanTree<-matcoph.utm[577:591,577:591]
dist_MinSpanTree<-dist_MinSpanTree/1000 #we turn the distance in kilometers
#we ordered this table to match the organization of the genetic table
dist_MinSpanTree<-dist_MinSpanTree[order(colnames(dist_MinSpanTree)),
                                   order(colnames(dist_MinSpanTree))]

#Distances between sampling point using roads. This distances has been 
#computed using google maps https://maps.google.fi/maps?hl=en&tab=wl in July 
#2013 by copying and pasting coordinates of the plantation
dist_roads<-as.matrix(read.table("dist_roads.txt",header=TRUE))

#in order to plot the tree on the map, we built another tree with WGS84 
#coordinates
xy<-rbind(coordinates(parcPoints),coordinates(sampPoints)[,c(1,2)])
row.names(xy)[577:591]<-c("BBU","DEN","DAQ1","DAQ2","DAQ3","NHO","PEM1","PEM2",
                          "PGA1","PGA2","RND","ROS","SRC1","SRC2","SOR")
colnames(xy)<-c("Longitude", "Latitude")
matdist<-vegdist(xy,method="euclidean")
#we use this function to find the minimum spanning tree
MinSpaTre<-spantree(matdist)

op<-par(pty="s")
image(alt,col=brewer.pal(9,"Greys"))
image(veget,col="darkgreen",add=T)
lines(MinSpaTre,xy,lwd=2)
plot(parcPoints,add=T, col="black",bg="blue",pch=21)
plot(parcPoly,add=T, col="yellow",lwd=0.1)
plot(sampPoints,add=T, col="black",bg="red",pch=21)
# #adding the ID of the sampled population
# scatterutil.eti(sampPoints$POINT_X,sampPoints$POINT_Y,
#                 as.character(sampPoints$name_id),0.5)
# #or a simplier way
# text(sampPoints,sampPoints$name_id)
box()
par(op)


###############################################################################
#Canonical Correspondence Analysis
###############################################################################

BRAt<-BRAcc
BRADE<-df2genind(BRAt[,14:27],ncode=3,ind.names=BRAt$sample_ID,
                 pop=BRAt$pop_ID,missing=0,ploidy=1)
BRADE@other$xy<-BRAt[,4:5]
BRADEpop<-genind2genpop(BRADE,process.other=T,missing="0")
cca1<-cca(as.data.frame(BRADEpop$tab),site_table[,-c(1:3)],scann=F)
cca1
plot(cca1)


###############################################################################
#clusterisation with DAPC method using adegenet package
###############################################################################

#the analysis here are performed for the clone-corrected dataset, but the same 
#analysis can be performed with the other dataset by replacing BRAcc by the 
#complete (BRA) or the conservative clone-corrected dataset (BRAcccons) in the 
#following line code, and then rerun other lines of code
BRAt<-BRAcc #name of the input file

#converting data to a genind format
BRADE<-df2genind(BRAt[,14:27],ncode=3,ind.names=BRAt$sample_ID, 
                 pop=BRAt$pop_ID,missing=0,ploidy=1)
BRADE@other$xy<-BRAt[,4:5]
#determination of the number of clusters
clustBRADE<- find.clusters(BRADE,max.n.clust=35)
#with 40 PCs, we lost nearly no information
clustBRADE<- find.clusters(BRADE,n.pca=40,max.n.clust=35) #chose 3 clusters
#which individuals in which clusters per population
table(pop(BRADE),clustBRADE$grp)
#DAPC by itself, first we try to optimized the number of principal component 
#(PCs) to retain to perform the analysis
dapcBRADE<-dapc(BRADE,clustBRADE$grp,n.da=5,n.pca=100)
temp<-optim.a.score(dapcBRADE)
dapcBRADE<-dapc(BRADE,clustBRADE$grp,n.da=5,n.pca=30)
temp<-optim.a.score(dapcBRADE) #based on this result, we finaly chose 7 PCs
dapcBRADE<-dapc(BRADE,clustBRADE$grp,n.da=7,n.pca=7)
#STRUCTURE-like graphic
compoplot(dapcBRADE,lab=NA)
scatter(dapcBRADE,xax=1, yax=2)

BRADEpop<-genind2genpop(BRADE,process.other=T,missing="0")

image(alt,col=brewer.pal(9,"Greys"))
stars(table(pop(BRADE),dapcBRADE$assign),draw.segment=TRUE,
      locations=BRADEpop@other$xy,
      #locations=cbind(jitter(BRADEpop@other$xy$longitude,200),
      #                jitter(BRADEpop@other$xy$latitude,200)),
      add=T,len=0.5)


###############################################################################
#clusterisation with GENELAND package
###############################################################################

MCMC(coordinates=BRAcc[,c(4,5)], geno.hap=BRAcc[,c(14:27)],varnpop=TRUE,
     npopmax=15,spatial=TRUE,freq.model="Uncorrelated", nit=100000,
     thinning=100, path.mcmc="c:/Users/bbarres/GENELAND/BRA/",
     filter.null.alleles=FALSE,delta.coord=0.003,npopinit=2)

PostProcessChain(coordinates=BRAcc[,c(4,5)],
                 path.mcmc="c:/Users/bbarres/GENELAND/BRA/",
                 nxdom=100,nydom=100, burnin=200)

Plotnpop(path.mcmc="c:/Users/bbarres/GENELAND/BRA/",burnin=200)

PlotTessellation(coordinates=BRAcc[,c(4,5)],
                 path.mcmc="c:/Users/bbarres/GENELAND/BRA/")

PosteriorMode(coordinates=BRAcc[,c(4,5)],
              path.mcmc="c:/Users/bbarres/GENELAND/BRA/")


###############################################################################
#Definition of functions to compute diversity indices
###############################################################################

#Allelic Richness computation
#data: a dataset at the 'genind' format from the package 'adegenet'
AllRich<-function(data)
{
  #Conversion from 'genind' object to 'genpop' object
  datapop<-genind2genpop(data, process.other=TRUE, other.action=mean,
                         quiet = TRUE)
  #First, determining the smaller number of allele across sampled population
  matloc<-t(matrix(data=datapop@loc.fac,nrow=(dim(datapop@tab)[2]), 
                   ncol=(dim(datapop@tab)[1])))
  matpop<-matrix(data=datapop@pop.names, nrow=(dim(datapop@tab)[1]), 
                 ncol=(dim(datapop@tab)[2]))
  conf<-list(matpop, matloc)
  effN<-(tapply(datapop@tab, conf, sum))
  effN<- effN[order(as.numeric(rownames(effN))),]
  colnames(effN)<-locNames(datapop)
  echMin<-min(effN)
  
  #Second, build of the matrix of total number of sampled allele 
  truc<-t(as.matrix(table(datapop@loc.fac)))
  x<-matrix(nrow=(dim(effN)[1]), ncol=(dim(effN)[2]), data=truc,byrow=TRUE)
  effTot<-matrix(rep(t(effN),t(x)), nrow=(dim(datapop@tab)[1]), 
                 ncol=(dim(datapop@tab)[2]), byrow=TRUE)
  
  #Third, compute the matrix of Ar for each population/loci combination
  #(see El Mousadik and Petit 1996 for details)
  CoMat<-matrix(nrow=(dim(datapop@tab)[1]),ncol=(dim(datapop@tab)[2]))
  for (i in 1:(dim(datapop@tab)[1])) {
    for (j in 1:(dim(datapop@tab)[2])) {
      CoMat[i,j]<-(1-(nCm(effTot[i,j]-datapop@tab[i,j],echMin)/
        nCm(effTot[i,j],echMin)))
    }
  }
  
  #Allelic richness in each population, for each LOCUS
  ArLOC<-(tapply(CoMat, conf, sum))
  ArLOC<-ArLOC[order(as.numeric(dimnames(ArLOC)[[1]])),]
  colnames(ArLOC)<-locNames(datapop)
  ##determining mean Allelic Richness across site and loci
  #determining mean Allelic Richness across loci
  #Ar<-(apply(ArLOC,1,mean))
  rez<-list("Minimum Sampling Size"=echMin,"Allelic Richness Matrix"=ArLOC)
  return(rez)
}

AllRich(BRADE)[[2]]
Ar<-apply(AllRich(BRADE)[[2]],1,mean)

#Private Allelic Richness computation
#data: a dataset at the 'genind' format from the package 'adegenet'
PrivAllRich<-function(data)
{
  #Conversion from 'genind' object to 'genpop' object
  datapop<-genind2genpop(data, process.other=TRUE,other.action=mean,quiet=TRUE)
  #First, determining the smaller number of allele across sampled population
  matloc<-t(matrix(data=datapop@loc.fac,nrow=(dim(datapop@tab)[2]), 
                   ncol=(dim(datapop@tab)[1])))
  matpop<-matrix(data=datapop@pop.names, nrow=(dim(datapop@tab)[1]), 
                 ncol=(dim(datapop@tab)[2]))
  conf<-list(matpop, matloc)
  effN<-(tapply(datapop@tab, conf, sum))
  effN<- effN[order(as.numeric(rownames(effN))),]
  colnames(effN)<-locNames(datapop)
  echMin<-min(effN)
  
  #Second, build of the matrix of total number of sampled allele 
  truc<-t(as.matrix(table(datapop@loc.fac)))
  x<-matrix(nrow=(dim(effN)[1]), ncol=(dim(effN)[2]), data=truc,byrow=TRUE)
  effTot<-matrix(rep(t(effN),t(x)), nrow=(dim(datapop@tab)[1]), 
                 ncol=(dim(datapop@tab)[2]), byrow=TRUE)
  
  #Third, compute the matrix of Pijg for each population/loci combination
  #(see Kalinowski 2004 for details)
  CoMat<-matrix(nrow=(dim(datapop@tab)[1]),ncol=(dim(datapop@tab)[2]))
  for (i in 1:(dim(datapop@tab)[1])) {
    for (j in 1:(dim(datapop@tab)[2])) {
      CoMat[i,j]<-(1-(nCm(effTot[i,j]-datapop@tab[i,j],echMin)/
                        nCm(effTot[i,j],echMin)))
    }
  }
  #fourth, compute the product of Qijg for each population/loci combination
  #(see Kalinowski 2004 for details)
  CoMat2<-matrix(nrow=(dim(datapop@tab)[1]),ncol=(dim(datapop@tab)[2]))
  for (i in 1:(dim(datapop@tab)[1])) {
    for (j in 1:(dim(datapop@tab)[2])) {
      CoMat2[i,j]<-(nCm(effTot[i,j]-datapop@tab[i,j],echMin)/
                      nCm(effTot[i,j],echMin))
    }
  }
  #fifth, compute the product of Kijg for each population/loci combination
  #(see Kalinowski 2004 for details)
  CoMat3<-matrix(nrow=(dim(datapop@tab)[1]),ncol=(dim(datapop@tab)[2]))
  temp<-c()
  for (i in 1:(dim(datapop@tab)[1])) {
    temp<-as.matrix(CoMat2[-i,])
    ifelse(dim(temp)[2]==1,CoMat3[i,]<-apply(temp,1,prod),
           CoMat3[i,]<-apply(temp,2,prod))
  }
  CoMat4<-CoMat*CoMat3
  #Private Allelic richness in each population, for each LOCUS
  PrivArLOC<-(tapply(CoMat4, conf, sum))
  PrivArLOC<-PrivArLOC[order(as.numeric(dimnames(PrivArLOC)[[1]])),]
  colnames(PrivArLOC)<-locNames(datapop)
  ##determining mean Allelic Richness across site and loci
  #determining mean Allelic Richness across loci
  #Ar<-(apply(ArLOC,1,mean))
  rez<-list("Minimum Sampling Size"=echMin,
            "Private Allelic Richness Matrix"=PrivArLOC)
  return(rez)
}

PrivAllRich(BRADE)
PrivAr<-apply(PrivAllRich(BRADE)[[2]],1,mean)

#Heterozygosity computation
#data: a dataset at the 'genind' format
HeterNei<-function(data)
{
  #Conversion from 'genind' object to 'genpop' object
  datapop<-genind2genpop(data, process.other=TRUE, other.action=mean,
                         missing="0",quiet = TRUE)
  #Heterozygosity (Nei 1987) in each population, for each LOCUS
  HsLOC<-matrix(nrow=(dim(datapop@tab)[1]),
                ncol=(length(data@loc.names)), byrow=TRUE)
  for (i in (1:(dim(datapop@tab)[1]))) {
    dataLOC<-genind2loci(data[data$pop==levels(data$pop)[i]])
    ss<-summary(dataLOC)
    HsLOC[i,]<-sapply(ss, function(x) H(x$allele))
  }
  #determining mean Heterozygosity across loci
  Hs<-(apply(HsLOC,1,mean))
  attr(Hs,"names")<-data@pop.names
  return(Hs)
}

HetNei<-HeterNei(BRADE)

#A function is directly implemented in adegenet: 'Hs'. It gives slightly 
#different results, even if they seem to be correlated. Investigation to be 
#conducted, why results are different (maybe different treatment of missing 
#data when computing allele frequencies)

Hs(BRADE)


###############################################################################
#Comparison of diversity between the two groups of populations
###############################################################################

endemic<-c("NHO","PGA2","PGA1","SOR","SRC2","SRC1")
invaded<-c("BBU","DEN","ROS","DAQ3","DAQ2","DAQ1","RND","PEM2","PEM1")

boxplot(Ar[endemic],Ar[invaded])
wilcox.test(Ar[endemic],Ar[invaded])
t.test(Ar[endemic],Ar[invaded])
boxplot(HetNei[endemic],HetNei[invaded])
wilcox.test(HetNei[endemic],HetNei[invaded])
boxplot(PrivAr[endemic],PrivAr[invaded])
wilcox.test(PrivAr[endemic],PrivAr[invaded])
t.test(PrivAr[endemic],PrivAr[invaded])


###############################################################################
#Comparison of diversity between regions
###############################################################################

#Resampling procedure to balance the different number of populations sampled 
#in the two regions. The aim here is to compare gene diversity between the two 
#groups of populations endemic vs invaded

#we want to compare the allelic richness (Ar), private allelic richness (PAr) 
#and the genetic diversity (He) in the two different regions: "endemic" (close 
#to the natural compartment) vs "invaded" (not in the direct vicinity of the 
#natural compartment). Because there are less "endemic" populations sampled 
#than "invaded" population (6 populations vs 10 populations), we are using a 
#nested-resampling procedure: first we are sampling 6 populations (the minimum 
#number of populations in the 2 different regions), and then we sampled 4 
#individuals (the minimum number of individuals in a population) in each of the 
#sampled populations. Then we compute Ar and He for each resampled datasets at 
#the region level; this procedure was repeated 100 times. A one-way ANOVA 
#followed by a post hoc Tuckey test was used to compare the levels of the 
#three different indices of diversity (Gladieux et al 2010)


resampleDIV<-function(data,grp1,grp2,nbsim) 
  #this function perform a balanced resampling of 6 populations in each defined 
  #group, and within each population it samples 4 individuals to create 
  #sub-dataset; then the function compute the allelic richness, the genetic 
  #diversity and the private allelic richness for each resampled dataset 
  
  #data: a data frame of the same format as "BRA"
  #grp1: first group of populations
  #grp2: second group of populations
  #nbsim: the number of resampling procedure
  
{
  ArDistr<-c()
  PrivArDistr<-c()
  HeDistr<-c()
  SimPop<-c()
  for (i in 1:nbsim) {
    popsamp<-sample(grp1, size=6, replace=FALSE)
    BRAinvad<-data[data$pop_ID %in% (popsamp),]
    resamp<-c()
    for (j in 1:6) {
      listind<-as.character(sample(BRAinvad[BRAinvad$pop_ID == popsamp[j],
                                            "sample_ID"],
                                   size=4, replace=FALSE))
      resamp<-rbind(resamp,BRAinvad[BRAinvad$sample_ID %in% (listind),])
    }
    popsamp<-sample(grp2, size=6, replace=FALSE)
    BRAendem<-data[data$pop_ID %in% (popsamp),]
    for (k in 1:6) {
      listind<-as.character(sample(BRAendem[BRAendem$pop_ID == popsamp[k],
                                            "sample_ID"],
                                   size=4, replace=FALSE))
      resamp<-rbind(resamp,BRAendem[BRAendem$sample_ID %in% (listind),])
    }
    resamp$pop_ID<-as.character(resamp$pop_ID)
    resamp$pop_ID[1:24]<-rep("BRAinvad",24)
    resamp$pop_ID[25:48]<-rep("BRAendem",24)
    resampADE<-df2genind(resamp[,14:27],ncode=3,ind.names=resamp$sample_ID,
                         pop=resamp$pop_ID,missing=0,ploidy=1)
    
    ArDistr<-rbind(ArDistr,apply(AllRich(resampADE)[[2]],1,mean))
    PrivArDistr<-rbind(PrivArDistr,apply(PrivAllRich(resampADE)[[2]],1,mean))
    HeDistr<-rbind(HeDistr,HeterNei(resampADE))
    SimPop<-c(SimPop,resampADE)
  }
  rez<-list("Ar distribution"=ArDistr,"PrivAr distribution"=PrivArDistr,
            "He distribution"=HeDistr,
            "Simulated Pop"=SimPop)
  return(rez)
}


endemic<-c("NHO","PGA2","PGA1","SOR","SRC2","SRC1")
invaded<-c("BBU","DEN","ROS","DAQ3","DAQ2","DAQ1","RND","PEM2","PEM1")

BRAt<-BRAcc #name of the input file
distri<-resampleDIV(BRAt,invaded,endemic,100)

#comparison between Allelic richness of endemic and invaded zones
boxplot(distri[[1]][,1], distri[[1]][,2],
        main="Comparison of Allelic Richness",
        names=c(colnames(distri[[1]])[1],colnames(distri[[1]])[2]))
wilcox.test(distri[[1]][,1], distri[[1]][,2],paired=TRUE)
t.test(distri[[1]][,1], distri[[1]][,2],paired=TRUE)
plot(density(distri[[1]][,1]))
lines(density(distri[[1]][,2]),col="red")

#comparison between Private Allelic richness of endemic and invaded zones
boxplot(distri[[2]][,1], distri[[2]][,2],
        main="Comparison of Private Allelic Richness",
        names=c(colnames(distri[[1]])[1],colnames(distri[[1]])[2]))
wilcox.test(distri[[2]][,1], distri[[2]][,2],paired=TRUE)
t.test(distri[[2]][,1], distri[[2]][,2],paired=TRUE)
plot(density(distri[[2]][,2]),col="red")
lines(density(distri[[2]][,1]))

#comparison between Nei diversity of endemic and invaded zones
boxplot(distri[[3]][,2], distri[[3]][,1],
        main="Comparison of Nei Diversity", 
        names=c(colnames(distri[[3]])[2],colnames(distri[[3]])[1]))
wilcox.test(distri[[3]][,2], distri[[3]][,1],paired=TRUE)
t.test(distri[[3]][,2], distri[[3]][,1],paired=TRUE)
plot(density(distri[[3]][,2]))
lines(density(distri[[3]][,1]),col="red")


###############################################################################
#Connectivity index computation
###############################################################################

#first we built the combined table of coordinates of plantations and sampling 
#sites
xy.utm<-rbind(coordinates(parcPoints.utm),coordinates(sampPoints.utm)[,c(1,2)])
row.names(xy.utm)[577:591]<-c("BBU","DEN","DAQ1","DAQ2","DAQ3","NHO","PEM1",
                              "PEM2","PGA1","PGA2","RND","ROS","SRC1","SRC2",
                              "SOR")
colnames(xy.utm)<-c("Longitude", "Latitude")
matdistb.utm<-vegdist(xy.utm,method="euclidean",diag=TRUE, upper=TRUE)/1000
#recovering the area information for each identified plantations in square km 
#(reason why we divided by 1000000)
plant_area<-(sapply(slot(parcPoly.utm, "polygons"), slot, "area"))/1000000
#we add the area of the sampling sites to the table of area since the surface 
#on which the samples were collected was limited, we fixe this area to 1 ha
plant_area<-c(plant_area,rep(0.01,15))
#'alpha' parameter of Laine and Hanski 2006, here we put 1km (take care that 
#the same unit is used in the Weidist matrix)
alpha<-1
racArea<-sqrt(plant_area) #conversion of area in square meter to the square 
#root of the area in square km
Weidist<-as.matrix(matdistb.utm)
Weidist<-exp(-alpha*Weidist)
connect<-matrix(nrow=(dim(Weidist)[1]),ncol=(dim(Weidist)[2]))
dimnames(connect)<-dimnames(Weidist)
for (i in 1:length(racArea)) {
  connect[,i]<-Weidist[,i]*racArea
}
for (i in 1:length(racArea)) {
  connect[i,i]<-0
}
connect<-apply(connect,2,sum)
connect<-cbind("PATCH_ID"=attr(connect,"names"),"connect"=as.numeric(connect))
connect[577:591][order(attr(connect[577:591],"names"))]
plot(Ar, connect[577:591][order(attr(connect[577:591],"names"))])



###############################################################################
#END
###############################################################################