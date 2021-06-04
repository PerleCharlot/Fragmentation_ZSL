################
# Fragmentation #
#################
#01/08/2018 - 10/08/2018

# ***** FIRST STEP: EXTRACTION SUITABLE MAPS
## DIRECTORY WHERE TO SAVE BIG OUTPUTS home/ucbtpvi/Scratch

Args<-commandArgs()
library(plyr)
library(raster)
library(rgeos)
library(data.table)
library(smoothr)

#setwd("C:/ZSL/Madagascar")
setwd("D:/Perle/R_rep/ZSL/Madagascar")
load('grids_Madagascar1995.RDATA')

#Load and read everything
load('spdatafragm_Madagascar1995_half.RDATA')
speciesprefESA <- fread('speciesprefESA.csv')
mean_min_pop <- read.table('average_population_size.csv.xls', sep = ',', header = TRUE)
mml <- read.table('M_ids_Zonation_params.csv', sep = ',', header = TRUE)

# function to identify removed pixel
pixel_removed <- function( data_polygons,nb_poly,min_patch_size){
  data=c()
  for (poly in 1:nb_poly) {
    colp = as.data.frame(data_polygons[poly])
    names(colp)=poly
    data = append(data, colp) }
  #Identify suitable areas too small according to minimum patch size and extract index of those areas
  list_pixel_removed = c()
  area_poly=matrix(0,nrow=nb_poly,ncol=1)
  for (poly in 1:nb_poly) {
    area_poly[poly] <- sum(unlist(data[poly]))
    if (area_poly [poly]< min_patch_size) {
      list_pixel_removed = append(list_pixel_removed, row.names(data_polygons[[poly]]))}}
  return(list_pixel_removed)
}

idspecies<-Args[5]
idspecies = taxid_list[7626]
idspecies =7626
LU2=merge(LUesa,fgrids, by='seq_grd_300m') #get LU at 300m scale
range=subset(speciesbysitesparse,taxid==idspecies)
#get species at 300 scale, with LU
rangeLU=merge(LU2,range,by.x='seq_grd_30min',by.y='site')#936 698 obs
preferences=subset(speciesprefESA,taxid==idspecies)
if(length(preferences)==0){
  print('Habitat preferencies not found')}
HSM=subset(rangeLU, LC %in% preferences$esa_code)
HSM2=ddply(HSM[,c(2,4,5,6)],.(long,lat,seq_grd_300m),summarize,suitarea=sum(LCArea))#lat,long, area, cellcd for considered species

print('HSM done')

#Set up
rm <- 1
bm <- mml[which(mml$tax_id == idspecies),11]
if(length(bm)==0){
  print('Body mass not found')}
dispersal <- mml[which(mml$tax_id == idspecies),5]*1000 # dispersal in meters
density <- mml[which(mml$tax_id == idspecies),14]

#Calculation of the minimum patch size
minpop <- source("linear_regression_body_mass.R")
min_patch_size <- minpop$value / density
min_stepstone_size <- 10/density #drop areas that can't support at least 10 individuals
area_tresh <- units::set_units(min_stepstone_size, km^2) # drop areas that can't support at least 10 individuals
#Extract coordinates and superficie of suitable areas for the selected species
suitable_areas <- HSM2[,c(1,2,4)]

#Convert suitables areas dataframe to raster
sp_raster <- rasterFromXYZ(suitable_areas,crs="+proj=longlat +datum=WGS84") #change res according to LU map, in km 
#project in Mollweide
sp_raster_proj <- projectRaster(sp_raster, crs="+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m",res=300)
dim_raster <- dim(sp_raster)

#Convert suitable areas raster to polygons and buffer without dropping small stepping stones
sp_polygons1 <- rasterToPolygons(sp_raster_proj, fun=NULL, n=4, na.rm=TRUE, digits=6, dissolve=TRUE)#takes long time
#sp_polygons1 <- rasterToPolygons(sp_raster, fun=NULL, n=4, na.rm=TRUE, digits=6, dissolve=TRUE)#takes long time
sp_polygons1<- disaggregate(sp_polygons1)
number_initial_area <- length(sp_polygons1$suitarea)
areas1<- area(sp_polygons1)/1000000
sp_polygons1$suitarea <- areas1
area_BF <- sum(areas1) # BEFORE FRAGMENTATION IS APPLIED area_BF and area_BSS should be the same. is all suitable areas

# ***** FRAGMENTATION STEP
# buffer areas without removing stepping stones
sp_buffer1 <- buffer(sp_polygons1, width = dispersal)#takes long time
sp_disag1 <- disaggregate(sp_buffer1)
data_polygons1 <- over(sp_disag1, sp_polygons1, returnList = TRUE)#takes time
nb_poly1 <- length(data_polygons1)
list_pixel_removed1 <- pixel_removed(data_polygons1,nb_poly1,min_patch_size)
sp_polygons_AF <- sp_polygons1[- which(row.names(sp_polygons1) %in% list_pixel_removed1),]
area_AF <- sum(area(sp_polygons_AF)/1000000) # AFTER FRAGMENTATION IS APPLIED

# ***** POPULATION THRESHOLD STEP
area_BSS <- sum(sp_polygons1$suitarea)
sp_polygons2 <- drop_crumbs(sp_polygons1, threshold = area_tresh) ## drop stuff smaller than minimum stepping stone size
area_ASS <- sum(sp_polygons2$suitarea)# AFTER STEPPING STONE REMOVAL IS APPLIED

## buffer areas after removing stepping stones
sp_buffer <- buffer(sp_polygons2, width = dispersal)#takes long time
sp_disag <- disaggregate(sp_buffer)
#Connect the new bigger polygons with their summed area
data_polygons2 <- over(sp_disag, sp_polygons2, returnList = TRUE)#takes time
nb_poly2 <- length(data_polygons2)
list_pixel_removed2 <- pixel_removed(data_polygons2,nb_poly2,min_patch_size)
## this is the final polygon
sp_polygons_AFS <- sp_polygons2[- which(row.names(sp_polygons2) %in% list_pixel_removed2),]
area_AFS <- sum(area(sp_polygons_AFS))/1000000 # AFTER FRAGMENTATION IS APPLIED ON TOP OF STEPPING STONE REMOVAL
loss_SS <- area_BSS-area_ASS # loss due to stepping stones removal:before SS (stepping stones) and BF
loss_FS <- area_ASS-area_AFS # loss due to fragmentation after removal of stepping stones
loss_AF <- area_BSS-area_AF # loss due to fragmentation without removal of stepping stones
loss_AFS <- area_BSS - area_AFS # loss due to both stepping stones removal and fragmentation

r <- raster(nrow = dim_raster[1], ncol = dim_raster[2]) 
extent(r) <- extent(sp_polygons_AFS)
sp_final_raster <-rasterize(sp_polygons_AFS, r, 0.09)#2/3 min
#Convert raster to dataframe
sp_final_dataframe <- as.data.frame(sp_final_raster, xy = TRUE)
names(sp_final_dataframe) <- c('long', 'lat', 'suitarea')

print('FRAGMENTATION DONE')

frag_effect <- cbind(idspecies, area_BSS, area_AF,area_ASS, area_AFS, loss_AF,loss_SS,loss_FS, loss_AFS)
print(frag_effect)

write.csv(frag_effect, file = paste(idspecies, 'Mada1992_habitat_loss.csv', sep='_'))
write.csv(new_dataframe, file = paste(idspecies, 'Mada1992_suitarea_frag.csv', sep='_'))

par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(sp_polygons1,border='red',col='red')#initial, BSS. RED
plot(sp_polygons_AF,add = TRUE, border='orange', col='orange')#after fragm only, AF. ORANGE
plot(sp_polygons2,add = TRUE, border='green', col='green')#after SS only, ASS?. GREEN
plot(sp_polygons_AFS,add = TRUE, border='blue', col='blue')#after both, AFS. BLUE
legend("topright", inset=c(-0.2,0),title = paste('Suitable Habitat Map for Species', idspecies, sep =' '),
       legend=c(paste('BFS - lost by 1st fragmentation', round(area_BF), 'km2', sep=' '),
                           paste('AF - stepping stones removed', round(area_ASS), 'km2', sep=' '), 
                           paste('AS - lost by second fragmentation', round(area_AF), 'km2', sep=' '),
                           paste(' Final map', round(area_AFS), 'km2', sep=' ')),
       col=c("red", "orange", "green","blue"),lwd=1:1,cex=0.5, bty="n")
