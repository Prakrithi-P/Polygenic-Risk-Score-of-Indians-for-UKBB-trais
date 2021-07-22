library(ggplot2)
library(forcats)
library(AICcmodavg)
library(ggpubr)
library(dplyr)

setwd("/home/prakrithi.p/scratch/PRS_cardio_met/panukbb/ukbb_summ_stats/processed")
#read summary stats
input_summary_stats<-list.files("../","*.gz$") 
#count files
n<-length(input_summary_stats)      
#extract each file and sort based on chr_pos
for(i in 1:n){
  system(paste("zcat ../",input_summary_stats[i]," | nawk '$1=($1\"_\"$2 FS )' | { IFS= read -r line; printf \"%s\n\" \"$line\"; sort;} > ",input_summary_stats[i],".txt_sorted", sep=""))
}     #sort
#overlap summary stat SNPs with IGV SNPs
sumstats_sorted<-list.files("./","*sorted$")
for(i in 1:n){
  system(paste("join IGV_SNPs_new ",sumstats_sorted[i]," | awk '/EUR/ || $8<1 {print }' > ",sumstats_sorted[i],"_inIGV", sep=""))
}  
#read files, sort based on pval and extract top 100 SNPs
inIGV<-list.files("./","*.txt_sorted_inIGV$")  
for(i in 1:n){ 
 df<-read.csv(paste("./",inIGV[i],sep=""),sep=" ", header=T)
 df_pval_sorted <-df[order(df$pval_EUR),]
 top100<-head(df_pval_sorted, 100)
 top100$OR<-exp(top100$beta_EUR)
 info<-top100 %>% select (SNP,alt,OR,beta_EUR,pval_EUR,se_EUR,chr_pos,ref)
 write.table(info,paste("./",inIGV[i],"_top100_info",sep=""), sep="\t", row.names = F, quote=F)
}
#calculate PRS with info from summary statistics
top100_info<-list.files("./","*_top100_info_n$") 
for(i in 1:n){
  system(paste("plink --bfile ~/scratch/IGV/IGV_7l --score ",top100_info[i],"  1 2 3 header  --out ~/scratch/PRS_cardio_met/panukbb/ukbb_summ_stats/processed/",top100_info[i],"_PRS",sep=""))
} 
#extract SNPs to flip strand and re-calculate PRS
nopred<-list.files("./","*nopred$") 
for(i in 1:n){
  system(paste("awk '{print $2}' ",nopred[i]," > ",nopred[i],"_to_flip", sep=""))
}
flip<-list.files("./","*flip$") 
for(i in 1:n){
  system(paste("../plink --bfile ../IGV --score ",top100_info[i],"  1 2 3 header --flip ",flip[i]," --out ",top100_info[i],"_PRS",sep=""))
}
#add population labels
profile<-list.files("./","*_PRS.profile$") 
for(i in 1:n){
  system(paste("paste IGV_485_pop ",profile[i]," | awk '{print $1\"\t\"$2\"\t\"$3\"\t\"$4\"\t\"$5\"\t\"$6\"\t\"$7}' > ",top100_info[i],"_PRS",sep=""))
}
#calculate median population PRS
ind_PRS<-list.files("./","*_PRS$") 
for(i in 1:n){
 p<-read.csv(paste("./",ind_PRS[i],sep=""), sep="\t", header=T)
 pp<-p %>%
  group_by(POP) %>%
  summarise(median(SCORE))
 write.table(pp,paste("./",ind_PRS[i],"_popwise_median",sep=""), sep="\t", row.names = F, quote=F) #pop_medianPRS
 #std<-sd(pp$`median(SCORE)`)
 #write.table(std,paste("./",ind_PRS[i],"_popwise_median_stdev",sep=""), sep="\t", row.names = F, quote=F)  #pop_median_Std_dev
 #pdf(paste(ind_PRS[i],"_bp.pdf")) #to save separate boxplots
 #boxplot <- p %>% ggplot(aes(x=fct_reorder(POP,SCORE,.desc = T), y=SCORE,fill="red")) + geom_boxplot()+theme(axis.text.x = element_text(angle=60,vjust = 0.7),plot.title = element_text(hjust = 0.5))+geom_jitter(width=0.2,alpha=0.4)+ggtitle("Distribution of PRS scores across IGV Populations")+xlab("Population")+ylab("PRS")
 #print(boxplot)  #Boxplot
 #dev.off()
}

sink("ANOVA.txt")
an<-list()
for(i in 1:n){
  p<-read.csv(paste("./",ind_PRS[i],sep=""), sep="\t", header=T)
  AN1<-aov(SCORE~POP, p)   #1W ANOVA
  print(trait[i])
  print(summary(AN1))
}
sink()

f<-read.csv("filename_traits",sep="\t", header = T)
trait<-as.character(f$description)
bp<-list()
for(i in 1:n){
  p<-read.csv(paste("./",ind_PRS[i],sep=""), sep="\t", header=T)
  bp[[i]]<-p %>% ggplot(aes(x=fct_reorder(POP,SCORE,.desc = T), y=SCORE,fill="red")) + geom_boxplot()+theme(axis.text.x = element_text(angle=60,vjust = 0.7),plot.title = element_text(hjust = 0.5))+geom_jitter(width=0.2,alpha=0.4)+ggtitle(paste(trait[i]," - Distribution of PRS scores across IGV Populations",sep=""))+xlab("Population")+ylab("PRS")+theme(legend.position = "none")
}
pdf("boxplots_n.pdf")
for (i in 1:n) {
  print(bp[[i]])
}
dev.off()

#rename with Trait name
old<-list.files("./","*_PRS$")
for(i in 1:n){
  system(paste("mv ",old[i],"_popwise_median ",trait[i],"_popwise_median",sep=""))
}

#combine PRS scores and district coordinates 
pop_median<-list.files("./","*_popwise_median$")
for(i in 1:n){
  system(paste("join District_coords ",pop_median[i], " | sed 's/median(SCORE)/PRS/' > ", pop_median[i],"_with_coords",sep=""))
}
#spatial plots
library(rgdal)
library(tmap)
library(maptools)
library(tmap)
library(spatstat)
library(gstat) # Use gstat's idw routine
library(sp)    # Used for the spsample function
library(raster)
library(rgeos)
sink()
lat_long<-list.files("./","*_with_coords$")
sp<-list()
for(i in 1:n){
d<-read.csv(paste(lat_long[i],sep=""),sep=" ", header=TRUE)


###
UTM32n <- CRS("+proj=utm +zone=32 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
# World Geographic System 1984 (lat/long) - mapping
WGS84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84") 
p <-  SpatialPointsDataFrame(coords = d[,c("Longitude", "Latitude")], 
                             data = d, 
                             proj4string = UTM32n)

# Load India boudary map
shp<-readOGR("~/Downloads/India Shape/IGISMAP/Indian_States.shp")

# Replace point boundary extent with that of India
p@bbox <- shp@bbox
th  <-  as(dirichlet(as.ppp(p)), "SpatialPolygons")

# The dirichlet function does not carry over projection information
# requiring that this information be added manually
proj4string(th) <- proj4string(p)

# The tessellated surface does not store attribute information
# from the point data layer. We'll use the over() function (from the sp
# package) to join the point attributes to the tesselated surface via
# a spatial join. The over() function creates a dataframe that will need to
# be added to the `th` object thus creating a SpatialPolygonsDataFrame object
th.z     <- over(th, p, fn=mean)
th.spdf  <-  SpatialPolygonsDataFrame(th, th.z)

# Finally, we'll clip the tessellated  surface to the India boundaries
set_RGEOS_CheckValidity(2L)
th.clp   <- intersect(shp,th.spdf)

grd              <- as.data.frame(spsample(p, "regular", n=50000))
names(grd)       <- c("X", "Y")
coordinates(grd) <- c("X", "Y")
gridded(grd)     <- TRUE  # Create SpatialPixel object
fullgrid(grd)    <- TRUE  # Create SpatialGrid object

# Add P's projection information to the empty grid
proj4string(p) <- proj4string(p) # Temp fix until new proj env is adopted
proj4string(grd) <- proj4string(p)

# Interpolate the grid cells using a power value of 2 (idp=2.0)
P.idw <- gstat::idw(PRS ~ 1, p, newdata=grd, idp=2.0)

# Convert to raster object then clip to Texas
r       <- raster(P.idw)
r.m     <- mask(r, shp)

# Plot
sp[[i]]<-tm_shape(r.m) + 
  tm_raster(n=10,palette = "YlOrRd", auto.palette.mapping = FALSE,
            title=paste(trait[i]," (PRS)",sep="")) + 
  tm_shape(p) + tm_dots(size=0.05) +
  tm_legend(legend.outside=TRUE)+tm_legend(legend.outside=TRUE)+tm_text("POP", just="top", xmod=0, size = 0.5,fontface = 2, overwrite.lines = FALSE, auto.placement = TRUE)
}
#pdf("spatial_plots.pdf")
#for (i in 1:n) {
#  print(bp[[i]])
#  print(sp[[i]])
#}
#dev.off()
pdf("box_spatial_plots.pdf", height=15, width=8)
for (i in 1:n) {
  print(grid.arrange(bp[[i]],tmap_grob(sp[[i]]),nrow=2, ncol = 1))
}
dev.off()
