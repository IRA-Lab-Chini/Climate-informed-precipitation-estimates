# Publication: Douglas C. Jaks (2023), "Climate-informed IDF Curve Updates and Impact on Design Policy at USAF Installations". Thesis: AFIT-ENV-MS-23-M-202
# Note: 1) This script is an example of generating change factor in design precipitation estimates compared to NOAA Atlas 14 estimates
#       and future depth duration frequency (DDF) curves for Tyndall Air Force Base.
#       2) This script can be implemented for any other stations.

setwd("~/Setup Working Directory/")

#. Install required packages from CRAN library

#install.packages("extRemes")
#install.packages("ncdf4")
#install.packages("raster")
#install.packages("rgdal")
library(extRemes)
library(ncdf4)
library(raster)
library(rgdal)

dataset = "~/Setup Working Directory/Extraction_precip.nc"


# GCMs : (53) mpi-esm-lr.1.rcp45, (54) mpi-esm-lr.1.rcp85,
#        (31) gfdl-esm2m.1.rcp45, (32) gfdl-esm2m.1.rcp85
#        (9) canesm2.1.rcp45, (10) canesm2.1.rcp85
#        (39) hadgem2-es.1.rcp45, (40) hadgem2-es.1.rcp85

# Historical: # get only 1/1/1950 to 12/31/1999
t1 = 18264 ; t2 = 36525

# Future: # get only 1/1/2023 to 12/31/2072
# t1 = 44927.00 ; t2 = 63189.00

# # Future: # get only 1/1/2050	to 12/31/2099
# t1 = 54789.00 ; t2 = 73050.00

lat = row #
lon = col #


get.depth<-function(dataset, GCM,lat,lon,t1,t2){
  # bring the netCDF file
data<-nc_open(dataset)
# Save the print(nc) dump to a text file
{
  sink('Extraction_precip_metadata.txt')
  print(data)
  sink()
}

# bring the lat and lon and time elements
#lon <- ncvar_get(data, "lon")
#lat <- ncvar_get(data, "lat", verbose = F)
t <- ncvar_get(data, "Time") 
#cut down the dimensions to the Beltsville station only
#283.0685 , lon deg west
#lon_s<-lon[c(283.0312,283.0938)]
#39.0302, lat deg north
pr.array <- ncvar_get(data,"precip")
dim(pr.array) 
#which(pr.array=="NA")
#nc_close(data) 

#t.r<-round(t,digits=0) - 36522 #get rid of .5 at end of each: convert #days from 1800 to days f/ 1900

#t.sub<-t.r[t.r>=t1 & t.r<=t2] # get only 1/1/950 to 12/31/1999
#t.sub<-t.r[t.r>=44927.00 & t.r<=63189.00] # get only 1/1/2023 to 12/31/2072
#t.sub<-t.r[t.r>=54789.00 & t.r<=73050.00] # get only 1/1/2050	to 12/31/2099

#pr.1<-as.data.frame(pr.array[lat,lon,t1:t2,GCM]) # apply to entire array for the first grid cell
pr.1<-as.data.frame(pr.array[t1:t2,GCM]) # apply to entire array for the first grid cell
colnames(pr.1)<-"Prcp.mm"
pr.1dates<-order(pr.1$Prcp.mm,decreasing = TRUE)
pr.1val<-sort(pr.1$Prcp.mm,decreasing = TRUE)
obs<-as.data.frame(cbind(pr.1dates,pr.1val))
colnames(obs)<-c("date","Val")
obs$i<-0

#establish 7 day independence
for(index in 2:(nrow(obs)-1)){
  UL<-obs$date[index+1]
  LL<-obs$date[index-1]
  test<-obs$date[index]
  if(abs(test - UL) > 7 & abs(test - LL) > 7){
    obs$i[index]<-0
  }else(obs$i[index]<-1)
}

top50.1<-obs$Val[1:50]*0.0393701 # convert to inches
top50<-obs[obs$i==0,]
top50.2<-top50$Val[1:50]*0.0393701 
closeAllConnections()

# assumes the input time step/bin unit is years, search 

GEV.MLE<-fevd(top50.2, type="GEV", units="inches",
              method = "MLE")
#round(mean(a.top50$Prcp.in), digits=2)

MLE.par<-GEV.MLE$results$par
GEV.GMLE<-fevd(top50.2, type="GEV", units="inches",
               method = "GMLE")
GMLE.par<-GEV.GMLE$results$par
GEV.LMOM<-fevd(top50.2, type="GEV", units="inches",
               method = "Lmoments")
LMOM.par<-GEV.LMOM$results

RP<-return.level(GEV.LMOM, c(2,5,10,25,50,100))

return(RP)
}


# get.depth<-function(dataset,GCM,lat,lon,t1,t2)


# dataset = "Extraction_precip.nc"
# GCMs : (53) mpi-esm-lr.1.rcp45, (54) mpi-esm-lr.1.rcp85,
#        (31) gfdl-esm2m.1.rcp45, (32) gfdl-esm2m.1.rcp85
#        (9) canesm2.1.rcp45, (10) canesm2.1.rcp85
#        (39) hadgem2-es.1.rcp45, (40) hadgem2-es.1.rcp85
# Historical: # get only 1/1/1950 to 12/31/1999
# t1 = 18264 ; t2 = 36525
# Future: # get only 1/1/2023 to 12/31/2072
# t1 = 44927.00 ; t2 = 63189.00
# Future: # get only 1/1/2050	to 12/31/2099
# t1 = 54789.00 ; t2 = 73050.00
# lat = row #
# lon = col #
# dile nomenclature: state.GCM.hist/fut1/fut2.RCP45/85.gridcoord(latlon)


summary.table<-data.frame(matrix(nrow=192,ncol=10))
colnames(summary.table)<-c("Index","GCM","RCP","Time Period","2-yr","5-yr","10-yr","25-yr","50-yr","100-yr")
summary.table[,1]<-seq(1,64,1)
summary.table[,2]<-rep(c("access1-0.1.rcp45",
                         "access1-0.1.rcp85",
                         "access1-3.1.rcp45",
                         "access1-3.1.rcp85",
                         "bcc-csm1-1.1.rcp45",
                         "bcc-csm1-1.1.rcp85",
                         "bcc-csm1-1-m.1.rcp45",
                         "bcc-csm1-1-m.1.rcp85",
                         "canesm2.1.rcp45",
                         "canesm2.1.rcp85",
                         "ccsm4.6.rcp45",
                         "ccsm4.6.rcp85",
                         "cesm1-bgc.1.rcp45",
                         "cesm1-bgc.1.rcp85",
                         "cesm1-cam5.1.rcp45",
                         "cesm1-cam5.1.rcp85",
                         "cmcc-cm.1.rcp45",
                         "cmcc-cm.1.rcp85",
                         "cnrm-cm5.1.rcp45",
                         "cnrm-cm5.1.rcp85",
                         "csiro-mk3-6-0.1.rcp45",
                         "csiro-mk3-6-0.1.rcp85",
                         "ec-earth.8.rcp45",
                         "ec-earth.2.rcp85",
                         "fgoals-g2.1.rcp45",
                         "fgoals-g2.1.rcp85",
                         "gfdl-cm3.1.rcp45",
                         "gfdl-cm3.1.rcp85",
                         "gfdl-esm2g.1.rcp45",
                         "gfdl-esm2g.1.rcp85",
                         "gfdl-esm2m.1.rcp45",
                         "gfdl-esm2m.1.rcp85",
                         "giss-e2-r.6.rcp45",
                         "giss-e2-r.2.rcp85",
                         "hadgem2-ao.1.rcp45",
                         "hadgem2-ao.1.rcp85",
                         "hadgem2-cc.1.rcp45",
                         "hadgem2-cc.1.rcp85",
                         "hadgem2-es.1.rcp45",
                         "hadgem2-es.1.rcp85",
                         "inmcm4.1.rcp45",
                         "inmcm4.1.rcp85",
                         "ipsl-cm5a-lr.1.rcp45",
                         "ipsl-cm5a-lr.1.rcp85",
                         "ipsl-cm5a-mr.1.rcp45",
                         "ipsl-cm5a-mr.1.rcp85",
                         "miroc-esm.1.rcp45",
                         "miroc-esm.1.rcp85",
                         "miroc-esm-chem.1.rcp45",
                         "miroc-esm-chem.1.rcp85",
                         "miroc5.1.rcp45",
                         "miroc5.1.rcp85",
                         "mpi-esm-lr.1.rcp45",
                         "mpi-esm-lr.1.rcp85",
                         "mpi-esm-mr.1.rcp45",
                         "mpi-esm-mr.1.rcp85",
                         "mri-cgcm3.1.rcp45",
                         "mri-cgcm3.1.rcp85",
                         "noresm1-m.1.rcp45",
                         "noresm1-m.1.rcp85",
                         "cmcc-cms.1.rcp45",
                         "cmcc-cms.1.rcp85",
                         "giss-e2-h.6.rcp45",
                         "giss-e2-h.2.rcp85"),3)
summary.table[,3]<-c(rep("Historical",64),rep(c("RCP 4.5", "RCP 8.5"),64))
summary.table[,4]<-c(rep("1950-1999",64),rep("2023-2072",64),rep("2050-2099",64))
# Historical: # get only 1/1/1950 to 12/31/1999
# use these row numbers 1:18262
 t1 = 1 ; t2 = 18262
lat=1
lon=1


for(index in 1:64){
    t1 = 1 ; t2 = 18262
    temp<-get.depth("Extraction_precip.nc", summary.table$Index[index], lat,lon,t1,t2)
    summary.table$`2-yr`[index]<-temp[1]
    summary.table$`5-yr`[index]<-temp[2]
    summary.table$`10-yr`[index]<-temp[3]
    summary.table$`25-yr`[index]<-temp[4]
    summary.table$`50-yr`[index]<-temp[5]
    summary.table$`100-yr`[index]<-temp[6]
}

for(index in 65:128){
  t1 = 26664 ; t2 = 44926
  temp<-get.depth("Extraction_precip.nc", summary.table$Index[index], lat,lon,t1,t2)
  summary.table$`2-yr`[index]<-temp[1]
  summary.table$`5-yr`[index]<-temp[2]
  summary.table$`10-yr`[index]<-temp[3]
  summary.table$`25-yr`[index]<-temp[4]
  summary.table$`50-yr`[index]<-temp[5]
  summary.table$`100-yr`[index]<-temp[6]
}

for(index in 129:192){
  t1 = 36526 ; t2 = 54787
  temp<-get.depth("Extraction_precip.nc", summary.table$Index[index], lat,lon,t1,t2)
  summary.table$`2-yr`[index]<-temp[1]
  summary.table$`5-yr`[index]<-temp[2]
  summary.table$`10-yr`[index]<-temp[3]
  summary.table$`25-yr`[index]<-temp[4]
  summary.table$`50-yr`[index]<-temp[5]
  summary.table$`100-yr`[index]<-temp[6]
}



 
 
 matplot(seq(1,24,1),summary.table[4:9],type="o",col=1:6,
         main="STA USC00180700, Beltsville, MD \n DDF Curve \n 24-hr / Dist = EXP",
         xlab="Duration (hrs)",
         ylab="Depth (inches)")
 legend("bottomright", colnames(summary.table)[4:9], col=1:6)
 
 
# calculate change factors
 Change.Factors<-data.frame(matrix(nrow=128,ncol=9))
 colnames(Change.Factors)<-c("GCM","RCP","Time Period","2-yr","5-yr","10-yr","25-yr","50-yr","100-yr")
 Change.Factors[,1]<-summary.table$GCM[65:192]
 Change.Factors[,2]<-rep(c("RCP 4.5", "RCP 8.5"),64)
 Change.Factors[,3]<-c(rep("2023-2072",64),rep("2050-2099",64))
 
 Change.Factors[1:64,4:9]<-c(summary.table[65:128,5:10]/summary.table[1:64,5:10])
 Change.Factors[65:128,4:9]<-c(summary.table[129:192,5:10]/summary.table[1:64,5:10])
   
# RCP 4.5 
 par(mfrow=c(2,2))
 RCP45.CF<-Change.Factors[Change.Factors['RCP']=="RCP 4.5",]
 boxplot(RCP45.CF[1:32,4:9],
         main="STA USC00082391, Dowling Park, FL \n 30.2497N, -83.2594E \n Time Period: 2023-2072 / 24-hr Storm / 32 GCMs / RCP 4.5",
         xlab="Return Period",
         ylab="Change Factor",
         ylim=c(0.5,2.25))
boxplot(RCP45.CF[33:64,4:9],
        main="STA USC00082391, Dowling Park, FL \n 30.2497N, -83.2594E \n Time Period: 2050-2099 / 24-hr Storm / 32 GCMs / RCP 4.5",
        xlab="Return Period",
        ylab="Change Factor",
        ylim=c(0.5,2.25))
 
 # RCP 8.5 
 RCP85.CF<-Change.Factors[Change.Factors['RCP']=="RCP 8.5",]
 boxplot(RCP85.CF[1:32,4:9],
         main="STA USC00082391, Dowling Park, FL \n 30.2497N, -83.2594E \n Time Period: 2023-2072 / 24-hr Storm / 32 GCMs / RCP 8.5",
         xlab="Return Period",
         ylab="Change Factor",
         ylim=c(0.5,2.25))
 boxplot(RCP85.CF[33:64,4:9],
         main="STA USC00082391, Dowling Park, FL \n 30.2497N, -83.2594E \n Time Period: 2050-2099 / 24-hr Storm / 32 GCMs / RCP 8.5",
         xlab="Return Period",
         ylab="Change Factor",
         ylim=c(0.5,2.25))
 
 
 
 
 ################################################################################
 ################################################################################
 ############################# NOAA ATLAS 14 Information#########################
 ################################################################################
 ################################################################################
 NOAA.pds<-read.csv("PF_Depth_English_PDS (1).csv", sep=",",header=TRUE)
 NOAA.pds<-NOAA.pds[1:10,]
 Durations<-c(5/60,10/60,0.25,0.5,1,2,3,6,12,24)
 
 
 par(mfrow=c(1,1))
 plot(Durations, NOAA.pds$X100*quantile(RCP85.CF$`100-yr`[33:64], c(0.05,0.95))[2], type="o", col="cyan", pch="*", lty=1,
      main="STA USC00180700, Beltsville, MD \n 39.0302N -76.9315E",
      xlab="Duration (hours)",
      ylab="Depth (inches)")
 points(Durations, NOAA.pds$X100*median(RCP85.CF$`100-yr`[33:64]),col="blue", pch="*")
 lines(Durations, NOAA.pds$X100*median(RCP85.CF$`100-yr`[33:64]),col="blue", lty=2,lwd=3)
 points(Durations, NOAA.pds$X100*quantile(RCP85.CF$`100-yr`[33:64], c(0.05,0.95))[1],col="cyan", pch="*")
 lines(Durations, NOAA.pds$X100*quantile(RCP85.CF$`100-yr`[33:64], c(0.05,0.95))[1],col="cyan", lty=3) 
 
 points(Durations, NOAA.pds$X100*quantile(RCP85.CF$`100-yr`[33:64], c(0.125,0.875))[2],col="light blue", pch="*")
 lines(Durations, NOAA.pds$X100*quantile(RCP85.CF$`100-yr`[33:64], c(0.125,0.875))[2],col="light blue", lty=4) 
 points(Durations, NOAA.pds$X100*quantile(RCP85.CF$`100-yr`[33:64], c(0.125,0.875))[1],col="light blue", pch="*")
 lines(Durations, NOAA.pds$X100*quantile(RCP85.CF$`100-yr`[33:64], c(0.125,0.875))[1],col="light blue", lty=5) 
 
 points(Durations, NOAA.pds$UL100, col="dark red", pch="*")
 lines(Durations, NOAA.pds$UL100, col="dark red",lty=6)
 points(Durations, NOAA.pds$X100, col="red", pch="*")
 lines(Durations, NOAA.pds$X100, col="red",lty=7,lwd=3)
 points(Durations, NOAA.pds$LL100, col="dark red",pch="*")
 lines(Durations, NOAA.pds$LL100, col="dark red", lty=8)
 text(19,5,"Return Period: 100-yr \n Time Period: 2050-2099 \n RCP 8.5", col="black")
 
 legend(15.5,4,legend=c("Projected","90% CI - Projected","75% CI - Projected","ATLAS 14 Value","90% CI - ATLAS 14"),
        col=c("blue","cyan","light blue","red","dark red"),
        pch="*",lty=c(2,1,4,7,6), ncol=1)
 
 
 ################################################################################
 ################################################################################
 ####################### depth vs duration IDF curve ############################
 ################################################################################
 ################################################################################
 
 plot(Durations,NOAA.pds$X100*median(RCP85.CF$`100-yr`[33:64]),ylim=range(c(NOAA.pds$X100*quantile(RCP85.CF$`100-yr`[33:64], c(0.05,0.95))[1],
                                                                            NOAA.pds$X100*quantile(RCP85.CF$`100-yr`[33:64], c(0.05,0.95))[2])),lty=1,
      main="STA USC00180700, Beltsville, MD \n 39.0302N -76.9315E",
      xlab="Duration (hours)",
      ylab="Depth (inches)"
 )
 #make polygon where coordinates start with lower limit and then upper limit in reverse order
 with(NOAA.pds,polygon(c(Durations,rev(Durations)),c(NOAA.pds$X100*quantile(RCP85.CF$`100-yr`[33:64], c(0.05,0.95))[1],
                                                     rev(NOAA.pds$X100*quantile(RCP85.CF$`100-yr`[33:64], c(0.05,0.95))[2])),col = "grey95", border = FALSE))
 lines(Durations, NOAA.pds$X100*quantile(RCP85.CF$`100-yr`[33:64], c(0.05,0.95))[2],col="red", lty=2)
 lines(Durations, NOAA.pds$X100*quantile(RCP85.CF$`100-yr`[33:64], c(0.05,0.95))[1],col="red", lty=2)
 lines(Durations, NOAA.pds$X100*median(RCP85.CF$`100-yr`[33:64]),col="black", lty=3,lwd=2)
 lines(Durations, NOAA.pds$X100,col="blue", lty=4,lwd=2)
 text(18,5,"Return Period: 100-yr \n Time Period: 2050-2099 \n RCP 8.5", col="black")
 
 legend(15.5,4,legend=c("Projected","90% CI - Projected","ATLAS 14 Value"),
        col=c("black","red","blue"),lty=c(1,2,4), ncol=1)
 
 
 
 ################################################################################
 ################################################################################
 ############################## depth vs RP IDF curve ###########################
 ################################################################################
 ################################################################################
 # subset the NOAA alues to only the averages
 NOAA.avgs<-NOAA.pds[,c(1:9)]
 # now melt the data to put a RP column
 library(reshape2)
 NOAA.melted<-melt(NOAA.avgs,id="by.duration.for.ARI..years..")
 NOAA.melted$Return.Period<-c(rep(1,10),rep(2,10),rep(5,10),rep(10,10),rep(25,10),
                              rep(50,10),rep(100,10),rep(200,10))
 RCP85.depthvRP<-NOAA.melted[11:70,]
 RCP85.depthvRP$CF.depth<-0
 RCP85.depthvRP$CF.depth[1:10]<-RCP85.depthvRP$value[1:10]*median(RCP85.CF$`2-yr`[33:64])
 RCP85.depthvRP$CF.depth[11:20]<-RCP85.depthvRP$value[11:20]*median(RCP85.CF$`5-yr`[33:64])
 RCP85.depthvRP$CF.depth[21:30]<-RCP85.depthvRP$value[21:30]*median(RCP85.CF$`10-yr`[33:64])
 RCP85.depthvRP$CF.depth[31:40]<-RCP85.depthvRP$value[31:40]*median(RCP85.CF$`25-yr`[33:64])
 RCP85.depthvRP$CF.depth[41:50]<-RCP85.depthvRP$value[41:50]*median(RCP85.CF$`50-yr`[33:64])
 RCP85.depthvRP$CF.depth[51:60]<-RCP85.depthvRP$value[51:60]*median(RCP85.CF$`100-yr`[33:64])
 
 
 grapher_IDF_DvRP<-data.frame(matrix(nrow=6,ncol=31))
 colnames(grapher_IDF_DvRP)<-c("Return Period",
                               "1_NOAA",
                               "1_45_2372",
                               "1_85_2372",
                               "1_45_5099",
                               "1_85_5099",
                               "2_NOAA",
                               "2_45_2372",
                               "2_85_2372",
                               "2_45_5099",
                               "2_85_5099",
                               "3_NOAA",
                               "3_45_2372",
                               "3_85_2372",
                               "3_45_5099",
                               "3_85_5099",
                               "6_NOAA",
                               "6_45_2372",
                               "6_85_2372",
                               "6_45_5099",
                               "6_85_5099",
                               "12_NOAA",
                               "12_45_2372",
                               "12_85_2372",
                               "12_45_5099",
                               "12_85_5099",
                               "24_NOAA",
                               "24_45_2372",
                               "24_85_2372",
                               "24_45_5099",
                               "24_85_5099"
 )
 grapher_IDF_DvRP[,1]<-c(2,5,10,25,50,100)
 grapher_IDF_DvRP[,2:31]<-c(RCP85.depthvRP$value[c(5,15,25,35,45,55)],
                            RCP85.depthvRP$value[c(5,15,25,35,45,55)]*median(RCP45.CF$`2-yr`[1:32]),
                            RCP85.depthvRP$value[c(5,15,25,35,45,55)]*median(RCP85.CF$`2-yr`[1:32]),
                            RCP85.depthvRP$value[c(5,15,25,35,45,55)]*median(RCP45.CF$`2-yr`[33:64]),
                            RCP85.depthvRP$value[c(5,15,25,35,45,55)]*median(RCP85.CF$`2-yr`[33:64]),
                            RCP85.depthvRP$value[c(6,16,26,36,46,56)],
                            RCP85.depthvRP$value[c(6,16,26,36,46,56)]*median(RCP45.CF$`5-yr`[1:32]),
                            RCP85.depthvRP$value[c(6,16,26,36,46,56)]*median(RCP85.CF$`5-yr`[1:32]),
                            RCP85.depthvRP$value[c(6,16,26,36,46,56)]*median(RCP45.CF$`5-yr`[33:64]),
                            RCP85.depthvRP$value[c(6,16,26,36,46,56)]*median(RCP85.CF$`5-yr`[33:64]),
                            RCP85.depthvRP$value[c(7,17,27,37,47,57)],
                            RCP85.depthvRP$value[c(7,17,27,37,47,57)]*median(RCP45.CF$`5-yr`[1:32]),
                            RCP85.depthvRP$value[c(7,17,27,37,47,57)]*median(RCP85.CF$`5-yr`[1:32]),
                            RCP85.depthvRP$value[c(7,17,27,37,47,57)]*median(RCP45.CF$`5-yr`[33:64]),
                            RCP85.depthvRP$value[c(7,17,27,37,47,57)]*median(RCP85.CF$`5-yr`[33:64]),
                            RCP85.depthvRP$value[c(8,18,28,38,48,58)],
                            RCP85.depthvRP$value[c(8,18,28,38,48,58)]*median(RCP45.CF$`5-yr`[1:32]),
                            RCP85.depthvRP$value[c(8,18,28,38,48,58)]*median(RCP85.CF$`5-yr`[1:32]),
                            RCP85.depthvRP$value[c(8,18,28,38,48,58)]*median(RCP45.CF$`5-yr`[33:64]),
                            RCP85.depthvRP$value[c(8,18,28,38,48,58)]*median(RCP85.CF$`5-yr`[33:64]),
                            RCP85.depthvRP$value[c(9,19,29,39,49,59)],
                            RCP85.depthvRP$value[c(9,19,29,39,49,59)]*median(RCP45.CF$`5-yr`[1:32]),
                            RCP85.depthvRP$value[c(9,19,29,39,49,59)]*median(RCP85.CF$`5-yr`[1:32]),
                            RCP85.depthvRP$value[c(9,19,29,39,49,59)]*median(RCP45.CF$`5-yr`[33:64]),
                            RCP85.depthvRP$value[c(9,19,29,39,49,59)]*median(RCP85.CF$`5-yr`[33:64]),
                            RCP85.depthvRP$value[c(10,20,30,40,50,60)],
                            RCP85.depthvRP$value[c(10,20,30,40,50,60)]*median(RCP45.CF$`5-yr`[1:32]),
                            RCP85.depthvRP$value[c(10,20,30,40,50,60)]*median(RCP85.CF$`5-yr`[1:32]),
                            RCP85.depthvRP$value[c(10,20,30,40,50,60)]*median(RCP45.CF$`5-yr`[33:64]),
                            RCP85.depthvRP$value[c(10,20,30,40,50,60)]*median(RCP85.CF$`5-yr`[33:64])
                            
 )
 
 
 
 
 
 #  write.csv(grapher_IDF_DvRP, file="florida_grapher_IDFdepthvRP.csv" )
 
 
 plot(RCP85.depthvRP$Return.Period[c(10,20,30,40,50,60)],RCP85.depthvRP$CF.depth[c(10,20,30,40,50,60)],      main="STA USC00180700, Beltsville, MD \n 39.0302N -76.9315E",
      xlab="Return Period (years)",
      ylab="Depth (inches)",log="xy",ylim=c(0.1,round(max(RCP85.depthvRP$CF.depth)))
 )
 lines(RCP85.depthvRP$Return.Period[c(10,20,30,40,50,60)],RCP85.depthvRP$CF.depth[c(10,20,30,40,50,60)])
 lines(RCP85.depthvRP$Return.Period[c(1,11,21,31,41,51)],RCP85.depthvRP$CF.depth[c(1,11,21,31,41,51)])
 lines(RCP85.depthvRP$Return.Period[c(2,12,22,32,42,52)],RCP85.depthvRP$CF.depth[c(2,12,22,32,42,52)])
 lines(RCP85.depthvRP$Return.Period[c(3,13,23,33,43,53)],RCP85.depthvRP$CF.depth[c(3,13,23,33,43,53)]) 
 lines(RCP85.depthvRP$Return.Period[c(4,14,24,34,44,54)],RCP85.depthvRP$CF.depth[c(4,14,24,34,44,54)])
 lines(RCP85.depthvRP$Return.Period[c(5,15,25,35,45,55)],RCP85.depthvRP$CF.depth[c(5,15,25,35,45,55)])
 lines(RCP85.depthvRP$Return.Period[c(6,16,26,36,46,56)],RCP85.depthvRP$CF.depth[c(6,16,26,36,46,56)])
 lines(RCP85.depthvRP$Return.Period[c(7,17,27,37,47,57)],RCP85.depthvRP$CF.depth[c(7,17,27,37,47,57)])
 lines(RCP85.depthvRP$Return.Period[c(8,18,28,38,48,58)],RCP85.depthvRP$CF.depth[c(8,18,28,38,48,58)])
 lines(RCP85.depthvRP$Return.Period[c(9,19,29,39,49,59)],RCP85.depthvRP$CF.depth[c(9,19,29,39,49,59)])
 
 # looking onnly at the 24 hour storm event
 plot(RCP85.depthvRP$Return.Period[c(10,20,30,40,50,60)],RCP85.depthvRP$CF.depth[c(10,20,30,40,50,60)],      main="STA USC00180700, Beltsville, MD \n 39.0302N -76.9315E",
      xlab="Return Period (years)",
      ylab="Depth (inches)",log="xy",ylim=c(2,round(max(RCP85.depthvRP$CF.depth)))
      
 )
 lines(RCP85.depthvRP$Return.Period[c(10,20,30,40,50,60)],RCP85.depthvRP$CF.depth[c(10,20,30,40,50,60)],col="black",lty=2,lwd=3)
 lines(RCP85.depthvRP$Return.Period[c(10,20,30,40,50,60)],RCP85.depthvRP$value[c(10,20,30,40,50,60)],col="blue",lty=4,lwd=3)
 text(25,4,"Duration: 24-hr \n Time Period: 2050-2099 \n RCP 8.5", col="black")
 
 legend(15,3.5,legend=c("Projected Value","ATLAS 14 Value"),
        col=c("black","blue"),lty=c(2,4),lwd=3, ncol=1)
 
 
 #make polygon where coordinates start with lower limit and then upper limit in reverse order
 with(NOAA.pds,polygon(c(Durations,rev(Durations)),c(NOAA.pds$X100*quantile(RCP85.CF$`100-yr`[33:64], c(0.05,0.95))[1],
                                                     rev(NOAA.pds$X100*quantile(RCP85.CF$`100-yr`[33:64], c(0.05,0.95))[2])),col = "grey95", border = FALSE))
 lines(Durations, NOAA.pds$X100*quantile(RCP85.CF$`100-yr`[33:64], c(0.05,0.95))[2],col="red", lty=2)
 lines(Durations, NOAA.pds$X100*quantile(RCP85.CF$`100-yr`[33:64], c(0.05,0.95))[1],col="red", lty=2)
 lines(Durations, NOAA.pds$X100*median(RCP85.CF$`100-yr`[33:64]),col="black", lty=3,lwd=2)
 lines(Durations, NOAA.pds$X100,col="blue", lty=4,lwd=2)
 text(18,5,"Return Period: 100-yr \n Time Period: 2050-2099 \n RCP 8.5", col="black")
 
 legend(15.5,4,legend=c("Projected","90% CI - Projected","ATLAS 14 Value"),
        col=c("black","red","blue"),lty=c(1,2,4), ncol=1)
 
 
 
 
 
 ################################################################################
 ################################################################################
 # create files for Grapher using the IQRs
 grapher_boxplots<-data.frame(matrix(nrow=2,ncol=72))
 colnames(grapher_boxplots)<-c("x2-yr/RCP 4.5/2023-2072","y2-yr/RCP 4.5/2023-2072","m2-yr/RCP 4.5/2023-2072",
                               "x2-yr/RCP 8.5/2023-2072","y2-yr/RCP 8.5/2023-2072","m2-yr/RCP 8.5/2023-2072",
                               "x2-yr/RCP 4.5/2050-2099","y2-yr/RCP 4.5/2050-2099","m2-yr/RCP 4.5/2050-2099",
                               "x2-yr/RCP 8.5/2050-2099","y2-yr/RCP 8.5/2050-2099","m2-yr/RCP 8.5/2050-2099",
                               "x5-yr/RCP 4.5/2023-2072","y5-yr/RCP 4.5/2023-2072","m5-yr/RCP 4.5/2023-2072",
                               "x5-yr/RCP 8.5/2023-2072","y5-yr/RCP 8.5/2023-2072","m5-yr/RCP 8.5/2023-2072",
                               "x5-yr/RCP 4.5/2050-2099","y5-yr/RCP 4.5/2050-2099","m5-yr/RCP 4.5/2050-2099",
                               "x5-yr/RCP 8.5/2050-2099","y5-yr/RCP 8.5/2050-2099","m5-yr/RCP 8.5/2050-2099",
                               "x10-yr/RCP 4.5/2023-2072","y10-yr/RCP 4.5/2023-2072","m10-yr/RCP 4.5/2023-2072",
                               "x10-yr/RCP 8.5/2023-2072","y10-yr/RCP 8.5/2023-2072","m10-yr/RCP 8.5/2023-2072",
                               "x10-yr/RCP 4.5/2050-2099","y10-yr/RCP 4.5/2050-2099","m10-yr/RCP 4.5/2050-2099",
                               "x10-yr/RCP 8.5/2050-2099","y10-yr/RCP 8.5/2050-2099","m10-yr/RCP 8.5/2050-2099",
                               "x25-yr/RCP 4.5/2023-2072","y25-yr/RCP 4.5/2023-2072","m25-yr/RCP 4.5/2023-2072",
                               "x25-yr/RCP 8.5/2023-2072","y25-yr/RCP 8.5/2023-2072","m25-yr/RCP 8.5/2023-2072",
                               "x25-yr/RCP 4.5/2050-2099","y25-yr/RCP 4.5/2050-2099","m25-yr/RCP 4.5/2050-2099",
                               "x25-yr/RCP 8.5/2050-2099","y25-yr/RCP 8.5/2050-2099","m25-yr/RCP 8.5/2050-2099",
                               "x50-yr/RCP 4.5/2023-2072","y50-yr/RCP 4.5/2023-2072","m50-yr/RCP 4.5/2023-2072",
                               "x50-yr/RCP 8.5/2023-2072","y50-yr/RCP 8.5/2023-2072","m50-yr/RCP 8.5/2023-2072",
                               "x50-yr/RCP 4.5/2050-2099","y50-yr/RCP 4.5/2050-2099","m50-yr/RCP 4.5/2050-2099",
                               "x50-yr/RCP 8.5/2050-2099","y50-yr/RCP 8.5/2050-2099","m50-yr/RCP 8.5/2050-2099",
                               "x100-yr/RCP 4.5/2023-2072","y100-yr/RCP 4.5/2023-2072","m100-yr/RCP 4.5/2023-2072",
                               "x100-yr/RCP 8.5/2023-2072","y100-yr/RCP 8.5/2023-2072","m100-yr/RCP 8.5/2023-2072",
                               "x100-yr/RCP 4.5/2050-2099","y100-yr/RCP 4.5/2050-2099","m100-yr/RCP 4.5/2050-2099",
                               "x100-yr/RCP 8.5/2050-2099","y100-yr/RCP 8.5/2050-2099","m100-yr/RCP 8.5/2050-2099")
 # grapher_boxplots[1,]<-0
 #grapher_boxplots[2,]<-0
 
 #  RCP 4.5/2023-2072
 grapher_boxplots[1,
                  c(2,14,26,38,50,62)
 ]<-c(quantile(RCP45.CF$`2-yr`[1:32], c(0.25,0.75))[1],
      quantile(RCP45.CF$`5-yr`[1:32], c(0.25,0.75))[1],
      quantile(RCP45.CF$`10-yr`[1:32], c(0.25,0.75))[1],
      quantile(RCP45.CF$`25-yr`[1:32], c(0.25,0.75))[1],
      quantile(RCP45.CF$`50-yr`[1:32], c(0.25,0.75))[1],
      quantile(RCP45.CF$`100-yr`[1:32], c(0.25,0.75))[1])
 grapher_boxplots[2,
                  c(2,14,26,38,50,62)
 ]<-c(quantile(RCP45.CF$`2-yr`[1:32], c(0.25,0.75))[2],
      quantile(RCP45.CF$`5-yr`[1:32], c(0.25,0.75))[2],
      quantile(RCP45.CF$`10-yr`[1:32], c(0.25,0.75))[2],
      quantile(RCP45.CF$`25-yr`[1:32], c(0.25,0.75))[2],
      quantile(RCP45.CF$`50-yr`[1:32], c(0.25,0.75))[2],
      quantile(RCP45.CF$`100-yr`[1:32], c(0.25,0.75))[2])
 grapher_boxplots[1,
                  c(3,15,27,39,51,63)
 ]<-c(median(RCP45.CF$`2-yr`[1:32]),
      median(RCP45.CF$`5-yr`[1:32]),
      median(RCP45.CF$`10-yr`[1:32]),
      median(RCP45.CF$`25-yr`[1:32]),
      median(RCP45.CF$`50-yr`[1:32]),
      median(RCP45.CF$`100-yr`[1:32]))
 
 #   RCP 8.5/2023-2072
 grapher_boxplots[1,
                  c(5,17,29,41,53,65)
 ]<-c(quantile(RCP85.CF$`2-yr`[1:32], c(0.25,0.75))[1],
      quantile(RCP85.CF$`5-yr`[1:32], c(0.25,0.75))[1],
      quantile(RCP85.CF$`10-yr`[1:32], c(0.25,0.75))[1],
      quantile(RCP85.CF$`25-yr`[1:32], c(0.25,0.75))[1],
      quantile(RCP85.CF$`50-yr`[1:32], c(0.25,0.75))[1],
      quantile(RCP85.CF$`100-yr`[1:32], c(0.25,0.75))[1])
 grapher_boxplots[2,
                  c(5,17,29,41,53,65)
 ]<-c(quantile(RCP85.CF$`2-yr`[1:32], c(0.25,0.75))[2],
      quantile(RCP85.CF$`5-yr`[1:32], c(0.25,0.75))[2],
      quantile(RCP85.CF$`10-yr`[1:32], c(0.25,0.75))[2],
      quantile(RCP85.CF$`25-yr`[1:32], c(0.25,0.75))[2],
      quantile(RCP85.CF$`50-yr`[1:32], c(0.25,0.75))[2],
      quantile(RCP85.CF$`100-yr`[1:32], c(0.25,0.75))[2])
 grapher_boxplots[1,
                  c(6,18,30,42,54,66)
 ]<-c(median(RCP85.CF$`2-yr`[1:32]),
      median(RCP85.CF$`5-yr`[1:32]),
      median(RCP85.CF$`10-yr`[1:32]),
      median(RCP85.CF$`25-yr`[1:32]),
      median(RCP85.CF$`50-yr`[1:32]),
      median(RCP85.CF$`100-yr`[1:32]))
 
 #   RCP 4.5/2050-2099
 grapher_boxplots[1,
                  c(8,20,32,44,56,68)
 ]<-c(quantile(RCP45.CF$`2-yr`[33:64], c(0.25,0.75))[1],
      quantile(RCP45.CF$`5-yr`[33:64], c(0.25,0.75))[1],
      quantile(RCP45.CF$`10-yr`[33:64], c(0.25,0.75))[1],
      quantile(RCP45.CF$`25-yr`[33:64], c(0.25,0.75))[1],
      quantile(RCP45.CF$`50-yr`[33:64], c(0.25,0.75))[1],
      quantile(RCP45.CF$`100-yr`[33:64], c(0.25,0.75))[1])
 grapher_boxplots[2,
                  c(8,20,32,44,56,68)
 ]<-c(quantile(RCP45.CF$`2-yr`[33:64], c(0.25,0.75))[2],
      quantile(RCP45.CF$`5-yr`[33:64], c(0.25,0.75))[2],
      quantile(RCP45.CF$`10-yr`[33:64], c(0.25,0.75))[2],
      quantile(RCP45.CF$`25-yr`[33:64], c(0.25,0.75))[2],
      quantile(RCP45.CF$`50-yr`[33:64], c(0.25,0.75))[2],
      quantile(RCP45.CF$`100-yr`[33:64], c(0.25,0.75))[2])
 grapher_boxplots[1,
                  c(9,21,33,45,57,69)
 ]<-c(median(RCP45.CF$`2-yr`[33:64]),
      median(RCP45.CF$`5-yr`[33:64]),
      median(RCP45.CF$`10-yr`[33:64]),
      median(RCP45.CF$`25-yr`[33:64]),
      median(RCP45.CF$`50-yr`[33:64]),
      median(RCP45.CF$`100-yr`[33:64]))
 
 #   RCP 8.5/2050-2099
 grapher_boxplots[1,
                  c(11,23,35,47,59,71)
 ]<-c(quantile(RCP85.CF$`2-yr`[33:64], c(0.25,0.75))[1],
      quantile(RCP85.CF$`5-yr`[33:64], c(0.25,0.75))[1],
      quantile(RCP85.CF$`10-yr`[33:64], c(0.25,0.75))[1],
      quantile(RCP85.CF$`25-yr`[33:64], c(0.25,0.75))[1],
      quantile(RCP85.CF$`50-yr`[33:64], c(0.25,0.75))[1],
      quantile(RCP85.CF$`100-yr`[33:64], c(0.25,0.75))[1])
 grapher_boxplots[2,
                  c(11,23,35,47,59,71)
 ]<-c(quantile(RCP85.CF$`2-yr`[33:64], c(0.25,0.75))[2],
      quantile(RCP85.CF$`5-yr`[33:64], c(0.25,0.75))[2],
      quantile(RCP85.CF$`10-yr`[33:64], c(0.25,0.75))[2],
      quantile(RCP85.CF$`25-yr`[33:64], c(0.25,0.75))[2],
      quantile(RCP85.CF$`50-yr`[33:64], c(0.25,0.75))[2],
      quantile(RCP85.CF$`100-yr`[33:64], c(0.25,0.75))[2])
 grapher_boxplots[1,
                  c(12,24,36,48,60,72)
 ]<-c(median(RCP85.CF$`2-yr`[33:64]),
      median(RCP85.CF$`5-yr`[33:64]),
      median(RCP85.CF$`10-yr`[33:64]),
      median(RCP85.CF$`25-yr`[33:64]),
      median(RCP85.CF$`50-yr`[33:64]),
      median(RCP85.CF$`100-yr`[33:64]))
 
 grapher_boxplots[1,
                  c(1,13,25,37,49,61)]<-c(0.94,1.34,1.74,2.14,2.54,2.94)
 grapher_boxplots[1,
                  c(4,16,28,40,52,64)]<-c(0.98,1.38,1.78,2.18,2.58,2.98)
 grapher_boxplots[1,
                  c(7,19,31,43,55,67)]<-c(1.02,1.42,1.82,2.22,2.62,3.02)
 grapher_boxplots[1,
                  c(10,22,34,46,58,70)]<-c(1.06,1.46,1.86,2.26,2.66,3.06)
 grapher_boxplots[2,
                  c(1,13,25,37,49,61)]<-c(0.94,1.34,1.74,2.14,2.54,2.94)
 grapher_boxplots[2,
                  c(4,16,28,40,52,64)]<-c(0.98,1.38,1.78,2.18,2.58,2.98)
 grapher_boxplots[2,
                  c(7,19,31,43,55,67)]<-c(1.02,1.42,1.82,2.22,2.62,3.02)
 grapher_boxplots[2,
                  c(10,22,34,46,58,70)]<-c(1.06,1.46,1.86,2.26,2.66,3.06) 
 
 #write.csv(grapher_boxplots, file="florida_grapher_boxplots3.csv" )
 
 ################################################################################
 ################################################################################
 # create files for Grapher for IDFs
 
 
 #DURATIONS == c(5/60,10/60,0.25,0.5,1,2,3,6,12,24)
 grapher_IDF<-data.frame(matrix(nrow=10,ncol=79))
 colnames(grapher_IDF)<-c("xDuration",
                          "2_NOAA",
                          "2m_45_2372","2L_45_2372","2H_45_2372",
                          "2m_85_2372","2L_85_2372","2H_85_2372",
                          "2m_45_5099","2L_45_5099","2H_45_5099",
                          "2m_85_5099","2L_85_5099","2H_85_5099",
                          "5_NOAA",
                          "5m_45_2372","5L_45_2372","5H_45_2372",
                          "5m_85_2372","5L_85_2372","5H_85_2372",
                          "5m_45_5099","5L_45_5099","5H_45_5099",
                          "5m_85_5099","5L_85_5099","5H_85_5099",
                          "10_NOAA",
                          "10m_45_2372","10L_45_2372","10H_45_2372",
                          "10m_85_2372","10L_85_2372","10H_85_2372",
                          "10m_45_5099","10L_45_5099","10H_45_5099",
                          "10m_85_5099","10L_85_5099","10H_85_5099",
                          "25_NOAA",
                          "25m_45_2372","25L_45_2372","25H_45_2372",
                          "25m_85_2372","25L_85_2372","25H_85_2372",
                          "25m_45_5099","25L_45_5099","25H_45_5099",
                          "25m_85_5099","25L_85_5099","25H_85_5099",
                          "50_NOAA",
                          "50m_45_2372","50L_45_2372","50H_45_2372",
                          "50m_85_2372","50L_85_2372","50H_85_2372",
                          "50m_45_5099","50L_45_5099","50H_45_5099",
                          "50m_85_5099","50L_85_5099","50H_85_5099",
                          "100_NOAA",
                          "100m_45_2372","100L_45_2372","100H_45_2372",
                          "100m_85_2372","100L_85_2372","100H_85_2372",
                          "100m_45_5099","100L_45_5099","100H_45_5099",
                          "100m_85_5099","100L_85_5099","100H_85_5099"
 )
 grapher_IDF[,1]<-Durations
 grapher_IDF[,2:79]<-
   c(NOAA.pds$X2,
     NOAA.pds$X2*median(RCP45.CF$`2-yr`[1:32]),
     NOAA.pds$X2*quantile(RCP45.CF$`2-yr`[1:32], c(0.05,0.95))[1],
     NOAA.pds$X2*quantile(RCP45.CF$`2-yr`[1:32], c(0.05,0.95))[2],
     NOAA.pds$X2*median(RCP85.CF$`2-yr`[1:32]),
     NOAA.pds$X2*quantile(RCP85.CF$`2-yr`[1:32], c(0.05,0.95))[1],
     NOAA.pds$X2*quantile(RCP85.CF$`2-yr`[1:32], c(0.05,0.95))[2],
     NOAA.pds$X2*median(RCP45.CF$`2-yr`[33:64]),
     NOAA.pds$X2*quantile(RCP45.CF$`2-yr`[33:64], c(0.05,0.95))[1],
     NOAA.pds$X2*quantile(RCP45.CF$`2-yr`[33:64], c(0.05,0.95))[2],
     NOAA.pds$X2*median(RCP85.CF$`2-yr`[33:64]),
     NOAA.pds$X2*quantile(RCP85.CF$`2-yr`[33:64], c(0.05,0.95))[1],
     NOAA.pds$X2*quantile(RCP85.CF$`2-yr`[33:64], c(0.05,0.95))[2],
     NOAA.pds$X5,
     NOAA.pds$X5*median(RCP45.CF$`5-yr`[1:32]),
     NOAA.pds$X5*quantile(RCP45.CF$`5-yr`[1:32], c(0.05,0.95))[1],
     NOAA.pds$X5*quantile(RCP45.CF$`5-yr`[1:32], c(0.05,0.95))[2],
     NOAA.pds$X5*median(RCP85.CF$`5-yr`[1:32]),
     NOAA.pds$X5*quantile(RCP85.CF$`5-yr`[1:32], c(0.05,0.95))[1],
     NOAA.pds$X5*quantile(RCP85.CF$`5-yr`[1:32], c(0.05,0.95))[2],
     NOAA.pds$X5*median(RCP45.CF$`5-yr`[33:64]),
     NOAA.pds$X5*quantile(RCP45.CF$`5-yr`[33:64], c(0.05,0.95))[1],
     NOAA.pds$X5*quantile(RCP45.CF$`5-yr`[33:64], c(0.05,0.95))[2],
     NOAA.pds$X5*median(RCP85.CF$`5-yr`[33:64]),
     NOAA.pds$X5*quantile(RCP85.CF$`5-yr`[33:64], c(0.05,0.95))[1],
     NOAA.pds$X5*quantile(RCP85.CF$`5-yr`[33:64], c(0.05,0.95))[2],
     NOAA.pds$X10,
     NOAA.pds$X10*median(RCP45.CF$`10-yr`[1:32]),
     NOAA.pds$X10*quantile(RCP45.CF$`10-yr`[1:32], c(0.05,0.95))[1],
     NOAA.pds$X10*quantile(RCP45.CF$`10-yr`[1:32], c(0.05,0.95))[2],
     NOAA.pds$X10*median(RCP85.CF$`10-yr`[1:32]),
     NOAA.pds$X10*quantile(RCP85.CF$`10-yr`[1:32], c(0.05,0.95))[1],
     NOAA.pds$X10*quantile(RCP85.CF$`10-yr`[1:32], c(0.05,0.95))[2],
     NOAA.pds$X10*median(RCP45.CF$`10-yr`[33:64]),
     NOAA.pds$X10*quantile(RCP45.CF$`10-yr`[33:64], c(0.05,0.95))[1],
     NOAA.pds$X10*quantile(RCP45.CF$`10-yr`[33:64], c(0.05,0.95))[2],
     NOAA.pds$X10*median(RCP85.CF$`10-yr`[33:64]),
     NOAA.pds$X10*quantile(RCP85.CF$`10-yr`[33:64], c(0.05,0.95))[1],
     NOAA.pds$X10*quantile(RCP85.CF$`10-yr`[33:64], c(0.05,0.95))[2],
     NOAA.pds$X25,
     NOAA.pds$X25*median(RCP45.CF$`25-yr`[1:32]),
     NOAA.pds$X25*quantile(RCP45.CF$`25-yr`[1:32], c(0.05,0.95))[1],
     NOAA.pds$X25*quantile(RCP45.CF$`25-yr`[1:32], c(0.05,0.95))[2],
     NOAA.pds$X25*median(RCP85.CF$`25-yr`[1:32]),
     NOAA.pds$X25*quantile(RCP85.CF$`25-yr`[1:32], c(0.05,0.95))[1],
     NOAA.pds$X25*quantile(RCP85.CF$`25-yr`[1:32], c(0.05,0.95))[2],
     NOAA.pds$X25*median(RCP45.CF$`25-yr`[33:64]),
     NOAA.pds$X25*quantile(RCP45.CF$`25-yr`[33:64], c(0.05,0.95))[1],
     NOAA.pds$X25*quantile(RCP45.CF$`25-yr`[33:64], c(0.05,0.95))[2],
     NOAA.pds$X25*median(RCP85.CF$`25-yr`[33:64]),
     NOAA.pds$X25*quantile(RCP85.CF$`25-yr`[33:64], c(0.05,0.95))[1],
     NOAA.pds$X25*quantile(RCP85.CF$`25-yr`[33:64], c(0.05,0.95))[2],
     NOAA.pds$X50,
     NOAA.pds$X50*median(RCP45.CF$`50-yr`[1:32]),
     NOAA.pds$X50*quantile(RCP45.CF$`50-yr`[1:32], c(0.05,0.95))[1],
     NOAA.pds$X50*quantile(RCP45.CF$`50-yr`[1:32], c(0.05,0.95))[2],
     NOAA.pds$X50*median(RCP85.CF$`50-yr`[1:32]),
     NOAA.pds$X50*quantile(RCP85.CF$`50-yr`[1:32], c(0.05,0.95))[1],
     NOAA.pds$X50*quantile(RCP85.CF$`50-yr`[1:32], c(0.05,0.95))[2],
     NOAA.pds$X50*median(RCP45.CF$`50-yr`[33:64]),
     NOAA.pds$X50*quantile(RCP45.CF$`50-yr`[33:64], c(0.05,0.95))[1],
     NOAA.pds$X50*quantile(RCP45.CF$`50-yr`[33:64], c(0.05,0.95))[2],
     NOAA.pds$X50*median(RCP85.CF$`50-yr`[33:64]),
     NOAA.pds$X50*quantile(RCP85.CF$`50-yr`[33:64], c(0.05,0.95))[1],
     NOAA.pds$X50*quantile(RCP85.CF$`50-yr`[33:64], c(0.05,0.95))[2],
     NOAA.pds$X100,
     NOAA.pds$X100*median(RCP45.CF$`100-yr`[1:32]),
     NOAA.pds$X100*quantile(RCP45.CF$`100-yr`[1:32], c(0.05,0.95))[1],
     NOAA.pds$X100*quantile(RCP45.CF$`100-yr`[1:32], c(0.05,0.95))[2],
     NOAA.pds$X100*median(RCP85.CF$`100-yr`[1:32]),
     NOAA.pds$X100*quantile(RCP85.CF$`100-yr`[1:32], c(0.05,0.95))[1],
     NOAA.pds$X100*quantile(RCP85.CF$`100-yr`[1:32], c(0.05,0.95))[2],
     NOAA.pds$X100*median(RCP45.CF$`100-yr`[33:64]),
     NOAA.pds$X100*quantile(RCP45.CF$`100-yr`[33:64], c(0.05,0.95))[1],
     NOAA.pds$X100*quantile(RCP45.CF$`100-yr`[33:64], c(0.05,0.95))[2],
     NOAA.pds$X100*median(RCP85.CF$`100-yr`[33:64]),
     NOAA.pds$X100*quantile(RCP85.CF$`100-yr`[33:64], c(0.05,0.95))[1],
     NOAA.pds$X100*quantile(RCP85.CF$`100-yr`[33:64], c(0.05,0.95))[2]
   )
 
 #write.csv(grapher_IDF, file="Florida_grapher_IDF.csv" )
 