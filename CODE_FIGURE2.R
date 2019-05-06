
library(rpca)
library(ggplot2)
library(maps)
library(gridExtra)
require(gridExtra)
library(grid)
library(R.matlab)
library(RColorBrewer)
library(scales)

#################################################################################### CROP COVERAGE ###############################################################################
setwd("C:/Users/Ehsan/Desktop/Network/Robast PCA/5-Correlation Maps Mix/NEW/input")
DATA2=readMat("crops_pdsi_pc_rmaping.mat") 
wheat_pdsi_pc_rmaping=DATA2$wheat.pdsi.pc.rmaping

################################################################################### jet.colors #######################################################################

jet.colors <- c("#00007F","#0007E7","#0057FF","#0094ff","#00c3ff","#e9f507","#F5CB07","#FF6200","#E70800","#7F0000");

PiYG <- c("#276419" ,"#4D9221",   "#7FBC41",  "#B8E186", "#E6F5D0" , "#FDE0EF" , "#F1B6DA" , "#DE77AE" , "#C51B7D",  "#8E0152")
####################################################################################################################################################################

sst_max = matrix(  c(0.65), nrow=1,  ncol=1) ;
sst_min = matrix(  c(-0.65), nrow=1,  ncol=1) ;
pdsi_max = matrix(  c(0.75), nrow=1,  ncol=1) ;
pdsi_min = matrix(  c(-0.75), nrow=1,  ncol=1) ;

setwd("C:/Users/Ehsan/Desktop/Network/Robast PCA/2-MIC computation in R and saveas mat/results")
DATA2=readMat("6CROPS_SSTNCDC_COR_rmapping.mat") 

cor_sst_ncdc_wheat=DATA2$cor.sst.ncdc.wheat.rmaping

cor_sst_ncdc_wheat[,2,]=cor_sst_ncdc_wheat[,2,]-1


#########################
##                     ##
##                     ##
##         PCs         ##
##                     ##
##                     ##
#########################     
setwd("C:/Users/Ehsan/Desktop/Network/Robast PCA/1-R-RPCA/input/detrended yield of 6 crops/wheat")
country <- read.table("countries.txt",header = FALSE, sep = "\t");
wheat  <- read.table("sorted_residual.txt",header = FALSE, sep = "\t");
wheat=as.matrix(wheat)
M=t(wheat);
#############################################
library(rpca)
a=rpca(M,lambda = 1/sqrt(max(dim(M))), mu = prod(dim(M))/(4 * sum(abs(M))),
       term.delta = 10^(-7), max.iter = 500000, trace = FALSE,
       thresh.nuclear.fun = thresh.nuclear, thresh.l1.fun = thresh.l1,
       F2norm.fun = F2norm)
S=(as.matrix(a$S));
L=as.matrix(a$L);
##############################################
setwd("C:/Users/Ehsan/Desktop/Network/Robast PCA/1-R-RPCA/results/New/wheat")
loadings=t(a$L.svd$vt);
scores=a$L.svd$u%*%diag(a$L.svd$d)
sdv=as.matrix(a$L.svd$d)
variance=sdv/sum(sdv);
score_frame=as.data.frame(scores)
write.table(loadings, "loadings.txt", sep="\t")
write.table(scores, "scores.txt", sep="\t")
############################################
#high and low loadings and scores based on the mean and standard deviation
my.sdv.loading<-array(NA,dim=c(1,nrow(sdv)))
my.sdv.scores<-array(NA,dim=c(1,nrow(sdv)))
my.mean.loading<-array(NA,dim=c(1,nrow(sdv)))
my.mean.scores<-array(NA,dim=c(1,nrow(sdv)))
for (i in 1:nrow(sdv)) {
  my.sdv.loading[1,i]<-sd(loadings[,i])
  my.sdv.scores[1,i]<-sd(scores[,i])
  my.mean.loading[1,i]<-mean(loadings[,i])
  my.mean.scores[1,i]<-mean(scores[,i])
}
upper.loading<-my.mean.loading+my.sdv.loading
lower.loading<-my.mean.loading-my.sdv.loading
upper.score<-my.mean.scores+my.sdv.scores
lower.score<-my.mean.scores-my.sdv.scores
loading.output<-array(NA,dim=c(ncol(M),nrow(sdv)))
scores.output<-array(NA,dim=c(nrow(M),nrow(sdv)))

for (i in 1:ncol(M)) {
  for (j in 1:nrow(sdv)){
    if (loadings[i,j]>upper.loading[1,j] ) {loading.output[i,j]<-1} 
  }
}
for (i in 1:ncol(M)) {
  for (j in 1:nrow(sdv)){
    if ( loadings[i,j]<lower.loading[1,j] ) {loading.output[i,j]<- -1} 
  }
}
for (i in 1:nrow(M)) {
  for (j in 1:nrow(sdv)){
    if (scores[i,j]>upper.score[1,j]) {scores.output[i,j]<-1} 
  }
}
for (i in 1:nrow(M)) {
  for (j in 1:nrow(sdv)){
    if (scores[i,j]<lower.score[1,j]) {scores.output[i,j]<- -1} 
  }
}
##############################################

new=cbind2(country,loadings)
new=cbind(new, low="NA")
new=cbind(new, high="NA")

score_frame=cbind(score_frame, low="NA")
score_frame=cbind(score_frame, high="NA")
score_frame=cbind(score_frame, time="NA")
########################################################################################################################
##                    
##                                resacling AT and SST in order to have same color schame for all the maps         ##
##                   
########################################################################################################################
n_of_row=dim(wheat_pdsi_pc_rmaping)[1];
n_of_3rd=dim(wheat_pdsi_pc_rmaping)[3];
###############################################
wheat_pdsi_pc_rmaping_rescale_sig=wheat_pdsi_pc_rmaping;
a01<-array(NA,dim=c((n_of_row*n_of_3rd),1))

k=0;
for (i in 1:n_of_3rd) { 
  for (j in 1:n_of_row) {
    k=k+1;
    a01[k,1]=wheat_pdsi_pc_rmaping_rescale_sig[j,3,i]
  }
}
a01[abs(a01) < 0.2858] <- NA
a01 = rescale(a01)+2 
k=0;
for (i in 1:n_of_3rd) { 
  for (j in 1:n_of_row) {
    k=k+1;
    wheat_pdsi_pc_rmaping_rescale_sig[j,3,i]=a01[k,1]
  }
}
##############
wheat_pdsi_pc_rmaping_rescale_all=wheat_pdsi_pc_rmaping;
a02<-array(NA,dim=c((n_of_row*n_of_3rd),1))

k=0;
for (i in 1:n_of_3rd) { 
  for (j in 1:n_of_row) {
    k=k+1;
    a02[k,1]=wheat_pdsi_pc_rmaping_rescale_all[j,3,i]
  }
}


a02=rbind(a02, pdsi_max,pdsi_min)
a02 = rescale(a02)+2 
a02 <- a02[-nrow(a02),] 
a02=as.data.frame(a02)
a02 <- a02[-nrow(a02),] 
a02=as.data.frame(a02)

k=0;
for (i in 1:n_of_3rd) { 
  for (j in 1:n_of_row) {
    k=k+1;
    wheat_pdsi_pc_rmaping_rescale_all[j,3,i]=a02[k,1]
  }
}
#################

n_of_row=dim(cor_sst_ncdc_wheat)[1];
n_of_3rd=dim(cor_sst_ncdc_wheat)[3];
###############################################
cor_sst_ncdc_wheat_rescale_sig=cor_sst_ncdc_wheat;
a03<-array(NA,dim=c((n_of_row*n_of_3rd),1))

k=0;
for (i in 1:n_of_3rd) { 
  for (j in 1:n_of_row) {
    k=k+1;
    a03[k,1]=cor_sst_ncdc_wheat_rescale_sig[j,3,i]
  }
}

a03[abs(a03) < 0.2858] <- NA
a03=rbind(a03, sst_max,sst_min)
a03 = rescale(a03)
a03 <- a03[-nrow(a03),] 
a03=as.data.frame(a03)
a03 <- a03[-nrow(a03),] 
a03=as.data.frame(a03)

k=0;
for (i in 1:n_of_3rd) { 
  for (j in 1:n_of_row) {
    k=k+1;
    cor_sst_ncdc_wheat_rescale_sig[j,3,i]=a03[k,1]
  }
}


#####################################       mapping correlation of AT and SST



country.index.high<-which(loading.output[,1]==1)
countries.high<-as.matrix((country[country.index.high,1]))
countries.high<-t(countries.high)
country.index.low<-which(loading.output[,1]==-1)
countries.low<-as.matrix((country[country.index.low,1]))
countries.low<-t(countries.low)
################################
countries.low_loading<-as.matrix((loadings[country.index.low,1]))
countries.high_loading<-as.matrix((loadings[country.index.high,1]))
low_loading=data.frame(cbind2(as.list(t(countries.low)),as.double(countries.low_loading)))
high_loading=data.frame(cbind2(as.list(t(countries.high),as.double(countries.high_loading))))
################################
new[,"low"] <- NA
new[,"high"] <- NA
#################################
score_frame[,"low"] <- NA
score_frame[,"high"] <- NA
score_frame[,"time"] <- NA
score_frame[,"low"] <- upper.score[1,1]
score_frame[,"high"] <- lower.score[1,1]
score_frame[,"time"] <- 1964:2010
################################
country.index.highh=as.matrix(country.index.high)
country.index.loww=as.matrix(country.index.low)
for (i in 1:nrow(country.index.loww)) {
  new[country.index.loww[i,1],"low"]=as.matrix(new[country.index.loww[i,1],1])
}
for (i in 1:nrow(country.index.highh)) {
  new[country.index.highh[i,1],"high"]=as.matrix(new[country.index.highh[i,1],1])
}
####################################################
a1=wheat_pdsi_pc_rmaping_rescale_sig[,,1];
a1=a1[complete.cases(a1), ];
a1=as.data.frame(a1)
a2=wheat_pdsi_pc_rmaping_rescale_all[,,1];
a2=a2[complete.cases(a2), ];
a2=as.data.frame(a2)
a3=cor_sst_ncdc_wheat_rescale_sig[,,1];
a3=a3[complete.cases(a3), ];
a3=as.data.frame(a3)
legend_title <- "PC1"
map.world <- map_data("world")
gg_sst1 <- ggplot()
gg_sst1 <- gg_sst1 + geom_raster(data=  a2, aes(V2,V1,fill=V3))
#gg_sst1 <- gg_sst1 + geom_raster(data=a3, aes(V2,V1,fill=V3))
gg_sst1 <- gg_sst1 + geom_point(data= a1, aes(V2,V1,fill=V3),shape=15,fill="black", size=0.01)+theme(axis.title.x=element_blank())+theme(axis.title.y=element_blank())
gg_sst1 <- gg_sst1  +scale_fill_gradientn(
  colours=c(jet.colors,PiYG),values = rescale(
    c(rescale(seq(from = min(cor_sst_ncdc_wheat[,3,],na.rm=TRUE), 
                  to = max(cor_sst_ncdc_wheat[,3,],na.rm=TRUE),
                  length.out = length(jet.colors))),
      rescale(seq(from = min(wheat_pdsi_pc_rmaping[,3,],na.rm=TRUE),
                  to = max(wheat_pdsi_pc_rmaping[,3,],na.rm=TRUE),
                  length.out = length(PiYG)))+2)),limits = c(0, 3))
gg_sst1 <- gg_sst1 + geom_map(dat=map.world, map = map.world, aes(map_id=region), fill="NA", color="gray",size=0.25)
gg_sst1 <- gg_sst1 + geom_map(map = map.world, aes(map_id = new[,"V1"]), fill="NA", color="black",size=0.25)
gg_sst1 <- gg_sst1 + geom_map(map = map.world, aes(map_id = new[,"low"]), color="black",size=0.25,alpha=0)
gg_sst1 <- gg_sst1 + geom_map(map = map.world, aes(map_id = new[,"high"]), color="black",size=0.25,alpha=0)
gg_sst1 <- gg_sst1 + expand_limits(x = map.world$long, y = map.world$lat)
gg_sst1 <- gg_sst1 + theme(panel.grid=element_blank(), panel.border=element_blank())
gg_sst1 <- gg_sst1+ coord_equal() + theme(legend.position = "none")
gg_sst1 <- gg_sst1+ coord_equal()  + theme_bw()+theme(legend.position = "none")+theme(axis.title.x=element_blank(),axis.text.x=element_blank())+theme(axis.title.y=element_blank(),axis.text.y=element_blank()) + theme(legend.position = "none")
gg_sst1 <- gg_sst1+ggtitle(paste(round(variance[1,1]*100,2),"% of variability"));
gg_sst1 <- gg_sst1+  theme(plot.title = element_text(hjust = 0.5))
gg_sst1



