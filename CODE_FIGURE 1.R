library(ggplot2)
library(maps)
library(diagram)  
library(plotrix)
library(gridExtra)
library(rJava)
library(maptools)
library(geosphere)
#################################

#####################################
wheat_negative <- read.csv(file="C:/Users/Ehsan/Desktop/Data Incubator/INPUT/wheat_net_negative.csv", header=TRUE, sep=",")
wheat_trade <- read.csv(file="C:/Users/Ehsan/Desktop/Data Incubator/INPUT/wheat_trade.csv", header=TRUE, sep=",")
PDSI <- read.csv(file="C:/Users/Ehsan/Desktop/Data Incubator/INPUT/PDSI_in_croplnads_2012.csv", header=TRUE, sep=",")
#################################################
dev.off()
map.world<-map_data("world")
map.world <- subset(map.world, region!="Antarctica")
gg <- ggplot()
gg <- gg + geom_map(dat=map.world, map = map.world, aes(map_id=region), fill="white", color="black",bg="balck")
#gg <- gg + geom_map(data=wheat_negative,map = map.world, aes(map_id =from),fill = "grey60", colour = "black")
#gg <- gg + geom_map(data=wheat_negative,map = map.world, aes(map_id =to),fill = "grey60", colour = "black")
gg <- gg + expand_limits(x = map.world$long, y = map.world$lat)
gg <- gg + geom_raster(data=PDSI, aes(PDSI$long, PDSI$lat,fill=PDSI$pdsi ),alpha=0.7)
gg <- gg + scale_fill_gradient(low="red", high="yellow",guide = "PDSI")
gg <- gg +  guides(fill = guide_legend(keywidth = 3, keyheight = 7))
gg <- gg + geom_point(data=wheat_trade, aes(wheat_trade$fromlongg, wheat_trade$fromlatt),shape=21,fill="purple", size=3.5)
gg <- gg + geom_point(data=wheat_trade, aes(wheat_trade$tolongg, wheat_trade$tolatt),shape=24,fill="purple", size=3.5)
gg <- gg +theme(panel.background = element_rect(fill = "black"))
gg <- gg + coord_equal()+xlab("Long")+ylab("Lat")
gg <- gg + xlab("Long")+ylab("Lat")
############################################################
col.1 <- adjustcolor("deeppink", alpha=0.4)
col.2 <- adjustcolor("deeppink", alpha=0.4)
edge.pal <- colorRampPalette(c(col.1, col.2), alpha = TRUE)
edge.col <- edge.pal(100)
####################################################################################
for (i in 1:89) {
  edge.ind <- round(100*wheat_negative[i,]$corr / min(wheat_negative$corr))
  p1<-c(wheat_negative[i,7],wheat_negative[i,6])#long lat
  p2<-c(wheat_negative[i,9],wheat_negative[i,8])#ong lat
  
  #################################################################################################################
  #########
  #       #
  #   1   #
  #       #
  #########
  bezier.curve <- function(p1, p2, p3) {
    n <- seq(0,1,length.out=50)
    bx <- (1-n)^2 * p1[[1]] +
      (1-n) * n * 2 * p3[[1]] +
      n^2 * p2[[1]]
    by <- (1-n)^2 * p1[[2]] +
      (1-n) * n * 2 * p3[[2]] +
      n^2 * p2[[2]]
    
    data.frame(lon=bx, lat=by)
  }
  #################################################################################################################
  #########
  #       #
  #   2   #
  #       #
  #########
  bezier.uv.arc <- function(p1, p2) {
    # Get unit vector from P1 to P2
    u <- p2 - p1
    u <- u / sqrt(sum(u*u))
    d <- sqrt(sum((p1-p2)^2))
    # Calculate third point for spline
    m <- d / 2
    h <- floor(d * .2)
    # Create new points in rotated space 
    pp1 <- c(0,0)
    pp2 <- c(d,0)
    pp3 <- c(m, h)
    mx <- as.matrix(bezier.curve(pp1, pp2, pp3))
    # Now translate back to original coordinate space
    theta <- acos(sum(u * c(1,0))) * sign(u[2])
    ct <- cos(theta)
    st <- sin(theta)
    tr <- matrix(c(ct,  -1 * st, st, ct),ncol=2)
    tt <- matrix(rep(p1,nrow(mx)),ncol=2,byrow=TRUE)
    points <- tt + (mx %*% tr)
    tmp.df <- data.frame(points)
    colnames(tmp.df) <- c("lon","lat")
    tmp.df
  }
  ######################################################################################################################
  #########
  #       #
  #   3   #
  #       #
  #########
  bezier.uv.merc.arc <- function(p1,p2) {
    pp1 <- p1
    pp2 <- p2
    pp1[2] <- asinh(tan(p1[2]/180 * pi))/pi * 180
    pp2[2] <- asinh(tan(p2[2]/180 * pi))/pi * 180
    arc <- bezier.uv.arc(pp1,pp2)
    arc$lat <-  atan(sinh(arc$lat/180 * pi))/pi * 180
    arc
  }
  
  arc5 <- bezier.uv.merc.arc(p1,p2)
  gg<-gg+ geom_path(data=as.data.frame(arc5), lwd=1 ,aes(x=lon, y=lat, group=NULL),col=edge.col[edge.ind])
  
  #####################################################################################################################
}

col.1 <- adjustcolor("purple", alpha=0.4)
col.2 <- adjustcolor("purple", alpha=0.4)
edge.pal <- colorRampPalette(c(col.1, col.2), alpha = TRUE)
edge.col <- edge.pal(100)
####################################################################################
for (i in 1:12) {
  edge.ind <- round(100*wheat_trade[i,]$dollar / max(wheat_trade$dollar))
  
  p1<-c(wheat_trade[i,6],wheat_trade[i,5])#long lat
  p2<-c(wheat_trade[i,8],wheat_trade[i,7])#ong lat
  
  #################################################################################################################
  #########
  #       #
  #   1   #
  #       #
  #########
  bezier.curve <- function(p1, p2, p3) {
    n <- seq(0,1,length.out=50)
    bx <- (1-n)^2 * p1[[1]] +
      (1-n) * n * 2 * p3[[1]] +
      n^2 * p2[[1]]
    by <- (1-n)^2 * p1[[2]] +
      (1-n) * n * 2 * p3[[2]] +
      n^2 * p2[[2]]
    
    data.frame(lon=bx, lat=by)
  }
  #################################################################################################################
  #########
  #       #
  #   2   #
  #       #
  #########
  bezier.uv.arc <- function(p1, p2) {
    # Get unit vector from P1 to P2
    u <- p2 - p1
    u <- u / sqrt(sum(u*u))
    d <- sqrt(sum((p1-p2)^2))
    # Calculate third point for spline
    m <- d / 2
    h <- floor(d * .2)
    # Create new points in rotated space 
    pp1 <- c(0,0)
    pp2 <- c(d,0)
    pp3 <- c(m, h)
    mx <- as.matrix(bezier.curve(pp1, pp2, pp3))
    # Now translate back to original coordinate space
    theta <- acos(sum(u * c(1,0))) * sign(u[2])
    ct <- cos(theta)
    st <- sin(theta)
    tr <- matrix(c(ct,  -1 * st, st, ct),ncol=2)
    tt <- matrix(rep(p1,nrow(mx)),ncol=2,byrow=TRUE)
    points <- tt + (mx %*% tr)
    tmp.df <- data.frame(points)
    colnames(tmp.df) <- c("lon","lat")
    tmp.df
  }
  ######################################################################################################################
  #########
  #       #
  #   3   #
  #       #
  #########
  bezier.uv.merc.arc <- function(p1,p2) {
    pp1 <- p1
    pp2 <- p2
    pp1[2] <- asinh(tan(p1[2]/180 * pi))/pi * 180
    pp2[2] <- asinh(tan(p2[2]/180 * pi))/pi * 180
    arc <- bezier.uv.arc(pp1,pp2)
    arc$lat <-  atan(sinh(arc$lat/180 * pi))/pi * 180
    arc
  }
  arc5 <- bezier.uv.merc.arc(p1,p2)
  gg<-gg+ geom_path(data=as.data.frame(arc5), lwd=edge.ind*.03 ,aes(x=lon, y=lat, group=NULL),col="purple")
  #####################################################################################################################
}
gg