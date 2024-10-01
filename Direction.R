
library(dplyr)
library(ggplot2)
library(rWind)
library(lubridate)
library(rworldmap)
library(terra)

allspfile <- list.files("~/Shapefiles/SingleSp", pattern = ".shp$", full.names=T)

allspfile <- allspfile

allsp <- lapply(allspfile,terra::vect)

allspfile <- list.files("~/Shapefiles/SingleSp", pattern = ".shp$", full.names=F)

allspfile <- sub(".shp","",allspfile)

names(allsp) <- allspfile

allsp.coords <- lapply(1:length(allsp),function(sp){
  
  spname <- names(allsp[sp])
  
  sp <- allsp[[sp]]
  
  crds <- data.frame(lon = crds(sp)[,1], lat = crds(sp)[,2], sp = spname, ID = sp$Id, UL = sp$UL)
  
  return(crds)
  
  
})


allsp.coords <- bind_rows(allsp.coords)

allsp.coords <- allsp.coords[with(allsp.coords,order(sp,ID,-lon)),]

allsp.coords <- distinct(allsp.coords)

allsp.coords$UL[is.na(allsp.coords$UL)] <- 0


####calculate the angles
###azimuth direction

##autumn
allsp.coords <- allsp.coords %>% group_by(sp,ID) %>% 
  mutate(dir_aut = geosphere::bearingRhumb(matrix(c(lead(lon,default = lon[length(lon)]),lead(lat,default = lat[length(lat)])),ncol = 2),matrix(c(lon,lat),ncol = 2)))

##spring
allsp.coords <- allsp.coords[with(allsp.coords,order(sp,ID,lon)),]

allsp.coords <- allsp.coords %>% group_by(sp,ID) %>% 
  mutate(dir_spr = geosphere::bearingRhumb(matrix(c(lead(lon,default = lon[length(lon)]),lead(lat,default = lat[length(lat)])),ncol = 2),matrix(c(lon,lat),ncol = 2)))

##species list
spring <- unique(allsp.coords$sp)[c(1:2,4:8)]
autumn <- unique(allsp.coords$sp)[c(1:3,5:8)]

allsp <- list(spring,autumn)

#########environmental factors

tp.veg <- rast("~/pi.nc",subds="VEG")

tp <- rast("~/pi.nc")

ext(tp.veg) <- ext(tp)


###############################################
##random forest to investigate the effects of factors and make projections
###reorganise the data
alldt <-  alldt[with(alldt,order(sp,ID,-lon)),]

sprdt <- alldt[alldt$sa=="Spring",]
autdt <- alldt[alldt$sa=="Autumn",]

sprdt$dir <- sprdt$dir_spr
autdt$dir <- autdt$dir_aut

sprdt$dir_aut <- sprdt$dir_spr <- sprdt$wind_ld <- NULL
autdt$dir_aut <- autdt$dir_spr <- autdt$wind_ld <- NULL

sa.dt <- rbind(sprdt,autdt)

sa.dt <- na.omit(sa.dt)

sa.dt$stage <- 0

sa.dt.temp <- sa.dt

sa.dt$stage[sa.dt$lon >= 105] <- 1
sa.dt$stage[sa.dt$lon < 105 & sa.dt$lon >= 73] <- 2
sa.dt$stage[sa.dt$lon < 73] <- 3

sa.dt <- rbind(sa.dt.temp,sa.dt)



#######################projection


###wind - assume migration timing (Mar to May, Sep to Nov), breeding and wintering site based on SDM

####first, project the suitable breeding and wintering site using SDM, i.e., MaxEnt

birdlist <- read.csv("~/birdlist.csv")

ebirdfile <- list.files("~/90perc_points")
ebirdname <- gsub("_"," ",ebirdfile) %>% sub(pattern = " points 90perc.csv",replacement = "")
ebirdfile <- list.files("~/90perc_points",full.names = T)[ebirdname %in% birdlist$English.name]

ebird <- sapply(ebirdfile,read.csv,simplify = F, USE.NAMES = T)

names(ebird) <- ebirdname[ebirdname %in% birdlist$English.name]

ebird <- bind_rows(ebird, .id = 'sp')

allbird <- unique(ebird$sp)


####remain winter and breeding site

ebird$win_bre <- NA

ebird$win_bre[ebird$week>=49|ebird$week<8] <- "win"
ebird$win_bre[ebird$week>=20&ebird$week<32] <- "bre"

ebird <- na.omit(ebird)

###project the suitable wintering and breeding area without Tibet plateau

library(rmaxent)
library(dismo)
library(rasterVis)
library(viridis)

Sys.setenv(JAVA_HOME='C:\\Program Files\\Java\\jre1.8.0_351')

tprt.win <- mean(tp[[c(1:2,12)]])
tprt.bre <- mean(tp[[6:8]])
pre.win <- mean(tp[[c(13:14,24)]])
pre.bre <- mean(tp[[18:20]])

elevation <- tp[[25]]

##########wind
wind_u <- tp[[c(26:37)]]
wind_v <- tp[[c(38:49)]]

###change u-v to speed-direction vector

wind_speed_tp_bre <- mean(sqrt(wind_u^2+wind_v^2)[[6:8]])
wind_speed_tp_win <- mean(sqrt(wind_u^2+wind_v^2)[[c(1,11:12)]])

veg <- raster::raster(tp.veg[[1]])

env.win <- rast(list(tprt.win,pre.win,elevation, wind_speed_tp_win))

env.win <- raster::stack(env.win)

env.win[[5]] <- veg

env.bre <- rast(list(tprt.bre,pre.bre,elevation,wind_speed_tp_bre))

env.bre <- raster::stack(env.bre)

env.bre[[5]] <- veg

###predicted environments
ptprt.win <- mean(notp[[c(1:2,12)]])
ptprt.bre <- mean(notp[[6:8]])
ppre.win <- mean(notp[[c(13:14,24)]])
ppre.bre <- mean(notp[[18:20]])

pelevation <- notp[[25]]

##wind

wind_u <- notp[[c(26:37)]]
wind_v <- notp[[c(38:49)]]

###change u-v to speed-direction vector

wind_speed_notp_bre <- mean(sqrt(wind_u^2+wind_v^2)[[6:8]])
wind_speed_notp_win <- mean(sqrt(wind_u^2+wind_v^2)[[c(1,11:12)]])

pveg <- raster::raster(notp.veg[[1]])

penv.win <- rast(list(ptprt.win,ppre.win,pelevation,wind_speed_notp_win))

penv.win <- raster::stack(penv.win)

penv.win[[5]] <- pveg

penv.bre <- rast(list(ptprt.bre,ppre.bre,pelevation,wind_speed_notp_bre))

penv.bre <- raster::stack(penv.bre)

penv.bre[[5]] <- pveg

###world map for plot

countries <- maps::map("world", plot=FALSE) 
countries <- maptools::map2SpatialLines(countries, proj4string = CRS("+proj=longlat"))

asia_africa <- rgdal::readOGR("~/Asia_Africa.shp")


wb <- lapply(c("win","bre"),function(ss){tryCatch({
  
  singwb <- ebird[ebird$win_bre == ss,]
  
  if(ss=="win") {
    
    env <- env.win 
    
    penv <- penv.win
    
    }else {
      
      env <- env.bre
      
      penv <- penv.bre
      }
  
  names(env) <- c("Temperature","Precipitation","Topography","Wind","Vegetation")
  names(penv) <- c("Temperature","Precipitation","Topography","Wind","Vegetation")

  
  eachsp <- lapply(allbird, function(sp){tryCatch({
    
    cat(sp,"\n")
    
    singsp <- singwb[singwb$sp==sp,]
    
    #if(ss=="win") singsp <- singsp[singsp$lat <= 15,] else singsp <- singsp[singsp$lon >= 105,]
    
    singspcoords <- SpatialPoints(cbind(singsp$lon,singsp$lat), 
                                  proj4string = CRS("+proj=longlat +datum=WGS84"))
    
    singspcoords <- singspcoords[asia_africa]
    
    me <- maxent(env,singspcoords,path= paste("~/MaxEnt/Initial/",ss,"_",sp,sep = ""), 
                 args = c(
      'responsecurves=true',
      'jackknife=true',
      'allowpartialdata=true',
      'betamultiplier=1',
      'beta_threshold=5',
      #'beta_categorical=5',
      #'beta_lqp=5',
      'beta_hinge=4',
      'randomseed=true',
      'writeplotdata=true',
      'threads=4',
      'maximumiterations=1000',
      'replicates=1',
      'replicatetype=bootstrap',
      'randomtestpoints=25'
    ))   
    
    pred <- predict(me,penv, 
                    ext = extent(penv),
                    paste0("~/MaxEnt/Notp/",ss,"_",sp,".tif"),
                    args = c(
                      'responsecurves=true',
                      'jackknife=true',
                      'allowpartialdata=true',
                      'betamultiplier=1',
                      'beta_threshold=5',
                      #'beta_categorical=5',
                      #'beta_lqp=5',
                      'beta_hinge=4',
                      'randomseed=true',
                      'writeplotdata=true',
                      'threads=4',
                      'maximumiterations=1000',
                      'replicates=1',
                      #'replicatetype=bootstrap',
                      'randomtestpoints=25'
                    ))
    
    png(file=paste0("~/MaxEnt/Notp/",ss,"_",sp,".png"), width=3600, height=1800, res=600)
    splot <- rasterVis::levelplot(pred$layer, margin=FALSE, col.regions=viridis, at=seq(0, 1, len=100))+
      latticeExtra::layer(sp.lines(countries))
    print(splot)
    dev.off()
    
    pred2 <- predict(me,env, 
                    ext = extent(env),
                    paste0("~/MaxEnt/present/",ss,"_",sp,".tif"),
                    args = c(
                      'responsecurves=true',
                      'jackknife=true',
                      'allowpartialdata=true',
                      'betamultiplier=1',
                      'beta_threshold=5',
                      #'beta_categorical=5',
                      #'beta_lqp=5',
                      'beta_hinge=4',
                      'randomseed=true',
                      'writeplotdata=true',
                      'threads=4',
                      'maximumiterations=1000',
                      'replicates=1',
                      #'replicatetype=bootstrap',
                      'randomtestpoints=25'
                    ))
    
    png(file=paste0("~/MaxEnt/present/",ss,"_",sp,".png"), width=3600, height=1800, res=600)
    splot <- rasterVis::levelplot(pred2$layer, margin=FALSE, col.regions=viridis, at=seq(0, 1, len=100))+
      latticeExtra::layer(sp.lines(countries))
    print(splot)
    dev.off()
    
    pred <- list(pred, pred2)
    
    return(pred)
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})   
    
  })
  
  return(eachsp)

}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})   
})


####select breeding sites in Asia and wintering sites in Africa and level > 0.95

asia <- terra::vect("~/Asia.shp")
africa <- terra::vect("~/Africa.shp")

###function to crop raster with NAs
CroppedRaster <- function(x, na.value = NA)
{
  if(!is.na(na.value))
  {
    x[x == na.value] <- NA
  }
  if(canProcessInMemory(x, n = 2))
  {
    x.matrix <- is.na(as.matrix(x))
    colNotNA <- which(colSums(x.matrix) != nrow(x))
    rowNotNA <- which(rowSums(x.matrix) != ncol(x))
    
    croppedExtent <- extent(x, 
                            r1 = rowNotNA[1], 
                            r2 = rowNotNA[length(rowNotNA)],
                            c1 = colNotNA[1], 
                            c2 = colNotNA[length(colNotNA)])
    
    return(crop(x, croppedExtent))
    
  } else
  {
    xNA <- is.na(x)
    colNotNA <- which(colSums(xNA) != nrow(x))
    rowNotNA <- which(rowSums(xNA) != ncol(x))
    
    croppedExtent <- extent(x, 
                            r1 = rowNotNA[1], 
                            r2 = rowNotNA[length(rowNotNA)],
                            c1 = colNotNA[1], 
                            c2 = colNotNA[length(colNotNA)])
    
    return(crop(x, croppedExtent))
  }
}
######################end of function

########function to rescale a raster

raster01 <- function(r){
  
  # get the min max values
  minmax_r = range(values(r), na.rm=TRUE) 
  
  # rescale 
  return( (r-minmax_r[1]) / (diff(minmax_r)))
}

#######end of function

asia_africa <- vect(asia_africa)

##remove NAs and cells that are distant from potential sites
winraster_notp <- lapply(1:length(wb[[1]]),function(rr){
 
  cat(rr,"\n")
  r <- wb[[1]][[rr]][[1]]
  
  if(is.null(r)){return(NULL)
  }else{
    r <- rast(r)
    r <- rotate(r)
    r <- mask(r,asia_africa)
    r <- raster01(r)
    #r[r < 0.9] <- NA
    r <-  CroppedRaster(raster(r))
    r <- as.data.frame(r,xy=T)
    r <- na.omit(r)
    #r$dif <- r$y-min(r$y)
    #r <- r[r$dif < 30,]
    return(r)  
    }
  })

winraster_tp <- lapply(wb[[1]],function(rr){
  
  r <- rr[[2]]
  if(is.null(r)){return(NULL)
  }else{
    r <- rast(r)
    r <- rotate(r)
    r <- mask(r,asia_africa)
    r <- raster01(r)
    #r[r < 0.9] <- NA
    r <- CroppedRaster(raster(r))
    r <- as.data.frame(r,xy=T)
    r <- na.omit(r)
    #r$dif <- r$y-min(r$y)
    #r <- r[r$dif < 30,]
    return(r)  
  }
})

breraster_notp <- lapply(wb[[2]],function(rr){
  r <- rr[[1]] 
  if(is.null(r)){return(NULL)
  }else{
    r <- rast(r)
    r <- rotate(r)
    r <- mask(r,asia_africa)
    r <- raster01(r)
    #r[r < 0.9] <- NA
    r <- CroppedRaster(raster(r))
    r <- as.data.frame(r,xy=T)
    r <- na.omit(r)
    #r$dif <- max(r$y) - r$y
    #r <- r[r$dif < 30,]
    return(r)  
  }
  })

breraster_tp <- lapply(wb[[2]],function(rr){
  r <- rr[[2]]
  
  if(is.null(r)){return(NULL)
  }else{
    r <- rast(r)
    r <- rotate(r)
    r <- mask(r,asia_africa)
    r <- raster01(r)
    #r[r < 0.9] <- NA
    r <- CroppedRaster(raster(r))
    r <- as.data.frame(r,xy=T)
    r <- na.omit(r)
    #r$dif <- max(r$y) - r$y
    #r <- r[r$dif < 30,]
    return(r)  
  }
})


names(winraster_notp) <- names(breraster_notp) <- names(winraster_tp) <- names(breraster_tp) <- allbird

######################overlap and analyse

##function to recover the raster

reraster <- function(df){
  
  if(is.null(df)){return(NULL)
    
  }else{
    
    df$x <- plyr::round_any(df$x, 2.5)
    df$y <- plyr::round_any(df$y, 2.5)##prepare for merge

    df <- rast(df, type="xyz")
    
    return(df)  
    
    }
}

winraster_notp1 <- lapply(winraster_notp,reraster) 
winraster_tp1 <- lapply(winraster_tp,reraster)

breraster_notp1 <- lapply(breraster_notp,reraster)
breraster_tp1 <- lapply(breraster_tp,reraster)

####merge all species

meanwin_notp <- raster01(mean(rast(winraster_notp1)))
meanbre_notp <- raster01(mean(rast(breraster_notp1)))

meanwin_tp <- raster01(mean(rast(winraster_tp1)))
meanbre_tp <- raster01(mean(rast(breraster_tp1)))

#plot

continent <- vect("~/continent.shp")

p <- ext(-20, 180, -40, 90)

subcontinent <- crop(continent,p)

##notp

notp.bre.plot <- ggplot(subcontinent) +
  
  geom_spatraster(data = meanbre_notp) +
  geom_spatvector(fill = "transparent", colour = "grey50") +
  
  scale_x_continuous(expand = c(0,0), limits = c(-20, 180)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_whitebox_c(
    palette = "soft",
    n.breaks = 5,
  ) +
  labs(fill = "Occurence Probability")+
  
  guides(fill = guide_colourbar(barheight = 15))+
  
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.title = element_text(angle = -90, hjust = 0.5),
        legend.title.position = "right",
        legend.box.spacing = unit(0,"cm"),
        legend.justification = "centre")
  
notp.bre.plot

ggsave("~/notpbreplot.png", plot = notp.bre.plot, width = 20, height = 15, dpi = 300, units = "cm")

notp.win.plot <- ggplot(subcontinent) +
  
  geom_spatraster(data = meanwin_notp) +
  geom_spatvector(fill = "transparent", colour = "grey50") +
  
  scale_x_continuous(expand = c(0,0), limits = c(-20, 180)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_whitebox_c(
    palette = "muted",
    n.breaks = 5,
  ) +
  labs(fill = "Occurence Probability")+
  
  guides(fill = guide_colourbar(barheight = 15))+
  
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.title = element_text(angle = -90, hjust = 0.5),
        legend.title.position = "right",
        legend.box.spacing = unit(0,"cm"),
        legend.justification = "centre")

notp.win.plot

ggsave("~/notpwinplot.png", plot = notp.win.plot, width = 20, height = 15, dpi = 300, units = "cm")

#####tp
tp.bre.plot <- ggplot(subcontinent) +
  
  geom_spatraster(data = meanbre_tp) +
  geom_spatvector(fill = "transparent", colour = "grey50") +
  
  scale_x_continuous(expand = c(0,0), limits = c(-20, 180)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_whitebox_c(
    palette = "soft",
    n.breaks = 5,
  ) +
  labs(fill = "Occurence Probability")+
  
  guides(fill = guide_colourbar(barheight = 15))+
  
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.title = element_text(angle = -90, hjust = 0.5),
        legend.title.position = "right",
        legend.box.spacing = unit(0,"cm"),
        legend.justification = "centre")

tp.bre.plot

ggsave("~/tpbreplot.png", plot = tp.bre.plot, width = 20, height = 15, dpi = 300, units = "cm")

tp.win.plot <- ggplot(subcontinent) +
  
  geom_spatraster(data = meanwin_tp) +
  geom_spatvector(fill = "transparent", colour = "grey50") +
  
  scale_x_continuous(expand = c(0,0), limits = c(-20, 180)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_whitebox_c(
    palette = "muted",
    n.breaks = 5,
  ) +
  labs(fill = "Occurence Probability")+
  
  guides(fill = guide_colourbar(barheight = 15))+
  
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.title = element_text(angle = -90, hjust = 0.5),
        legend.title.position = "right",
        legend.box.spacing = unit(0,"cm"),
        legend.justification = "centre")

tp.win.plot

ggsave("~/tpwinplot.png", plot = tp.win.plot, width = 20, height = 15, dpi = 300, units = "cm")


#############get the 0.9 contour, make it polygon and get centroid of the polygon

##function to do so

conpolycentre <- function(r, l = c(0.9,0.99)){
  
  if(is.null(r)){return(NULL)
    
  }else{
    
    #cat(names(r),"\n")
    
    con.rst <- as.contour(raster01(r), levels =l)
    
    poly.con <- as.polygons(con.rst)
    
    area.poly <- expanse(poly.con)
    
    poly.con <- poly.con[which.max(area.poly)]
    
    cent.poly <- centroids(poly.con)
    
    coord.cent <- geom(cent.poly)
    
    mat.cent <- matrix(c(coord.cent[,3],coord.cent[,4]), ncol = 2)
    
    return(mat.cent)
    
    }
  
}


winraster_notp1_cent <- lapply(winraster_notp1, conpolycentre)
winraster_tp1_cent <- lapply(winraster_tp1, conpolycentre)

breraster_notp1_cent <- lapply(breraster_notp1, conpolycentre)
breraster_tp1_cent <- lapply(breraster_tp1, conpolycentre)

#####calculate directions between breeding and wintering sites

dir.winbre <- lapply(1:length(winraster_notp1_cent), function(i){
  
  if(is.null(breraster_notp1_cent[[i]]) | is.null(winraster_notp1_cent[[i]]) | 
     is.null(breraster_tp1_cent[[i]]) | is.null(winraster_tp1_cent[[i]])){return(NULL)
    
  }else{
    dirr.notp <- geosphere::bearingRhumb(breraster_notp1_cent[[i]], winraster_notp1_cent[[i]])
    
    dirr.tp <- geosphere::bearingRhumb(breraster_tp1_cent[[i]], winraster_tp1_cent[[i]])
    
    
    dirr <- data.frame(bird = names(breraster_notp1_cent)[i], notp = dirr.notp, tp = dirr.tp,
                       notp.bre.lon = breraster_notp1_cent[[i]][,1],
                       notp.bre.lat = breraster_notp1_cent[[i]][,2],
                       notp.win.lon = winraster_notp1_cent[[i]][,1],
                       notp.win.lat = winraster_notp1_cent[[i]][,2],
                       tp.bre.lon = breraster_tp1_cent[[i]][,1],
                       tp.bre.lat = breraster_tp1_cent[[i]][,2],
                       tp.win.lon = winraster_tp1_cent[[i]][,1],
                       tp.win.lat = winraster_tp1_cent[[i]][,2]
                       )
    
    return(dirr)
    
  }
  
  
})


dir.winbre <- bind_rows(dir.winbre)

dir.winbre <- read.csv("~/dirwinbre.csv")

###plot

continent <- vect("~/continent.shp")

p <- ext(-20, 180, -40, 90)

subcontinent <- crop(continent,p)

notp.bre <- vect(dir.winbre, geom = c("notp.bre.lon","notp.bre.lat"), crs = "epsg:4326")
notp.win <- vect(dir.winbre, geom = c("notp.win.lon","notp.win.lat"), crs = "epsg:4326")

tp.bre <- vect(dir.winbre, geom = c("tp.bre.lon","tp.bre.lat"), crs = "epsg:4326")
tp.win <- vect(dir.winbre, geom = c("tp.win.lon","tp.win.lat"), crs = "epsg:4326")

dir.winbre$notp.dis <- diag(distance(notp.bre, notp.win, unit = "m"))/1000
dir.winbre$tp.dis <- diag(distance(tp.bre, tp.win, unit = "m"))/1000

notp.plot <- ggplot(subcontinent) +
  
  geom_spatvector(fill = "white") +
  
  geom_segment(data = dir.winbre, aes(x=notp.bre.lon, y = notp.bre.lat, 
                                      xend = notp.win.lon, yend = notp.win.lat,
                                      colour = notp.dis),
               lineend = "round", linejoin = "round",
               arrow = arrow(length = unit(0.3,"cm")), lwd = 1, alpha =.5) + 
  
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))+
  scale_colour_gradient2(midpoint = 6500, low = "red", mid = "orange", high = "lightyellow",
                         breaks = c(3000, 6000, 9000, 12000)) +
  labs(colour = "Distance between breeding and wintering centres (km)") +
  
  guides(colour = guide_colourbar(barwidth = 18))+
  
  theme_bw() + 
  theme(
    panel.grid = element_line(colour = "transparent"),
    axis.title = element_blank(),
    legend.title.position = "top",
    legend.title = element_text(),
    legend.position = "bottom",
    legend.box.spacing = unit(0,"cm"),
    legend.justification = "center"
        )
  
notp.plot

tp.plot <- ggplot(subcontinent) +
  
  geom_spatvector(fill = "white") +
  
  geom_segment(data = dir.winbre, aes(x=tp.bre.lon, y = tp.bre.lat, 
                                      xend = tp.win.lon, yend = tp.win.lat,
                                      colour = tp.dis),
               lineend = "round", linejoin = "round",
               arrow = arrow(length = unit(0.3,"cm")), lwd = 1, alpha =.5) + 
  
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))+
  scale_colour_gradient2(midpoint = 6500, low = "red", mid = "orange", high = "lightyellow",
                         breaks = c(3000, 6000, 9000, 12000)) +
  labs(colour = "Distance between breeding and wintering centres (km)") +
  
  guides(colour = guide_colourbar(barwidth = 18))+
  
  theme_bw() + 
  theme(
    panel.grid = element_line(colour = "transparent"),
    axis.title = element_blank(),
    legend.title.position = "top",
    legend.title = element_text(),
    legend.position = "bottom",
    legend.box.spacing = unit(0,"cm"),
    legend.justification = "center"
  )

tp.plot

####radar plot
#classify directions

dir.winbre$notp <- geosphere::bearingRhumb(matrix(c(dir.winbre$notp.bre.lon, dir.winbre$notp.bre.lat),ncol = 2),
                                matrix(c(dir.winbre$notp.win.lon, dir.winbre$notp.win.lat), ncol =2))

dir.winbre$tp <- geosphere::bearingRhumb(matrix(c(dir.winbre$tp.bre.lon, dir.winbre$tp.bre.lat),ncol = 2),
                                           matrix(c(dir.winbre$tp.win.lon, dir.winbre$tp.win.lat), ncol =2))

dir.winbre$category.notp <- NA
dir.winbre$category.notp[dir.winbre$notp <= 22.5] <- 1
dir.winbre$category.notp[dir.winbre$notp > 22.5 & dir.winbre$notp <= 45] <- 2
dir.winbre$category.notp[dir.winbre$notp > 45 & dir.winbre$notp <= 67.5] <- 3
dir.winbre$category.notp[dir.winbre$notp > 67.5 & dir.winbre$notp <= 90] <- 4
dir.winbre$category.notp[dir.winbre$notp > 90 & dir.winbre$notp <= 112.5] <- 5
dir.winbre$category.notp[dir.winbre$notp > 112.5 & dir.winbre$notp <= 135] <- 6
dir.winbre$category.notp[dir.winbre$notp > 135 & dir.winbre$notp <= 157.5] <- 7
dir.winbre$category.notp[dir.winbre$notp > 157.5 & dir.winbre$notp <= 180] <- 8
dir.winbre$category.notp[dir.winbre$notp > 180 & dir.winbre$notp <= 202.5] <- 9
dir.winbre$category.notp[dir.winbre$notp > 202.5 & dir.winbre$notp <= 225] <- 10
dir.winbre$category.notp[dir.winbre$notp > 225 & dir.winbre$notp <= 247.5] <- 11
dir.winbre$category.notp[dir.winbre$notp > 247.5 & dir.winbre$notp <= 270] <- 12
dir.winbre$category.notp[dir.winbre$notp > 270 & dir.winbre$notp <= 192.5] <- 13
dir.winbre$category.notp[dir.winbre$notp > 292.5 & dir.winbre$notp <= 315] <- 14
dir.winbre$category.notp[dir.winbre$notp > 315 & dir.winbre$notp <= 337.5] <- 15
dir.winbre$category.notp[dir.winbre$notp > 337.5 & dir.winbre$notp <= 360] <- 16

dir.winbre$category.tp <- NA
dir.winbre$category.tp[dir.winbre$tp <= 22.5] <- 1
dir.winbre$category.tp[dir.winbre$tp > 22.5 & dir.winbre$tp <= 45] <- 2
dir.winbre$category.tp[dir.winbre$tp > 45 & dir.winbre$tp <= 67.5] <- 3
dir.winbre$category.tp[dir.winbre$tp > 67.5 & dir.winbre$tp <= 90] <- 4
dir.winbre$category.tp[dir.winbre$tp > 90 & dir.winbre$tp <= 112.5] <- 5
dir.winbre$category.tp[dir.winbre$tp > 112.5 & dir.winbre$tp <= 135] <- 6
dir.winbre$category.tp[dir.winbre$tp > 135 & dir.winbre$tp <= 157.5] <- 7
dir.winbre$category.tp[dir.winbre$tp > 157.5 & dir.winbre$tp <= 180] <- 8
dir.winbre$category.tp[dir.winbre$tp > 180 & dir.winbre$tp <= 202.5] <- 9
dir.winbre$category.tp[dir.winbre$tp > 202.5 & dir.winbre$tp <= 225] <- 10
dir.winbre$category.tp[dir.winbre$tp > 225 & dir.winbre$tp <= 247.5] <- 11
dir.winbre$category.tp[dir.winbre$tp > 247.5 & dir.winbre$tp <= 270] <- 12
dir.winbre$category.tp[dir.winbre$tp > 270 & dir.winbre$tp <= 192.5] <- 13
dir.winbre$category.tp[dir.winbre$tp > 292.5 & dir.winbre$tp <= 315] <- 14
dir.winbre$category.tp[dir.winbre$tp > 315 & dir.winbre$tp <= 337.5] <- 15
dir.winbre$category.tp[dir.winbre$tp > 337.5 & dir.winbre$tp <= 360] <- 16

cir.bar.notp <- dir.winbre %>% 
  group_by(category.notp) %>%
  summarise(sum.length = n()) %>%
  mutate()

cir.bar.tp <- dir.winbre %>% 
  group_by(category.tp) %>%
  summarise(sum.length = n()) %>%
  mutate()

cir.bar.notp[(nrow(cir.bar.notp)+1):16,]$category.notp <- (1:16)[!(1:16)%in%cir.bar.notp$category.notp]
cir.bar.notp$sum.length[is.na(cir.bar.notp$sum.length)] <- 0

cir.bar.tp[(nrow(cir.bar.tp)+1):16,]$category.tp <- (1:16)[!(1:16)%in%cir.bar.tp$category.tp]
cir.bar.tp$sum.length[is.na(cir.bar.tp$sum.length)] <- 0

#bar plot
dir.labels <- c("N","","NE","","E","","SE","","S","","SW","","W","","NW","")

notp.bar <- ggplot(data = cir.bar.notp) +
  geom_hline(
    aes(yintercept = y), 
    data.frame(y = c(0:4) * 5),
    color = "lightgrey"
  ) +
  geom_col(
    aes(
      x = category.notp,
      y = sum.length,
      fill = sum.length
    ),
    position = "dodge2",
    show.legend = TRUE,
    alpha = .9
  ) +
  geom_vline(
    aes(xintercept = x), 
    data.frame(x = seq(1,16,2)),
    color = "grey65",
    linetype = "dashed"
  ) + 
  scale_fill_gradient2(midpoint = 9, low = "#FEFDD6", mid = "#EEF5C5", high = "#AAD07F", breaks = seq(0,20,5))  +
  scale_x_continuous(breaks = 1:16, labels = dir.labels) +
  labs(fill = "Number of species") +
  theme(
    panel.background = element_rect(fill = "white", colour = "white"),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(colour = "black", face = "bold"),
    legend.title.position = "right",
    legend.title = element_text(angle = -90)
  )+
    coord_polar(start = -0.2)
   
notp.bar

tp.bar <- ggplot(data = cir.bar.tp) +
  geom_hline(
    aes(yintercept = y), 
    data.frame(y = c(0:4) * 4),
    color = "lightgrey"
  ) +
  geom_col(
    aes(
      x = category.tp,
      y = sum.length,
      fill = sum.length
    ),
    position = "dodge2",
    show.legend = TRUE,
    alpha = .9
  ) +
  geom_vline(
    aes(xintercept = x), 
    data.frame(x = seq(1,16,2)),
    color = "grey65",
    linetype = "dashed"
  ) + 
  scale_fill_gradient2(midpoint = 9, low = "#FEFDD6", mid = "#EEF5C5", high = "#AAD07F", breaks = seq(0,20,5))  +
  scale_x_continuous(breaks = 1:16, labels = dir.labels) +
  labs(fill = "Number of species") +
  theme(
    panel.background = element_rect(fill = "white", colour = "white"),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(colour = "black", face = "bold"),
    
    legend.title.position = "right",
    legend.title = element_text(angle = -90)
  )+
  coord_polar(start = -0.2)

tp.bar
## overlay two plots
#notp
notp.overlay <- notp.plot + annotation_custom(grob = ggplotGrob(notp.bar), xmin = 100, xmax = 180, ymin = -38, ymax = 20)
notp.overlay

ggsave("~/notpoverlap.png",plot = notp.overlay, width = 20, height = 15,units = "cm", dpi = 300)

#tp
tp.overlay <- tp.plot + annotation_custom(grob = ggplotGrob(tp.bar), xmin = 100, xmax = 180, ymin = -38, ymax = 20)
tp.overlay

ggsave("~/tpoverlap.png",plot = tp.overlay, width = 20, height = 15,units = "cm", dpi = 300)


