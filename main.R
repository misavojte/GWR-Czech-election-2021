# Package managment ====
# spwgr ... GWR
# rgdal, sf ... spatial data manipulation
# ggplot2, corrplot ... plotting
# MASS ... stepwise regression

if (!require("pacman")) install.packages("pacman")
pacman::p_load("rgdal","ggplot2", "spgwr", "MASS", "sf", "corrplot")



# Import spatial data ====
#read shp data (obceGen.shp) from directory (/data)
#parameters encoding & use_iconv for correct encoding of Czech language
data <- readOGR(dsn="data",layer="obceGen", encoding = "UTF-8", use_iconv = T)

#data <- data[data@data$OBY_CEL > 500,]



# Simplified EDA ====
#create boxplot from CSSD gain and save outliers etc.
bp <- boxplot(data@data$CSSD)

#define reverse in operator
`%!in%` <- Negate(`%in%`)

#create data2 without outliers (not used)
data2 <- data[data@data$CSSD %!in% bp$out,]

#create corrplot (with outliers)
cr <- cor(data@data[,c(13:34)],method = "spearman")
corrplot(cr)

#show significant collinearity (|R| > .8)
subset(as.data.frame.table(cr), Freq < 1 & abs(Freq) > 0.8)

#remove variables
subIndep <- subset(data@data[,13:34], select=-c(MUZI, VEK65))

#corrplot 2
corrplot(cor(subIndep))



# Regression (not spatial) ====
#create base reg. model
model <- lm(data@data$CSSD ~ ., data = subIndep)
io <- lm(data@data$CSSD ~ 1, data= subIndep)

#create stepwise regression (bidirectional - forward & backward selection of variables)
resModel <- step(io, direction = "both", scope=formula(model), trace=0)

#compare base & stepwise
anova(model)
summary(model)

anova(resModel)
summary(resModel)

#prepare data for plotting coefs
bpdata <- data.frame(promenna = names(resModel$coefficients), koeficient = resModel$coefficients)
bpdata$promenna <- factor(bpdata$promenna, levels = bpdata$promenna)

#plot coefs
ggplot(bpdata, aes(x=promenna, y=koeficient)) +geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust=1))

#plot model parameters (heteroscedasticity...)
plot(resModel)

#just checking some characteristics
data.frame(data@data$NEZAMEST, data@data$ZAMESTEL, data@data$ZAMESTNA)

#adding residuals and predicted values from regression
data@data$residsCSSD <- residuals(resModel)
data@data$predictedCSSD <- resModel$fitted.values



# Making maps in R ====
#transforming spatial data for future plotting (UTM 33N for Czech Republic)
stDat <- st_transform(st_as_sf(data), st_crs(32633))

#creating polygon of Czech Rep. from municipal data
#it will be used for having gray color in NoData areas
#could be imported instead...
uni <- st_union(stDat)
tx <- st_as_text(uni) #to near-WKT format
ex <- sub("(.*\\(\\((.*?)\\).*)", replacement = '\\2', x = tx) #removing all interior ring from near-WKT format
uni <- st_as_sfc(paste0("POLYGON ((",ex,"))"), crs = st_crs(32633)) #making a proper WKT data and converting them to spatial object (utm 33n)

#plot resids (global) (as a map)
ggplot() + geom_sf(data = uni, fill = "grey50", colour = "grey50") + geom_sf(data = stDat, aes(fill = residsCSSD), colour = NA) + 
  scale_fill_gradient2(low = "#a50026", 
                    mid = "white", 
                    high = "#006837",
                    name = "Rezidua (ČSSD)")

#plot vysledek socani (as a map)
ggplot() + geom_sf(data = uni, fill = "grey50", colour = "grey50") + geom_sf(data = stDat, aes(fill = CSSD * 100), colour = NA) + 
  scale_fill_gradient2(low = "#a50026", 
                       mid = "white", 
                       high = "#006837",
                       name = "Volební zisk [%]")



# Geographically Weighted Regression ====
#see: https://gdsl-ul.github.io/san/gwr.html 

# create optimal ADAPTIVE bandwidth (based on cross-validation)
# this can take a lot of time (MacBook Air 2020 = 10 mins)
GWRbandwidth <- gwr.sel(formula(resModel), data = data@data, coords = coordinates(data), adapt = T)

# create GWR model
# using stepwise reg. model as base and previously calculated bandwidth
# it's adaptive = taking in account proportion of k-nearest neighbours
# this can take even more time (9.42)
GWRmodel <- gwr(formula(resModel), data = data@data, coords = coordinates(data), adapt = GWRbandwidth, hatmatrix=TRUE, se.fit=TRUE) 

#Time could be significantly reduced if using fixed bandwidth!

#see how big the bandwidth was
#plot bandwidths (in meters)
plot(GWRmodel$bandwidth)
#used bandwidths (in meters)
summary(GWRmodel$bandwidth)


#create function for plotting a map from GWR data
plotKoef <- function(koef) {
  stDat[,paste0("k",koef)] <- GWRmodel$SDF[[koef]]
  print(ggplot() + geom_sf(data = uni, fill = "grey50", colour = "grey50") + geom_sf(data = stDat, aes(fill = stDat[,paste0("k",koef)][[1]]), colour = NA) + 
    scale_fill_gradient2(low = "#a50026", 
                         mid = "white", 
                         high = "#006837",
                         name = koef))
}

#testing the function
plotKoef("V_VYSOKO")

#create and save GWR plots to /output directory
#and write data to spatial object
for (i in names(GWRmodel$SDF)) {
  stDat[,paste0("k",i)] <- GWRmodel$SDF[[i]]
  png(paste0("output/plots/",i,".png"))
  plotKoef(i)
  dev.off()
}

#add resids data to spatial object
stDat$residsGWR <- GWRmodel$SDF$gwr.e

changeRes <- abs(GWRmodel$SDF$gwr.e) - abs(residuals(resModel))
stDat$changeRes <- changeRes

#export shp with resids, coeffs etc.
st_write(stDat, dsn="output", driver="ESRI Shapefile")



#Now just plotting some maps in R Studio ====
#create GWR resids map
ggplot() + geom_sf(data = uni, fill = "grey50", colour = "grey50") + geom_sf(data = stDat, aes(fill = GWRmodel$SDF$gwr.e), colour = NA) + 
  scale_fill_gradient2(low = "#a50026", 
                       mid = "white", 
                       high = "#006837",
                       name = "Rezidua (ČSSD)")
                                 
#map GWRresids-Regresids
ggplot() + geom_sf(data = uni, fill = "grey50", colour = "grey50") + geom_sf(data = stDat, aes(fill = changeRes), colour = NA) + 
  scale_fill_gradient2(low = "#a50026", 
                       mid = "white", 
                       high = "#006837",
                       name = "Změna")          

#plot local R2
plotKoef("localR2")

#show predictor..
ggplot() + geom_sf(data = uni, fill = "grey50", colour = "grey50") + geom_sf(data = stDat, aes(fill = VERICI), colour = NA) + 
  scale_fill_gradient2(low = "#a50026", 
                       mid = "white", 
                       high = "#006837",
                       name = "Věřící [%]")
