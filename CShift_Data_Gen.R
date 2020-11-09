start.time <- Sys.time()

library(haven)
library(raster)
library(ncdf4)
library(maptools)
library(rgdal)
library(parallel)
library(doSNOW)
library(snow)
library(reshape2)
library(sqldf)
library(tidyverse)
library(caret)
library(RSQLite)
library(DBI)
library(sf)
library(tiff)
library(stringr)
library(sf)
library(rgeos)
library(data.table)
library(fastDummies)
library(foreach)
library(VIM)
library(velox)

## Functions

bin_fun = function(x){
  x[is.na(x)] = 0
  x = ifelse(x == 1,1,0)
}

na.i = function(A,B){
  setDT(cs.df)[is.na(A), A:=B] 
}

my.mean = function(x){ mean(x,na.rm=TRUE) }

is.nan.data.frame <- function(x){
  do.call(cbind, lapply(x, is.nan))
}

####

setwd("C:/Users/fc3/Box Sync/CShift/Data")

numParallelCores <- max(1, detectCores()-1)
cl <- makeCluster(rep("localhost", numParallelCores), type = "SOCK")
registerDoSNOW(cl)

gfsf_sql <- dbConnect(RSQLite::SQLite(), dbname = "CShift_db.sqlite3")

################################################ Agricultural census###################################################################

files = files <- list.files("C:/Users/fc3/Box Sync/CShift/Data/Total_nacional", pattern = "*.sav$")
nm.dt = c("prod","crop","hhd","ppl","home")

ag.cs = list()
for (i in 1:length(files)) {
  path = file.path(paste0("C:/Users/fc3/Box Sync/CShift/Data/Total_nacional/",files[i]))
  a = read_sav(path)
  if(i==2){
    a$P_S6P46 = str_replace(a$P_S6P46, "^0+" ,"")
    a$broad_crop = substr(a$P_S6P46, 1, 7)
    a = a %>% group_by(ENCUESTA,P_S6P45B) %>% mutate(id = row_number())
  }
  dbWriteTable(gfsf_sql,nm.dt[i], a, overwrite=T)
}

################################################## Veredas shapefile##################################################################

poly <- readOGR("C:/Users/fc3/Box Sync/GFSF/Maps/Veredas_de_Colombia.shp")
c.poly = gCentroid(poly, byid = T)
c.poly = cbind(poly@data[,c("DPTOMPIO","CODIGO_VER")],c.poly@coords)
c.poly = as.data.frame(sapply(c.poly,function(x) gsub("(^|[^0-9])0+", "",as.character(x))))
colnames(c.poly)[2] = "COD_VEREDA"

######################################## CHELSA temperature and precipitation#########################################################

nm.dt = c("prec","tmax","tmin")

s = list()
for (i in 1:length(nm.dt)) {
  grids <- list.files(paste0("C:/Users/fc3/Box Sync/CShift/Maps/",nm.dt[i]), pattern = "*.tif$")
  s[[i]] <- stack(paste0("C:/Users/fc3/Box Sync/CShift/Maps/",nm.dt[i],"/", grids))
}

pb <- txtProgressBar(min = 0, max = length(nm.dt)*12, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
x <- as.data.frame(foreach(i=1:3, .combine='cbind', .packages = c("raster","velox", "Matrix")) %:%
                     foreach(j=1:12, .combine='cbind',
                             .options.snow = opts) %dopar% {
                               ras_crop = crop(s[[i]]@layers[[j]],poly)
                               vx.ras = velox(ras_crop)
                               ex = vx.ras$extract(poly, fun=my.mean,small=T)
                               return(ex)
                             })
close(pb)

w <- foreach(i=1:3, .combine='c', .packages = c("raster")) %:%
  foreach(j=1:12, .combine='c',
          .options.snow = opts) %dopar% {
            name = gsub(".*CHELSAcruts_","",s[[i]]@layers[[j]]@file@name)
            name = gsub("_2014.*","",name)
            return(name)
          }

colnames(x) = w
x = cbind(poly@data$CODIGO_VER,x)
weather = left_join(x,c.poly[,2:4],by=c("poly@data$CODIGO_VER"="COD_VEREDA"))
colnames(weather)[1] = "COD_VEREDA"
weather$COD_VEREDA = gsub("(^|[^0-9])0+","",as.character(weather$COD_VEREDA))

##Impute small polygons with neighbors
na.nm = colnames(weather)[colSums(is.na(weather)) > 0]
weather[is.nan(weather)] <- NA
weather = VIM::kNN(weather,variable = na.nm,numFun = my.mean, dist_var = c("x","y"),imp_var=F)

close(pb)

for (j in 14:ncol(weather)) {
  dd = function(n){(n/10)}
  weather[,j] = dd(weather[,j])
}

############################################################ Altitude##############################################################

dem = raster("C:/Users/fc3/Box Sync/CShift/Maps/dem.tif")
vx.dem = velox(dem)
altitude = as.data.frame(foreach(fn = c(mean,sd), .combine='cbind', .packages = c("raster","velox", "Matrix")) %dopar% {
  a = vx.dem$extract(poly, fun=fn,small=T)
  return(a)
})
altitude = cbind(poly@data$CODIGO_VER,altitude)
altitude = left_join(altitude,c.poly, by = c("poly@data$CODIGO_VER"="COD_VEREDA"))
colnames(altitude) = c("COD_VEREDA","alt_mn","alt_sd","nd","lat","lon")
altitude$COD_VEREDA = gsub("(^|[^0-9])0+", "",as.character(altitude$COD_VEREDA))

# Impute missing values with mean of neighbors
na.nm = colnames(altitude)[colSums(is.na(altitude)) > 0]
altitude[is.nan(altitude)] <- NA
altitude = VIM::kNN(altitude,variable = na.nm,numFun = my.mean, dist_var = c("lat","lon"),imp_var=F)

################################################################## Prices########################################################

price = read.csv("price.csv")
price$DPTOMPIO = as.character(price$DPTOMPIO)

## Create Spatial Poits Datafrme of Markets
## Barranquilla, Bogota, Bucaramanga, Cali, Cartagena, Cucuta, Medellin, Pereira

mkt.pts = data.frame(LONGITUDE=c(-74.8070, -74.0721, -73.1227, -76.5320, -75.4832, -72.4967, -75.5742, -75.6906),
                     LATITUDE=c(11.0041, 4.7110, 7.1193, 3.4516, 10.3932, 7.8891, 6.2486, 4.8087),
                     MYID=c("Barranquilla", "Bogota", "Bucaramanga", "Cali", "Cartagena", "Cucuta", "Medellin", "Pereira"),
                     INDEX=1:8)
coordinates(mkt.pts) <- c("LONGITUDE","LATITUDE")

# parallel calculation of closest city
{## Export the environment variables to each cluster
  clusterExport(cl,ls())
  
  ## Load the library "rgeos" to each cluster
  clusterEvalQ(cl, library(rgeos))
  
  ## Split the data
  ID.Split<-clusterSplit(cl,unique(poly$OBJECTID))
  
  ## Define a function that calculates the distance of one ID in relation to the poly2
  a<-function(x) gDistance(spgeom1 = poly[x,], spgeom2 = mkt.pts, byid=TRUE)
  
  ## Run the function in each cluster
  system.time(m<-clusterApply(cl, x=ID.Split, fun=a))
  
  ## Merge the results
  output<- do.call("cbind", m)
  colnames(output) = poly@data$CODIGO_VER}

cl_min = function(x){which(grepl(min(x),x))}
close_city = apply(output,2,cl_min)

close_city = ifelse(close_city == "1","8001",
                    ifelse(close_city == "2","11001",
                           ifelse(close_city == "3", "68001",
                                  ifelse(close_city == "4","13001",
                                         ifelse(close_city == "5","76001",
                                                ifelse(close_city == "6","54001",
                                                       ifelse(close_city == "7","5001","66001")))))))
close_city = as.data.frame(cbind(poly@data[,c("CODIGO_VER","DPTOMPIO")],close_city))
price = left_join(close_city,price,by=c("close_city"="DPTOMPIO"))

### Change to wide format
month = c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")

wide_price = reshape(price,idvar = "CODIGO_VER",timevar = "broad_crop", direction = "wide")
cod_ver = wide_price[,c("close_city.Avocado","DPTOMPIO.Avocado","nm_mpio.Avocado" )]
wide_price = wide_price[, -grep("DPTOMPIO", colnames(wide_price))]
wide_price = wide_price[, -grep("close_city", colnames(wide_price))]
wide_price = wide_price[, -grep("nm.mpio", colnames(wide_price))]
wide_price = cbind(cod_ver,wide_price)

colnames(wide_price)[4] = "COD_VEREDA"
wide_price$COD_VEREDA = gsub("(^|[^0-9])0+", "",as.character(wide_price$COD_VEREDA))

######################################################### Colombian Soils############################################################

poly2 = st_read("C:/Users/fc3/Box Sync/GFSF/Maps/Veredas_de_Colombia.shp")
crs = crs(poly)
soil <- st_read("C:/Users/fc3/Box Sync/GFSF/Maps/soils.shp")
soil = st_transform(soil[,14],crs)
pi <- st_intersection(soil, poly2)

## add in areas in m2
attArea <- pi %>%
  mutate(area = st_area(.) %>% as.numeric())

## for each field, get area per soil type
ucs_ver = attArea %>%
  as_tibble() %>%
  dplyr::group_by(CODIGO_VER, COMPONENTE) %>%
  dplyr::summarize(area = sum(area))

ucs_ver$COMPONENTE = as.character(ucs_ver$COMPONENTE)
soil_type = data.frame(ucs_ver$COMPONENTE, do.call(rbind, strsplit(ucs_ver$COMPONENTE, split = ":", fixed = TRUE)))

soil_ver <- cbind(data.frame(ucs_ver),soil_type[,2:3])
wide_soil = dcast(soil_ver, CODIGO_VER ~ X1, value.var = "area", drop = FALSE)

wide_soil$`Asociación Campoalegre, San Juan` <- as.integer(wide_soil$`Asociación Campoalegre, San Juan`|wide_soil$`Asociación Campoalegre - San Juan`)
wide_soil$`Asociación Campoalegre - San Juan` = NULL
wide_soil$`Asociación Guarinó - Samaná` <- as.integer(wide_soil$`Asociación Guarinó - Samaná`|wide_soil$`Asociación Guarinó, Samaná`)
wide_soil$`Asociación Guarinó, Samaná` = NULL
wide_soil$`Asociación Iberia - San Lorenzo` <- as.integer(wide_soil$`Asociación Iberia - San Lorenzo`|wide_soil$`Asociación Iberia, San Lorenzo`)
wide_soil$`Asociación Iberia, San Lorenzo` = NULL
wide_soil$`Asociación Peñas, Ventanas` <- as.integer(wide_soil$`Asociación Peñas, Ventanas`|wide_soil$`Asociación Peñas - Ventanas`)
wide_soil$`Asociación Peñas - Ventanas` = NULL
wide_soil$`Asociación Santa Isabel - Pensilvania` <- as.integer(wide_soil$`Asociación Santa Isabel - Pensilvania`|wide_soil$`Asociación Santa Isabel, Pensilvania`)
wide_soil$`Asociación Santa Isabel, Pensilvania` = NULL
wide_soil$`Asociación Taudia - Chinchiná` <- as.integer(wide_soil$`Asociación Taudia - Chinchiná`|wide_soil$`Asociación Taudia, Chinchiná`)
wide_soil$`Asociación Taudia, Chinchiná` = NULL
wide_soil$`Asociación Chinchiná - Azufrado` <- as.integer(wide_soil$`Asociación Chinchiná - Azufrado`|wide_soil$`Asociación Chinchimá, Azufrado`)
wide_soil$`Asociación Chinchimá, Azufrado` = NULL
wide_soil$`Asociación Río Arma - Castillo` <- as.integer(wide_soil$`Asociación Río Arma - Castillo`|wide_soil$`Asociación Rioarma, Castilla`)
wide_soil$`Asociación Rioarma, Castilla` = NULL
wide_soil$Consociación <- as.integer(wide_soil$Consociación|wide_soil$Consicuación)
wide_soil$Consicuación = NULL
wide_soil$`Grupo indiferenciado` <- as.integer(wide_soil$`Grupo indiferenciado`|wide_soil$`Grupo Indiferenciado`)
wide_soil$`Grupo Indiferenciado` = NULL
wide_soil$`Misceláneo rocoso` <- as.integer(wide_soil$`Misceláneo rocoso`|wide_soil$`Misceláneos Rocosos`)
wide_soil$`Misceláneos Rocosos` = NULL
wide_soil$`NA` <- as.integer(wide_soil$`NA`|wide_soil$`N/A`)
wide_soil$`N/A` = NULL
wide_soil$`Asociación Fluventic Dystrudepts; Oxic Dystrudepts` = NULL
wide_soil$`Asociación; Typic Hapludox; Typic Kandiudox; Typic Humaquepts; Typic Dystrudepts` = NULL
wide_soil$Complejo = NULL
wide_soil$`NA` = NULL

colnames(wide_soil)[1] = "COD_VEREDA"
wide_soil$COD_VEREDA = as.character(wide_soil$COD_VEREDA)
wide_soil$COD_VEREDA = gsub("(^|[^0-9])0+", "",as.character(wide_soil$COD_VEREDA))

############################### Merge: weather, soil, altitude, prices ###########################################################################

ls = c("weather","wide_soil","altitude","wide_price")
ls2 = list()
for (k in 1:length(ls)) {
  a = as.data.frame(get(ls[k]))
  a = sapply(a, as.numeric)
  ls2[[k]] = as.data.table(a)
}

ver_df = Reduce(function(x,y) merge.data.table(x = x, y = y, by = "COD_VEREDA", all=T),
                ls2)
ver_df$COD_VEREDA = gsub("(^|[^0-9])0+", "",as.character(ver_df$COD_VEREDA))
cult = read.csv("crop.csv")
crop = dcast(cult, COD_VEREDA ~ cult, value.var = "AREA_COSECHADA", drop = FALSE)
coffee_ver = crop$COD_VEREDA[which(crop$`1610050` != 0 | crop$`1610010` != 0 | crop$`1610060` != 0 | crop$`1610070` != 0 )]
ver_df = ver_df[which(ver_df$COD_VEREDA %in% coffee_ver),]

# dbWriteTable(gfsf_sql, "ver_dfcl", ver_df, overwrite=T)

###########################################
############SQL############################
###########################################

gc()

## gfsf_sql database created in sql with tables: prod, cult, ppl, hhd, home, weather, soil ##
## subset tables by coffee producing veredas
## subset prod by UPA
## drop unfiltererd tables

co_ver <- paste0("'", coffee_ver, "'", collapse=", ")
whereIn <- paste0("(", co_ver, ")")

myTables <- c("crop","prod","ppl","hhd","home") 
for (i in 1:length(myTables)) {
  if(i == 1){
    f = dbReadTable(gfsf_sql,"crop")
    f = f[which(f$P_S6P45A==1),]
    dbWriteTable(gfsf_sql, "crop", f, overwrite=T)
  }
  if(i == 2){
    g = dbReadTable(gfsf_sql,"prod")
    g = g[which(g$TIPO_UC==1),]
    dbWriteTable(gfsf_sql, "prod", g, overwrite=T)
  }
  sqlStatement <- paste("SELECT *",
                         "FROM", myTables[i],
                         "WHERE COD_VEREDA",
                         "IN",whereIn)
  a = dbGetQuery(gfsf_sql, sqlStatement)
  if(myTables[[i]] == "crop"){
    b = a$ENCUESTA
  }
  if(myTables[[i]] != "crop"){
    a = a[which(a$ENCUESTA %in% b),]
  }
  dbWriteTable(gfsf_sql, paste0(myTables[i],"cl"), a, overwrite=T)
  dbRemoveTable(gfsf_sql,myTables[i])
}

nm.db = c(paste0(myTables,"cl"),"ver_dfcl")
d = list()
for (i in 1:length(nm.db)) {
  if(i %in% 1:5){
    a = as.data.table(dbReadTable(gfsf_sql, nm.db[i]))
    setkey(a,"ENCUESTA")
    d[[i]] = a
  } else {
    a = as.data.table(dbReadTable(gfsf_sql, nm.db[i]))
    a <- a[, COD_VEREDA:=as.character(COD_VEREDA)]
    setkey(a,"COD_VEREDA")
    d[[i]] = a
  }
}

cs.df = merge.data.table(d[[1]], unique(d[[6]]), all.x=TRUE, by='COD_VEREDA')
e = Reduce(function(x, y) x[unique(y, by="ENCUESTA"), by = "ENCUESTA", all = T], d[2:5])
cs.df = merge.data.table(cs.df, unique(e,by="ENCUESTA"), all.x=TRUE, by='ENCUESTA')
cs.df = cs.df %>% select(-contains("i."))
cs.df = cs.df %>% select(-contains(".y"))
colnames(cs.df) = gsub(".x", "",as.character(colnames(cs.df)))

## Impute NAs in binary columns
## Create dummies for cat variables
## transform cat variables to factor

cat = c("P_S6P71","P_S7P95","P_S8P109","P_S8P113A","P_S4P15","S05_TENENCIA",
        "P_S6P47A","P_S6P50","P_S6P60","P_S15P160","P_S15P162",
        "P_S15P163","P_S15P178","P_S15P170","P_S15P172","P_S15P175A","P_S15P176")

cs.df = dummy_cols(cs.df, select_columns = cat, remove_selected_columns = T)
cs.df = cs.df[, !duplicated(as.list(cs.df)), with = FALSE]

v = sapply(cs.df, function(x) { length(unique(x[!is.na(x)]))<4 })
bin = colnames(cs.df)[v]
cs.df[, (bin) :=  lapply(.SD, bin_fun), .SDcols = bin]

## Remove plot observations with no productivity per hectare (P_S6P59_UNIF)
## Remove columns with more that 40% NA and impossible to impute
## Impute NAs of livestock columns depending on the base question

cs.df = na.omit(cs.df,cols = "P_S6P59_UNIF")

drop = c("Asociación","Inasociación","close_city.Avocado","DPTOMPIO.Avocado","nm_mpio.Avocado","COD_PARQUE",
         "PRED_ETNICA","P_S4P16","SNH","SNM","SN9","P_S5PAUTOS","P_S6P66","P_S6P68","P_S11P138A","P_S11P138B",
          "P_S11P139","P_S11P139A","P_S11P139B","P_S11P140","P_S12P142","P_S12P143","P_S12P144","P_S12P145","P_S12P146",
         "P_S12P147","P_S12P148","P_S12P149","P_S15P158","P_S15P158A","P_S15P158B","P_S15P165","P_S15P165A",
         "P_S15P169","P_S15P171","P_S15P175B","P_S15P175C","ITER_HG","ID_PROD","P_S15P179A","P_S15P179B","P_S15P179C",
         "TOT_PROD_HOGAR","P_S15P161","P_S11P138")
cs.df[, (drop) := NULL]

### Cattle
im_na = c("P_S7P83A","P_S7P83B","P_S7P83C","P_S7P83D",
          "P_S7P83E","P_S7P83F","P_S7P84A","P_S7P84B",
          "P_S7P84C","P_S7P84D","P_S7P84E","P_S7P84F","P_S7P85B")
cs.df[, (im_na) :=  lapply(.SD, function(x) ifelse(is.na(x),cs.df$P_S7P82[which(is.na(x))], x)), .SDcols = im_na]

### Hogs
im_na = c("P_S7P87_SP1","P_S7P87_SP2","P_S7P87_SP3",
          "P_S7P87_SP4","P_S7P88","P_S7P89A","P_S7P89B",   
          "P_S7P89C","P_S7P89D","P_S7P89E","P_S7P89F")
cs.df[, (im_na) :=  lapply(.SD, function(x) ifelse(is.na(x),cs.df$P_S7P86[which(is.na(x))], x)), .SDcols = im_na]

### Poultry
im_na = c("P_S7P92A","P_S7P92B","P_S7P93A","P_S7P93B")
cs.df[, (im_na) :=  lapply(.SD, function(x) ifelse(is.na(x),cs.df$P_S7P90[which(is.na(x))], x)), .SDcols = im_na]

### Other major species
im_na = c("P_S7P102A","P_S7P102B","P_S7P102C",
          "P_S7P102D","P_S7P102E","P_S7P102F",
          "P_S7P102G","P_S7P102H","P_S7P102I",
          "P_S7P102J","P_S7P102K","P_S7P102L")
cs.df[, (im_na) :=  lapply(.SD, function(x) ifelse(is.na(x),cs.df$P_S7P101[which(is.na(x))], x)), .SDcols = im_na]

### Other minor species
im_na = c("P_S7P106A","P_S7P106B","P_S7P106C",
          "P_S7P106D","P_S7P106E","P_S7P106F",
          "P_S7P106G","P_S7P106H","P_S7P106I",
          "P_S7P106J","P_S7P106K","P_S7P106L",
          "P_S7P106M","P_S7P106N","P_S7P106O",
          "P_S7P106P","P_S7P106Q")
cs.df[, (im_na) :=  lapply(.SD, function(x) ifelse(is.na(x),cs.df$P_S7P105[which(is.na(x))], x)), .SDcols = im_na]

smp = sample_n(cs.df,20000)
vc.crop = c("1610010","1313010","1802010","1312010","1640010","1311010","1318010","1919990")
vc.citr = c("1323010","1324010","1322010")
smp$broad_crop = ifelse((smp$broad_crop %in% vc.citr), "1323010",
                        ifelse(!(smp$broad_crop %in% vc.crop), "9999999",smp$broad_crop))
smp = sapply(smp,as.numeric)
# write.csv(smp,"sample_df.csv",row.names = F)
# dbWriteTable(gfsf_sql, "cs.df", cs.df, overwrite=T)

stopCluster(cl)

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken
