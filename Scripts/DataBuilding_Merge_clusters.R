# IMPORTANT: Script to merge clusters to data with values 
## 
library(dplyr)
library(tidyr)
library(broom)
library(rgdal)
# gen and den data from Databuliding_note.Rmd

# data <- unzip("/home/amira.azizan/Documents/Y320/GenDivDatabase/Data/TerrestrialEcos.zip")
# biogeo <- rgdal::readOGR(data[7])
# levels(biogeo$REALM)[levels(biogeo$REALM) == "NA"] <-  "NE"
# biogeo$BIOME <- as.factor(biogeo$BIOME)
# saveRDS(biogeo, "Data/biogeo.rds")


##################
# Gendata + biogeo
getREALMS <- function(df){
  biogeo <- readRDS("Data/biogeo.rds")
  points <- sp::SpatialPoints(coords=df[,c("long","lat")], proj4string = sp::CRS("+proj=longlat +datum=WGS84 +no_defs"))
  biogeo.point <- raster::extract(x=biogeo[,c("ECO_NAME","REALM","BIOME")], y=points)
  df <- cbind(df, biogeo.point[,c("ECO_NAME","REALM","BIOME")])
  df
}

gendata <- getREALMS(gen[!is.na(gen$long),])
dendata <- getREALMS(den[!is.na(den$long),])

###################################
## Merge gendata with clus.res to identify the clusters
# Initialise df
gen_clus <- data.frame()
xl <- list()
yl <- list()

for (j in sort(unique(gendata$sp))){
    
    gendata_sp <- subset(gen, sp %in% j)
    dendata_sp <- subset(den, sp %in% j)
    clus_sp <- subset(clus.res, sp %in% j)
    clus_sp <- clus_sp[which(clus_sp$distance.lowerthanDI == TRUE),]
    
    if(nrow(clus_sp) <= 1) next
    
    #if(j == "Panthera leo") {
    #  # subset to TD using SECR only
    #  dendata_sp <- dendata_sp[grep("SECR", dendata_sp$td),]
    #}
    
    # Separate data based on source
    clus_sp_gen <- clus_sp[which(clus_sp$source == "Popgen"),c("country","long","lat","cluster_center","distance.lowerthanDI")]
    clus_sp_gen <- unique(clus_sp_gen)
    
    clus_sp_den <- clus_sp[which(clus_sp$source == "Popden"),c("ref","country","long","lat","cluster_center","distance.lowerthanDI")]
    clus_sp_den <- unique(clus_sp_den)

    x <- dplyr::left_join(gendata_sp,clus_sp_gen, by=c("country","long","lat"))
    if(any(is.na(x$cluster_center))) x <- x[!is.na(x$cluster_center),]
    xl[[j]] <- x
    
    y <- dplyr::left_join(clus_sp_den, dendata_sp)
    if(any(is.na(y$cluster_center))) y <- y[!is.na(y$cluster_center),]
    y$cluster_center <- droplevels(as.factor(y$cluster_center))
    yl[[j]] <- y

# # Calculate the DP average across studies in each cluster_center before merging with the GD
    DP_mean <- aggregate(dp ~ cluster_center, data = y, FUN=mean, na.action=na.omit)
    DP_max <- aggregate(dp ~ cluster_center, data = y, FUN=max, na.action=na.omit)
    DP_sd <- aggregate(dp ~ cluster_center, data = y, FUN=sd, na.action=na.omit)
    DP_n <- aggregate(dp ~ cluster_center, data = y, FUN=n_distinct, na.action=na.omit)
    DP_median <- aggregate(dp ~ cluster_center, data = y, FUN=median, na.action=na.omit)
    
    DP_mean <- data.frame(cluster_center=DP_max[,1], DP_max=DP_max[,2], DP_mean=DP_mean[,2], DP_sd=DP_sd[,2], DP_n=DP_n[,2], DP_median=DP_median[,2])
    
    DP_ref <- tapply(y$ref, y$cluster_center, unique)
    DP_ref <- lapply(DP_ref, paste, collapse=", ")
    DP_ref <- do.call(rbind,DP_ref)
    DP_ref <- data.frame(ref=DP_ref, cluster_center=rownames(DP_ref))
    
    DP_mean <- merge(DP_mean, DP_ref, all.x=T)
    
    DP_pop <- tapply(y$population, y$cluster_center, unique)
    DP_pop <- lapply(DP_pop, paste, collapse=", ")
    DP_pop <- do.call(rbind,DP_pop)
    DP_pop <- data.frame(population=DP_pop, cluster_center=rownames(DP_pop))
    
    
    DP_mean <- merge(DP_mean, DP_pop, all.x=T)
    
    gen_dp <- dplyr::left_join(x, DP_mean, by="cluster_center", suffix=c(".gen",".den"))

    gen_clus <- rbind(gen_clus, gen_dp)



}

gen_clus

# remove rows without DP
gen_clus <- gen_clus[!is.na(gen_clus$DP_mean),]



##################################
# Aggregate data per country
## Merge gendata with popgen based on country name


gen$SxC <- paste(gen$sp, gen$country, sep="_")
den$dp <- as.numeric(den$dp)
den$SxC <- paste(den$sp, den$country, sep="_")
gendata <- gen[order(gen$sp),]

dendata.country <- den[which(den$SxC %in% gen$SxC),]
table(dendata.country$sp)

dendata.c <- dendata.country %>%
    group_by(sp,country) %>%
    summarise(
        DP_n = n(),
        N.SS = n_distinct(population), 
        DP_mean = mean(dp, na.rm = TRUE), 
        DP_sd = sd(dp, na.rm = TRUE), 
        DP_ref = paste(unique(ref), collapse=", "),
        DP_median = median(dp, na.rm = TRUE)
    )


gen_count <- merge(gendata, dendata.c, by=c("sp","country"), all.x=TRUE)
gen_count <- gen_count[!is.na(gen_count$DP_mean),]

colnames(gen_count) <- gsub("_mean", "", colnames(gen_count))


##################################
# Aggregate data per econame
# data source: https://academic.oup.com/bioscience/article/51/11/933/227116?login=false
## Merge gendata with popgen based on econame name


gendata <- getREALMS(gen[!is.na(gen$long),])
dendata <- getREALMS(den[!is.na(den$long),])

gendata$SxE <- paste(gendata$sp, gendata$ECO_NAME, sep="_")
dendata$dp <- as.numeric(dendata$dp)
dendata$SxE <- paste(dendata$sp, dendata$ECO_NAME, sep="_")

dendata.eco <- dendata[which(dendata$SxE %in% gendata$SxE),]
table(dendata.eco$sp)

dendata.e <- dendata.eco %>%
  group_by(sp,ECO_NAME) %>%
  summarise(
    DP_n = n(),
    N.SS = n_distinct(population), 
    DP_mean = mean(dp, na.rm = TRUE), 
    DP_sd = sd(dp, na.rm = TRUE), 
    DP_ref = paste(unique(ref), collapse=", "),
    DP_median = median(dp, na.rm = TRUE), 
    
  )


gen_eco <- merge(gendata, dendata.e, by=c("sp","ECO_NAME"), all.x=TRUE)
gen_eco <- gen_eco[!is.na(gen_eco$DP_mean),]





## Add body size (from Johnson's data)
pdenIV <- read.csv("/home/amira.azizan/Documents/Y320/GenDivDatabase/Data/etc/pdenIV.csv")
pdenIV$sp <- pdenIV$Sname
pdenIV$sp <- reorder(pdenIV$sp, pdenIV$BS)
BS_df <- unique(pdenIV[,c("sp","BS")])
BS_DF <- BS_df[BS_df$sp %in% gen_clus$sp,]
gen_clus$BS <- NULL


gen_clus <- left_join(gen_clus, BS_df)
gen_count <- left_join(gen_count, BS_df)
gen_eco <- left_join(gen_eco, BS_df)

# Allocate generation time (from Pacifici et al.)
gentime <- read.csv("Data/etc/GenerationtimeFelidae.csv", sep=",")
gentime$GL <- round(gentime$GenerationLength_d/365,0)
names(gentime)[[1]] <- "sp"

gen_clus <- left_join(gen_clus, gentime[,c("sp","GL","GenerationLength_d")])
gen_count <- left_join(gen_count, gentime[,c("sp","GL","GenerationLength_d")])
gen_eco <- left_join(gen_eco, gentime[,c("sp","GL","GenerationLength_d")])


## Add DP x home range
## Home range
sigma <- read.csv("Data/etc/sigma.csv")
sigma$Sname
names(sigma)[1] <- "sp"
gen_clus <- left_join(gen_clus, sigma[,c("sp","HR","DI")])
gen_count <- left_join(gen_count, sigma[,c("sp","HR", "DI")])
gen_eco <- left_join(gen_eco, sigma[,c("sp","HR", "DI")])


## Add conservation status (28/04/2020)
# Only otocolobus manul's status has changed in 2019/2020
iucn <- read.csv("Data/etc/IUCN/redlist_Felidae/simple_summary.csv")
iucn <- iucn[,c("sp","genusName","redlistCategory","populationTrend")]
iucn$redlistCategory <- factor(iucn$redlistCategory, levels=c("Endangered","Vulnerable","Near Threatened","Least Concern"))
iucn$iucn <- iucn$redlistCategory
levels(iucn$iucn) <- c("T","T", "NT", "NT") 


gen_clus <- left_join(gen_clus, iucn)
gen_count <- left_join(gen_count, iucn)
gen_eco <- left_join(gen_eco, iucn)

#gen_clus$population <- gsub("Whole Bhutan", "Bhutan", gen_clus$population)
gen_count$population <- gsub("Whole Bhutan", "Bhutan", gen_count$population)

#gen_clus$population <- gsub("Tanzania", "United Republic of Tanzania", gen_clus$population)
gen_count$population <- gsub("Tanzania", "United Republic of Tanzania", gen_count$population)




gen_clus <- unique(gen_clus)
gen_count <- unique(gen_count)
gen_eco <- unique(gen_eco)
rownames(gen_clus) <- 1:nrow(gen_clus)

gen_clus <- gen_clus[-which(gen_clus$min_time < 1970),]
summary(gen_clus$min_time)

gen_count <- gen_count[-which(gen_count$min_time < 1970),]

gen_eco <- gen_eco[-which(gen_eco$min_time < 1970),]

gen_clus <- droplevels(gen_clus)
gen_count <- droplevels(gen_count)
gen_eco <- droplevels(gen_eco)
# Arrange species by size
gen_clus$sp <- factor(gen_clus$sp, levels=as.character(unique(gen_clus$sp[order(gen_clus$BS)])))
gen_count$sp <- factor(gen_count$sp, levels=as.character(unique(gen_count$sp[order(gen_count$BS)])))
gen_eco$sp <- factor(gen_eco$sp, levels=as.character(unique(gen_eco$sp[order(gen_eco$BS)])))
# Cluster data ------------------------------------------------------------
# Number of cluster per species
tapply(gen_clus$cluster_center, gen_clus$sp, n_distinct)


gen_clus$genus
gen_clus$label <- paste0("atop(italic(",substr(gen_clus$genus, start = 1, stop = 1),".", gsub(".*. ", "", gen_clus$sp),"),", gsub(" ", "~", trimws(gen_clus$population)),")")

# Transform variables

gen_clus$DP_log <- log10(as.numeric(gen_clus$DP_mean))
gen_clus$BS_log <- log10(gen_clus$BS)
# 
# Transform variables
gen_count$DP_log <- log10(gen_count$DP)
gen_count$BS_log <- log10(gen_count$BS)

gen_eco$DP_log <- log10(as.numeric(gen_eco$DP_mean))
gen_eco$BS_log <- log10(gen_eco$BS)


clus.den <- clus.res[clus.res$source == "Popden",]
dendata <- merge(dendata, clus.den)

clus.gen <- clus.res[clus.res$source == "Popgen",]
gendata <- merge(gendata, clus.gen)

gen_count.l <- gen_count

#####################################
#save(gendata, dendata, file="Data/GenDen_Rawdata_V4.rda")
save(gen_clus, gen_count, file="Data/Popgenden_biogeo_V4.rda")
#save(gen_clus.l, gen_count.l, file="Data/Popgenden_biogeo_V4_LION.rda")
#####################################

