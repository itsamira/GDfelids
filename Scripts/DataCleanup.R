library(ggplot2)
library(dplyr)

##
# Genetic data -----------------------------------------------------------
# Load pooled data 
dfll <- read.csv("/home/amira.azizan/Documents/Y320/GenDivDatabase/Data/Rawdata/Link to Microsat_rawpaper_V10.csv", sep="\t", stringsAsFactors=FALSE)

# Fix locus name
fixlocusname <- function(x){  
  # remove space
  x <- gsub(" ","",x)
  # toupper
  x <- toupper(x)
  # if FCA leading zero 
  if (grepl("^FC",x)) x <- paste0("FCA",sprintf("%03d", as.numeric(gsub("[A-z]","", x))))
  if (grepl("^F[0-9]",x)) x <- paste0("F",sprintf("%03d", as.numeric(gsub("[A-z]","", x))))
  return(x)
}
dfll$locus <- sapply(dfll$locus, fixlocusname)


## Remove nonFCA loci
CleanDF <- function(x){
  # Where there is a mix of locus type (FCA + other species), keep the most locus
  loc <- gsub("[0-9]", "", substr(x$locus,1,3))
  if (any(grep("F$",loc))) loc[grep("F$",loc)] <- 'FCA'
  sl <- which(max(table(loc)) == table(loc))
  if (length(sl) == 1 & !all(grepl("^F", names(sl)))){
    x <- x[grep(names(sl), x$locus),]
  } else x <- x[grep("^F", x$locus),]
  return(x)
}

ncol(dfll)
nrow(dfll)
# Assign to another vector
df.1 <- dfll
# Clean columns /!\
names(df.1)
df.1 <- df.1[,-c(3,4,15, 27:30)]

# Load estimates calculated straight from genotype data
file_raw <- "/home/amira.azizan/Documents/Y320/GenDivDatabase/Data/Rawdata/Rawdata_V3.csv"
# rawdata <- read.csv("Data/Rawdata_V1.csv", sep="\t", stringsAsFactors=FALSE)
rawdata <- read.csv(file_raw, sep="\t", stringsAsFactors=FALSE)
rawdata <- rawdata[,c(1:3,8,10,18:21)]
names(rawdata)
# sp : species name
# source : either from raw data or reported in study
# ref : study reference (format = first author_first title word_year)
# population name. Global = includes all population in genotype data.
# locus : locus name. If AVERAGE, No locus specific estimates reported.

#####################
# Fix data #####
# Filter 1 : Population with multiple countries
df.filter1 <- df.1 %>% filter(grepl("\\,",country) | is.na(country))

#' Multiple countries, No raw data
df.filter1.noraw <- df.filter1[!df.filter1$ref %in% rawdata$ref,]

# Filter 2 : Population with mixed loci (FCA + other species loci)
df.filter2 <- df.1 %>% filter(grepl(" ",Species.specific))
df.filter2.noraw <- df.filter2 %>% filter(!ref %in% rawdata$ref)


# Filter 1 + 2
df.2 <- df.1 %>% filter(grepl("\\,",country) | is.na(country) | grepl(" ",Species.specific) | grepl("\\+",population))

refl <- unique(df.2$ref)

# raw data with different pop in dfll
k.1 <- intersect(rawdata$ref, dfll$ref)
df.a <- dfll[dfll$ref %in% k.1,]
pop.ref.raw <- tapply(rawdata$pop[rawdata$ref %in% k.1],rawdata$ref[rawdata$ref %in% k.1], function(x) length(unique(x)))
pop.ref.pap <- tapply(df.a$pop, df.a$ref, function(x) length(unique(x)))
pop.df <- data.frame(raw=pop.ref.raw, paper=pop.ref.pap)
pop.df[,"diff"] <- pop.df[,1] - pop.df[,2]
refl <- c(refl, rownames(pop.df[pop.df$diff > 1,]))
refl <- unique(refl)

# raw data from paper with filter 1+2 : population with multiple countries and mixed loci
raw.1 <- rawdata[which(rawdata$ref %in% refl),]

# raw data references not in dfll
k <- setdiff(rawdata$ref, dfll$ref)
raw.1 <- rbind(raw.1, rawdata[rawdata$ref %in% k,])

# Fixed raw data information to be identical as in reported paper
subs <- names(raw.1) %in% names(df.1) #subset columns
raw.1 <- raw.1[,subs]
raw.1$country <- as.character(NA)

# raw.1 <- raw.1[raw.1$ref %in% k,]
## Some manual modification had to be done so that the population names are +/- standardised (spellings, localities in raw genotype data =/= locality reported in the paper etc.) as well as geographic coordinates
file_pop <- "/home/amira.azizan/Documents/Y320/GenDivDatabase/Data/PopName_V2.csv"
popname <- read.csv(file_pop, sep="\t", stringsAsFactors=FALSE)
popnamelong <- tidyr::gather(popname, popsource, population, Pop_report, Pop_raw)
popnamelong <- popnamelong[!is.na(popnamelong$population) & !is.na(popnamelong$country),]
popnamelong$popsource <- NULL
popnamelong <- unique(popnamelong)

changepopName <- function(x){
  
  s <- which(popname$ref %in% x[["ref"]] & popname$sp %in% x[["sp"]] & popname$Pop_raw %in% x[["population"]])
  if (!length(s)==0){
    if(!is.na(popname$Pop_report[s])){
      x[["population"]] <- paste(as.character(popname[s,"Pop_report"]))
    }
    x[["country"]] <- paste(as.character(popname[s,"country"]))
  }
  
  as.data.frame(t(x), stringsAsFactors=F)
}

repairpopname <- function(rawdf){
  rawdf$population <- trimws(rawdf$population)
  df.1 <- apply(rawdf, 1, changepopName)
  dfl <- do.call(plyr::rbind.fill,df.1)
  cl <- sapply(rawdf, class)
  if(!any(names(cl) %in% "country")){
    cl <- c(cl, "character")
    names(cl)[[length(cl)]] <- "country"
  }
  cl[which(cl == "factor")] <- "character"
  funs <- sapply(paste0("as.",cl), match.fun)
  dfl[] <- Map(function(dd,f) f(as.character(dd)), dfl, funs)
  dfl
}



raw.2 <- apply(raw.1, 1, changepopName)

dfl <- do.call(plyr::rbind.fill,raw.2)
cl <- sapply(raw.1, class)
funs <- sapply(paste0("as.",cl), match.fun)
dfl[] <- Map(function(dd,f) f(as.character(dd)), dfl, funs)
raw.3 <- dfl

## 

lsref <- unique(raw.3$ref)
raw.4 <- data.frame()
for (i in lsref) {
    foo <- raw.3[raw.3$ref %in% i,]
    # Remove Global population
    if(any(foo$population == "Global")) foo <- foo[-which(foo$population == "Global"),]
    pref <- popnamelong[popnamelong$ref %in% foo$ref,]
    pref <- pref[,c("ref","sp","country","min_time","max_time","Latitude","Longitude","population","Check_GD")]

    if (all(is.na(foo$country))) {
        foo$country <- NULL
        foo <- merge(foo, pref, by=c("sp","population","ref"), all.x=TRUE)
    } else { foo <- merge(foo, pref, by=c("sp","population","country","ref"), all.x=TRUE) }

    # Remove unwanted population (checked manually)
    foo <- foo[which(foo$Check_GD == 1),]
    raw.4 <- rbind(raw.4,foo)
}

# Raw data replacing reported estimates 

df.3 <- plyr::rbind.fill(
        df.1 %>% filter(!ref %in% unique(raw.4$ref)), 
        # remove data from df.1 which has raw data)
        raw.4)
df.3[,c(1:5,22:24)] <- sapply(df.3[,c(1:5,22:24)], as.character)
df.3 <- df.3[!is.na(df.3$ref),]


## Tag population with average or locus sp estimates
df.3$AVERAGE <- ifelse(df.3$locus == "AVERAGE","AVE","LOCUS")
df.3$locus[df.3$locus =="AVERAGE"] <- NA

df.3a <- df.3[,!grepl("_SD|_SE|per_polylocus|Check_GD|No.of.Loci", colnames(df.3))]
df.3a$ref <- as.factor(df.3a$ref)
df.3a$locus <- sapply(df.3a$locus, fixlocusname)





write.table(df.3a, "Data/Gendiv_clean_V2.csv", sep="\t")

# Function to calculate average estimates per population
sum_index <- list(mean=~mean(.x, na.rm=TRUE), SD=~sd(.x, na.rm=TRUE))

# Locus-specific estimates with low sample size was removed
df.locus <- df.3a %>% filter(AVERAGE == "LOCUS", sample.size > 5 ) %>% 
    group_by(sp,ref,country,population,Latitude,Longitude, min_time, max_time) %>% 
    mutate(mono=n.allele <= 1) %>% # Set monomorphic locus
    group_modify(~CleanDF(.x)) %>%
    summarise(across(where(is.numeric), sum_index), 
        No.of.Loci=n_distinct(locus), 
        per_polylocus=(No.of.Loci-sum(mono))/No.of.Loci,
        AVERAGE=paste("LOCUS"),
        Species.specific = unique(ifelse(is.na(Species.specific), 
    toupper(paste0(unique(gsub("[0-9]","",substr(locus,1,3))), collapse=";")), paste(unique(Species.specific)))),
        ) %>% 
arrange(ref)
names(df.locus) <- gsub("_mean$","", names(df.locus))

# Add estimates with averaged estimates (exclude data with raw genotype)
df.4 <- df.3 %>% filter(AVERAGE == "AVE", !ref %in% raw.3$ref)


df.final <- plyr::rbind.fill(df.locus, df.4) %>% arrange(ref)
df.final$locus <- NULL
df.final$source <- NULL

# Add variance based on number of loci
s <- is.na(df.final$Ho_SD) & !is.na(df.final$Ho_SE)
df.final[s,]$Ho_SD <- df.final[s,]$Ho_SE*sqrt(as.numeric(df.final[s,]$No.of.Loci))
s <- is.na(df.final$He_SD) & !is.na(df.final$He_SE)
df.final[s,]$He_SD <- as.numeric(df.final[s,]$He_SE)*sqrt(as.numeric(df.final[s,]$No.of.Loci))
s <- is.na(df.final$n.allele_SD) & !is.na(df.final$n.allele_SE)
df.final[s,]$n.allele_SD <- df.final[s,]$n.allele_SE*sqrt(as.numeric(df.final[s,]$No.of.Loci))

####################
df.final$No.of.Loci <- as.numeric(df.final$No.of.Loci)
df.final$per_polylocus <- as.numeric(df.final$per_polylocus)

# Filter data
df.lowss <- df.final %>% filter(No.of.Loci >= 5, sample.size >= 5)
print(list(
    "Some populations from these references were removed due to low No.of.Loci and Sample.size:", 
    setdiff(df.final$ref, df.lowss$ref)))

df.nocap <- df.lowss %>% filter(!grepl("captive",country,ignore.case=TRUE),!grepl("captive",population,ignore.case=TRUE), !grepl(",",country,ignore.case=TRUE),!is.na(country))
print(list("Captive populations and multiple countries were removed:", setdiff(df.lowss$ref, df.nocap$ref)))

df.nomix <- df.nocap %>%filter(!grepl("\\s", Species.specific))
print(list("Populations using a combination of FCA and other species derived loci were removed:", setdiff(df.nocap$ref, df.nomix$ref)))

df.lowper <- df.nomix[-which(df.nomix$per_polylocus < 0.5),]

print(list("Some populations with low % polymorphic loci were removed:", setdiff(df.nomix$ref, df.lowper$ref)))

print(paste("From", n_distinct(dfll$ref),"ref to", n_distinct(df.nomix$ref)))

# df.poly <- df.nomix %>% filter(per_polylocus > 0.5, is.na(per_polylocus))
df.lowper$Species.specific <- recode(df.lowper$Species.specific, 
    "F;FCA" = "FCA",
    "FCA;F" = "FCA")


gen <- df.lowper
colnames(gen) <- tolower(colnames(gen))
colnames(gen)[c(5,6)] <- c("lat","long")

gen <- gen[-grep("intro|hist", gen$population), ]

gen$ho[which(gen$ho == 1)] <- NA


filename_filtered <- "/home/amira.azizan/Documents/Y320/GenDivDatabase/Data/Gendiv_filtered_pop_averaged_v2.csv"
write.table(gen, filename_filtered, sep="\t", row.names=FALSE)




#
#gen <- read.csv("Data/Gendiv_filtered_pop_averaged_v2.csv", stringsAsFactors = F, sep="\t")
gen <- gen[-grep("zoo", gen$population, ignore.case = T),]
#summary(gen)
#135
gen$refy <- as.numeric(sapply(strsplit(gen$ref, "_"), function(x) paste(gsub("[A-z]","", x[3]))))

# Density data -----------------------------------------------------------
den <- read.csv("Data/FromStef/alldendata.csv")
# Add lion's data
lion_df <- read.csv("Data/FromStef/Lion Dataset.csv")
lion_df$Sname <- "Panthera leo"
names(lion_df)
lion_df <- rename(lion_df, "FirstAuthor"='Author', 
                  "YearP" = "yearP",
                  "SS"="Study.Site",
                  "N.Y.lat"="Lat",
                  "E.X.long"="Long",
                  "Years.Study"="yearStudy",
                  "TD"="method",
                  "Area.SS"="Area_SS"
)
den <- den[,names(den) %in% c("Sname","FirstAuthor","Journal","YearP", "Years.Study", "Country","SS","N.Y.lat","E.X.long","Area.SS","DP","Period",
                              "TD")]
den <- rbind(den, lion_df[,names(lion_df) %in% names(den)])
names(den)
colnames(den)[c(1,6:8)] <- c("sp","population","lat","long")

colnames(den) <- tolower(colnames(den))
den$ref <- paste(den$firstauthor,den$journal,den$yearp, sep="_")

den$dp <- as.numeric(den$dp)

# check country name coherence
sort(setdiff(den$country, gen$country))
sort(setdiff(gen$country, den$country))
den$country <- den$country %>% 
  recode("Myanmar (Burma)" = "Myanmar", 
         "Tanzania. United Republic of" = "United Republic of Tanzania",
         "United States" = "USA",
         "South Africa" = "Republic of South Africa",
         "Cameron" = "Cameroon",
         "Croatia-Slovakia" = "Slovakia",
         "United Kingdom" = "Scotland")

#saveRDS(den,"Data/Dendata_corrected.rds") #-> to be used in InteractiveMap


# Simplified df
colneeds <- c("sp","ref","country","lat","long","population")
g_df <- gen[,names(gen) %in% colneeds]
g_df$source <- "Popgen"


p_df <- den[,colneeds]
p_df$source <- "Popden"

names(g_df)

gp_df <- rbind(g_df, p_df)
gp_df$ID <- rownames(gp_df)

# Cluster dataset -----------------------------------------------------------
# 1. Remove pop without geographic coordinates
gp_df <- gp_df[!is.na(gp_df$long),]

# 2. Get biogeographic regions 
df <- gp_df
biogeo <- readRDS("Data/biogeo.rds")
points <- sp::SpatialPoints(coords=df[,c("long","lat")], proj4string = sp::CRS("+proj=longlat +datum=WGS84 +no_defs"))
biogeo.point <- raster::extract(x=biogeo[,c("ECO_NAME","REALM","BIOME")], y=points)

df <- cbind(df, biogeo.point[,c("ECO_NAME","REALM","BIOME")])

gp_df <- df

save(gen, den, gp_df,file="Data/GPall_data.rda")
