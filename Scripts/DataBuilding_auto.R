################ Title: Data building script
library(tidyr)
library(dplyr)

# load from DataBuilding_note.Rmd
sp_name <- sort(unique(gp_df$sp))
gp_full <- data.frame()

for (j in sp_name){
  sp_df.i <- droplevels(gp_df[gp_df$sp %in% j,])
  print(j)
  
  # - remove duplicated coordinates as it requires one point per localisation for clustering
  
  dups <- duplicated(sp_df.i[,c("long","lat","source")])
  if (any(dups)) {
    print(paste(length((dups[dups])), "Dups detected"))
    sp_df <- sp_df.i[!dups,]
    sp_df.ii <- sp_df.i[dups,]
  } else sp_df <- sp_df.i
  
  
  # - initialise column
  sp_df$cluster_den <- NA
  sp_df$cluster_gen <- NA
  sp_df$closest_ID <- NA
  sp_df$distance <- NA
  sp_df$cluster_center <- NA
  
  nbiogeo <- table(sp_df$REALM)
  nbiogeo
  
  # Apply clustering within each biogeo -----------------------------------------------------------------------
  table(sp_df$REALM, sp_df$source)
  
  # stop if there's only one kind of data (Popden/Popgen)
  if(length(unique(sp_df$source)) == 1) next
  
  # Which biogeo has both GD and PD
  kk <- tapply(sp_df$source, sp_df$REALM, n_distinct)
  sp_df <- subset(sp_df, REALM %in% names(which(kk == 2)))
  sp_df$REALM <- droplevels(sp_df$REALM)
  
  for (k in unique(sp_df$REALM)){
    
    bg_i <- which(sp_df$REALM == k)
    
    # Check number of pop in each dataset
    foo <- table(sp_df$source[bg_i])
    print(k)
    print(foo)
    # using Popden coordinates as centers
    cls_df <- sp_df[bg_i,]
    den_cls <- kmeans(cls_df[,c("long","lat")], centers=cls_df[which(cls_df$source == "Popden"),c("long","lat")], trace=TRUE)
    
    cls_df[,"cluster_den"] <- den_cls$cluster
    cls_df$cluster_den <- paste0(cls_df$cluster_den,"_", cls_df$REALM)
    
    # using Popgen coordinates as centers
    gen_cls <- kmeans(cls_df[,c("long","lat")], centers=cls_df[which(cls_df$source == "Popgen"),c("long","lat")], trace=TRUE)
    cls_df[,"cluster_gen"] <- gen_cls$cluster
    cls_df$cluster_gen <- paste0(cls_df$cluster_gen,"_",cls_df$REALM)
    
    sp_df[bg_i, "cluster_den"] <- cls_df$cluster_den
    sp_df[bg_i, "cluster_gen"] <- cls_df$cluster_gen
    
    # join unassigned data points to the nearest unassigned point
    ## based on gen
    bg_df <- sp_df[bg_i, ]
    # - calculate distances between the nearest points -Popgen <-> Popden, (in meters)
    xx <- bg_df[which(bg_df$source == "Popgen"),c("long","lat")]
    yy <- bg_df[which(bg_df$source == "Popden"),c("long","lat")]
    mat_gen2den <- geosphere::distm(xx, 
                                    yy)
    # get distance
    bg_df[which(bg_df$source == "Popgen"),"distance"] <- apply(mat_gen2den, 1, min)
    bg_df[which(bg_df$source == "Popden"),"distance"] <- apply(mat_gen2den, 2, min)
    
    # get id of the closest and store it in bg_df$closest_ID
    bg_df$closest_ID[which(bg_df$source == "Popden")] <- bg_df$ID[which(bg_df$source == "Popgen")][apply(mat_gen2den, 2, which.min)]
    bg_df$closest_ID[which(bg_df$source == "Popgen")] <- bg_df$ID[which(bg_df$source == "Popden")][apply(mat_gen2den, 1, which.min)]
    
    # - assign to the nearest cluster (not given by kmeans)
    # sort df in ascending order of the distance
    bg_df <- bg_df[order(bg_df$distance),]
    # by row operation
    # modify clus_den
    for (i in 1:nrow(bg_df)){
      x <- bg_df[i, "ID"]
      # refresh cluster info
      b <- xtabs(~cluster_den+source, data=bg_df)
      nclus.d <- which(rowSums(b > 0) == 2)
      
      clos_id <- bg_df[bg_df$ID == x,"closest_ID"]
      if(!bg_df$cluster_den[bg_df$ID == clos_id] %in% names(nclus.d)){
        
        bg_df$cluster_den[bg_df$ID == x] <- bg_df[bg_df$ID == clos_id,"cluster_den"]
      }
    }
    # modify clus_gen
    for (i in 1:nrow(bg_df)){
      x <- bg_df[i, "ID"]
      # refresh cluster info
      b <- xtabs(~cluster_gen+source, data=bg_df)
      nclus.g <- which(rowSums(b > 0) == 2)
      clos_id <- bg_df[bg_df$ID == x,"closest_ID"]
      if(!bg_df$cluster_gen[bg_df$ID == clos_id] %in% names(nclus.g)){
        bg_df$cluster_gen[bg_df$ID == x] <- bg_df[bg_df$ID == clos_id,"cluster_gen"]
      }
    }
    
    sp_df[bg_i,] <- bg_df
  }
  
  
  
  # - insert back dups (based on long/lat)
  if(any(dups)){
    poo <- merge(sp_df.ii[,c("long","lat","source")], 
                 sp_df[,c("long","lat","source","cluster_den","cluster_gen","closest_ID","distance","cluster_center")], all.x=T, all.y=F)
    poop <- merge(sp_df.ii, unique(poo))
    sp_df <- rbind(sp_df, poop)
  }
  
  if(nrow(sp_df) == 0) next
  # Number of data in each cluster
    # based on genetic pop
  if (any(is.na(sp_df$cluster_gen))) sp_df <- sp_df[!is.na(sp_df$cluster_gen),]
  a <- xtabs(~cluster_gen+source, data=sp_df)
  nclus.g <- which(rowSums(a > 0) == 2)
  a[rownames(a) %in% names(nclus.g),]
  
  clus_data <- sp_df[sp_df$cluster_gen %in% names(nclus.g),]
  # average distance across clusters
  summary(clus_data$distance)
  
  # based on density pop
  b <- xtabs(~cluster_den + source, data=sp_df)
  nclus.d <- which(rowSums(b > 0) == 2)
  b[rownames(b) %in% names(nclus.d),]
  clus_data <- sp_df[sp_df$cluster_den %in% names(nclus.d),]
  
  
  # average distance across clusters
  summary(clus_data$distance)
  
  if(length(nclus.g) > length(nclus.d))
  {
    sp_df$cluster_center <- sp_df$cluster_gen
  } else if(abs(length(nclus.g) - length(nclus.d)) <= 1) {
    bb <- table(sp_df$source)
    if(names(which.max(bb))=="Popgen") sp_df$cluster_center <- sp_df$cluster_gen else sp_df$cluster_center <- sp_df$cluster_den
  } else sp_df$cluster_center <- sp_df$cluster_den
  
  gp_full <- rbind(gp_full, sp_df)
  
}


