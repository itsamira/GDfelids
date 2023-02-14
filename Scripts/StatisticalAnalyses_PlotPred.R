library(ggeffects)
library(ggrepel)
library(ggpubr)
library(dplyr)
set.seed(1703)

no_green = c("#00007F", 
             "#0000B2", "#0000E5", "#0019FF", "#004DFF", "#007FFF", "#00B2FF", 
             "#00E5FF", "#7FFFFF", "#FFFFA5", "#FFFF54", 
             "#FFE500", "#FFB300", "#FF7F00", "#FF4C00", "#FF1900", 
             "#B20000", "#660000")


scales::show_col(no_green)
# reorder based on generation length
boo <- unique(gen_count[,c("sp","GenerationLength_d","BS")])
#boo <- boo[order(boo$GenerationLength_d),]
#gen_count$sp <- factor(gen_count$sp, levels = boo$sp)

no_green <- setNames(no_green, levels(gen_count$sp))
theme_set(theme_classic())
# source("03_ModelGenDen)

# Cluster data set
pred <- ggpredict(Ho.lmer.out, type="fixed", terms="DP_log[all]")
df_gen <- df_Ho[which(!df_Ho$ho_out & !df_Ho$dp_out),]
df_gen <- df_Ho
df_mean <-  df_gen %>% group_by(sp) %>% 
  dplyr::summarise(y = median(ho), x = log10(median(DP_mean)))

miny <- -0.8
maxy <- 2

df_gen$outliers <- df_gen$dp_out | df_gen$ho_out

Ho.cl <- ggplot() +
  geom_point(data=df_gen, aes(x= DP_log, y=ho, size=log10(sample.size), shape=outliers), 
             col="grey", show.legend = F) + #all data
  geom_point(data=df_mean, aes(x= x, y=y, fill=sp), size=6, pch=22) + #median/sp 
  scale_fill_manual(values = no_green ) +
  scale_shape_manual(values=c(1,8)) + 
  geom_line(data = pred, aes(x = x, y = predicted), col="black", size=1) +
  scale_y_continuous(limits=c(0,1)) + scale_x_continuous(limits=c(miny, maxy)) + 
  theme(legend.position = "none") +  
  labs(x = "", y = "Observed heterozygosity", title = "A)") 
  
# extreme values = low GD and DP
ext <- df_gen[df_gen$DP_log %in% head(sort(df_gen$DP_log),5),]
ext <- ext[!ext$outliers,]

ext$label_1 <- c("c", "d", "b", "e", "a")

Ho.cl <- Ho.cl + geom_point(data=ext, aes(x=DP_log, y=ho, col=sp), size=4) + scale_color_manual(values=no_green) + 
   geom_text_repel(data = ext, aes(x = DP_log, y = ho, label=label_1), 
                   size=3.5, direction="both", nudge_x=-0.1, nudge_y=0, 
                   hjust=0, parse=F, box.padding = 0, max.overlaps = Inf, min.segment.length=0)
Ho.cl

# 
# 
# add coeff
mod <- Ho.lmer.out
b <- round(unique(coef(mod)$sp$DP_log),2)
ci <- confint(mod, method="boot", n=1000)
ci <- as.data.frame(ci)[4,]
pv <- anova(mod)
pv <- paste(round(pv[6], 3))

cof <- paste0("==",b, "~(",round(ci[1], 2), "-",round(ci[2], 2),")")
Ho.cl.p <- Ho.cl + 
  geom_text(aes(x =1.25, y = 1), 
            label = paste(expression(beta[PD]), as.character(cof), "\n~ test"), 
            parse = TRUE,
            col="black",
            size=4) +
  
  geom_text(aes(x =1.25, y = 0.95), 
            label=paste0("p-value = ",pv),
            parse = F,
            col="black",
            size=4)

## He

pred <- ggpredict(He.lmer.out, type="fixed", terms="DP_log[all]")
df_gen <- df_He[!df_He$he_out & !df_He$dp_out,]
df_gen <- df_He

df_mean <-  df_gen %>% group_by(sp) %>% 
  dplyr::summarise(y = median(he), x = log10(median(DP_mean)))

df_gen$outliers <- df_gen$dp_out | df_gen$he_out

He.cl <- ggplot() +
  geom_point(data=df_gen, aes(x= DP_log, y=he, size=log10(sample.size), shape=outliers), 
             col="grey", show.legend = F) + #all data
  geom_point(data=df_mean, aes(x= x, y=y, fill=sp), size=6, pch=22) + #median/sp 
  scale_fill_manual(values = no_green ) +
  scale_shape_manual(values=c(1,8)) + 
  geom_line(data = pred, aes(x = x, y = predicted), col="black", size=1) +
  scale_y_continuous(limits=c(0,1)) + scale_x_continuous(limits=c(miny, maxy)) + 
  theme(legend.position = "none") +  labs(x = "", y = "Expected heterozygosity", title = "B)") 

# extreme values = low GD and DP
ext <- df_gen[df_gen$DP_log %in% head(sort(df_gen$DP_log),5),]
ext <- ext[!ext$outliers,]
ext$population.gen[grep("Rio", ext$population.gen)] <- "South Brazil"
ext$label_1 <- c("c", "d", "b", "e", "a")

He.cl <- He.cl + geom_point(data=ext, aes(x=DP_log, y=he, col=sp),size=4) + scale_color_manual(values=no_green) + 
   geom_text_repel(data = ext, aes(x = DP_log, y = he, label=label_1), 
                   size=3.5, direction="both", nudge_x=-0.1, nudge_y=0, 
                   hjust=0, parse=F, box.padding = 0, max.overlaps = Inf, min.segment.length=0)


# add coeff
mod <- He.lmer.out
b <- round(unique(coef(mod)$sp$DP_log),2)
ci <- confint(mod, method="boot", n=1000)
ci <- as.data.frame(ci)[4,]
pv <- anova(mod)
pv <- if(paste(round(pv[6], 1)) == "0") paste("< 0.01") else paste("=",round(pv[6], 1))

cof <- paste0("==",b, "~(",round(ci[1], 2), "-",round(ci[2], 2),")")
He.cl.p <- He.cl + 
  geom_text(aes(x =1.25, y = 1), 
            label = paste(expression(beta[PD]), as.character(cof), "\n~ test"), 
            parse = TRUE,
            col="black",
            size=4) +
  
  geom_text(aes(x =1.25, y = 0.95), 
            label=paste0("p-value ",pv),
            parse = F,
            col="black",
            size=4)

He.cl
## AR
pred <- ggpredict(mna.lmer.out, type="fixed", terms="DP_log[all]")
df_gen <- df_mna[!df_mna$ar_out & !df_mna$dp_out,]



df_gen <- df_mna

df_gen$outliers <- df_gen$dp_out | df_gen$ar_out

df_mean <-  df_gen %>% group_by(sp) %>% 
  dplyr::summarise(y = log10(median(n.allele)), x = log10(median(DP_mean)))

mna.cl <- ggplot() +
  geom_point(data=df_gen, aes(x= DP_log, y=n.allele_log, size=log10(sample.size), shape=outliers), col="grey", show.legend = F) + #all data
  geom_point(data=df_mean, aes(x= x, y=y, fill=sp), size=6, pch=22) + #median/sp 
  scale_fill_manual(values = no_green ) +
  scale_shape_manual(values=c(1,8)) + 
  geom_line(data = pred, aes(x = x, y = predicted), col="black", size=1) +
  scale_y_continuous(limits=c(0,1.3)) + scale_x_continuous(limits=c(miny, maxy)) + 
  theme(legend.position = "none") +  labs(x=expression(log[10]~Population~density~(ind/100~km^2)), y = expression(log[10]~Allelic~richness), title = "C)") 

# extreme values = low GD and DP
ext <- df_gen[df_gen$DP_log %in% head(sort(df_gen$DP_log),5),]
ext <- ext[!ext$outliers,]
ext$population.gen[grep("Rio", ext$population.gen)] <- "South Brazil"
ext$label_1 <- c("c", "d", "e", "a")

mna.cl <- mna.cl + geom_point(data=ext, aes(x=DP_log, y=n.allele_log, col=sp),size=4) + scale_color_manual(values=no_green) + 
   geom_text_repel(data = ext, aes(x = DP_log, y = n.allele_log, label=label_1), 
                   size=3.5, direction="both", nudge_x=-0.2, nudge_y=0.01, 
                   hjust=0, parse=F, box.padding = 0, max.overlaps = Inf, min.segment.length=0)
mna.cl
# 

# add coeff
mod <- mna.lmer.out
b <- round(unique(coef(mod)$sp$DP_log),2)
ci <- confint(mod, method="boot", n=1000)
ci <- as.data.frame(ci)[4,]
pv <- anova(mod)
pv <- paste(round(pv[6], 3))

cof <- paste0("==",b, "~(",round(ci[1], 2), "-",round(ci[2], 2),")")

mna.cl.p <- mna.cl + 
  geom_text(aes(x =1.25, y = 1.3), 
            label = paste(expression(beta[PD]), as.character(cof), "\n~ test"), 
            parse = TRUE,
            col="black",
            size=4) +
  
  geom_text(aes(x =1.25, y = 1.2), 
            label=paste0("p-value = ",pv),
            parse = F,
            col="black",
            size=4)


###
library(gridExtra)
library(grid)
get_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

############"
#' Country dataset

# Ho
df_gen <- df_Ho.c[!df_Ho.c$ho_out & !df_Ho.c$dp_out,]
df_gen <- df_Ho.c
df_gen$outliers <- df_gen$dp_out | df_gen$ho_out

df_mean <- df_gen %>% group_by(sp) %>% dplyr::summarise(y = median(ho), x = log10(median(DP)))

df_gen$labels <- gsub("_","\n", df_gen$SxC)

df_gen$spp <- paste0(substr(df_gen$sp, 1, 1),".", gsub("^.* ", "", df_gen$sp))

pred <- ggpredict(Ho.c.lmer.out, type="fixed", terms="DP_log[all]")

Ho.ct <- ggplot() +
  geom_point(data=df_gen, aes(x= DP_log, y=ho, size=log10(sample.size), shape=outliers), col="grey", show.legend = F) + #all data
  geom_point(data=df_mean, aes(x= x, y=y, fill=sp), size=6, pch=22) + #median/sp 
  scale_fill_manual(values = no_green ) +
  scale_shape_manual(values=c(1,8)) + 
  geom_line(data = pred, aes(x = x, y = predicted), col="black", size=1) +
  scale_y_continuous(limits=c(0,1)) + scale_x_continuous(limits=c(miny, maxy)) + 
  theme(legend.position = "none") +  labs(x = "", y = "Observed heterozygosity", title = "D)") 


# add coeff
mod <- Ho.c.lmer.out
b <- round(unique(coef(mod)$sp$DP_log),2)
ci <- confint(mod, method="boot", n=1000)
ci <- as.data.frame(ci)[4,]
pv <- anova(mod)
pv <- paste(round(pv[6], 4))

cof <- paste0("==",b, "~(",round(ci[1], 2), "-",round(ci[2], 2),")")
Ho.ct.p <- Ho.ct + 
  geom_text(aes(x =1.25, y = 1), 
            label = paste(expression(beta[PD]), as.character(cof), "\n~ test"), 
            parse = TRUE,
            col="black",
            size=4) +
  
  geom_text(aes(x =1.25, y = 0.95), 
            label=paste0("p-value = ", pv),
            parse = F,
            col="black",
            size=4)

##He
df.gen <- df_He.c[!df_He.c$he_out & !df_He.c$dp_out,]

df_gen <- df_He.c
df_mean <- df_gen %>% group_by(sp) %>% dplyr::summarise(y = median(he, na.rm = T), x = log10(median(DP)))
mod <- He.c.lmer.out

df_gen$outliers <- df_gen$dp_out | df_gen$he_out

df_gen$labels <- gsub("_","\n", df_gen$SxC)

df_gen$spp <- paste0(substr(df_gen$sp, 1, 1),".", gsub("^.* ", "", df_gen$sp))

pred <- ggpredict(He.c.lmer.out, type="fixed", terms="DP_log[all]")

He.ct <- ggplot() +
  geom_point(data=df_gen, aes(x= DP_log, y=he, size=log10(sample.size), shape=outliers), col="grey", show.legend = F) + #all data
  geom_point(data=df_mean, aes(x= x, y=y, fill=sp), size=6, pch=22) + #median/sp 
  scale_fill_manual(values = no_green, name="species") +
  scale_shape_manual(values=c(1,8)) + 
  guides (fill= guide_legend (nrow=3, byrow= TRUE )) +
  geom_line(data = pred, aes(x = x, y = predicted), col="black", size=1) +
  scale_y_continuous(limits=c(0,1)) + scale_x_continuous(limits=c(miny, maxy)) + 
  theme(legend.position="bottom", legend.text=element_text(face="italic", size=11)) +  labs(x = "", y = "Expected heterozygosity", title = "E)") 
# add coeff

b <- round(unique(coef(mod)$sp$DP_log),2)
ci <- confint(mod, method="boot", n=1000)
ci <- as.data.frame(ci)[4,]
pv <- anova(mod)
pv <- paste(round(pv[6], 4))

cof <- paste0("==",b, "~(",round(ci[1], 2), "-",round(ci[2], 2),")")
He.ct.p <- He.ct + 
  geom_text(aes(x =1.25, y = 1), 
            label = paste(expression(beta[PD]), as.character(cof), "\n~ test"), 
            parse = TRUE,
            col="black",
            size=4) +
  
  geom_text(aes(x =1.25, y = 0.95), 
            label=paste0("p-value = ",pv),
            parse = F,
            col="black",
            size=4)



#AR
df_gen <- df_mna.c[!df_mna.c$ar_out & !df_mna.c$dp_out,]

df_gen <- df_mna.c
df_mean <- df_gen %>% group_by(sp) %>% dplyr::summarise(y = log10(median(n.allele, na.rm=T)), x = log10(median(DP)))
mod <- mna.c.lmer.out
df_gen$outliers <- df_gen$dp_out | df_gen$ar_out

df_gen$labels <- gsub("_","\n", df_gen$SxC)

df_gen$spp <- paste0(substr(df_gen$sp, 1, 1),".", gsub("^.* ", "", df_gen$sp))

pred <- ggpredict(mna.c.lmer.out, type="fixed", terms="DP_log[all]")

mna.ct <- ggplot() +
  geom_point(data=df_gen, aes(x= DP_log, y=n.allele_log, size=log10(sample.size), shape=outliers), col="grey", 
             show.legend = F) + #all data
  geom_point(data=df_mean, aes(x= x, y=y, fill=sp), size=6, pch=22) + #median/sp 
  scale_fill_manual(values = no_green ) +
  scale_shape_manual(values=c(1,8)) + 
  geom_line(data = pred, aes(x = x, y = predicted), col="black", size=1) +
  scale_y_continuous(limits=c(0,1.3)) + scale_x_continuous(limits=c(miny, maxy)) + 
  theme(legend.position="none") + 
  labs(x = expression(log[10]~Population~density~(ind/100~km^2)), y = expression(log[10]~Allelic~richness), title = "F)")


# add coeff

b <- round(unique(coef(mod)$sp$DP_log),2)
ci <- confint(mod, method="boot", n=1000)
ci <- as.data.frame(ci)[4,]
pv <- anova(mod)
pv <- paste(round(pv[6], 4))

cof <- paste0("==",b, "~(",round(ci[1], 2), "-",round(ci[2], 2),")")
mna.ct.p <- mna.ct + 
  geom_text(aes(x =1.25, y = 1.3), 
            label = paste(expression(beta[PD]), as.character(cof), "\n~ test"), 
            parse = TRUE,
            col="black",
            size=4) +
  
  geom_text(aes(x =1.25, y = 1.2), 
            label=paste0("p-value = ", pv),
            parse = F,
            col="black",
            size=4)


#
leg2 <- get_legend(He.ct.p)
tit <- 
pred2 <- grid.arrange(arrangeGrob(Ho.cl.p + theme(plot.title = element_text(size=15, face="bold"), 
                                                  axis.title = element_text(size=15), axis.text = element_text(size=10)), 
                                  Ho.ct.p + theme(plot.title = element_text(size=15, face="bold"), 
                                                  axis.title = element_text(size=15), axis.text = element_text(size=10)), 
                                  He.cl.p + theme(plot.title = element_text(size=15, face="bold"),
                                                  axis.title = element_text(size=15), axis.text = element_text(size=10)), 
                                  He.ct.p+ theme(legend.position="hide", plot.title = element_text(size=15, face="bold"),
                                                 axis.title = element_text(size=15), axis.text = element_text(size=10)) ,
                                  mna.cl.p + theme(plot.title = element_text(size=15, face="bold"),
                                                   axis.title = element_text(size=15), axis.text = element_text(size=10)), 
                                  mna.ct.p + theme(plot.title = element_text(size=15, face="bold"), 
                                                   axis.title = element_text(size=15), axis.text = element_text(size=10)), nrow=3, ncol=2, top="top"), 
                                            heights=c(15))


ggsave("Figures/fig2_GDandPD.png", plot=pred2, width=13.5, height=17)



grid.arrange(arrangeGrob(Ho.cl.p + theme(plot.title = element_text(size=15, face="bold")), 
                         He.cl.p + theme(plot.title = element_text(size=15, face="bold")), 
                         mna.cl.p + theme(plot.title = element_text(size=15, face="bold")), 
                         nrow=3, ncol=1, top="Cluster data set"), 
             arrangeGrob(
                         Ho.ct.p + theme(plot.title = element_text(size=15, face="bold")), 
                        
                         He.ct.p+ theme(legend.position="hide", 
                                        plot.title = element_text(size=15, face="bold")) ,
                        
                         mna.ct.p + theme(plot.title = element_text(size=15, face="bold")), 
                         nrow=3, ncol=1, top="Country data set"),
             heights=c(15))




















