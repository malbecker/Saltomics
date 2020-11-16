##################################################################
##                                                              ##
##                 Code and Analyses for                        ##
## Molecular mechanisms of local adaptation for salt-tolerance  ##
##                                                              ##
## Authors: M. Albecker, A Stuckert, C. Balakrishnan, M. McCoy  ##
##                                                              ##
##################################################################

## Load packages

library(sleuth)
library(blme)
library(nlme)
library(optimx)
library(reshape2)
library(brglm)
library(lme4)
library(ggplot2)
library(lattice)
library(data.table)
library(pcaExplorer)
library(DESeq2)
require(RColorBrewer)
library(tidyverse)
library(ggthemes)
library(apeglm)
require(survminer)
require(survival)

################################
##    Plotting Aesthetics     ##
################################

eggshape = c("0" = 19, "4" = 1)
col.loc = cols = c("inland"="darkseagreen","coastal" = "dodgerblue4")
col.loc_cap = cols1 = c("Inland"="darkseagreen", "Coastal" = "dodgerblue4")
PCAcols = c("Inland" = "darkseagreen4", "Coastal" = "dodgerblue4")
pop_shapes <- c("BL" = 21,"LO" = 22,"PWL" = 23,"WHF" = 24,"BOD" = 21, "CSI" = 24,"DQ" = 22, "LH" = 23)
loc_shapes <- c("inland"= 21,"coastal" = 22)
line_type <- c("4_inland" = "dashed",  "0_inland"= "solid", "4_coastal"="dashed", "0_coastal"="solid")
line_type_cap <- c("4_Inland" = "dashed",  "0_Inland"= "solid", "4_Coastal"="dashed", "0_Coastal"="solid")
salinity <- c("0" = 15, "4" = 16,"6"=17)
eggpop_shape <- c("4_WHF" =2,"0_WHF"=17, "4_PWL" =5,"0_PWL"=18, "4_CSI"=2, "0_CSI"=17, "4_LH" =5, "0_LH" =18, 
                  "4_BOD"=1, "0_BOD"=16, "4_DQ"  =0, "0_DQ"=15,  "4_BL" =1, "0_BL"=16 , "4_LO"=0,  "0_LO"= 15)
eggtadx = c("0_0" = "0 ppt", "0_6" = "6 ppt", "4_4" = "4 ppt", "4_6" = "6 ppt")
shapesite = c("0_Inland_BL"=21,"4_Inland_BL"=21,
              "0_Coastal_BOD"=21, "4_Coastal_BOD"=21,
              "4_Coastal_DQ"=22,"0_Coastal_DQ"=22,
              "0_Coastal_LH"=23,"4_Coastal_LH"=23,
              "0_Inland_LO"=22,"4_Inland_LO" = 22,
              "0_Inland_PWL" =23,"4_Inland_PWL" =23,
              "4_Inland_WHF"=24,"0_Inland_WHF"=24,
              "0_Coastal_CSI"=24,"4_Coastal_CSI"=24)
colsite = c("0_Inland_BL"="darkseagreen" ,"4_Inland_BL"="darkseagreen" ,
            "0_Coastal_BOD"="dodgerblue4" , "4_Coastal_BOD"="dodgerblue4" ,
            "4_Coastal_DQ"="dodgerblue4" ,"0_Coastal_DQ"="dodgerblue4" ,
            "0_Coastal_LH"="dodgerblue4" ,"4_Coastal_LH"="dodgerblue4" ,
            "0_Inland_LO"="darkseagreen" ,"4_Inland_LO" = "darkseagreen" ,
            "0_Inland_PWL" ="darkseagreen" ,"4_Inland_PWL" ="darkseagreen" ,
            "4_Inland_WHF"="darkseagreen" ,"0_Inland_WHF"="darkseagreen" ,
            "0_Coastal_CSI"="dodgerblue4" ,"4_Coastal_CSI"="dodgerblue4")

################################
## Tadpole Phenotype Analyses ##
################################

######## Load Datasets #########

# Set working directory
setwd("~/Documents/GitHub/Saltomics/data/")

# Egg developmental data 
egg = read.csv("Egg_devel.csv")

# Egg counts and hatching counts
hatch = read.csv("Egg_Hatch.csv")

# Physiology data on Day 6 of acclimation
tadphysio = read.csv("HC_physio.csv")

# Tadpole survivorship during acclimation period
tadacc1 = read.csv("Tadpole Acclim_2017.csv")

# Survival (tadacc in long form)
survdat2 = read.csv("survdat2.csv")


######## Embryo Development #########

# Data wrangling
egg1<-filter(egg, Age_adj<50) # Some eggs did not hatch - this removes all non-hatching eggs.
egg$ppt <- as.factor(egg$ppt)

# Is there a difference in egg development according to salinity or location?
a1 <- lmer(log(stage) ~ Age_adj * ppt * Loc + (1|Pop), data = egg1)
a2 <- lmer(log(stage) ~ Age_adj + ppt + Loc + (1|Pop), data = egg1)
a3 <- lmer(log(stage) ~ 1 + ppt + Loc + (1|Pop), data = egg1)
a4 <- lmer(log(stage) ~ Age_adj + 1 + Loc + (1|Pop), data = egg1)
a5 <- lmer(log(stage) ~ Age_adj + ppt + 1 + (1|Pop), data = egg1)
anova(a1,a2) # Additive model is best
anova(a2,a3) # Yes, differences according to age
anova(a2,a4) # No difference according to salinity
anova(a2,a5) # No differences according to location

# Plot (not included in paper)
ggplot(egg1,aes(x=Age_adj, y=stage, group= factor(Loc), colour =factor(Loc)))+
  xlab("Age (hr)") +
  ylab("Gosner Stage") +
  geom_point()+
  stat_summary(position = position_dodge(0.3), geom = "line")+
  scale_fill_manual(values = cols1) +
  scale_colour_manual(values = cols1) +
  theme_hc(base_size = 24, base_family = "Times")+ 
  theme(axis.line = element_line(colour = "black"))+
  theme(axis.text= element_text(colour = "black")) +
  theme(legend.position="none")

######## Embryo Hatching #########

hatch$prop.hatch = hatch$Hatched/hatch$Eggs # Divide the number hatched by total laid in clutch

# How does salinity and location affect the proportion of eggs that hatch?
b1 = glmer(prop.hatch ~ factor(Egg_env) * Loc * (1|Pop), weights = Eggs, data = hatch, family = binomial(), na.action="na.exclude")
b2 = glmer(prop.hatch ~ factor(Egg_env) + Loc + (1|Pop),  weights = Eggs, data = hatch, family = binomial(), na.action="na.exclude")
b3 = glmer(prop.hatch ~ 1 * Loc + (1|Pop),  weights = Eggs, data = hatch, family = binomial(), na.action="na.exclude")
b4 = glmer(prop.hatch ~ factor(Egg_env) * 1 + (1|Pop), weights = Eggs, data = hatch, family = binomial(), na.action="na.exclude")

anova(b1,b2) # Significant interaction with salinity and location
anova(b1,b3) # Yes, significant effect of salinity
anova(b1,b4) # Yes, significant effect of location

# Plot (Figure 2a)
hatch$predicted = plogis(predict(b1))
hatch$eggloc = paste(hatch$Egg_env,hatch$Loc, sep="_")
hatch$eggpop = paste(hatch$Egg_env,hatch$Pop, sep="_")

ggplot(hatch,aes(x=Egg_env, y=prop.hatch, group= factor(eggloc)))+
  xlab("Salinity of Embryonic Environment (ppt)") +
  ylab("Proportion Embryos Hatched") +
  geom_boxplot(aes(fill=Loc),alpha = 0.4)+
  geom_jitter(aes(shape = Loc, fill = Loc),size = 4,  width = 0.35)+
  scale_fill_manual(values = col.loc) +
  scale_colour_manual(values = col.loc) +
  scale_shape_manual(values = loc_shapes) +
  theme_hc(base_size = 24, base_family = "Times")+ 
  theme(axis.line = element_line(colour = "black"))+
  theme(axis.text= element_text(colour = "black"))+
  theme(legend.position="none")

######## Tadpole Survival on Day 6 #########

# Data wrangling
tadacc = tadacc1[-c(475:480),] #removes suspicious outlier inland clutch (uncertain of cause)
tadacc$eggsal = as.factor(tadacc$eggsal)
tadacc$target_sal = as.factor(tadacc$target_sal)
tadacc$prop.surv = (tadacc$survive/50) # Divides number in cup by starting density of 50 tads
tadacc6 = filter(tadacc, day == "6")
tadacc6$weight = rep(50, nrow(tadacc6))

# Are there differences in overall tadpole survival on day 6? 

c1 = glmer(prop.surv ~ eggsal * loc * target_sal + (1|pop), # Changed optimizer to improve model convergence
            data=tadacc6, family = binomial(), weights = weight, glmerControl(optimizer = "bobyqa"),na.action="na.exclude") 
c2 = glmer(prop.surv ~ eggsal + loc + target_sal + (1|pop), 
             data=tadacc6, family = binomial(), weights = weight, glmerControl(optimizer = "bobyqa"),na.action="na.exclude")
c3 = glmer(prop.surv ~ 1 * loc * target_sal + (1|pop), 
            data=tadacc6, family = binomial(),  weights = weight, glmerControl(optimizer = "bobyqa"),na.action="na.exclude")
c4 = glmer(prop.surv ~ eggsal * 1 * target_sal + (1|pop), 
            data=tadacc6, family = binomial(),  weights = weight, glmerControl(optimizer = "bobyqa"),na.action="na.exclude")
c5 = glmer(prop.surv ~ eggsal * loc * 1 + (1|pop), 
            data=tadacc6, family = binomial(),  weights = weight, glmerControl(optimizer = "bobyqa"),na.action="na.exclude")
anova(c1,c2,test="Chisq") # Yes, significant interaction
anova(c1,c3,test="Chisq") # Significant impact of egg salinity
anova(c1,c4,test="Chisq") # Significant impact of location
anova(c1,c5,test="Chisq") # Significant impact of tadpole salinity

# Plot (Figure 2b)
tadacc6$predicted = plogis(predict(c1))
tadacc6$eggloc = paste(tadacc6$eggsal,tadacc6$loc, sep="_")
tadacc6$eggpop = paste(tadacc6$eggsal,tadacc6$pop, sep="_")

ggplot(tadacc6,aes(x=factor(target_sal), y=prop.surv, group= factor(eggloc)))+
  xlab("Target Salinity of Tadpole Environment (ppt)") +
  ylab("Proportion Survived") +
  geom_smooth(aes(y = predicted, colour = loc, linetype = eggloc), se = FALSE)+
  geom_jitter(aes(shape = eggpop,colour=loc, size = 2),width = 0.25)+
  #stat_summary(position = position_dodge(0.35),size=1.5)+
  scale_fill_manual(values = col.loc) +
  scale_colour_manual(values = col.loc) +
  scale_shape_manual(values = eggpop_shape) +
  scale_linetype_manual(values = line_type)+
  theme_hc(base_size = 24, base_family = "Times")+ 
  theme(axis.line = element_line(colour = "black"))+
  theme(axis.text= element_text(colour = "black"))+
  theme(legend.position="none")


######## Tadpole Survival Through Time #########

survdat2$eggsal = as.factor(survdat2$eggsal)
survdat2$target_sal = as.factor(survdat2$target_sal)
survdat2$dead = 0
for(i in 1:nrow(survdat2)){ #Create binary column - 1 = alive at MM, 0 = dead at MM
  if (survdat2$alive[i] == "1"){survdat2$dead[i] = 0}
  else (survdat2$dead[i] = 1)
}

# Kaplan-Meier Survival Estimates
surv1 <- survfit(Surv(day,dead) ~ eggsal + loc + target_sal, data = survdat2) 
sumdat = surv_summary(surv1)
sumdat$eggloc = paste(sumdat$eggsal,sumdat$loc,sep="_")

#Cox Regression (nonparametric hazard estimation)
cx1 <- coxph(Surv(day,dead) ~ eggsal + loc + target_sal, data = survdat2)
summary(cx1)
sumdat1 = filter(sumdat,target_sal=="12")

ggplot(data = sumdat1, aes(x=factor(time),y=surv,group = eggloc,colour = loc,shape=factor(eggsal)))+
  geom_line(position = position_dodge(0.3)) +
  geom_point(size=6,position = position_dodge(0.3)) +
  geom_errorbar(aes(ymin=lower,ymax=upper,colour = loc),width=0,position = position_dodge(0.3))+
  scale_fill_manual(values = cols) +
  scale_colour_manual(values = cols) +
  scale_shape_manual(values = eggshape) +
  #scale_linetype_manual(values = lineshift2)+ 
  theme_hc(base_size = 24, base_family = "Times")+ 
  theme(axis.line = element_line(colour = "black"))+
  theme(axis.text= element_text(colour = "black")) +
  xlab("Acclimation period (day)") +
  ylab("Survival in 12ppt (Target Salinity)")+
  theme(legend.position="none")

(summs = sumdat1 %>%
    group_by(eggsal) %>%
    summarize("mean" = mean(surv, na.rm = TRUE)))


######## Tadpole Length #########

tadphysio$EggEnv. = factor(tadphysio$EggEnv)
d1 = lmer(log(Length) ~ EggEnv. * TadEnv * Loc + (1|Cup) + (1|Pop), data = tadphysio, na.action="na.exclude")
d2 = lmer(log(Length) ~ EggEnv. + TadEnv + Loc + (1|Cup) + (1|Pop), data = tadphysio, na.action="na.exclude")
d3 = lmer(log(Length) ~ 1 + TadEnv + Loc + (1|Cup) + (1|Pop), data = tadphysio, na.action="na.exclude")
d4 = lmer(log(Length) ~ EggEnv. + TadEnv + 1 + (1|Cup) + (1|Pop), data = tadphysio, na.action="na.exclude")
d5 = lmer(log(Length) ~ EggEnv. + 1 + Loc + (1|Cup) + (1|Pop), data = tadphysio, na.action="na.exclude")

anova(d1,d2) # No interaction
anova(d2,d3) # No Effect of Embryonic Salinity 
anova(d2,d4) # Marginal effect of location (p = 0.06)
anova(d2,d5) # Effect of tadpole environment

# Plot (Figure 3a)
tadphysio$predicted = exp(predict(d2))
tadphysio$eggloc = paste(tadphysio$EggEnv,tadphysio$Loc, sep="_")
tadphysio$eggpop = paste(tadphysio$EggEnv,tadphysio$Pop, sep="_")

ggplot(tadphysio,aes(x=factor(TadEnv), y=Length, group= factor(eggloc)))+
  xlab("Target Salinity of Tadpole Environment (ppt)") +
  ylab("Total Length (mm)") +
  geom_jitter(aes(shape = eggpop, colour=Loc), size = 1, width = 0.25)+
  geom_smooth(aes(y = predicted, fill= Loc,colour = Loc, linetype = eggloc), se = TRUE)+
  scale_fill_manual(values = col.loc_cap) +
  scale_colour_manual(values = col.loc_cap) +
  scale_shape_manual(values = eggpop_shape) +
  scale_linetype_manual(values = line_type_cap)+
  theme_hc(base_size = 24, base_family = "Times")+ 
  theme(axis.line = element_line(colour = "black"))+
  theme(axis.text= element_text(colour = "black"))+
  theme(legend.position="none")


######## Tadpole Plasma Osmolality #########

# I had an issue with osmometer calibration, so I excluded populations that are unreliable 
goodpops = c("BOD","CSI","LH","PWL","WHF")
tadphysio1 = filter(tadphysio, Pop %in% goodpops) %>%
  droplevels()

# Are there differences in the osmolality of tadpoles?

tadphysio1$EggEnv = as.factor(tadphysio1$EggEnv)
f1 = glmer(Osmo ~ EggEnv * TadEnv * Loc + (1|Pop) + (1|Cup), data = tadphysio1, family = poisson(),
            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000),tol=0.001), na.action="na.exclude") # Difficulty converging
f2 = glmer(Osmo ~ EggEnv + TadEnv + Loc + (1|Pop) + (1|Cup), data = tadphysio1, family = poisson(), na.action="na.exclude")
f3 = glmer(Osmo ~ EggEnv + TadEnv + 1 + (1|Pop) + (1|Cup), data = tadphysio1, family = poisson(),
            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000),tol=0.001), na.action="na.exclude")
f4 = glmer(Osmo ~ TadEnv + 1 + Loc + (1|Pop) + (1|Cup), data = tadphysio1, family = poisson(), na.action="na.exclude")
f5 = glmer(Osmo ~ EggEnv + 1 + Loc + (1|Pop) + (1|Cup), data = tadphysio1, family = poisson(), na.action="na.exclude")
anova(f1,f2) # No significant interaction
anova(f2,f3) # Location effect  
anova(f2,f4) # Marginal egg environment effect (p = 0.09)
anova(f2,f5) # Effect of tadpole salinity (p > 0.01)

# Plot (Figure 3b)
tadphysio1$predicted = exp(predict(f2))
tadphysio1$eggloc = paste(tadphysio1$EggEnv,tadphysio1$Loc, sep="_")
tadphysio1$eggpop = paste(tadphysio1$EggEnv,tadphysio1$Pop, sep="_")

ggplot(tadphysio1,aes(x=factor(TadEnv), y=Osmo, group= factor(eggloc)))+
  xlab("Salinity of Tadpole Environment (ppt)") +
  ylab("Plasma Osmolality (mOsm/L)") +
  geom_jitter(aes(shape = eggpop, colour=Loc), size = 1, width = 0.25)+
  geom_smooth(aes(y = predicted, fill= Loc,colour = Loc, linetype = eggloc), se = TRUE)+
  scale_fill_manual(values = col.loc_cap) +
  scale_colour_manual(values = col.loc_cap) +
  scale_shape_manual(values = eggpop_shape) +
  scale_linetype_manual(values = line_type_cap)+
  theme_hc(base_size = 24, base_family = "Times")+ 
  theme(axis.line = element_line(colour = "black"))+
  theme(axis.text= element_text(colour = "black"))+
  theme(legend.position="none")


######################################
## Tadpole Gene Expression Analyses ##
######################################

## Load and format transcriptome data

# Establish the pathway to pull from
base_dir <- "~/Desktop/Work/DataSets/Saltomics" 

# Attach identifying information to folders and samples
sample_id <- dir(file.path(base_dir, "kallisto_quants"))

kal_dirs <- paste0("~/Desktop/Work/DataSets/Saltomics/kallisto_quants/",sample_id)
samples <- read.csv("~/Desktop/Work/DataSets/Saltomics/omic_samples.csv", header = TRUE) # Sample meta-data
samples$location = sub("inland","Inland",as.character(samples$location)) 
samples$location = sub("coastal","Coastal", as.character(samples$location))

samples <- samples[order(samples$sample),] #Aligns samples data with sample_id data

# Combine into a single combined factor
samples$path <- kal_dirs

# Creates the table from kallisto
so <- sleuth_prep(samples) 

#Pull out data and format for downstream applications
so1 = so[["obs_raw"]] # I want the not-normalized data
so2 = full_join(so1,samples,by = "sample")

# Prepare by building the design
samples$treatment <- paste(samples$egg, samples$tad, samples$location, sep = "_")
samples$eggtad <- paste(samples$egg,samples$tad, sep="_")
samples$sample. <- paste(samples$egg,samples$tad,samples$population,samples$ind,sep="_")

# Add information to dataset
so2$treatment <- paste(so2$egg, so2$tad, so2$location, sep = "_")
so2$eggtad <- paste(so2$egg,so2$tad, sep="_")

## Pull in Annovations

# Uniref Annotations
u1g <- fread("~/Desktop/Work/DataSets/Saltomics/uniprot-reviewed_yes.tab", header = TRUE,sep="\t")
u2g <- u1g[,c(1,2,5,4)]
colnames(u2g) <- c("entry", "entry_name", "gene_name","protein_id")

# Trim the gene name column 
new_genes = separate(u2g, gene_name, into = c("gene_name", "trash"), sep = "\\s")
new_genes$gene_name <- tolower(new_genes$gene_name) # lowercase

# merge 
a2g <- new_genes[,-4]

# Add full annotation first
ann1 <- read.table("~/Desktop/Work/DataSets/Saltomics/cinerea.annotation.txt", header = FALSE, fill = TRUE)
ann <- ann1[,c(1,3)]
colnames(ann) <- c("target_id","entry")

ann <- dplyr::left_join(ann, a2g, by = "entry")
ann$gene_name <- tolower(ann$gene_name)

# Full Datset with all annotations and identifiers
so_anno1 <- dplyr::full_join(ann,so2,by="target_id") 

# Restrict to just annotated genes
so_anno <- filter(so_anno1, gene_name != "NA") # Cuts down to 36% of dataset that is annotated

## FOR PCA -- ALL TRANSCRIPTS ##

## Create dataframe with counts
so_anno1$sample. = paste(so_anno1$egg,so_anno1$tad,so_anno1$population,so_anno1$ind,sep="_")
dat1 <- data.frame("id" = so_anno1$sample.,"target_id" = so_anno1$target_id,"count" = so_anno1$est_counts) # Transcript level

dat2 <- dat1 %>%
  group_by(id, target_id) %>%
  summarise(count_sum = sum(count, na.rm = TRUE))

dat2$count_sum = as.integer(dat2$count_sum)
dat_transcript = data.frame(spread(dat2, id, count_sum))

# For downstream idenfication
gene_names <- so_anno[,c(1,4)]
gene_names <- gene_names %>% unique()
dat3 = left_join(dat_transcript,gene_names, by = "target_id")
dat3$uniquetarget = paste(dat3$target_id,dat3$gene_name,sep="_")
dat4 = dat3 %>% separate(uniquetarget,
                         into = c("firsthalf", "secondhalf"), 
                         remove = FALSE)


dat3_transcript = dat3
rownames(dat3_transcript) <- dat3_transcript$uniquetarget
dat3_transcript = dat3_transcript[,c(-1,-33,-34)]
colnames(dat3_transcript) <- sub("X", "", colnames(dat3_transcript))
dim(dat3_transcript) # 71,016

## FOR PCA -- JUST ANNOTATED GENES ##

## Create dataframe with counts
so_anno$sample. = paste(so_anno$egg,so_anno$tad,so_anno$population,so_anno$ind,sep="_")
dat1. <- data.frame("id" = so_anno$sample.,"gene" = so_anno$target_id,"count" = so_anno$est_counts) # Transcript level
dat2. <- dat1. %>%
  group_by(id, gene) %>%
  summarise(count_sum = sum(count, na.rm = TRUE))
dat2.$count_sum = as.integer(dat2.$count_sum)
dat3_gene = data.frame(spread(dat2., id, count_sum))
rownames(dat3_gene) <- dat3_gene[,1]
dat3_gene = dat3_gene[,-1]
colnames(dat3_gene) <- sub("X", "", colnames(dat3_gene))

dim(dat3_gene) # 25,418

# To nest clutch in location - Values are assigned manually because its messy
samples$family = as.factor(c(1,1,2,3,2,3,4,1,1,2,4,3,2,3,4,1,1,2,3,4,2,3,4,1,1,2,3,4,1,2,3)) # 1:4 for each clutch for each location

#To make test variable dataset 
coldata1 <- data.frame("id" = samples$sample., 
                       "location" = samples$location, 
                       "egg" = factor(samples$egg),
                       "tad" = factor(samples$tad),
                       "family"= samples$family) 

rownames(coldata1) <- coldata1[,1]
coldata1 = coldata1[,-1]

# Test for matching
all(rownames(coldata1) == colnames(dat3_transcript)) # Usually False
all(rownames(coldata1) == colnames(dat3_gene)) # Usually False

#If false, use below code. 
dat3_transcript <- dat3_transcript[,rownames(coldata1)]
dat3_gene <- dat3_gene[,rownames(coldata1)]

#Effect of location on all transcripts (For PCA)
dds_pca <- DESeqDataSetFromMatrix(countData = dat3_transcript, 
                                  colData = coldata1,
                                  design = ~ location + egg + tad + location:family + location:egg + location:tad) 


#Effect of location on just annotated genes (For DeSeq2)
dds_loc <- DESeqDataSetFromMatrix(countData = dat3_gene, 
                                  colData = coldata1,
                                  design = ~ location + egg + tad + location:family + location:egg + location:tad) 

# For WGCNA (only on annotated genes)
dds_loc_tr = dds_loc

## Transform Data 
rl_loc <- DESeq2::vst(dds_pca)  

## Launch PCA Explorer
pcaExplorer(dds = dds_pca, rl_loc) #Will launch app

## Extract high loadings
pcaobj <- prcomp(t(SummarizedExperiment::assay(rl_loc)))
hi_loadings(pcaobj, whichpc = 2, topN = 10,exprTable=counts(dds_pca)) 

## Extract PC1 and PC2 data

pcaData <- plotPCA(rl_loc, intgroup = c("location","egg","tad"), returnData = TRUE)
pcaData$Population <- so_anno$population[match(pcaData$name,so_anno$sample.)]
percentVar <- round(100 * attr(pcaData, "percentVar"),digits=2)

## Plot Specs
pcaData$group = paste(pcaData$egg, pcaData$tad, pcaData$location,sep=":")
pcaData$trt = paste(pcaData$egg, pcaData$tad, sep = "_")

# Plot PCA (Figure 4)
ggplot(pcaData, aes(x = PC1, y = PC2, shape = Population)) +
  geom_point(aes(fill = location),size = 4) +
  stat_ellipse(aes(group = location, colour = location),size = 1,alpha = 1)+
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  scale_colour_manual(name = "Location",values = PCAcols)+
  scale_linetype_manual(name = "Salinity Exposure", 
                        labels = c("0ppt egg: 0ppt tad","0ppt egg: 6ppt tad","4ppt egg: 4ppt tad", "4ppt egg: 6ppt tad"), 
                        values = trt_lines)+
  scale_fill_manual(values= PCAcols)+ guides(fill = FALSE) +
  scale_shape_manual(values = pop_shapes,
                     labels = c("Inland 1 (BL)", "Coastal 1 (BOD)", "Coastal 2 (CSI)","Coastal 3 (DQ)","Coastal 4 (LH)","Inland 2 (LO)","Inland 3 (PWL)","Inland 4 (WHF)")) + #guides(shape = FALSE) +
  theme_bw(base_size = 24, base_family = "Times")+
  theme(axis.text.x = element_text(size=20,colour = "black",angle = 45, hjust=1),
        axis.title.x = element_text(size=24,face="bold")) +
  theme(axis.text.y = element_text(size=20,colour = "black"),
        axis.title.y = element_text(size=24,face="bold")) +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 2)) + theme(legend.position="none")

# To identify DE genes
gene_names <- so_anno[,c(1,4)]
gene_names <- gene_names %>% unique()
gene_names[gene_names$target_id == "NODE_1_length_16363_cov_4496.792504_g0_i0",] # To identify unnamed transcripts, use this.

## DESEQ ANALYSIS
levels(dds_loc$location) <- c("Coastal", "Inland")
dds_loc$location <- relevel(dds_loc$location, ref = "Inland") 
dds_loc$egg <- relevel(dds_loc$egg, ref = "0") 
dds_loc$tad <- relevel(dds_loc$tad, ref = "0")

de_loc <- DESeq(dds_loc) 
resultsNames(de_loc) # For Comparisons

# Location Comparison - Coastal vs. Inland
loc_shrink <- lfcShrink(de_loc, coef="location_Coastal_vs_Inland", type="apeglm")
sum(loc_shrink$padj < 0.05, na.rm=TRUE) # 633 DE genes
loc_go <- data.frame(loc_shrink[,-1])
loc_go <- rownames_to_column(loc_go, "target_id")
loc_go1 <- left_join(loc_go, gene_names, by = "target_id")
loc_go2 <- filter(loc_go1, padj <= 0.05)
loc_go3 <- arrange(loc_go2, -log2FoldChange) # 2nd row is HBB3, use above code to find other missing values

# Egg Salinity Comparison - 4ppt vs 0ppt
de_sal <- lfcShrink(de_loc, coef="egg_4_vs_0", type="apeglm")
sum(de_sal$padj < 0.05, na.rm=TRUE) # 9
de_sal <- data.frame(de_sal[,-1])
de_sal <- rownames_to_column(de_sal, "target_id")
de_sal <- left_join(de_sal, gene_names, by = "target_id")
de_sal <- filter(de_sal, padj <= 0.05)
de_sal <- arrange(de_sal, -log2FoldChange) # Use above code to fill in missing gene names

# Tadpole Salinity Comparison - 6ppt vs 0ppt
de_sal2 <- lfcShrink(de_loc, coef="tad_6_vs_0",type="apeglm")
sum(de_sal2$padj < 0.05, na.rm=TRUE) # 44
de_sal2 <- data.frame(de_sal2[,-1])
de_sal2 <- rownames_to_column(de_sal2, "target_id")
de_sal2 <- left_join(de_sal2, gene_names, by = "target_id")
de_sal2 <- filter(de_sal2, padj <= 0.05)
de_sal2 <- arrange(de_sal2, -log2FoldChange)

# Interaction Comparison - Location and Egg salinity
de_interaction <- lfcShrink(de_loc, coef="locationCoastal.egg4",type="apeglm")
sum(de_interaction$padj < 0.05, na.rm=TRUE) # 5
de_interaction <- data.frame(de_interaction[,-1])
de_interaction <- rownames_to_column(de_interaction, "target_id")
de_interaction <- left_join(de_interaction, gene_names, by = "target_id") 
de_interaction <- filter(de_interaction, padj <= 0.05)
de_interaction <- arrange(de_interaction, -log2FoldChange)

# Interaction Comparison - Location and Tadpole salinity
de_interaction2 <- lfcShrink(de_loc, coef="locationCoastal.tad6",type = "apeglm")
sum(de_interaction2$padj < 0.05, na.rm=TRUE) # 7
de_interaction2 <- data.frame(de_interaction2[,-1])
de_interaction2 <- rownames_to_column(de_interaction2, "target_id")
de_interaction2 <- left_join(de_interaction2, gene_names, by = "target_id")
de_interaction2 <- filter(de_interaction2, padj <= 0.05)
de_interaction2 <- arrange(de_interaction2, -log2FoldChange)

# Combine all into single dataset (Supplemental Material 1)
all.genes = rbind(loc_go3,de_sal, de_sal2,de_interaction,de_interaction2,fill=TRUE)
all.genes[all.genes$gene_name == "tjp3",] # To look at patterns of single genes

# Plot genes of interest 

genePlot =function(transcript_id,refdata){
  
  target_gene = transcript_id
  gene_name = unique(refdata$gene_name[refdata$target_id == target_gene])
  
  # Arrange dataset for plot
  b_temp <- plotCounts(de_loc, gene = target_gene, 
                       intgroup = c("location"), returnData = TRUE, xlab="treatment")
  setDT(b_temp, keep.rownames = TRUE)
  colnames(b_temp)[1]<-"sample."
  c_temp = left_join(b_temp,samples,by ="sample.")
  c_temp$eggtadloc = paste0(c_temp$egg,c_temp$tad,c_temp$location.x)
  c_temp$eggloc = paste(c_temp$egg,c_temp$location.x,sep = "_")
  c_temp$egglocpop = paste(c_temp$egg,c_temp$location.x,c_temp$population, sep = "_")
  
  # Plot Candidates
  (geneplot_candidate = ggplot(c_temp,  aes(x = eggtad, y = count)) + 
    ggtitle(gene_name,  subtitle = "") +
    ylab("Normalized Counts") + 
    xlab("Salinity of Tadpole Environment") +
    geom_vline(aes(xintercept = 2.5))+
    geom_boxplot(aes(fill= location.x),size = 0.75,colour = "black",alpha = 0.6)+
    geom_point(aes(group = egglocpop, shape = egglocpop, fill = location.x),position = position_dodge(0.4), size = 3)+
    scale_fill_manual(values = PCAcols)+ # Use colsite with fill = egglocpop when you want to recreate legend.
    scale_shape_manual(name = "Population",
                       values = shapesite,
                       labels = c("0_Inland_BL"="Inland 1 (BL)","4_Inland_BL"="Inland 1 (BL)",
                                  "0_Coastal_BOD"="Coastal 1 (BOD)", "4_Coastal_BOD"="Coastal 1 (BOD)",
                                  "4_Coastal_DQ"="Coastal 3 (DQ)","0_Coastal_DQ"="Coastal 3 (DQ)",
                                  "0_Coastal_LH"="Coastal 4 (LH)","4_Coastal_LH"="Coastal 4 (LH)",
                                  "0_Inland_LO"="Inland 2 (LO)","4_Inland_LO" = "Inland 2 (LO)",
                                  "0_Inland_PWL" ="Inland 3 (PWL)","4_Inland_PWL" ="Inland 3 (PWL)",
                                  "4_Inland_WHF"="Inland 4 (WHF)","0_Inland_WHF"="Inland 4 (WHF)",
                                  "0_Coastal_CSI"="Coastal 2 (CSI)","4_Coastal_CSI"="Coastal 2 (CSI)"))+
    scale_x_discrete(labels = eggtadx)+
    theme_hc(base_size = 22, base_family = "Times")+ 
    theme(axis.line = element_line(colour = "black"))+
    theme(axis.text= element_text(colour = "black")) + theme(legend.position="none")
  ) # Saved as 6x8 inch PDF
  return(geneplot_candidate)
}

# Ion Transporters
(atp1b1 = genePlot("S668054",so_anno))
(atp6v1g2 = genePlot("R666300",so_anno))
(nalcn = genePlot("NODE_36985_length_332_cov_2.131313_g33311_i0",so_anno))

# Cell Adhesion 
(cdh17 = genePlot("NODE_11566_length_1263_cov_5.210098_g9470_i0",so_anno))
(cdh23 = genePlot("TRINITY_DN10889_c0_g2_i1",so_anno))
(pkhd1l1 = genePlot("TRINITY_DN4134_c0_g4_i1",so_anno)) # Most downregulated transcript of several
(cdh26 = genePlot("TRINITY_DN4588_c0_g4_i1",so_anno))
(cldn1 = genePlot("NODE_7127_length_1796_cov_8.010789_g5786_i0",so_anno))
(ocln = genePlot("NODE_2513_length_3000_cov_4.461046_g2012_i0",so_anno))
(pkp3 = genePlot("TRINITY_DN20493_c0_g4_i1",so_anno))
(gjb3 = genePlot("NODE_35414_length_259_cov_2.225490_g29797_i0",so_anno))
(ppl = genePlot("NODE_9920_length_1231_cov_2.914966_g6925_i1",so_anno))
(eppk = genePlot("NODE_1284_length_3806_cov_9.718377_g1008_i0",so_anno))
(tjp3 = genePlot("NODE_8282_length_1405_cov_2.746667_g5806_i0",so_anno))

# Cytoskeleton
(odc1a = genePlot("NODE_13734_length_1071_cov_78.253861_g4759_i1",so_anno))
(oaz = genePlot("R671092",so_anno))
(tgm3 = genePlot("TRINITY_DN7644_c0_g1_i1",so_anno))

(col14a = genePlot("NODE_23809_length_460_cov_1.851852_g18355_i0",so_anno))
(actn4 = genePlot("TRINITY_DN1213_c0_g1_i2",so_anno))
(myo1e = genePlot("TRINITY_DN1969_c0_g1_i1",so_anno))
(lama3 =  genePlot("NODE_31240_length_416_cov_2.039370_g27584_i0",so_anno))
(arhgap8 = genePlot("TRINITY_DN11822_c0_g2_i1",so_anno))
(pkhd1l1 = genePlot("TRINITY_DN4134_c0_g4_i1",so_anno))
(cadhf = genePlot("NODE_1237_length_3404_cov_7.922365_g826_i0",so_anno))
(scel = genePlot("NODE_3453_length_2621_cov_5.753674_g2788_i0",so_anno))

# Glycerol
(gpd1 = genePlot("NODE_10087_length_1412_cov_1023.620915_g8225_i0",so_anno))

# Environmental exposure
(slc26a4 = genePlot("NODE_14770_length_999_cov_1.945021_g12203_i0",so_anno))
(aqp5 = genePlot("TRINITY_DN16791_c0_g1_i1",so_anno))
(atp6v1g3 = genePlot("NODE_18031_length_656_cov_21.474210_g13194_i0",so_anno))

###############
##   WGCNA   ##
###############

library(WGCNA)
library(flashClust)
library(nlme)

# VSD transformation on deseq object (above code for "rl")
rl_ge <- DESeq2::vst(dds_loc_tr) 
rl1 <- assay(rl_ge)
rl1 <- as.data.frame(rl1)
setDT(rl1, keep.rownames = TRUE)[]
colnames(rl1)[1] <- "xen_gene_name"

# Match up gene names with transformed dataset and trim
rl2 = as.data.frame(t(rl1[,-1]))
colnames(rl2) <- t(rl1[,1])

# Exclude genes with low or no variance
gsg = goodSamplesGenes(rl2, verbose = 3)
gsg$allOK #False means there are genes to be excluded. 

# If False, run function below to remove 
if (!gsg$allOK) 
{if (sum(!gsg$goodGenes)>0)
  printFlush(paste("Removing genes:", paste(names(rl2)[!gsg$goodGenes], collapse= ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(rl2)[!gsg$goodSamples], collapse=", ")))
  rl9= rl2[gsg$goodSamples, gsg$goodGenes] # Removed 109 genes
}

gsg = goodSamplesGenes(rl9, verbose = 3)
gsg$allOK # Should be true now


# Load trait data for correlations
traitdat <- data.frame("sample" = samples$sample., 
                       "location" = samples$location,
                       "egg" = samples$egg,
                       "tad" = samples$tad)
rownames(traitdat) <- traitdat[,1]
traitdat = traitdat[,-1]

traitdat$location = as.numeric(traitdat$location)
traitdat$egg = as.numeric(traitdat$egg)
traitdat$tad = as.numeric(traitdat$tad)
table(rownames(traitdat)==rownames(rl9)) #True, so both datasets line up correctly.

# Physiology data on Day 6 of acclimation
tadphysio$group = paste(tadphysio$EggEnv,tadphysio$TadEnv,tadphysio$Pop,tadphysio$Individual,sep= "_")
# Recall issue with calibration, so use dataset with reliable data (tadphysio1)

## Osmolality
tadphysio1$group = paste(tadphysio1$EggEnv,tadphysio1$TadEnv,tadphysio1$Pop,tadphysio1$Individual,sep= "_")
osmo_df <- tadphysio1 %>% 
  group_by(group) %>%
  summarize("osmo" = mean(Osmo,na.rm = TRUE))

## Length
length_df <- tadphysio %>% 
  group_by(group) %>%
  summarize("length" = mean(Length,na.rm = TRUE))

## Load phenotype data for correlations
osmo_df1 = data.frame(osmo_df)
length_df1 = data.frame(length_df)
phen_dat = merge(osmo_df1, length_df1)
traitdat$group <- rownames(traitdat)
trait_phen_dat = left_join(traitdat,phen_dat,by = "group") # reduces only to the samples that were used for WGCNA
rownames(trait_phen_dat) = trait_phen_dat$group
trait_phen_dat = trait_phen_dat[,-4]

## Calculate network connectivity
A = adjacency(t(rl9),type="signed") 
k = as.numeric(apply(A,2,sum))-1 # standardized connectivity
Z.k = scale(k)
thresholdZ.k = -2.5 # standard practice is -2.5
outlierColor = ifelse(Z.k<thresholdZ.k,"red","white") # establish colors
sampleTree = flashClust(as.dist(1-A), method = "average")

## Convert traits to a color representation where red indicates high values
traitColors = data.frame(numbers2colors(trait_phen_dat,signed=TRUE))
dimnames(trait_phen_dat)[[2]] = paste(names(trait_phen_dat))
datColors = data.frame(outlier = outlierColor,traitColors)

## Plot dendrogram
plotDendroAndColors(sampleTree,groupLabels=names(datColors),
                    colors=datColors,main="Sample Dendrogram and Trait Heatmap")

## Choose set of soft-threshold power
powers = c(c(1:10), seq(from =10, to=30, by=1)) #choosing a set of soft-thresholding powers

## Call network topology analysis function
sft = pickSoftThreshold(rl9, powerVector=powers, verbose = 5, 
                        networkType="signed",allowWGCNAThreads(nThreads = 3)) 

## Plot 
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], 
     xlab= "Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed R^2", 
     type= "n", main= paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], 
     labels=powers, col="red")
abline(h=0.90, col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab= "Soft Threshold (power)", ylab="Mean Connectivity", 
     type="n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, col="red")
dev.off()

## Define soft power and adjacency
softPower = 8 # first value that crosses soft threshold line.
adj = adjacency(rl9, power = softPower, type = "signed") #specify network type

## Translate adjacency into topological overlap matrix and calculate the corresponding dissimilarity

TOM = TOMsimilarity(adj, TOMType="signed") #Topological Overlap Matrix
dissTOM = 1-TOM 

## Generate a clustered gene tree

geneTree = flashClust(as.dist(dissTOM), method="average") 
#plot(geneTree, xlab="", sub="", main= "Gene Clustering on TOM-based dissimilarity", labels= FALSE, hang=0.04)

## Set the minimum number of genes to cluster into a module and cut the trees
minModuleSize = 30 
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize)
#load("~/Desktop/Work/Manuscripts/Saltomics/Evolution/WGCNA/dynamicMods.Rdata") # Load from Cluster
dynamicColors= labels2colors(dynamicMods) # Color assignments

## Assign genes eigenvalues to be sorted into different modules

MEList= moduleEigengenes(rl9, colors= dynamicColors, softPower = 8)
MEs= MEList$eigengenes
MEDiss= 1-cor(MEs)
METree= flashClust(as.dist(MEDiss), method= "average")

## Plot tree showing how the eigengenes cluster together
plot(METree, main= "Clustering of module eigengenes", xlab= "", sub= "") # Merge? 

## Set threhold for merging modules and merge
MEDissThres = 0.2
merge = mergeCloseModules(rl9, dynamicColors, cutHeight= MEDissThres, verbose =3)

## Establish new colors and module IDs
mergedColors = merge$colors
mergedMEs = merge$newMEs
dev.off()

## Plot original dendrogram with new, merged module colors below 
#plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), 
#                    dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang=0.05)

## Establish new color schemes and data 
moduleColors = mergedColors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs

## Define number of genes and samples
nGenes = ncol(rl9)
nSamples = nrow(rl9)

## Recalculate MEs with new color labels
MEs0 = moduleEigengenes(rl9, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

## Run correlations on different treatments (location, egg, tadpole)
moduleTraitCor = cor(MEs, trait_phen_dat, use= "p") # p = pearson
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

## Print correlation heatmap between modules and traits
textMatrix= paste(signif(moduleTraitCor, 2), "\n(", 
                  "p = ",signif(moduleTraitPvalue, 1), ")", sep= "")
dim(textMatrix)= dim(moduleTraitCor)
par(mar= c(6, 8.5, 3, 3))

## Display the correlation values with a heatmap plot
#pdf(file="heatmap.pdf")
labeledHeatmap(Matrix= moduleTraitCor, # all modules, not just significant ones
               xLabels= names(trait_phen_dat), 
               yLabels= names(MEs), 
               ySymbols= names(MEs), 
               colorLabels= FALSE, 
               colors= blueWhiteRed(50), 
               textMatrix= textMatrix, 
               setStdMargins= FALSE, 
               cex.text= 0.5, 
               zlim= c(-1,1), 
               main= paste("Module-trait relationships"))
dev.off()

## for-loop specs
rownum = nrow(moduleTraitCor)
colnum = ncol(moduleTraitCor)

## Exclude non-significant modules from heatmap
MEcor <- as.data.frame(matrix("NA", nrow = rownum, ncol = colnum), row.names = rownames(moduleTraitPvalue))
colnames(MEcor) <- colnames(moduleTraitCor)
MEcor[,1] <- as.character(unlist(MEcor[,1]))
MEcor[,2] <- as.character(unlist(MEcor[,2]))
MEcor[,3] <- as.character(unlist(MEcor[,3]))
MEcor[,4] <- as.character(unlist(MEcor[,4]))
MEcor[,5] <- as.character(unlist(MEcor[,5]))

for (i in 1:rownum){
  for (j in 1:colnum){
    MEcor[i,j] <- if(moduleTraitPvalue[i,j] <= 0.055) {signif(moduleTraitCor[i,j], 2)}else{as.character("0")}
  }
}
MEcor = data.matrix(MEcor)

## Heatmap with just significant values 
#pdf(file="significant_heatmap.pdf")
labeledHeatmap(Matrix = MEcor,
               xLabels = names(trait_phen_dat),
               plotLegend = FALSE, 
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(100),
               textMatrix = textMatrix ,
               setStdMargins = FALSE,
               cex.text = 0.8,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

dev.off()
sig_mod = moduleColors 

## Module Membership pvalues
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(rl9, MEs, use = "p"));
MMPvalue=as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM.", modNames, sep="");
geneModuleMembership = geneModuleMembership %>%
  rownames_to_column() %>%
  rename("rowname" = "target_id")
names(MMPvalue) = paste("p.MM", modNames, sep="");

preGO <- data.frame("target_id" = so_anno$target_id, "gene_name" = so_anno$gene_name)
preGO1 <- preGO[!duplicated(preGO),] # Get rid of duplicates

for(i in 1:ncol(geneModuleMembership)){
  tempdat = data.frame("target_id" = geneModuleMembership$target_id,
                       "Connectivity" = geneModuleMembership[,i])
  tempdat2 = arrange(tempdat, Connectivity)
  tempdat3 = left_join(tempdat2,preGO1, by = "target_id")
  tempdat3$module = rep(colnames(geneModuleMembership[i]),nrow(tempdat3))
  #write.csv(tempdat3, paste0("ModuleMembership_",colnames(geneModuleMembership[i]),".csv"))
}

## Assign genes to module membership
datME = moduleEigengenes(rl9,mergedColors)$eigengenes
datKME = signedKME(rl9, datME, outputColumnName="MM.") # Calculates module membership
genes = names(rl9)
geneInfo0 = data.frame(gene=genes,moduleColor=moduleColors) # Gene assignments
color = data.frame(geneInfo0, datKME, MMPvalue) # From original WGCNA analysis 

## Define the modules of interest
whichModule = unique(sig_mod)
modcol = paste("kME",whichModule,sep="")

## Establish site order for plotting tadpole salinity
xnames1 = c("0_0_BL_3","0_0_LO_3","0_0_PWL_3","0_0_WHF_1","4_4_BL_6","4_4_LO_2","4_4_PWL_6","4_4_WHF_2",
            "0_6_BL_4","0_6_LO_2","0_6_PWL_3","0_6_WHF_2","4_6_BL_2","4_6_LO_2","4_6_PWL_2","4_6_WHF_3",
            "0_0_BOD_3","0_0_DQ_5", "0_0_LH_5","0_6_BOD_3","4_4_BOD_5","4_4_DQ_6","4_4_LH_6","4_4_CSI_2",
            "0_6_DQ_2","0_6_LH_2","0_6_CSI_2","4_6_BOD_3","4_6_DQ_2","4_6_LH_2","4_6_CSI_1")
xnames1 = as.factor(xnames1)
xnames2 = as.numeric(xnames1)

## Re-run datasets 
vsd = rl9 
vsd %>% drop_na()
a.vsd = data.frame(t(vsd)) 
allkME = as.data.frame(signedKME(t(a.vsd), MEs)) # Calculates eigengene-based connectivity, aka module membership

## For-loop to extract genes in each significant module
for(i in 1:length(modcol)){
  tempdat = order(allkME[,modcol[i]],decreasing=T)
  tempdat2 = a.vsd[tempdat,] # Arrange genes in module according to connectivity
  setDT(tempdat2, keep.rownames = TRUE)[]
  colnames(tempdat2)[1] = "target_id"
  #write.csv(tempdat2, paste("ModuleMembership",modcol[i],".csv",sep="_"), na="", row.names=FALSE)
}

## Attach to dataframe for GOrilla
GoDat <- left_join(Ranked_ModMemb,preGO1, by = "target_id")

## Make CSV for GOrilla (use this at http://cbl-gorilla.cs.technion.ac.il/)
write.csv(GoDat, paste("RankedModMembership_GOrilla.csv",sep="_"), na="", row.names=FALSE)

