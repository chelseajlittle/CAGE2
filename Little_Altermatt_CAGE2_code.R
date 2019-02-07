## This code is for the analysis of the results presented in:
## Little C & Altermatt F. Differential Resource Consumption in Leaf Litter Mixtures by Native and Non-Native Amphipods
## Submitted to Aquatic Ecology
## Version 20181026
## contact with any questions: chelsea.jean.little@gmail.com
## Eawag: Swiss Federal Institute of Aquatic Science and Technology,
## Department of Aquatic Ecology, Überlandstrasse 133, CH-8600 Dübendorf, Switzerland

#### required packages for this analysis #####

# load packages
library(dplyr)
library(ggplot2)
library(FSA)
library(gridExtra)
library(multcomp)
library(multcompView)
library(tidyr)

##### load the data #####

## set a path "path.data" to wherever you downloaded the data to

## read in data 

CAGE2data <- read.csv(paste(path.data, "/Experimental_data_CAGE2.csv",
                            sep=""), 
                      sep=";", dec=".")


#### data processing for the leaf mass loss part ####

## to do this we use allometric relationships we previously developed
## source: Little C.J. & Altermatt F. (2018) 
## Species turnover and invasion of dominant freshwater invertebrates alter biodiversity-ecosystem-function relationship. 
## Ecological Monographs 88, 461-480.

# convert oak area loss to mass loss
CAGE2data$Mass_loss_Oak_mg <- (CAGE2data$Area_loss_Oak_mm2*0.0762) - 16.34

# convert beech area loss to mass loss
CAGE2data$Mass_loss_Beech_mg <- (CAGE2data$Area_loss_Beech_mm2*0.0194) + 54.84

# convert alder area loss to mass loss
CAGE2data$Mass_loss_Alder_mg <- (CAGE2data$Area_loss_Alder_mm2*0.0552) + 10.89

#### Before doing the analyses, the mass loss must be converted to a few different metrics

## calculate the total mass loss in each enclosure
## (because some enclosures have the three-species mix, must sum the loss of different species)

CAGE2data$Total_mass_loss_mg <- rowSums(CAGE2data[,c("Mass_loss_Oak_mg",
                                                     "Mass_loss_Beech_mg",
                                                     "Mass_loss_Alder_mg")],
                                        na.rm = TRUE)
  

## the mass loss should also be corrected for mass loss in the controls

## first calculate this mass loss in the controls only
## make a dataset for this

CAGE2_control_data <- CAGE2data[which(CAGE2data$Amphi.TRT=="CTR"),-c(3:7)]

## calculate the mean mass loss in controls

CAGE2_control_summary <- CAGE2_control_data %>%
  group_by(Leaf.TRT) %>%
  summarise(mean_loss = mean(Total_mass_loss_mg), 
            sd_loss = sd(Total_mass_loss_mg))

## there is one outlier, which we will exclude

CAGE2_control_data_outexclude <- 
  CAGE2_control_data[which(CAGE2_control_data$Total_mass_loss_mg < 300),]

CAGE2_control_summary_2 <- CAGE2_control_data_outexclude %>%
  group_by(Leaf.TRT) %>%
  summarise(mean_loss = mean(Total_mass_loss_mg), 
            sd_loss = sd(Total_mass_loss_mg))

## store these values as variables

CTRmeans <- CAGE2_control_summary_2 %>% 
  dplyr::select(mean_loss) %>% unlist(use.names = FALSE)

Oak_CTR_loss <- CTRmeans[4]
Beech_CTR_loss <- CTRmeans[3]
Alder_CTR_loss <- CTRmeans[2]
SP3_CTR_loss <- CTRmeans[1]

## make a dataset which is only the non-control data, i.e. the enclosures which had amphipods

CAGE2_exp_data <- CAGE2data[which(CAGE2data$Amphi.TRT!="CTR"),]

## subtract the mass loss in controls from the total mass loss in amphipod enclosures

CAGE2_exp_data$Total_mass_loss_mg_adjusted <- NA

# for Oak
CAGE2_exp_data$Total_mass_loss_mg_adjusted[CAGE2_exp_data$Leaf.TRT=="Oak"] <- 
  CAGE2_exp_data$Total_mass_loss_mg[CAGE2_exp_data$Leaf.TRT=="Oak"] - Oak_CTR_loss

# for Beech
CAGE2_exp_data$Total_mass_loss_mg_adjusted[CAGE2_exp_data$Leaf.TRT=="Beech"] <- 
  CAGE2_exp_data$Total_mass_loss_mg[CAGE2_exp_data$Leaf.TRT=="Beech"] - 
  Beech_CTR_loss

# for Alder
CAGE2_exp_data$Total_mass_loss_mg_adjusted[CAGE2_exp_data$Leaf.TRT=="Alder"] <- 
  CAGE2_exp_data$Total_mass_loss_mg[CAGE2_exp_data$Leaf.TRT=="Alder"] -
  Alder_CTR_loss

# for 3SP
CAGE2_exp_data$Total_mass_loss_mg_adjusted[CAGE2_exp_data$Leaf.TRT=="3SP"] <- 
  CAGE2_exp_data$Total_mass_loss_mg[CAGE2_exp_data$Leaf.TRT=="3SP"] - SP3_CTR_loss

#### now the density of amphipods must be processed/calculated ####

## to estimate average density, we just averaged the initial density (12) and final density
## for a few enclosures, there were >12 amphipods at the end...
## whether by experimenter error or somehow they came in during experiment, we don't know.
## in that case, we used final density as the average (because don't know how it changed during experiment)

CAGE2_exp_data$avg_density <- ifelse(CAGE2_exp_data$Amphis.Final > 12, 
                                CAGE2_exp_data$Amphis.Final, 
                                ((CAGE2_exp_data$Amphis.Final+12)/2))

## calculate the biomass of an enclosurse by multiplying this average density
## by the average weight of the amphipods weighed from that enclosure

CAGE2_exp_data$amphi_biomass <-
  CAGE2_exp_data$avg_density*CAGE2_exp_data$Weight.per.amphi..mg.

#### calculate density-specific consumption rates ####

## for any enclosure where adjusting the mass loss by the mass loss in controls made it negative
## we re-set it to simply be zero mass loss

CAGE2_exp_data$Total_mass_loss_mg_zeroadjusted <- 
  ifelse(CAGE2_exp_data$Total_mass_loss_mg_adjusted < 0, 
         0,
         CAGE2_exp_data$Total_mass_loss_mg_adjusted)

## calculate daily mass loss by dividing by length of experiment (27 days)
## do this per capita (just divide by density) and per biomass (divide by biomass)

CAGE2_exp_data$mg_loss_per_amphi <- 
  (CAGE2_exp_data$Total_mass_loss_mg_zeroadjusted/CAGE2_exp_data$avg_density)/27

CAGE2_exp_data$mg_loss_per_biomass <-
  (CAGE2_exp_data$Total_mass_loss_mg_zeroadjusted/CAGE2_exp_data$amphi_biomass)/27

#### data wrangling to permit plotting ####

## for any plotting, etc., the treatment codes will have to be combined
CAGE2_exp_data$TRT <- paste(CAGE2_exp_data$Amphi.TRT, 
                            CAGE2_exp_data$Leaf.TRT, sep="-")

CAGE2_exp_data$TRT <- factor(CAGE2_exp_data$TRT, 
                             levels = c("GF-Alder", "GR-Alder", 
                                        "GF-Oak", "GR-Oak",
                                        "GF-Beech", "GR-Beech",
                                        "GF-3SP", "GR-3SP"), ordered=TRUE)

#### examine biomass-specific total consumption per enclosure ####

## summary statistics

totconsum_perbio_summary <- CAGE2_exp_data %>%
  group_by(TRT, Amphi.TRT) %>%
  summarise(mean_consum = mean(mg_loss_per_biomass),
            sd_consum = sd(mg_loss_per_biomass),
            se_consum =
              sd(mg_loss_per_biomass)/sqrt(length(mg_loss_per_biomass)))

## plot of total consumption rates (Figure 2a)

amphicols <- c("GF" = "#d8daeb", "GR" = "#f1a340")

(totconsum_plot <- ggplot(data=totconsum_perbio_summary, 
                          aes(x=TRT, y=mean_consum, fill=Amphi.TRT)) +
  xlab("Treatment") + ylab ("Consumption (per biomass)")+
  geom_bar(position=position_dodge(), stat="identity", 
           color="black")+
   scale_fill_manual("Species", values=amphicols)+
  geom_point(data=CAGE2_exp_data, 
             mapping=aes(x=TRT, y=mg_loss_per_biomass), shape=21,
             fill="darkgray", alpha=0.5, size=4)+
  geom_errorbar(aes(ymin=mean_consum-se_consum, ymax=mean_consum+se_consum),
                width=0.2)+
  theme_classic()+
  theme(axis.text.x=element_text(size=10), 
        axis.text.y=element_text(size=10),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=12),
        legend.position=c(0.8,0.8)))

## test a linear model for examining this

totconsum_lm_fact <- lm(mg_loss_per_biomass ~ Amphi.TRT * Leaf.TRT,
                        CAGE2_exp_data)

shapiro.test(resid(totconsum_lm_fact)) 

## the shapiro test shows that the residuals are highly non-normal, 
## this is not included, but we tested some data transformations, and they did not help much with residuals

## therefore, we used non-parametric tests

## a kruskal-wallis factorial model

kruskal.test(mg_loss_per_biomass ~ TRT, data = CAGE2_exp_data)

## perform a Dunn test on this full model 

pairwise_perbio_test <- dunnTest(mg_loss_per_biomass ~ TRT, 
                                 data = CAGE2_exp_data,
                                 method="holm") 
pairwise_perbio_test_table <- pairwise_perbio_test$res

## in fact, we want to just compare the two amphipod species for each leaf type

comparisons <- c("GF-Alder - GR-Alder", 
                 "GF-Beech - GR-Beech",
                 "GF-Oak - GR-Oak",
                 "GF-3SP - GR-3SP")

pairwise_perbio_test_table_subset <- 
  pairwise_perbio_test_table[which(pairwise_perbio_test_table$Comparison 
                                   %in% comparisons),]

## apply the holm-bonferroni correction to this

pairwise_perbio_test_table_subset$P.adj.2 <- 
  p.adjust(pairwise_perbio_test_table_subset$P.unadj, "holm")

## this indicates that amphipod species are not different, so redo the kruskal-wallis test on just leaf type

kruskal.test(mg_loss_per_biomass ~ Leaf.TRT, 
             data = CAGE2_exp_data) # very significant

## apply a Dunn test

dunnTest(mg_loss_per_biomass ~ Leaf.TRT, 
                                 data = CAGE2_exp_data,
                                 method="holm") 

dunntest_perbio_leaf <- dunnTest(mg_loss_per_biomass ~ Leaf.TRT, 
                                 data = CAGE2_exp_data,
                                 method="holm")$res

##  generate letters for display output 

perbio_leaf_pvals_adj <- matrix(NA,4,4)
colnames(perbio_leaf_pvals_adj) <- levels(CAGE2_exp_data$Leaf.TRT)
rownames(perbio_leaf_pvals_adj) <- levels(CAGE2_exp_data$Leaf.TRT)

## we will do this manually

diag(perbio_leaf_pvals_adj) <- 1
perbio_leaf_pvals_adj[2,1] <- dunntest_perbio_leaf[1,4]
perbio_leaf_pvals_adj[3,1] <- dunntest_perbio_leaf[2,4]
perbio_leaf_pvals_adj[4,1] <- dunntest_perbio_leaf[4,4]
perbio_leaf_pvals_adj[3,2] <- dunntest_perbio_leaf[3,4]
perbio_leaf_pvals_adj[4,2] <- dunntest_perbio_leaf[5,4]
perbio_leaf_pvals_adj[4,3] <- dunntest_perbio_leaf[6,4]

multcompLetters(perbio_leaf_pvals_adj)

#### examine the per-capita total consumption rates ####

## summary statistics

totconsum_peramphi_summary <- CAGE2_exp_data %>%
  group_by(TRT, Amphi.TRT) %>%
  summarise(mean_consum = mean(mg_loss_per_amphi),
            sd_consum = sd(mg_loss_per_amphi),
            se_consum =
              sd(mg_loss_per_amphi)/sqrt(length(mg_loss_per_amphi)))

## plot it (figure 2b)

(totconsum_plot_peramphi <- ggplot(data=totconsum_peramphi_summary, 
                          aes(x=TRT, y=mean_consum, fill=Amphi.TRT)) +
  xlab("Treatment") + ylab ("Consumption (per capita)")+
  geom_bar(position=position_dodge(), stat="identity", 
           color="black")+
   scale_fill_manual("Species", values=amphicols)+
  geom_point(data=CAGE2_exp_data, 
             mapping=aes(x=TRT, y=mg_loss_per_amphi), shape=21,
             fill="darkgray", alpha=0.5, size=4)+
  geom_errorbar(aes(ymin=mean_consum-se_consum, ymax=mean_consum+se_consum),
                width=0.2)+
  theme_classic()+
  theme(axis.text.x=element_text(size=10), axis.text.y=element_text(size=10),
        axis.title.x=element_text(size=12), axis.title.y=element_text(size=12),
        legend.position=c(0.8,0.8)))

## apply the kruskal-wallis test on the factorial model 

kruskal.test(mg_loss_per_amphi ~ TRT, data = CAGE2_exp_data)

## run a posthoc Dunn test

pairwise_peramphi_test <- dunnTest(mg_loss_per_amphi ~ TRT, 
                                 data = CAGE2_exp_data,
                                 method="holm") 
pairwise_peramphi_test_table <- pairwise_peramphi_test$res

## extract only the relevant comparisons from table

pairwise_peramphi_test_table_subset <- 
  pairwise_peramphi_test_table[which(pairwise_peramphi_test_table$Comparison 
                                   %in% comparisons),]

## apply the holm-bonferroni correction to this
pairwise_peramphi_test_table_subset$P.adj.2 <- 
  p.adjust(pairwise_peramphi_test_table_subset$P.unadj, "holm")

## generate letters for the big comparison

peramphi_pvals_adj <- matrix(NA,8,8)
colnames(peramphi_pvals_adj) <- levels(CAGE2_exp_data$TRT)
rownames(peramphi_pvals_adj) <- levels(CAGE2_exp_data$TRT)

## we will do this manually

diag(peramphi_pvals_adj) <- 1
peramphi_pvals_adj[2,1] <- pairwise_peramphi_test_table[12,4] 
peramphi_pvals_adj[3,1] <- pairwise_peramphi_test_table[5,4]
peramphi_pvals_adj[4,1] <- pairwise_peramphi_test_table[23,4]
peramphi_pvals_adj[5,1] <- pairwise_peramphi_test_table[3,4]
peramphi_pvals_adj[6,1] <- pairwise_peramphi_test_table[17,4]
peramphi_pvals_adj[7,1] <- pairwise_peramphi_test_table[1,4]
peramphi_pvals_adj[8,1] <- pairwise_peramphi_test_table[8,4]
peramphi_pvals_adj[3,2] <- pairwise_peramphi_test_table[14,4] 
peramphi_pvals_adj[4,2] <- pairwise_peramphi_test_table[27,4] 
peramphi_pvals_adj[5,2] <- pairwise_peramphi_test_table[13,4] 
peramphi_pvals_adj[6,2] <- pairwise_peramphi_test_table[21,4] 
peramphi_pvals_adj[7,2] <- pairwise_peramphi_test_table[11,4] 
peramphi_pvals_adj[8,2] <- pairwise_peramphi_test_table[15,4] 
peramphi_pvals_adj[4,3] <- pairwise_peramphi_test_table[25,4] 
peramphi_pvals_adj[5,3] <- pairwise_peramphi_test_table[6,4] 
peramphi_pvals_adj[6,3] <- pairwise_peramphi_test_table[19,4] 
peramphi_pvals_adj[7,3] <- pairwise_peramphi_test_table[4,4] 
peramphi_pvals_adj[8,3] <- pairwise_peramphi_test_table[10,4] 
peramphi_pvals_adj[5,4] <- pairwise_peramphi_test_table[24,4] 
peramphi_pvals_adj[6,4] <- pairwise_peramphi_test_table[28,4] 
peramphi_pvals_adj[7,4] <- pairwise_peramphi_test_table[22,4] 
peramphi_pvals_adj[8,4] <- pairwise_peramphi_test_table[26,4] 
peramphi_pvals_adj[6,5] <- pairwise_peramphi_test_table[18,4] 
peramphi_pvals_adj[7,5] <- pairwise_peramphi_test_table[2,4] 
peramphi_pvals_adj[8,5] <- pairwise_peramphi_test_table[9,4]
peramphi_pvals_adj[7,6] <- pairwise_peramphi_test_table[16,4]
peramphi_pvals_adj[8,6] <- pairwise_peramphi_test_table[20,4] 
peramphi_pvals_adj[8,7] <- pairwise_peramphi_test_table[7,4] 

multcompLetters(peramphi_pvals_adj)

## note: significance letters were added to plots manually in Adobe Illustrator

#### Examine percent loss in monocultures vs. mixtures ####

## calculate percent loss of leaves of each type in each enclosure

CAGE2data$percent_loss_alder <- 
  (CAGE2data$Area_Alder_before - CAGE2data$Area_Alder_after)/
  CAGE2data$Area_Alder_before

CAGE2data$percent_loss_beech <- 
  (CAGE2data$Area_Beech_before - CAGE2data$Area_Beech_after)/
  CAGE2data$Area_Beech_before

CAGE2data$percent_loss_oak <- 
  (CAGE2data$Area_Oak_before - CAGE2data$Area_Oak_after)/
  CAGE2data$Area_Oak_before

## replace any negative values with zero

CAGE2data$percent_loss_alder <- ifelse(CAGE2data$percent_loss_alder < 0, 0,
                                       CAGE2data$percent_loss_alder)

CAGE2data$percent_loss_beech <- ifelse(CAGE2data$percent_loss_beech < 0, 0,
                                       CAGE2data$percent_loss_beech)

CAGE2data$percent_loss_oak <- ifelse(CAGE2data$percent_loss_oak < 0, 0,
                                       CAGE2data$percent_loss_oak)

## extract relevant columns and convert the data from wide to long

CAGE2_forlong_2 <- CAGE2data[,c(1:3,26:29)]

CAGE2_long_2 <- gather(CAGE2_forlong_2, species, measurement,
                     percent_loss_alder:percent_loss_oak, factor_key=TRUE)

## remove the rows with no measurement

CAGE2_long_2 <- subset(CAGE2_long_2, (!is.na(CAGE2_long_2[,6])))

## summary statistics of percent mass loss

percentloss_diversity_summary <- CAGE2_long_2 %>%
  group_by(diversity_level,Amphi.TRT,species) %>%
  summarise(mean_loss = mean(measurement),
            sd_loss = sd(measurement),
            se_loss = sd(measurement)/sqrt(length(measurement)))

## set the order for factor levels to be displayed in plots

percentloss_diversity_summary$diversity_level <-
  factor(percentloss_diversity_summary$diversity_level,
         levels=c("mono", "mix"), ordered=TRUE)

percentloss_diversity_summary$species <-
  factor(percentloss_diversity_summary$species,
         levels=c("percent_loss_alder", "percent_loss_beech",
                  "percent_loss_oak"), 
         ordered=TRUE)

levels(percentloss_diversity_summary$species)[levels(percentloss_diversity_summary$species)=="percent_loss_alder"] <- "alder"
levels(percentloss_diversity_summary$species)[levels(percentloss_diversity_summary$species)=="percent_loss_beech"] <- "beech"
levels(percentloss_diversity_summary$species)[levels(percentloss_diversity_summary$species)=="percent_loss_oak"] <- "oak"

CAGE2_long_2$diversity_level <-
  factor(CAGE2_long_2$diversity_level,
         levels=c("mono", "mix"), ordered=TRUE)

CAGE2_long_2$species <-
  factor(CAGE2_long_2$species,
         levels=c("percent_loss_alder", "percent_loss_beech",
                  "percent_loss_oak"), 
         ordered=TRUE)

levels(CAGE2_long_2$species)[levels(CAGE2_long_2$species)=="percent_loss_alder"] <- "alder"
levels(CAGE2_long_2$species)[levels(CAGE2_long_2$species)=="percent_loss_beech"] <- "beech"
levels(CAGE2_long_2$species)[levels(CAGE2_long_2$species)=="percent_loss_oak"] <- "oak"

## plot this data (Figure 3)

(percent_loss_plot <- ggplot(aes(y=mean_loss, x=Amphi.TRT, fill=diversity_level),
                   data=percentloss_diversity_summary) +
  xlab("Treatment") + ylab ("Percent Area Loss")+
  geom_bar(position=position_dodge(), stat="identity", 
           color="black")+
  scale_fill_manual("Diversity", values=monomixcols)+
  geom_point(data=CAGE2_long_2, 
             mapping=aes(y=measurement, x=Amphi.TRT, 
                         group=diversity_level),
             position=position_dodge(width=0.9),
             shape=21, fill="darkgray", alpha=0.5, size=4)+
  geom_errorbar(aes(ymin=mean_loss-se_loss, ymax=mean_loss+se_loss),
               position=position_dodge(width=0.9), width=0.2)+
  facet_wrap(vars(species), nrow = 1)+
  theme_classic()+
  theme(axis.text.x=element_text(size=10), 
        axis.text.y=element_text(size=10),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=12),
        legend.position=c(0.8,0.8)))

## model percent loss as a linear function of leaf type, amphipod species, and diversity level

percentlosslm <- lm(measurement ~ species + Amphi.TRT + diversity_level, 
                data=CAGE2_long_2)

shapiro.test(resid(percentlosslm)) 

## the residuals are not normal, however, if we square-root transform, they become normal

percentlosslm2 <- lm(sqrt(measurement) ~ 
                       species + Amphi.TRT + diversity_level, 
                data=CAGE2_long_2)

shapiro.test(resid(percentlosslm2)) 

## so, we do a model selection round to see whether a factorial or additive model is best
## and if all factors are included

percentloss_fullfact_lm <- lm(sqrt(measurement) ~ 
                       species * Amphi.TRT * diversity_level, 
                data=CAGE2_long_2)

step_percent_loss <- step(percentloss_fullfact_lm)

summary(step_percent_loss)

## because the interaction is significant, in order to do posthoc testing,
## we must create treatment codes that show all three factors combined

CAGE2_long_2$TRTcombo <- paste(CAGE2_long_2$species,
                               CAGE2_long_2$diversity_level, 
                               CAGE2_long_2$Amphi.TRT, 
                               sep="_")

## order this factor for display in plots later

CAGE2_long_2$TRTcombo <- factor(CAGE2_long_2$TRTcombo, 
                                levels=c("alder_mono_CTR", "alder_mix_CTR",
                                         "alder_mono_GF", "alder_mix_GF",
                                         "alder_mono_GR", "alder_mix_GR",
                                         "beech_mono_CTR", "beech_mix_CTR",
                                         "beech_mono_GF", "beech_mix_GF",
                                         "beech_mono_GR", "beech_mix_GR",
                                         "oak_mono_CTR", "oak_mix_CTR",
                                         "oak_mono_GF", "oak_mix_GF",
                                         "oak_mono_GR", "oak_mix_GR"),
                                ordered=TRUE)

## for posthoc testing, redo the linear model using this factor rather than the three-factor interaction

monomix_combinedtreatment_lm <- 
  lm(sqrt(measurement) ~ TRTcombo, 
     data = CAGE2_long_2)

## define specific contrasts, the relevant ones
## we only compared treatment combination with two factor levels in common,
## and did not make all possible comparisons (in order to reduce number of tests)

tukeytest_chosencombos <- 
  glht(monomix_combinedtreatment_lm, 
       linfct = mcp(TRTcombo = c("alder_mix_GF - alder_mix_CTR == 0 ",
                                 "alder_mix_GR - alder_mix_CTR == 0 ",
                                 "alder_mono_CTR - alder_mix_CTR == 0",
                                 "beech_mix_CTR - alder_mix_CTR == 0",
                                 "oak_mix_CTR - alder_mix_CTR == 0",
                                 "alder_mono_GF - alder_mono_CTR == 0",
                                 "alder_mono_GR - alder_mono_CTR == 0",
                                 "beech_mono_CTR - alder_mono_CTR == 0",
                                 "oak_mono_CTR - alder_mono_CTR == 0",
                                 "alder_mono_GF - alder_mix_GF == 0",
                                 "alder_mix_GR - alder_mix_GF == 0",
                                 "beech_mix_GF - alder_mix_GF == 0",
                                 "oak_mix_GF - alder_mix_GF == 0",
                                 "alder_mono_GR - alder_mono_GF == 0",
                                 "beech_mono_GF - alder_mono_GF == 0",
                                 "oak_mono_GF - alder_mono_GF == 0",
                                 "alder_mono_GR - alder_mix_GR == 0",
                                 "beech_mix_GR - alder_mix_GR == 0",
                                 "oak_mix_GR - alder_mix_GR == 0",
                                 "beech_mono_GR - alder_mono_GR == 0",
                                 "oak_mono_GR - alder_mono_GR == 0",
                                 "beech_mono_CTR - beech_mix_CTR == 0",
                                 "beech_mix_GF - beech_mix_CTR == 0",
                                 "beech_mix_GR - beech_mix_CTR == 0",
                                 "oak_mix_CTR - beech_mix_CTR == 0",
                                 "beech_mono_GF - beech_mono_CTR == 0",
                                 "beech_mono_GR - beech_mono_CTR == 0",
                                 "oak_mono_CTR - beech_mono_CTR == 0",
                                 "beech_mono_GF - beech_mix_GF == 0",
                                 "beech_mix_GR - beech_mix_GF == 0",
                                 "oak_mix_GF - beech_mix_GF == 0",
                                 "beech_mono_GR - beech_mono_GF == 0",
                                 "oak_mono_GF - beech_mono_GF == 0",
                                 "beech_mono_GR - beech_mix_GR == 0",
                                 "oak_mix_GR - beech_mix_GR == 0",
                                 "oak_mono_GR - beech_mono_GR == 0",
                                 "oak_mono_CTR - oak_mix_CTR == 0",
                                 "oak_mix_GF - oak_mix_CTR == 0",
                                 "oak_mix_GR - oak_mix_CTR == 0",
                                 "oak_mono_GF - oak_mono_CTR == 0",
                                 "oak_mono_GR - oak_mono_CTR == 0",
                                 "oak_mono_GF - oak_mix_GF == 0",
                                 "oak_mix_GR - oak_mix_GF == 0",
                                 "oak_mono_GR - oak_mono_GF == 0",
                                 "oak_mono_GR - oak_mix_GR == 0")))

summary(tukeytest_chosencombos)          

#### visualize the relative leaf loss of different species in the mixtures (Figure 4) ####

## subset only the data from the mixed mesocosms

CAGE2_mixture_data <- CAGE2_long[which(CAGE2_long$Leaf.TRT == "3SP"),]

## make the plot

leafcols <- c("Mass_loss_Alder_mg" = "#c2e699", 
              "Mass_loss_Oak_mg" = "#78c679", 
              "Mass_loss_Beech_mg" = "#31a354")

(stackedbar <- ggplot(CAGE2_mixture_data, 
                     aes(x=ID, y=measurement, fill=species))+
  ylab("Mass lost (mg)")+ ylim(c(0,600))+
  geom_bar(stat="identity", 
           color="black")+
  scale_fill_manual("Leaf Species", values = leafcols, 
                    labels = c("Alder", "Beech", "Oak"))+
    geom_segment(x=58, xend=65, y=333, yend=333)+
    geom_segment(x=66, xend=73, y=580, yend=580)+
    geom_segment(x=74, xend=76, y=290, yend=290)+
  theme_classic()+
  theme(axis.text.x=element_blank(), 
        axis.text.y=element_text(size=10),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=12),
        axis.ticks.x=element_blank(), 
        legend.position=c(0.1,0.8)))


#### some summary information about the data ####  

## size of amphipods of the two species

weightsumtable <- CAGE2_exp_data %>%
  group_by(Amphi.TRT) %>%
  summarise(total_amphis_weighed = sum(Amphis.Weighed),
            total_weight = sum(Weight.total..mg.))

weightsumtable <- as.data.frame(weightsumtable)

weightsumtable$avg_weight <- weightsumtable$total_weight /
  weightsumtable$total_amphis_weighed

## summarize the survival rates of the two species

CAGE2_exp_data$survival <- ifelse(CAGE2_exp_data$Amphis.Final > 12, 1,
                                  CAGE2_exp_data$Amphis.Final/12)

CAGE2_exp_data %>% group_by(Amphi.TRT) %>%
  summarise(mean_survival = mean(survival))