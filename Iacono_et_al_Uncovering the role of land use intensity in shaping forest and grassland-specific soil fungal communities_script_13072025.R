#' ---
#' title: "Analysis of IMS data from the Biodiversity Exploratories collected in 300 experimental plots (EPS) in forest (150) and grassland (150)"
#' output:
#'   pdf_document:
#'     keep_tex: true
#' ---


#=== Loading packages ====

#+ 
library(readxl) # import data from excel file
library(magrittr) # piping 
library(tidyverse) # opinionated collection of R packages designed for data science to keep data tidy
library(scales) # used for scaling
library(vegan) # contains the functions for all the species diversity measures
library(kableExtra) # generates tables 
library(ggplot2) # a grammar of graphics package for plotting
library(ggrepel) #necessary to avoid label overlap in figures
library(ggpubr) # based on ggplot to create publication ready figures
library(dismo) # package for the creation of bioclimatic variables
library(vegetarian) # functions for the calculation of biodiversity indices
library(vegan) # contains function for the analysis of ecological data
library(lmerTest) # contains functions for testing the fit of  mixed linear models and define variance percentages
library(lme4) # contains functions for testing the fit of  mixed linear models
library(emmeans) # package to perform estimated marginaL means (EMM) as aposthoc for mixed linear models
library(multcomp) # contains functions to extract compact letter displays from EMM analysis
library(multcompView) # contains functions to extract compact letter displays from EMM analysis
library(phyloseq) # contains functions for the analysis of microbiome data 
library(DESeq2) # contains functions for differential analysis
library(mia) # contains functions for microbial biodiversity analysis
library(miaViz) #visualization for Mia package
library(corrplot) # necessary for the creation of correlation plots
library(indicspecies) #used to implement the calculation for the indicator species
library(lavaan)
library(stringr) #to work with strings
library(EcolUtils) #for permutational rarefication
library(psych) #for correlation with bonferroni correction

#=== Color code ======
pal_lui <- c("steelblue", "grey70", "coral2")
names(pal_lui) <- levels(site_data$lui)


#=== Loading the data ====
##=== Community data ====
#' Importing dataframe containing the ASV (here named OTUs) abundances
spec_abund<- read_excel("./data/ITS2021.xlsx", sheet = "Matrix_ITS21")

#' Importing the taxonomy table
tax_mat<- read_excel("./data/ITS2021.xlsx", sheet = "Taxonomy_ITS21")%>%
  column_to_rownames("OTU")

##===Site data ====
###=== Basic information ====
sites <- read_excel("./data/EPS_details.xlsx")%>%
  column_to_rownames("EPS") #naming this coulmn consistently between dataset allows further merging

###=== Land Use Intensity information ====

#'Import Land Use Index information for the grassland plots
LUI_global <- read_excel("data/LUI_global_2017_2020.xlsx")
LUI_global<-LUI_global%>% 
  mutate(EPS=str_sub(plotid, 1, 5))

#' Imported the dataset containing information about the land use index for forests SMI
SMI <- read_excel("data/31217_9_data_SMI.xlsx",
                  col_types = c("skip", "skip", "text",
                                "text", "text", "skip", "numeric",
                                "numeric", "numeric"))

SMI<-SMI%>%
  filter(year == 2017 | year ==2018 | year == 2019 | year ==2020)%>%
  group_by(plotid)%>%
  summarize(across(SMId:SMI, mean, na.rm = T))%>% 
  mutate(EPS=str_sub(plotid, 1, 5))

SMI[is.na(SMI)] <- 0

##=== Environment data ====

#' Soil environment
#' importing the environmental data table. Columns with light information are skipped
env <- read_csv("data/soil_env.csv",
                   col_types = cols(LWDR_300 = col_skip(),
                                    LWUR_300 = col_skip(), 
                                    PAR_200 = col_skip(),
                                    SWDR_300 = col_skip(), 
                                    SWUR_300 = col_skip()))

#' After the death of the stands in plot HEW02, it was replaced by EP HEW51. 
#' For simplicity HEW02 row is removed from the dataset and row HEW51 is renamed HEW02
env<-env%>%
  filter(!(plotID == "HEW02"))%>%
  mutate(plotID = sub("HEW51", "HEW02", plotID))

#' Infinite and NAs values originating from the previous summary are replaced by 0s
env[sapply(env, is.infinite)] <- 0
env[sapply(env, is.na)] <- 0
env[sapply(env, is.nan)] <- 0

env <- env %>% merge(sites[c(4:10)], .,  by.x = 0, by.y = "plotID", all.y = TRUE) 

env <- env %>% 
  rename(Row.names = "EPS", `Landuse type` = "Type"  )

##=== Soil composition data ====
#importing the dataset with soil chemistry and soil texture information
#requires
#library(readxl)

chem <- read_excel("data/soil_abiotic_data.xlsx",
                   col_types = c("text", "text", "text",
                                 "text", "numeric", "numeric", "numeric",
                                 "numeric", "numeric", "numeric",
                                 "numeric", "numeric", "numeric",
                                 "numeric", "numeric", "numeric",
                                 "numeric", "numeric", "numeric",
                                 "numeric"))

chem<-chem%>%
  rename(plotid = "EPS") %>% 
  mutate(EPS=str_sub(EPS, 1, 5))

#=== Cleaning and preparing data ====
##=== Rarefaction of community data ====
#' Turning the dataframe to have the sites as rows and the names of the sites as rownames

spec_abund<- spec_abund%>% 
  as.data.frame()%>%
  column_to_rownames("OTU")%>%
  t()%>%
  as.data.frame()%>%
  rownames_to_column(var = "EPS")%>%
  mutate(EPS = str_sub(EPS, 1, 5))%>%
  column_to_rownames("EPS")

#' Calculating the minimum number of seq  uences for the rearefaction
min(rowSums(spec_abund))

spec_abund_rar <- spec_abund%>%
  rrarefy.perm(., sample = 42138, n = 999)

#' Verifying rarefaction happened correctly
mean(rowSums(spec_abund_rar))

##=== Normalization using Deseq =====


spec_ds<-DESeqDataSetFromMatrix(countData = t(spec_abund),
                                colData = site_data[c(1:4, 6:58)],
                                design=~1
)

spec_ds$group<-factor(paste(spec_ds$Type, spec_ds$lui))

design(spec_ds)<-~group

spec_ds<-DESeq(spec_ds, test="Wald", fitType="parametric", sfType = "poscounts")

spec_abund_norm<-counts(spec_ds, normalized = TRUE)%>%
  t()

spec_ds_f<-DESeqDataSetFromMatrix(countData = t(spec_abund[c(51:100, 151:200, 251:300), ]),
                                  colData = site_data %>% filter(Type == "Forest"),
                                  design= ~ lui)

spec_ds_f<-DESeq(spec_ds_f, test="Wald", fitType="parametric", sfType = "poscounts")

spec_ds_g<-DESeqDataSetFromMatrix(countData = t(spec_abund[c(1:50, 101:150, 201:250), ]),
                                  colData = site_data %>% filter(Type == "Grassland"),
                                  design= ~ lui)

spec_ds_g<-DESeq(spec_ds_g, test="Wald", fitType="parametric", sfType = "poscounts")

EPS_norm <- row.names(spec_abund_norm)
n_spec_norm <- specnumber(spec_abund_norm)
shannon_norm <- vegan::diversity(spec_abund_norm)
simpson_norm <- vegan::diversity(spec_abund_norm, index = "simpson")
invsimpson_norm <- vegan::diversity(spec_abund_norm, index = "inv")
N1_norm <- exp(shannon_norm)
N2_norm <- exp(invsimpson_norm)
E1_norm <- N1_norm/n_spec_norm 
E2_norm <- N2_norm/n_spec_norm

spec_div_norm<-data.frame(EPS_norm, n_spec_norm, shannon_norm, simpson_norm, invsimpson_norm, N1_norm, N2_norm,E1_norm, E2_norm)%>%
  as_tibble(rownames = "EPS")

spec_div_rar<-merge(spec_div_rar, site_data[c(1, 5, 6, 14)], by = "EPS")

##=== Calculation of diversity indices ====
#biodiversity indices
EPS <- row.names(spec_abund_rar)
n_spec <- specnumber(spec_abund_rar)
shannon <- vegan::diversity(spec_abund_rar)
simpson <- vegan::diversity(spec_abund_rar, index = "simpson")
N1 <- exp(shannon)
N2 <- vegan::diversity(spec_abund_rar, index = "invsimpson")
E1 <- N1/n_spec 
E2 <- N2/n_spec
beta_w <- betadiver(spec_abund_rar, method = NA, order = FALSE)

spec_div_rar<-data.frame(n_spec,
                         shannon, 
                         simpson,
                         N1,
                         N2,
                         E1,
                         E2
                         )%>%
  as_tibble(rownames = "EPS")



##=== Categorization of land use intensity =====
LUI_global<-LUI_global%>% #create a categorical variable
  mutate(lui = factor(ntile(LUI, 3), 
                      levels = c("1", "2", "3"), labels = c("Low", "Medium", "High")),
         lui_norm = ((LUI-min(LUI))/(max(LUI)-min(LUI)))*3)

SMI<-SMI%>%
  mutate(lui=factor(ntile(SMI, 3),
                    levels = c("1", "2", "3"), labels = c("Low", "Medium", "High")),
         lui_norm=((SMI-min(SMI))/(max(SMI)-min(SMI)))*3)

lui<-rbind(SMI%>%
             dplyr::select(EPS, lui, lui_norm), 
           LUI_global%>%
             dplyr::select(EPS, lui, lui_norm)) %>% 
  mutate(lui_norm = round(lui_norm, 2))

lui <- lui %>% 
  merge(., sites[c(4:10)], by.x = "EPS", by.y = 0, all.x = TRUE) %>% 
  rename(`Landuse type` = "Type")

##== Creating differential  light integrals and summarising the data to obtain monthly values =====
#' This step is necessary to create bioclimatic variables as the function dismo::biovars() works on a 12 monthly values.

env<-env%>%
  filter(datetime>"2017-01-01" & datetime<"2021-12-31")%>%
  mutate(hour = hour(datetime), day = day(datetime), month = month(datetime), year = year(datetime))%>%
  mutate(par_dli = PAR*SD_Olivieri*0.0036, 
         lwdr_dli = LWDR * SD_Olivieri,
         lwur_dli = LWUR * SD_Olivieri,
         swdr_dli = SWDR * SD_Olivieri,
         swur_dli = SWUR * SD_Olivieri)%>%
  group_by(plotID, year, month)%>%
  summarise(mois = mean(SM_10, na.rm = T), # we use 10 cm because there is no moisture data for 5 cm
            Ts_10_max = max(Ts_10, na.rm = T),
            Ts_10_min = min(Ts_10, na.rm = T),
            par_dli = mean(par_dli, na.rm = T),
            par_dli = sum(par_dli),
            lwdr_dli = sum(lwdr_dli),
            lwur_dli = sum (lwur_dli),
            swdr_dli = sum (swdr_dli),
            swur_dli = sum (swur_dli)
  )%>%
  ungroup()%>%
  group_by(plotID, month)%>% #the monthly values used for the analysis are the average of the 6 year prior the sampling
  summarise(mois = mean(mois, na.rm = T), # we use 10 cm because there is no moisture data for 5 cm
            Ts_10_max = mean(Ts_10_max, na.rm = T),
            Ts_10_min = mean(Ts_10_min, na.rm = T),
            par_dli = mean(par_dli, na.rm = T),
            par_dli = mean(par_dli),
            lwdr_dli = mean(lwdr_dli),
            lwur_dli = mean(lwur_dli),
            swdr_dli = mean(swdr_dli),
            swur_dli = mean(swur_dli))

##=== Calculation of bioclimatic variables ====
biovars = data.frame()
names = levels(as.factor(env$EPS))
source("./functions/soil_biovars.R") # see additional R file for the function

for (i in 1:length(levels(as.factor(env$EPS)))) {
  plot = names [i]
  print(plot)
  df = env%>%
    filter(EPS == plot)
  biovalues = soil_biovars(df$mois, df$Ts_10_min, df$Ts_10_max, df$par_dli)
  biovalues = cbind(EPS = plot, round(biovalues, 3))
  biovars=rbind(biovars, biovalues)
}

biovars[2:25] <- sapply(biovars[2:25],as.numeric)

#biovars <- biovars%>%
#mutate(EPS=paste(EPS, "2021", sep = "-"))%>% #the data represent an average from 2016 to 2021
#mutate(EPS=str_sub(plotid, 1, 5))

##=== Merging datasets together =====

env_all <- env %>% 
  group_by(EPS, Exploratorium, Type) %>% 
  summarise(across(mois:swur_dli, ~mean(.x, na.rm = TRUE)))%>%
  merge(. , biovars, 
        by= "EPS", 
        )

rm(i, names, plot) #removes objects created by the loop

#' Creating a matrix containing the sites information
#' 
site_data <- merge(lui, env_all, by = c("EPS", "Exploratorium", "Type")) %>% 
  merge(., chem, by = c("EPS")) %>%
  rename (Exploratorium.x = "Exploratorium", Type.x = "Type") %>% 
  dplyr::select(-c(Exploratorium.y, Type.y, Plotid))

site_data<-site_data%>%
  mutate(Type_lui= paste(Type, lui, sep = "_"), .before = "Elevation")

#' Scaling the site data

site_data_scaled <- site_data%>%
  dplyr::select(EPS,
                Exploratorium,
                Type,
                lui,
                lui_norm,
                Elevation:Coarse_Sand)%>%
  mutate(across(Elevation:Coarse_Sand, ~ c(scale(., center = TRUE, scale = TRUE)))) #scaling to account for the different units of the variables

#=== Statistics ====
##=== Exploratory statistics ====
###=== Site data ====
site_summary <- site_data |>
  group_by(Type, lui) |>
  summarise(across(mois:Coarse_Sand, list(mean = ~mean(., na.rm = T), sd = ~sd(., na.rm = T)), .names = 
                     "{.col}.{.fn}")) |>
  mutate_if(is.numeric, round, digits = 2)|>
  pivot_longer(cols = c(mois.mean:Coarse_Sand.sd), values_to = "value", names_to = "variable")%>%
  separate(variable, c("variable", "stat"), sep = "\\.")%>%
  pivot_wider(names_from = c(Type, lui, stat), values_from = value)%>%
  mutate(variable = factor(variable, levels = c("mois","Ts_10_max",
                                                "Ts_10_min",
                                                "par_dli",
                                                "lwdr_dli",
                                                "lwur_dli",
                                                "swdr_dli",
                                                "swur_dli",
                                                "bio1",
                                                "bio2",
                                                "bio3",
                                                "bio4",
                                                "bio5",
                                                "bio6",
                                                "bio7",
                                                "bio8",
                                                "bio9",
                                                "bio10",
                                                "bio11",
                                                "bio12",
                                                "bio13",
                                                "bio14",
                                                "bio15",
                                                "bio16",
                                                "bio17",
                                                "bio18",
                                                "bio19",
                                                "bio20",
                                                "bio21",
                                                "bio22",
                                                "bio23",
                                                "bio24",
                                                "Total_C",
                                                "Inorganic_C",
                                                "Organic_C",
                                                "Total_N",
                                                "CN_ratio",
                                                "Total_S",
                                                "CS_ratio",
                                                "Olsen_P",
                                                "pH",
                                                "Clay",
                                                "Fine_Silt",
                                                "Medium_Silt",
                                                "Coarse_Silt",
                                                "Fine_Sand",
                                                "Medium_Sand",
                                                "Coarse_Sand" )
  )
  )

#' Exporting summary of site data
site_summary%>%
  write.csv("./tables/env_summary.csv", row.names = FALSE)

###=== Diversity indices ====
spec_div_rar<-spec_div_rar%>%
  group_by(EPS)%>%
  summarize(n_spec = mean(n_spec),
            shannon = mean(shannon),
            simpson = mean (simpson),
            N1 = mean(N1),
            N2 = mean(N2),
            E1 = mean(E1),
            E2 = mean(E2))

spec_div_rar <- spec_div_rar%>%
  merge(lui, ., by = "EPS")

#' Summarising the diversity indices dataset
div_summary <- spec_div_rar |>
  group_by(Type, lui) |>
  summarise(across(c(n_spec, N1:E2), list(mean = ~mean(., na.rm = T), sd = ~sd(., na.rm = T)), .names = 
                     "{.col}.{.fn}")) |>
  mutate_if(is.numeric, round, digits = 2)|>
  pivot_longer(cols = c(n_spec.mean:E2.sd), values_to = "value", names_to = "variable")%>%
  separate(variable, c("variable", "stat"), sep = "\\.")%>%
  pivot_wider(names_from = c(Type, lui, stat), values_from = value)%>%
  mutate(variable = factor(variable, levels = c("n_spec",
                                                "N1",
                                                "N2",
                                                "E1",
                                                "E2")
  )
  )

div_summary%>%
  write.csv("./tables/div_summary.csv", row.names = FALSE)

##=== Mixed Linear Model test====
###=== Site variables =====
####=== Effect of ecosystem type and LUI ==== 

#' Exploratory is considered as a fixed effect in the model
lme_res<-data.frame()
lme_post<-data.frame()
plots <-list()
for (i in 1:length(names(site_data)[c(12:59)])){
  names = names(site_data)[c(12:59)]
  model<-lmerTest::lmer(get(names[i])~Type*lui+(1|Exploratorium), data = site_data)
  print(names[i])
  anova_model<-as.data.frame(anova(model))
  anova_model$var<-paste(names[i])
  lme_res<-rbind(lme_res, anova_model)
  posthoc<-emmeans(model, c("Type", "lui"),
                   adjust = "mvt",
                   lmer.df = "satterthwaite"
  )
  outREML<-cld(posthoc,
               level = 0.05,
               Letters = letters,
               decreasing=TRUE)
  outREML.df<-as.data.frame(outREML)
  outREML.df$var<-paste(names[i])
  lme_post<-rbind(lme_post, outREML.df[ ,c(9, 1:8)])
  posthoc.f<-as.data.frame(posthoc)
  p<-ggplot(posthoc.f, aes(x=Type, y=emmean, group=lui, color=lui)) + theme_bw()+
    geom_pointrange(aes(ymin=lower.CL, ymax=upper.CL), size=1.25,
                    position=position_dodge(0.7))+
    labs(x = "Ecosystem type", y = paste(names[i], "(emmean)"), color = "Land use intensity")+
    scale_color_manual(values = c("tomato", "springgreen4","blue3"))
  plots[[i]]<- p
}

#' Exporting the MLM results for site
lme_res%>%
  mutate_if(is.numeric, round, digits = 3)%>%
  rownames_to_column("Factor")%>%
  dplyr::select("var", "Factor",  "Sum Sq",  "Mean Sq", "NumDF",   "DenDF",   "F value", "Pr(>F)")%>%
  write.csv("./tables/mlm_res_env.csv", row.names = FALSE)

####=== Table with estimated marginal means (emm) ====
lme_post%>%
  dplyr::select(!(lower.CL:upper.CL))%>%
  arrange(var, Type, lui)%>%
  pivot_wider(names_from = c(Type, lui), values_from = c(emmean:.group), names_sort = TRUE)%>%
  dplyr::select(var, emmean_Forest_Low, 
                SE_Forest_Low, 
                df_Forest_Low, 
                .group_Forest_Low, 
                emmean_Forest_Medium, 
                SE_Forest_Medium, 
                df_Forest_Medium, 
                .group_Forest_Medium,
                emmean_Forest_High, 
                SE_Forest_High, 
                df_Forest_High, 
                .group_Forest_High,
                emmean_Grassland_Low, 
                SE_Grassland_Low, 
                df_Grassland_Low, 
                .group_Grassland_Low,
                emmean_Grassland_Medium, 
                SE_Grassland_Medium, 
                df_Grassland_Medium, 
                .group_Grassland_Medium,
                emmean_Grassland_High, 
                SE_Grassland_High, 
                df_Grassland_High, 
                .group_Grassland_High)%>%
  mutate(var = factor(var, levels = c("mois","Ts_10_max",
                                                "Ts_10_min",
                                                "par_dli",
                                                "lwdr_dli",
                                                "lwur_dli",
                                                "swdr_dli",
                                                "swur_dli",
                                                "bio1",
                                                "bio2",
                                                "bio3",
                                                "bio4",
                                                "bio5",
                                                "bio6",
                                                "bio7",
                                                "bio8",
                                                "bio9",
                                                "bio10",
                                                "bio11",
                                                "bio12",
                                                "bio13",
                                                "bio14",
                                                "bio15",
                                                "bio16",
                                                "bio17",
                                                "bio18",
                                                "bio19",
                                                "bio20",
                                                "bio21",
                                                "bio22",
                                                "bio23",
                                                "bio24",
                                                "Total_C",
                                                "Inorganic_C",
                                                "Organic_C",
                                                "Total_N",
                                                "CN_ratio",
                                                "Total_S",
                                                "CS_ratio",
                                                "Olsen_P",
                                                "pH",
                                                "Clay",
                                                "Fine_Silt",
                                                "Medium_Silt",
                                                "Coarse_Silt",
                                                "Fine_Sand",
                                                "Medium_Sand",
                                                "Coarse_Sand" )
                           )
  )%>%
  arrange(var)%>%
  write.csv("./tables/mlm_posthoc_env.csv", row.names = FALSE)

###===Diversity indices====
####===Effect of ecosystem and lui =====
lme_res_div<-data.frame()
lme_post_div<-data.frame()
plots_div <-list()
for (i in 1:length(names(spec_div_rar)[c(11, 14:17)])){
  names = names(spec_div_rar)[c(11, 14:17)]
  model_div<-lmerTest::lmer(get(names[i])~Type*lui+(1|Exploratorium), data = spec_div_rar)
  print(names[i])
  anova_model_div<-as.data.frame(anova(model_div))
  anova_model_div$var<-paste(names[i])
  lme_res_div<-rbind(lme_res_div, anova_model_div)
  posthoc_div<-emmeans(model_div, c("Type", "lui"),
                   adjust = "mvt",
                   lmer.df = "satterthwaite"
  )
  outREML_div<-cld(posthoc_div,
               level = 0.05,
               Letters = letters,
               decreasing=TRUE)
  outREML_div.df<-as.data.frame(outREML_div)
  outREML_div.df$var<-paste(names[i])
  lme_post_div<-rbind(lme_post_div, outREML_div.df[ ,c(9, 1:8)])
  posthoc_div.f<-as.data.frame(posthoc_div)
  p_div<-ggplot(posthoc_div.f, aes(x=Type, y=emmean, group=lui, color=lui)) + 
    theme_bw()+
    geom_pointrange(aes(ymin=lower.CL, ymax=upper.CL), size=1.25,
                    position=position_dodge(0.7))+
    ylim(0, NA)+
    labs(x = "Ecosystem type", y = paste(names[i], "(emmean)"), color = "Land use intensity")+
    scale_color_manual(values = pal_lui)
  plots_div[[i]]<- p_div
}

lme_res_div%>%
  mutate_if(is.numeric, round, digits = 3)%>%
  rownames_to_column("Factor")%>%
  dplyr::select("var", "Factor",  "Sum Sq",  "Mean Sq", "NumDF",   "DenDF",   "F value", "Pr(>F)")%>%
  write.csv("./tables/mlm_res_div.csv", row.names = FALSE)

####=== Table with estimated marginam means (emm)====
lme_post_div%>%
  mutate_if(is.numeric, round, digits = 3)%>%
  dplyr::select(!(lower.CL:upper.CL))%>%
  mutate(lui = factor(lui, levels = c("Low", "Medium", "High")))%>%
  arrange(var, Type, lui)%>%
  pivot_wider(names_from = c(Type, lui), values_from = c(emmean:.group), names_sort = TRUE)%>%
  dplyr::select(var, emmean_Forest_Low, 
                SE_Forest_Low, 
                df_Forest_Low, 
                .group_Forest_Low, 
                emmean_Forest_Medium, 
                SE_Forest_Medium, 
                df_Forest_Medium, 
                .group_Forest_Medium,
                emmean_Forest_High, 
                SE_Forest_High, 
                df_Forest_High, 
                .group_Forest_High,
                emmean_Grassland_Low, 
                SE_Grassland_Low, 
                df_Grassland_Low, 
                .group_Grassland_Low,
                emmean_Grassland_Medium, 
                SE_Grassland_Medium, 
                df_Grassland_Medium, 
                .group_Grassland_Medium,
                emmean_Grassland_High, 
                SE_Grassland_High, 
                df_Grassland_High, 
                .group_Grassland_High)%>%
  mutate(var = factor(var, levels = c("n_spec",
                                      "shannon",
                                      "simpson",
                                      "N1",
                                      "N2",
                                      "E1",
                                      "E2")
                      )
         )%>%
  arrange(var)%>%
  write.csv("./tables/mlm_posthoc_div.csv", row.names = FALSE)





library("RColorBrewer")
col <- colorRampPalette(brewer.pal(10, "RdYlBu"))(256)


###=== Visualization ====
####====Plots for diversity =====
#' Code for Figure 1 in manuscript
p_div_1 <- plots_div[[4]]+labs(x = "Ecosystem", y = "Exp(Shannon)")
p_div_2 <- plots_div[[5]]+labs(x = "Ecosystem", y = "Simpson Diversity")
p_div_3 <- plots_div[[6]]+labs(x = "Ecosystem", y = "Evenness (q=1)")
p_div_4 <- plots_div[[7]]+labs(x = "Ecosystem", y = "Evenness (q=2)")

(p_div_1 + p_div_2 + p_div_3 + p_div_4) +
  plot_layout(
    guides = 'collect',
    axes = 'collect',
    axis_titles = 'collect'
  ) &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",  # Align legend items horizontally
    legend.margin = margin(t = -10)  # Reduce space above legend if needed
  )

ggsave("./plots/Figure_1.svg",
width = 160,
height = 160,
units = "mm")

ggsave("./plots/Figure_1.eps",
       bg = "white",
       width = 160, 
       height = 160,  # Wider than tall
       units = "mm",
       device = cairo_ps,  # Specify EPS format
       dpi = 1000)      # Adjust DPI if needed (default is 300)


##=== Correlation test =====

###=== Environment variables correlation ====
### To account for the effect of the region on the correlations, we calculated the correlation on the 
### residuals of fitted mixed effect models for each environmental variable using Spearman correlation index.
### On the full dataset
# Fit mixed models for each predictor separately
env_resid<-data.frame(row.names = rownames(site_data))


#lme_post<-data.frame()
#plots <-list()
for (i in 1:length(names(site_data)[c(12:59)])){
  names = names(site_data)[c(12:59)]
  colname<-names[i]
  model<-lmerTest::lmer(get(names[i])~(1|Exploratorium), data = site_data, na.action = na.exclude)
  print(colname)
  resid <- resid(model)
  env_resid[[colname]]<-resid
}



env_resid<-env_resid%>%
  mutate(EPS = site_data$EPS, 
         Exploratorium = site_data$Exploratorium, 
         Type = site_data$Type, 
         lui = site_data$lui, 
         lui_norm = site_data$lui_norm,
         .before = mois)

## Check for normality of residuals
names <- env_resid%>%
  dplyr::filter(Type == "Forest" & lui == "High")%>%
  dplyr::select(6:53)%>%
  colnames()

env_resid%>%
  dplyr::filter(Type == "Forest" & lui == "High")%>%
  dplyr::select(6:53)%>%
  lapply(., shapiro.test)

for (i in names){
  print(i)
  print(env_resid%>%
          ggdensity(i)
  )
  
  print(env_resid%>%
          ggqqplot(i))
  
  #print(env_resid%>%
  #dplyr::select(i)%>%
  #shapiro.test())
  
  #hist(var, col='steelblue')
}
cor_all <-env_resid%>%
  dplyr::select(c(mois:Ts_10_min, bio1:bio19, Total_C:Coarse_Sand))%>%
  scale(. , center = F, scale = T)%>%
  psych::corr.test(. ,
                   use = "pairwise.complete.obs",
                   method="spearman",
                   adjust="bonferroni",
                   alpha=.05,
                   ci=TRUE,
                   minlength=5)

tiff("./plots/testFIGS1.tif", width = 16.9, height = 16.9, units = "cm", res = 300)
corrplot(cor_all$r,
         diag = FALSE,
         type="upper",
         order="original",
         p.mat = cor_all$p,
         sig.level = 0.01,
         na.label = 'square',
         na.label.col = 'white',
         insig = "blank",
         addgrid.col = "lightgrey",
         tl.cex = .8,
         tl.srt=90,
         tl.col = "black",
         cl.cex = .8,
         title = "Overall",
         col = COL2('PuOr', 10))
dev.off()

# saved in svg (600x600)

## Forest 
### Only for forest in high lui
### High
### Only for forest in high lui
cor_for_H <-env_resid%>%
  filter(Type == "Forest" & lui == "High")%>%
  dplyr::select(c(mois:Ts_10_min, bio1:bio19, Total_C:Coarse_Sand))%>%
  scale(. , center = F, scale = T)%>%
  psych::corr.test(. ,
                   use = "pairwise.complete.obs",
                   method="spearman",
                   adjust="bonferroni",
                   alpha=.05,
                   ci=TRUE,
                   minlength=5)




## Low
### Only for forest in low lui
cor_for_L <-env_resid%>%
  filter(Type == "Forest" & lui == "Low")%>%
  dplyr::select(c(mois:Ts_10_min, bio1:bio19, Total_C:Coarse_Sand))%>%
  scale(. , center = F, scale = T)%>%
  psych::corr.test(. ,
                   use = "complete",
                   method="spearman",
                   adjust="bonferroni",
                   alpha=.05,
                   ci=TRUE,
                   minlength=5)




## Grassland
### High
### Only for grassland in high lui
cor_grs_H <-env_resid%>%
  filter(Type == "Grassland" & lui == "High")%>%
  dplyr::select(c(mois:Ts_10_min, bio1:bio19, Total_C:Coarse_Sand))%>%
  scale(. , center = F, scale = T)%>%
  psych::corr.test(. ,
                   use = "complete",
                   method="spearman",
                   adjust="bonferroni",
                   alpha=.05,
                   ci=TRUE,
                   minlength=5)

corrplot(cor_grs_H$r,
         diag = FALSE,
         type="upper",
         order="hclust",
         p.mat = cor_grs_H$p,
         sig.level = 0.01,
         na.label = 'square',
         na.label.col = 'white',
         insig = "blank",
         addgrid.col = "lightgrey",
         tl.cex = .8,
         tl.srt=90,
         tl.col = "black",
         cl.cex = .8,
         title = "Grassland High",
         col = COL2('PuOr', 10))


### Low
### Only for grassland in low lui
cor_grs_L <-env_resid%>%
  filter(Type == "Grassland" & lui == "Low")%>%
  dplyr::select(c(mois:Ts_10_min, bio1:bio19, Total_C:Coarse_Sand))%>%
  scale(. , center = F, scale = T)%>%
  psych::corr.test(. ,
                   use = "complete",
                   method="spearman",
                   adjust="bonferroni",
                   alpha=.05,
                   ci=TRUE,
                   minlength=5)

corrplot(cor_grs_L$r,
         diag = FALSE,
         type="upper",
         order="hclust",
         p.mat = cor_grs_L$p,
         sig.level = 0.01,
         na.label = 'square',
         na.label.col = 'white',
         insig = "blank",
         addgrid.col = "lightgrey",
         tl.cex = .8,
         tl.srt=90,
         tl.col = "black",
         cl.cex = .8,
         title = "Grassland Low",
         col = COL2('PuOr', 10))

# Correlation with diversity indices
## Preparing dataset for correlation
envdiv_data<-merge(spec_div_rar, site_data, by = "EPS")%>%
  dplyr::rename(lui = lui.x,
                lui_norm = lui_norm.x,
                Exploratorium = Exploratorium.x,
                Type = Type.x, 
                Elevation = Elevation.x, 
                Elevation.sd = `Elevation sd.x`, 
                Slope.d = `Slope d.x`, 
                Slope.sd = `Slope sd.x`, 
                Slope.p = `Slope p.x`)%>%
  dplyr::select(EPS:Type, n_spec:E2, mois:Coarse_Sand)

# Fit mixed models for each predictor separately
div_resid<-data.frame(row.names = rownames(envdiv_data))
#lme_post<-data.frame()
#plots <-list()
for (i in 1:length(names(envdiv_data)[c(6:60)])){
  names = names(envdiv_data)[c(6:60)]
  colname<-names[i]
  model<-lmerTest::lmer(get(names[i])~(1|Exploratorium), data = envdiv_data, na.action = na.exclude)
  print(colname)
  resid <- resid(model)
  div_resid[[colname]]<-resid
}

div_resid<-div_resid%>%
  mutate(EPS = envdiv_data$EPS, 
         Exploratorium = envdiv_data$Exploratorium, 
         Type = envdiv_data$Type, 
         lui = envdiv_data$lui, 
         lui_norm = envdiv_data$lui_norm,
         .before = n_spec)

#Correlation on the full dataset
corr_div<-div_resid%>%
  #filter(Type == "Grassland" & lui == "High")%>%
  dplyr::select(c(n_spec:Coarse_Sand))%>%
  scale(. , center = F, scale = T)%>%
  psych::corr.test(. ,
                   use = "complete",
                   method="spearman",
                   adjust="bonferroni",
                   alpha=.05,
                   ci=TRUE,
                   minlength=5)


#Correlation in Forest:High (FH)

corr_div_FH<-div_resid%>%
  filter(Type == "Forest" & lui == "High")%>%
  dplyr::select(c(n_spec:Coarse_Sand))%>%
  scale(. , center = F, scale = T)%>%
  psych::corr.test(. ,
                   use = "complete",
                   method="spearman",
                   adjust="bonferroni",
                   alpha=.05,
                   ci=TRUE,
                   minlength=5)


#Correlation in Forest:Low (FL)
corr_div_FL<-div_resid%>%
  filter(Type == "Forest" & lui == "Low")%>%
  dplyr::select(c(n_spec:Coarse_Sand))%>%
  scale(. , center = F, scale = T)%>%
  psych::corr.test(. ,
                   use = "complete",
                   method="spearman",
                   adjust="bonferroni",
                   alpha=.05,
                   ci=TRUE,
                   minlength=5)


#Correlation in Grassland:High (GH)

corr_div_GH<-div_resid%>%
  filter(Type == "Grassland" & lui == "High")%>%
  dplyr::select(c(n_spec:Coarse_Sand))%>%
  scale(. , center = F, scale = T)%>%
  psych::corr.test(. ,
                   use = "complete",
                   method="spearman",
                   adjust="bonferroni",
                   alpha=.05,
                   ci=TRUE,
                   minlength=5)


#Correlation in Grassland:Low (FL)
corr_div_GL<-div_resid%>%
  filter(Type == "Grassland" & lui == "Low")%>%
  dplyr::select(c(n_spec:Coarse_Sand))%>%
  scale(. , center = F, scale = T)%>%
  psych::corr.test(. ,
                   use = "complete",
                   method="spearman",
                   adjust="bonferroni",
                   alpha=.05,
                   ci=TRUE,
                   minlength=5)


corrplot(cbind(corr_div$r[ -c(1:7),c(1,2,6)], corr_div_FH$r[ -c(1:7),c(1,2,6)], corr_div_FL$r[ -c(1:7),c(1,2,6)], corr_div_GH$r[ -c(1:7),c(1,2,6)], corr_div_GL$r[ -c(1:7),c(1,2,6)]),
         diag = TRUE,
         is.corr = TRUE,
         type="full",
         order="original",
         p.mat = cbind(corr_div$p[ -c(1:7),c(1,2,6)], corr_div_FH$p[ -c(1:7),c(1,2,6)],corr_div_FL$p[ -c(1:7),c(1,2,6)], corr_div_GH$p[ -c(1:7),c(1,2,6)],corr_div_GL$p[ -c(1:7),c(1,2,6)]),
         sig.level = 0.01,
         na.label = 'square',
         na.label.col = 'white',
         insig = "blank",
         addgrid.col = "lightgrey",
         tl.cex = .8,
         tl.srt=90,
         tl.col = "black",
         cl.cex = .8,
         col = COL2('PuOr', 10))


###=== Correlation with diversity indices====
## Preparing dataset for correlation
envdiv_data<-merge(spec_div_rar, site_data, by = "EPS")%>%
  dplyr::rename(lui = lui,
                lui_norm = lui_norm,
                Exploratorium = Exploratorium,
                Type = Type, 
                Elevation = Elevation, 
                Elevation.sd = `Elevation sd`, 
                Slope.d = `Slope d`, 
                Slope.sd = `Slope sd`, 
                Slope.p = `Slope p`)%>%
  dplyr::select(EPS, Exploratorium, Type, lui, lui_norm, Type_lui, n_spec:E2, mois:Ts_10_min, bio1:bio20, Total_C:Coarse_Sand)

# Fit mixed models for each predictor separately
div_resid<-data.frame(row.names = rownames(envdiv_data))

for (i in 1:length(names(envdiv_data)[c(7:52)])){
  names = names(envdiv_data)[c(7:52)]
  colname<-names[i]
  model<-lmerTest::lmer(get(names[i])~(1|Exploratorium), data = envdiv_data, na.action = na.exclude)
  print(colname)
  resid <- resid(model)
  div_resid[[colname]]<-resid
}

div_resid<-div_resid%>%
  mutate(EPS = envdiv_data$EPS, 
         Exploratorium = envdiv_data$Exploratorium, 
         Type = envdiv_data$Type, 
         lui = envdiv_data$lui, 
         lui_norm = envdiv_data$lui_norm,
         .before = n_spec)

#' Correlation on the full dataset
corr_div<-div_resid%>%
  #filter(Type == "Grassland" & lui == "High")%>%
  dplyr::select(c(n_spec:Coarse_Sand))%>%
  scale(. , center = F, scale = T)%>%
  psych::corr.test(. ,
                   use = "complete",
                   method="spearman",
                   adjust="bonferroni",
                   alpha=.05,
                   ci=TRUE,
                   minlength=5)


#Correlation in Forest:High (FH)

corr_div_FH<-div_resid%>%
  filter(Type == "Forest" & lui == "High")%>%
  dplyr::select(c(n_spec:Coarse_Sand))%>%
  scale(. , center = F, scale = T)%>%
  psych::corr.test(. ,
                   use = "complete",
                   method="spearman",
                   adjust="bonferroni",
                   alpha=.05,
                   ci=TRUE,
                   minlength=5)


#Correlation in Forest:Low (FL)
corr_div_FL<-div_resid%>%
  filter(Type == "Forest" & lui == "Low")%>%
  dplyr::select(c(n_spec:Coarse_Sand))%>%
  scale(. , center = F, scale = T)%>%
  psych::corr.test(. ,
                   use = "complete",
                   method="spearman",
                   adjust="bonferroni",
                   alpha=.05,
                   ci=TRUE,
                   minlength=5)


#Correlation in Grassland:High (GH)

corr_div_GH<-div_resid%>%
  filter(Type == "Grassland" & lui == "High")%>%
  dplyr::select(c(n_spec:Coarse_Sand))%>%
  scale(. , center = F, scale = T)%>%
  psych::corr.test(. ,
                   use = "complete",
                   method="spearman",
                   adjust="bonferroni",
                   alpha=.05,
                   ci=TRUE,
                   minlength=5)


#Correlation in Grassland:Low (FL)
corr_div_GL<-div_resid%>%
  filter(Type == "Grassland" & lui == "Low")%>%
  dplyr::select(c(n_spec:Coarse_Sand))%>%
  scale(. , center = F, scale = T)%>%
  psych::corr.test(. ,
                   use = "complete",
                   method="spearman",
                   adjust="bonferroni",
                   alpha=.05,
                   ci=TRUE,
                   minlength=5)

###=== Visualization =====
####=== Evironmental variables by ecosystem ====
#' General correlation plot between environmental variables
ggcorrplot(cor_all$r,
           p.mat = cor_all$p,
           sig.level = 0.01,
           insig = 'blank',
           method = "circle",
           type = "lower",
           outline.color = "black",
           lab_size = 6) +
  labs(x = NULL, y = NULL) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave("./plots/Figure_S1.eps",
       bg = "white",
       width = 400,
       height = 400,  # Wider than tall
       units = "mm",
       device = cairo_ps,  # Specify EPS format
       dpi = 1000)      # Adjust DPI if needed (default is 300)

#' Plots for correlations in the single combinations of ecosystem and land use intensity
cor_1 <- ggcorrplot(cor_for_H$r,
                    p.mat = cor_for_H$p,
                    sig.level = 0.01,
                    insig = 'blank',
                    method = "circle",
                    type = "lower",
                    outline.color = "black",
                    lab_size = 6,
                  title = "Forest \nHigh LUI") +
  labs(x = NULL, y = NULL) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

cor_2 <- ggcorrplot(cor_for_L$r,
                    p.mat = cor_for_L$p,
                    sig.level = 0.01,
                    insig = 'blank',
                    method = "circle",
                    type = "lower",
                    outline.color = "black",
                    lab_size = 6,
                    title = "Forest \nLow LUI") +
  labs(x = NULL, y = NULL) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

cor_3 <- ggcorrplot(cor_grs_H$r,
                    p.mat = cor_grs_H$p,
                    sig.level = 0.01,
                    insig = 'blank',
                    method = "circle",
                    type = "lower",
                    outline.color = "black",
                    lab_size = 6,
                    title = "Grassland \nHigh LUI") +
  labs(x = NULL, y = NULL) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

cor_4 <- ggcorrplot(cor_grs_L$r,
                    p.mat = cor_grs_L$p,
                    sig.level = 0.01,
                    insig = 'blank',
                    method = "circle",
                    type = "lower",
                    outline.color = "black",
                    lab_size = 6,
                    title = "Grassland \nLow LUI") +
  labs(x = NULL, y = NULL) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

cor_1+cor_2+cor_3+cor_4+plot_layout(guides = 'collect',axis_titles = 'collect')

ggsave("./plots/Figure_2.eps",
bg = "white",
width = 700,
height = 700,  # Wider than tall
units = "mm",
device = cairo_ps,  # Specify EPS format
dpi = 1000)      # Adjust DPI if needed (default is 300)


cor_div_1<- ggcorrplot(corr_div$r[ -c(1:6),c(4,5,6,7)],
                       p.mat = corr_div$p[ -c(1:6),c(4,5,6,7)],
                       sig.level = 0.01,
                       insig = 'blank',
                       method = "circle",
                       type = "lower",
                       outline.color = "black",
                       lab_size = 6,
                       title = "Overall") +
  labs(x = NULL, y = NULL) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

cor_div_2<- ggcorrplot(corr_div_FH$r[ -c(1:6),c(4,5,6,7)],
                       p.mat = corr_div_FH$p[ -c(1:6),c(4,5,6,7)],
                       sig.level = 0.01,
                       insig = 'blank',
                       method = "circle",
                       type = "lower",
                       outline.color = "black",
                       lab_size = 6,
                       title = "Forest - High LUI") +
  labs(x = NULL, y = NULL) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

cor_div_3<- ggcorrplot(corr_div_FL$r[ -c(1:6),c(4,5,6,7)],
                       p.mat = corr_div_FL$p[ -c(1:6),c(4,5,6,7)],
                       sig.level = 0.01,
                       insig = 'blank',
                       method = "circle",
                       type = "lower",
                       outline.color = "black",
                       lab_size = 6,
                       title = "Forest - Low LUI") +
  labs(x = NULL, y = NULL) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

cor_div_4<- ggcorrplot(corr_div_GH$r[ -c(1:6),c(4,5,6,7)],
                       p.mat = corr_div_GH$p[ -c(1:6),c(4,5,6,7)],
                       sig.level = 0.01,
                       insig = 'blank',
                       method = "circle",
                       type = "lower",
                       outline.color = "black",
                       lab_size = 6,
                       title = "Grassland - High LUI") +
  labs(x = NULL, y = NULL) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

cor_div_5<- ggcorrplot(corr_div_GL$r[ -c(1:6),c(4,5,6,7)],
                       p.mat = corr_div_GL$p[ -c(1:6),c(4,5,6,7)],
                       sig.level = 0.01,
                       insig = 'blank',
                       method = "circle",
                       type = "lower",
                       outline.color = "black",
                       lab_size = 6,
                       title = "Grassland - Low LUI") +
  labs(x = NULL, y = NULL) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

cor_div_5/cor_div_4/cor_div_3/cor_div_2/cor_div_1+plot_layout(guides = 'collect', axis_titles = 'collect')+ plot_annotation(tag_levels = 'A')

ggsave("./plots/Figure_5.eps",
       bg = "white",
       width = 700,
       height = 700,  # Wider than tall
       units = "mm",
       device = cairo_ps,  # Specify EPS format
       dpi = 1000)      # Adjust DPI if needed (default is 300)


corrplot(cbind(corr_div$r[ -c(1:7),c(4,5,6,7)], 
               corr_div_FH$r[ -c(1:7),c(4,5,6,7)], 
               corr_div_FL$r[ -c(1:7),c(4,5,6,7)], 
               corr_div_GH$r[ -c(1:7),c(4,5,6,7)], 
               corr_div_GL$r[ -c(1:7),c(4,5,6,7)]),
         diag = TRUE,
         is.corr = TRUE,
         type="full",
         order="original",
         p.mat = cbind(corr_div$p[ -c(1:7),c(4,5,6,7)], 
                       corr_div_FH$p[ -c(1:7),c(4,5,6,7)],
                       corr_div_FL$p[ -c(1:7),c(4,5,6,7)], 
                       corr_div_GH$p[ -c(1:7),c(4,5,6,7)],
                       corr_div_GL$p[ -c(1:7),c(4,5,6,7)]),
         sig.level = 0.01,
         na.label = 'square',
         na.label.col = 'white',
         insig = "blank",
         addgrid.col = "lightgrey",
         tl.cex = .8,
         tl.srt=90,
         tl.col = "black",
         cl.cex = .8,
         col = COL2('PuOr', 10))

##=== Differential abundance analysis ====
spec_ds<-DESeqDataSetFromMatrix(countData = t(spec_abund),
                                colData = site_data[c(1:4, 6:58)],
                                design=~1
)

spec_ds$group<-factor(paste(spec_ds$Type, spec_ds$lui))

design(spec_ds)<-~group 

spec_ds<-DESeq(spec_ds, test="Wald", fitType="parametric", sfType = "poscounts")

res_for<- results(spec_ds, 
                  contrast = c("group", "Forest High", "Forest Low"  ), 
                  #name = "group_Grassland.Low_vs_Forest.Low", 
                  tidy = TRUE)


res_for_n<-merge(res_for, tax_mat, by.x="row", by.y=0)%>%
  dplyr::filter(padj <0.01)

res_grs<- results(spec_ds, 
                  contrast = c("group", "Grassland High", "Grassland Low"  ), 
                  #name = "group_Grassland.Low_vs_Forest.Low", 
                  tidy = TRUE)

res_grs_n<-merge(res_grs, tax_mat, by.x="row", by.y=0)%>%
  dplyr::filter(padj <0.01)

###==== Exporting tables summaries ====
res_n%>%
  dplyr::select(Kingdom, Phylum, Class, Order, Family, Genus, Species, Guild, Guild, Secondary_lifestyle, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj)%>%
  group_by(Order, Family, Genus)%>%
  summarize(baseMean_avg = mean(baseMean, na.rm = T), across(log2FoldChange:padj, max))%>%
  filter(abs(log2FoldChange)>10 & padj<0.01)%>%
  arrange(., desc(log2FoldChange))%>%
  write.csv("./tables/differential_ecosystem.csv", row.names = FALSE)

res_for_n%>%
  dplyr::select(Kingdom, Phylum, Class, Order, Family, Genus, Species, Guild, Guild, Secondary_lifestyle, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj)%>%
  #group_by(Order, Family, Genus)%>%
  #summarize(baseMean_avg = mean(baseMean, na.rm = T), across(log2FoldChange:padj, max))%>%
  filter(abs(log2FoldChange)>20 & padj<0.01)%>%
  arrange(., desc(log2FoldChange))%>%
  write.csv("./tables/differential_forest_land_use.csv", row.names = FALSE)

res_grs_n%>%
  dplyr::select(Kingdom, Phylum, Class, Order, Family, Genus, Species, Guild, Guild, Secondary_lifestyle, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj)%>%
  #group_by(Order, Family, Genus)%>%
  #summarize(baseMean_avg = mean(baseMean, na.rm = T), across(log2FoldChange:padj, max))%>%
  filter(abs(log2FoldChange)>20 & padj<0.01)%>%
  arrange(., desc(log2FoldChange))%>%
  write.csv("./tables/differential_grassland_land_use.csv", row.names = FALSE)

###=== Visualization ======

ggmaplot(res_n, main = expression("Forest High" %->% "Forest Low"),
         fdr = 0.05, fc = 20, size = 1.5,
         genenames = as.vector(res_n$Species),
         select.top.method = c("fc"),
         legend = "top", top = 20,
         font.label = c("bold", 11, "grey75"),
         font.legend = "bold",
         font.main = "bold",
         ggtheme = ggplot2::theme_minimal())+
  scale_color_manual(labels = c("Forest High", "Forest Low", "Not different"), 
                     values = c("palegreen4", "tan4", "grey75"))+ labs(color = "Ecosystem type")

res_for_n%>%
  ggmaplot(., main = expression("Forest:High" %->% "Forest:Low"),
           fdr = 0.05, fc = 20, size = 1.5,
           palette = c("#BC3C29FF", "#E18727FF", "#9498A0FF"),
           genenames = as.vector(res_for_n$Species),
           select.top.method = c("fc"),
           legend = "top", top = 40,
           font.label = c("bold", 11, "#9498A0FF"),
           font.legend = "bold",
           font.main = "bold",
           ggtheme = ggplot2::theme_minimal())+
  xlim(0,NA)+
  scale_color_manual(labels = c("Forest High", "Forest Low", "Not different"), values = c("#BC3C29FF", "#E18727FF", "#9498A0FF"))+
  labs(color = "Land use intensity")

ggmaplot(res_grs_n, main = expression("Grassland:High" %->% "Grassland:Low"),
         fdr = 0.05, fc = 20, size = 1.5,
         #palette = c("red2", "steelblue", "grey50"),
         genenames = as.vector(res_grs_n$Species),
         select.top.method = c("fc"),
         legend = "top", top = 40,
         font.label = c("bold", 11, "grey75"),
         font.legend = "bold",
         font.main = "bold",
         ggtheme = ggplot2::theme_minimal())+
  xlim(0,NA)+
  scale_color_manual(labels = c("High", "Low", "Not different"), values = c("#20854EFF", "#3CA951FF", "#9498A0FF"))+
  labs(color = "Land use intensity")

res_n%>%
  group_by(Order)%>%
  summarize(across(log2FoldChange, max))%>%
  arrange(., desc(log2FoldChange))%>%
  mutate(color = case_when(log2FoldChange>0 ~ "Forest", 
                           log2FoldChange<0 ~ "Grassland"))%>%
  ggbarplot(x = "Order", y = "log2FoldChange", orientation = "horiz", position = position_dodge(), fill = "color", palette =c("#BC3C29FF", "#EFB118FF"))+
  labs(y = "Fold Change (log2)", x = "Order", fill = "Ecosystem type")+
  ylim(-35, 35)+
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))+
  theme_bw()+
  theme(legend.position = "bottom")

a<-res_for_n%>%
  group_by(Genus)%>%
  summarize(across(log2FoldChange, max))%>%
  mutate(color = case_when(log2FoldChange>0 ~ "Forest:High", 
                           log2FoldChange<0 ~ "Forest:Low"))%>%
  merge(.,
        tax_mat[c(7, 9)], 
        by = "Genus",
        all.x = FALSE,
        all.y = FALSE)%>%
  unique()%>%
  group_by(Guild)%>%
  arrange(., desc(log2FoldChange), desc(Guild), .by_group=TRUE)%>%
  filter(Guild!="unknown")%>%
  ggbarplot(x = "Genus", y = "log2FoldChange", orientation = "horiz", position = position_dodge(), fill = "Guild")+
  labs(y = "Fold Change (log2)", x = "Genus")+
  ylim(-30, 30)+
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))+
  scale_fill_manual(name = "Trophic Guild", values = pal_guild, aesthetics = c("color", "fill"))+
  theme_minimal()+
  theme(legend.position = "bottom")

b<-res_grs_n%>%
  group_by(Genus)%>%
  summarize(across(log2FoldChange, max))%>%
  mutate(color = case_when(log2FoldChange>0 ~ "Grassland:High", 
                           log2FoldChange<0 ~ "Grassland:Low"))%>%
  merge(.,
        tax_mat[c(7, 9)], 
        by = "Genus",
        all.x = FALSE,
        all.y = FALSE)%>%
  unique()%>%
  group_by(Guild)%>%
  arrange(., desc(log2FoldChange), Guild, .by_group=TRUE)%>%
  filter(Guild!="unknown")%>%
  ggbarplot(x = "Genus", y = "log2FoldChange", orientation = "horiz", position = position_dodge(), fill = "Guild")+
  labs(y = "Fold Change (log2)", x = "Genus", fill = "Use intensity")+
  ylim(-40, 40)+
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))+
  scale_fill_manual(name = "Trophic Guild", values = pal_guild, aesthetics = c("color", "fill"))+
  theme_minimal()+
  theme(legend.position = "bottom")

ggarrange(a, b,
          ncol = 2,
          common.legend = TRUE,
          legend = "bottom",
          labels = c("A", "B"))


##=== Distance based redundancy analysis (db-RDA) ====

db_rda<-capscale(formula= spec_abund_rar ~ Type + lui_norm + Exploratorium +  bio1+bio2+bio4+bio5+bio6+bio7+bio8+bio9+bio10+bio11+bio12+bio13+bio14+bio15+bio16+bio17+bio18+bio19+bio20+Total_C+ Inorganic_C+Organic_C+Total_N+CN_ratio+Total_S+CS_ratio+Olsen_P+pH+Clay+Fine_Silt+Medium_Silt+Coarse_Silt+Fine_Sand+Medium_Sand+Coarse_Sand,
                 data = site_data_scaled,
                 distance ="bray", 
                 sqrt.dist = FALSE,
                 na.action=na.exclude)

db_rda_f<-capscale(formula= spec_abund_rar[c(51:100, 151:200, 251:300), ] ~ lui_norm + Exploratorium + bio1+bio2+bio4+bio5+bio6+bio7+bio8+bio9+bio10+bio11+bio12+bio13+bio14+bio15+bio16+bio17+bio18+bio19+bio20+Total_C+ Inorganic_C+Organic_C+Total_N+CN_ratio+Total_S+CS_ratio+Olsen_P+pH+Clay+Fine_Silt+Medium_Silt+Coarse_Silt+Fine_Sand+Medium_Sand+Coarse_Sand,
                   data = site_data_scaled %>% filter(Type == "Forest"),
                   dist ="bray", 
                   sqrt.dist = FALSE,
                   na.action=na.exclude)

db_rda_g<-capscale(formula= spec_abund_rar[c(1:50, 101:150, 201:250), ] ~ lui_norm + Exploratorium + bio1+bio2+bio4+bio5+bio6+bio7+bio8+bio9+bio10+bio11+bio12+bio13+bio14+bio15+bio16+bio17+bio18+bio19+bio20+Total_C+ Inorganic_C+Organic_C+Total_N+CN_ratio+Total_S+CS_ratio+Olsen_P+pH+Clay+Fine_Silt+Medium_Silt+Coarse_Silt+Fine_Sand+Medium_Sand+Coarse_Sand,
                   data = site_data_scaled %>% filter(Type == "Grassland"),
                   dist ="bray", 
                   sqrt.dist = FALSE,
                   na.action=na.exclude)

###=== Selecting variables ====

fwd.sel <- ordiR2step(capscale(formula= spec_abund_rar ~ 1,
                                   data = site_data_scaled[c(2, 3, 5, 11:20, 22:58)],
                                   distance ="bray", sqrt.dist = FALSE,
                                   na.action=na.exclude), # lower model limit (simple!)
                          scope = formula(db_rda), # upper model limit (the "full" model)
                          direction = "forward",
                          R2scope = TRUE, # can't surpass the "full" model's R2
                          pstep = 1000,
                          trace = TRUE) # change to TRUE to see the selection process!

fwd.sel_f <- ordiR2step(capscale(formula= spec_abund_rar[c(51:100, 151:200, 251:300), ] ~ 1,
                                 data = site_data_scaled[c(2, 3, 5, 11:20, 22:58)]%>% filter(Type == "Forest"),
                                 distance ="bray", 
                                 sqrt.dist = FALSE,
                                 na.action=na.exclude), # lower model limit (simple!)
                        scope = formula(db_rda_f), # upper model limit (the "full" model)
                        direction = "forward",
                        R2scope = TRUE, # can't surpass the "full" model's R2
                        pstep = 1000,
                        trace = TRUE) # change to TRUE to see the selection process!

fwd.sel_g <- ordiR2step(capscale(formula= spec_abund_rar[c(1:50, 101:150, 201:250), ] ~ 1,
                                 data = site_data_scaled[c(2, 3, 5, 11:20, 22:58)] %>% filter(Type == "Grassland"),
                                 distance ="bray", 
                                 sqrt.dist = FALSE,
                                 na.action=na.exclude), # lower model limit (simple!)
                        scope = formula(db_rda_g), # upper model limit (the "full" model)
                        direction = "forward",
                        R2scope = TRUE, # can't surpass the "full" model's R2
                        pstep = 1000,
                        trace = TRUE) # change to TRUE to see the selection process!

### The RDA model is rewritten using only the selected variables

fwd.sel$call

db_rda_sel<-capscale(formula = spec_abund_rar ~ Type + Medium_Sand + Inorganic_C + 
                       pH + CN_ratio + Total_N + Elevation + lui_norm + Total_C + 
                       Coarse_Silt + Olsen_P + bio13 + Fine_Sand + Fine_Silt + bio16 + 
                       Total_S + CS_ratio + bio12 + bio18 + bio6 + bio19 + bio10 + 
                       Medium_Silt + Clay + bio17, 
                     data = site_data_scaled[-21],
                     distance = "bray", 
                     sqrt.dist = FALSE, 
                     na.action = na.exclude)

db_rda_sel<-capscale(formula = spec_abund_rar ~ Type + Exploratorium + pH + 
           CN_ratio + Inorganic_C + Total_N + lui_norm + Coarse_Silt + 
           Total_C + Olsen_P + bio13 + Fine_Sand + bio16 + Total_S + 
           bio6 + bio12 + bio18 + Medium_Sand + CS_ratio + Fine_Silt + 
           bio19 + bio10 + bio15, data = site_data_scaled[c(2, 3, 5, 11:20, 22:58)], 
           distance = "bray", 
           sqrt.dist = FALSE, 
           na.action = na.exclude)

fwd.sel_f$call

db_rda_sel_f<-capscale(formula = spec_abund_rar[c(51:100, 151:200, 251:300), ] ~ Medium_Sand + pH + lui_norm + CN_ratio + Elevation + 
  Fine_Sand + Coarse_Silt + Olsen_P + Fine_Silt + Total_S +  Total_C + Total_N + CS_ratio + Inorganic_C + Medium_Silt + 
  bio16, data = site_data_scaled[-21] %>% filter(Type == "Forest"), 
distance = "bray", sqrt.dist = FALSE, na.action = na.exclude)

db_rda_sel_f<-capscale(formula = spec_abund_rar[c(51:100, 151:200, 251:300),] ~ Exploratorium + pH + CN_ratio + lui_norm + Coarse_Silt + 
  Olsen_P + Total_S + Fine_Silt + Total_C + Total_N + Coarse_Sand + 
  Inorganic_C + Medium_Silt, data = site_data_scaled[c(2, 3,5, 11:20, 22:58)] %>% filter(Type == "Forest"), 
distance = "bray", 
sqrt.dist = FALSE, 
na.action = na.exclude)


fwd.sel_g$call

db_rda_sel_g<-capscale(formula = spec_abund_rar[c(1:50, 101:150, 201:250), ] ~ Inorganic_C + pH + Medium_Sand + Total_N + lui_norm + 
  Elevation + Organic_C + Olsen_P + Total_S + CN_ratio + bio12 +  Fine_Sand + Coarse_Silt, data = site_data_scaled %>% filter(Type =="Grassland") %>% dplyr::select(-c(bio3)), distance = "bray", 
sqrt.dist = FALSE, na.action = na.exclude)

db_rda_sel_g<-capscale(formula = spec_abund_rar[c(1:50, 101:150, 201:250),] ~ Exploratorium + Inorganic_C + lui_norm + pH + Total_N + 
  CN_ratio + Olsen_P + bio13 + Total_S + Coarse_Silt + Medium_Sand + 
  Fine_Sand + Organic_C, data = site_data_scaled[c(2, 3, 5, 11:20, 22:58)] %>% filter(Type == "Grassland"), distance = "bray", 
sqrt.dist = FALSE, na.action = na.exclude)


RsquareAdj(db_rda)
RsquareAdj(db_rda_sel)

RsquareAdj(db_rda_f)
RsquareAdj(db_rda_sel_f)


RsquareAdj(db_rda_g)
RsquareAdj(db_rda_sel_g)


###==== Partitioning of the variance =====
env.soil <- subset(site_data_scaled, select = c(Total_C, Inorganic_C, Organic_C, Total_N, CN_ratio, Total_S, CS_ratio, Olsen_P, pH, Clay, Fine_Silt, Medium_Silt, Coarse_Silt, Fine_Sand, Medium_Sand, Coarse_Sand))
env.temp <- subset(site_data_scaled, select = c(bio1, bio2,  bio4, bio5, bio6, bio7, bio8, bio9, bio10))
env.mois <- subset(site_data_scaled, select = c(bio11, bio12, bio13, bio14, bio15, bio16, bio17, bio18, bio19, bio20))
env.lui <- subset(site_data_scaled, select = lui_norm)
# Partition the variation in community composition
spec.part.all <- varpart(as.matrix(vegdist(spec_abund_rar, "bray")), env.soil, env.temp, env.mois, env.lui)
spec.part.all$part  # access results!

### Partitioning of the variance
spec.part.for <- subset( site_data_scaled %>% filter(Type =="Forest"), select = c(Total_C, Inorganic_C, Organic_C, Total_N, CN_ratio, Total_S, CS_ratio, Olsen_P, pH, Clay, Fine_Silt, Medium_Silt, Coarse_Silt, Fine_Sand, Medium_Sand, Coarse_Sand))
env.temp.for <- subset(site_data_scaled%>% filter(Type =="Forest"), select = c(bio1, bio2,  bio4, bio5, bio6, bio7, bio8, bio9, bio10))
env.mois.for <- subset(site_data_scaled%>% filter(Type =="Forest"), select = c(bio11, bio12, bio13, bio14, bio15, bio16, bio17, bio18, bio19, bio20))
env.lui.for <- subset(site_data_scaled%>% filter(Type =="Forest"), select = lui_norm)
# Partition the variation in community composition
spec.part.for <- varpart(as.matrix(vegdist(spec_abund_rar[c(51:100, 151:200, 251:300), ], "bray")), spec.part.for, env.temp.for, env.mois.for, env.lui.for)
spec.part.for$part  # access results!

### Partitining of the variance
spec.part.grs <- subset( site_data_scaled %>% filter(Type =="Grassland"), select = c(Total_C, Inorganic_C, Organic_C, Total_N, CN_ratio, Total_S, CS_ratio, Olsen_P, pH, Clay, Fine_Silt, Medium_Silt, Coarse_Silt, Fine_Sand, Medium_Sand, Coarse_Sand))
env.temp.grs <- subset(site_data_scaled%>% filter(Type =="Grassland"), select = c(bio1, bio2,  bio4, bio5, bio6, bio7, bio8, bio9, bio10))
env.mois.grs <- subset(site_data_scaled%>% filter(Type =="Grassland"), select = c(bio11, bio12, bio13, bio14, bio15, bio16, bio17, bio18, bio19, bio20))
env.lui.grs <- subset(site_data_scaled%>% filter(Type =="Grassland"), select = lui_norm)
# Partition the variation in community composition
spec.part.grs <- varpart(as.matrix(vegdist(spec_abund_rar[c(1:50, 101:150, 201:250), ], "bray")), spec.part.grs, env.temp.grs, env.mois.grs, env.lui.grs)
spec.part.grs$part  # access results!

###=== ANOVA of the RDA with selected variables====
####=== ANOVA TEST on the general RDA =====
anova.cca(db_rda_sel) # overall test of the significant of the analysis
anova.cca(db_rda_sel, by="axis", perm.max=999) # test axes for significance
anova.cca(db_rda_sel, by="terms", perm.max=999) %>% # test for sign. environ. variables
  kbl(caption = "ANOVA of dbRDA environmental terms", digits = 3) %>%
  kable_classic(full_width = F, html_font = "Calibri")%>%
  save_kable(. , "./tables/anova_rda_env_terms.docx")
anova.cca(db_rda_sel,
          by="margin",
          perm.max=999)%>%
  as.data.frame()%>%
  write.csv("./tables/anovaIII_rda_env_terms.csv", row.names = TRUE)


####=== ANOVA TEST on the RDA of the forest data ====
anova.cca(db_rda_sel_f) # overall test of the significant of the analysis
anova.cca(db_rda_sel_f, by="axis", perm.max=999) # test axes for significance
anova.cca(db_rda_sel_f, by="terms", perm.max=999) %>% # test for sign. environ. variables
  kbl(caption = "ANOVA of dbRDA environmental terms in forest", digits = 3) %>%
  kable_classic(full_width = F, html_font = "Calibri")%>%
  save_kable(. , "./tables/anova_rda_env_terms_forest.docx")
anova.cca(db_rda_sel_f, 
          by="margin", 
          perm.max=999)%>%
  as.data.frame()%>%
  write.csv("./tables/anovaIII_rda_env_terms_forest.csv")

####=== ANOVA TEST on the RDA of the grassland data ====
anova.cca(db_rda_sel_g) # overall test of the significant of the analysis
anova.cca(db_rda_sel_g, by="axis", perm.max=999) # test axes for significance
anova.cca(db_rda_sel_g, by="terms", perm.max=999) %>% # test for sign. environ. variables
  kbl(caption = "ANOVA of dbRDA environmental terms in grassland", digits = 3) %>%
  kable_classic(full_width = F, html_font = "Calibri")%>%
  save_kable(. , "./tables/anova_rda_env_terms_grassland.docx")
anova.cca(db_rda_sel_g, 
          by="margin", 
          perm.max=999)%>% as.data.frame()%>%
  write.csv(. , "./tables/anovaIII_rda_env_terms_grassland.csv")

###=== Visualization =====
#' Plotting variance partitioning
plot(spec.part.all,
     Xnames = c("Soil", "Temp", "Mois", "LUI"), # name the partitions
     bg = c("chartreuse3", "mediumpurple", "darkorange", "darkgoldenrod2"), alpha = 80, # colour the circles
     digits = 2, # only show 2 digits
     cex = 1.5)

plot(spec.part.for,
     Xnames = c("Soil", "Temp", "Mois", "LUI"), # name the partitions
     bg = c("chartreuse3", "mediumpurple", "darkorange", "darkgoldenrod2"), alpha = 80, # colour the circles
     digits = 2, # only show 2 digits
     cex = 1.5)

plot(spec.part.grs,
     Xnames = c("Soil", "Temp", "Mois", "LUI"), # name the partitions
     bg = c("chartreuse3", "mediumpurple", "darkorange", "darkgoldenrod2"), alpha = 80, # colour the circles
     digits = 2, # only show 2 digits
     cex = 1.5)

#' Plotting the result of the dbRDA with ggpubr
plot<-ordiplot(db_rda_sel, scaling = 2) #scaling 2 is used to have angles between variables reflect their correlation

perc <- round(100*(summary(db_rda_sel)$cont$importance[2, 1:2]), 2)

uscores = as.data.frame(plot$sites)
vscores = as.data.frame(plot$species)
escores = as.data.frame(plot$biplot)
uscores1 <- dplyr::inner_join(site_data_scaled, rownames_to_column(data.frame(uscores), var = "EPS"), by = "EPS")

ordi_oa<-ggplot(uscores1) + 
  geom_point(aes(x = CAP1, y = CAP2, fill = lui, 
                 shape = Type), color = 'black', size = 3) +
  geom_point(data = vscores, aes(x = CAP1, y = CAP2), col = "darkgrey", shape = 3)+
  #scale_color_brewer(palette = "Set1") +
  scale_color_manual(values = pal_lui, 
                     aesthetics = c("colour", "fill"),
                     breaks = c("Low", "Medium", "High"),
                     labels = c("Low", "Medium", "High")) +
  scale_shape_manual(values = c(22, 24, 21)) +
  geom_segment(data = escores, aes(x = 0, y = 0, xend = CAP1, yend = CAP2), arrow=arrow(length=unit(0.2,"cm")),
               alpha = 0.75, color = 'black', linetype = 'dashed')+
  geom_label_repel(data = escores, aes(x = CAP1, y = CAP2, label = rownames(escores)), size = 6) +
  guides(fill = guide_legend(override.aes = list(shape=21)))+
  theme_bw() +
  labs(x = paste("CAP 1 (", perc[1], "%)", sep = ""), y = paste("CAP 2 (", perc[2], "%)", sep = ""), fill = "LUI", shape = "Ecosystem")+
  xlim(-1.5, 1.5)+
  ylim(-1.5, 1.5)+
  coord_fixed(ratio = 1)+
  theme(strip.text.y = element_text(angle = 0), legend.position = "bottom")

ordi_oa+wrap_elements(~plot(spec.part.all,
             Xnames = c("Soil", "Temp", "Mois", "LUI"), # name the partitions
             bg = c("chartreuse3", "darkorange","mediumpurple", "darkgoldenrod2"), alpha = 80, # colour the circles
             digits = 2, # only show 2 digits
             cex = 1.5))

ggsave("./plots/Figure_6.pdf",
       bg = "white",
       width = 700,
       height = 700,  # Wider than tall
       units = "mm",
       fallback_resolution = 1000,  # Helps with raster elements
       device = cairo_pdf,  # Specify EPS format
       dpi = 1000)      # Adjust DPI if needed (default is 300)

## Ordination plot for forest
plot_f<-ordiplot(db_rda_sel_f, scaling = 2) #scaling 2 is used to have angles between variables reflect their correlation

perc_f <- round(100*(summary(db_rda_sel_f)$cont$importance[2, 1:2]), 2)

uscores_f = as.data.frame(plot_f$sites)
vscores_f = as.data.frame(plot_f$species)
escores_f = as.data.frame(plot_f$biplot)
uscores1_f <- dplyr::inner_join(site_data_scaled, rownames_to_column(data.frame(uscores_f), var = "EPS"), by = "EPS")

ordi_f<-ggplot(uscores1_f) + 
  geom_point(aes(x = CAP1, y = CAP2,  fill = lui, 
                 shape = Exploratorium), color = 'black', size = 3) +
  geom_point(data = vscores_f, aes(x = CAP1, y = CAP2), col = "darkgrey", shape = 3)+
  scale_color_manual(values = pal_lui, 
                     aesthetics = c("colour", "fill"),
                     breaks = c("Low", "Medium", "High"),
                     labels = c("Low", "Medium", "High"))+
  scale_shape_manual(values = c(22, 24, 21)) +
  geom_label_repel(data = escores_f, aes(x = CAP1, y = CAP2, label = rownames(escores_f)), size = 6) +
  geom_segment(data = escores_f, aes(x = 0, y = 0, xend = CAP1, yend = CAP2), arrow=arrow(length=unit(0.2,"cm")),
               alpha = 0.75, color = 'black', linetype = 'dashed')+
  guides(fill = guide_legend(override.aes = list(shape=21)))+
  theme_bw() +
  labs(x = paste("CAP 1 (", perc_f[1], "%)", sep = ""), y = paste("CAP 2 (", perc_f[2], "%)", sep = ""))+
  xlim(-1.5, 1.5)+
  ylim(-1.5, 1.5)+ 
  coord_fixed(ratio = 1)+
  theme(strip.text.y = element_text(angle = 0), legend.position = "bottom")


## Ordination plot for grassland
plot_g<-ordiplot(db_rda_sel_g, scaling = 2) #scaling 2 is used to have angles between variables reflect their correlation

perc_g <- round(100*(summary(db_rda_sel_g)$cont$importance[2, 1:2]), 2)

uscores_g = as.data.frame(plot_g$sites)
vscores_g = as.data.frame(plot_g$species)
escores_g = as.data.frame(plot_g$biplot)
uscores1_g <- dplyr::inner_join(site_data_scaled, rownames_to_column(data.frame(uscores_g), var = "EPS"), by = "EPS")

ordi_g<-ggplot(uscores1_g) + 
  geom_point(aes(x = CAP1, y = CAP2,  fill = lui, 
                 shape = Exploratorium), color = 'black', size = 3) +
  geom_point(data = vscores_g, aes(x = CAP1, y = CAP2), col = "darkgrey", shape = 3)+
  scale_color_manual(values = pal_lui, 
                     aesthetics = c("colour", "fill"),
                     breaks = c("Low", "Medium", "High"),
                     labels = c("Low", "Medium", "High"))+
  scale_shape_manual(values = c(21:25)) +
  geom_label_repel(data = escores_g, aes(x = CAP1, y = CAP2, label = rownames(escores_g)), size = 6) +
  geom_segment(data = escores_g, aes(x = 0, y = 0, xend = CAP1, yend = CAP2), arrow=arrow(length=unit(0.2,"cm")),
               alpha = 0.75, color = 'black', linetype = 'dashed')+
  guides(fill = guide_legend(override.aes = list(shape=21)))+
  theme_bw() +
  labs(x = paste("CAP 1 (", perc_g[1], "%)", sep = ""), y = paste("CAP 2 (", perc_g[2], "%)", sep = ""))+
  xlim(-1.5, 1.5)+
  ylim(-1.5, 1.5)+ 
  coord_fixed(ratio = 1)+
  theme(strip.text.y = element_text(angle = 0), legend.position = "bottom")

(ordi_f + wrap_elements(~ plot(spec.part.grs,
                               Xnames = c("Soil", "Temp", "Mois", "LUI"), # name the partitions
                               bg = c("chartreuse3", "darkorange", "mediumpurple", "darkgoldenrod2"), 
                               alpha = 80, # colour the circles
                               digits = 2, # only show 2 digits
                               cex = 1.5))) / 
  (ordi_g + wrap_elements(~ plot(spec.part.for,
                                 Xnames = c("Soil", "Temp", "Mois", "LUI"), # name the partitions
                                 bg = c("chartreuse3", "darkorange", "mediumpurple", "darkgoldenrod2"), 
                                 alpha = 80, # colour the circles
                                 digits = 2, # only show 2 digits
                                 cex = 1.5))) + 
  plot_annotation(tag_levels = 'A')

ggsave("./plots/Figure_S2.pdf",
       bg = "white",
       width = 350,
       height = 350,  # Wider than tall
       units = "mm",
       fallback_resolution = 1000,  # Helps with raster elements
       device = cairo_pdf,  # Specify EPS format
       dpi = 1000)      # Adjust DPI if needed (default is 300)


##=== Identifying indicator species =====

library(indicspecies)

otu_df_Type.lui<-merge(site_data[c(1,6)],spec_abund_rar%>%as.data.frame()%>%rownames_to_column("EPS"),by= "EPS")

indval_Type.lui <- multipatt(otu_df_Type.lui[c(3:10413)], otu_df_Type.lui$Type_lui,
                         duleg = TRUE,
                         control = how(nperm=999)) 

summary(indval_Type.lui, alpha=0.001, indvalcomp=TRUE)

indval_type_lui.sign<-as.data.frame(indval_Type.lui$sign)%>%
  mutate(p.adj = p.adjust(p.value, "BH"))%>%
  filter(p.adj<0.05)%>%
  merge(., tax_mat, by = 0)%>%
  dplyr::select(!Sequence)%>%
  arrange(index, Phylum)

###=== Export table ==== 
indval_type_lui.sign%>%
  mutate_if(is.numeric, round, digits = 3)%>%
  dplyr::filter(stat>0.5)%>%
  mutate(index = factor(index, levels = c("1", "3", "2",  "4", "6", "5"), labels = c("Forest High", "Forest Medium", "Forest Low", "Grassland High", "Grassland Medium", "Grassland Low")))%>%
  separate(index, into = c("Type", "LUI"))%>%
  dplyr::arrange(Type, LUI)%>%
  dplyr::select(Type, LUI, Kingdom:Secondary_lifestyle, stat:p.adj)%>%
  write.csv("./tables/indval_ecosystem_landuse.csv", row.names = FALSE)

###=== Visualization ======
spec_1<-spec_abund_rar%>%
  merge(., 
        site_data[,c(1:6)],
        by.x = 0, by.y = "EPS")%>%
  column_to_rownames("Row.names")%>%
  as.data.frame()%>%
  pivot_longer(cols = c(OTU00001:OTU10411), 
               names_to = "OTU", 
               values_to="abundance")%>%
  merge(., 
        tax_mat[,c(2:11)], 
        by.x = "OTU", 
        by.y = 0,
        all.x = T)%>%
  ungroup()%>%
  group_by(Type, lui, Phylum)%>%
  summarise(abundance = sum(abundance))%>%
  ungroup()%>%
  group_by(Type, lui)%>%
  mutate(abundance_tot = sum(abundance), abundance_perc = 100*(abundance/abundance_tot))%>%
  mutate(Phylum = ifelse(abundance_perc < 0.5, "Others", Phylum))%>%
  mutate(Phylum = factor(Phylum, levels = c("Ascomycota", "Basidiomycota", "Chytridiomycota", "Mortierellomycota", "Mucoromycota", "Glomeromycota", "Rozellomycota", "unclassified", "Others"),
                         labels = c("Ascomycota", "Basidiomycota", "Chytridiomycota", "Mortierellomycota", "Mucoromycota", "Glomeromycota", "Rozellomycota",  "Unclassified", "Others")
  ))%>%
  #mutate(Phylum = fct_reorder(as.factor(Phylum), abundance_perc, .desc = TRUE))%>%
  ggbarplot(x = "lui", y = "abundance_perc", fill = "Phylum", color = "white", facet.by = "Type", ggtheme = ggplot2::theme_minimal())+
  scale_fill_manual(name = "Phylum", values = pal_phyla, aesthetics = "fill")+
  labs(y = "Fungal relative abundance (%)", x = "Land Use Intensity (LUI) Level", fill = "Phylum")


spec_2<-spec_abund_rar%>%
  merge(., 
        site_data[,c(1:6)],
        by.x = 0, by.y = "EPS")%>%
  column_to_rownames("Row.names")%>%
  as.data.frame()%>%
  pivot_longer(cols = c(OTU00001:OTU10411), 
               names_to = "OTU", 
               values_to="abundance")%>%
  merge(., 
        tax_mat[,c(2:11)], 
        by.x = "OTU", 
        by.y = 0,
        all.x = T)%>%
  ungroup()%>%
  group_by(Type, lui, Guild)%>%
  summarise(abundance = sum(abundance))%>%
  ungroup()%>%
  group_by(Type, lui)%>%
  mutate(abundance_tot = sum(abundance), abundance_perc = 100*(abundance/abundance_tot))%>%
  mutate(Guild = ifelse(abundance_perc < 0.15, "Others", Guild))%>%
  mutate(Guild = factor(Guild, levels = c("pathotroph", "saprotroph", "symbiotroph", "unknown"), labels = c("Pathotroph", "Saprotroph", "Symbiotroph", "Unknown")))%>%
  #mutate(Guild = fct_reorder(as.factor(Guild), abundance_perc, .desc = TRUE))%>%
  ggbarplot(x = "lui", y = "abundance_perc", fill = "Guild", color = "white", facet.by = "Type", ggtheme = ggplot2::theme_minimal())+
  scale_fill_manual(name = "Trophic Guild", values = pal_guild, aesthetics = "fill")+
  labs(y = "Fungal relative abundance (%)", x = "Land Use Intensity (LUI) Level", fill = "Trophic Guild")


ind_1<-indval_type_lui.sign %>%
  filter(index %in% c(1, 2, 3, 4, 5, 6)) %>%
  group_by(index) %>%
  dplyr::count(Phylum) %>%
  mutate(abundance_perc = n / sum(n) * 100) %>%
  mutate(index = factor(index, levels = c(2,
                                          3,
                                          1, 
                                          5,
                                          6,
                                          4),
                        labels = c("Forest Low", 
                                   "Forest Medium",
                                   "Forest High",
                                   "Grassland Low",
                                   "Grassland Medium",
                                   "Grassland High")
  )
  )%>%
  mutate(Phylum = ifelse(abundance_perc < 1, "Others", Phylum))%>%
  mutate(Phylum = factor(Phylum, levels = c("Ascomycota", "Basidiomycota", "Chytridiomycota", "Mortierellomycota", "Mucoromycota", "Glomeromycota", "Rozellomycota",  "unclassified", "Others"),
                         labels = c("Ascomycota", "Basidiomycota", "Chytridiomycota", "Mortierellomycota", "Mucoromycota", "Glomeromycota", "Rozellomycota",  "Unclassified", "Others")
  ))%>%
  separate(index, into=c("Type", "lui"))%>%
  mutate(lui = factor(lui, levels = c("Low", "Medium", "High")))%>%
  ggbarplot(x = "lui", 
            y = "abundance_perc", 
            fill = "Phylum", 
            color = "white", 
            label = FALSE,
            facet.by = "Type",
            ggtheme = ggplot2::theme_minimal())+
  scale_fill_manual(name = "Phylum", values = pal_phyla, aesthetics = "fill")+
  labs(y = "Indicator relative abundance (%)", x = "Land Use Intensity (LUI) Level", fill = "Phylum")

ind_2<-indval_type_lui.sign %>%
  filter(index %in% c(1, 2, 3, 4, 5, 6)) %>%
  group_by(index) %>%
  dplyr::count(Guild) %>%
  mutate(abundance_perc = n / sum(n) * 100) %>%
  ungroup()%>%
  mutate(Guild = factor(Guild, levels = c("pathotroph", "saprotroph", "symbiotroph", "unknown"), labels = c("Pathotroph", "Saprotroph", "Symbiotroph", "Unknown")))%>%
  mutate(index = factor(index, levels = c(2,
                                          3,
                                          1, 
                                          5,
                                          6,
                                          4),
                        labels = c("Forest Low", 
                                   "Forest Medium",
                                   "Forest High",
                                   "Grassland Low",
                                   "Grassland Medium",
                                   "Grassland High")
  )
  )%>%
  separate(index, into=c("Type", "lui"))%>%
  mutate(lui = factor(lui, levels = c("Low", "Medium", "High")))%>%
  ggbarplot(x = "lui", 
            y = "abundance_perc", 
            fill = "Guild", 
            color = "white", 
            label = FALSE,
            facet.by = "Type",
            ggtheme = ggplot2::theme_minimal())+
  scale_fill_manual(name = "Trophic Guild", values = pal_guild, aesthetics = "fill")+
  labs(y = "Indicator relative abundance (%)", x = "Land Use Intensity (LUI) Level", fill = "Trophic Guild")

#ggsave("./plots/indicator_species_guild.png", height = 36, units = "cm")

(spec_1+ind_1)/(spec_2+ind_2)+
  plot_layout(guides = 'collect', axes = 'collect_x', axis_titles = 'collect_x')+
  plot_annotation(tag_levels = 'A')

ggsave("./plots/Figure_3.eps",
       bg = "white",
       width = 240, 
       height = 240,  # Wider than tall
       units = "mm",
       device = cairo_ps,  # Specify EPS format
       dpi = 1000)      # Adjust DPI if needed (default is 300)
