#rm(list = ls())
####_####
####___Settings __________________####
####_####
#install.packages("remotes")
#install.packages("readr")
#install.packages("tidyr")
#install.packages("vegan")
#install.packages("reshape")
#install.packages("reshape2")
#install.packages("plyr")
#install.packages("ggplot2")
#install.packages("stringr")
#install.packages("Hmisc")
#install.packages("ggpubr")
#install.packages("fitdistrplus")
#install.packages("RColorBrewer")
#install.packages("dplyr")
#install.packages("car")
#install.packages("MuMIn")
#install.packages("MASS")
#install.packages("pscl")
#install.packages("visreg")
#install.packages ("brms")
#install.packages("loo")
#install.packages("nlme")
#install.packages("mgcv")
#install.packages("tidymv")
#install.packages("tidybayes")
#install.packages("splines")
#install.packages("DHARMa")
require(DHARMa, quietly = TRUE) ## may be missing ... library(broom)
#install.packages("forecast")
#install.packages("nnest")
#install.packages("strucchange")
# remove.packages(c("StanHeaders", "rstan", "TMB"))
#install.packages('TMB')
#install.packages("rstan", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
#install.packages("Matrix", type = 'source')
#install.packages("bcp")
#install.packages("scales")
library(rstan)
library(Matrix)
library(TMB)
library(devtools)
library(scales)
#install.packages("devtools")
#devtools::session_info()
#install.packages("nnet")
#install.packages("patchwork")
#install.packages("knitr")
#install.packages("kableExtra")



#### Library #####
library(patchwork)
library(remotes)
library(readr)
library(tidyr)
library(vegan)
library(reshape)
library(reshape2)
library(plyr)
library(ggplot2)
library(stringr)
library(Hmisc)
library(ggpubr)
library(fitdistrplus)
library(RColorBrewer)
library(dplyr)
library(car)
library(MuMIn)
library(MASS)
library(pscl)
library(visreg)
library (brms)
library(loo)
library(nlme)
library(mgcv)
library(tidymv)
library(tidybayes)
library(splines)
library(DHARMa)
require(DHARMa, quietly = TRUE) ## may be missing ... library(broom)
library(forecast)
library(strucchange)
library(bcp)
library(nnet)
library(knitr)
library(kableExtra)
library(tibble)



#### Remove Scientific notation ####

options(scipen = 999)


#### Replace accent function ####
rmaccent <- function(data){
  data %>%
    mutate(across(.cols = everything(),
                  .fns = ~ stringi::stri_trans_general(., id = "Latin-ASCII")))
}

# - - - - - - - - - - - - - - - - 
# - - - - - - - - - - - - - - - - 
# - - - - - - - - - - - - - - - - 
# - - - - - - - - - - - - - - - - 
####_####
####___Functions__________________ ####
####_####

# Function to generate zeros

##### Data Manipulation functions ####

# Create zeros and fill NAs with equivalent data
CreateZeros <- function(data, location_data, taxonomic_data, control_column, variable_to_zero, value_to_zero){
  
  # Check number of zeros to be added 
  sppresent <- data %>% group_by(!!sym(control_column)) %>% dplyr::summarise(!!variable_to_zero := length(unique(!!sym(variable_to_zero)))) #number of species per transect
  spmissing <- length(unique(data[[variable_to_zero]]))-sppresent[[variable_to_zero]]  #calculate number of species missing in each transect  
  length(unique(data[[control_column]]))
  
  # Create zeros sheet
  makezero <- subset(data, select=c(control_column, variable_to_zero, value_to_zero)) #subset database: necessary variables
  makezeros <- pivot_wider(makezero, names_from=variable_to_zero, values_from=value_to_zero, values_fn=sum) # Transform each species to a column
  makezeros[, !(names(makezeros) %in% c(control_column))] <- makezeros[, !(names(makezeros) %in% c(control_column))] %>% 
    replace( is.na(.), 0) #replace abundance NAs with zeros
  withzeros <- pivot_longer(makezeros, cols = -c(control_column), names_to = variable_to_zero, values_to = value_to_zero)
  length(unique(withzeros[[control_column]]))
  dim(withzeros) 
  ZerosSheet <- withzeros %>% dplyr::filter(!!sym(value_to_zero) == 0)
  dim(ZerosSheet) 
  
  #Create row id column to deal with future duplicates
  ZerosSheet1 <- ZerosSheet %>% 
    mutate(row_id=row_number()) %>%
    relocate(row_id) 
  
  data2 <- data %>% 
    mutate(row_id=row_number()) %>%
    relocate(row_id) 
  
  #Fill ZerosSheet with taxonomic information
  ZerosSheet1 <- merge(ZerosSheet1, taxonomic_data)
  ZerosSheet1 <- ZerosSheet1[!duplicated(ZerosSheet1$row_id), ] #remove duplicates
  nrow(ZerosSheet) - nrow(ZerosSheet1) #must be zero
  
  #Fill ZerosSheet with location information
  ZerosSheet1 <- merge(ZerosSheet1, location_data)
  ZerosSheet1 <- ZerosSheet1[!duplicated(ZerosSheet1$row_id), ] #remove duplicates
  nrow(ZerosSheet) - (nrow(ZerosSheet1)) #must be zero
  
  # Guarantee all columns in ZeroSheet are equal to those in data
  ZerosSheet1 <- suppressMessages(left_join(ZerosSheet1, data2))
  ZerosSheet1 <- ZerosSheet1[!duplicated(ZerosSheet1$row_id), ] #remove duplicates
  nrow(ZerosSheet) - (nrow(ZerosSheet1)) #must be zero
  
  # Add zeros to the original non-zero data frame
  data2 <- rbind(ZerosSheet1, data2)
  
  # Check if the created dataframe has the expected dimensions and number of transects
  if ((sum(spmissing)+nrow(data)) - nrow(data2) == 0 & length(unique(data2[[control_column]])) == length(unique(data[[control_column]]))) {
    cat("Zeros were successfully added!") 
    
  } else {
    cat("Something is wrong! The final dataframe doesn`t show the expected dimensions")
    
  }
  
  # Remove row_id column created to deal with duplicates
  census_all <- data2[, !names(data2) %in% "row_id"]
  
  # return the desired dataframe
  census_all
  
}

##### Stats functions ####

# Function generate standard error 
StandardError <- function(data){
  sd(data, na.rm=TRUE)/sqrt(sum(!is.na(data)))
}

# - Function generates sd, quantile, mean and sum
MeanSdQuantile <- function(data, columns = c("island" ,"year", "mpa"), variables = c("density")){
  data <- data %>% group_by(across(all_of(columns))) %>% 
    dplyr::summarize(across(variables,
                            list(mean = ~mean(.,na.rm=TRUE), 
                                 sum = ~sum(.,na.rm=TRUE),
                                 sd = ~sd(.,na.rm=TRUE), 
                                 quantile25 = ~quantile(., probs=0.25, na.rm=TRUE), 
                                 quantile95 = ~quantile(., probs=0.95, na.rm=TRUE),
                                 se = ~StandardError(.)),
                            .names="{.col}_{.fn}"),                                                                          
                     .groups = "drop")
}

##### Graph Functions ####

root_trans <- function() {
  trans_new(
    name = 'root_trans',
    transform = function(x) x^(1/10),
    inverse = function(x) x^10,
    breaks = function(x) {  # Define a custom function to return desired breaks
      c(180, 420)  # Custom break values
    },
    format = scales::label_number()  # Use standard number formatting
  )
}

##### Resampling functions ####


# Function: minimal sampling area by zone
choosing_transec <- function (DATA){
  # Convert the character vector into a list of symbols
  
  transect_table <- DATA %>% 
    group_by(Id_locality) %>% 
    dplyr:: summarise(ntransect = length(unique(transect_id)),
                      area = mean(area))  %>%
    mutate(area = ntransect * area) 
  
  minSampleArea <- min(transect_table$area)
  
  resu<- unlist (lapply (split(DATA,DATA$Id_locality), function(x){
    x   <- droplevels(x)
    y   <- split (x,x$Id_locality)
    tab <- rep(1:length(y),1000)
    transec <- numeric();  Area=0; i=1
    while (Area<minSampleArea){ ## Define the limit based on smallest location area
      pool <- y[[ tab[i] ]][,"transect_id"]; pool
      pool <- pool[!pool%in%transec];pool
      if (length(pool)==0){
        i <- i+1
      } else{   
        transec[i] <- sample(pool,1,replace=F)
        rm(pool)
        i <- i+1
        Area <- sum (x$area [x$transect_id%in%transec]);
      }
    }
    transec
  }))
  resu=resu[!is.na(resu)]
  names(resu)=NULL
  return(resu)
}


# - Bootstrap Loop Function 
# Resample transect with replace based on minimal sample area (in this case 120)
Bootstrap_ys <- function(data, transectsId, N){
  #create Id_locality column
  data$Id_locality <- paste(data$location, data$sampling_season_year, data$mpa, sep="_")
  ## create output object
  listOfTables <- list()
  ### bootstrap loop
  for (m in 1:N) {
    chosen_transec <- choosing_transec(transectsId)
    sub2 <- droplevels (transectsId[transectsId$transect_id %in%
                                      chosen_transec, ])
    # tapply(sub2$area, list (sub2$sampling_season_year, sub2$Id_locality), sum) # check area selected
    database_red <- droplevels (data [data$transect_id %in% 
                                        sub2$transect_id, ])
    ##select epi_mar
    grouperfilter <- database_red %>%
      filter( species_name %in% 
                "epinephelus_marginatus" )
    length(unique(grouperfilter$transect_id)) # check census number
    
    #select target variables
    grouper2 <-  grouperfilter[,c("Id_locality", "biomass", "density", "transect_id", "total_length")]
    #do not remove total NAs! They account for mean, replace with 0 instead. (Except for total length where NA is necessary)
    grouper2$biomass[is.na(grouper2$biomass)] <- 0
    
    # total per transect (except total_length that is the mean)
    grouper2 <- grouper2 %>% group_by(Id_locality, transect_id) %>% summarise(density = sum(density),
                                                                              biomass = sum(biomass),
                                                                              total_length = mean(total_length, na.rm=TRUE))
    #calculate mean statistics
    grouper2 <- ddply(grouper2,.(Id_locality), summarise, 
                      density = mean(density),
                      biomass = mean(biomass), 
                      total_length = mean(total_length, na.rm=TRUE)
    )
    #Add lap tracker
    grouper2$lap <- m
    #Output table
    listOfTables[[m]] <- separate(grouper2, Id_locality, c("island", "year", "mpa"), sep = "_")
    #Boot loop control 
    print(m)
    
  }
  #True output
  listOfTables
}

# - Bootstrap Loop Function Zone and Period
# Resample transect with replace based on minimal sample area (in this case 120)
Bootstrap_pz_class <- function(data, transectsId, N){
  #create Id_locality column
  data$Id_locality <- paste(data$yearcat, data$mpa, sep="_")
  ## create output object
  listOfTables <- list()
  ### bootstrap loop
  for (m in 1:N) {
    chosen_transec <- choosing_transec(transectsId)
    sub2 <- droplevels (transectsId[transectsId$transect_id %in%
                                      chosen_transec, ])
    # tapply(sub2$area, list (sub2$sampling_season_year, sub2$Id_locality), sum) # check area selected
    database_red <- droplevels (data [data$transect_id %in% 
                                        sub2$transect_id, ])
    ##select epi_mar
    grouperfilter <- database_red %>%
      filter( species_name %in% 
                "epinephelus_marginatus" )
    length(unique(grouperfilter$transect_id)) # check census number
    
    #select target variables
    grouper2 <-  grouperfilter[,c("Id_locality", "biomass", "density", "transect_id", "total_length","class")]
    #do not remove total NAs! They account for mean, replace with 0 instead. (Except for total length where NA is necessary)
    grouper2$biomass[is.na(grouper2$biomass)] <- 0
    
    # total per transect (except total_length that is the mean)
    grouper2 <- grouper2 %>% group_by(Id_locality, transect_id, class) %>% summarise(density = sum(density),
                                                                                     biomass = sum(biomass),
                                                                                     total_length = mean(total_length, na.rm=TRUE))
    #calculate mean statistics
    grouper2 <- ddply(grouper2,.(Id_locality, class), summarise, 
                      density = mean(density),
                      biomass = mean(biomass), 
                      total_length = mean(total_length, na.rm=TRUE), 
                      biomass_sum = sum(biomass) )
    #Add lap tracker
    grouper2$lap <- m
    #Output table
    listOfTables[[m]] <- separate(grouper2, Id_locality, c("yearcat", "mpa"), sep = "_")
    #Boot loop control 
    print(m)
    
  }
  #True output
  listOfTables
}

# - - - - - - - - - - - - - - - - 
# - - - - - - - - - - - - - - - - 
# - - - - - - - - - - - - - - - - 
# - - - - - - - - - - - - - - - - 
####_####
####___Data__________________ ####
####_####
##### Import data sets  ####

#combined sheet
census <- read_csv("censussheets/census_all.csv")
dim(census)
length(unique(census$species_code))
length(unique(census$species_name))
length(unique(census$transect_id))

#location info
location_data <- read_csv("censussheets/Darwincorefiles/location_information_combined.csv")
dim(location_data)
length(unique(location_data$transect_id))
length(unique(location_data$location))

#taxonomic info
taxonomic_data <- read_csv("censussheets/Darwincorefiles/taxonomic_information_combined.csv")
dim(taxonomic_data)
length(unique(taxonomic_data$species_name))
length(unique(taxonomic_data$species_code))

##### Create zeros #####
census_all <- CreateZeros(data = census, location_data = location_data, taxonomic_data = taxonomic_data, control_column ="transect_id", variable_to_zero ="species_name", value_to_zero = "abundance")
length(unique(census_all$species_name))
length(unique(census_all$transect_id))

# Check data
rm(location_data, taxonomic_data)

##### Create useful columns #####

censusSheet <- census_all

#Create year column by season
censusSheet$sampling_season_year <- censusSheet$sampling_season 
censusSheet$sampling_season_year <- gsub("summer_","",as.character(censusSheet$sampling_season_year))

#area column
censusSheet$area <- 40

#control column
censusSheet$control_column <- 1

##### Select data #####

#Remove winter samples
censusSheet <- censusSheet[censusSheet$sampling_season != "winter_2012", ] 
dim(censusSheet) #28474 lines 36 columns 


#Select target islands
censusSheet <- censusSheet %>%
  filter(location %in% c("arvoredo_island","deserta_island","gale_island","xavier_island")) #24843 lines 37 columns 

#Rename islands names
censusSheet$location <- recode_factor(censusSheet$location, arvoredo_island = "arvoredo",
                                      deserta_island = "deserta",
                                      gale_island = "gale",
                                      xavier_island = "xavier")

#Rename mpa column
censusSheet$mpa <- replace(censusSheet$no_take_zone, censusSheet$no_take_zone=="yes", "protected")
censusSheet$mpa <- replace(censusSheet$mpa, censusSheet$mpa=="no", "unprotected")

#Id_locality column
censusSheet$Id_locality <- paste(censusSheet$location, censusSheet$sampling_season_year, censusSheet$mpa, sep="_")

##### Add variables #####

#Calculate selected rows biomass (grams/m2)
class(censusSheet$abundance)
censusSheet$abundance <- as.numeric(censusSheet$abundance)
class(censusSheet$total_length)
censusSheet$total_length <- as.numeric(censusSheet$total_length)
censusSheet$weigth2  <- censusSheet$abundance*(censusSheet$a*(censusSheet$total_length^censusSheet$b))
censusSheet$weigth2[is.na(censusSheet$weigth2)] <- 0
censusSheet$biomass <- censusSheet$weigth2/censusSheet$area ### Biomass 

#Turn biomass Nans into 0
censusSheet <- censusSheet %>% 
  mutate(biomass = ifelse(is.na(biomass), 0, biomass))

#Check if biomass is transformed and if total_length NAs remain present
summary(is.na(censusSheet$biomass))
summary(is.na(censusSheet$total_length))

#Calculate selected rows density (abundance/m2)
censusSheet$density <- censusSheet$abundance/censusSheet$area 

#Create year categories column
censusSheet$yearcat <- cut(as.numeric(censusSheet$sampling_season_year), breaks=c(2007, 2010, 2013, 2016, 2019, 2023),
                           labels=c('2008-2010', '2011-2013', '2014-2016', '2017-2019','2020-2023'))
sum(is.na(censusSheet$yearcat)) 

# Check data
length(unique(censusSheet$species_name)) #n species
length(unique(censusSheet$species_code)) #n codes
length(unique(censusSheet$transect_id)) #n transect_id

##### Final Sheets #####


# - - - - - - - - - - - - - - - - 

#Add grouper size classes
censusSheet <- censusSheet %>%
  mutate(class=case_when(
    species_name == "epinephelus_marginatus" & total_length <= 10 ~ "Recruit",
    species_name == "epinephelus_marginatus" & total_length > 10 & total_length <= 30 ~ "Juvenile_1",
    species_name == "epinephelus_marginatus" & total_length > 30 & total_length <= 46 ~ "Juvenile_2",
    species_name == "epinephelus_marginatus" & total_length > 46 & total_length <= 73 ~ "Adult",
    species_name == "epinephelus_marginatus" & total_length > 73 ~ "Reproductive_Matrix"
  ))

###### • Dataframe - all species ####
head(censusSheet) 

#Filter grouper data only
grouperSheet <- censusSheet %>%
  filter( species_name %in% "epinephelus_marginatus" )

##### • Dataframe - dusky grouper ####
head(grouperSheet) 

##### Explore dusky groper df ####

x11() # use windows() for windows system
hist(grouperSheet$total_length)
hist(grouperSheet$abundance) 
hist(grouperSheet$density) 
hist(grouperSheet$biomass) 
hist(log(grouperSheet$biomass+1)) 
# Data is zero inflated

# - - - - - - - - - - - - - - - - 
# - - - - - - - - - - - - - - - - 
# - - - - - - - - - - - - - - - - 
# - - - - - - - - - - - - - - - - 
####_####
####___ Bootstraps_____________ ####
####_####
# Analysis of data using the time series as a continuous.

##### 1. Hurdle-gamma #####
# Bootstrap using all years and sites (mpa+island [protenction,location]).

###### Sampling effort table ####
transect_table_ys <- censusSheet %>% 
  group_by(sampling_season_year, mpa, location, area) %>% 
  dplyr:: summarise(ntransect = length(unique(transect_id)))  %>%
  mutate(area = ntransect * area)

# Minimal, maximal sample area: 
hist(transect_table_ys$ntransect)
min(transect_table_ys$area)
max(transect_table_ys$area)
# Data has uneven sampling efforts

###### Resampling ####

# Create unique transect data with number of registers per transect
transectUnique_ys <- data.frame(cast(censusSheet, formula=transect_id+sampling_season_year+mpa+location+area ~.,
                                     value = "control_column", fun.aggregate = sum, na.rm=TRUE)) 
# Remove number of registers per transect
transectUnique_ys <- droplevels (transectUnique_ys[,-ncol(transectUnique_ys)])
transectUnique_ys$Id_locality <- paste(transectUnique_ys$sampling_season_year, 
                                       transectUnique_ys$mpa, transectUnique_ys$location, sep="_")

#Bootstrap with replacement   
grouperdataboot_ys <- Bootstrap_ys(data=censusSheet, transectsId = transectUnique_ys, N=1000)


#Binding Boot output
grouperdataboot_ys <- bind_rows(grouperdataboot_ys)
dim(grouperdataboot_ys)
head(grouperdataboot_ys)
###### Statistics #####

# Boot histograms

hist(grouperdataboot_ys$density)
hist(grouperdataboot_ys$biomass)
hist(grouperdataboot_ys$total_length)

# Calculate statistics

bootmean_ys <- MeanSdQuantile(grouperdataboot_ys, columns = c("year", "mpa", "island"), variables = c("total_length", "density", "biomass")) %>% 
  mutate_if(is.numeric, ~replace_na(., 0))

hist(bootmean_ys$density_mean)
hist(bootmean_ys$biomass_mean)
hist(bootmean_ys$total_length_mean)

###### Tables Stats #####

#Create data frame to generate table and calculate mean statistics per period/site
bootmean_ys_selected <- subset(bootmean_ys, select = c("island", "mpa", "year", "density_mean", "density_sd", "biomass_mean", "biomass_sd"))

#Rename columns
names(bootmean_ys_selected)[names(bootmean_ys_selected) == "mpa"] <- "protection"
names(bootmean_ys_selected)[names(bootmean_ys_selected) == "density_mean"] <- "density mean"
names(bootmean_ys_selected)[names(bootmean_ys_selected) == "density_sd"] <- "density sd"
names(bootmean_ys_selected)[names(bootmean_ys_selected) == "biomass_mean"] <- "biomass mean"
names(bootmean_ys_selected)[names(bootmean_ys_selected) == "biomass_sd"] <- "biomass sd"
names(bootmean_ys_selected)[names(bootmean_ys_selected) == "site"] <- "site"

#Rename site names
bootmean_ys_selected$island <- recode_factor(bootmean_ys_selected$island, arvoredo_protected = "ArvoredoMPA",
                                             deserta_protected = "Deserta",
                                             gale_protected = "Galé",
                                             arvoredo_unprotected = "ArvoredoNPZ",
                                             xavier_unprotected = "Xavier")

#Create table

table <- kable(bootmean_ys_selected, format = "latex", booktabs = TRUE, digits = 3) %>%
  kable_stysing(latex_options = c("striped", "scale_down"))

writeLines(as.character(table), "my_table.tex")

# Convert the table to a character string
latex_table_string <- as.character(table)

# macOS version
write.table(latex_table_string, pipe("pbcopy"), sep = "\n", row.names = FALSE, col.names = FALSE, quote = FALSE)


# Calculate mean per year 2019 forward
bootmean_ys_selected_2019 <- subset(bootmean_ys,  year >= 2019)
bootmean_ys_selected_2019 %>% group_by(year) %>% summarise(year_density = sum(density_mean),
                                                           year_biomass = sum(biomass_mean)) %>% summarise(mean_year_density = mean(year_density),
                                                                                                           mean_year_biomass = mean(year_biomass))
# Calculate mean per year 2019 backwards
bootmean_ys_selected_2008 <- subset(bootmean_ys,  year <= 2019)
bootmean_ys_selected_2008 %>% group_by(year) %>% summarise(year_density = sum(density_mean),
                                                           year_biomass = sum(biomass_mean)) %>% summarise(mean_year_density = mean(year_density),
                                                                                                           mean_year_biomass = mean(year_biomass))

##### 2. Life Stage Histogram #####
# Bootstrap using time periods and protection (mpa), with life stage structuring.

# Create unique transect data with number of registers per transect
transectUnique_pz <- data.frame(cast(censusSheet, formula=transect_id+yearcat+mpa+area ~.,
                                     value = "control_column", fun.aggregate = sum, na.rm=TRUE)) 
# Remove number of registers per transect
transectUnique_pz <- droplevels(transectUnique_pz[,-ncol(transectUnique_pz)])
transectUnique_pz$Id_locality <- paste(transectUnique_pz$yearcat, 
                                       transectUnique_pz$mpa, sep="_")

###### Resampling ####

#Bootstrap with replacement   
Bootstrap_pz_class <- Bootstrap_pz_class(data=censusSheet, transectsId = transectUnique_pz, N=1000)
View(Bootstrap_pz_class)
#Binding Boot output
Bootstrap_pz_class <- bind_rows(Bootstrap_pz_class)
dim(Bootstrap_pz_class)


###### Statistics #####

# Boot histograms

hist(Bootstrap_pz_class$density)
hist(Bootstrap_pz_class$biomass)
hist(Bootstrap_pz_class$total_length)

# Calculate statistics

bootmean_pz_class_o <- MeanSdQuantile(Bootstrap_pz_class, columns = c("yearcat", "mpa", "class"), variables = c("total_length", "density", "biomass", "biomass_sum")) %>% 
  mutate_if(is.numeric, ~replace_na(., 0))

hist(bootmean_pz_class_o$density_mean)
hist(bootmean_pz_class_o$biomass_mean)
hist(bootmean_pz_class_o$total_length_mean)

# - - - - - - - - - - - - - - - - 
# - - - - - - - - - - - - - - - - 
# - - - - - - - - - - - - - - - - 
# - - - - - - - - - - - - - - - - 
####_####
####___ Hurdle-models___________ ####
####_####

#####  Data inspection ####
# Data distribution
dim(bootmean_ys)

# log biomass
bootmean_ys$logbiomass <- log(bootmean_ys$biomass_mean+1)

hist(bootmean_ys$density_mean) #zero-inflated
hist(bootmean_ys$biomass_mean) #zero-inflated
hist(bootmean_ys$logbiomass) #zero-inflated

#####  Breakpoint inspection ####
# Perform breakpoint analysis biomass
bp_lm_b <- breakpoints(biomass_mean ~ year, data = bootmean_ys)
summary(bp_lm_b)
plot(bp_lm_b)
# Extract the breakpoint indices
break_indices_b <- bp_lm_b$breakpoints

# Look up the corresponding years in your dataset
bootmean_ys$year[break_indices_b] #Sugest 3 breakpoints (2013 2016 2018)


# Perform breakpoint analysis density
bp_lm_d <- breakpoints(density_mean ~ year, data = bootmean_ys)
summary(bp_lm_d)
plot(bp_lm_d)
# Extract the breakpoint indices
break_indices_d <- bp_lm_d$breakpoints

# Look up the corresponding years in your dataset
bootmean_ys$year[break_indices_d] #Sugest 3 breakpoints (2013 2016 2018)


# Isolate variables acording to zone
protected_bootmean_ys <- filter(bootmean_ys, mpa == "protected")
unprotected_bootmean_ys <- filter(bootmean_ys, mpa == "unprotected")

# Breakpoints
bcp_biomass_l <- bcp(bootmean_ys$biomass_mean)
bcp_biomass_p_l <- bcp(protected_bootmean_ys$biomass_mean)
bcp_biomass_u_l <- bcp(unprotected_bootmean_ys$biomass_mean)
bcp_density_l <- bcp(bootmean_ys$density_mean)
bcp_density_p_l <- bcp(protected_bootmean_ys$density_mean)
bcp_density_u_l <- bcp(unprotected_bootmean_ys$density_mean)
# Plot the posterior probabilities of change points
plot(bcp_biomass_l)
plot(bcp_biomass_p_l)
plot(bcp_biomass_u_l)
plot(bcp_density_u_l)
plot(bcp_density_p_l)
plot(bcp_density_u_l)

# Create datasheet with density and biomass posterior probabilities
bootmean_ys_breaks <- data.frame(
  bootmean_ys,
  breaks_biomass = bcp_biomass_l$posterior.prob,
  breaks_density = bcp_density_l$posterior.prob
)

# Check most significant years
# Posterior probabilities
post_prob_b_l <- bcp_biomass_l$posterior.prob
post_prob_d_l <- bcp_density_l$posterior.prob
# Set a threshold for determining significant change points
threshold <- 0.9  
# Identify the indices where the posterior probability exceeds the threshold
significant_bio_indices_l <- which(post_prob_b_l > threshold)
significant_den_indices_l <- which(post_prob_d_l > threshold)

# Find the corresponding years for these significant change points
bootmean_ys$year[significant_bio_indices_l]
bootmean_ys$year[significant_den_indices_l]

#####  Raw graphs ####

###### • Bubbleplot #####

###Summary table#
summarytable <- bootmean_ys[, c("island", "year", "mpa", "density_mean", "density_sd", "biomass_mean", "biomass_sd")]
clipr::write_clip(summarytable) 

#Create site column
summarytable$site <- paste(summarytable$island, summarytable$mpa, sep = "_")
bootmean_ys$site <- paste(bootmean_ys$island, bootmean_ys$mpa, sep = "_")

#Create binary column
summarytable$binary <- (summarytable$density_mean>0)+0

#Summarize mean abundance
GrouperSampleYearsBoot <- bootmean_ys %>% 
  group_by(site, island, year, mpa) %>% 
  dplyr::summarize(density=mean(density_mean),
                   biomass=mean(biomass_mean),
  )

#Create binary column
GrouperSampleYearsBoot$binary <- (GrouperSampleYearsBoot$density>0)+0

#Mean values and SD per period
grp2008<- filter(GrouperSampleYearsBoot, year < 2019) 
mean(grp2008$density)
sd(grp2008$density)
mean(grp2008$biomass)
sd(grp2008$biomass)


#Mean values and SD per period
grp2023<- filter(GrouperSampleYearsBoot, year > 2018) 
mean(grp2023$density)
sd(grp2023$density)
mean(grp2023$biomass)
sd(grp2023$biomass)

# T-test
t.test(grp2008$density,grp2023$density)
t.test(grp2008$binary,grp2023$binary)

#Percentages
#Biomass
(mean(grp2023$biomass)/mean(grp2008$biomass)*100)-100
#density
(mean(grp2023$density)/mean(grp2008$density)*100)-100

#Create binary column
summarytable$binary <- (summarytable$density_mean>0)+0

# Bubble plot 
bubble_plot <- ggplot(summarytable, aes(x=year, y=factor(site, levels = c("xavier_unprotected",
                                                                          "arvoredo_unprotected",
                                                                          "deserta_protected",
                                                                          "arvoredo_protected",
                                                                          "gale_protected")), color= factor(binary), fill = factor(binary))) +
  geom_point(size = ((summarytable$density_mean)))+
  
  scale_color_manual(values = c("1" = "#9DA29B",
                                "0"="white")) +
  
  geom_rect(fill = "orange", alpha = .005, color = "white",
            aes(ymin = 0, ymax = 2.5,
                xmin = -Inf, xmax = Inf)) +
  geom_rect(fill = "turquoise", alpha = .005, color = "white",
            aes(ymin = 2.5, ymax = 6,
                xmin = -Inf, xmax = Inf)) +
  geom_point(size = (((summarytable$density_mean))+0.02)*170)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 60, hjust=1, size=8)) +
  
  ggtitle(expression(paste("Occurence of ",italic("Epinephelus marginatus")))) +
  labs(y = "Islands", x = "Year") 


# Bubble plot with SD
bubble_plot_sd <- ggplot(summarytable, aes(x=year, y=factor(site, levels = c("xavier_unprotected",
                                                                             "arvoredo_unprotected",
                                                                             "deserta_protected",
                                                                             "arvoredo_protected",
                                                                             "gale_protected")), color= factor(binary), fill = factor(binary))) +
  geom_point(size = ((summarytable$density_mean)))+
  
  geom_rect(fill = "orange", alpha = .005, color = "white",
            aes(ymin = 0, ymax = 2.5,
                xmin = -Inf, xmax = Inf)) +
  geom_rect(fill = "turquoise", alpha = .005, color = "white",
            aes(ymin = 2.5, ymax = 6,
                xmin = -Inf, xmax = Inf)) +
  geom_point(size = (((summarytable$density_mean))+0.02)*160, alpha=1, color="darkgrey", stroke = 0)+
  geom_point(size = (((summarytable$density_mean+summarytable$density_sd))+0.02)*160, alpha=0.13, color="black", stroke = 0)+
  
  # Add white points for density_mean == 3.4
  #geom_point(data = subset(summarytable, density_mean == 0), 
  #aes(size = (density_mean + 0.02) * 160), 
  #shape = 21, color = "white", fill = "white", stroke = 0) +
  
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 60, hjust=1, size=8)) +
  
  ggtitle(expression(paste("Occurence of ",italic("Epinephelus marginatus")))) +
  labs(y = "Islands", x = "Year") 

# Convert year to a factor for discrete x-axis plotting
grouperdataboot_ys$year <- factor(grouperdataboot_ys$year)

# Mean density for each 'year' and 'zone' combination
grouperdataboot_ys_means <- grouperdataboot_ys %>%
  group_by(year, mpa) %>%
  summarize(mean_density = mean(density),
            mean_biomass = mean(biomass), 
            mean_sd_density = sd(density),
            mean_sd_biomass = sd(biomass),
            .groups = 'drop')

##### . Density Plot ####

#Line plot
ggplot(grouperdataboot_ys, aes(x = year, y = density)) +
  geom_jitter(alpha=0.4, size = 0.4, width = 0.1, height = 0, shape = 1, color = "black") + 
  geom_line(data = grouperdataboot_ys_means, aes(x = year, y = mean_density, group = mpa, color = mpa), linewidth=1) +
  facet_wrap(~mpa, scales = "free_y") +
  theme_bw() +
  theme(
    legend.position = "none",  # Remove legend if not needed
    axis.text.x = element_text(angle = 45, hjust = 1)  # Incline the x-axis text
  ) +
  labs(y = "Density (individuals/m²)", x = "Year")  # Add labels

#line plot
bufferlineplot_densraw <- ggplot(grouperdataboot_ys_means, aes(x = year, y = mean_density)) +
  geom_ribbon(aes(ymin = pmax(mean_density - mean_sd_density, 0), ymax = mean_density + mean_sd_density, group = mpa, fill= mpa), alpha = 0.3) +
  geom_line(data = grouperdataboot_ys_means, aes(x = year, y = mean_density, group = mpa, color = mpa), linewidth=1) +
  facet_wrap(~mpa) +
  theme_minimal() +
  theme(
    legend.position = "none",  
    axis.text.x = element_text(angle = 45, hjust = 1)  # Incline the x-axis text
  ) +
  labs(y = "Density (individuals/m²)", x = "Year")  # Add labels


#Boxplot
ggplot(grouperdataboot_ys, aes(x = year, y = density)) +
  geom_boxplot() +  # Add boxplots
  facet_wrap(~mpa) +
  theme_minimal() +
  theme(
    legend.position = "none",  # Remove legend if not needed
    axis.text.x = element_text(angle = 45, hjust = 1)  # Incline the x-axis text
  ) +
  labs(y = "Density (individuals/m²)", x = "Year")  # Add labels


#####. Biomass Plots ####

#line plot
ggplot(grouperdataboot_ys, aes(x = year, y = biomass)) +
  geom_point(alpha = 0.2) + 
  geom_line(data = grouperdataboot_ys_means, aes(x = year, y = mean_biomass, group = mpa, color = mpa), linewidth=1) +
  facet_wrap(~mpa) +
  theme_minimal() +
  theme(
    legend.position = "none",  
    axis.text.x = element_text(angle = 45, hjust = 1)  # Incline the x-axis text
  ) +
  labs(y = "Biomass (g/m²)", x = "Year")  # Add labels

#line plot
ggplot(grouperdataboot_ys, aes(x = year, y = biomass)) +
  geom_point(alpha = 0.2) + 
  geom_smooth(data = grouperdataboot_ys_means, aes(x = year, y = mean_biomass, group = mpa, color = mpa), linewidth=1) +
  facet_wrap(~mpa) +
  theme_minimal() +
  theme(
    legend.position = "none",  
    axis.text.x = element_text(angle = 45, hjust = 1)  # Incline the x-axis text
  ) +
  labs(y = "Biomass (g/m²)", x = "Year")  # Add labels

#line plot
bufferlineplot_bioraw <- ggplot(grouperdataboot_ys_means, aes(x = year, y = mean_biomass)) +
  geom_ribbon(aes(ymin = pmax(mean_biomass - mean_sd_biomass, 0), ymax = mean_biomass + mean_sd_biomass, group = mpa, fill= mpa), alpha = 0.3) +
  geom_line(data = grouperdataboot_ys_means, aes(x = year, y = mean_biomass, group = mpa, color = mpa), linewidth=1) +
  facet_wrap(~mpa) +
  theme_minimal() +
  theme(
    legend.position = "none",  
    axis.text.x = element_text(angle = 45, hjust = 1)  # Incline the x-axis text
  ) +
  labs(y = "Biomass (g/m²)", x = "Year")  # Add labels


#Boxplot
ggplot(grouperdataboot_ys, aes(x = year, y = biomass)) +
  geom_boxplot() +  # Add boxplots
  facet_wrap(~mpa) +
  theme_minimal() +
  theme(
    legend.position = "none",  # Remove legend if not needed
    axis.text.x = element_text(angle = 45, hjust = 1)  # Incline the x-axis text
  ) +
  labs(y = "Biomass (g/m²)", x = "Year")  # Add labels

#line plot zoomed
#0 to 60
ys_bio_plot_normals_raw <- ggplot(grouperdataboot_ys, aes(x = year, y = biomass)) +
  geom_jitter(alpha=0.4, size = 0.4, width = 0.1, height = 0, shape = 1, color = "black") + 
  geom_line(data = grouperdataboot_ys_means, aes(x = year, y = mean_biomass, group = mpa, color = mpa), size = 1.2) +
  facet_wrap(~mpa) +
  coord_cartesian(ysim = c(0, 60)) +  # Set Y-axis limits to zoom in
  theme_bw() +
  theme(
    legend.position = "none",  
    axis.text.x = element_text(angle = 45, hjust = 1),  # Incline the x-axis text
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_blank(),   # Remove panel background
    strip.background = element_blank(),
    strip.text.x = element_blank()
  ) +
  labs(y = "Biomass (g/m²)", x = "Year")  # Add labels

#60 to 420 (xˆ10)
ys_bio_plot_extremes_raw <- ggplot(grouperdataboot_ys, aes(x = year, y = biomass)) +
  geom_jitter(alpha=0.4,size = 0.4, width = 0.1, height = 0, shape = 1, color = "black") + 
  facet_wrap(~mpa) +
  coord_cartesian(ysim = c(60.1, 450)) +
  scale_y_continuous(trans = root_trans(), labels = scales::label_number()) +  # Apply square root transformation  # Adjust the limits for the transformed scale
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_blank(),  # Remove panel background
    axis.text.x.bottom = element_blank(),  # Remove y-axis text labels
    axis.ticks.x.bottom = element_blank(),  # Remove y-axis ticks
    axis.title.x.bottom = element_blank(),
    axis.title.y = element_blank()
  ) 

# Combine the plots
combined_plot_raw <- ys_bio_plot_extremes_raw / ys_bio_plot_normals_raw+ plot_layout(heights = c(0.5, 2))


##### Models ######

###### 1.Density hurdle-gamma #####

# Checking intercept : in this case the reference levels for the intercept are > protected and > 2008
levels(as.factor(bootmean_ys$mpa))
min(bootmean_ys$year)
# setting unprotected as reference level
bootmean_ys$mpa <- factor(bootmean_ys$mpa, levels = c("unprotected", "protected"))
bootmean_ys$year <- as.numeric(as.character(bootmean_ys$year))


# Checking priors
priors <- get_prior(bf(density_mean ~ year + mpa, hu ~ year + mpa), data = bootmean_ys,
                    family = hurdle_gamma(link = "log", link_shape = "log", link_hu = "logit"))

# Fit the model with the specified priors
ys_den_hurdlegam_priors <- brm(bf(density_mean ~ year + mpa, 
                                  hu ~ year + mpa), data = bootmean_ys,
                               family = hurdle_gamma(link = "log", link_shape = "log", link_hu = "logit"),
                               prior = priors,
                               chains = 4, iter = 20000, warmup = 10000, control = list(adapt_delta = 0.95))

summary(ys_den_hurdlegam_priors) # Good model fit: Rhat=1 and Bulk_ESS + Tail_ESS >1000
conditional_effects(ys_den_hurdlegam_priors)
plot(ys_den_hurdlegam_priors)
pp_check(ys_den_hurdlegam_priors, type = "scatter_avg") 
pp_check(ys_den_hurdlegam_priors, type = "dens_overlay")

# Loo test
plot(loo(ys_den_hurdlegam_priors, cores = getOption("mc.cores", 1)))
loo_ys_priors <- loo(ys_den_hurdlegam_priors)

######. Plot predicted values ####
(dens_predict <- bootmean_ys %>%
   tidybayes::add_predicted_draws(ys_den_hurdlegam_priors) %>%  # adding the posterior distribution
   ggplot(aes(x = year, y = density_mean)) +  
   facet_wrap(~mpa) +
   stat_lineribbon(aes(y = .prediction), .width = c(.95, .75, .50),  # regression line and CI
                   alpha = 0.5, colour = "black") +
   geom_point(data = bootmean_ys, colour = "darkseagreen4", size = 2) +   # raw data
   scale_fill_brewer(palette = "Greys") +
   ysab("Predicted density (m2)") +  # latin name for red knot
   xlab("Year") +
   theme_bw() +
   theme(legend.title = element_blank()))

######. Plot predicted values with mean bootstrap values ####
grouperdataboot_ys_means$year <- as.numeric(as.character(grouperdataboot_ys_means$year))
(dens_predict_comb <- bootmean_ys %>%
    tidybayes::add_predicted_draws(ys_den_hurdlegam_priors) %>%  # adding the posterior distribution
    ggplot(aes(x = year, y = density_mean)) +  
    facet_wrap(~mpa) +
    stat_lineribbon(aes(y = .prediction), .width = c(.95, .75, .50),  # regression line and CI
                    alpha = 0.5, colour = "black") +
    geom_line(data = grouperdataboot_ys_means, aes(x = year, y = mean_density, group = mpa, color = mpa), linewidth=1) +   # raw data
    scale_fill_brewer(palette = "Greys") +
    ysab("Predicted density (m2)") +  # latin name for red knot
    xlab("Year") +
    theme_bw() +
    theme(legend.title = element_blank()))

######. Plot predicted values with mean bootstrap values and sd ####
bootmean_ys$year <- as.numeric(as.character(bootmean_ys$year))
grouperdataboot_ys_means$year <- as.numeric(as.character(grouperdataboot_ys_means$year))

graph_density_2023 <- ggplot() +
  # Elements from dens_predict_comb (base layers)
  stat_lineribbon(data = bootmean_ys %>% tidybayes::add_predicted_draws(ys_den_hurdlegam_priors),
                  aes(x = year, y = .prediction, group = mpa),
                  .width = c(.95, .75, .50), alpha = 0.5, colour = "black") +
  
  # Separate geom_ribbon layers for each 'mpa' group
  geom_ribbon(data = subset(grouperdataboot_ys_means, mpa == "protected"), 
              aes(x = year, ymin = pmax(mean_density - mean_sd_density, 0), ymax = mean_density + mean_sd_density),
              fill = "turquoise", alpha = 0.15) +
  geom_ribbon(data = subset(grouperdataboot_ys_means, mpa == "unprotected"), 
              aes(x = year, ymin = pmax(mean_density - mean_sd_density, 0), ymax = mean_density + mean_sd_density),
              fill = "orange", alpha = 0.15) +
  
  # Add dashed lines for the upper boundary, setting custom colors for color based on 'mpa'
  
  geom_line(data = subset(grouperdataboot_ys_means, mpa == "protected"), 
            aes(x = year, y = mean_density + mean_sd_density, group = mpa),
            linetype = "dashed", size = 0.3, color = "turquoise", alpha=0.8) +
  geom_line(data = subset(grouperdataboot_ys_means, mpa == "unprotected"), 
            aes(x = year, y = mean_density + mean_sd_density, group = mpa),
            linetype = "dashed", size = 0.3, color = "orange", alpha=0.8) +
  
  # Add dashed lines for the lower boundary, setting custom colors for color based on 'mpa'
  geom_line(data = subset(grouperdataboot_ys_means, mpa == "protected"), 
            aes(x = year, y = pmax(mean_density - mean_sd_density, 0), group = mpa),
            linetype = "dashed", size = 0.3, color = "turquoise", alpha=0.8) +
  geom_line(data = subset(grouperdataboot_ys_means, mpa == "unprotected"), 
            aes(x = year, y = pmax(mean_density - mean_sd_density, 0), group = mpa),
            linetype = "dashed", size = 0.3, color = "orange", alpha=0.8) +
  
  # Add solid lines for mean density, setting custom colors for color based on 'mpa'
  geom_line(data = subset(grouperdataboot_ys_means, mpa == "protected"), 
            aes(x = year, y = mean_density, group = mpa),
            color = "turquoise", 
            size = 0.8, alpha = 0.9) +
  geom_line(data = subset(grouperdataboot_ys_means, mpa == "unprotected"), 
            aes(x = year, y = mean_density, group = mpa),
            color = "orange", 
            size = 0.8, alpha = 0.9) +
  
  # Common elements
  facet_wrap(~mpa) +
  coord_cartesian(ysim = c(0.10, 0.25)) +
  scale_fill_brewer(palette = "Greys") +
  ysab("Density (individuals/m²)") +
  xlab("Year") +
  theme_bw() +
  scale_x_continuous(breaks = 2008:2023)+
  theme(legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1), 
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        strip.background = element_blank())

# Combine plots
combined_denplot_model <- graph_density_2023_2 / graph_density_2023_1 + plot_layout(heights = c(0.5, 2))


###### 2.Biomass hurdle-gamma #####

# Checking intercept : in this case the reference levels for the intercept are > protected and > 2008
levels(as.factor(bootmean_ys$mpa))
min(bootmean_ys$year)
# setting unprotected as reference level
bootmean_ys$mpa <- factor(bootmean_ys$mpa, levels = c("unprotected", "protected"))
bootmean_ys$year <- as.numeric(as.character(bootmean_ys$year))

# Fit the model with the specified priors
ys_bio_hurdlegam_priors <- brm(bf(biomass_mean ~ year + mpa, 
                                  hu ~ year + mpa), data = bootmean_ys,
                               family = hurdle_gamma(link = "log", link_shape = "log", link_hu = "logit"),
                               chains = 4, iter = 20000, warmup = 10000, control = list(adapt_delta = 0.9))

summary(ys_bio_hurdlegam_priors) # Good model fit: Rhat=1 and Bulk_ESS + Tail_ESS >1000
conditional_effects(ys_bio_hurdlegam_priors)
plot(ys_bio_hurdlegam_priors)
pp_check(ys_bio_hurdlegam_priors, type = "scatter_avg") 
pp_check(ys_bio_hurdlegam_priors, type = "dens_overlay")

# Loo test
plot(loo(ys_bio_hurdlegam_priors, cores = getOption("mc.cores", 1)))
loo_yz_bio_priors <- loo(ys_bio_hurdlegam_priors)


######. Plot predicted values #####
(model_fit <- bootmean_ys %>%
   tidybayes::add_predicted_draws(ys_bio_hurdlegam_priors) %>%  # adding the posterior distribution
   ggplot(aes(x = year, y = biomass_mean)) +  
   facet_wrap(~mpa) +
   stat_lineribbon(aes(y = .prediction), .width = c(.95, .75, .50),  # regression line and CI
                   alpha = 0.5, colour = "black") +
   geom_point(data = bootmean_ys, colour = "darkseagreen4", size = 2) +   # raw data
   scale_fill_brewer(palette = "Greys") +
   ysab("Predicted biomass (g/m2)") +  # latin name for red knot
   xlab("Year") +
   theme_bw() +
   theme(legend.title = element_blank()))

######. Plot predicted values with mean bootstrap values ####
grouperdataboot_ys_means$year <- as.numeric(as.character(grouperdataboot_ys_means$year))
(bio_predict_comb <- bootmean_ys %>%
    tidybayes::add_predicted_draws(ys_bio_hurdlegam_priors) %>%  # adding the posterior distribution
    ggplot(aes(x = year, y = biomass_mean)) +  
    facet_wrap(~mpa) +
    stat_lineribbon(aes(y = .prediction), .width = c(.95, .75, .50),  # regression line and CI
                    alpha = 0.5, colour = "black") +
    geom_line(data = grouperdataboot_ys_means, aes(x = year, y = mean_biomass, group = mpa, color = mpa), linewidth=1) +   # raw data
    scale_fill_brewer(palette = "Greys") +
    ysab("Predicted biomass (g/m2)") +  # latin name for red knot
    xlab("Year") +
    theme_bw() +
    theme(legend.title = element_blank()))

# Plot predicted values without extreme observation
bootmean_ys$year <- as.numeric(as.character(bootmean_ys$year))
(ys_bio_plot_normals <- bootmean_ys %>%
    tidybayes::add_predicted_draws(ys_bio_hurdlegam_priors) %>%  # adding the posterior distribution
    ggplot(aes(x = year, y = biomass_mean)) +  
    facet_wrap(~mpa) +
    coord_cartesian(ysim = c(0, 60)) +
    stat_lineribbon(aes(y = .prediction), .width = c(.95, .75, .50),  # regression line and CI
                    alpha = 0.5, colour = "black") +
    geom_point(data = bootmean_ys, colour = "darkseagreen4", size = 2) +   # raw data
    scale_fill_brewer(palette = "Greys") +
    ysab("Predicted biomass (g/m)") +  # latin name for red knot
    xlab("Year") +
    theme_bw() +
    theme(legend.title = element_blank(),
          legend.position = "none",  
          axis.text.x = element_text(angle = 45, hjust = 1),  # Incline the x-axis text
          panel.grid.major = element_blank(),  # Remove major grid lines
          panel.grid.minor = element_blank(),  # Remove minor grid lines
          panel.background = element_blank(),   # Remove panel background
          strip.background = element_blank()))

# Plot predicted values extreme observations only
(ys_bio_plot_extremes <- bootmean_ys %>%
    tidybayes::add_predicted_draws(ys_bio_hurdlegam_priors) %>%  # adding the posterior distribution
    ggplot(aes(x = year, y = biomass_mean)) +  
    facet_wrap(~mpa) +
    coord_cartesian(ysim = c(60.1, 180)) +
    stat_lineribbon(aes(y = .prediction), .width = c(.95, .75, .50),  # regression line and CI
                    alpha = 0.5, colour = "black") +
    geom_point(data = bootmean_ys, colour = "darkseagreen4", size = 2) +   # raw data
    scale_fill_brewer(palette = "Greys") +
    ysab("Predicted biomass (g/m") +  # latin name for red knot
    xlab("Year") +
    scale_y_continuous(breaks = c(180)) +
    theme_bw() +
    theme(legend.title = element_blank(),
          panel.grid.major = element_blank(),  # Remove major grid lines
          panel.grid.minor = element_blank(),  # Remove minor grid lines
          panel.background = element_blank(),  # Remove panel background
          axis.text.x.bottom = element_blank(),  # Remove y-axis text labels
          axis.ticks.x.bottom = element_blank(),  # Remove y-axis ticks
          axis.title.x.bottom = element_blank(),
          axis.title.y = element_blank()))

# Combine the plots
combined_plot <- ys_bio_plot_extremes / ys_bio_plot_normals+ plot_layout(heights = c(0.3, 2))

# Plot predicted values line

# Mean density for each 'year' and 'zone' combination
bootmean_ys_raw <- bootmean_ys %>%
  group_by(year, mpa) %>%
  summarize(mean_density = mean(density_mean),
            mean_biomass = mean(biomass_mean), 
            .groups = 'drop')


(model_fit <- bootmean_ys %>%
    tidybayes::add_predicted_draws(ys_bio_hurdlegam_priors) %>%  # adding the posterior distribution
    ggplot(aes(x = year, y = biomass_mean)) +  
    facet_wrap(~mpa) +
    stat_lineribbon(aes(y = .prediction), .width = c(.95, .75, .50),  # regression line and CI
                    alpha = 0.5, colour = "black") +
    geom_line(data = bootmean_ys_raw, aes(x = year, y = mean_biomass), colour = "darkseagreen4", size = 1) +
    scale_fill_brewer(palette = "Greys") +
    ysab("Predicted biomass") +  # latin name for red knot
    xlab("Year") +
    theme_bw() +
    theme(legend.title = element_blank()))

######. Plot predicted values with mean bootstrap values and sd ####

bootmean_ys$year <- as.numeric(as.character(bootmean_ys$year))
grouperdataboot_ys_means$year <- as.numeric(as.character(grouperdataboot_ys_means$year))

graph_biomass_2023 <- ggplot() +
  # Elements from dens_predict_comb (base layers)
  stat_lineribbon(data = bootmean_ys %>% tidybayes::add_predicted_draws(ys_bio_hurdlegam_priors),
                  aes(x = year, y = .prediction, group = mpa),
                  .width = c(.95, .75, .50), alpha = 0.5, colour = "black") +
  
  # Separate geom_ribbon layers for each 'mpa' group
  geom_ribbon(data = subset(grouperdataboot_ys_means, mpa == "protected"), 
              aes(x = year, ymin = pmax(mean_biomass - mean_sd_biomass, 0), ymax = mean_biomass + mean_sd_biomass),
              fill = "turquoise", alpha = 0.15) +
  geom_ribbon(data = subset(grouperdataboot_ys_means, mpa == "unprotected"), 
              aes(x = year, ymin = pmax(mean_biomass - mean_sd_biomass, 0), ymax = mean_biomass + mean_sd_biomass),
              fill = "orange", alpha = 0.15) +
  
  # Add dashed lines for the upper boundary, setting custom colors for color based on 'mpa'
  
  geom_line(data = subset(grouperdataboot_ys_means, mpa == "protected"), 
            aes(x = year, y = mean_biomass + mean_sd_biomass, group = mpa),
            linetype = "dashed", size = 0.3, color = "turquoise", alpha=0.8) +
  geom_line(data = subset(grouperdataboot_ys_means, mpa == "unprotected"), 
            aes(x = year, y = mean_biomass + mean_sd_biomass, group = mpa),
            linetype = "dashed", size = 0.3, color = "orange", alpha=0.8) +
  
  # Add dashed lines for the lower boundary, setting custom colors for color based on 'mpa'
  geom_line(data = subset(grouperdataboot_ys_means, mpa == "protected"), 
            aes(x = year, y = pmax(mean_biomass - mean_sd_biomass, 0), group = mpa),
            linetype = "dashed", size = 0.3, color = "turquoise", alpha=0.8) +
  geom_line(data = subset(grouperdataboot_ys_means, mpa == "unprotected"), 
            aes(x = year, y = pmax(mean_biomass - mean_sd_biomass, 0), group = mpa),
            linetype = "dashed", size = 0.3, color = "orange", alpha=0.8) +
  
  # Add solid lines for mean density, setting custom colors for color based on 'mpa'
  geom_line(data = subset(grouperdataboot_ys_means, mpa == "protected"), 
            aes(x = year, y = mean_biomass, group = mpa),
            color = "turquoise", 
            size = 0.8, alpha = 0.9) +
  geom_line(data = subset(grouperdataboot_ys_means, mpa == "unprotected"), 
            aes(x = year, y = mean_biomass, group = mpa),
            color = "orange", 
            size = 0.8, alpha = 0.9) +
  
  # Common elements
  facet_wrap(~mpa) +
  scale_fill_brewer(palette = "Greys") +
  ysab("Biomass (g/m²)") +
  xlab("Year") +
  theme_bw() +
  scale_x_continuous(breaks = 2008:2023)+
  theme(legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1), 
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        strip.background = element_blank())

# Combine the plots
combined_bioplot_model <- graph_biomass_2023_2 / graph_biomass_2023_1 + plot_layout(heights = c(0.3, 2))

# - - - - - - - - - - - - - - - - 
# - - - - - - - - - - - - - - - - 
# - - - - - - - - - - - - - - - - 
# - - - - - - - - - - - - - - - - 
####_####
####__ Life Stage Histogram_________ ####
####_####

#Rename life-stage classes
grouperdataboot_ys_class$class <- recode_factor(grouperdataboot_ys_class$class, Juvenile_1 = "Juvenile",
                                               Juvenile_2 = "Subadult",
                                               Reproductive_Matrix = "Mega-spaw")
#Add yearcat column
grouperdataboot_ys_class$yearcat <- cut(as.numeric(grouperdataboot_ys_class$year), breaks=c(2007, 2010, 2013, 2016, 2019, 2023),
                                       labels=c('2008-2010', '2011-2013', '2014-2016', '2017-2019','2020-2023'))

grouperdataboot_ys_class_mean <- MeanSdQuantile(grouperdataboot_ys_class, columns = c("yearcat", "mpa", "class"), variables = c("total_length", "density", "biomass", "biomass_sum")) %>% 
  mutate_if(is.numeric, ~replace_na(., 0))

# Remove NAs
grouperdataboot_ys_class <- grouperdataboot_ys_class %>% filter(!is.na(class))


###### Density histogram #####

# Density :: Year and mpa relative density for each loop plot #
bootstrap_pz_class <- Bootstrap_pz_class
# Remove NAs
bootstrap_pz_class <- bootstrap_pz_class %>% filter(!is.na(class))
bootstrap_pz_class$class <- factor(bootstrap_pz_class$class, levels = c("Recruit", "Juvenile_1", "Juvenile_2", "Adult", "Reproductive_Matrix"))    

##### . Histogram #####
#Rename life-stage names
bootstrap_pz_class$class <- recode_factor(bootstrap_pz_class$class, Juvenile_1 = "Juvenile",
                                          Juvenile_2 = "Subadult",
                                          Reproductive_Matrix = "Mega-spa")

# Relevel class to ensure their order
bootstrap_pz_class$class <- factor(bootstrap_pz_class$class, levels = c("Recruit", "Juvenile", "Subadult", "Adult", "Mega-spa"))    


# Relevel 'mpa' to ensure 'protected' comes before 'unprotected'
bootstrap_pz_class$mpa <- factor(bootstrap_pz_class$mpa, levels = c("protected", "unprotected"))

# Use facet_grid to arrange the facets in two dimensions
hist_dens <- ggplot(bootstrap_pz_class, aes(x = class, y = density, fill = as.factor(mpa))) +
  geom_bar(stat = "summary", fun = "mean", width = 0.95) +
  geom_point(size = 0.1, alpha = 0.2, colour = "black") +
  facet_grid(mpa~yearcat) + # Rows by 'mpa', columns by 'yearcat'
  scale_fill_manual(values = c("protected" = "darkturquoise", "unprotected" = "darkorange"),
                    name = "MPA Status") +
  labs(title = "Density of Class by Year",
       x = "Size Class", y = "Density") +
  theme_bw() +
  theme(legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 10),
        plot.title = element_text(hjust = 0.5))

# Print the plot
print(hist_dens)

# - - - - - - - - - - - - - - - - 
# - - - - - - - - - - - - - - - - 
# - - - - - - - - - - - - - - - - 
# - - - - - - - - - - - - - - - - 
####_####
####__ Life Stage models _________ #####
####_####

#### Data prep ####
# data: grouperSheet <- subset(grouperSheet, abundance == 0)
class_mult <- grouperSheet %>% filter(!is.na(class))
class_mult_mean <- class_mult %>% select(yearcat, mpa, location, class) %>%  distinct()

#Rename life-stage names
class_mult_mean$class <- recode_factor(class_mult_mean$class, Juvenile_1 = "Juvenile",
                                               Juvenile_2 = "Subadult",
                                               Reproductive_Matrix = "Mega-spa")

# Define base levels
class_mult_mean$mpa <- as.factor(class_mult_mean$mpa)
class_mult_mean$mpa <- relevel(class_mult_mean$mpa, ref = "unprotected")
class_mult_mean$class <- relevel(class_mult_mean$class, ref = "Juvenile")
class_mult_mean$yearcat <- relevel(class_mult_mean$yearcat, ref = "2014-2016")

#Intercept is: Adults, mpaunprotected and yearcat 2014-2016/year2008 
multinom_model_pz <- multinom(class ~ yearcat + mpa, data = class_mult_mean)
model_summary<-summary(multinom_model_pz)
coef_df<-as.data.frame(model_summary$standard.errors)
clipr::write_clip(coef_df)
# Calculate type II Wald chi-square tests for the coefficients of the model with the bes AIC:
anova_results_pz <- Anova(multinom_model_pz, type="II")
summary(anova_results_pz) # p-values close to zero indicate that the effects of the predictors in the response variable are statistically significant

# Z-values and p-values
z_values <- model_summary$coefficients / model_summary$standard.errors
p_values <- (1 - pnorm(abs(z_values))) * 2
clipr::write_clip(p_values)

#### Pearson residuals ####
# Predicted probabilities
predicted_probs <- predict(multinom_model_pz, type = "probs")

# Dataframe with original and the predicted probabilities
data_with_probs <- cbind(class_mult_mean, predicted_probs)
dim(predicted_probs)
# Convert to a one-hot encoded matrix
actual_classes <- model.matrix(~ class - 1, data = data_with_probs)

# Calculate the residuals
residuals <- actual_classes - predicted_probs

# Since the residuals are in matrix form, we can convert them to a more manageable data frame
residuals_df <- as.data.frame(residuals)

# View the structure of the residuals data frame
str(residuals_df)

#Plot residuals
residuals_graphs <- list()

for (i in colnames(residuals_df)) {
  # Use a local function to force immediate evaluation
  residuals_graphs[[i]] <- local({
    column_data <- as.numeric(residuals_df[[i]])
    ggplot() +
      geom_histogram(aes(x = column_data), binwidth = 0.1, fill = 'blue', color = 'black') +
      labs(x = paste0("Residuals for Class ", i), y = "Count") +
      theme_minimal()
  })
}
wrap_plots(residuals_graphs)

##### Coefficients plot #####
# Model coefficients and tidy data frame
coef_df <- as.data.frame(coef(multinom_model_pz))
coef_df$Variable <- rownames(coef_df)
coef_df <- tidyr::pivot_longer(coef_df, cols = -Variable, names_to = "Class", values_to = "Estimate")

# plot
ggplot(coef_df, aes(x = Variable, y = Estimate)) +
  geom_point() +
  facet_wrap(~Class, scales = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Variable", y = "Coefficient Estimate", title = "Coefficients of Multinomial Logistic Regression Model")


##### Coefficients probability plot #####
# Convert the log-odds to probabilities
coef_df$Probability <- 1 / (1 + exp(-coef_df$Estimate))

# Plot the probabilities instead of the log-odds
ggplot(coef_df, aes(x = Variable, y = Probability)) +
  geom_point() +
  facet_wrap(~Class, scales = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Variable", y = "Probability", title = "Probabilities from Multinomial Logistic Regression Model")


##### Barplot probabilities #####

# Convert log-odds to probabilities for all coefficients
log_odds_matrix <- coef(multinom_model_pz)
probabilities_matrix <- 1 / (1 + exp(-log_odds_matrix))
# Probabilities plots
barplot(as.matrix(probabilities_matrix), beside = TRUE, legend = rownames(probabilities_matrix))

##### Predicted probabilities vs original #####

# Calculate predicted probabilities
pred_probs <- predict(multinom_model_pz, type = "probs")
# Convert to data frame and add actual class for comparison
pred_probs_df <- data.frame(pred_probs, ActualClass = class_mult_mean$class)
# Melt datarame
pred_probs_df_long <- reshape2::melt(pred_probs_df, id.vars = "ActualClass")

# Plot predicted probabilities
ggplot2::ggplot(pred_probs_df_long, ggplot2::aes(x = ActualClass, y = value, fill = variable)) +
  ggplot2::geom_bar(stat = "identity", position = "dodge") +
  ggplot2::theme_minimal() +
  ggplot2::labs(x = "Actual Class", y = "Predicted Probability", fill = "Predicted Class")

##### Each class tendency ####
# Data frame for prediction
predict_dt <- expand.grid(yearcat=unique(class_mult_mean$yearcat),
                          mpa=unique(class_mult_mean$mpa),
                          class=levels(as.factor(class_mult_mean$class)))

# Predict probabilities
predict_probs <- predict(multinom_model_pz, newdata=predict_dt, type="probs")

# bind data
predict_dt <- cbind(predict_dt[, -which(names(predict_dt) == "class")], predict_probs)

#transform to ggplot format
long_predict_data <- tidyr::pivot_longer(
  predict_dt,
  cols = c("Adult", "Juvenile", "Subadult", "Recruit", "Mega-spa"),  
  names_to = "class",
  values_to = "probability"
)

# Sample code to create a numeric yearcat representation
long_predict_data$yearcat_numeric <- as.numeric(factor(long_predict_data$yearcat, 
                                                       levels = c("2008-2010", "2011-2013", "2014-2016", "2017-2019", "2020-2023"),
                                                       labels = c(2009, 2012, 2015, 2018, 2021)))


# Use ggplot2 to plot the values
ggplot(long_predict_data, aes(x = yearcat_numeric, y = probability, color = class)) +
  facet_wrap(~mpa) +
  geom_point(size=0.7) +
  geom_line(size=1) + # Add a trend line
  scale_x_continuous(name = "Year period", breaks = c(1, 2, 3, 4, 5), labels = c("2008-2010", "2011-2013", "2014-2016", "2017-2019", "2020-2023")) +
  labs(y = "Predicted Probabilities", color = "Class") +
  theme_bw() +
  theme(legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),   
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        strip.background = element_blank())

# Use ggplot2 to plot trends 
ggplot(long_predict_data, aes(x = yearcat_numeric, y = probability, color = class)) +
  facet_wrap(~mpa) +
  geom_point(size=0.7) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 3)) +  # Add a trend line
  scale_x_continuous(name = "Year period", breaks = c(1, 2, 3, 4, 5), labels = c("2008-2010", "2011-2013", "2014-2016", "2017-2019", "2020-2023")) +
  labs(y = "Predicted Probabilities", color = "Class") +
  theme_bw() +
  theme(legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),   
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        strip.background = element_blank())


##### GAM residuals visualization ####

# GAM model for each class and mpa combination
gam_models <- list()
for (class_name in unique(long_predict_data$class)) {
  for (mpa_status in unique(long_predict_data$mpa)) {
    # Subset data for each class and mpa
    subset_data <- subset(long_predict_data, class == class_name & mpa == mpa_status)
    
    # Fit the GAM model
    gam_model <- gam(probability ~ s(yearcat_numeric, k = 3), data = subset_data)
    
    # Calculate residuals and add them to the dataframe
    subset_data$residuals <- residuals(gam_model)
    
    # Store the model and modified data
    gam_models[[paste(class_name, mpa_status, sep = "_")]] <- list(model = gam_model, data = subset_data)
  }
}

# Combine residuals data into a single dataframe
all_residuals_data <- do.call(rbind, lapply(gam_models, function(x) x$data))

#Plot residuals
ggplot(all_residuals_data, aes(x = yearcat_numeric, y = residuals, color = class)) +
  facet_wrap(~mpa) +
  geom_point(size=0.7) +
  geom_hline(yintercept = 0, linetype = "dashed") + # Add a horizontal line at y=0
  scale_x_continuous(name = "Year period", breaks = c(1, 2, 3, 4, 5), labels = c("2008-2010", "2011-2013", "2014-2016", "2017-2019", "2020-2023")) +
  labs(y = "Residuals", color = "Class") +
  theme_bw() +
  theme(legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),   
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        strip.background = element_blank())

