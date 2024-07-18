#loading in the necessary libraries

library(tidyverse)
library(sf)
library(tmap)
library(janitor)
library(spatstat)
library(spdep)
library(rstan)
library(geostan)
library(SpatialEpi)
library(tidybayes)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

#loading in the shapefile
full_sf<- read_sf("Shape/LSOA_2011_EW_BFC.shp")

#creating seperate shapefiles for westminister and camden to demarcate the areas on a map
west_sf <- dplyr::filter(full_sf, grepl("Westminster", LSOA11NM))
cam_sf <- dplyr::filter(full_sf, grepl("Camden", LSOA11NM))

#filtering the shapefile in such a manner that only wards from the borough of Camden and Westminster remain
clean_sf <- dplyr::filter(full_sf, grepl("Camden|Westminster", LSOA11NM))

#plotting the shapefiles together and highlighting the wards of the two boroughs
tm_shape(clean_sf) +
  tm_polygons(col = NA, alpha = 0.5) +
tm_shape(west_sf) +
  tm_polygons(col = 'red', alpha = 0.5) +
  tm_shape(cam_sf) +
  tm_polygons(col = 'blue', alpha = 0.5)

#loading in the population data
population <- read_csv("population.csv", show_col_types = FALSE)

#loading in the IMD Score
imd <- read_csv("IMD.csv", show_col_types = FALSE)

#Loading in the average house price value (updated till beginning of 2018)
house_price <- read_csv("house_price.csv", show_col_types = FALSE)

#Merging the aforementioned CSVs with the cleaned shapefile
clean_sf <- merge(clean_sf, population, by.x = c("LSOA11CD", "LSOA11NM"), by.y = c("lsoa_code", "lsoa_name"), all.x = TRUE)
clean_sf <- merge(clean_sf, imd, by.x = c("LSOA11CD", "LSOA11NM"), by.y = c("lsoa_code", "lsoa_name"), all.x = TRUE)
clean_sf <- merge(clean_sf, house_price, by.x = c("LSOA11CD", "LSOA11NM"), by.y = c("lsoa_code", "lsoa_name"), all.x = TRUE)

#Adding in crime data with crime commiteed over the last month
data <- read_csv("2023-03-metropolitan-street.csv", show_col_types = FALSE)

#cleaning data names in order to have a seamless merge with the shapefile
clean_data<-data%>%
  clean_names()

#filering the data to create a dataframe with crimes that contain theft from a person
mug <- clean_data%>%
  filter(str_detect(crime_type, 'Theft from the person'))

#creating a count of all the incidents in an area by counting the repetiton of an LSOA code
#creating a new table with the values
counts <- table(mug$lsoa_code)
incident_table <- data.frame(lsoa_code = names(counts), total_incidents = as.numeric(counts))

#Merging the incidents of theft with the cleaned shapefile
clean_sf <- merge(clean_sf, incident_table, by.x = ("LSOA11CD"), by.y = ("lsoa_code"), all.x = TRUE)

#As opposed to dropping NA, the missing values have been filled with 0 to avoid the shapefile being incomplete and not rendering a map
clean_sf <- replace(clean_sf, is.na(clean_sf),0)

#Creating a new column in the shapefile with expected incidents of theft using the expected function
clean_sf$expinc <- round(expected(population = clean_sf$population, cases = clean_sf$total_incidents, n.strata = 1), 0)

#Creating a spatial object to extract nodes and edges from 
sp.object <- as(clean_sf, "Spatial")

# converting the spatial object into an adjacency matrix and extracting nodes and edges
adjacencyMatrix <- shape2mat(sp.object)

# extracting the requried components for the ICAR model
extractComponents <- prep_icar_data(adjacencyMatrix)

#assigning features of the matrix as variables to use in the Stan code
n <- as.numeric(extractComponents$group_size)
nod1 <- extractComponents$node1
nod2 <- extractComponents$node2
n_edges <- as.numeric(extractComponents$n_edges)

#assigning shapefile columns to be used in the Stan code
y <- clean_sf$total_incidents
x <- st_drop_geometry(clean_sf[c(3:4)])
e <- clean_sf$expinc

#creating a spatial dataset for Stan
clean_sf.dataset <- list(N=n, N_edges=n_edges, node1=nod1, node2=nod2, Y=y, X=x, E=e, k=2)

icar_poisson_fit = stan("mug_two.stan", data=clean_sf.dataset, iter=20000, chains=6, verbose = FALSE)

#preventing the exponential notation of integers
options(scipen = 999)

#summarising the alpha, beta and sigma values after running the Stan code on the dataset
summary(icar_poisson_fit, pars=c("alpha", "beta", "sigma"), probs=c(0.025, 0.975))$summary

#getting a shorter summary of the phi values  using the head function
head(summary(icar_poisson_fit, pars=c("phi"), probs=c(0.025, 0.975))$summary)

#creating a new dataframe to ensure that the rHAT is below 1.1 for all values
diagnostic.checks <- as.data.frame(summary(icar_poisson_fit, pars=c("alpha", "beta", "sigma", "phi", "lp__"), probs=c(0.025, 0.5, 0.975))$summary)
# creating a  binary variable
diagnostic.checks$valid <- ifelse(diagnostic.checks$Rhat < 1.1, 1, 0)
# changing the variable into a table where 1 means the value is under the rHat of 1.1
table(diagnostic.checks$valid)

# summarising the mu of the stanfit
head(summary(icar_poisson_fit, pars=c("mu"), probs=c(0.025, 0.975))$summary)

#removing posterior results to create a new dataframe
relativeRisk.results <- as.data.frame(summary(icar_poisson_fit, pars=c("mu"), probs=c(0.025, 0.975))$summary)
#adding row numbers to a new dataframe and selecting the 4 required columns
row.names(relativeRisk.results) <- 1:nrow(relativeRisk.results)
relativeRisk.results <- relativeRisk.results[, c(1,4,5,7)]
# third, rename the columns appropriately
colnames(relativeRisk.results)[1] <- "rr"
colnames(relativeRisk.results)[2] <- "rrlower"
colnames(relativeRisk.results)[3] <- "rrupper"
colnames(relativeRisk.results)[4] <- "rHAT"

#creating the risk map by adding the results to the shapefile
clean_sf$rr <- relativeRisk.results[, "rr"]
clean_sf$rrlower <- relativeRisk.results[, "rrlower"]
clean_sf$rrupper <- relativeRisk.results[, "rrupper"]

#creating categories to define if an area has if there will be a ssignificant risk, if at all any, in an area
clean_sf$Significance <- NA
clean_sf$Significance[clean_sf$rrlower<1 & clean_sf$rrupper>1] <- 0    # NOT SIGNIFICANT
clean_sf$Significance[clean_sf$rrlower==1 | clean_sf$rrupper==1] <- 0  # NOT SIGNIFICANT
clean_sf$Significance[clean_sf$rrlower>1 & clean_sf$rrupper>1] <- 1    # SIGNIFICANT INCREASE
clean_sf$Significance[clean_sf$rrlower<1 & clean_sf$rrupper<1] <- -1   # SIGNIFICANT DECREASE

#creating labels to add to the map
RiskCategorylist <- c(">0.0 to 0.25", "0.26 to 0.50", "0.51 to 0.75", "0.76 to 0.99", "1.00 & <1.01",
                      "1.01 to 1.10", "1.11 to 1.25", "1.26 to 1.50", "1.51 to 1.75", "1.76 to 2.00", "2.01 to 3.00")

#creating a colour palatte from colourbrewer

RRPalette <- c("#1a1a1a","#4d4d4d","#878787","#bababa","white","#e0e0e0","#fddbc7","#f4a582","#d6604d","#b2182b","#67001f")

#assigning the risk values to the labelling in the RiskCategorylist object
clean_sf$RelativeRiskCat <- NA
clean_sf$RelativeRiskCat[clean_sf$rr>= 0 & clean_sf$rr <= 0.25] <- -4
clean_sf$RelativeRiskCat[clean_sf$rr> 0.25 & clean_sf$rr <= 0.50] <- -3
clean_sf$RelativeRiskCat[clean_sf$rr> 0.50 & clean_sf$rr <= 0.75] <- -2
clean_sf$RelativeRiskCat[clean_sf$rr> 0.75 & clean_sf$rr < 1] <- -1
clean_sf$RelativeRiskCat[clean_sf$rr>= 1.00 & clean_sf$rr < 1.01] <- 0
clean_sf$RelativeRiskCat[clean_sf$rr>= 1.01 & clean_sf$rr <= 1.10] <- 1
clean_sf$RelativeRiskCat[clean_sf$rr> 1.10 & clean_sf$rr <= 1.25] <- 2
clean_sf$RelativeRiskCat[clean_sf$rr> 1.25 & clean_sf$rr <= 1.50] <- 3
clean_sf$RelativeRiskCat[clean_sf$rr> 1.50 & clean_sf$rr <= 1.75] <- 4
clean_sf$RelativeRiskCat[clean_sf$rr> 1.75 & clean_sf$rr <= 2.00] <- 5
clean_sf$RelativeRiskCat[clean_sf$rr> 2.00 & clean_sf$rr <= 10] <- 6

#creating a plot relative risk
rr_map <- tm_shape(clean_sf) + 
  tm_fill("RelativeRiskCat", style = "cat", title = "Relavtive Risk of Theft", palette = RRPalette, labels = RiskCategorylist) +
  tm_polygons(alpha = 0, border.alpha = 1, border.col = "orange") +
  tm_layout(frame = FALSE, legend.outside = TRUE, legend.title.size = 0.8, legend.text.size = 0.7) +
  tm_compass(position = c("right", "top")) + tm_scale_bar(position = c("right", "bottom"))

#creating a plot of significance regions
sg_map <- tm_shape(clean_sf) + 
  tm_fill("Significance", style = "cat", title = "Significance of Threat", 
          palette = c("#33a6fe", "white", "#fe0000"), labels = c("Significantly low", "Not Significant", "Significantly high")) +
  tm_polygons(alpha = 0, border.alpha = 1, border.col = "orange") +
  tm_layout(frame = FALSE, legend.outside = TRUE, legend.title.size = 0.8, legend.text.size = 0.7) +
  tm_compass(position = c("right", "top")) + tm_scale_bar(position = c("right", "bottom"))

#generating the plots together
tmap_arrange(rr_map, sg_map, ncol = 2, nrow = 1)

#extracting  exceedence probabilities from the icar_possion_fit object
#calculating the probability of areas having a relative risk ratio > 1.0
threshold <- function(x){mean(x > 1.00)}
excProbrr <- icar_poisson_fit %>% spread_draws(mu[i]) %>% 
  group_by(i) %>% summarise(mu=threshold(mu)) %>%
  pull(mu)

#adding  exceedance values to the clean_sf
clean_sf$excProb <- excProbrr

#generating labels as done earlier
ProbCategorylist <- c("<0.01", "0.01-0.09", "0.10-0.19", "0.20-0.29", "0.30-0.39", "0.40-0.49","0.50-0.59", "0.60-0.69", "0.70-0.79", "0.80-0.89", "0.90-0.99", "1.00")

clean_sf$ProbCat <- NA
clean_sf$ProbCat[clean_sf$excProb>=0 & clean_sf$excProb< 0.01] <- 1
clean_sf$ProbCat[clean_sf$excProb>=0.01 & clean_sf$excProb< 0.10] <- 2
clean_sf$ProbCat[clean_sf$excProb>=0.10 & clean_sf$excProb< 0.20] <- 3
clean_sf$ProbCat[clean_sf$excProb>=0.20 & clean_sf$excProb< 0.30] <- 4
clean_sf$ProbCat[clean_sf$excProb>=0.30 & clean_sf$excProb< 0.40] <- 5
clean_sf$ProbCat[clean_sf$excProb>=0.40 & clean_sf$excProb< 0.50] <- 6
clean_sf$ProbCat[clean_sf$excProb>=0.50 & clean_sf$excProb< 0.60] <- 7
clean_sf$ProbCat[clean_sf$excProb>=0.60 & clean_sf$excProb< 0.70] <- 8
clean_sf$ProbCat[clean_sf$excProb>=0.70 & clean_sf$excProb< 0.80] <- 9
clean_sf$ProbCat[clean_sf$excProb>=0.80 & clean_sf$excProb< 0.90] <- 10
clean_sf$ProbCat[clean_sf$excProb>=0.90 & clean_sf$excProb< 1.00] <- 11
clean_sf$ProbCat[clean_sf$excProb == 1.00] <- 12


#mapping exceedance probabilities
tm_shape(clean_sf) + 
  tm_fill("ProbCat", style = "cat", title = "Probability", palette = "OrRd", labels = ProbCategorylist) +
  tm_polygons(alpha = 0, border.alpha = 1, border.col = "black") +
  tm_layout(frame = FALSE, legend.outside = TRUE, legend.title.size = 0.8, legend.text.size = 0.7) +
  tm_compass(position = c("right", "top")) + tm_scale_bar(position = c("right", "bottom"))

