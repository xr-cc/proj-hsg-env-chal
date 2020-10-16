### Load libraries
library(tidyverse)
library(googledrive)
library(openxlsx)
library(Amelia)
library(xtable)
library(sp)
library(spdep)
library(rgdal)
library(spatialreg)
library(tmap)
library(GWmodel)

############### Retrieving Data ###############

# Link to Google Drive data file
drive_link <- "https://drive.google.com/file/d/1bLlhX3VWOcODsOobgNLWkDSlw8BILQS_/view?usp=sharing"
drive_file <- drive_get(as_id(drive_link))

# Download Files
drive_download(as_id(drive_file$id), overwrite=TRUE)

# Read Files into R
ces3 <- read.xlsx("CES3 with per capita income - Bay Area.xlsx")
hatlas <- read.csv(file.path(path,"Health_Atlas_Data.csv"))
shp_file <- st_read(file.path(path,"CESJune2018Update_SHP/CES3June2018Update.shp")) 

############### Data Cleaning and Processing ###############

# CES3.0 
ces3_cleaned <- ces3 %>% 
  dplyr::mutate(SB.535.Disadvantaged.Community=ifelse(SB.535.Disadvantaged.Community=='No',0,1)) %>%
  dplyr::rename(Per_capita_income = `Per.Capita.Income.(In.2018.Inflation.Adjusted.Dollars)`) %>%
  dplyr::mutate(Bracket = ifelse(Per_capita_income < 50000, 'Per Capita Income <$50,000', ifelse((Per_capita_income >= 50000) & (Per_capita_income < 100000), 'Per Capita Income $50,000-$99,999', 'Per Capita Income $100,000+'))) %>%
  dplyr::select(-contains(c("Pctl","Percentile",".to.","$25","$100","$50","$60","Tract.2","Tract.3",".County","County."))) %>%
  dplyr::mutate(Majority_less_75k = ifelse(`Less.than.$75,000` > 50.0, 'Majority Household Income <$75,000', 'Majority Household Income >$75,000'))

# Health Atlas
hatlas_cleaned <- hatlas %>%
  dplyr::select(contains(c("GEO.display.label","geoid","arthritis_crudeprev","bphigh_crudeprev","casthma_crudeprev","chd_crudeprev","copd_crudeprev","diabetes_crudeprev","highchol_crudeprev","kidney_crudeprev","mhlth_crudeprev","obesity_crudeprev","phlth_crudeprev","stroke_crudeprev"))) %>%
  dplyr::filter(!(geoid=="06")) %>%
  dplyr::rename(Arthritis = arthritis_crudeprev,
                High_blood_pressure= bphigh_crudeprev,
                Current_asthma = casthma_crudeprev,
                Coronary_heart_disease = chd_crudeprev,
                COPD = copd_crudeprev,
                Diabetes = diabetes_crudeprev,
                High_cholesterol = highchol_crudeprev,
                Chronic_kidney_disease = kidney_crudeprev,
                Poor_mental_health = mhlth_crudeprev,
                Obesity = obesity_crudeprev,
                Poor_physical_health = phlth_crudeprev,
                Stroke = stroke_crudeprev) %>%
  dplyr::mutate(geoid=as.numeric(sub('.', '', geoid)))

# Replace "N" with NA and convert character numerical columns to true numeric in Health Atlas
hatlas_cleaned[hatlas_cleaned == "N"] <- NA
hatlas_cleaned[,3:14] <- sapply(sapply(hatlas_cleaned[,3:14],as.character), as.numeric)

# Merge CES3.0 with Health Atlas
ces_health <- ces3_cleaned %>% inner_join(hatlas_cleaned, by = c("Census.ID" = "geoid"))

############### Using PCA to create health index / Disease Burden ###############
require(missMDA)
require(FactoMineR)

# extract health outcomes and perform PCA
crudeprev <- ces_health[,41:52]
nb <- estim_ncpPCA(crudeprev, ncp.min=0, ncp.max=5, method.cv="Kfold", nbsim=50)
imputed <- imputePCA(crudeprev, ncp=nb$ncp)
res.pca <- PCA(imputed$completeObs, scale.unit = TRUE)

# set PC1 as Disease Burden index
ces_health$Disease_burden <- res.pca$ind$coord[,1]

############### Preprocessing CES3.0 predictors before Modelling ###############
model_subset <- ces_health[,c(4:5, 10:11, 14:25, 31:35, 39, 53)]
preprocessing <- preProcess(model_subset[,-c(1:4,22:23)], method = c("center", "scale", "knnImpute","nzv")) 
model_subset_prep <- predict(preprocessing, model_subset, model_subset)

### Merging preprocessed dataset with CES3.0 shapefile
ces3_shp <- model_subset_prep %>% inner_join(shp_file, by = c("Census.ID" = "tract"))
spd <- sf::as_Spatial(st_geometry(ces3_shp$geometry), IDs = as.character(1:nrow(ces3_shp)))

df <- ces3_shp
df$geometry <- NULL
df <- as.data.frame(ces3_shp)

### create the SpatialPolygonsDataFrame and compute Queen continuity distance metric
spd <- sp::SpatialPolygonsDataFrame(spd, data = df)
nb <- poly2nb(spd, queen=TRUE)
lw <- nb2listw(nb, style="B", zero.policy=TRUE)
Wij <- as.matrix( as(lw, "symmetricMatrix") )

ces3_shp <- st_as_sf(ces3_shp)

# map of PC1/Disease Burden
tm_shape(ces3_shp) +
  tm_polygons(col = "Disease_burden", style = "cont", n=8, palette = "RdPu", border.alpha = 0.1, title = "", midpoint = 0, legend.reverse=TRUE) +
  tm_layout(main.title = "Disease Burden Index",  main.title.size = 1.1, frame = FALSE, legend.outside = FALSE) 


############### Models: OLS, Spatial Lag, GWR ###############

### initialize formula object
model <- Disease_burden ~ Ozone + PM2.5 + Diesel.PM + Drinking.Water + Pesticides + Tox..Release +
  Traffic + Cleanup.Sites + Groundwater.Threats + Haz..Waste + Imp..Water.Bodies + 
  Solid.Waste + Education + Linguistic.Isolation + Unemployment + Housing.Burden + Poverty

### OLS Model - inappropriate but try to see how much residual spatial distribution there is and check VIF
ols.fit <- lm(model, data=ces3_shp)
summary(ols.fit)
car::vif(ols.fit)

ces3_shp <- mutate(ces3_shp, olsresid = resid(ols.fit))

tm_shape(ces3_shp) +
  tm_polygons(col = c("olsresid"), style = "cont", palette = "Reds", 
              border.alpha = 0, title = "", midpoint = 0) +
  tm_scale_bar(breaks = c(0, 2.5, 5), text.size = 1, position = c("right", "bottom")) +
  tm_layout(main.title = "Residuals from linear regression",  main.title.size = 0.95, frame = FALSE, legend.outside = TRUE) +
  tm_facets(sync = TRUE, ncol = 2)

# Estimate Moran's I for Disease Burden
bay.sp <- as(ces3_shp, "Spatial")
bayb <- poly2nb(bay.sp, queen=TRUE)
bayw <- nb2listw(bayb, style="W", zero.policy = TRUE)

moran.plot(bay.sp$Disease_burden, listw=bayw, xlab="PC1", ylab="Standardized Lagged PC1",
           main=c("Moran Scatterplot for PC1", "in Bay Area"), zero.policy = TRUE )

moran.mc(bay.sp$Disease_burden, listw = bayw, nsim=999, zero.policy = TRUE)
lm.morantest(ols.fit, listw = bayw, zero.policy=TRUE)

### Spatial Lag Model / Global Model
lag.fit <- lagsarlm(model, data = bay.sp, listw = bayw, zero.policy=TRUE) 
summary(lag.fit)
moran.mc(resid(lag.fit), bayw, nsim=999, zero.policy=TRUE)

### Geographically Weighted Regression / Local Model
# Prep distance matrix based on Euclidean
bay.sp@data <- bay.sp@data[, -c(24:50)]
DM <- gw.dist(dp.locat=coordinates(bay.sp))

# Cross-validation to estimate kernal bandwith
bw.rd <- bw.gwr.lcr(model, data=bay.sp, kernel="gaussian",
                    lambda=0,lambda.adjust=TRUE,cn.thresh=30,
                    adaptive=TRUE, p=2, theta=0, longlat=F,dMat=DM)

# Fit GWR model 
gwr.ridge <- gwr.lcr(model, data=bay.sp, bw=60, dMat=DM, kernel='gaussian',
                     lambda=0, lambda.adjust=TRUE, cn.thresh=30,
                     adaptive=TRUE, p=2, theta=0, longlat=F, cv=T)

summary(gw.ridge)

# Plotting GWR coefficients
m1 <- tm_shape(gwr.ridge$SDF) +
  tm_polygons(col = c("Poverty"), style = "pretty", n=6, palette = "Purples", interactive=TRUE,
              border.alpha = 0, title = "", midpoint = 0, legend.reverse=TRUE)  +
  tm_layout(main.title = "Poverty",  main.title.size = 1.1, frame = FALSE, legend.outside = FALSE) 

m2 <- tm_shape(gwr.ridge$SDF) +
  tm_polygons(col = c("Imp..Water.Bodies"), style = "pretty", n=4, palette = "Blues", interactive=TRUE,
              border.alpha = 0, title = "", midpoint = 0, legend.reverse=TRUE)+
  tm_layout(main.title = "Impaired Water Bodies",  main.title.size = 1.1, frame = FALSE, legend.outside = FALSE) 

m3 <- tm_shape(gwr.ridge$SDF) +
  tm_polygons(col = c("Unemployment"), style = "pretty", n=3, palette = "Greens", interactive=TRUE,
              border.alpha = 0, title = "", midpoint = 0, legend.reverse=TRUE)  +
  tm_layout(main.title = "Unemployment",  main.title.size = 1.1, frame = FALSE, legend.outside = FALSE) 

png("gwr.png", width = 1100, height = 480)
tmap_arrange(m1,m2,m3, ncol = 3)
dev.off()





