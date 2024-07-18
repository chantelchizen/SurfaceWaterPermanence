#### Model Wetland Pond Classification
#### Use surface water permanence to model Stewart and Kantrud (1971) wetland classes

#### Last updated: July 9, 2024
#### Author: Chantel Chizen


# Packages ----

library(caret)
library(dplyr)
library(gstat)
library(ROSE)
library(sf)
library(sp)

# Constants ----

# define wetland class transformation lists based on Steward & Kantrud (1971) classification system
months_list <- list(0, 1, 2, 3, 4, 5, 6, 7)
class_list <- list(1, 2, 2, 3, 3, 4, 4, 5)

# set up train control for the random forest model
fitControl <- trainControl(method = "cv",
                           number = 10,
                           sampling = "smote")

# Functions ----

# convert pond permanence to wetland class
transform_pond_perm <- function(value) {
  if (value %in% unlist(months_list)) {
    return(unlist(class_list)[match(value, unlist(months_list))])
  } else {
    return(NA)  # Or some default value
  }
}

# Load wetland data ----

## provide csv file path here ====

wetland_class_summary <- read.csv("wetland_total_pondpermanence.csv")
#wetland_class_summary <- read.csv("wetland_consecutive_pondpermanence.csv")

wetland_class_summary <- wetland_class_summary[complete.cases(wetland_class_summary[, "year_2017"]), ] %>%
  subset(select = -c(X, gee, class_1))

wetland_class_summary$class <- factor(wetland_class_summary$class)

## provide shape file path here ====
site_polygons <- st_read("site_polygons.shp") %>%
  mutate(wetland = as.numeric(wetland)) %>% # convert wetland # to numeric
  arrange(wetland) %>% # reorder list in wetland order
  st_make_valid()

site_polygons$centroids <- st_coordinates(st_centroid(site_polygons$geometry))

site_polygons$lat <- site_polygons$ centroids[,2]
site_polygons$long <- site_polygons$ centroids[,1]

site_polygons <- site_polygons %>%
  subset(select = c(wetland, lat, long)) %>%
  data.frame()

wetland_class_summary <- left_join(wetland_class_summary, site_polygons, by = "wetland")

# Calculate mean pond permanence over all years  ----

wetland_class_summary <- wetland_class_summary  %>%
  rowwise() %>%
  mutate(pond_perm_mean = round(mean(c_across(starts_with("year")))))

wetland_class_summary$pond_perm_mean <- factor(wetland_class_summary$pond_perm_mean)


# Convert pond permanence to wetland class ----

# Create new columns with transformations
wetland_class_summary <- wetland_class_summary %>%
  mutate(across(starts_with("pond_perm"), ~ transform_pond_perm(.), .names = "class_{col}"))

wetland_class_summary$class_pond_perm_mean <- factor(wetland_class_summary$class_pond_perm_mean, levels = c(1:5))
wetland_class_summary$class <- factor(wetland_class_summary$class, levels = c(1:5))

# assign soil zone
wetland_class_summary$soilzone <- plyr::revalue(wetland_class_summary$layer, c("kantrud_armriver" = "dark_brown", "kantrud_bauche" = "black", "kantrud_hebert"="black"))

wetland_class_summary_m <- wetland_class_summary %>%
  subset(select = -c(wetland, layer, geometry)) %>%
  mutate(class = as.numeric(class))


# Model wetland class ----

## pond permanence only ====

confusionMatrix(data = wetland_class_summary$class_pond_perm_mean,
                reference = wetland_class_summary$class)

## predicted wetland class based on mean ====

# Linear model ====

class_lm <- lm(class ~ class_pond_perm_mean, data = wetland_class_summary_m)

summary(class_lm)
sigma(class_lm)
summary(class_lm)$r.squared


# Linear model with soil zone ====

class_lm <- lm(class ~ class_pond_perm_mean * soilzone * area, data = wetland_class_summary_m)

summary(class_lm)
sigma(class_lm)
summary(class_lm)$r.squared


# Machine learning models -----

wetland_class_summary_m$soilzone <-
  plyr::revalue(wetland_class_summary_m$soilzone, c("black" = "1", "dark_brown" = "2"))

wetland_class_summary_m$impact <-
  plyr::revalue(wetland_class_summary_m$impact, c("cultivated" = "1", "none" = "2"))

wetland_class_summary_m <- wetland_class_summary_m %>%
  mutate(class_pond_perm_mean = as.numeric(class_pond_perm_mean),
         pond_perm_mean = as.numeric(pond_perm_mean),
         class = factor(class, levels = c(1:5)),
         soilzone = as.numeric(soilzone),
         impact = as.numeric(impact)
         )

# Set seed for reproducibility
set.seed(123)

trainIndex <- createDataPartition(wetland_class_summary_m$class, p = 0.7,
                                  list = FALSE,
                                  times = 1)

train_data_sp <- wetland_class_summary_m[trainIndex,]
test_data_sp  <- wetland_class_summary_m[-trainIndex,]

table(train_data_sp$class)

train_data <- train_data_sp %>%
  subset(select = -c(lat, long))

test_data <- test_data_sp %>%
  subset(select = -c(lat, long))


# variogram for spatial autocorrelations

coordinates(train_data_sp) <- ~ long+lat
coordinates(test_data_sp) <- ~ long+lat

variogram_train <- variogram(class ~ 1, train_data_sp)
variogram_test <- variogram(class ~ 1, test_data_sp)

# Plot the variograms
plot(variogram_train, main="Variogram - Training Data")
plot(variogram_test, main="Variogram - Validation Data")

## random forest ====

set.seed(123)

rf.model = train(class ~ ., data = train_data, method = "rf",
                 trControl = fitControl,
                 na.action = na.omit)

rf.pred <- predict(rf.model,test_data)

## CART ====

set.seed(123)

cart.model = train(class ~ ., data = train_data, method = "rpart",
                   trControl = fitControl,
                   na.action = na.omit,)

cart.pred <- predict(cart.model,test_data)


## compare models ====

confusionMatrix(rf.pred, test_data$class) # random forest performed the best
confusionMatrix(cart.pred, test_data$class)

gbmImp <- varImp(rf.model, scale = FALSE)
gbmImp$importance$relative <- format(round(gbmImp$importance$Overall / sum(gbmImp$importance$Overall), 2), nsmall = 2)

gbmImp$importance[order(gbmImp$importance$relative, decreasing = T),]


save(rf.model, file = "~/git/WetlandClassification/processed_data/rfmodel.RData")

# Remove ephemeral wetland class ----

wetland_class_summary_ne <- wetland_class_summary_m %>%
  subset(class != 1) %>%
  mutate(class = factor(class, levels = c(2,3,4,5)))

set.seed(123)

trainIndex <- createDataPartition(wetland_class_summary_ne$class, p = 0.7,
                                  list = FALSE,
                                  times = 1)

train_data <- wetland_class_summary_ne[trainIndex,]
test_data  <- wetland_class_summary_ne[-trainIndex,]

table(train_data$class)


set.seed(2000)

## random forest ====
rf.model.ne = train(class ~ ., data = train_data, method = "rf",
                 trControl = fitControl,
                 na.action = na.omit)

rf.pred.ne <- predict(rf.model.ne,test_data)

## compare models ====

confusionMatrix(rf.pred.ne, test_data$class) # random forest performed the best


# Remove soil zone ----

wetland_class_summary_ns <- wetland_class_summary_m %>%
  subset(select = -c(soilzone))

set.seed(123)

trainIndex <- createDataPartition(wetland_class_summary_ns$class, p = 0.7,
                                  list = FALSE,
                                  times = 1)

train_data <- wetland_class_summary_ns[trainIndex,]
test_data  <- wetland_class_summary_ns[-trainIndex,]

table(train_data$class)


set.seed(123)

## random forest ====
rf.model.ns = train(class ~ ., data = train_data, method = "rf",
                    trControl = fitControl,
                    na.action = na.omit)

rf.pred.ns <- predict(rf.model.ns,test_data)


## compare models ====

confusionMatrix(rf.pred.ns, test_data$class) # random forest performed the best

save(rf.model.ns, file = "~/git/WetlandClassification/processed_data/rfmodel_nosoilzone.RData")
