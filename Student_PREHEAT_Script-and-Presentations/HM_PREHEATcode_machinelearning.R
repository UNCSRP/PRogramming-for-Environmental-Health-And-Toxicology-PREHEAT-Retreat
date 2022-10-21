## PREHEAT Retreat
# Authored by Hannah Matthews
# 10/21/2022
# Group 5: Asthma 

# Install libraries
# We will need need tidyverse, randomForest, janitor, and xgboost
# Remember you can pass a vector of library names to install multiple libraries at once
install.packages("tidyverse")
install.packages("randomForest")
install.packages("janitor")
install.packages("xgboost")


# Load libraries
library(tidyverse)
library(randomForest)
library(janitor)
library(xgboost)
...


# Load data

data_demo <- read.csv("./3. Data Analysis/Data/Demographic_Data.csv")
data_proteomics <- read.csv("./3. Data Analysis/Data/Processed_Proteomics_Data.csv")

# Change to tidy format 

data_tidy <-
  data_proteomics %>% 
  pivot_longer(cols=c(2:ncol(data_proteomics)), names_to = "SampleID", values_to = 'Value') %>% 
  pivot_wider(names_from = 'Gene_Name', values_from='Value')

# In the demographics dataframe, let's keep only our SampleID and variable to predict on

data_demo <-
  data_demo %>% 
  select("SampleID", "Asthma_Status")

# Join the two dataframes on SampleID, preserving all the rows in proteomics (what kind of join is this? Type ?join in your console to see the documentation)
# Drop the SampleID column afterwards

data_predict <-
  data_tidy %>% 
  left_join(data_demo, by='SampleID') %>% 
  select(c(2:ncol(data_tidy)), "Asthma_Status") %>% 
  clean_names()

# Note that clean_names might change the column name of your prediction value of interest 
# Make sure you take a look at your resulting dataframe to see the new colume name

# Build a linear model. Remember that if you are doing a binomial regression you will need to convert your outcome to a binary output
# Multinomial regression will use the multinom() function from nnet package

data_predict <-
  data_predict %>% 
  mutate(asthma_status = ifelse(asthma_status == 'NAS',0,1))

# Shuffle the data set and split into a test and prediction data

shuffle_idx = sample(nrow(data_predict), nrow(data_predict))
data_predict_shuffled <- data_predict[shuffle_idx,]

split_idx <- sample(c(T,F), size=nrow(data_predict), replace=T, prob=c(0.80,0.20))

data_train <- data_predict_shuffled[split_idx,]
data_test <- data_predict_shuffled[!split_idx,]

#Generalized linear model (predicting on asthma status across all data points)

mdl <- glm(formula=asthma_status~., family=binomial, data=data_train)
summary(mdl)

# Predict on your test data and plot accuracy

logit_pred <- predict(mdl, data_test %>% select(-asthma_status), type='response')

# Convert your coefficients into a dataframe so you can plot feature importances
# Coefficients are the slopes of the log odds (the bigger the abs of the slope the bigger the influence of that gene/protein predictive capacity)

mdl_coeff <- as.data.frame(coefficients(mdl))

mdl_coeff <- abs(mdl_coeff) %>% 
  mutate(gene = rownames(.))

mdl_coeff <- as.numeric(as.character(mdl_coeff))

  clean_names() %>% 
  dplyr::filter(gene != "(Intercept)") %>% 
  arrange(desc(coefficients_mdl)) %>% 
  
  dplyr::slice(c(2:11)) %>% 


mdl_coeff= mdl_coeff[order(mdl_coeff$coefficients(mdl)),]

ggplot(mdl_coeff, aes(x=gene, y=coefficients_mdl)) + 
  geom_col()

# The last step is selecting the top 10 to graph
# Slice is selecting the first 10 rows

## Making a bar graph

# Hint: you will need to manipulate this dataframe, sort it, and take the first 10 rows to plot

# Next step, repeat but for a Random Forest classifier

rf <- randomForest(formula=??, data=??, importance=T)

rf_predict <- predict(??, ??)

rf_importances <- as.data.frame(importance(rf_model))

# Now you can plot feature importance after sorting and taking the top 10 features

# Lastly, use xgboost to do the same

xgb <- xgboost(data=??, label=??, nrounds=20)

xgb_predict(??, ??)

# Don't worry about feature importances for xgboost

