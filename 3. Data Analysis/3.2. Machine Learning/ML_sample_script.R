# Install libraries
# We will need need tidyverse, randomForest, janitor, and xgboost
# Remember you can pass a vector of library names to install multiple libraries at once
install.packages(??)


# Load libraries
library(??)
library(??)
...


# Load data

data_proteomics <- read.csv(??)
data_demo <- read.csv(??)

# Change to tidy format 

data_tidy <-
  data_proteomics %>% 
  pivot_longer(cols=??, names_to=??, values_to=??) %>% 
  picot_wider(names_from=??, values_from=??) 

# In the demographics dataframe, let's keep only our SampleID and variable to predict on

data_demo <-
  data_demo %>% 
  select(?)

# Join the two dataframes on SampleID, preserving all the rows in proteomics (what kind of join is this? Type ?join in your console to see the documentation)
# Drop the SampleID column afterwards

data_predict <-
  first_dataframte %>% 
  join_function(second_dataframe, by=??) %>% 
  select(??) %>% 
  clean_names()

# Note that clean_names might change the column name of your prediction value of interest 
# Make sure you take a look at your resulting dataframe to see the new colume name

# Build a linear model. Remember that if you are doing a binomial regression you will need to convert your outcome to a binary output
# Multinomial regression will use the multinom() function from nnet package

data_predict <-
  data_predict %>% 
  mutate(prediction = ifelse(??,??,??))

# Shuffle the data set and split into a test and prediction data

shuffle_idx = sample(??, ??)
data_predict_shuffled <- data_predict[??,]

split_idx <- sample((T,F), size=??, replace=T, prob=??)

data_train <- data_predict_shuffled[??,]
data_test <- data_predict_shuffled[??,]

mdl <- glm(formula=??, data=??)

# Predict on your test data and plot accuracy

logit_pred <- predict(??, ??)

# Convert your coefficients into a dataframe so you can plot feature importances

mdl_coeff <- as.data.frame(coefficients(??))

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

