##################################################################################################
##################################################################################################
##   PREHEAT Retreat - 2 Group Analyses - CEYM - draft
##
##   Kyle Roell
##   October 20-21
##
##   These analyses will detail how to do (simple) 2-group analyses in R using our
##   sample proteomics dataset.
##
##################################################################################################
##################################################################################################

##############################
##############################
##
##  Install and Load Packages
##
##############################
##############################
install.packages("tidyverse")
install.packages("ggplot2")
install.packages("reshape2") 
install.packages("Hmisc") 
install.packages("superheat") 
install.packages("gtsummary")


library(tidyverse)
library(ggplot2)
library(reshape2)
library(Hmisc)
library(superheat)
library(gtsummary)

##############################
##############################
##
##  Load Data
##
##############################
##############################
# Demographic Data
demo.data = read.csv("/cloud/project/3. Data Analysis/Data/Demographic_Data.csv")

# Proteomics Data
proteomics.data = read.csv("/cloud/project/3. Data Analysis/Data/Processed_Proteomics_Data.csv")

# Make sure everyone has two samples
table(demo.data$Subject_Number)
table(demo.data$PrevPost)


##############################
##############################
##
##  Data Exploration
##
##############################
##############################
### Let's just take a look at our data before doing anything analyses, maybe should do some of this before manipulation
# I like to use description functions like describe from Hmisc package, others exist
describe(demo.data)


### Can plot big data, that you can't just visually inspect manually
# I want to plot boxplots of all genes per sample, same as what we did in pre-processing script
# Need to reformat data to tidy format ;) 
proto.long = proteomics.data %>% pivot_longer(cols=-c("Gene_Name"))
proto.long$PrevPost = gsub("_.*", "", proto.long$name)
ggplot(data=proto.long, aes(x=name, y=log10(value), col=PrevPost)) + geom_boxplot(outlier.shape = NA)


##############################
##############################
##
##  Dataset Manipulation
##
##############################
##############################
# Combine data
proteomics.data = proteomics.data %>% filter(!is.na(Gene_Name)) %>%   # filter let's us filter by some condition, here Gene_Name that is NOT NA (aka not missing)
  mutate(Gene_Name = gsub("-", "_", Gene_Name))  # mutate allows us to create a new variable or modify existing one, here I'm just changing - to _ in gene names

# First transpose the data so subjects are on rows
proto.tr.data = data.frame(t(proteomics.data)) # "t" transposes dataset (columns to rows, rows to columns)
colnames(proto.tr.data) = proto.tr.data[1,]
proto.tr.data = proto.tr.data[-1,]

# Create ID column
proto.tr.data$SampleID = rownames(proto.tr.data)


# Merge two datasets so everything is together
combined.data = merge(demo.data, proto.tr.data, by = "SampleID")
combined.data = data.frame(combined.data)



######### Long vs Wide format
# Here we want to put our data into long format, once for Pre exposure and once for Post exposure
# We then can merge those two datasets into one by Subject ID

long.data.pre = combined.data %>% filter(PrevPost == "Pre") %>%     # filter to only "Pre" exposure
  dplyr::select(c("Subject_Number", proteomics.data$Gene_Name)) %>%        # select alows you to just use the specified columns
  pivot_longer(cols = -c(Subject_Number), names_to = "Gene", values_to = "Value_Pre") # pivot longer turns the data into long form
                                                                                       # cols says which columns to pivot longer -- i am telling it to use everything EXCEPT Subject_Number
                                                                                       # names_to  says where to store info from column names for data specified by "cols"
                                                                                       # values_to  says where to store the actual data

### So basically to create this long dataset, I am saying:
# 1. Take the combined dataset (both demographic and proteomics datasets) and filter for only those with "Pre" for PrevPost
# 2. Next, select ONLY the columns corresponding to "Subject_Name" or any gene name
# 3. Finally, create a new dataset where I want to just use the data in the gene columns, by Subject ID, and create a new column "Gene" to tell me which gene name that value corresponds to. 



# Do the same thing for post exposure data
long.data.post = combined.data %>% filter(PrevPost == "Post") %>% 
  dplyr::select(c("Subject_Number", proteomics.data$Gene_Name)) %>%
  pivot_longer(cols = -c(Subject_Number), names_to = "Gene", values_to = "Value_Post")

# Combine pre and post by subject number
long.data = merge(long.data.pre, long.data.post, by = c("Gene", "Subject_Number"))



##############################
##############################
##
##  Paired T-Test Example
##
##############################
##############################
# Now we can set this through to the t test function
# Make sure value are numberic
long.data$Value_Pre = as.numeric(long.data$Value_Pre)
long.data$Value_Post = as.numeric(long.data$Value_Post)

# Here I am running a BUNCH of t-tests, 1 per Gene 
stat.test = long.data %>%
  group_by(Gene) %>%     # group_by Gene, so run through each gene 1 by 1
  summarise(pval = t.test(Value_Pre, Value_Post, paired=T)$p.value,  # extract p-value from t-test
            estimate = t.test(Value_Pre, Value_Post, paired=T)$estimate,  # extract estimates (mean difference between groups) from t-test
            statistic = t.test(Value_Pre, Value_Post, paired=T)$statistic) %>%  # extract test statistic from t-test
  mutate(pval.BH = p.adjust(pval, method="BH"))   # multiple test corrections
 


# This is doing the same thing as above, but in a loop (which may be more intuitive to understand but, generally, much slower run time)
results_matrix = matrix(NA, nrow=nrow(proteomics.data), ncol=7)
colnames(results_matrix) = c("Estimate", "Statistic", "Pval", "Pval.BH", "Mean.Pre", "Mean.Post", "Log2FC")
rownames(results_matrix) = proteomics.data$Gene_Name

for(i in proteomics.data$Gene_Name) {

  pre.data = as.numeric(unlist(combined.data %>% filter(PrevPost == "Pre") %>% arrange(Subject_Number) %>% dplyr::select(i)))
  post.data = as.numeric(unlist(combined.data %>% filter(PrevPost == "Post") %>% arrange(Subject_Number) %>% dplyr::select(i)))
  
  t.test.results = t.test(pre.data, post.data, paired=T)
  
  results_matrix[i, "Estimate"] = t.test.results$estimate[1]
  results_matrix[i, "Statistic"] = t.test.results$statistic
  results_matrix[i, "Pval"] = t.test.results$p.value

  results_matrix[i, "Mean.Pre"] = mean(pre.data)
  results_matrix[i, "Mean.Post"] = mean(post.data)
  results_matrix[i, "Log2FC"] = log2(mean(post.data)/mean(pre.data))
  
}

results_matrix[,"Pval.BH"] = p.adjust(results_matrix[,"Pval"], method="BH")

# Convert to dataframe
results_matrix = data.frame(results_matrix)


# We can easily sort these lists from small to largest p-value using order() 
results_matrix_order = results_matrix[order(results_matrix$Pval),]

head(results_matrix)
head(results_matrix_order)



##############################
##############################
##
##  Plotting - Boxplot
##
##############################
##############################
# A boxplot allows us to see differences across groups
# Generally, the x-axis is used to display groups, and the y-axis is used to display actual continous values

# Let's use the top gene from the paired t-test analysis
# Get our ordered matrix, and select first gene (genes are rownames)
rownames(results_matrix_order[1,]);   # "TGM2"


# Now we can use that gene to plot the data for that specific gene
combined.data.tgm2 = combined.data[,c("Subject_Number", "PrevPost", "TGM2")];
combined.data.tgm2$PrevPost = as.factor(combined.data.tgm2$PrevPost);
combined.data.tgm2$TGM2 = as.numeric(combined.data.tgm2$TGM2);

ggplot(data=combined.data.tgm2, aes(x=PrevPost, y=TGM2, col=PrevPost)) + geom_boxplot();    ### If your data look weird, might be because values aren't numeric (or factor or whatever) when they should be



##############################
##############################
##
##  Plotting - Volcano
##
##############################
##############################
# A volcano plot shows how p values relate to fold change (fc)

# Fold change is generally reported as log2(fc), where fold change describes the change between pre and post (in this study) values
# Here we are calculating log2 fold change as: log2(mean(post.data)/mean(pre.data))
# We take log2 because this allows us to say, log2fc=1 corresponds to the meaning the post value is double the pre value

# Set rownames
results_matrix$Gene = rownames(results_matrix)


# Here we are creating categories for what values are significantly up or down regulated, based on p value and fc

results_matrix$diffexp = "Not Significant"
results_matrix$diffexp[results_matrix$Log2FC > log2(1.5) & results_matrix$Pval < 0.05] = "Up Regulated, Significant"
results_matrix$diffexp[results_matrix$Log2FC < -log2(1.5) & results_matrix$Pval < 0.05] = "Down Regulated, Significant"


# ggplot is a great plotting package in r, that is super customizeable 
# Here, I am basically just saying put Log2FC on the x axis, use -log10(Pval) on the y-axis, and color it based on the diffexp value I just defined above
ggplot(data=results_matrix, aes(x=Log2FC, y=-log10(Pval), col=diffexp)) + xlim(c(-2,2)) +
  geom_point() + theme_minimal()  + scale_color_manual(values=c("dodgerblue", "grey", "firebrick1"), name="")



##############################
##############################
##
##  Further Analyses
##
##############################
##############################
#### These are all samples of what you will be doing tomorrow, I have done this here for PRE exposure, and you will do it for POST exposure.
# Please use pieces from this code, or code anywhere in the document, if necessary and modify to suit your analysis.
# And of course, you are more than welcome to write your own code to do these analyses.


# For all these questions we are looking at PRE exposure, so let's just create a long dataset with only pre exposure across the various columns
long.data.pre = combined.data %>% filter(PrevPost == "Pre") %>% 
  dplyr::select(c("Subject_Number", proteomics.data$Gene_Name, "Responder_Status", "Sex", "Race", "Age", "BMI", "Asthma_Status")) %>%
  pivot_longer(cols = -c(Subject_Number, Responder_Status, Sex, Race, Age, BMI, Asthma_Status), names_to = "Gene", values_to = "Value_Pre") %>%
  arrange(Gene)
long.data.pre$Value_Pre = as.numeric(long.data$Value_Pre)

long.data.pre.loop = combined.data[which(combined.data$PrevPost=="Pre"),]



# Prepare data set for POST exposure analysis, create long dataset with only post exposure across the various columns 
long.data.post = combined.data %>% filter(PrevPost == "Post") %>% 
  dplyr::select(c("Subject_Number", proteomics.data$Gene_Name, "Responder_Status", "Sex", "Race", "Age", "BMI", "Asthma_Status")) %>%
  pivot_longer(cols = -c(Subject_Number, Responder_Status, Sex, Race, Age, BMI, Asthma_Status), names_to = "Gene", values_to = "Value_Post") %>%
  arrange(Gene)
long.data.post$Value_Post = as.numeric(long.data$Value_Post)
long.data.post.loop = combined.data[which(combined.data$PrevPost=="Post"),]

##############################################
##############################################
## 
##  4. BMI-associated differences
##
##############################################
##############################################

### Group 4: BMI-associated differences ### 
### At baseline (pre-exposure): What proteins differ by BMI? ### 
stat.test.BMI.pre = long.data.pre %>%
  group_by(Gene) %>%
  summarise(as_tibble(as.list(summary(lm(Value_Pre ~ BMI))$coefficients["BMI",]))) %>%
  rename(pval=`Pr(>|t|)`) %>%
  mutate(pval.BH = p.adjust(pval, method="BH"))

stat.test.BMI.pre.significant <- subset(stat.test.BMI.pre, pval < 0.05)
# We found that 14 genes were statistically significant at the 0.05 alpha level. 

stat.test.BMI.pre.significant_adjusted <- subset(stat.test.BMI.pre, pval.BH < 0.05)
# When we ran an Benjamini & Hochberg adjusted t test, which controlled for 
# false positives or Type 1 errors, we did not find that any of the 14 originally 
# identified genes or any other ones were statistically significant. 


####  Post-exposure: What proteins differ by BMI? ### 
stat.test.BMI.post = long.data.post %>%
  group_by(Gene) %>%
  summarise(as_tibble(as.list(summary(lm(Value_Post ~ BMI))$coefficients["BMI",]))) %>%
  rename(pval=`Pr(>|t|)`) %>%
  mutate(pval.BH = p.adjust(pval, method="BH"))

stat.test.BMI.post.significant <- subset(stat.test.BMI.post, pval < 0.05)
# We found that 302 genes were statistically significant at the 0.05 alpha level. 

stat.test.BMI.post.significant_adjusted <- subset(stat.test.BMI.post, pval.BH < 0.05)
# When we ran an Benjamini & Hochberg adjusted t test, which controlled for 
# false positives or Type 1 errors, we did not find that any of the 302 originally 
# identified genes or any other ones were statistically significant. 


### Compare these two findings: are they the same? ### 
stat.test.BMI.similarities <- merge(stat.test.BMI.post.significant, stat.test.BMI.pre.significant, by = "Gene") 
# 3 of the genes are the same  
# IGHV5_10_1 
# IGLC7 
# MUC13 

stat.test.BMI.differences <- setdiff(stat.test.BMI.post.significant$Gene, stat.test.BMI.pre.significant$Gene)
# 3 genes difference confirmed 
view(stat.test.BMI.differences)

check <- stat.test.BMI.post.significant %>% filter(Gene %in% c("IGHV5_10_1", "IGLC7", "MUC13"))
# confirmed 

unique(check$Gene) # unique names of distinct variables (excluding duplicates)
n_distinct(check$Gene) # number of distinct variables (excluding duplicates)
nrow(check) # number of rows 


# Review demographics associated with the 3 genes 
stat.test.demographics.similarities <- merge(stat.test.BMI.similarities, long.data.post, by = "Gene")
# stat.test.demographics.similarities2 <- stat.test.demographics.similarities %>% filter(Gene %in% c("IGHV5_10_1", "IGLC7", "MUC13"))
# already done 
stat.test.demographics.similarities %>%  
  select(Responder_Status, Sex, Asthma_Status, BMI) %>% 
  tbl_summary()

stat.test.demographics.similarities %>%  
  select(Responder_Status, Sex, Asthma_Status, BMI, Gene) %>% 
  tbl_summary(by = "Gene")
# Same 3 people 


table(stat.test.demographics.similarities$BMI, stat.test.demographics.similarities$Responder_Status)

table(stat.test.demographics.similarities$Asthma_Status, stat.test.demographics.similarities$Responder_Status)

