install.packages("tidyverse")
install.packages("ggplot")
install.packages("reshape2")
install.packages("Hmisc")
install.packages("superheat")

library(tidyverse)
library(ggplot2)
library(reshape2)
library(Hmisc)
library(superheat)

## Load Data

demo.data = read.csv("./3. Data Analysis/Data/Demographic_Data.csv")
proteomics.data = read.csv("./3. Data Analysis/Data/Processed_Proteomics_Data.csv")

table(demo.data$Subject_Number)
table(demo.data$PrevPost)

## Data Exploration 

describe(demo.data)

proto.long = proteomics.data %>% pivot_longer(cols=-c("Gene_Name"))
proto.long$PrevPost = gsub("_.*", "", proto.long$name)
ggplot(data=proto.long, aes(x=name, y=log10(value), col=PrevPost)) + geom_boxplot(outlier.shape = NA)

superheat(log(proteomics.data[,2:ncol(proteomics.data)]), pretty.order.rows = T, 
          bottom.label.text.angle = 90)

## Dataset Manipulation 

# Combine data
proteomics.data = proteomics.data %>% filter(!is.na(Gene_Name)) %>%   
  mutate(Gene_Name = gsub("-", "_", Gene_Name))

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
  pivot_longer(cols = -c(Subject_Number), names_to = "Gene", values_to = "Value_Pre")

# Do the same thing for post exposure data

long.data.post = combined.data %>% filter(PrevPost == "Post") %>% 
  dplyr::select(c("Subject_Number", proteomics.data$Gene_Name)) %>%
  pivot_longer(cols = -c(Subject_Number), names_to = "Gene", values_to = "Value_Post")


# Combine pre and post by subject number

long.data = merge(long.data.pre, long.data.post, by = c("Gene", "Subject_Number"))


##  Paired T-Test Example

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

##  Plotting - Boxplot

rownames(results_matrix_order[1,]);   # "TGM2"


# Now we can use that gene to plot the data for that specific gene

combined.data.tgm2 = combined.data[,c("Subject_Number", "PrevPost", "TGM2")];
combined.data.tgm2$PrevPost = as.factor(combined.data.tgm2$PrevPost);
combined.data.tgm2$TGM2 = as.numeric(combined.data.tgm2$TGM2);

ggplot(data=combined.data.tgm2, aes(x=PrevPost, y=TGM2, col=PrevPost)) + geom_boxplot();    ### If your data look weird, might be because values aren't numeric (or factor or whatever) when they should be



##  Plotting - Volcano

results_matrix$Gene = rownames(results_matrix)

# Here we are creating categories for what values are significantly up or down regulated, based on p value and fc

results_matrix$diffexp = "Not Significant"
results_matrix$diffexp[results_matrix$Log2FC > log2(1.5) & results_matrix$Pval < 0.05] = "Up Regulated, Significant"
results_matrix$diffexp[results_matrix$Log2FC < -log2(1.5) & results_matrix$Pval < 0.05] = "Down Regulated, Significant"



# ggplot is a great plotting package in r, that is super customizeable 
# Here, I am basically just saying put Log2FC on the x axis, use -log10(Pval) on the y-axis, and color it based on the diffexp value I just defined above

ggplot(data=results_matrix, aes(x=Log2FC, y=-log10(Pval), col=diffexp)) + xlim(c(-2,2)) +
  geom_point() + theme_minimal()  + scale_color_manual(values=c("dodgerblue", "grey", "firebrick1"), name="")


##  Further Analyses

long.data.pre = combined.data %>% filter(PrevPost == "Pre") %>% 
  dplyr::select(c("Subject_Number", proteomics.data$Gene_Name, "Responder_Status", "Sex", "Race", "Age", "BMI", "Asthma_Status")) %>%
  pivot_longer(cols = -c(Subject_Number, Responder_Status, Sex, Race, Age, BMI, Asthma_Status), names_to = "Gene", values_to = "Value_Pre") %>%
  arrange(Gene)
long.data.pre$Value_Pre = as.numeric(long.data$Value_Pre)


long.data.pre.loop = combined.data[which(combined.data$PrevPost=="Pre"),]

# Group 5: Asthma-associated proteins PRE
# At baseline (pre-exposure): What proteins differ by asthma status?

stat.test.asthma.pre = long.data.pre %>%
  group_by(Gene) %>%
  summarise(pval = t.test(Value_Pre ~ Asthma_Status)$p.value,
            as_tibble(as.list(t.test(Value_Pre ~ Asthma_Status)$estimate)),
            statistic = t.test(Value_Pre ~ Asthma_Status)$statistic) %>%
  mutate(pval.BH = p.adjust(pval, method="BH"))


##  Further Analyses

long.data.post = combined.data %>% filter(PrevPost == "Post") %>% 
  dplyr::select(c("Subject_Number", proteomics.data$Gene_Name, "Responder_Status", "Sex", "Race", "Age", "BMI", "Asthma_Status")) %>%
  pivot_longer(cols = -c(Subject_Number, Responder_Status, Sex, Race, Age, BMI, Asthma_Status), names_to = "Gene", values_to = "Value_Post") %>%
  arrange(Gene)
long.data.post$Value_Post = as.numeric(long.data$Value_Post)


long.data.post.loop = combined.data[which(combined.data$PrevPost=="Post"),]

# Group 5: Asthma-associated proteins
# Post-exposure: What proteins differ by asthma status?

stat.test.asthma.post = long.data.post %>%
  group_by(Gene) %>%
  summarise(pval = t.test(Value_Post ~ Asthma_Status)$p.value,
            as_tibble(as.list(t.test(Value_Post ~ Asthma_Status)$estimate)),
            statistic = t.test(Value_Post ~ Asthma_Status)$statistic) %>%
  mutate(pval.BH = p.adjust(pval, method="BH"))



