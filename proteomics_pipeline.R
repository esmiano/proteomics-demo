# ==============================================================================
# PROTEOMICS PIPELINE
# ==============================================================================
# Date created: 07-OCT-2026
# 
# Description:
# A basic pipeline for processing and analysis of proteomics data.

# ==============================================================================
# PREPARATION
# ==============================================================================

# Remove artefacts
rm(list=ls())

# Load packages
library(tidyverse)
library(vroom)
library(Hmisc)

# Read in data
df_data <- vroom('')

# ==============================================================================
# MASTER DATAFRAME
# ==============================================================================

# Define new master protein dataframe
df_data <-
  # Take complete protein dataframe as input and pipe
  df_data %>%
  # Filter to only include Master Protein entries and exclude low FDR confidence entries
  filter(Master == "IsMasterProtein" & Protein.FDR.Confidence..Combined != "Low")

# ==============================================================================
# SAMPLE DATAFRAME
# ==============================================================================

# define new object for sample information
df_sample <- 
  # create dataframe with single column containing Raw Abundance column numbers
  # identify which columns, by number, have the pattern "Abundance..F" in their name
  data.frame(Raw.column.number = str_which(string = colnames(df_data),
                                           pattern = "Abundance..F"))

# add column with associated raw column names
df_sample$Raw.column.name <- 
  # use the raw column numbers previously identified to index the column names
  colnames(df_data)[df_sample$Raw.column.number]

# add new column with column numbers for normalised values
df_sample$Normalised.column.number <- 
  # identify which columns, by number, have the pattern "Abundances..N" in their name
  str_which(string = colnames(df_data), pattern = "Abundances..N")

# add column with associated normalised column names
df_sample$Normalised.column.name <- 
  # use the normalised column numbers previously identified to index the column names
  colnames(df_data)[df_sample$Normalised.column.number]

# add column containing replicate names
df_sample$Replicate <-
  # use string extract to match either "F1", "F2", or "F3" in "Raw.column.name" column
  # replicates are the same for raw and normalised values
  str_extract(string = df_sample$Raw.column.name,
              pattern = "(F1|F2|F3)")

# add column containing condition names
df_sample$Condition <- 
  # use string extract to match either "Control.b" or "Exp.B" in "Raw.column.name" column
  str_extract(string=df_sample$Raw.column.name,
              pattern= "Control.b|Exp.B")

# add column containing concatenated Condition and Replicate information
df_sample$SampleID <- 
  # use paste0 to concatenate Condition and Replicate column values together
  paste0(df_sample$Condition, sep="_", df_sample$Replicate)

# ==============================================================================
# TASK 4: LOG2 TRANSFORMATION
# ==============================================================================

# define new log2 transformed dataframe
df_log2_data <- 
  # take master protein dataframe as input
  df_data %>%
  # log2 transform all values in the normalised and raw columns
  mutate(across(df_sample$Normalised.column.number, log2),
         across(df_sample$Raw.column.number, log2))

# ==============================================================================
# TASK 5: PCA
# ==============================================================================

# copy dataframe before modifying
df_pca.in <-
  df_log2_data

# filter columns to just the now-log2 transformed column numbers
df_pca.in <-
  df_pca.in[,df_sample$Normalised.column.number]

# perform PCA analysis
df_pca.out <- 
  prcomp(t(df_pca.in))

# index out of the full output using $
df_pca.out <- 
  df_pca.out$x

# merge the data with df_sample to pull in individual sample info
df_pca.out <-
  merge(x = df_sample,
        by.x = "Normalised.column.name",
        y = df_pca.out,
        by.y = 0)

# define pca results object (without converting to dataframe)
pca_res <- 
  prcomp(t(df_pca.in))

# define variance explained object
var_exp <-
  # calculate the variance explained by each principle component as a percentage
  # round to 2 decimal places
  round(((pca_res$sdev^2 / sum(pca_res$sdev^2)) * 100), 2)

# plotting PC1 & PC2
ggplot(data = df_pca.out,
       # select x and y data
       aes(x = PC1, y = PC2,
           # define point shape by condition
           shape = Condition,
           # define point shape by replicate
           colour = Replicate)) +
  # add points and set size to 2.5
  geom_point(size = 2.5) +
  # add horizontal line at y intercept 0
  geom_hline(yintercept = 0) +
  # add vertical line at x intercept 0
  geom_vline(xintercept = 0) +
  # add title
  ggtitle("Principal Component Analysis: PC1 & PC2") +
  # change legend label names
  scale_shape_discrete(labels = c("Control B", "Experiment B")) +
  # add variance explained percentage to x axis
  xlab(paste0("PC1 (", var_exp[1], "%)")) +
  # add variance explained percentage to y axis
  ylab(paste0("PC2 (", var_exp[2], "%)")) +
  # set theme to classic (blank background) for clarity
  theme_classic()

# plotting PC3 & PC4
ggplot(data = df_pca.out,
       # select x and y data
       aes(x = PC3, y = PC4,
           # define point shape by condition
           shape = Condition,
           # define point shape by replicate
           colour = Replicate)) +
  # add points and set size to 2.5
  geom_point(size = 2.5) +
  # add horizontal line at y intercept 0
  geom_hline(yintercept = 0) +
  # add vertical line at x intercept 0
  geom_vline(xintercept = 0) +
  # add title
  ggtitle("Principal Component Analysis: PC3 & PC4") +
  # change legend label names
  scale_shape_discrete(labels = c("Control B", "Experiment B")) +
  # add variance explained percentage to x axis
  xlab(paste0("PC3 (", var_exp[3], "%)")) +
  # add variance explained percentage to y axis
  ylab(paste0("PC4 (", var_exp[4], "%)")) +
  # set theme to classic (blank background) for clarity
  theme_classic()

# ==============================================================================
# TASK 6: T-TEST
# ==============================================================================

# add t-test p-value column to log2 dataframe
df_log2_data <- 
  # take log2 dataframe as input
  df_log2_data %>% 
  # perform for each row
  rowwise() %>% 
  mutate(
    Ttest_ConB.v.ExpB_p.value = t.test(x = c(Abundances..Normalized...F1..127..Sample..Exp.B,
                                             Abundances..Normalized...F2..127..Sample..Exp.B, 
                                             Abundances..Normalized...F3..127..Sample..Exp.B),
                                       y = c(Abundances..Normalized...F1..126..Sample..Control.b,
                                             Abundances..Normalized...F2..126..Sample..Control.b,
                                             Abundances..Normalized...F3..126..Sample..Control.b),
                                       paired = T)$p.value)

# add t-test LogFC column to log2 dataframe
df_log2_data <- 
  # take log2 dataframe as input
  df_log2_data %>%
  # perform for each row
  rowwise() %>% 
  mutate(
    Ttest_ConB.v.ExpB_LogFC = t.test(x = c(Abundances..Normalized...F1..127..Sample..Exp.B,
                                           Abundances..Normalized...F2..127..Sample..Exp.B, 
                                           Abundances..Normalized...F3..127..Sample..Exp.B),
                                     y = c(Abundances..Normalized...F1..126..Sample..Control.b,
                                           Abundances..Normalized...F2..126..Sample..Control.b,
                                           Abundances..Normalized...F3..126..Sample..Control.b),
                                     paired = T)$estimate)

# task 6 output
median(df_log2_data$Ttest_ConB.v.ExpB_p.value, na.rm=T)
median(df_log2_data$Ttest_ConB.v.ExpB_LogFC, na.rm=T)

# ------------------------------------------------------------------------------
# FDR ADJUSTMENT
# ------------------------------------------------------------------------------

# add new column for fdr-adjusted p-values
df_log2_data$FDR_adj.p.values <- p.adjust(p = df_log2_data$Ttest_ConB.v.ExpB_p.value,
                                          method = "fdr")

# task 7 output
median(df_log2_data$FDR_adj.p.values)

# ------------------------------------------------------------------------------
# VOLCANO PLOT
# ------------------------------------------------------------------------------

# call ggplot function and take log2 dataset as input
ggplot(data = df_log2_data,
       # select t-test logfc values for the x values
       aes(x = Ttest_ConB.v.ExpB_LogFC,
           # select -log10 t-test p-values for the y values
           y = -log10(Ttest_ConB.v.ExpB_p.value),
           # define point shape by FDR-adjusted p-value significance
           shape = FDR_adj.p.values < 0.05)) +
  # add points
  geom_point() +
  # add title
  ggtitle("T-Test Results for Experiment B vs Control B") +
  # add x-axis label
  xlab("LogFC") +
  # add y-axis label, making p italic
  ylab(expression(paste("- log"[10]*italic(" p"), "-value"))) +
  # rename figure legend
  scale_shape_discrete(name = expression(paste("FDR-adjusted ", italic("p"), "-values")),
                       # change label names, false = not siginificant (p>0.05)
                       labels = c(expression(paste("Not significant (", italic("p")>0.05, ")")),
                                  # true = significant (p<0.05)
                                  expression(paste("Significant (", italic("p")<0.05, ")")))) +
  # set theme to classic (blank background) for clarity
  theme_classic()
