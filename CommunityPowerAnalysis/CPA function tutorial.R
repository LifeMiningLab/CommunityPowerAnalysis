library(tidyverse)
library(viridis)
library(dplyr)
library("ggsci")
library("ggplot2")
library("gridExtra")
library(gtools)
library(data.table)
library(ggpubr)

result_dt <- KeggResultProcess(inputDir = '/Users/jwbaek9506/Documents/antibiotic project paper/result_bak', filePattern = '.fa.faa.annotation.table.txt')


# Example usage
# 
communityPower <- UniqueKeggCPA(result_dt)

ggboxplot(communityPower, x="Community",y="UniqueKeggCount", fill = "Community", legend = 'none', outlier.size=0.6, size=0.1, ylab = 'Number of uniqkegg') +
  rotate_x_text(angle=45) +
  ggtitle("Number of uniqkegg") +
  rremove('x.text') +
  rremove('xlab') +
  # font("xlab", size = 6)+
  font("ylab", size = 6)+
  font("y.text", size=6) +
  # font("xy.text", size = 6)+
  font("legend.title", size = 6)+
  font("legend.text",size = 6) +
  font("title", size = 6, face = "bold") +
  theme(legend.key.size = unit(0.3, 'cm')) +
  guides(color = guide_legend(title.position = "top",
                              title.hjust = 0.5))


# 
CPAResult <- CommunityPowerAnalysis(result_dt)

test <- CPAResult %>% select(Community, Sample, CP_Score) %>% unique()

ggboxplot(test, x="Community",y="CP_Score", fill = "Community", legend = 'none', outlier.size=0.6, size=0.1, ylab = 'Community power score') +
  rotate_x_text(angle=45) +
  ggtitle("Community power analysis") +
  rremove('x.text') +
  rremove('xlab') +
  # font("xlab", size = 6)+
  font("ylab", size = 6)+
  font("y.text", size=6) +
  # font("xy.text", size = 6)+
  font("legend.title", size = 6)+
  font("legend.text",size = 6) +
  font("title", size = 6, face = "bold") +
  theme(legend.key.size = unit(0.3, 'cm')) +
  guides(color = guide_legend(title.position = "top",
                              title.hjust = 0.5))
