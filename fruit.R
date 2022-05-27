
#### SNP data and Nucleotide diversity analysis ####

### Import necessary module ####
#install packages required
#install.packages('tidyverse','emmeans','lme4')
library(tidyverse)
library(emmeans)
library(lme4)


### Working directory and reading data table ####
setwd('~/Desktop/Genomics/ag1000G3/All_country/R/')
##
div_tab <- read.table('div_tab.csv', sep = ",", header = T, quote = "\"",
                      stringsAsFactors = TRUE)
#View(div_tab)
Hap_tab <- read.table('Hap_stat.csv', sep = ",", header = T, quote = "\"",
                      stringsAsFactors = TRUE)
#View(Hap_tab)

#### basic data manipulation #####
div_tab$pops <- as.factor(div_tab$pops)
Hap_tab$country <- as.factor(Hap_tab$country)

### Statistical test - verify data distribution #### 
# Verifying  normality
hist(div_tab$pi, breaks = 'sturges')
qqnorm(div_tab$pi)
qqline(div_tab$pi)

hist(div_tab$tajimaD, breaks = 'sturges')
qqnorm(div_tab$tajimaD)
qqline(div_tab$tajimaD)

hist(Hap_tab$H, breaks = 'sturges')
qqnorm(Hap_tab$H)
qqline(Hap_tab$H)

# Conc: The data doesn't follow the normal distribution #

#### Nucleotide diversity statistics #####
#### Nucleotide diversity varies between mosquito populations ?
kruskal.test(pi~pops, data = div_tab)

## Calculate the median and the mean 
mean(div_tab$pi)
median(div_tab$pi)

# fit a model
fit1 <- glm(pi~pops, data = div_tab, family = quasipoisson())
deviance(fit1)/df.residual(fit1)
summary(fit1)
car::Anova(fit1, type='III')

# Multiple comparison 
nuc_div_emm <- emmeans(fit1, pairwise ~ pops, type="response", adjust="none")
nuc_div_contrast_df <- summary(nuc_div_emm$contrasts, infer=TRUE)%>%as.data.frame()
nuc_div_emm_df <- summary(nuc_div_emm$emmeans, infer=TRUE)%>%as.data.frame()


##### Tadjima statistics ####
# Nucleotide diversity varies between mosquito populations ?
kruskal.test(tajimaD~pops, data = div_tab)

## Calculate the median and the mean 
mean(div_tab$tajimaD, na.rm=TRUE)
median(div_tab$tajimaD, na.rm=TRUE)

####
fit2 <- glm(tajimaD~pops, data = div_tab, family = gaussian())
deviance(fit2)/df.residual(fit2)
summary(fit2)
car::Anova(fit2, type='III')

# Multiple comparison 
tajima_emm <- emmeans(fit2, pairwise ~ pops, type="response", adjust="none")
tajima_contrast_df <- summary(tajima_emm$contrasts, infer=TRUE)%>%as.data.frame()
tajima_emm_df <- summary(tajima_emm$emmeans, infer=TRUE)%>%as.data.frame()

# data summary as dataFrame
sum_tab <-  div_tab%>%
  dplyr::group_by(pops)%>%
  dplyr::summarise(me_pi = median(pi), 
                   min_pi = min(pi),
                   max_pi = max(pi),
                   me_D = median(tajimaD, na.rm = T),
                   min_D = min(tajimaD, na.rm = T),
                   max_D = max(tajimaD, na.rm = T),
                   me_theta = median(watt_theta), 
                   min_theta = min(watt_theta),
                   max_theta = max(watt_theta)
                   )
Hap_tab1 <- Hap_tab%>%
  dplyr::group_by(country)%>%
  dplyr::summarise(me_h = median(H),
                   min_h = min(H),
                   max_h = max(H))

max(sum_tab$max_pi)
min(sum_tab$me_pi)
median(sum_tab$me_pi)

max(sum_tab$max_D)
min(sum_tab$min_D)
max(sum_tab$me_D)
min(sum_tab$me_D)

max(div_tab$pi)
min(div_tab$pi)

max(sum_tab$max_theta)
min(sum_tab$min_theta)
max(sum_tab$me_theta)
min(sum_tab$me_theta)

max(div_tab$tajimaD, na.rm = T)
min(div_tab$tajimaD, na.rm = T)

max(Hap_tab$H, na.rm = T)
min(Hap_tab$H, na.rm = T)

max(Hap_tab1$me_h, na.rm = T)
min(Hap_tab1$me_h, na.rm = T)

# Save df to csv
write.csv(sum_tab, 'R_sumtab.csv')
write.csv(Hap_tab1, 'R_Hap_tab1.csv')

#### Data plotting ####
# Nucleotide diversity plot 
ggplot(data = div_tab, aes(x = pops, y = pi))+
  geom_boxplot(outlier.shape = NA, width=0.4, na.rm = T)+
  stat_summary(fun = mean, geom = "point",
               shape=20,size=1,color="red",fill='red')+
  theme_bw() +ylim(0,0.05)+
  labs(x ="", y = "Nucleotide diversity")+
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        panel.grid = element_blank(),
        legend.position = 'none')

# Tajima's D plot
ggplot(data = div_tab, aes(x = pops, y = tajimaD))+
  geom_boxplot(outlier.shape = NA, width=0.4, na.rm = T)+
  stat_summary(fun = mean, geom = "point",
               shape=20,size=1,color="red",fill='red')+
  geom_hline(yintercept = 0, size=.1)+ theme_bw() + ylim(-4,5)+
  labs(x ="", y = "Tajima's D")+
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        panel.grid = element_blank(),
        legend.position = 'none')

# Watterson theta plot
ggplot(data = div_tab, aes(x = pops, y = watt_theta))+
  geom_boxplot(outlier.shape = NA, width=0.4, na.rm = T)+
  stat_summary(fun = mean, geom = "point",
               shape=20,size=1,color="red",fill='red')+
  theme_bw() +
  labs(x ="", y = "Watterson’s theta")+
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        panel.grid = element_blank(),
        legend.position = 'none')

## Other plot : violon plot
ggplot(data = div_tab, aes(x = pops, y = pi))+
  geom_violin(width=0.5, na.rm = T)+
  stat_summary(fun = mean, geom = "point",
               shape=20,size=1,color="red",fill='red')+
  theme_bw() +
  labs(x ="", y = "Nucleotide diversity")+
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        panel.grid = element_blank(),
        legend.position = 'none')

ggplot(data = div_tab, aes(x = pops, y = tajimaD))+
  geom_violin(width=0.5, na.rm = T)+
  stat_summary(fun = mean, geom = "point",
               shape=20,size=1,color="red",fill='red')+
  theme_bw() +
  labs(x ="", y = "Tajima's D")+
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        panel.grid = element_blank(),
        legend.position = 'none')

### Haplotype diversity plot
ggplot(data = Hap_tab, aes(x = country, y = H))+
  geom_boxplot(outlier.shape = NA, width=0.4, na.rm = T)+
  stat_summary(fun = mean, geom = "point",
               shape=20,size=2,color="red",fill='red')+
  theme_bw() +
  labs(x ="", y = "Haplotype diversity")+
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        panel.grid = element_blank(), legend.position = 'none')

ggplot(data = Hap_tab, aes(x = country, y = H))+
  geom_violin(width=0.4, na.rm = T)+
  stat_summary(fun = mean, geom = "point",
               shape=20,size=2,color="red",fill='red')+
  theme_bw() +
  labs(x ="", y = "Haplotype diversity")+
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        panel.grid = element_blank(), legend.position = 'none')


