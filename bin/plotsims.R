#######
#plotsims
#plotting simulation output
#I hate plotting in Python so I import the data into R for final plotting
#tristan dennis 08.03.21 tristanpwdennis@gmail.com
######

library(tidyverse)
library(cowplot)
library(MASS)
library(viridis)

setwd('/Users/tristanpwdennis/OneDrive - University of Glasgow/MOVE/manuscripts/move_review_2021/move_r/')

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}


#import data
metadata_terrible30 = read.csv("simoutput/everythingisterribleafter30.csv")
geo_terrible30 = read.csv("simoutput/everythingisterribleafter30_geodf.csv")
metadata_onefocusat30 = read.csv("simoutput/onefocusat30.csv")
geo_onefocusat30 = read.csv("simoutput/onefocusat30_geodf.csv")
metadata_twofociatthirty = read.csv("simoutput/twofocusat30.csv")
geo_twofociatthirty = read.csv("simoutput/twofociafter30_geodf.csv")


#add condition labels and bind metadata together
metadata_terrible30 = metadata_terrible30 %>% mutate(condition = 'thinning')
metadata_onefocusat30 = metadata_onefocusat30 %>% mutate(condition = 'singlefocus')
metadata_twofociatthirty = metadata_twofociatthirty %>% mutate(condition = 'twofoci')

metafull = rbind(metadata_terrible30, metadata_onefocusat30, metadata_twofociatthirty)
#repl silly time labels
metafull$time <- gsub('then', 'past', metafull$time)
metafull$time <- gsub('now', 'present', metafull$time)
#get local density
metafull$density <- get_density(metafull$x, metafull$y, n = 1000)


pretreat = metafull %>% filter(condition == 'twofoci' & time == 'past') %>% 
  ggplot(aes(x=x, y=y, colour=density))+
  geom_density_2d() +
  theme_minimal()+
  theme(legend.position = "none")+
  theme(axis.title = element_blank())+
  labs(x="Eastings", y="Northings")


thinning = metafull %>% filter(time == 'present' & condition == 'thinning') %>% 
  ggplot(aes(x=x, y=y, colour=density))+
  geom_density_2d() +
  theme_minimal()+
  scale_color_viridis() +
  theme(legend.position = "none")+
  theme(axis.title = element_blank())+
  labs(x="Eastings", y="Northings")

single_locus = metafull %>% filter(time == 'present' & condition == 'singlefocus') %>% 
  ggplot(aes(x=x, y=y, colour=density))+
  geom_density_2d() +
  theme_minimal()+
  scale_color_viridis() +
  theme(legend.position = "none")+
  theme(axis.title = element_blank())+
  labs(x="Eastings", y="Northings")

twofoci = metafull %>% filter(time == 'present' & condition == 'twofoci') %>% 
  ggplot(aes(x=x, y=y, colour=density))+
  geom_density_2d() +
  theme_minimal()+
  scale_color_viridis() +
  theme(legend.position = "none")+
  theme(axis.title = element_blank())+
  labs(x="Eastings", y="Northings")


plot_grid(pretreat, thinning, single_locus, twofoci, labels = c('A', 'B', 'C', 'D'), align = 'v', ncol =2)

ptpretreat = metafull %>% filter(condition == 'twofoci' & time == 'past') %>% 
  ggplot(aes(x=x, y=y, colour=density))+
  geom_density_2d(colour='grey', size=1) +
  scale_colour_viridis() +
  geom_point(size=3, alpha=0.5)+
  theme_minimal()+
  theme(legend.position = "none")+
  theme(axis.title = element_blank())+
  labs(x="Eastings", y="Northings")

ptthinning = metafull %>% filter(time == 'present' & condition == 'thinning') %>% 
  ggplot(aes(x=x, y=y, colour=density))+
  geom_density_2d(colour='grey', size=1) +
  scale_colour_viridis() +
  geom_point(size=3, alpha=0.5)+
  theme_minimal()+
  theme(legend.position = "none")+
  theme(axis.title = element_blank())+
  labs(x="Eastings", y="Northings")

ptsingle_locus = metafull %>% filter(time == 'present' & condition == 'singlefocus') %>% 
  ggplot(aes(x=x, y=y, colour=density))+
  geom_density_2d(colour='grey', size=1) +
  scale_colour_viridis() +
  geom_point(size=3, alpha=0.5)+
  theme_minimal()+
  theme(legend.position = "none")+
  theme(axis.title = element_blank())+
  labs(x="Eastings", y="Northings")

pttwofoci = metafull %>% filter(time == 'present' & condition == 'twofoci') %>% 
  ggplot(aes(x=x, y=y, colour=density))+
  geom_density_2d(colour='grey', size=1) +
  scale_colour_viridis() +
  geom_point(size=3, alpha=0.5)+
  theme_minimal()+
  theme(legend.position = "none")+
  theme(axis.title = element_blank())+
  labs(x="Eastings", y="Northings")

s = plot_grid(ptpretreat, ptthinning, ptsingle_locus, pttwofoci, labels = c('A', 'B', 'C', 'D'), align = 'v', ncol =2)


ggsave(file="test.svg", plot=s, width=10, height=8)
