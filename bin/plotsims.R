#######
#plotsims
#plotting simulation output
#tristan dennis 08.03.21 tristanpwdennis@gmail.com
######

library(tidyverse)

setwd('OneDrive - University of Glasgow/MOVE/manuscripts/move_review_2021/move_r/')

ibddata <- read.csv('metadata/isolation_by_distance.txt')
metadata <- read.delim('metadata/spatial_sim_individuals.txt')

#label group data based on generations ago born
metadata <- metadata %>% mutate(group = case_when(birth_time_ago < 10 ~ "alive",
                                                  birth_time_ago > 20 ~ "ancient"))
#remove nonsensical infinity values
ibddata <- ibddata %>% filter(geog_dist != 'Inf') 

#add metadata to each individual in the divergence data
ibddata <- left_join(ibddata, metadata, by = c('ind1' = 'vcf_label')) %>%  #join metadata to tskit data for ind1
  left_join(., metadata, by = c('ind2' = 'vcf_label')) 



#plot geography
metadata %>% ggplot(aes(x=x,y=y)) +
  geom_point()+
  facet_wrap(~group) +
  theme_minimal()

ibddata %>% ggplot(aes(x=geog_dist, y=ind_div, colour=group.x))+
  geom_point()
  

