library(ggplot2)
library(dplyr)
library(viridis)
library(ggpointdensity)
library(raster)
theme_set(theme_bw())


#two populations around two foci (settlements)

df1popxmean = 6
df1popymean = 3

#pop1
#dispersed subpop (higher sd)
df1dispersed <-  tibble(x = rnorm(50, mean=pop1xmean,sd = 1),
              y = rnorm(50, mean=pop1ymean,sd = 1))
#our concentrated subpop pop 
df1conc <- tibble(x = rnorm(50, mean=pop1xmean,sd = 0.5),
                  y = rnorm(50, mean=pop1ymean,sd = 0.5))
#combine into one pop
df1 <- rbind(df1conc, df1dispersed)

#add popinfo label
df1$pop <- paste0("pop1")
#add focal point distance
df1$focalx = pop1xmean
df1$focaly = pop1ymean


df2popxmean = 9
df2popymean = 3

df2 <- tibble(x = rnorm(250, mean=df2popxmean,sd = 0.5),
              y = rnorm(250, mean=df2popymean,sd = 0.5))
df2$pop <- paste0("pop2")


df2$focalx = df2popxmean
df2$focaly = df2popymean

c <- rbind(df1, df2)

for(i in 1:nrow(c)) {
  s<-raster::pointDistance(c(c$x[i], c$y[i]), c(c$focalx[i], c$focaly[i]), lonlat = F)
  c[i, 6] <- s
}

c$dist <- 1- c$...6
c$ltwodist <- log2(c$dist)


c %>% 
  ggplot( mapping = aes(x = x, y = y, colour=pop, size=dist)) +
  geom_point(alpha = 0.5) +
  xlim(2,11) +
  xlab("Easting") +
  ylab("Northing") +
  annotate(geom = "text", x=df1popxmean, y=df1popymean, label="Site 1", color = 'black') +
  annotate(geom = "text", x=df2popxmean, y=df2popymean, label="Site 2", color = 'black')











