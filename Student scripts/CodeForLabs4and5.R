
library(tmap)
library(raster)
library(gstat)
library(sf)
library(tidyverse)
library(tigris)
library(plotly)
library(terra)
library(spdep)
library(sfdep)
rm(list=ls())

#--------------------------------------------------
# Geographic borders
#--------------------------------------------------
state.border.LonLat  <- states(cb = TRUE)               %>% 
                           select(NAME, geometry)       %>%
                           filter(NAME == "California")  
  
county.border.LonLat  <- counties("CA", cb = TRUE) %>% select(GEOID, geometry)  


#--------------------------------------------------
# SOVI Data
#--------------------------------------------------
county.sovi <- read.csv("SOVI2020_CA_County.csv") %>%
  mutate(FIPS = str_pad(STCNTY, width = 5, pad = "0")) %>%
  rename(GEOID = FIPS) %>%
  mutate_all( ~ ifelse(. == -999, NA, .)) %>%
  select(-matches("M_|MP_|SPL_|EPL_|E_|F_")) 


# Make a version for the AEA Map projection
state.border.AEA   <- st_transform(state.border.LonLat ,3310)
county.border.AEA  <- st_transform(county.border.LonLat,3310)


#--------------------------------------------------
# Ozone Data
#--------------------------------------------------
point.ozone        <- read.csv("./Data/CA_OzonePopulation.csv")
point.ozone <- point.ozone[-1]
point.ozone.LonLat <- st_as_sf(point.ozone, coords=c("LONGITUDE","LATITUDE"),crs=4326)
point.ozone.AEA    <- st_transform(point.ozone.LonLat,3310)

#--------------------------------------------------
# Merge them all together
# This is hacked together because I couldn't get my head
# around the tidyverse version. 
#--------------------------------------------------
point.ozone.AEA <- st_join(point.ozone.AEA,county.border.AEA)

stations_in_county_split <- split(point.ozone.AEA$OZONE_1000PPB,point.ozone.AEA$GEOID)

Newdata <- data.frame(GEOID = names(stations_in_county_split),
                      OZONE_1000PPB = unlist(lapply(stations_in_county_split,mean,na.rm=TRUE)))


county.AllData        <- merge(county.sovi,Newdata,by="GEOID",all.x=TRUE,all.y=FALSE)
county.AllData.sf     <- full_join(county.border.AEA,county.AllData,by="GEOID")
county.AllData.LonLat <- st_transform(county.AllData.sf,4326)
county.AllData.AEA    <- st_transform(county.AllData.sf,3310)
names(county.AllData.AEA)


file.remove("./Data/county.AllData.AEA.geojson")
st_write(county.AllData.AEA, "./Data/county.AllData.AEA.geojson", driver = "GeoJSON")

point.ozone.AEA <- merge(point.ozone.AEA,county.sovi[,c(2,5,6)],by="GEOID",all.x=TRUE,all.y=FALSE)
file.remove("./Data/point.ozone.AEA.geojson")
st_write(point.ozone.AEA, "./Data/point.ozone.AEA.geojson", driver = "GeoJSON",append=FALSE)





#######################################################################
#######################################################################
#######################################################################
rm(list=ls())

ozone.sf <- st_read("./Data/point.ozone.AEA.geojson")
county.data.sf <- st_read("./Data/county.AllData.AEA.geojson")

state.border.sf  <- states(cb = TRUE) %>%  
  select(NAME, geometry) %>% 
  filter(NAME == "California") %>%
  st_transform(3310)

county.border.sf <- counties("CA", cb = TRUE) %>% 
  select(GEOID, geometry)   %>%
  st_transform(3310)

state.border.terra <- vect(state.border.sf)


#--------------------------------------------------
# Part 1, explain what this code is doing
#--------------------------------------------------
ozone.terra      <- vect(ozone.sf)
ozone.voronoi.sf <- voronoi(ozone.terra)
ozone.voronoi.sf <- st_as_sf(crop(ozone.voronoi.sf,state.border.terra))

map.voronoi <- tm_shape(ozone.voronoi.sf) +
               tm_polygons("OZONE_1000PPB",palette="YlGnBu")+
               tm_legend(position = c("left", "bottom"))+
               tm_layout(title="Voroni Tesselation of Ozone in CA (1000 PPB)",
                         title.size=.8,title.fontface=2)

map.county  <- tm_shape(county.data.sf) +
               tm_polygons("OZONE_1000PPB",palette="YlGnBu",breaks=seq(0,100,by=20))+
               tm_legend(position = c("left", "bottom"))+
               tm_layout(title="County Averages of Ozone in CA (1000 PPB)",
                         title.size=.8,title.fontface=2)

map.point   <- tm_shape(county.data.sf) +
               tm_borders()+
               tm_shape(ozone.sf) +
               tm_dots(col="OZONE_1000PPB",palette="YlGnBu",size=.4,alpha=.6)+
               tm_legend(position = c("left", "bottom"))+
               tm_layout(title="Point values of Ozone in CA (1000 PPB)",
                         title.size=.8,title.fontface=2)

tmap_arrange(map.point,map.voronoi,map.county)

# Explain the three plots
# What are the advantages and disadvantages of each way of summarising this data?

moran.data    <- county.data.sf[county.data.sf$OZONE_1000PPB >  0 , ]
tracts_empty  <- st_is_empty(moran.data)
moran.data    <- moran.data[which(tracts_empty==FALSE), ]

matrix.queen  <- poly2nb(moran.data, queen=TRUE)
weights.queen <- nb2listw(matrix.queen, style='B', zero.policy = T)
plot(weights.queen, st_coordinates(st_centroid(moran.data)), col='black', lwd=1, cex=.2)

moran.plot(moran.data$EP_POV150, 
           listw= weights.queen,
           xlab = "Value of Ozone per county (1000 PPB)",
           ylab = "Average Ozone in neighbouring counties (1000 PPB)",
           zero.policy = T)



# Choose your data
moran.data <- ozone.voronoi.sf
critical_threshold <- 0.1

s <- moran.data[moran.data$OZONE_1000PPB >  0 , ]
tracts_empty  <- st_is_empty(s)
s    <- s[which(tracts_empty==FALSE), ]



#https://mgimond.github.io/simple_moransI_example/

s <- moran.data[moran.data$OZONE_1000PPB >  0 , ]
tracts_empty  <- st_is_empty(s)
s    <- s[which(tracts_empty==FALSE), ]

names(s)
s$OZONE_1000PPB
hist(s$OZONE_1000PPB, main=NULL)
boxplot(s$OZONE_1000PPB, horizontal = TRUE)
tm_shape(s) + tm_fill(col="OZONE_1000PPB", style="quantile", n=8, palette="Greens") +
  tm_legend(outside=TRUE)


nb <- poly2nb(s, queen=TRUE)
nb[1]
lw <- nb2listw(nb, style="W", zero.policy=TRUE)
lw$weights[1]
inc.lag <- lag.listw(lw, s$Income)
inc.lag


# try adding in scale
plot(inc.lag ~ s$Income, pch=16, asp=1)
M1 <- lm((inc.lag) ~ (s$Income))
abline(M1, col="blue")
coef(M1)[2]
I <- moran(s$OZONE_1000PPB, lw, length(nb), Szero(lw))[1]
I
moran.test(s$Income,lw, alternative="greater")



MC<- moran.mc(s$OZONE_1000PPB, lw, nsim=999, alternative="greater")
MC
plot(MC)
set.seed(131)
s$rand1 <- sample(s$OZONE_1000PPB, length(s$OZONE_1000PPB), replace = FALSE)
s$rand2 <- sample(s$OZONE_1000PPB, length(s$OZONE_1000PPB), replace = FALSE)
s$rand3 <- sample(s$OZONE_1000PPB, length(s$OZONE_1000PPB), replace = FALSE)

tm_shape(s) + tm_fill(col=c("OZONE_1000PPB", "rand1", "rand2", "rand3"),
                      style="quantile", n=8, palette="Greens", legend.show = FALSE) + tm_facets( nrow=1)







#---------------------------------------------------
# Note, you can make the Moran's plot more easily these days
#---------------------------------------------------

moran.plot(s$OZONE_1000PPB, 
           listw= lw,
           xlab = "Value of Ozone per polygon (1000 PPB)",
           ylab = "Average Ozone in neighbouring polygons (1000 PPB)",
           zero.policy = T)

moran.test(s$OZONE_1000PPB, listw= lw)
plot(moran.mc(s$OZONE_1000PPB, lw, nsim=999, alternative="greater"))


#---------------------------------------------------
# We can either create the Raw LISA values using the theoretical approach or Monte Carlo
# 
#  - Padjust allows us to account for mulitple testing issues.
#  - adjust.x allows us to omit polygons with no neighbours as needed
#  - zero.policy allows us to keep but ignore polygons with no neighbours
#  - nsim: number of simulations, If this crashes the cloud try reducing to say 50
#---------------------------------------------------
LocalMoran_Output <- localmoran_perm(x = s$OZONE_1000PPB, 
                                     listw = lw, 
                                     nsim=499,
                                     adjust.x = TRUE,
                                     zero.policy = T)

#---------------------------------------------------
# The column names are less intuitive to newbies, so let's rename them 
# (to see the originals, type ?localmoran into the console)
# The skew and kurtosis are because the histogram of IRP I might not look normal
# To do this I force the output to be a standard dataframe
#---------------------------------------------------
LocalMoran_Table <- as.data.frame(LocalMoran_Output)
names(LocalMoran_Table) <- c("Observed_LocalI", "IRP_estimate_LocalI",   "Variance_I",
                             "ZScore_I",  "PValue_Theoretical",  "PValue_MonteCarlo",
                             "Skew_I","Kurtosis_I")

#---------------------------------------------------
# Merge with our data
#---------------------------------------------------
s <- cbind(s,LocalMoran_Table[,4:6])

#---------------------------------------------------
# It turns out that the quadrants are secretly stored in the output, 
# so we don't even have to manually calculate them
#---------------------------------------------------
s$LISA_Quadrant <- attr(LocalMoran_Output,"quadr")$mean

#---------------------------------------------------
# Remove polygons that are unlikely to be significant
# YOU GET TO CHOOSE THE THRESHOLD HERE
# Find the rows where your p-value is > your threshold
#---------------------------------------------------
RowsOverThreshold <- which(s$PValue_MonteCarlo > critical_threshold) 

#---------------------------------------------------
# Make the column a character to make life easy
# Then rename the quadrant in those rows p>0.05, or 
# whatever your threshold is
#---------------------------------------------------
s$LISA_Quadrant_Plot <- as.character(s$LISA_Quadrant)
s$LISA_Quadrant_Plot[RowsOverThreshold] <- paste("P>",critical_threshold,sep="")

#---------------------------------------------------
# Make a factor again
#---------------------------------------------------
s$LISA_Quadrant_Plot <- factor(s$LISA_Quadrant_Plot,
                               levels= c("High-High",  
                                                         "High-Low",
                                                         "Low-High",
                                                         "Low-Low", paste("P>",critical_threshold,sep="")))

#---------------------------------------------------
# And make a plot
#---------------------------------------------------
mapLISA <- 
  tm_shape(s)  +
  tm_fill( "LISA_Quadrant_Plot",id="NAME",  alpha=.6,
           palette=  c("#ca0020","#f4a582","#92c5de","#0571b0","white"), title="") +
  tm_borders(alpha=.5) +
  tm_legend(position = c("left", "bottom"))+
  tm_layout(title = paste("LISA,sig= ",critical_threshold))


#---------------------------------------------------
# plot against the raw values 
#---------------------------------------------------
mapRaw <- tm_shape(s) + 
  tm_fill(col = "OZONE_1000PPB", 
              style = "pretty",
              palette = "YlGnBu", 
              border.alpha = 0, 
              title = "",alpha=.6) +  
  tm_layout(title = "Mean Ozone",  main.title.size = 0.95)+
  tm_borders(alpha=.5)+
  tm_legend(position = c("left", "bottom"))
  


tmap_options(check.and.fix = TRUE)
tmap_mode("plot")
tmap_arrange(mapRaw,mapLISA)



list_nb <- poly2nb(s, queen = TRUE)
empty_nb <- which(card(list_nb) == 0)
empty_nb 
tes_subset <- s[-empty_nb, ]
empty_polygons <- s[empty_nb, ]
empty_polygons$nghbrhd  # print neighborhood names

tes_subset <- s

# Now that we removed empty neighbor sets (tes_subset)
# Identify neighbors with queen contiguity (edge/vertex touching)
tes_nb <- poly2nb(tes_subset, queen = TRUE)

# Binary weighting assigns a weight of 1 to all neighboring features 
# and a weight of 0 to all other features
tes_w_binary <- nb2listw(tes_nb, style="B")

# Calculate spatial lag of TreEqty
tes_lag <- lag.listw(tes_w_binary, tes_subset$OZONE_1000PPB)


globalG.test(tes_subset$OZONE_1000PPB, tes_w_binary)

# Identify neighbors, create weights, calculate spatial lag
tes_nbs <- s |> 
  mutate(
    nb = st_contiguity(geometry),        # neighbors share border/vertex
    wt = st_weights(nb),                 # row-standardized weights
    tes_lag = st_lag(OZONE_1000PPB, nb, wt)    # calculate spatial lag of TreEqty
  ) 

# Calculate the Gi using local_g_perm
tes_hot_spots <- tes_nbs |> 
  mutate(
    Gi = local_g_perm(OZONE_1000PPB, nb, wt, nsim = 999)
    # nsim = number of Monte Carlo simulations (999 is default)
  ) |> 
  # The new 'Gi' column itself contains a dataframe 
  # We can't work with that, so we need to 'unnest' it
  unnest(Gi) 


# Cursory visualization
# Plot looks at gi values for all locations
tes_hot_spots |> 
  ggplot((aes(fill = gi))) +
  geom_sf(color = "black", lwd = 0.15) +
  scale_fill_gradient2() # makes the value 0 (random) be the middle

