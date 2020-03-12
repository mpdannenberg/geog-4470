# Lab 4

# Load model functions
source('eco_model.r')
library(pracma)

# Load data and add a header
dat <- read.table('./data/TdrTsoilWindLaiPrecParTaVpdMay2001-Daytime.txt', sep=' ', header=FALSE)
vars <- c('Month','Day','Hour','Minute','Tdr','Tsoil','Wind','Lai','Prec','Par','Ta','Vpd')
names(dat) <- vars

# Add year and date to the data frame
dat$Year <- 2001
dat$Date <- ISOdate(dat$Year, dat$Month, dat$Day, dat$Hour, dat$Minute, tz='EST')
dat$Vpd <- dat$Vpd * 1000

# Define global variables and model parameters
define_global_variables()

# Calculate solar zenith, elevation, azimuth, and declination angles
sun <- compute_sun_angles(LATI, LONGI, dat$Year, dat$Month, dat$Day, dat$Hour, dat$Minute)

# Calculate top of canopy diffuse/direct components of PAR and NIR
rad <- compute_toc_down_rad(dat$Par, sun$zenith, sun$jday)

# Leaf area index
lai <- dat$Lai

# Leaf orientation
# x=1 (the default) means that plants have a "spherical" leaf angle distribution (leaves are oriented in all directions)
# x<1 means that plants have more horizontally oriented leaves
# x>1 means that plants have more vertically oriented leaves
x <- 1

# calculate sunlit and shaded components of LAI based on the leaf angle distribution (x)
i <- seq(from=0, to=89, by=1) # Zenith angles from zero to 89 degrees
kb <- sqrt(x^2 + tan(i*pi/180)^2 / (x + 1.774 * (x+1.182)^-0.733))
beer <- function(lai, kb, omega) exp(-kb * omega * lai) # Beer's law
trans_D <- outer(lai, kb, beer, omega=omega) # apply function across all LAIs and angles
trans_d <- apply(trans_D, 1, function(trans_D, i) trapz(i, 2*trans_D*sin(i*pi/180.0)*cos(i*pi/180.0)*(pi/180)), i=i) # integrate across all angles

Kd <- -log(trans_d)/lai
Kb <- omega * sqrt(x^2 + tan(sun$zenith)^2) / (x + 1.774*(x+1.182)^-0.733)

sunlit_lai <- (1 - exp(-Kb*lai)) / Kb
shaded_lai <- lai - sunlit_lai

# Compute canopy shortwave radiation
canopy_shortwave <- compute_canopy_intercepted_shortwave_rad(rad$par_D,
                                                             rad$par_d, 
                                                             rad$nir_D, 
                                                             rad$nir_d, 
                                                             lai, 
                                                             Kb, 
                                                             Kd)

# Compute canopy longwave radiation
canopy_longwave <- compute_longwave_rad(dat$Ta, dat$Tsoil, dat$Vpd, rad$cloud)

# Compute net radiation
canopy_net_rad <- compute_canopy_net_rad(canopy_shortwave, 
                                         canopy_longwave,
                                         rad, 
                                         sunlit_lai, 
                                         shaded_lai, 
                                         Kb, 
                                         Kd)

##### NEW FOR LAB 4 #####

# Estimate aerodynamic resistance and conductance
ra <- canopy_aerodynamic_resistance(dat$Wind, hs, h2)
ra[dat$Wind <= 0] <- 1.0 / g1_blayer
ga <- 1/ra

# Compute soil resistance
rs_u <- compute_soil_surface_resistance(dat$Tdr)
ra_u <- ra

# Calculate evapotranspiration using Penman-Monteith model
gs_sunlit <- 0.025
LE_sunlit <- penman_monteith(dat$Ta, 
                             canopy_net_rad$Rnet_sunlit, 
                             gs_sunlit,
                             ga, 
                             dat$Vpd)

gs_shaded <- 0.01
LE_shaded <- penman_monteith(dat$Ta, 
                             canopy_net_rad$Rnet_shaded, 
                             gs_shaded,
                             ga, 
                             dat$Vpd)

LE_understory <- penman_monteith(dat$Ta, 
                                 canopy_net_rad$Rnet_floor, 
                                 1/rs_u,
                                 1/ra_u, 
                                 dat$Vpd)

# Measured heat fluxes
flx <- read.table('./data/GppLeHMay2001-Daytime.txt', sep=' ', header=T)


# Question 1: make separate time series of wind and aerodynamic RESISTANCE, and a scatterplot of wind vs. ra. Interpret plots: when is resistance highest? Diurnal patterns? What does this mean for fluxes?
plot(dat$Date, dat$Wind, 
     type='l', 
     xlab='Date',
     ylab='Wind speed (m/s)',
     frame=F)

plot(dat$Date, ra, 
     type='l', 
     xlab='Date',
     ylab='Aerodynamic resistance (s/m)',
     frame=F)

plot(dat$Wind, ra, 
     xlab='Wind speed (m/s)', 
     ylab='Aerodynamic resistance (s/m)',
     frame=F)

# Question 2: make a time series of soil surface resistance (rs_u) and a scatterplot of rs_u vs. soil moisture (Tdr). How does rs_u vary with soil moisture?
plot(dat$Date, dat$Tdr, 
     type='l', 
     xlab='Date',
     ylab='Soil moisture (%)',
     frame=F)

plot(dat$Date, rs_u, 
     type='l', 
     xlab='Date',
     ylab='Soil resistance (s/m)',
     frame=F)

plot(dat$Tdr, rs_u, 
     xlab='Soil moisture (%)', 
     ylab='Soil resistance (s/m)',
     frame=F)


# Question 3: calculate total LE of the ecosystem (LE_total). Plot LE_total and its three components as a time series

# Question 4: Make a scatterplot of estimated vs. measured LE. How do they compare? What do you think might cause some of the differences?

# Quesiton 5: Do something fun.
