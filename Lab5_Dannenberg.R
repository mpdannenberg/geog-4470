# Lab 5

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

# Estimate aerodynamic resistance and conductance
ra <- canopy_aerodynamic_resistance(dat$Wind, hs, h2)
ra[dat$Wind <= 0] <- 1.0 / g1_blayer
ga <- 1/ra

# Compute soil resistance
rs_u <- compute_soil_surface_resistance(dat$Tdr)
ra_u <- ra

#########################
##### New for Lab 5 #####
#########################

# Estimate photosynthesis parameters
psn_para <- compute_psn_parameters(dat$Ta, Kb, sunlit_lai, shaded_lai)
fsoil_theta <- 1.0

# Calcuate photosynthesis and stomatal conductance for sunlit leaves
sunlit <- 1;
psn_sunlit <- compute_bwb_farq_psn(psn_para, 
                                   canopy_shortwave$sunlit_apar,
                                   dat$Ta, 
                                   dat$Vpd, 
                                   fsoil_theta, 
                                   sunlit)

out <- canopy_bwb_stomatal_conductance(psn_sunlit, dat$Ta, m, dat$Vpd)
gs_sunlit <- out$gs # Conductance for water (m/s)
gc_sunlit <- out$gc # Conductance for CO2 (m/s)

rs_sunlit <- 1/gs_sunlit # Resistance for water (s/m)
rs_sunlit[gs_sunlit<=0] <- 1e6 

# Calculate photosynthesis and stomatal conductance for shaded leaves
sunlit <- 0
psn_shaded <- compute_bwb_farq_psn(psn_para, 
                                   canopy_shortwave$shaded_apar,
                                   dat$Ta, 
                                   dat$Vpd, 
                                   fsoil_theta, 
                                   sunlit)

out <- canopy_bwb_stomatal_conductance(psn_shaded, dat$Ta, m, dat$Vpd)
gs_shaded <- out$gs # Conductance for water
gc_shaded <- out$gc # Conductance for CO2

rs_shaded <- 1/gs_shaded
rs_shaded[gs_shaded<=0] <- 1e6

########################################
##### Now return to Lab 4 material #####
########################################

# Calculate evapotranspiration using Penman-Monteith model
LE_sunlit <- penman_monteith(dat$Ta, 
                             canopy_net_rad$Rnet_sunlit, 
                             gs_sunlit,
                             ga, 
                             dat$Vpd)

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

# Measured carbon, water, and heat fluxes
flx <- read.table('./data/GppLeHMay2001-Daytime.txt', sep=' ', header=T)

# Question 1: Make time series plots of gs_sunlit and gs_shaded on the same graph. How do they vary diurnally and seasonally? Why do you think?

# Question 2: Make scatterplots of gs_sunlit vs. PAR, temperature, VPD, and soil moisture (Tdr). Compare these to the same plots in the Bonan book. How are they similar/different? Why do you think they might be different and how could you make them more comparable?

# Question 3: Make a scatterplot of estimated vs. measured total LE. How do they compare? What do you think might cause some of the differences? How does it compare to LE in Lab 4 with static gs?

# Quesiton 4: Do something fun.
