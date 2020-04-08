# Lab 6

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

# Calculate photosynthesis and stomatal conductance for shaded leaves
sunlit <- 0
psn_shaded <- compute_bwb_farq_psn(psn_para, 
                                   canopy_shortwave$shaded_apar,
                                   dat$Ta, 
                                   dat$Vpd, 
                                   fsoil_theta, 
                                   sunlit)

Anet <- (psn_sunlit-psn_para$Rd_sunlit)*sunlit_lai + (psn_shaded-psn_para$Rd_shaded)*shaded_lai

# Measured carbon, water, and heat fluxes
flx <- read.table('./data/GppLeHMay2001-Daytime.txt', sep=' ', header=T)

##### Question 1: Make time series plots of psn_sunlit and psn_shaded on the same graph. How do they vary diurnally and seasonally? Why do you think?

##### Question 2: Make scatterplots of Anet vs. PAR, psn_sunlit vs. canopy_shortwave$sunlit_apar, and psn_shaded vs. canopy_shortwave$shaded_apar. Describe the relationships. Based on these plots, if you wanted to make a simplified photosynthesis model, what could you do? What would you need to know and what assumptions would you need to make? 

##### Question 3: Make a scatterplot of estimated vs. measured Anet/GPP. How do they compare? What do you think might cause some of the differences?

##### Question 4: CO2 fertilization
CA <- c(280, 400, 560, 800, 1120)
psn_para <- compute_psn_parameters(25, 0.5, max(sunlit_lai), min(shaded_lai))

# CO2 fertilization of sunlit leaves
sunlit <- 1
psn_sunlit_co2 <- compute_bwb_farq_psn(psn_para, 
                                   max(canopy_shortwave$sunlit_apar),
                                   25, 
                                   1000, 
                                   fsoil_theta, 
                                   sunlit)

# CO2 fertilization of shaded leaves
sunlit <- 0
# Modify from lines 117-122 for shaded leaves!

##### Quesiton 5: Do something fun.

