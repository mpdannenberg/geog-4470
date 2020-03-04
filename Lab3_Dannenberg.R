# Lab 3

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

# Question 1: plot sunlit and shaded leaf area
plot(dat$Date, sunlit_lai, type='l', ylim=c(0,5))
lines(dat$Date, shaded_lai, col=2)
# Question 2: plot direct and diffuse PAR and NIR for the canopy
# Question 3: plot net radiation of the canopy, soil, and stand
# Question 4: do something fun
