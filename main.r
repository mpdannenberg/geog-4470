#########################################################
#####  Main file: estimate radiation, conductance,  #####
#####  photosynthesis, and evapotranspiration with  #####
#####             process based model               #####
#########################################################


# Load model functions
source('eco_model.r')
library(pracma)

# Load data
dat <- read.table('./data/TdrTsoilWindLaiPrecParTaVpdMay2001-Daytime.txt')
names(dat) <- c('Month', 'Day','Hour','Minute','Tdr','Tsoil','Wind','Lai','Prec','Par','Ta','Vpd')
dat$Year <- 2001

# Convert VPD from kPa to Pa
dat$Vpd <- dat$Vpd*1000

# Define global variables and model parameters
define_global_variables()

# Calculate solar zenith, elevation, azimuth, and declination angles
sun <- compute_sun_angles(LATI, LONGI, dat$Year, dat$Month, dat$Day, dat$Hour, dat$Minute)

# Calculate top of canopy diffuse/direct components of PAR and NIR
rad <- compute_toc_down_rad(dat$Par, sun$zenith, sun$jday)

# calculate sunlit and shaded components of LAI
i <- seq(from=0, to=89, by=1)
kb <- sqrt(x^2 + tan(i*pi/180)^2 / (x + 1.774 * (x+1.182)^-0.733))
beer <- function(lai, kb, omega) exp(-kb * omega * lai) # Beer's law
trans_D <- outer(dat$Lai, kb, beer, omega=omega) # apply function across all LAIs and angles
trans_d <- apply(trans_D, 1, function(trans_D, i) trapz(i, 2*trans_D*sin(i*pi/180.0)*cos(i*pi/180.0)*(pi/180)), i=i) # integrate across all angles

Kd <- -log(trans_d)/dat$Lai
Kb <- omega * sqrt(x^2 + tan(sun$zenith)^2) / (x + 1.774*(x+1.182)^-0.733)

sunlit_lai <- (1 - exp(-Kb*dat$Lai)) / Kb
shaded_lai <- dat$Lai - sunlit_lai

# Compute canopy shortwave radiation
canopy_shortwave <- compute_canopy_intercepted_shortwave_rad(rad$par_D,
                                                             rad$par_d, 
                                                             rad$nir_D, 
                                                             rad$nir_d, 
                                                             dat$Lai, 
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

out <- canopy_bwb_stomatal_conductance(psn_sunlit, dat$Ta, m, dat$Vpd)
gs_sunlit <- out$gs
gc_sunlit <- out$gc

rs_sunlit <- 1/gs_sunlit
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
gs_shaded <- out$gs
gc_shaded <- out$gc

rs_shaded <- 1/gs_shaded
rs_shaded[gs_shaded<=0] <- 1e6

# Calculate evapotranspiration using Penman-Monteith model
ET_sunlit <- penman_monteith(dat$Ta, 
                             canopy_net_rad$Rnet_sunlit, 
                             gs_sunlit,
                             ga, 
                             dat$Vpd)

ET_shaded <- penman_monteith(dat$Ta, 
                             canopy_net_rad$Rnet_shaded, 
                             gs_shaded,
                             ga, 
                             dat$Vpd)

ET_understory <- penman_monteith(dat$Ta, 
                                 canopy_net_rad$Rnet_floor, 
                                 1/rs_u,
                                 1/ra_u, 
                                 dat$Vpd)

S <- dat$Tdr / porosity

Exfil_soil <- compute_potential_exfiltration(S, 
                                             soil_depth, 
                                             Ksat, 
                                             m_z, 
                                             p1,
                                             p2, 
                                             p_0, 
                                             porosity)

ET_understory <- pmin(Exfil_soil, ET_understory)

# Add outputs to data frame
dat$apar <- canopy_shortwave$sunlit_apar*sunlit_lai + 
  canopy_shortwave$shaded_apar*shaded_lai
dat$rnet <- canopy_net_rad$Rnet_stand
dat$psn <- (psn_sunlit*sunlit_lai + psn_shaded*shaded_lai)*0.012
dat$A <- psn_sunlit*sunlit_lai + psn_shaded*shaded_lai # um/m^2/s
dat$gc <- 1000000*(gc_sunlit*sunlit_lai + gc_shaded*shaded_lai) / (22.4/1000.0) # um/m^2/s
dat$et <- ET_sunlit*sunlit_lai + ET_shaded*shaded_lai + ET_understory


