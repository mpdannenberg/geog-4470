## Functions for ecological model

date2doy <- function(yr, mo, dy){
  yr <- as.character(yr)
  mo <- as.character(mo)
  dy <- as.character(dy)
  
  dt <- paste(c(yr,'-',mo,'-',dy), collapse='')
  doy <- strftime(dt, format='%j')
  return(doy)
}

define_global_variables <- function(){
  
  # Site constants
  PAR2SWIR <<- 0.43 # PAR to shortwave radiation ratio
  STD_MERAD <<- -75.0 # Time zone longitude
  LATI <<- 35.9736 # Site latitude (CHANGE TO MMS)
  LONGI <<- -79.1004 # Site longitude (CHANGE TO MMS)
  hs <<- 20.0 # Screen height
  h2 <<- 16.2 # Default canopy height
  PA <<- 101325.0 # Air pressure in Pa
  CP <<- 1010.0 # Specific heat of air (J kg-1 K-1)
  CA <<- 410.0 # Atmospheric CO2 concentration (ppm)
  
  # Physical constants
  SO <<- 1367.0 # Solar constant (W m-2)
  
  # Leaf constants
  alpha_par <<- 0.80
  alpha_nir <<- 0.20
  rho_s <<- 0.10
  x <<- 1.0 # Leaf angle distribution parameter
  omega <<- 1.0 # Clumping index
  m <<- 8.0 # Ball-Berry model slope
  g1_blayer <<- 0.005 # Leaf boundary layer conductance in m/s
  g0 <<- 0.010 # Parameter in Leuning's gs function in mol m-2 s-1
  
  # Soil constants
  Ksat <<- 0.6 # Saturated conductance (m day-1)
  wilting_p <<- 0.08 # Wilting poit soil water content (unitless)
  porosity <<- 0.54 # Porosity at depth zero (unitless)
  p1 <<- 0.478 # Air entry pressure in meters of water
  p2 <<- 0.186 # Pore size index (unitless)
  p_0 <<- 4000.0 # Porosity decay in m-1
  m_z <<- 0.24 # Saturated conductivity decay parameter (in m-1)
  soil_depth <<- 0.325 # Soil depth in rooting zone
  
}

canopy_aerodynamic_resistance <- function(u_h, h, h_o){
  h_u <- 2.0
  cn <- 0.002 # wind attenuation coefficient
  
  # Compute zero plane displacement
  d_o <- 0.7 * h_o
  
  # Compute the roughness length
  zo_o <- 0.1 * h_o
  u_o <- compute_toc_wind(u_h, h, h_o)
  u_m <- (u_o * exp(cn * (0.5 * h_o)/h_o - 1))
  u_m <- sapply(u_m, function(x) max(x,0))
  
  if(h<d_o) h <- d_o + 0.0001
  
  ra <- (log((h_o-d_o)/zo_o)/0.41)^2 / u_m
  ra[u_m==0] <- 9999999
  
  return(ra)
}

canopy_bwb_stomatal_conductance <- function(psn, Tair, a1, vpda){
  svp <- 610.7 * exp(17.38 * Tair / (239.0+Tair))
  rh <- 1.0 - vpda/svp
  gc <- g0 + a1 * psn * rh/CA
  
  # convert mol m-2 s-1 to m s-1
  gc <- gc * (22.4 / 1000.0)
  
  # convert gc (Co2) to gs (H2O)
  gs <- 1.56 * gc
  
  return(data.frame(gc=gc, gs=gs))
}

compute_bwb_farq_psn <- function(psnp, apar_wm2, Tair, vpda, fstheta, leaftype){
  svp <- 610.7 * exp(17.38 * Tair / (239.0 + Tair))
  rh <- 1.0 - vpda/svp
  apar <- apar_wm2*4.55
  Ca <- CA
  
  # calculate atmospheric O2 in Pa, assuming 21% O2 by volume
  O2 <- psnp$O2
  Ko <- psnp$Ko
  
  # Convert Pa to ppm
  Kc <- psnp$Kc * 1e6/PA
  gamma <- psnp$gamma * 1e6/PA
  
  ppe <- psnp$ppe
  if(leaftype==1){
    Vmax <- psnp$Vmax_sunlit
    Rd <- psnp$Rd_sunlit
  } else {
    Vmax <- psnp$Vmax_shaded
    Rd <- psnp$Rd_shaded
  }
  
  Jmax <- 2.1 * Vmax
  
  g1 <- 0.01 + m * rh/Ca
  K <- Kc * (1.0 + O2/Ko)
  
  # Calculate rubisco-limited assimilation rate
  aa <- Ca*g1 + g1*K - 1.0
  bb <- Ca*g0 + g0*K - Vmax*Ca*g1 + Vmax + Vmax*g1*gamma + Rd*Ca*g1 - Rd + Rd*g1*K
  cc = -Vmax*Ca*g0 + Vmax*gamma*g0 + Rd*Ca*g0 + Rd*K*g0
  
  det <- bb*bb - 4.0*aa*cc
  Av <- (-bb + sqrt(det)) / (2.0*aa)
  Av[det<0] <- -1.0
  
  # Calculate RuBP regeneration limited assimilation rate
  aa <- 0.7
  bb <- -Jmax - (apar/ppe)
  cc <- Jmax * apar/ppe
  J <- (-bb - sqrt(bb*bb - 4*aa*cc)) / (2.0/aa)
  
  aa <- 4.5*Ca*g1 - 4.5 + 10.5*gamma*g1
  bb <- 4.5*Ca*g0 + 10.5*gamma*g0 - J*g1*Ca + J + gamma*g1*J + 4.5*Rd*Ca*g1-4.5*Rd+10.5*gamma*g1*Rd
  cc <- -J*Ca*g0+4.5*Rd*Ca*g0+10.5*gamma*g0*Rd
  
  det <- bb*bb - 4.0*aa*cc
  Aj <- (-bb+sqrt(det))/(2.0*aa)
  Aj[det<0] <- -1.0
  
  # Calculate net assimilation rate
  A = pmin(Aj,Av)
  A[Av == -1.0 & Aj == -1.0] <- 0.0
  A[Av == -1.0 & Aj != -1.0] <- Aj[Av == -1.0 & Aj != -1.0]
  A[Aj == -1.0 & Av != -1.0] <- Av[Aj == -1.0 & Av != -1.0]
  
  return(A)
}

compute_canopy_intercepted_shortwave_rad <- function(par_D, par_d, nir_D, nir_d, lai, Kb, Kd){
  Q_sc_par <- par_D*(exp(-Kb*sqrt(alpha_par)*lai)-exp(-Kb*lai))/2.0
  Q_d_par <- par_d*(1.0-exp(-Kd*sqrt(alpha_par)*lai))/(Kd*sqrt(alpha_par)*lai)
  Q_par_soil <- par_D*exp(-Kb*sqrt(alpha_par)*lai)+par_d*exp(-Kd*sqrt(alpha_par)*lai)
  Q_nir_soil <- nir_D*exp(-Kb*sqrt(alpha_nir)*lai)+nir_d*exp(-Kd*sqrt(alpha_nir)*lai)
  
  sunlit_apar <- alpha_par*(par_D*Kb+ Q_sc_par + Q_d_par)
  shaded_apar <- alpha_par*(Q_sc_par + Q_d_par)
  
  Q_sc_nir <- nir_D*(exp(-Kb*sqrt(alpha_nir)*lai)-exp(-Kb*lai))/2.0
  Q_d_nir <- par_d*(1.0-exp(-Kd*sqrt(alpha_par)*lai))/(Kd*sqrt(alpha_par)*lai)
  
  sunlit_anir <- alpha_nir*(nir_D*Kb+Q_sc_nir+Q_d_nir)
  shaded_anir <- alpha_nir*(Q_sc_nir+Q_d_nir)
  
  rad <- data.frame(sunlit_apar=sunlit_apar,
              shaded_apar=shaded_apar,
              sunlit_anir=sunlit_anir,
              shaded_anir=shaded_anir)
  return(rad)
}

compute_canopy_net_rad <- function(rad_short, rad_long, toc_rad, sunlit_lai, shaded_lai, Kb, Kd){
  
  lai <- sunlit_lai + shaded_lai
  
  Rnet_sunlit <- rad_short$sunlit_apar + rad_short$sunlit_anir +
    rad_long$floor * (1.0 - exp(-Kd*lai)) / (Kd*lai) +
    rad_long$air * (1.0 - exp(-Kd*lai)) / (Kd*lai) +
    -2.0 * rad_long$canopy * (1.0 - exp(-Kd*lai)) / (Kd*lai)
  
  Rnet_shaded <- rad_short$shaded_apar + rad_short$shaded_anir + 
    rad_long$floor * (1.0 - exp(-Kd*lai)) / (Kd*lai) + 
    rad_long$air * (1.0 - exp(-Kd*lai)) / (Kd*lai) + 
    -2.0 * rad_long$canopy * (1.0 - exp(-Kd*lai)) / (Kd*lai)
  
  Rnet_floor <- (1.0 - rho_s) * (toc_rad$par_D * exp(-Kb*sqrt(alpha_par)*lai)+ 
                                  toc_rad$par_d * exp(-Kd*sqrt(alpha_par)*lai)+ 
                                  toc_rad$nir_D * exp(-Kb*sqrt(alpha_nir)*lai)+
                                  toc_rad$nir_d * exp(-Kd*sqrt(alpha_nir)*lai))+
                                  rad_long$air * exp(-Kd*lai)+ 
                                  rad_long$canopy * (1.0-exp(-Kd*lai))-
                                  rad_long$floor
  
  Rnet_canopy <- sunlit_lai*Rnet_sunlit + shaded_lai*Rnet_shaded
  Rnet_stand <- Rnet_canopy + Rnet_floor
  
  out <- data.frame(Rnet_sunlit = Rnet_sunlit, 
                    Rnet_shaded = Rnet_shaded, 
                    Rnet_floor = Rnet_floor, 
                    Rnet_canopy = Rnet_canopy, 
                    Rnet_stand = Rnet_stand)
  return(out)
  
}

compute_longwave_rad <- function(Tair, Tsoil, vpd, cloud){
  # Compute emitted longwave radiation from the canopy, floor, and air
  
  N = cloud
  e_s = 610.7 * exp(17.38 * Tair / ( 239.0 + Tair))
  e_a = (e_s-vpd)/100.0
  
  epsilon_air = 1.24 * (e_a/(Tair+273.13))^(1.0/7.0)
  epsilon_air =(1.0-0.84*N)*epsilon_air + 0.84*N
  
  Ln_air = epsilon_air*(5.6704e-8) * (Tair+273.13)^4.0
  Ln_canopy = 1.0*(5.6704e-8) * (Tair+273.13)^4.0
  Ln_floor = 0.97*(5.6704e-8) * (Tsoil+273.13)^4.0
  
  emitted_rad <- data.frame(air=Ln_air, canopy=Ln_canopy, floor=Ln_floor, air_emissivity=epsilon_air)
  return(emitted_rad)
  
}

compute_potential_exfiltration <- function(Sr, depth_s, K_sat, Ksat_decay, psi_air_entry, pore_size_index, pore_decay, pore0){
  # Estimate mean porosity
  porosity_average <- pore0 * pore_decay * (1.0-exp(-1.0*depth_s/pore_decay))
  
  # Estimate mean saturated conductivity
  if (Ksat_decay > 0){
    Ksat_average <- Ksat_decay* Ksat *(1.0-exp(-1.0*depth_s/Ksat_decay))
  } else Ksat_average <- Ksat
  
  S <- max(0,min(Sr,1))
  
  # Plug everything into the equation for max infiltration
  potential_exfiltration <- S ^ ((1.0 / (2.0*pore_size_index))+2.0) * 
                            sqrt((8.0 * porosity_average * 
                            Ksat_average * psi_air_entry) / 
                            (3.0 * (1.0 + 3.0 * pore_size_index) * 
                            (1.0 + 4.0 * pore_size_index)))
  
  potential_exfiltration = min(0.001, potential_exfiltration)
  pe_wm2 = potential_exfiltration * 1e6 * 597.3 * 4.18 / (24.0*3600.0)
  return(pe_wm2)
}

compute_psn_parameters <- function(t, Kb, L_sunlit, L_shaded){
  # Local static variables
  Kc25 <- 404.0  # (ubar) MM const carboxylase, 25 deg C
  q10Kc <- 2.1 # (DIM) Q_10 for kc
  Ko25 <- 248.0 # (mbar) MM const oxygenase, 25 deg C
  q10Ko <- 1.2 # (DIM) Q_10 for ko
  act25 <- 3.6 # (umol/mgRubisco/min) Rubisco activity
  q10act <- 2.4 # (DIM) Q_10 for Rubisco activity
  Kn <- 0.52 # extinction coefficient for Nitrogen distribution
  cica <- 0.67 # Ci to Ca ratio 
  alpha <- 0.10 # quantum yield efficiency
  
  L <- L_sunlit + L_shaded
  
  # Calculate atmospheric O2 in Pa, assumes 21% O2 by volume
  O2 <- 0.21 * PA
  
  # Correct kinetic constants for temperature, and do unit conversions
  Ko <- Ko25 * q10Ko^((t-25.0)/10.0)
  Ko <- Ko * 100.0 # mbar --> Pa
  
  Kc <- Kc25 * (1.8*q10Kc)^((t-25.0)/10.0) / q10Kc
  Kc[t>15.0] <- Kc25 * q10Kc^((t[t>15.0]-25.0)/10.0)
  act <- act25 * (1.8*q10act)^((t-25.0)/10.0) / q10act
  act[t>15.0] <- act25 * q10act^((t[t>15.0]-25.0)/10.0)
  
  Kc <- Kc * 0.10 # ubar --> Pa
  act <- act * 1e6/60.0 # umol/mg/min --> umol/kg/s
  
  # Calculate gamma (Pa), assumes Vomax/Vcmax = 0.21
  gamma <- 0.5 * 0.21 * Kc * O2 / Ko
  
  # Calculate Vmax from leaf nitrogen data and Rubisco activity
  Vmax25 <- 59.0 # umol/m2leaf/s at the top of the canopy (Lai et al, 2002, PCE, 25:1095-1119)
  
  Vmax25_canopy <- L * Vmax25 * (1.0 - exp(-L*Kn)) / (Kn*L)
  Vmax25_sunlit <- L * Vmax25 * (1.0 - exp(-Kn-Kb*L)) / (Kn+Kb*L)
  Vmax25_shaded <- (Vmax25_canopy - Vmax25_sunlit) / L_shaded
  Vmax25_sunlit <- Vmax25_sunlit / L_sunlit
  
  Vmax_sunlit <- Vmax25_sunlit * exp(0.051*(t-25.0))/(1.0+exp(0.205*(t-41.0)))
  Vmax_shaded <- Vmax25_shaded * exp(0.051*(t-25.0))/(1.0+exp(0.205*(t-41.0)))
  
  Rd_sunlit <- 0.015 * Vmax_sunlit
  Rd_shaded <- 0.015 * Vmax_shaded
  
  ppe <- 1.0 / (alpha*(4.0*cica*CA+2.0*gamma)/(cica*CA-gamma))
  
  psn_para <- data.frame(gamma=gamma,
                   Vmax_sunlit=Vmax_sunlit,
                   Vmax_shaded=Vmax_shaded,
                   Ko=Ko,
                   Kc=Kc,
                   O2=O2,
                   Rd_sunlit=Rd_sunlit,
                   Rd_shaded=Rd_shaded,
                   ppe=ppe)
  
  return(psn_para)
}

compute_soil_surface_resistance <- function(theta){
  gs <- 1.0 / (-83000.0 * theta + 16100.0)
  gs[theta > 0.185] <- 0.001429
  
  rs <- 1.0 / gs
  return(rs)
}

compute_sun_angles <- function(lat, lon, yr, mo, dt, hr, mi){
  solar_time <- hr + mi/60.0 + (lon-STD_MERAD)/15.0
  hangle <- (12.0-solar_time)*15.0*pi/180.0
  
  # Compute Julian day
  jday <- as.numeric(mapply(date2doy, yr=yr, mo=mo, dy=dt))
  
  # Sun declination angle
  declangle <- 23.45 * sin(2.0*pi*(284.0+jday)/365.0) * pi/180.0
  
  # sin(sun_elevation)
  sinh <- sin(lat*pi/180.0)*sin(declangle)+cos(lat*pi/180.0)*cos(declangle)*cos(hangle)
  
  zenith <- acos(sinh)
  azimuth <- asin(cos(declangle)*sin(hangle)/cos(pi/2.0-zenith))
  
  out <- data.frame(jday = jday, zenith = zenith, sun_decl = declangle, hourangle = hangle)
  return(out)
}

compute_toc_down_rad <- function(par, sun_zenith, jday){
  
  # Constants from Weiss & Norman 1985, Eq. 12
  A <- 0.9
  B <- 0.7
  C <- 0.88
  D <- 0.68
  
  # Cosine of the solar zenith angle
  cos_theta <- cos(sun_zenith)
  cos_theta[cos_theta<0] <- 0 
  
  # Potential visible (V) and NIR (N) radiation at top of atmosphere
  S_V <- 0.473*1367.0*(1.0+0.033*cos(2.0*pi*(jday-10)/365))
  S_N <- 0.527*1367.0*(1.0+0.033*cos(2.0*pi*(jday-10)/365))
  
  # Potential direct visible on a horizontal surface
  R_DV <- S_V*exp(-0.185/cos_theta)*cos_theta
  
  # Potential diffuse visible on a horizontal surface
  R_dV <- 0.4*(S_V-R_DV)*cos_theta
  
  # Total potential visible radiation on a horizontal surface
  R_V <- R_DV + R_dV
  
  # Absorbed NIR by atmospheric water vapor
  cos_theta_noZero <- cos_theta
  cos_theta_noZero[cos_theta==0] <- 0.0000000000000001 # avoid division by zero
  R_aN <- S_N * 10^(-1.195+0.4459*log10(1.0/cos_theta_noZero)-0.0345*(log10(1.0/cos_theta_noZero)^2.0))
  
  # Potential direct NIR on a horizontal surface
  R_DN <- (S_N*exp(-0.06/cos_theta)-R_aN)*cos_theta
  R_DN[R_DN<0] <- 0 # Can't have negative direct radiation
  
  # Potetial diffuse NIR on a horizontal surface
  R_dN <- 0.6*(S_N-R_DN-R_aN)*cos_theta
  
  # Total potential NIR on a horizontal surface
  R_N <- R_DN + R_dN
  
  # Ratio of actual total radiation to potential total radiation on the ground
  RATIO <- (par/4.55)*(1.0/PAR2SWIR) / (R_V+R_N)
  
  # Actual fraction of direct visible radiation
  RATIO[RATIO>A] <- A
  f_DV <- (R_DV/(R_V*0.928)) * (1.0 - ((A-RATIO)/B) ^ (2.0/3.0))
  f_DV[R_V==0] <- 0
  
  # Actual fraction of direct NIR radiation
  RATIO[RATIO>C] <- C
  f_DN <- (R_DN/R_N) * (1.0 - ((C-RATIO)/D) ^ (2.0/3.0))
  f_DN[R_N==0] <- 0
  
  rad <- data.frame(par_D = f_DV * par/4.55,
                    par_d = par/4.55 - f_DV * par/4.55,
                    nir_D = f_DN * (par/4.55) * ((1-PAR2SWIR)/PAR2SWIR),
                    nir_d = (par/4.55) * ((1-PAR2SWIR)/PAR2SWIR) - f_DN * (par/4.55) * ((1-PAR2SWIR)/PAR2SWIR))
                    
  return(rad)
  
}

compute_toc_wind <- function(u_h, h, z){
  # Returns wind speed at height z in stratum
  #   from RHESSys, which was from BGC 4.11.
  
  # Compute winds
  d_o <- 0.63 * z
  z_o <- 0.1 * z
  
  # Wind speed at toc log decrease from screen height to toc
  u_z <- u_h * log((z-d_o)/z_o) / log((h-d_o)/z_o)
  return(u_z)
}

penman_monteith <- function(Tair, Rnet, gs, ga, vpda){
  # Calculate forest canopy ET using penman-monteith equation
  
    # dt:     delta t, a small temperature increase
    # rho:    air density (kg/m2)
    # lhvap:  latent heat of water (J/kg)
    # s:      Slope of saturated vapor pressure - temp
    # es:     Saturated vapor pressure at t2 in Pa
    # ET:     Evapotranspiration in W/m2
    # gamma:  Psychrometric constant
    # z0:     Surface roughness
    # d0:     Zero plane displacement
  
  # Density of air (rho) as a f'n of air temp
  rho <- 1.292 - (0.00428*Tair)

  # Latent heat of vaporization as a f'n of air temp
  lhvap = 2.5023e6 - 2430.54*Tair
  
  # Saturation vapor pressures at Tair (Pa)
  es <- 610.7 * exp(17.38 * Tair / ( 239.0 + Tair))
  
  # Slope of es-T curve at Tair (Pa/deg C) (from Campbell & Norman 1998)
  s <- 17.38*239.0*es / (239.0+Tair)^2
  
  # Calculate gamma
  gamma <- CP * PA / ( 0.622*lhvap )
  
  # Evaporation in W/m2
  ET <- ((s*Rnet) + (rho*CP*vpda*ga)) / (gamma*(1.0 + ga/gs) +s)
}

