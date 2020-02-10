# Lab 2

# Load in model functions
source('eco_model.r')

# a vector of strings
vars <- c('Month','Day','Hour','Minute','Tdr','Tsoil','Wind','Lai','Prec','Par','Ta','Vpd')

# Load data and add a header
dat <- read.table('./data/TdrTsoilWindLaiPrecParTaVpdMay2001-Daytime.txt', sep=' ', header=FALSE)
vars <- c('Month','Day','Hour','Minute','Tdr','Tsoil','Wind','Lai','Prec','Par','Ta','Vpd')
names(dat) <- vars

# Add year and date to the data frame
dat$Year <- 2001
dat$Date <- ISOdate(dat$Year, dat$Month, dat$Day, dat$Hour, dat$Minute, tz='EST')

# Load all site variables
define_global_variables()

# Calculate sun angles
sun <- compute_sun_angles(LATI, LONGI, dat$Year, dat$Month, dat$Day, dat$Hour, dat$Minute)

# Calculate top of canopy radiation
rad <- compute_toc_down_rad(dat$Par, sun$zenith, sun$jday)

# Question 1: plot radiation
plot(dat$Date, rad$par_D, type='l', xlab='Time', ylab='Radiation (W m-2)', frame=FALSE)
lines(dat$Date, rad$par_d, col=2)
lines(dat$Date, rad$nir_D, col=3)
lines(dat$Date, rad$nir_d, col=4)
legend(dat$Date[1], 320, 
       legend=c('PAR, direct','PAR, diffuse', 'NIR, direct', 'NIR, diffuse'), 
       col=c(1,2,3,4), 
       lty=1)

# Question 2: Modify compute_toc_down_rad to return fractions of PAR and NIR that are diffuse and graph

# Question 3: Run compute_sun_angles and compute_toc_down_rad for +10 and -10 degrees latitude










