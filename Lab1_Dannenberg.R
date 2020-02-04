# Lab 1 practice session

# 1) Defining variables
lat <- 35.9736
-79.1004 -> lon

vars <- c('Month','Day','Hour','Minute','Tdr','Tsoil','Wind','Lai','Prec','Par','Ta','Vpd')
vars[1]

# 2) Load some data
dat <- read.table('./data/TdrTsoilWindLaiPrecParTaVpdMay2001-Daytime.txt', sep=' ', header=FALSE)
names(dat) <- vars

# 3) Working with data frames
dat$Year <- 2001
Ta <- dat$Ta
dat$Vpd <- dat$Vpd * 1000 # kPa -> Pa

# 4) Calculate
mean(dat$Tdr)
min(dat$Tsoil)
max(dat$Ta)
sum(dat$Prec)

# 5) Calculate dates
dat$Date <- ISOdate(dat$Year, dat$Month, dat$Day, dat$Hour, dat$Minute, tz='EST')

# 6) Plot a time series
plot(dat$Date, dat$Ta,
     type='l',
     xlab='Time',
     ylab='Air temperature (degrees C)',
     frame=FALSE)
lines(dat$Date, dat$Tsoil, col='red')







