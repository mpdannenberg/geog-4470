# Lab 1 practice session

# When you use R for the first time in this class:
  # Create a "Labs" folder on your computer
  # Create "data" and "output" folders within the "Labs" folder
  # In RStudio, create a new project (save as Labs.Rproj)
  # Create a new script for each lab ("Lab1_<your last name here>.R")

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
legend(dat$Date[1], 10, legend=c('Air','Soil'), col=c('black','red'), lty=1)

# 7) Plot a scatterplot
plot(dat$Par, dat$Ta,
     xlab='PAR',
     ylab='Temperature',
     frame=FALSE)

# 8) Linear regression
mdl <- lm(dat$Ta~dat$Par)
mdl
summary(mdl)
abline(mdl)

# 9) Write a function to convert celcius to kelvin
c2k <- function(x){
  y <- x + 273.15
  return(y)
}
dat$Ta_kelvin <- c2k(dat$Ta)

# 10) Indexing
noon <- dat[dat$Hour==12 & dat$Minute==0, ]












