# Lab 7

# Load data and add a header
dat <- read.table('./data/GppParTminVpdRedNir-daily.csv', 
                  sep=',', 
                  header=TRUE,
                  na.strings = c('NaN',NaN))
dat$Date <- as.Date(seq(from=0, to=dim(dat)[1]-1, by=1), origin = "2001-01-01")
dat$VPD <- dat$VPD * 100 # hPa -> Pa

# Question 1: calculate NDVI and FPAR
dat$NDVI <- (dat$NIR - dat$Red) / (dat$NIR + dat$Red)
# Replace this line with your FPAR calculation (dat$FPAR <- ???)
# Replace this line with your APAR calculation (dat$APAR <- ???)

# Question 2: environmental stress
T1 <- -6.0
T2 <- 10.0
D1 <- 650.0
D2 <- 1650.0

dat$fT <- (dat$Tmin - T1) / (T2 - T1)
dat$fT[dat$Tmin < T1] <- 0 # f(T)=0 when Tmin is less than T1
dat$fT[dat$Tmin >= T2] <- 1 # f(T)=1 when Tmin is greater than T2

dat$fD <- (D2 - dat$VPD) / (D2 - D1)
dat$fD[dat$VPD < D1] <- 1
dat$fD[dat$VPD >= D2] <- 0

dat$fE <- apply(dat[,c('fT','fD')], 1, min)

# Question 3: GPP
lue <- 1.165 # g C MJ-1
# replace this line with your GPP calculation (dat$GPP_est <- ???)

cor(dat$GPP, dat$GPP_est, use = "complete.obs")^2 # R^2
mean(abs(dat$GPP - dat$GPP_est), na.rm=TRUE) # mean absolute error
mean((dat$GPP - dat$GPP_est), na.rm=TRUE) # mean bias

