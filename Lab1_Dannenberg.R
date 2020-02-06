# Lab 1 practice session

# When you use R for the first time in this class:
  # Create a "Labs" folder on your computer
  # Create "data" and "output" folders within the "Labs" folder
  # In RStudio, create a new project (save as Labs.Rproj)
  # Create a new script for each lab ("Lab1_<your last name here>.R")

##### 1) Defining variables
# single numeric values
lat <- 35.9736
-79.1004 -> lon

# a vector of strings
vars <- c('Month','Day','Hour','Minute','Tdr','Tsoil','Wind','Lai','Prec','Par','Ta','Vpd')
vars[1]

##### 2) Load some data
# read.table function loads data into a data frame
# sep: value separating the columns, usually either a comma (','), a space (' '), or a tab ('\t')
# header: does the first row of the file include the names of the variables? (If yes, header=TRUE; if no, header=FALSE)
dat <- read.table('./data/TdrTsoilWindLaiPrecParTaVpdMay2001-Daytime.txt', sep=' ', header=FALSE)

# Name the variables of the data frame
names(dat) <- vars

##### 3) Working with data frames
# select (or create) individual columns of the data frame using $
dat$Year <- 2001
Ta <- dat$Ta
dat$Vpd <- dat$Vpd * 1000 # kPa -> Pa

##### 4) Calculate simple univariate statistics
# HINT: you need this for Lab 1, question 3
mean(dat$Tdr) # mean/average value
min(dat$Tsoil) # minimum value
max(dat$Ta) # maximum value
sum(dat$Prec) # sum of entire vector/column

##### 5) Calculate dates
# ISOdate function: convert individual date elements (yr, month, etc) into R-formatted dates
# HINT: you need this for Lab 1, question 1
dat$Date <- ISOdate(dat$Year, dat$Month, dat$Day, dat$Hour, dat$Minute, tz='EST')

# 6) Plot a time series
# create a plot of time on the x-axis (dat$Date) and a variable on the y-axis (e.g., dat$Ta)
# xlab: set the x-axis label
# ylab: set the y-axis label
# type: set the type of plot (by default it is plotted as points, type='l' plots as a line)
# HINT: you need this for Lab 1, question 1
plot(dat$Date, dat$Ta,
     type='l',
     xlab='Time',
     ylab='Air temperature (degrees C)',
     frame=FALSE)

# Add another line on top of the previous one
# col: set the color of the line
lines(dat$Date, dat$Tsoil, col='red')

# Create a legend
# First two values are the x and y coordinates for the upper-left corner of the legend
# legend: set strings to use for the legend
# col: set the color for each legend element
# lty: line type (1: straight line, 2: dashed line, etc)
legend(dat$Date[1], 10, legend=c('Air','Soil'), col=c('black','red'), lty=1)

##### 7) Plot a scatterplot
# Make a scatterplot of two variables plotted against each other
# Same basic structure as time series (above), but leave the default type
# HINT: you need this for Lab 1, question 2
plot(dat$Par, dat$Ta,
     xlab='PAR',
     ylab='Temperature',
     frame=FALSE)

##### 8) Linear regression
# Fit an ordinary least squares linear regression model (lm function) between two variables
# Use Wilkinson notation to define the model: lm(<response variable> ~ <predictor variable(s)>)
# HINT: you need this for Lab 1, question 2
mdl <- lm(dat$Ta~dat$Par)

# Display model in the console
mdl

# Display summary statistics of the model in the console
summary(mdl)

# Add the fitted regression line to the scatterplot
abline(mdl)

##### 9) Write a function to convert celcius to kelvin
# Create a customized function
# Functions take input, perform actions on that input, and then 'return' output
# HINT: you need this for Lab 1, question 4
c2k <- function(x){
  y <- x + 273.15
  return(y)
}

# Use the function to convert air temperature from celcius to kelvin
dat$Ta_kelvin <- c2k(dat$Ta)

##### 10) Indexing
# Select a subset of observations or variables using Boolean logic (AND/OR)

# Select all observations for a single day
# HINT: you need to modify this for Lab 1, question 3e
day31 <- dat[dat$Day==31, ]

# Select all observations for exactly noon
noon <- dat[dat$Hour==12 & dat$Minute==0, ]

# Select all observations of the two temperature variables
temperatures <- dat[ , c('Ta', 'Tsoil')]









