#==========================================================#
# Name:        Basic luminescence data analysis            #
#                   and summary statistics                 #
#                                                          #
#                                                          #
# Author:  Sourav Saha, May 2020 (Kreutzer et al., 2017)   #
#==========================================================#

##########################################################################################
##########################################################################################

#============================================================#
#  Name:    Check_Luminescence_Package                       #
#                                                            #
#  Purpose: Script to check whether Luminescence package     #
#           has been installed on this computer. If the      #
#           package is not found, then the script attempts   #
#           to download it from CRAN and install it          #
#                                                            #
# Author:   GAT Duller, November 2018                        #
#============================================================#


if(!("Luminescence"%in%installed.packages())){
  message("Luminescence Package not installed - now installing...")
  
  #First, create the library path in the Users Documents folder where they will have write access
  dir.create(Sys.getenv("R_LIBS_USER"), recursive = TRUE, showWarnings=FALSE)
  
  #Then set the repository to use to download the package
  r = getOption("repos")
  r["CRAN"] = "https://cloud.r-project.org"
  options(repos = r)
  
  #Finally undertake the installation
  install.packages("Luminescence", lib=Sys.getenv("R_LIBS_USER"))
  
  .libPaths(c(Sys.getenv("R_LIBS_USER"),.libPaths()))
  #Then activate the library
  library("Luminescence")
  
} else {
  
  message("Luminescence Package already installed. No action taken.")
  message("Checking Luminescence Package version number...")
  message("")
  library("Luminescence")
  
}


##########################################################################################
##########################################################################################
#load the data

# Setting up the directory and file input
# "Directory need to change if using in other computer"

setwd("C:/Users/path")
# Remove any incomplete rows of data where ED or ED_Err are not present
LumData<-read.csv(file.choose()) # read the csv file; # e.g., use Example1.csv

LumData <- LumData[is.finite(LumData$ED) & is.finite(LumData$ED_Err),c("ED", "ED_Err")]
print(LumData)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#             MUST CHECK for CONVERSION from Seconds to Gray            #
#               defult values may not be applicable to you              #
Al_SS <- 0.0951393245047504 # must be corrected for each month          #
SG_fact <- 1.09 # constant for the current reader                       #
dose.rate <- (Al_SS*SG_fact) # (Gy/s)                                   #
LumData <- Second2Gray(LumData, dose.rate, error.propagation = "omit")  #
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

print(LumData) # for quick look of the file to see the variable names
# Save the De values in Gy
write.csv(LumData, file = "xxxxx(Gy).csv") # change the name of the file

##########################################################################################
##########################################################################################
#==========================================================#
# Name:    Plot the kernel density estimate with           #
#          summary statistics                              #
#                                                          #
# Purpose: To see the initial distribution                 #
#                                                          #
# Authos:  Sourav Saha, May 2020 (Kreutzer et al., 2017)   #
#==========================================================#
if(!require(Luminescence)){
  library("Luminescence")
}

# Attach this data file to make accessing the fields easier
attach(LumData) 
# Creating a data frame
AnData <- data.frame(De, De.error) 

#run calculation
results <-plot_KDE(AnData,
                   cex = 1.1,
                   main = "Equivalent dose distribution",
                   xlab = "Equivalent dose [Gy]",
                   ylab = c("Density",
                            "Cumulative age counts"),
                   boxplot = TRUE,
                   summary = c("n",
                               "mean",
                               "se.abs",
                               "median",
                               "skewness",
                               "kurtosis"),
                   summary.pos = "topright",
                   col = c("blue", "orange"),
                   output = TRUE) # xlim = c(0, 200) # set the x limit for zoom

##########################################################################################
##########################################################################################

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#                              Age models:                              #
#                      1. Central Age Model (CAM)                       #
#                      2. Minimum Age Model (MAM)                       #
#                      3. Average Dose Model                            #
#                      4. Common Age Model                              #
#                      5. Finite Mixture Model (FMM)                    #
#                      6. Maximum Age Model (MAXAM)                     #
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

#===============================================================#
# Name:    Calculate logged Central Age Model (CAM)             #
#                                                               #
# Purpose: Calculates the CAM dose as proposed by               #
#          Galbraith et al. (1999) and show the results         #
# Note:    The sigmab value needs to be obtained from separate  #
#          experiments.                                         #
#                                                               #
# Authos:  Sourav Saha, May 2020 (Kreutzer et al., 2017)        #
#===============================================================#
if(!require(Luminescence)){
  library("Luminescence")
}

# Attach this data file to make accessing the fields easier
attach(LumData) 
# Creating a data frame
AnData <- data.frame(De, De.error)

#run calculation
results <- calc_CentralDose(AnData,
                            sigmab=0.15,
                            log = TRUE,
                            plot=TRUE)

#show summary
get_RLum(results)

##########################################################################################
##########################################################################################
#===============================================================#
# Name:    Calculate logged Minimum Age Model (MAM)             #
#                                                               #
# Purpose: Calculates the MAM dose as proposed by               #
#          Galbraith et al. (1999) and show the results         #
# Note:    The sigmab value needs to be obtained from separate  #
#          experiments.                                         #
#          par is the parameters (use 3 or 4 to check)          #
#                                                               #
# Authos:  Sourav Saha, May 2020 (Kreutzer et al., 2017)        #
#===============================================================#
if(!require(Luminescence)){
  library("Luminescence")
}

# Attach this data file to make accessing the fields easier
attach(LumData) 
# Creating a data frame
AnData <- data.frame(De, De.error)

#run calculation
results <- calc_MinDose(AnData,
                        sigmab=0.15,
                        log = TRUE,
                        par = 3,
                        bootstrap = FALSE,
                        level = 0.95,
                        plot = TRUE,
                        multicore = FALSE)

#show summary
get_RLum(results)

##########################################################################################
##########################################################################################
#===============================================================#
# Name:    Calculate Average Dose Model                         #
#                                                               #
# Purpose: Calculates the averages dose as proposed by          #
#          Guerin et al., 2017 and show the results             #
# Note:    The sigma_m value needs to be obtained from separate #
#          experiments.                                         #
#          Nb_BE is the #of bootstrapping (by default is 500)   #
#                                                               #
# Authos:  Sebastian Kreutzer, September 2017                   #
#===============================================================#
if(!require(Luminescence)){
  library("Luminescence")
}

# Attach this data file to make accessing the fields easier
attach(LumData) 
# Creating a data frame
AnData <- data.frame(De, De.error)

#run calculation
results <- calc_AverageDose(AnData,
                            sigma_m = 0.15,
                            Nb_BE = 10000,
                            na.rm = TRUE,
                            plot = TRUE,
                            verbose = TRUE)

#show summary
get_RLum(results)

##########################################################################################
##########################################################################################
#===============================================================#
# Name:    Calculate (un-)logged Common Age Model               #
#                                                               #
# Purpose: Calculates the Common Age dose as proposed by        #
#          Galbraith et al. (1999) and show the results         #
# Note:    The sigmab value needs to be obtained from separate  #
#          experiments.                                         #
#                                                               #
# Authos:  Sourav Saha, May 2020 (Kreutzer et al., 2017)        #
#===============================================================#
if(!require(Luminescence)){
  library("Luminescence")
}

# Attach this data file to make accessing the fields easier
attach(LumData) 
# Creating a data frame
AnData <- data.frame(De, De.error)

#run calculation
results <- calc_CommonDose(AnData,
                           sigmab=0.15,
                           log = TRUE)

#show summary
get_RLum(results)

##########################################################################################
##########################################################################################
#===============================================================#
# Name:    Calculate Finite Mixture Model (FMM)                 #
#                                                               #
# Purpose: Calculates the FMM doses as proposed by              #
#          Galbraith et al. (1999) and show the results         #
# Note:    Change the "dose.scale" to have your appropriate     #
#           dose range; here 0 to 500 years.                    #
#          Change "sigmab" for each experiment; low value for   #
#           finar resolution and high for coarse resolution.    #
#          An warning massage mean your "sigmab" value is two   #
#           low or too high for the given distribution.         #
#          n.components should be varied, e.g., 2:5 or 3:5      #
#          Lowest BIC score for k & proportion should be used   #
#                                                               #
# Authos:  Sourav Saha, May 2020 (Kreutzer et al., 2017)        #
#===============================================================#
if(!require(Luminescence)){
  library("Luminescence")
}

# Attach this data file to make accessing the fields easier
attach(LumData) 
# Creating a data frame
AnData <- data.frame(De, De.error)

#run calculation
results <- calc_FiniteMixture(AnData,
                              sigmab=0.15,
                              n.components=c(2:7),
                              grain.probability = TRUE,
                              dose.scale=c(0, 500),
                              pdf.weight = TRUE)  #init.values = list(gamma=log(4),sigma=0.2,p0=0.5,mu=log(5))

#show summary
get_RLum(results) #This is based on the max k in n.components (change it to the Lowest BIC k)

##########################################################################################
##########################################################################################
#===============================================================#
# Name:    Calculate Maximum Age Model (MAXAM)                  #
#                                                               #
# Purpose: Calculates the MAM dose as proposed by               #
#          Jon Olley and show the results                       #
# Note:    The sigmab value needs to be obtained from separate  #
#          experiments.                                         #
#          par is the parameters (use 3 or 4 to check)          #
#                                                               #
# Authos:  Sourav Saha, May 2020 (Kreutzer et al., 2017)        #
#===============================================================#
if(!require(Luminescence)){
  library("Luminescence")
}

# Attach this data file to make accessing the fields easier
attach(LumData) 
# Creating a data frame
AnData <- data.frame(De, De.error)

#run calculation
results <- calc_MaxDose(AnData,
                        sigmab=0.15,
                        log = TRUE,
                        par = 3,
                        bootstrap = FALSE,
                        level = 0.95,
                        plot = TRUE,
                        multicore = FALSE)

#show summary
get_RLum(results)

##########################################################################################
##########################################################################################

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#                  Visualisation of De distributions:                   #
#                      1. Abanico Plot (CAM)                            #
#                      2. Abanico Plot (MAM)                            #
#                      3. Radial Plot (CAM)                             #
#                      4. Radial Plot (MAM)                             #
#                      5. KDE plot                                      #
#                      6. ViolinPlot                                    #
#                      7. Histogram                                     #
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

#=========================================================#
# Name:    Plot_Abanico with CAM                          #
#                                                         #
# Purpose: To create an Abanico plot of the De data from  #
#          Analyst. A number of parameters can be changed #
#          using variables defined in this script         #
#                                                         #
# Authos:  Sourav Saha, May 2020 (Kreutzer et al., 2017)  #
#=========================================================#
if(!require(Luminescence)){
  library("Luminescence")
}

# Attach this data file to make accessing the fields easier
attach(LumData) 
# Creating a data frame
AnData <- data.frame(De, De.error)

# Calculate the CAM De for z (either logged or unlogged model)
CAM <- calc_CentralDose(AnData, log = TRUE, plot=TRUE)
# Plot now
results <- plot_AbanicoPlot(AnData,
                            log.z=TRUE,
                            z.0 = CAM$summary$de,
                            zlim = c(0, 500),
                            ylim = c(-100, 100),
                            dispersion = "qr",
                            zlab = c("Equivalent dose (Gy)"),
                            main = "Abanico plot with CAM",
                            summary = c("n",
                                        "in.2s",
                                        "mean",
                                        "se.abs",
                                        "kurtosis",
                                        "skewness"),
                            summary.pos = "topleft",
                            rotate = FALSE,
                            rug = TRUE,
                            grid.col = FALSE,
                            y.axis = FALSE,
                            line = CAM,
                            lwd=2,
                            line.col = "red",
                            line.label = "CAM",
                            cex = 1.1) # bw = 0.1 # you can add bandwith
# bar.col = FALSE # for the small bar
# error.bars = TRUE

##########################################################################################
##########################################################################################

#=========================================================#
# Name:    Plot_Abanico with MAM                          #
#                                                         #
# Purpose: To create an Abanico plot of the De data from  #
#          Analyst. A number of parameters can be changed #
#          using variables defined in this script         #
#                                                         #
# Authos:  Sourav Saha, May 2020 (Kreutzer et al., 2017)  #
#=========================================================#
if(!require(Luminescence)){
  library("Luminescence")
}

# Attach this data file to make accessing the fields easier
attach(LumData) 
# Creating a data frame
AnData <- data.frame(De, De.error)

# Calculate the MAM De for z (either logged or unlogged model)
MAM <- calc_MinDose(AnData, sigmab=0.15, log=TRUE, par=3)
# Plot now
results <- plot_AbanicoPlot(AnData,
                            log.z=TRUE,
                            z.0 = MAM$summary$de,
                            zlim = c(0, 500),
                            ylim = c(-100, 100),
                            dispersion = "qr",
                            zlab = c("Equivalent dose (Gy)"),
                            main = "Abanico plot with MAM",
                            summary = c("n",
                                        "in.2s",
                                        "mean",
                                        "se.abs",
                                        "kurtosis",
                                        "skewness"),
                            summary.pos = "topleft",
                            rotate = FALSE,
                            rug = TRUE,
                            grid.col = FALSE,
                            y.axis = FALSE,
                            line = MAM,
                            lwd=2,
                            line.col = "red",
                            line.label = "MAM",
                            cex = 1.1,
                            bw = 0.1) # error.bars = TRUE

##########################################################################################
##########################################################################################

#=========================================================#
# Name:    Radial plotter with CAM                        #
#                                                         #
# Purpose: To create a radial plot of the De data from    #
#          Analyst.                                       #
#                                                         #
# Authos:  Sourav Saha, May 2020 (Kreutzer et al., 2017)  #
#=========================================================#
if(!require(Luminescence)){
  library("Luminescence")
}

# Attach this data file to make accessing the fields easier
attach(LumData) 
# Creating a data frame
AnData <- data.frame(De, De.error)

# Calculate the CAM De for z (either logged or unlogged model)
CAM <- calc_CentralDose(AnData, log = TRUE, plot=TRUE)

# Plot now
results <- plot_RadialPlot(AnData, 
                           na.rm = TRUE,
                           log.z = TRUE,
                           #log.z = FALSE, zlim = c(0, 500), ## if limit z-axis
                           centrality = CAM$summary$de, 
                           grid.col = "none", 
                           y.ticks = TRUE,
                           output = TRUE, 
                           z.0 = CAM$summary$de, 
                           lwd = 1,
                           xlab = c("Relative standard error (%)", "Precision"),
                           ylab = "Standard estimate (%)",
                           zlab = "Equivalent dose [Gy]", 
                           rug = FALSE,
                           summary = c("n", "in.2s"),
                           cex = 1.1)

##########################################################################################
##########################################################################################

#=========================================================#
# Name:    Radial plotter with MAM                        #
#                                                         #
# Purpose: To create a radial plot of the De data from    #
#          Analyst.                                       #
#                                                         #
# Authos:  Sourav Saha, May 2020 (Kreutzer et al., 2017)  #
#=========================================================#
if(!require(Luminescence)){
  library("Luminescence")
}

# Attach this data file to make accessing the fields easier
attach(LumData) 
# Creating a data frame
AnData <- data.frame(De, De.error)

# Calculate the MAM De for z (either logged or unlogged model)
MAM <- calc_MinDose(AnData, sigmab=0.15, log=TRUE, par=3)

# Plot now
results <- plot_RadialPlot(AnData, 
                           na.rm = TRUE,
                           log.z = TRUE,
                           #log.z = FALSE, zlim = c(0, 500), ## if limit z-axis
                           centrality = MAM$summary$de, #(e.g. centrality = "median.weighted")
                           grid.col = "none", 
                           y.ticks = TRUE,
                           output = TRUE, 
                           z.0 = MAM$summary$de, 
                           lwd = 1,
                           xlab = c("Relative standard error (%)", "Precision"),
                           ylab = "Standard estimate (%)",
                           zlab = "Equivalent dose [Gy]", 
                           rug = FALSE,
                           summary = c("n", "in.2s"),
                           cex = 1.1)

##########################################################################################
##########################################################################################

#==========================================================#
# Name:    Violin plot with summary statistics             #
#                                                          #
# Purpose: To see the overall distribution                 #
#                                                          #
# Authos:  Sourav Saha, May 2020 (Kreutzer et al., 2017)   #
#==========================================================#
if(!require(Luminescence)){
  library("Luminescence")
}

# Attach this data file to make accessing the fields easier
attach(LumData) 
# Creating a data frame
AnData <- data.frame(De, De.error) 

#run calculation
results <- plot_ViolinPlot(AnData,
                           main = "Violin plot",
                           summary = c("n", "mean", "se.abs", "kurtosis", "skewness"),
                           summary.pos = "topright",
                           xlab= c("De (Gy)"),
                           cex = 1.1) 
# xlim = c(0, 200) # set the x limit for zoom

##########################################################################################
##########################################################################################

#==========================================================#
# Name:    Plot the kernel density estimate with           #
#          summary statistics                              #
#                                                          #
# Purpose: To see the overall distribution                 #
#                                                          #
# Authos:  Sourav Saha, May 2020 (Kreutzer et al., 2017)   #
#==========================================================#
if(!require(Luminescence)){
  library("Luminescence")
}

# Attach this data file to make accessing the fields easier
attach(LumData) 
# Creating a data frame
AnData <- data.frame(De, De.error) 

#run calculation
results <-plot_KDE(AnData,
                   na.rm = TRUE, 
                   values.cumulative = TRUE, 
                   order = TRUE,
                   rug = TRUE,
                   cex = 1.1,
                   main = "Equivalent dose distribution",
                   xlab = "Equivalent dose [Gy]",
                   ylab = c("Density",
                            "Cumulative age counts"),
                   boxplot = TRUE,
                   summary = c("n",
                               "mean",
                               "se.abs",
                               "median",
                               "skewness",
                               "kurtosis"),
                   summary.pos = "topright",
                   bw = 10,
                   output = TRUE) # xlim = c(0, 200) # set the x limit for zoom #bw = "nrd0"

##########################################################################################
##########################################################################################

#==========================================================#
# Name:    Plot the Histogram with   summary statistics    #
#                                                          #
# Purpose: To see the overall distribution                 #
#                                                          #
# Authos:  Sourav Saha, May 2020 (Kreutzer et al., 2017)   #
#==========================================================#
if(!require(Luminescence)){
  library("Luminescence")
}

# Attach this data file to make accessing the fields easier
attach(LumData) 
# Creating a data frame
AnData <- data.frame(De, De.error) 

# use custom breaks
#breaks <- hist(AnData$ED, breaks = 15, plot = FALSE)$breaks

#run calculation
results <- plot_Histogram(AnData,
                          rug = TRUE,
                          normal_curve = TRUE,
                          cex = 1.1,
                          pch = 1,
                          colour = c("grey", "black", "blue", "green"),
                          summary = c("n","mean.weighted","sdabs.weighted"),
                          summary.pos = "topright",
                          main = "Histogram of De",
                          mtext = "De (Gy)",
                          xlab= c("De (Gy)"),
                          ylab = c("Relative density",
                                   "Standard error"))
#breaks = breaks)

##########################################################################################
##########################################################################################

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#              Athermal Fading Correction for Feldspar grains           #
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

#==========================================================#
# Name: Calculate g-value (%/decade) and age correction    #
#          Calculate g-values for multiple-aliquots        #
# Purpose: correct ages for athermal fading                #
#                                                          #
# Authos:  Sourav Saha, May 2020 (Kreutzer et al., 2017)   #
#==========================================================#
if(!require(Luminescence)){
  library("Luminescence")
}

#==========================================================#
# keep the BIN/BINX files in the same project folder       #
#==========================================================#

## set up the directory 
setwd("C:/Users/soura/Box/Current_Working_Files/Post-Doc_UCLA/Everything in OSL/Measured IRSL samples/Mission Creek USGS Project 1/Analyzed")

## Read the Bin/BINX file
bin_fading <- choose.files(default=paste0(getwd(), "/*.*")) # e.g., SA-J1518-19_21-Suk-AA-225-BG-Fading
setwd(dirname(bin_fading))

## Read the BINX file
bin <- read_BIN2R(bin_fading)

## Select all the aliquots in the script
aliquot.position <- unique(bin@METADATA[, "POSITION"])

## Run the fading analysis (NOTE: change the data to measure either fading_data50 or fading_data225)
## Read the data for the 50  and 225 degree C
bin_De_IRSL_50 <- subset(bin, LTYPE == "IRSL"  & TEMPERATURE == 50)
bin_De_IRSL_225 <- subset(bin, LTYPE == "IRSL"  & TEMPERATURE == 225)


## For all the aliquotes in the BIN/BINX file
for(i in aliquot.position) {
  #fading_data <- Risoe.BINfileData2RLum.Analysis(bin_De_IRSL_50, pos=i)
  fading_data <- Risoe.BINfileData2RLum.Analysis(bin_De_IRSL_225, pos=i)
  
  g_value <- analyse_FadingMeasurement(
    fading_data,
    plot = TRUE,
    verbose = TRUE,
    n.MC = 1000,
    structure = c("Lx", "Tx"),
    signal.integral = c(1:5), background.integral = c(200:250))
  
}



##########################################################################################
##########################################################################################
#==========================================================#
# Name: Calculate g-value (%/decade) and age correction    #
#          Calculate g-values for single-aliquot           #
# Purpose: correct ages for athermal fading                #
#                                                          #
# Authos:  Sourav Saha, May 2020 (Kreutzer et al., 2017)   #
#==========================================================#
if(!require(Luminescence)){
  library("Luminescence")
}

#==========================================================#
# keep the BIN/BINX files in the same project folder       #
#==========================================================#

## set up the directory 
setwd("C:/Users/saha/Desktop/Fading update")

## Read the Bin/BINX file 
bin_fading <- choose.files(default=paste0(getwd(), "/*.*"), caption = "fading",)
setwd(dirname(bin_fading))

## Read the BINX file
bin <- read_BIN2R(bin_fading)

## Read the data for the 50  and 225 degree C
bin_De_IRSL_50 <- subset(bin, LTYPE == "IRSL"  & TEMPERATURE == 50)
bin_De_IRSL_225 <- subset(bin, LTYPE == "IRSL"  & TEMPERATURE == 225)

## Define the aliquot for each run 
Aliquot <- 17

## Assigning data frame of the variables for fading
fading_data50 <- Risoe.BINfileData2RLum.Analysis(bin_De_IRSL_50, pos=Aliquot)
fading_data225 <- Risoe.BINfileData2RLum.Analysis(bin_De_IRSL_225, pos=Aliquot)

## Run the fading analysis (NOTE: change the data to measure either fading_data50 or fading_data225)
g_value <- analyse_FadingMeasurement(
  fading_data225,
  plot = TRUE,
  verbose = TRUE,
  n.MC = 10000,
  structure = c("Lx", "Tx"),
  signal.integral = c(1:5),
  background.integral = c(200:250))

#g_value <- analyse_FadingMeasurement(
#fading_data50,
#plot = TRUE,
#verbose = TRUE,
#n.MC = 10000,
#structure = c("Lx", "Tx"),
#signal.integral = c(1:5),
#background.integral = c(200:250))

## Age correction using the g-value after Huntley & Lamothe, 2001

## Manually assign the ages and err
Age <- 40
Err <- 5

results <- calc_FadingCorr(age.faded = c(Age,Err),
                           g_value = g_value, 
                           n.MC = 10000,
                           txtProgressBar = TRUE,
                           seed = NULL, 
                           interval = c(0.01, 500))
#show summary
get_RLum(results)
##########################################################################################
##########################################################################################

#==========================================================#
# Name: Calculate g-value (%/decade) and age correction    #
#          Calculate g-values manually                     #
# Purpose: correct ages for athermal fading                #
#                                                          #
# Authos:  Sourav Saha, May 2020 (Kreutzer et al., 2017)   #
#==========================================================#
if(!require(Luminescence)){
  library("Luminescence")
}

## Install the lubridate package if you do not alreday have
#install.packages('lubridate')

if(!require(lubridate)){
  library("lubridate")
}

## set up the directory
setwd("C:/Users/saha/Desktop/Fading update")
## Remove any incomplete rows of data where ED or ED_Err are not present
FedData<-read.csv(file.choose()) # read the csv file; # e.g., Fadingdata1.csv

## Set the variables
TimeStart <- FedData$StartDate # when you start the measurements for each position
TimeEnd <- FedData$EndDate # when you end the measurements for each position
lxtx <- FedData$LxTx # Lx/Tx read from the analyst
lxtx_err <- FedData$LxTx_err # Lx/Tx error
Age <- FedData$age[1] # insert the DRAC age
Err <- FedData$err[1] # insert the DRAc err

## Normalize Lx/Tx and relative (Lx/Tx)error
lxtx2 <- lxtx/lxtx[1]
lxtx_err2 <- (lxtx_err/lxtx)

## Estimate the effective time difference and delay
t1 <- parse_date_time(TimeStart, '%m/%d/%y %I:%M:%S %p')
t2 <- parse_date_time(TimeEnd, '%m/%d/%y %I:%M:%S %p')
t1 <- (t1+100)-(100/exp(1))

delay <- difftime(t2, t1, units = "secs")
delay <- (as.numeric(delay))


## Assigning data frame of the variables for fading
fading_data <- data.frame(lxtx2,lxtx_err2,delay)
print(fading_data)
## use this one instead if you have >1 position
#fading_data <- data.frame(FedData[, 5:6],delay)

# Run the fading analysis
g_value <- analyse_FadingMeasurement(fading_data, 
                                     plot = TRUE,
                                     verbose = TRUE,
                                     n.MC = 10000,
                                     plot.single = FALSE)


# to correct the age according to Huntley & Lamothe, 2001
results <- calc_FadingCorr(age.faded = c(Age,Err),
                           g_value = g_value, 
                           n.MC = 10000,
                           txtProgressBar = TRUE,
                           seed = NULL, 
                           interval = c(0.01, 500))

#show summary
get_RLum(results)























##########################################################################################
##########################################################################################
# Extra unnecessary scripts (in preparation)
##########################################################################################
##########################################################################################

#==========================================================#
# Name:    CAM_KDE                                         #
#                                                          #
# Purpose: To fit the Central Age Model to the De dataset  #
#          and then create a Kernel Density plot of the    #
#          data                                            #
#                                                          #
# GAT Duller, September 2018                               #
#==========================================================#

attach(LumData)
AnData <- data.frame(De, De.error)

results <- calc_CentralDose(
  AnData,
  log = TRUE,
  plot=FALSE)

#Finally, plot a Kernel Density Estimate plot of the data
plot_KDE(AnData, cex = 1.5)

#show summary
print(as.list(results$summary), digits = 2)


##########################################################################################
#==========================================================#
# Name:    Radial plotter showing two sub-groups           #
#                                                          #
# Purpose: If you want to show two subpopulation           #
#          side-by-side                                    #
#                                                          #
#==========================================================#

attach(LumData)
AnData <- data.frame(De, De.error)
order.pop <- order(AnData$De)
AnData$De [order.pop]
AnData$De.error[order.pop]
AnData <- data.frame (AnData$De[order.pop], AnData$De.error[order.pop])
print(AnData)


## now the data set is split into sub-groups, one is manipulated
data.1 <- AnData[1:5,]
data.2 <- AnData[6:12,]
## now a common dataset is created from the two subgroups
data.3 <- list(data.1, data.2)
## now the two data sets are plotted in one plot
plot_RadialPlot(data = data.3)
## now with some graphical modification
plot_RadialPlot(data = data.3,
                na.rm = TRUE,
                log.z = TRUE,
                grid.col = "none", 
                y.ticks = TRUE,
                output = TRUE, 
                lwd = 1,
                xlab = c("Relative standard error (%)", "Precision"),
                ylab = "Standard estimate (%)",
                zlab = "Equivalent dose [Gy]", 
                rug = FALSE,
                summary = c("n", "in.2s"),
                cex = 1.1,
                col = c("darkblue", "darkgreen"),
                bar.col = c("lightblue", "lightgreen"),
                pch = c(2, 6),
                summary.pos = "sub",
                legend = c("Sample 1", "Sample 2"))




#=========================================================#
# Name:    Plot_Abanico with MAM                          #
#                                                         #
# Purpose: To create an Abanico plot of the De data from  #
#          Analyst. A number of parameters can be changed #
#          using variables defined in this script         #
#                                                         #
# Authos:  Sourav Saha, May 2020 (Kreutzer et al., 2017)  #
#=========================================================#
if(!require(Luminescence)){
  library("Luminescence")
}

# Attach this data file to make accessing the fields easier
attach(LumData) 
# Creating a data frame
AnData <- data.frame(De, De.error)

# Calculate the MAM De for z (either logged or unlogged model)
MAM <- calc_MinDose(AnData, sigmab=0.15, log=TRUE, par=3)
# Calculate FMM
FMM <- calc_FiniteMixture(AnData,
                          sigmab=0.15,
                          n.components=c(2:5),
                          grain.probability = TRUE,
                          dose.scale=c(0, 200),
                          pdf.weight = TRUE)  #init.values = list(gamma=log(4),sigma=0.2,p0=0.5,mu=log(5))

#show summary
get_RLum(FMM)
# Plot now
results <- plot_AbanicoPlot(AnData,
                            log.z=TRUE,
                            z.0 = FMM$summary$de,
                            zlim = c(0, 250),
                            ylim = c(-100, 100),
                            dispersion = "qr",
                            zlab = c("Equivalent dose (Gy)"),
                            main = "Abanico plot with MAM",
                            summary = c("n",
                                        "in.2s",
                                        "mean",
                                        "se.abs",
                                        "kurtosis",
                                        "skewness"),
                            summary.pos = "topleft",
                            rotate = FALSE,
                            rug = TRUE,
                            grid.col = FALSE,
                            y.axis = FALSE,
                            line = FMM,
                            lwd=2,
                            line.col = "red",
                            line.label = "FMM",
                            cex = 1.1,
                            bw = 0.1) # error.bars = TRUE

##########################################################################################
##########################################################################################
