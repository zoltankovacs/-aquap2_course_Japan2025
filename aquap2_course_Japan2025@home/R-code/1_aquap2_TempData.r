
# Copy the file "externalFiles/TempCalib_XDS" into the settings-home directory

# Create a home folder for the experimental analyses 
# In this created folder create an R-Studio project
# Within this R-Studio project, launch the shared codes 

# This script is intended to be worked through "manually": 
# The standard analysis procedure file ("anproc.r") is constantly modified, and the basic call to "gdmm" is always the same.
#
# Hashtag (#) is to make comments for ourself in the code



# 0.) Launch "aquap2" package in R 
library(aquap2) # load aqap2 package so we can use its functions

genFolderStr() # create the necessary folder structure ONCE!


# 1.) Creat sample list based on metadata file 
esl(form = "xls", rnd = FALSE, showFirstRows = TRUE, timeEstimate = TRUE) 

exportSampleList() # create randomised measurement order

# Use the sample list created above (or your own sample list) when recording the spectra
# After data acquisition, move the rawdata file containing the spectra and a (possible) logger file containing e.g., temperature data into the "rawdata" folder
# Rename the file containing the spectra to the experiment name as given in the metadata (see ?gfd for further information)



# 2.) Import data after experiment and copying all files on place 
fullData <- gfd(ttl = FALSE) # every detail of the spectral data is defined in the "metadata.R" file

# The complete dataset is saved in the "R-data" folder and can be read in simply by running:
fullData <- gfd()



# 3.) Plot raw spectra from imported dataset 
plot(fullData, colorBy = "C_waterTemp") 
plot_spectra(fullData, colorBy = "C_waterTemp", pg.where = "")


# 4.) Remove unwanted data from our dataset 
dataReduced <- ssc(fullData, C_Water == "AIR", include=FALSE) # use the "sub-select class" (ssc) function 
plot_spectra(dataReduced, colorBy = "C_waterTemp", pg.where = "")


# 5.) Set analysis procedure file for the range of 1300-1600 nm and plot using cube
## set in "matedata/anproc.r": 
### spl.wl <- c("1300-to-1600")
cu <- gdmm(dataReduced) # calculate PCA on the reduced dataset, then plot the results 
plot_spectra(cu, colorBy = "C_waterTemp", pg.where = "")
plot(cu, pg.fns = "_noAir") # "pg.fns" is to change the name of PDF not to overwrite previous ones 


# 6.) Do PCA on the dateset 
## set in "matedata/anproc.r": 
### do.pca <- TRUE
### pca.colorBy <- c("C_Water", "C_waterTemp")
cu2 <- gdmm(dataReduced)  # calculate PCA
plot(cu2)   # plot the results of PCA 


# 7.) Analyses in different wavelength ranges and pretreatments 
## set in "metadata/anproc.r": 
### spl.wl <- c("680-to-800", "900-to-1040", "1300-to-1600")
### do.pca <- FALSE
cu3 <- gdmm(dataReduced)
plot_spectra(cu3, colorBy = "C_waterTemp", pg.fns = "_byRange") # plot spectra

# 7a.) Set analysis procedure to do pre-treatment 
## set in "metadata/anproc.r": 
### dpt.pre <- "sgol@2-21-0"
cu3 <- gdmm(dataReduced)  
plot_spectra(cu3, colorBy = "C_waterTemp", pg.fns = "_byRange_sgol") # plot the results of PCA calculated 


# 7b.) Set analysis procedure to do pre-treatment  and calculate PCA
## set in "metadata/anproc.r": 
### dpt.post <- "msc"
### do.pca <- TRUE
cu4 <- gdmm(dataReduced)  
plot_spectra(cu4, colorBy = "C_waterTemp", pg.fns = "_byRange_sgol+msc") # plot the spectra 
plot(cu4, pg.fns = "_byRange_sgol+msc") # plot the results of PCA calculated 


# 8.) Look at on the structure of data and cube (cu) - not required for analysis 
dataReduced  # or View(dataReduced)
cu4 # or View(cu4)


# 9.) Split data by water type within the range 1300-1600 nm and run PCA 
## set in "metadata/anproc.r":
### spl.var <- "C_Water"
### spl.wl <- c("1300-to-1600")
### do.pca <- TRUE
### pca.colorBy <- c("C_waterTemp", "C_conSNr")  # colouring accordinng to diff. variables 
cu5 <- gdmm(dataReduced)
plot(cu5, pg.fns = "_byWater_sgol+msc") # plot the results of PCA 
View(cu5)  # look at the structure of the generated results 


# 10.) Do regression on water temperature with PLSR 
## set in "metadata/anproc.r":
### do.pca <- FALSE
### do.pls <- TRUE
### pls.ncomp <- 2 # set the number of latent variables (LV) to be used in the modelling 
### pls.regOn <- c("Y_waterTemp")  # variable to be predicted 
### pls.valid <- "C_SampleNr"  # default: 10-fold cross-validation but it can be changed to active class  
### pls.colorBy <- "C_waterTemp" # colouring of the points 
cu6 <- gdmm(dataReduced) # calculate PLSR 
plot(cu6, pg.fns = "_byWater_sgol+msc") # plot the results of PLSR 


# 11.) Create classic aquagrams 
## set in "metadata/anproc.r":
### do.pls <- FALSE
### do.aqg <- TRUE
### aqg.vars <- "C_waterTemp"


# 11a.) Calculate classic aquagram 
### aqg.mod <- "classic"  
cu7 <- gdmm(dataReduced)
plot(cu7, pg.fns = "_byWater") # plot "classic" aquagram 


# 11b.) set anproc to do classic aquagram by subtracting lowest temperature 
## set in "metadata/anproc.r":
### aqg.spectra <- "all"
### aqg.minus <- "20"
### aqg.mod <- "classic-diff"  
cu8 <- gdmm(dataReduced)
plot(cu8, pg.fns = "_byWater")  # plot "classic-difference" aquagram 

