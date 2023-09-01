

#### prepare experiment ####

# create a home folder for each experiment
# there, in this folder created above, create an R-Studio project
# from within this R-Studio project, launch the aquap2 package
library(aquap2)

# create the necessary folder structure
# only needs to be called once:
# genFolderStr() 

# provide the required information in the metadata file. (metadata/metadata.r)

# aquap2 can help in creating a randomized sample list that can be used during 
#data acquisition
exportSampleList()

#### data acquisition ####
# 1) use the sample list created above (or your own sample list) when recording the 
# spectra
# 2) after data acquisition, move the rawdata file containing the spectra and a 
# possible logger file containing e.g. temperature data into the "rawdata" folder
# 3) Rename the file containing the spectra to the experiment name as given in
# the metadata (see ?gfd for further information)



#### importing raw data, plotting ####
# we assume that we already did our experiment, and we already have the acquired 
# spectra copied into the folder "rawdata" and renamed

# import raw spectra and data from the temperature and humidity data logger
fd <- gfd(slType = "xls", trhLog = "ESPEC")
# The complete dataset is saved in the folder "R-data" 
# and can, in the future, be read in simply by calling (again) "gfd":
fd <- gfd()


# we can look at the dataset:
fd
str(fd, max.level = 3)
str(fd, max.level = 4)

# plot the spectra
plot(fd, pg.where="")

# now we focus on the first overtone only
fdCut <- selectWls(fd, 1300, 1600)
plot(fdCut, pg.where="")


#### function gdmm ####
# use the function "gdmm" to produce the so called "cube", an object containing
# datasets and models according to the analysis procedure
cube <- gdmm(fd, ap=getap(spl.wl=c("800-to-2400", "1300-to-1600")))
# look at the generated cube
cube
str(cube, max.level = 4)
# get dataset #1 from the cube:
ds1 <- getcd(cube, 1)
ds1
# plot all raw spectra from a cube
plot_spectra(cube, pg.where="")

# now split the dataset by waterNames and also keep averaged consecutive scans
# instead of overwriting the (default) analysis procedure as above, we now
# use a different analysis procedure instead ("ap1.r")
cube <- gdmm(fd, ap=getap("ap1.r"))
cube
# manually add additional data pre-treatment
cube <- gdmm(fd, ap=getap("ap1.r", dpt.pre=c("snv", "sgol@2-21-1")))
cube


#### subtracting spectra & plotting ####
# from the last cube, lets subtract the averaged MilliQ from the "Ob" and the 
# "StU" water
# and plot the resulting spectra
# first get the MilliQ samples, keep only the control, and average all together
MQsFull <- getcd(cube, 1)
MQsFull
MQsFull <- ssc(MQsFull, C_Group == "Cont")
MQsFull
MQsFull_avg <- do_avg(MQsFull)
MQsFull_avg
# the same for the first overtone only, but all in one line:
MQsCut_avg <- do_avg(ssc(getcd(cube, 4), C_Group == "Cont"))
MQsCut_avg

#
# extract the "Ob" and "StU" dataset, subtract the averaged MQs samples and plot
ObFull <- getcd(cube, 8)
plot(ObFull - MQsFull_avg, pg.where="", pg.main="|Subtr: Ob-MQ|", colorBy="C_Group")
StUCut <- getcd(cube, 12)
plot(StUCut - MQsCut_avg, pg.where="", pg.main="|Subtr: StU-MQ|", colorBy="C_Group")
MQsCut <- getcd(cube, 10)
plot(MQsCut - MQsCut_avg, pg.where="", pg.main="|Subtr: MQs-MQ|", colorBy="C_Group")
ObCut <- getcd(cube, 11)
plot(ObCut - MQsCut_avg, pg.where="", pg.main="|Subtr: Ob-MQ|", colorBy="C_Group")

#
# now look at the subtraction spectra of the "Ob" water again, but this time 
# without pre-treatment:
fdCut <- selectWls(fd, 1300, 1600)
fdOb <- ssc(fdCut, C_waterNames == "Ob")
fdOb
fdMQ <- ssc(fdCut, C_waterNames=="MQs" & C_Group=="Cont")
fdMQ
fdMQ_avg <- do_avg(fdMQ)
plot(fdOb - fdMQ_avg, pg.where="", pg.main="|Subtr: Ob-MQ|", 
     pg.sub="no data pre-treatment; ", colorBy="C_Group")
plot(do_avg(fdOb, clv = "C_SampleNr") - fdMQ_avg, pg.where="", 
     pg.main="|Subtr: Ob-MQ| ConsScan avg", pg.sub="no data pre-treatment; ", 
     colorBy="C_Group")
# as a reminder, the raw spectra of the "Ob" water without subtracting the 
# averaged MQs samples
plot(fdOb, pg.where="", pg.main="|Ob|", colorBy="C_Group")


#### Make Models: PCA & PLSR ####
# we are here employing a new analysis procedure, "ap2.r", where we switch PCA 
# and PLSR to 'TRUE'.
cube <- gdmm(fd, getap("ap2.r", do.pls=FALSE))
cube
plot(cube, pg.fns="_allWaters") # is plotting to PDF in the folder "results"
# XXX crash here
#
# lets investigate consecutive scans, relative humidity and absolute time via PLSR
# 5 fold crossvalidation
cube_val5 <- gdmm(ssc(fd, C_waterNames=="Ob"), getap("ap2.r", pls.valid=5)) 
# (we overwrite "pls.valid")
cube_val5
plot_pca(cube_val5, pg.fns="_Ob")
plot_pls(cube_val5, pg.fns="_Ob_val5") # is plotting to PDF in the folder "results"
#
# class-variable based crossvalidation
cube_valCl <- gdmm(ssc(fd, C_waterNames=="Ob"), getap("ap2.r")) # no overwriting 
# of "pls.valid"
plot_pls(cube_valCl, pg.fns="_Ob_valCl")
# --> to regress on absolute time probably gives the best results



#### focus on water types ####
# lets ignore the treatment for a moment and focus on trying to differentiate 
# the different water types
# we use a new analysis procedure "ap3.r"
cube <- gdmm(ssc(fd, C_ECRM=="RM"), getap("ap3.r", aqg.mod="classic", 
                                          do.pca=T, do.sim=T))
cube
plot(cube, pg.fns="_wType")
# -> PCA shows separation between water types in comp2 vs comp3 (*)
#
# now only repeat the aquagram
tf <- "TempCalib_XDS"
cube <- gdmm(ssc(fd, C_ECRM=="RM"), getap("ap3.r", aqg.mod="aucs.dce"), 
                                  tempFile=tf)
plot_aqg(cube, pg.fns="_wType")
# subtract the MQs
cube <- gdmm(ssc(fd, C_ECRM=="RM"), getap("ap3.r",aqg.mod="aucs.dce-diff"), 
                                  tempFile=tf)
plot_aqg(cube, pg.fns="_wType")
# -> looks like this separation between water types (* above) can not be seen 
# in the aquagram



#### focus on treatment (C_Group) within each water type ####
# using "ap4,r"
cube <- gdmm(ssc(fd, C_ECRM=="RM"), getap("ap4.r", 
              aqg.mod="classic", aqg.vars="C_Group",aqg.minus="Cont",
              do.pca=T, do.sim=T))
cube
plot(cube, pg.fns="_wt_group")
#
# now only the aquagram
# please note the values provided at parameters "aqg.vars" and "aqg.minus" in 
# the file "ap4.r"
cube <- gdmm(ssc(fd, C_ECRM=="RM"), getap("ap4.r", aqg.mod="aucs.dce",
                    aqg.vars="C_Group",aqg.minus="Cont"), tempFile=tf)
plot_aqg(cube, pg.fns="_wt_group")
# subtract the MQs
cube <- gdmm(ssc(fd, C_ECRM=="RM"), getap("ap4.r", aqg.mod="aucs.dce-diff"), 
                                      tempFile=tf)
plot_aqg(cube, pg.fns="_wt_group")
# -> the "Ob" water seems to show differences between groups; 
# seen in  aquagram, but not so in PCA and SIMCA



#### focus on the "Ob" water - using EMSC ####
# we now focus only on the "Ob" water, trying to improve the separation between 
# groups (i.e. treatment) by applying EMSC, using 
# a) loading vectors and 
# b) the regression vector 
# from the regression on absolute time as done above.
# EMSC will be applied to the complete dataset
# Finally, the resulting two datasets (approach a and b) are then used in PCA, 
# SIMCA and Aquagram to differentiate, again, between groups.

#
## Approach a: using loading vector for EMSC ##
# first, organize the loading vector
cube <- gdmm(fd, getap("ap2.r", spl.var="C_waterNames", do.pls=FALSE)) 
cube
dsOb <- getcd(cube, 2) # the dataset
pcaMod_Ob <- getcm(cube, 2, "pca")
str(pcaMod_Ob)
ldOb <- pcaMod_Ob$model$loadings[, c(1,2)]
ldOb[1:10,]
matplot(ldOb, type="l") # plot to check, looks good
#
# apply EMSC
dsObEmsc_pca <- do_emsc(dsOb, ldOb)
plot(dsOb-dsObEmsc_pca, pg.where="", pg.main="|Subtr: Ob - do_emsc(Ob)", 
     pg.sub="Visualizing the effect of EMSC via PCA loading vectors") 
# we see, there is a difference between original and emsc-treated data
#
# now make PCA, SIMCA, Aquagram
cube <- gdmm(dsObEmsc_pca, getap("ap4.r", aqg.mod="classic", 
                             aqg.vars="C_Group",aqg.minus="Cont", 
                             do.pca=T, do.sim=T, spl.do.csAvg=FALSE))
plot(cube, pg.fns="_ObGrp_emscPCA")
#
# now only the aquagram
cube <- gdmm(dsObEmsc_pca, getap("ap4.r", aqg.mod="aucs.dce",
                             aqg.vars="C_Group",aqg.minus="Cont", 
                             spl.do.csAvg=FALSE), tempFile=tf)
plot_aqg(cube, pg.fns="_ObGrp_emscPCA")
# subtract the MQs
cube <- gdmm(dsObEmsc_pca, getap("ap4.r", aqg.mod="aucs.dce-diff", 
                             spl.do.csAvg=FALSE), tempFile=tf)
plot_aqg(cube, pg.fns="_ObGrp_emscPCA")


#
## Approach b: using regression vector for EMSC ##
# first, organize the regression vector (using "ap5.r")
cube <- gdmm(ssc(fd, C_waterNames=="Ob"), getap("ap5.r"))
plot(cube, pg.fns="_ObAll")
dsOb <- getcd(cube, 1) # the dataset)
# -> the probably best model is to regress on humidity (Y_RelHum) and 
# to crossvalidate based on consecutive scan
allMods <- getcm(cube, 1, "plsr")
str(allMods, max.level = 2)
unlist(allMods$regrOn); unlist(allMods$valid) # -> we want the model number 3
plsMod <- allMods$model[[3]]
str(plsMod, max.level = 2)
coefs <- plsMod$coefficients
str(coefs)
plsRegVec <- as.data.frame(plsMod$coefficients[,1,9])
matplot(plsRegVec, type="l") # check for plausibility
#
# apply EMSC
dsObEmsc_pls <- do_emsc(dsOb, plsRegVec)
plot(dsOb-dsObEmsc_pls, pg.where="", pg.main="|Subtr: Ob - do_emsc(Ob)", 
     pg.sub="Visualizing the effect of EMSC via PLSR regression vector") 
# we see, there is a difference between original and emsc-treated data
#
# now make PCA, SIMCA, Aquagram
cube <- gdmm(dsObEmsc_pls, getap("ap4.r", aqg.mod="classic", 
                                 aqg.vars="C_Group",aqg.minus="Cont", 
                                 do.pca=T, do.sim=T, spl.do.csAvg=FALSE))
plot(cube, pg.fns="_ObGrp_emscPLSR")
#
# now only the aquagram
cube <- gdmm(dsObEmsc_pls, getap("ap4.r", aqg.mod="aucs.dce",
                                 aqg.vars="C_Group",aqg.minus="Cont", 
                                 spl.do.csAvg=FALSE), tempFile=tf)
plot_aqg(cube, pg.fns="_ObGrp_emscPLSR")
# subtract the MQs
cube <- gdmm(dsObEmsc_pls, getap("ap4.r", aqg.mod="aucs.dce-diff", 
                                 spl.do.csAvg=FALSE), tempFile=tf)
plot_aqg(cube, pg.fns="_ObGrp_emscPLSR")

# -> applying EMSC via the regression vector does not help, while
# doing it via the PCA loadings did increase separation between groups



#### Leave out consecutive scans ####
# By far the best separation between groups in the "Ob" water was when 
# consecutive scans were averaged, and no EMSC was applied
# Lets look at the outcome when analyzing each consecutive scan
# in each water separately
cube <- gdmm(ssc(fd, C_ECRM=="RM"), getap("ap4.r", 
                          spl.var=c("C_conSNr", "C_waterNames"), 
                          spl.do.csAvg = FALSE, 
                          aqg.mod="aucs.dce",  
                          aqg.vars="C_Group", aqg.minus="Cont"))
cube
# for this plot, it might be good to set "aqg_alwaysPlotAvgAqg" in the 
# "settings.r" file to "FALSE".
plot_aqg(cube, pg.fns="_wt_csnr")

# -> quite good separation between groups in almost all consecutive scans
# -> how come that leaving the consScans not averaged together is so much
# -> worse then averaging them together?

# see also the presentation "Highlights", summarizing some observaations




