
# This script is intended to be worked through without modifying the standard analysis procedure file ("anproc.r") 
#
# Hashtag (#) is to make comments for ourself in the code


# 1.) Importing data from different file formats 

# fullData <- gfd(filetype = "tabDelim.txt", trhLog = F, slType = "xls", multiplyRows = F) # from " .txt" format with existing sample list 

fullData <- gfd(getmd(expName = "LacAceNaClKCl"), filetype = "Pirouette.pir", trhLog = F, slType = NULL, multiplyRows = F) # from " .pir" format without sample list 



# 2.) Creating extra columns and renaming variables 

levels(fullData$header$C_Group) <- c("Ace", "Lac", "NaCl", "KCl")

# fullData <- reColor(fullData) # creating colour values for the whole dataset



# 3.) Preliminary investigation of the spectral data  

pca.colorBy <- c("C_ECRM", "C_Group", "C_Range", "C_conSNr") # defining colouring variables 

cu0 <- gdmm(fullData, getap(do.pca = T, pca.colorBy = pca.colorBy)) # calculating PCA 

plot_spectra(cu0, colorBy = pca.colorBy, pg.fns = "_(fullData)") # spectra plotting from "cu"
plot_pca(cu0, pg.fns = "_(fullData)") # plotting PCA from "cu"



# 4.) Reducing dataset 

fullData2 <- ssc(dataset = fullData, criteria = C_Group == "KCl", include = FALSE) # removing "KCl" data 


# 4.1. PCA on the reduced dataset 
cu01 <- gdmm(fullData2, getap(do.pca = T, pca.colorBy = pca.colorBy)) # calculating PCA by solutions 

plot_pca(cu01, pg.fns = "_(fullData2)") # plotting PCA from "cu" 


# 4.2. PCA on the reduced dataset by solutions 
pca.colorBy <- c("C_Range", "C_real_cc", "C_conSNr") # defining colouring variables 

cu02 <- gdmm(fullData2, getap(spl.var = "C_Group", do.pca = T, pca.colorBy = pca.colorBy)) # calculating PCA by solutions 

plot_pca(cu02, pg.fns = "_(fullData2_bySolution)") # plotting PCA from "cu" 



# 5.) Reducing dataset for comon concentration range 

fullData3 <- ssc(dataset = fullData2, criteria = C_Range == "2", include = TRUE) # removing "3" data 
fullData3 <- reColor(fullData3)


# 5.1. PCA on the reduced dataset 
pca.colorBy <- c("C_Group", "C_real_cc", "C_conSNr") # defining colouring variables 

cu03 <- gdmm(fullData3, getap(do.pca = T, pca.colorBy = pca.colorBy)) # calculating PCA by solutions 

plot_pca(cu03, pg.fns = "_(fullData3)") # plotting PCA from "cu" 


# 5.2. PCA on the reduced dataset by solutions 
pca.colorBy <- c("C_real_cc","C_conSNr") # defining colouring variables 

cu04 <- gdmm(fullData3, getap(spl.var = "C_Group", do.pca = T, pca.colorBy = pca.colorBy)) # calculating PCA by solutions 

plot_pca(cu04, pg.fns = "_(fullData3)_bySolution") # plotting PCA from "cu" 



# 6.) PLSR on the reduced datasets by solutions with defined nr LV
pls.ncomp <- 6 # defining number of latent variables for PLSR 

cu05 <- gdmm(fullData3, getap(spl.var = "C_Group", do.pls = T, pls.regOn = "Y_real_cc", pls.ncomp = pls.ncomp, pls.colorBy = "C_real_cc")) # calculating PLSR by solution
plot_pls(cu05, pg.fns = "_(fullData2)_bySolution_v2") # plotting PLSR from "cu" 



# 7.) Creating aquagrams 

# 7.1. Creating classic aquagrams 
cu06 <- gdmm(fullData3, getap(spl.var = "C_Group", do.aqg = T, aqg.vars = "C_real_cc", aqg.mod = "classic")) # calculating aquagrams 
plot_aqg(cu06, pg.fns = "_(fullData2)") # plotting Aquagram from "cu" 


# 7.2. Creating classic difference aquagrams 
cu07 <- gdmm(fullData3, getap(spl.var = "C_Group", do.aqg = T, aqg.vars = "C_real_cc", aqg.spectra = "all", 
                              aqg.minus = "2", aqg.mod = "classic-diff")) # calculating difference aquagrams  

plot_aqg(cu07, pg.fns = "_(fullData2)") # plotting Aquagram from "cu" 