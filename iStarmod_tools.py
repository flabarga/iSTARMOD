# -*- coding: utf-8 -*-
#=======================================================================
## Name:     STARMOD-Tools
# Filename: starmod_tools.py
# Type:     auxiliary functions for starmod
################################################################
# Language: Python 3.5
# Purpose:  Auxliary functions for the component iSTARMOD:
#            Read config files
#            Read FITS files
#            Building data extracted from the read files
#           
################################################################

# from astropy.io import fits
# import math
# import numpy as np
# import scipy.interpolate as interpl
# import matplotlib.pyplot as plt
# import copy as cpy
from scipy.optimize import*
from iStarmod_fits import *
   



def readFitsFile(inputFITS, order):

    objectFITS = FITSObject(inputFITS)
    if ".fits" in inputFITS:
        
        #Retrieving basic data of the fits format file, the headers and the proper data
        # objectFITS.fits_hdu_in.info()
        header = fits.getheader(inputFITS)
        # headerKeys = header.keys()
        count = 0
        for iterh in header.keys():
            count += 1
        numberOfCards = count #len(headerKeys)
        numberOfAxis = int(objectFITS.fits_hdu_in[0].header['NAXIS'])
        try:
            objectFITS.object      = header["OBJECT"]
        except KeyError:
            print ("KEY NOT FOUND - OBJECT")
        try:
            objectFITS.date        = float(header["MJD-OBS"])
        except KeyError:
            print ("KEY NOT FOUND - MJD-OBS")
            objectFITS.date        = float(header["BJD"])
            print ("USED BJD: ", objectFITS.date)
        try:
            numberOfPointPerSpec = int(objectFITS.fits_hdu_in[0].header['NAXIS1'])
        except KeyError:
            print ("KEY NOT FOUND - NAXIS1")
            numberOfPointPerSpec = int(-1)
        if numberOfAxis >= 2:
            try:
                numberOfOrders = int(objectFITS.fits_hdu_in[0].header['NAXIS2'])
            except KeyError:
                print ("KEY NOT FOUND - NAXIS2")
        else:
            numberOfOrders = 1
        ##FOR DEBUGGING PURPOSES
        ##----------------------------------------------------------------------------------
        #for iterh in range(numberOfCards):
        #    try:
        #        print ("[ ", headerKeys[iterh], ",",header[iterh]," ]")
        #    except fits.VerifyError:
        #        print ("error in format")
        ##----------------------------------------------------------------------------------
        
        try:
            ctype1 = objectFITS.fits_hdu_in[0].header['CTYPE1']
        except KeyError:
            print ("KEY NOT FOUND - CTYPE1")
            ctype1 = -1
        inlist = objectFITS.fits_hdu_in[0].data
        # print ("inlist: ", inlist)
        #Building the Array/s of spectrum's ordinate of wavelenghts
        if numberOfAxis == 2 and ctype1 == "MULTISPE":
            #In a multispec format, the wavelenghts are specified by means of
            #a series of parameters provided in the headers 
            # specsLambdaValues, inlist2, objectFITS.specs = 
            objectFITS.readLambdaDataMultiSpec(header,numberOfCards,numberOfPointPerSpec,numberOfOrders, inlist)
            #the way the spectrum is plotted depends on the format
            #for it in range(numberOfOrders):
            #    plt.plot(specsLambdaValues[it],inlist2[it])
            # objectFITS.specsLambdaValues = specsLambdaValues
            # objectFITS.dataValues  = inlist2 
            # objectFITS.numberOfOrders = numberOfOrders
            return objectFITS
        elif numberOfAxis == 1 and (ctype1 == "LINEAR" or ctype1 == "Linear" or ctype1 == "WAVELENGTH"):
            #in a 'onedspec'-linear format, the wavelenghts are specified by means of 
            #only one set of parameters of the linear formulae
            objectFITS.initialLambda = float(objectFITS.fits_hdu_in[0].header['CRVAL1'])
            objectFITS.deltaStep = float(objectFITS.fits_hdu_in[0].header['CDELT1'])
            print("CDELT1 = ",  objectFITS.deltaStep)
            specsLambdaValues = []
            for i in range(numberOfPointPerSpec):
                specsLambdaValues.append(objectFITS.initialLambda+i*objectFITS.deltaStep) 
            ##FOR DEBUGGING PURPOSES
            #print len(specsLambdaValues), len(inlist)
            ##----------------------------------------------------------------------------------
            objectFITS.specsLambdaValues = specsLambdaValues
            objectFITS.dataValues  = inlist 
            objectFITS.numberOfOrders = numberOfOrders
            return objectFITS
        elif numberOfAxis == 0:
            
            #######################################This is the case for a CARMENES FITS file
            objectFITS.bCarmenes = True
            objectFITS.carmenesType = len(objectFITS.fits_hdu_in)
            # print ("carmenesType = ", objectFITS.carmenesType)
            objectFITS.object      = header["OBJECT"]
            
            ########################ORIGINAL or TELLURIC CORRECTED#######################
            if objectFITS.carmenesType == 5 or objectFITS.carmenesType == 6: 
                objectFITS.specsLambdaValues = objectFITS.fits_hdu_in[4].data
                objectFITS.dataValues  = objectFITS.fits_hdu_in[1].data
                objectFITS.errorValues = objectFITS.fits_hdu_in[3].data
                objectFITS.date        = float(header["MJD-OBS"])
                try:
                    objectFITS.RVfromFITS  = header["HIERARCH CARACAL SERVAL RV"]
                except KeyError:
                    print ("KEY NOT FOUND - HIERARCH CARACAL SERVAL RV")
                try:
                    objectFITS.BaryCorr  = header["HIERARCH CARACAL BERV"]
                except KeyError:
                    print ("KEY NOT FOUND - HIERARCH CARACAL BERV")
            #####################STACKED#########################
            elif objectFITS.carmenesType == 4: 
                objectFITS.specsLambdaValues = objectFITS.fits_hdu_in[3].data
                objectFITS.dataValues  = objectFITS.fits_hdu_in[1].data
                
                for i in range(len(objectFITS.specsLambdaValues)):
                    for j in range(len(objectFITS.specsLambdaValues[i])):
                        objectFITS.specsLambdaValues[i][j] = np.exp(objectFITS.specsLambdaValues[i][j])

            objectFITS.numberOfOrders = objectFITS.specsLambdaValues.shape[0]
            # for i in range(len(objectFITS.specsLambdaValues)):
            #         print (str(i) + " = [" + str(objectFITS.specsLambdaValues[i][0]) + "," + str(objectFITS.specsLambdaValues[i][len(objectFITS.specsLambdaValues[i])-1]) + "]")
            
            # ################Normalization###############
            # norm = 0.0
            # count = 0
            # for i in range(objectFITS.numberOfOrders):
            #     for j in range(len(objectFITS.dataValues[i])):
            #         if not math.isnan(objectFITS.dataValues[i][j]) and objectFITS.dataValues[i][j] != 0.0:
            #             norm += objectFITS.dataValues[i][j]
            #             count += 1
            # norm = norm/count
            # if norm != 0.0:
            #     for i in range(objectFITS.numberOfOrders):
            #         for j in range(len(objectFITS.dataValues[i])):
            #             if not math.isnan(objectFITS.dataValues[i][j]):                        
            #                 objectFITS.dataValues[i][j] = objectFITS.dataValues[i][j]/norm
            # ############################################
            try:
                objectFITS.SNR = float(header["HIERARCH CARACAL FOX SNR " + str(order)])
            except KeyError:
                    print ("KEY NOT FOUND - " + "HIERARCH CARACAL FOX SNR " + str(order))
                    objectFITS.SNR = 167.
            return objectFITS
            
    elif ".dat" in inputFITS:
        i = 0
        with open(inputFITS, 'r') as inFile:
            for line in iter(inFile.readline,''):
                inputDataStr = line.split()
                objectFITS.specsLambdaValues.append(float(inputDataStr[0]))
                if inputDataStr[1] != "nan":
                    objectFITS.dataValues.append(float(inputDataStr[1]))
                objectFITS.numberOfOrders = 1
        objectFITS.initialLambda = objectFITS.specsLambdaValues[0]
        objectFITS.deltaStep     = objectFITS.specsLambdaValues[1] - objectFITS.specsLambdaValues[0]
        # print ("= [",objectFITS.specsLambdaValues[0],",",objectFITS.specsLambdaValues[len(objectFITS.specsLambdaValues)-1],"]")

        return objectFITS
    return
        
class spcfParams(object):
    def __init__(self, spcfFile):
        self.spcfFile = spcfFile
        self.workpath = ""
        self.objectInputFits0 = ""
        self.objectInputFits = ""
        self.writeOutputSpec = False
        self.synthSpecName = ""
        self.substractedSpecName = ""
        self.numIter =  8
        self.pixelRange = []
        self.pixelExcl = []
        self.primaryStarfits = []
        self.PRMRad = [-1.0,'fix']
        self.PRMRot = [-1.0,'fix']
        self.PRMWeight = [-1.0,'fix']
        self.secondaryStarfits = []
        self.SECRad = [-1.0,'fix']
        self.SECRot = [-1.0,'fix']
        self.SECWeight = [-1.0,'fix']
        self.tertiaryStarfits = []
        self.TERRad = [-1.0,'fix']
        self.TERRot = [-1.0,'fix']
        self.TERWeight = [-1.0,'fix']
        self.spectraFormatMode = ""
        self.spectraFormatAperture = -1
        self.spectraFormatBand = -1.0
        self.spectraFormatAperturePrimary = -1
        self.diffFormatPrimary = False
        self.cnt = 1
        self.initialLambda = 0.0
        self.finalLambda   = 0.0
        self.deltaStep     = 0.0
        self.numberOfWorkingDataValues = 0
        self.linespcf = "Halpha"
        self.refpath = ""
        self.secRefpath = ""
        self.MaxFluxDisplayed_Obj = 5.0
        self.MaxFluxDisplayed_Sub = 4.0
        self.FindMax_Tolerance = 0.001
        self.deltaStep = -1.0
        self.const = -1.0
        self.nrpoints = 0
        self.teff = 0.0
        
    def readConfigFile(self):
        
        #self.workpath = workpath
        pixelExcl = []
        for i in range(5):
            pixelExcl.append([0,0])
        #print pixelExcl
        writeOutputSpec = False
        spec = False
        i = 0
        with open(self.spcfFile, 'r') as inFile:
            for line in iter(inFile.readline,''):
                if 'IM_PATH' in line:
                    start = line.find('=') + 1
                    end   = line.find('#') - 1
                    impathStr = line[start:end].split()
                    self.workpath = impathStr[0]
                    print ("IM_PATH = ", impathStr, '\\')
                if 'OBJ_NAME' in line:
                    start = line.find('=') + 1
                    end   = line.find('#') - 1 
                    objInStr = line[start:end].split()
                    objectInputFits0 = objInStr[0]
                    objectInputFits = self.workpath + objectInputFits0
                    print("objectinputfits : ", objectInputFits)
                    self.objectInputFits0 = objectInputFits0
                    self.objectInputFits = objectInputFits 
                    #print objectInputFits
                if 'T_EFF' in line:
                    start = line.find('=') + 1
                    end   = line.find(' #')  - 1
                    teff = float(line[start:end])
                    self.teff = teff
                    print("teff = ", str(self.teff))
                if ('SYN_SPEC' in line) and spec == False: #Have we read this line yet?
                    start = line.find('=') + 1
                    end   = line.find('#') - 1 
                    if ('YES' in line[start:end]):
                        writeOutputSpec = True
                        self.writeOutputSpec = writeOutputSpec
                    spec = True
                    #print writeOutputSpec
                if 'SYN_NAME' in line:
                    start = line.find('=') + 1
                    end   = line.find('#')  - 1
                    synNameStr = line[start:end].split()
                    synthSpecName = synNameStr[0]
                    synthSpecName = self.workpath + synthSpecName
                    self.synthSpecName = synthSpecName
                    #print synthSpecName 
                if 'SUB_NAME' in line and writeOutputSpec == True:
                    start = line.find('=') + 1
                    end   = line.find('#')  - 1
                    subsSpecNameStr = line[start:end].split()
                    substractedSpecName = subsSpecNameStr[0]
                    substractedSpecName = self.workpath + substractedSpecName
                    self.substractedSpecName = substractedSpecName
                    # print (subsSpecNameStr )
                elif 'SUB_NAME' in line:
                    substractedSpecName = ''
                    self.substractedSpecName = substractedSpecName
                    #print substractedSpecName 
                if 'N_ITER' in line:
                    start = line.find('=') + 1
                    end   = line.find(' #')  - 1
                    numIter = int(line[start:end])
                    self.numIter = numIter
                    #print numIter 
                if 'PIX_ZONE' in line:
                    start = line.find('=') + 1
                    end   = line.find(' #')  - 1
                    line[start:end].split()
                    pixelRangeStr = line[start:end].split()
                    if pixelRangeStr != [] and len(pixelRangeStr)==3:
                        if 'pxl' in pixelRangeStr[2]:
                            keyw = True
                        else:   
                            keyw = False                 
                        pixelRange = [int(pixelRangeStr[0]),int(pixelRangeStr[1]),keyw]
                    elif len(pixelRangeStr) >= 2:
                        pixelRange = [int(pixelRangeStr[0]),int(pixelRangeStr[1]), True] 
                    self.pixelRange = pixelRange
                    # print ("PIX_ZONE = ", pixelRange)
                if 'PIX_EXCL' in line:
                    start = line.find('=') + 1
                    end   = line.find('#') - 1
                    pixelRangeStr = line[start:end].split()
                    if pixelRangeStr != [] and len(pixelRangeStr)==3:
                        if pixelRangeStr[2] == 'pxl':
                            keyw = True
                        else:   
                            keyw = False
                        pixelExcl[i] = [int(pixelRangeStr[0]),int(pixelRangeStr[1]), keyw]
                    elif len(pixelRangeStr) == 2:
                        pixelExcl[i] = [int(pixelRangeStr[0]),int(pixelRangeStr[1]), True]
                    # print (str(i), str(pixelExcl[i]))
                    if i<4: i+=1
                if 'DELTA_STEP' in line:
                    start = line.find('=') + 1
                    end   = line.find('#')  - 1   
                    line[start:end].split()
                    delta_step_str_ = line[start:end].split()
                    delta_step_str = delta_step_str_[0]
                    delta_step = float(delta_step_str)
                    self.deltaStep = delta_step
                    print ('Parameter delta_step = ', str(self.deltaStep)) 
                if 'NR_POINTS' in line:
                    start = line.find('=') + 1
                    end   = line.find('#')  - 1   
                    line[start:end].split()
                    nrpoints_str_ = line[start:end].split()
                    nrpoints_str = nrpoints_str_[0]
                    nrpoints = int(nrpoints_str)
                    self.nrpoints = nrpoints
                    print ('Parameter nr_points = ', str(self.nrpoints)) 
                if 'CONST_' in line:
                    start = line.find('=') + 1
                    end   = line.find('#')  - 1   
                    line[start:end].split()
                    const_str_ = line[start:end].split()
                    const_str = const_str_[0]
                    const = float(const_str)
                    self.const = const
                    print ('Parameter const_ = ', str(self.const))
                if 'PRM_NAME' in line:
                    start = line.find('=') + 1
                    end   = line.find('#')  - 1
                    prNameStr = line[start:end].split()
                    primaryStarfits = prNameStr[0]
                    if 'NONE' not in primaryStarfits and 'NONE' not in self.refpath:
                        primaryStarfits = self.refpath + primaryStarfits
                    elif 'NONE' not in primaryStarfits:
                        primaryStarfits = self.workpath + primaryStarfits
                    self.primaryStarfits = primaryStarfits
                    print ("primaryStarfits : ", primaryStarfits)
                if 'PRM_RAD' in line and primaryStarfits != '':
                    start = line.find('=') + 1
                    end   = line.find(' #')  - 1   
                    line[start:end].split()
                    PRMRadStr = line[start:end].split()
                    if PRMRadStr != []:
                        if PRMRadStr[1] == 'fix':
                            keyw = True
                        else:   
                            keyw = False                 
                        PRMRad = [float(PRMRadStr[0]),keyw]
                    self.PRMRad = PRMRad 
                    # print ('PRMRad = ',str(PRMRad))
                if 'PRM_ROT' in line and primaryStarfits != '':
                    start = line.find('=') + 1
                    end   = line.find('#')  - 1   
                    line[start:end].split()
                    # print (line[start:end])
                    PRMRotStr = line[start:end].split()
                    if PRMRotStr != []:
                        if PRMRotStr[1] == 'fix':
                            keyw = True
                        else:   
                            keyw = False                 
                        PRMRot = [float(PRMRotStr[0]),keyw]
                    self.PRMRot = PRMRot
                    # print ('PRMRot = ', str(PRMRot))
                if 'PRM_WGT' in line and primaryStarfits != '':
                    start = line.find('=') + 1
                    end   = line.find(' #')  - 1   
                    line[start:end].split()
                    PRMWeightStr = line[start:end].split()
                    if PRMWeightStr != []:
                        if PRMWeightStr[1] == 'fix':
                            keyw = True
                        else:   
                            keyw = False                 
                        PRMWeight = [float(PRMWeightStr[0]),keyw]
                    self.PRMWeight= PRMWeight
                    # print (str(PRMWeight))
                if 'SEC_NAME' in line:
                    start = line.find('=') + 1
                    end   = line.find(' #')  - 1   
                    secNameStr = line[start:end].split()
                    secondaryStarfits = secNameStr[0]
                    if 'NONE' not in secondaryStarfits and 'NONE' not in self.secRefpath:
                        secondaryStarfits = self.secRefpath + secondaryStarfits
                    elif 'NONE' not in secondaryStarfits:
                        secondaryStarfits = self.workpath + secondaryStarfits
                    self.secondaryStarfits = secondaryStarfits
                    # print (secondaryStarfits)
                if 'SEC_RAD' in line:
                    if (secondaryStarfits != '' and 'NONE' not in secondaryStarfits):
                        start = line.find('=') + 1
                        end   = line.find(' #')  - 1   
                        line[start:end].split()
                        SECRadStr = line[start:end].split()
                        if SECRadStr != []:
                            if SECRadStr[1] == 'fix':
                                keyw = True
                            else:   
                                keyw = False                 
                            SECRad = [float(SECRadStr[0]),keyw]
                        self.SECRad = SECRad 
                    # print ('SECRad = ',str(SECRad))
                    else:
                        self.SECRad = [0.0, False]
                if 'SEC_ROT' in line:
                    if (secondaryStarfits != '' and 'NONE' not in secondaryStarfits):
                        start = line.find('=') + 1
                        end   = line.find(' #')  - 1   
                        line[start:end].split()
                        SECRotStr = line[start:end].split()
                        if SECRotStr != []:
                            if SECRotStr[1] == 'fix':
                                keyw = True
                            else:   
                                keyw = False                 
                            SECRot = [float(SECRotStr[0]),keyw]
                        self.SECRot = SECRot
                        # print ('SECRot = ',str(SECRot))
                    else:
                        self.SECRot = [-1.0, False]
                if 'SEC_WGT' in line:
                    if (secondaryStarfits != '' and 'NONE' not in secondaryStarfits ):
                        start = line.find('=') + 1
                        end   = line.find(' #')  - 1   
                        line[start:end].split()
                        SECWeightStr = line[start:end].split()
                        if SECWeightStr != []:
                            if SECWeightStr[1] == 'fix':
                                keyw = True
                            else:   
                                keyw = False                 
                            SECWeight = [float(SECWeightStr[0]),keyw]
                        self.SECWeight = SECWeight
                        #print SECWeight
                    else:
                        self.SECWeight = [0.0, True]
                if 'TER_NAME' in line:
                    start = line.find('=') + 1
                    end   = line.find(' #')  - 1   
                    terNameStr = line[start:end].split()
                    tertiaryStarfits = terNameStr[0]
                    if 'NONE' not in tertiaryStarfits:
                         tertiaryStarfits = self.workpath + tertiaryStarfits
                    self.tertiaryStarfits = tertiaryStarfits
                    #print tertiaryStarfits
                if 'TER_RAD' in line:
                    if (tertiaryStarfits != '' and 'NONE' not in tertiaryStarfits):
                        start = line.find('=') + 1
                        end   = line.find(' #')  - 1   
                        line[start:end].split()
                        TERRadStr = line[start:end].split()
                        if TERRadStr != []:
                            if TERRadStr[1] == 'fix':
                                keyw = True
                            else:   
                                keyw = False                 
                            TERRad = [float(TERRadStr[0]),keyw]
                        self.TERRad = TERRad
                        #print TERRad
                    else:
                        self.TERRad = [0.0, False]
                if 'TER_ROT' in line:
                    if (tertiaryStarfits != '' and 'NONE' not in tertiaryStarfits):
                        start = line.find('=') + 1
                        end   = line.find(' #')  - 1   
                        line[start:end].split()
                        TERRotStr = line[start:end].split()
                        if TERRotStr != []:
                            if TERRotStr[1] == 'fix':
                                keyw = True
                            else:   
                                keyw = False                 
                            TERRot = [float(TERRotStr[0]),keyw]
                        self.TERRot = TERRot
                        #print TERRot
                    else:
                        self.TERRot = [-1.0, True]
                if 'TER_WGT' in line:
                    if (tertiaryStarfits != '' and 'NONE' not in tertiaryStarfits):
                        start = line.find('=') + 1
                        end   = line.find(' #')  - 1   
                        line[start:end].split()
                        TERWeightStr = line[start:end].split()
                        if TERWeightStr != []:
                            if TERWeightStr[1] == 'fix':
                                keyw = True
                            else:   
                                keyw = False                 
                            TERWeight = [float(TERWeightStr[0]),keyw]
                        self.TERWeight = TERWeight 
                        #print TERWeight 
                    else:
                        self.TERWeight = [-1.0, True]
                if 'MODE' in line:
                    start = line.find('= ') + 2
                    end   = line.find(' #')  - 1
                    if 'ech' in line[start:end]:
                        spectraFormatMode = 'echelon'
                    elif 'mult' in line[start:end]:
                        spectraFormatMode = 'multispec'
                    self.spectraFormatMode = spectraFormatMode
                    #print spectraFormatMode
               
                if 'APERTURE' in line:
                    start = line.find('=') + 1
                    end   = line.find(' #')  - 1   
                    spectraFormatAperture = int(line[start:end])
                    self.spectraFormatAperture = spectraFormatAperture
                    #print spectraFormatAperture
                if 'APERT_PRIMARY' in line:
                    start = line.find('=') + 1
                    end   = line.find(' #')  - 1   
                    if ('NONE' not in line[start:end]):
                        strSpectraFormatAperturePrimary = int(line[start:end])
                        self.spectraFormatAperturePrimary = strSpectraFormatAperturePrimary
                        self.diffFormatPrimary = True
                    #print spectraFormatAperture
                if 'BAND' in line:
                    start = line.find('=') + 1
                    end   = line.find(' #')  - 1
                    if line[start:end].isdigit():
                        spectraFormatBand = int(line[start:end])
                    else:
                        spectraFormatBand = -1
                    self.spectraFormatBand = spectraFormatBand
                    #print spectraFormatBand
                if 'LINE' in line:
                    start = line.find('=') + 1
                    end   = line.find('#')  - 1
                    linespcfStr = line[start:end].split()
                    print (start, end, linespcfStr)
                    self.linespcf = linespcfStr[0]
                if 'REF_PATH' in line:
                    start = line.find('=') + 1
                    end   = line.find('#') - 1
                    refpathStr = line[start:end].split()
                    self.refpath = refpathStr[0]
                    # print ("REF_PATH = ", refpathStr, '\\' )
                if 'SEC_PATH' in line:
                    start = line.find('=') + 1
                    end   = line.find('#') - 1
                    secRefpathStr = line[start:end].split()
                    self.secRefpath = secRefpathStr[0]
                
                if 'MAXFLUXDISP_OBJ' in line:
                    start = line.find('=') + 1
                    end   = line.find('#')  - 1   
                    line[start:end].split()
                    # print (line[start:end])
                    MaxFluxDisplayed_Obj_Str_ = line[start:end].split()
                    MaxFluxDisplayed_Obj_Str = MaxFluxDisplayed_Obj_Str_[0]
                    if MaxFluxDisplayed_Obj_Str_ != []:
                        MaxFluxDisplayed_Obj = float(MaxFluxDisplayed_Obj_Str)
                    self.MaxFluxDisplayed_Obj = MaxFluxDisplayed_Obj
                    print ('Max. Flux to Display = ', str(self.MaxFluxDisplayed_Obj))
                if 'MAXFLUXDISP_SUB' in line:
                    start = line.find('=') + 1
                    end   = line.find('#')  - 1   
                    line[start:end].split()
                    # print (line[start:end])
                    MaxFluxDisplayed_Sub_Str_ = line[start:end].split()
                    MaxFluxDisplayed_Sub_Str = MaxFluxDisplayed_Sub_Str_[0]
                    if MaxFluxDisplayed_Sub_Str_ != []:
                        MaxFluxDisplayed_Sub = float(MaxFluxDisplayed_Sub_Str)
                    self.MaxFluxDisplayed_Sub = MaxFluxDisplayed_Sub
                    print ('Max. Flux to Display (Subtracted) = ', str(self.MaxFluxDisplayed_Sub))
                if 'FINDMAX_TOLERANCE' in line:
                    start = line.find('=') + 1
                    end   = line.find('#')  - 1   
                    line[start:end].split()
                    # print (line[start:end])
                    FindMax_Tolerance_Str_ = line[start:end].split()
                    FindMax_Tolerance_Str = FindMax_Tolerance_Str_[0]
                    if FindMax_Tolerance_Str_ != []:
                        FindMax_Tolerance = float(FindMax_Tolerance_Str)
                    self.FindMax_Tolerance = FindMax_Tolerance
                    print ('Parameter Range for finding Max = ', str(self.FindMax_Tolerance))

        self.pixelExcl = pixelExcl
        # print (str(pixelExcl))
        
    def readConfigFileIter(self):
        with open(self.spcfFile, 'r') as inFile:
            for line in iter(inFile.readline,''):
                if 'PRM_NAME' in line:
                    start = line.find('=') + 1
                    end   = line.find('#')  - 1
                    prNameStr = line[start:end].split()
                    primaryStarfits = prNameStr[0]
                    if 'NONE' not in primaryStarfits and 'NONE' not in self.refpath:
                        primaryStarfits = self.refpath + primaryStarfits
                    elif 'NONE' not in primaryStarfits:
                        primaryStarfits = self.workpath + primaryStarfits
                    self.primaryStarfits = primaryStarfits
                    print ("primaryStarfits : ", primaryStarfits)
                if 'PRM_RAD' in line and primaryStarfits != '':
                    start = line.find('=') + 1
                    end   = line.find(' #')  - 1   
                    line[start:end].split()
                    PRMRadStr = line[start:end].split()
                    if PRMRadStr != []:
                        if PRMRadStr[1] == 'fix':
                            keyw = True
                        else:   
                            keyw = False                 
                        PRMRad = [float(PRMRadStr[0]),keyw]
                    self.PRMRad = PRMRad 
                if 'PRM_ROT' in line and primaryStarfits != '':
                    start = line.find('=') + 1
                    end   = line.find('#')  - 1   
                    line[start:end].split()
                    PRMRotStr = line[start:end].split()
                    if PRMRotStr != []:
                        if PRMRotStr[1] == 'fix':
                            keyw = True
                        else:   
                            keyw = False                 
                        PRMRot = [float(PRMRotStr[0]),keyw]
                    self.PRMRot = PRMRot
                if 'PRM_WGT' in line and primaryStarfits != '':
                    start = line.find('=') + 1
                    end   = line.find(' #')  - 1   
                    line[start:end].split()
                    PRMWeightStr = line[start:end].split()
                    if PRMWeightStr != []:
                        if PRMWeightStr[1] == 'fix':
                            keyw = True
                        else:   
                            keyw = False                 
                        PRMWeight = [float(PRMWeightStr[0]),keyw]
                    self.PRMWeight= PRMWeight
                if 'SEC_NAME' in line:
                    start = line.find('=') + 1
                    end   = line.find(' #')  - 1   
                    secNameStr = line[start:end].split()
                    secondaryStarfits = secNameStr[0]
                    if 'NONE' not in secondaryStarfits:
                        secondaryStarfits = self.workpath + secondaryStarfits
                    self.secondaryStarfits = secondaryStarfits
                if 'SEC_RAD' in line and (secondaryStarfits != '' and 'NONE' not in secondaryStarfits):
                    start = line.find('=') + 1
                    end   = line.find(' #')  - 1   
                    line[start:end].split()
                    SECRadStr = line[start:end].split()
                    if SECRadStr != []:
                        if SECRadStr[1] == 'fix':
                            keyw = True
                        else:   
                            keyw = False                 
                        SECRad = [float(SECRadStr[0]),keyw]
                    self.SECRad = SECRad 
                if 'SEC_ROT' in line and (secondaryStarfits != '' and 'NONE' not in secondaryStarfits):
                    start = line.find('=') + 1
                    end   = line.find(' #')  - 1   
                    line[start:end].split()
                    SECRotStr = line[start:end].split()
                    if SECRotStr != []:
                        if SECRotStr[1] == 'fix':
                            keyw = True
                        else:   
                            keyw = False                 
                        SECRot = [float(SECRotStr[0]),keyw]
                    self.SECRot = SECRot
                if 'SEC_WGT' in line and (secondaryStarfits != '' and 'NONE' not in secondaryStarfits ):
                     start = line.find('=') + 1
                     end   = line.find(' #')  - 1   
                     line[start:end].split()
                     SECWeightStr = line[start:end].split()
                     if SECWeightStr != []:
                         if SECWeightStr[1] == 'fix':
                             keyw = True
                         else:   
                             keyw = False                 
                         SECWeight = [float(SECWeightStr[0]),keyw]
                     self.SECWeight = SECWeight
                if 'TER_NAME' in line:
                    start = line.find('=') + 1
                    end   = line.find(' #')  - 1   
                    terNameStr = line[start:end].split()
                    tertiaryStarfits = terNameStr[0]
                    if 'NONE' not in tertiaryStarfits:
                         tertiaryStarfits = self.workpath + tertiaryStarfits
                    self.tertiaryStarfits = tertiaryStarfits
                if 'TER_RAD' in line and (tertiaryStarfits != '' and 'NONE' not in tertiaryStarfits):
                    start = line.find('=') + 1
                    end   = line.find(' #')  - 1   
                    line[start:end].split()
                    TERRadStr = line[start:end].split()
                    if TERRadStr != []:
                        if TERRadStr[1] == 'fix':
                            keyw = True
                        else:   
                            keyw = False                 
                        TERRad = [float(TERRadStr[0]),keyw]
                    self.TERRad = TERRad
                if 'TER_ROT' in line and (tertiaryStarfits != '' and 'NONE' not in tertiaryStarfits):
                    start = line.find('=') + 1
                    end   = line.find(' #')  - 1   
                    line[start:end].split()
                    TERRotStr = line[start:end].split()
                    if TERRotStr != []:
                        if TERRotStr[1] == 'fix':
                            keyw = True
                        else:   
                            keyw = False                 
                        TERRot = [float(TERRotStr[0]),keyw]
                    self.TERRot = TERRot
                if 'TER_WGT' in line and (tertiaryStarfits != '' and 'NONE' not in tertiaryStarfits):
                    start = line.find('=') + 1
                    end   = line.find(' #')  - 1   
                    line[start:end].split()
                    TERWeightStr = line[start:end].split()
                    if TERWeightStr != []:
                        if TERWeightStr[1] == 'fix':
                            keyw = True
                        else:   
                            keyw = False                 
                        TERWeight = [float(TERWeightStr[0]),keyw]
                    self.TERWeight = TERWeight 
                if 'APERTURE' in line:
                    start = line.find('=') + 1
                    end   = line.find(' #')  - 1   
                    spectraFormatAperture = int(line[start:end])
                    self.spectraFormatAperture = spectraFormatAperture
                if 'LINE' in line:
                    start = line.find('=') + 1
                    end   = line.find('#')  - 1
                    linespcfStr = line[start:end].split()
                    self.linespcf = linespcfStr[0]

    def checkWeights(self):
        if (self.PRMWeight[0] > 1.0):
            self.PRMWeight[0] = 1.0
        self.TERWeight[0] = 1.0 - self.PRMWeight[0]  - self.SECWeight[0] 
        self.TERWeight[0] = float('{:5.3f}'.format(self.TERWeight[0]))
        if (self.TERWeight[0] < 0.0):
            self.TERWeight[0] = 0.0
        self.SECWeight[0] = 1.00 - self.PRMWeight[0]
        self.SECWeight[0] = float('{:5.3f}'.format(self.SECWeight[0]))
        return 

def readConfigFile(spcfFile):
    configParams = spcfParams(spcfFile)
    configParams.readConfigFile()
    return configParams

class lambdaData(object):
    def __init__(self,spcfFile):
        self.spcfFile = spcfFile
        self.speedOfLight = 299792.458
        self.vRad = 0.0
        self.ldoLineStart = 0.0
        self.ldoLineEnd = 0.0
        self.ldoObs = 0.0
        self.ldoObs2 = 0.0
        self.ldoObs3 = 0.0
        self.ldoLineInput = []
        self.start = 0
        self.end = 0
        self.line_input = 'Halpha'
        self.lambdaDataDict = self.readConfigFile()    
        return
        
    def readConfigFile(self):
        self.lambdaDataDict = {}
        with open(self.spcfFile, 'r') as inFile:
            for line in iter(inFile.readline,''):
                end = line.find(':=') -1
                key = line[:end]
                start1 = line.find("lambda = ") + 9
                end1   = line.find("lambda_left")
                # data1 = float(line[start1:end1])
                # start2 = line.find("deltalambda =") + 13
                if end1 == -1:
                    end1 = line.find("deltalambda =")
                    data1 = float(line[start1:end1-1])
                    data2 = float(line[end1+13:])
                else:
                    start2 = line.find("deltalambda =") + 13
                    data1 = float(line[start1:end1])
                    data2 = float(line[start2:])
                self.lambdaDataDict.update( {key :[data1,data2]})
            # print (str(self.lambdaDataDict))
        return self.lambdaDataDict
        
    def calculateDopplerDisplacement(self, vRad, nstar2 = False, vRad2 = 0.0):

        displacementFactor = 1 + float(vRad/self.speedOfLight)
        self.ldoObs = self.ldoLineInput[0] * displacementFactor
        if nstar2:
            displacementFactor = 1 + float(vRad2/self.speedOfLight)
            self.ldoObs2 = self.ldoLineInput[0] * displacementFactor
        return self.ldoObs, self.ldoObs2
        
    def calculateRadialVelocityfromDopplerDisplacement(self):

        self.vRad = (((self.ldoObs/self.ldoLineInput[0])-1)*self.speedOfLight)
        return self.vRad
        
    def setLineStart(self, fLambda):

        self.ldoLineStart = fLambda
        return self.ldoLineStart
        
    def setLineEnd(self, fLambda):

        self.ldoLineEnd = fLambda
        return self.ldoLineEnd
        
    def setLineInputbyKey(self, key):
        self.line_input = key
        self.ldoLineInput = self.lambdaDataDict[key]
        return self.ldoLineInput
    
    def get_line_input(self):
        return self.line_input
    
    def setLineInputbyValue(self,value):
        self.ldoLineInput[0] = float(value)
        return self.ldoLineInput
             
    def selectLambdaRanges(self,workingLambdaValues, deltaStep):
        startPassed = False
        endPassed = False
        for iterLdo in range(len(workingLambdaValues)):
            if (self.ldoLineStart - workingLambdaValues[iterLdo]) < 0.1*deltaStep and startPassed == False:
                self.start = iterLdo
                startPassed = True
                # print (self.ldoLineStart)
            if (self.ldoLineEnd - workingLambdaValues[iterLdo] ) < 0.1*deltaStep and endPassed == False:
                self.end = iterLdo
                endPassed = True
                # print (self.ldoLineEnd)
            if(endPassed == False):
                self.end = iterLdo
        return 
    
    

    def find_local_max(self, x, y, x_center, delta, order = 1, abs_max = 0.0):

        from scipy.signal import find_peaks
        xx = np.array(x,float)
        yy = np.array(y,float)

        mask = (xx >= x_center - delta) & (xx <= x_center + delta)
        x_region = xx[mask]
        y_region = yy[mask]
        
        if len(y_region) == 0:
            return None, None, None

        peaks, _ = find_peaks(y_region, distance = 1.5)
        if len(peaks) == 0:
            return None, None, None

        # Find the peak with the highest y value
        if (order == 1):
            max_peak_idx = peaks[np.argmax(y_region[peaks])]
            return x_region[max_peak_idx], y_region[max_peak_idx], x_region
        ## Or the second highest...
        if (order == 2):
            second_max = abs_max
            if self.get_line_input()  == 'Halpha':
                second_max_distance_parameter = 1.5
            elif self.get_line_input() == 'CaIRT-a' or self.get_line_input() == 'HeD3':
                second_max_distance_parameter = 0.5
            second_max_distance = math.fabs(second_max - abs_max)
            while second_max_distance < second_max_distance_parameter: ## FREE PARAMETER second_max_distance. Depends on the specific set of star spectra
                xx_region = np.array([], dtype='float64')
                yy_region = np.array([], dtype='float64')
                for it in peaks:
                    if (x_region[it] != second_max):
                        xx_region = np.append(xx_region, x_region[it])
                        yy_region = np.append(yy_region,y_region[it])
                max_peak_idx = np.argmax(yy_region)
                second_max = xx_region[max_peak_idx]
                second_max_distance = math.fabs(second_max - abs_max)
                peaks = tuple(it for it in peaks if x_region[it] != second_max)
            max_peak_idx = np.argmax(yy_region)

            return xx_region[max_peak_idx], yy_region[max_peak_idx], x_region
    
    def getMaximumValue(self,fitsObject, linespcf, findMaxTolerance, nstars = 1, inputValue2 = 0.0):
        ## Once the maximum value has been calculated, the initial value of self.ldoObs is modified 
        ## with the value of the lambda for this maximum value

        max = -1.0
        # if (nstars != 1):
        #     findMaxTolerance = findMaxTolerance/10
        factorStart = 1-findMaxTolerance
        factorEnd   = 1+findMaxTolerance
        ##############################################################################
        debugging_value = 0.0
        workingLambdas = np.array(fitsObject.workingLambdaValues, dtype = 'float64')
        workingfluxes  = np.array(fitsObject.workingDataValues3, dtype = 'float64')
        _start_ = factorStart*self.ldoObs
        _end_   = factorEnd*self.ldoObs
        delta1   = self.ldoObs - _start_ 
        delta2  = _end_ - self.ldoObs
        
        if nstars == 1:

            print("self.ldoObs= ", self.ldoObs, "lambda_factorStart = ", _start_, "lambda_factorEnd= ", _end_)
           ######################################################################################
            mask = (workingLambdas >= self.ldoObs - delta1) & (workingLambdas <= self.ldoObs + delta2)
            lambdas_region = workingLambdas[mask]
            fluxes_region = workingfluxes[mask]
            abs_max_idx    = np.argmax(fluxes_region)
            self.ldoObs = lambdas_region[abs_max_idx]
            abs_max_flux = fluxes_region[abs_max_idx]
            ####################################################################################y
            # iterMax = len(fitsObject.workingDataValues3)-1
            # iterMax2 = len(fitsObject.workingDataValues3)-1
            # iterM = -1
            # for iterM in range(len(fitsObject.workingDataValues3)):
            #     debugging_lambda = fitsObject.workingLambdaValues[iterM]
            #     debugging_value  = math.fabs(fitsObject.workingDataValues3[iterM])
            #     if debugging_lambda < _start_:
            #         continue
            #     if debugging_value > max:
            #         max = math.fabs(fitsObject.workingDataValues3[iterM]) 
            #         iterMax = iterM
            #     if ((nstars == 1 and debugging_lambda > _end_) or (nstars > 1 and debugging_lambda > (self.ldoObs+self.ldoObs2)/2)):
            #         break
            # debugging_lambda = fitsObject.workingLambdaValues[iterMax]
            # debugging_value  = fitsObject.workingDataValues3[iterMax]
            
            # max = -1

        elif nstars == 2:
            # for debugging purposes
            mask = (workingLambdas >= self.ldoObs - delta1) & (workingLambdas <= self.ldoObs + delta1)
            lambdas_region = workingLambdas[mask]
            fluxes_region  = workingfluxes[mask]
            abs_max_idx    = np.argmax(fluxes_region)
            abs_max_lambda = lambdas_region[abs_max_idx]
            abs_max_flux   = fluxes_region[abs_max_idx]
            debugging_lambda, debugging_value, _ = self.find_local_max(workingLambdas, workingfluxes, self.ldoObs, delta1, order = 2, abs_max = abs_max_lambda)
            
            
            # _start_ = factorStart*self.ldoObs2
            # _end_   = factorEnd*self.ldoObs2
            # delta2 = _end_ -self.ldoObs2
            # debugging_lambda_arr, debugging_value_arr, _ = self.find_local_max(fitsObject.workingLambdaValues, fitsObject.workingDataValues3, self.ldoObs2, math.fabs(_start_-self.ldoObs2))
            
            self.ldoObs    = abs_max_lambda
            self.ldoObs2   = debugging_lambda

            # for iterM2 in range(len(fitsObject.workingDataValues3)):
            #     if fitsObject.workingLambdaValues[iterM2] < _start_:
            #         continue
            #     if math.fabs(fitsObject.workingDataValues3[iterM2]) > max:
            #         # print(fitsObject.workingLambdaValues[iterM])
            #         max = math.fabs(fitsObject.workingDataValues3[iterM2]) 
            #         iterMax2 = iterM2
            #     if fitsObject.workingLambdaValues[iterM2] > _end_:
            #         break
            #     debugging_lambda = fitsObject.workingLambdaValues[iterM2]
            # self.ldoObs  = fitsObject.workingLambdaValues[iterMax]
            # self.ldoObs2 = fitsObject.workingLambdaValues[iterMax2]


        return self.ldoObs, abs_max_flux, _start_, _end_, self.ldoObs2, debugging_value

    def lorentzian(self, x, amp, mean, gamma):
        return amp * gamma**2 / ((x - mean)**2 + gamma**2)
    
    def combined_model_lorentzian(self, x, amp1, mean1, gamma1, amp2, mean2, gamma2):
        return (self.lorentzian(x, amp1, mean1, gamma1) + self.lorentzian(x, amp2, mean2, gamma2))
    
    def gaussian(self, x, amp, mean, sigma):
        return (amp * np.exp(-(x - mean) ** 2 / (2 * sigma ** 2)))

    def combined_model_gaussian(self, x, amp1, mean1, sigma1, amp2, mean2, sigma2):
        return (self.gaussian(x, amp1, mean1, sigma1) + self.gaussian(x, amp2, mean2, sigma2))
    
    def paint_diff(self, xdata, ydata, summ_up, popt):
        
        from scipy.integrate import simps, quad
        # Datos y modelo
        amp1_fitted, mean1_fitted, gamma1_fitted, amp2_fitted, mean2_fitted, gamma2_fitted = popt
        x_fit = np.linspace(xdata[0], xdata[-1], len(summ_up))
        y_fit = self.combined_model_lorentzian(x_fit, *popt) + 1.0

        # Integral de los datos
        int_data = simps(ydata, xdata)

        # Integral del modelo
        int_model = quad(lambda x: self.combined_model_lorentzian(x, *popt) + 1, xdata[0], xdata[-1])[0]

        # Gráfico
        g = plt.figure(figsize=(10, 5), dpi=150)

        # Área bajo la curva de los datos
        plt.fill_between(xdata, ydata, alpha=0.4, label=f"Datos (Área = {int_data:.2f})")

        # Área bajo la curva del modelo
        plt.fill_between(x_fit, y_fit, alpha=0.4, label=f"Modelo (Área = {int_model:.2f})")

        # Curvas
        plt.plot(xdata, ydata, 'o', label="Datos", color='black')
        plt.plot(x_fit, y_fit, '-', label="Modelo ajustado", color='red')

        plt.xlabel("x")
        plt.ylabel("y")
        plt.title("Comparación de áreas bajo la curva")
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.show()


    def CalculateEW(self,fitsObject, fitsobject2 = None, initial_weight = 0.0, initial_weight2 = 0.0, nstar = 1):
        
        from scipy.integrate import simps, quad

        fDataValues   = cpy.copy(fitsObject.workingDataValues3[self.start:self.end])
        fLambdaValues = cpy.copy(fitsObject.workingLambdaValues[self.start:self.end])
        deltaLambda = fLambdaValues[len(fDataValues)-1] - fLambdaValues[0]
        print (len(fDataValues), self.start, self.end)
        # print("len(fLambdaValues)=", len(fLambdaValues))
        if nstar == 1:
            for it in range(len(fDataValues)):
                fDataValues[it] = fDataValues[it] + 1.0
                
            EW = - (np.trapz(fDataValues,fLambdaValues) - 
                (fDataValues[0] + fDataValues[len(fDataValues)-1]) * deltaLambda / 2) #Area of trapezium below pushed-up continuum 
            return EW, None, None, None, None, None, None, None
        
        if nstar == 2:
            # initial guess = [weight1, _ldo_guess_1, sigma1, weight2, _ldo_guess_2, sigma2]
            if(self.ldoObs>=self.ldoObs2):
                if self.get_line_input() == 'Halpha':
                    initial_guess =  [initial_weight, self.ldoObs, 2.0, initial_weight2, self.ldoObs2, 1.0]
                elif self.get_line_input() == 'HeD3':
                    initial_guess =  [initial_weight, self.ldoObs, 0.5, initial_weight2, self.ldoObs2, 0.1]
                elif self.get_line_input() == 'CaIRT-a':
                    initial_guess =  [initial_weight, self.ldoObs, 1.0, initial_weight2, self.ldoObs2, 0.5]
            else:
                if self.get_line_input() == 'Halpha':
                    initial_guess =  [initial_weight, self.ldoObs, 2.5, initial_weight2, self.ldoObs2, 0.10]
                elif self.get_line_input() == 'CaIRT-a':
                    initial_guess =  [initial_weight, self.ldoObs, 1.0, initial_weight2, self.ldoObs2, 1.0]
            # Fit the model to the data
            try:
                popt, _ = curve_fit(self.combined_model_lorentzian, fLambdaValues, fDataValues, p0=initial_guess)
                # popt, _ = curve_fit(self.combined_model_gaussian, fLambdaValues, fDataValues, p0=initial_guess)
                # perr = np.sqrt(np.diag(pcov))
                # print (infodict('fvec'))

                # Extract fitted parameters
                amp1_fitted, mean1_fitted, gamma1_fitted, amp2_fitted, mean2_fitted, gamma2_fitted = popt

                # Generate the fitted model
                # fitted_spectrum = self.combined_model_lorentzian(fLambdaValues, *popt)

                # Generate individual spectra
                star1_spectrum = self.lorentzian(fLambdaValues, amp1_fitted, mean1_fitted, gamma1_fitted)
                # star1_spectrum = self.gaussian(fLambdaValues, amp1_fitted, mean1_fitted, gamma1_fitted)
                # print("len(star1_spectrum) = ", len(star1_spectrum))
                star2_spectrum = self.lorentzian(fLambdaValues, amp2_fitted, mean2_fitted, gamma2_fitted)
                # star2_spectrum = self.gaussian(fLambdaValues, amp2_fitted, mean2_fitted, gamma2_fitted)
                # print("len(star2_spectrum) = ", len(star1_spectrum))
                
                star1_spectrum_trapz  = []
                star2_spectrum_trapz  = []
                summ_up               = []
                summ_up_trapz         = []
                for it in range(len(star1_spectrum)):
                    star1_spectrum_trapz.append(star1_spectrum[it] + 1.0)
                    star2_spectrum_trapz.append(star2_spectrum[it] + 1.0)
                    summ_up.append(star1_spectrum[it] + star2_spectrum[it])
                    summ_up_trapz.append(star1_spectrum[it] + star2_spectrum[it] + 1.0)
                    fDataValues[it] = fDataValues[it] + 1.0

                ### Calculation of the error ##############################################################
                # self.paint_diff (fLambdaValues, fDataValues, summ_up, popt)
                int_data = simps(fDataValues, fLambdaValues) # Integral of the original data 
                int_model1 = quad(lambda x: self.combined_model_lorentzian(x, *popt) + 1.0, fLambdaValues[0], fLambdaValues[-1])[0]# Integral of the adjusted model
                int_model2 = simps(summ_up_trapz,fLambdaValues)
                abs_error = int_model2 - int_data # Absolute error
                rel_error = 100 * abs(abs_error) / abs(int_data) # Relative error 
                print("int_data = ", int_data, ";  int_model1 = ", int_model1, ";  int_model2 = ", int_model2)
                print("abs_error = ", abs_error, "  relative: ", '{:2.1f}'.format(rel_error) + '%') 
                #############################################################################################
                
                EW1 = - (np.trapz(star1_spectrum_trapz,fLambdaValues) - 
                    (star1_spectrum_trapz[0] + star1_spectrum_trapz[len(star1_spectrum_trapz)-1]) * deltaLambda / 2) #Area of trapezium below pushed-up continuum 
                EW2 = - (np.trapz(star2_spectrum_trapz,fLambdaValues) - 
                    (star2_spectrum_trapz[0] + star2_spectrum_trapz[len(star2_spectrum_trapz)-1]) * deltaLambda / 2) #Area of trapezium below pushed-up continuum 

                return EW1, EW2, star1_spectrum, star2_spectrum, fLambdaValues,summ_up, abs_error, rel_error
            except RuntimeError:
                print("Optimal parameters not found: Number of calls to function has reached maxfev = 1400.")
                fDataValues2 = cpy.copy(fDataValues)
                for it in range(len(fDataValues2)):
                    fDataValues2[it] = fDataValues2[it] + 1.0
                
                EW = - (np.trapz(fDataValues2,fLambdaValues) - 
                    (fDataValues2[0] + fDataValues2[len(fDataValues)-1]) * deltaLambda / 2) 
                #Area of trapezium below pushed-up continuum 
                return EW, None, fDataValues, fDataValues, fLambdaValues, None, None, None
    
    
    def getFWHM(self,fitsObject):
        startFound = False
        maxIter = -1.0
        minIter = 0.0
        fLambdaValues, fDataValues = self.resampleFWHMDataInterval(fitsObject)
        endLambdaValue = 0.0
        startLambdaValue = -1.0
        iterMax = 1.0
        for iter in range(len(fLambdaValues)):
            itemToCompare = fDataValues[iter]
            if (itemToCompare > maxIter):
                maxIter =  itemToCompare
                iterMax = iter
            if(itemToCompare <= minIter):
                minIter = itemToCompare
        halfMaximum = (maxIter-minIter)/2
        halfMaximumValue = minIter + halfMaximum
        print(minIter, maxIter, halfMaximum, halfMaximumValue)
        for iter in range(len(fLambdaValues)):
            FWHMIter = math.fabs(math.fabs(fDataValues[iter])- halfMaximumValue)
            # print("############################# ", math.fabs(fDataValues[iter]) , FWHMIter)
            if (FWHMIter <= 0.1 and startFound == False):
                startLambdaValue = fLambdaValues[iter]
                startFound = True
            elif(FWHMIter <= 0.1):
                endLambdaValue = fLambdaValues[iter]
        if startLambdaValue == -1.0 and iterMax - 1 > 0 :
            iterMax = iterMax -1
            startLambdaValue = cpy.copy(fLambdaValues[iterMax])
        if endLambdaValue == 0.0 or (endLambdaValue - startLambdaValue) <0.0:
            endLambdaValue = fLambdaValues[-1]
        FWHM = endLambdaValue - startLambdaValue
        return FWHM
    
    def resampleFWHMDataInterval(self, fitsObject):
        fInitialDataValues = []
        fInitialLambdaValues = []
        for it in range(self.start,self.end):
            fInitialDataValues.append(cpy.copy(fitsObject.workingDataValues3[it]))
            fInitialLambdaValues.append(cpy.copy(fitsObject.workingLambdaValues[it]))
        interpolate = interpl.interp1d(fInitialLambdaValues, fInitialDataValues, "linear")
        step = (fInitialLambdaValues[1] - fInitialLambdaValues[0])/4
        finalStep = fInitialLambdaValues[len(fInitialLambdaValues)-1]
        lambdaValue = fInitialLambdaValues[0]
        fResampledLambdas = []
        fResampledValues  = []
        while(lambdaValue < finalStep):
            fResampledValues.append(interpolate(lambdaValue))
            fResampledLambdas.append(lambdaValue)
            lambdaValue += step
        return fResampledLambdas, fResampledValues
        
    def calculateErrorCayrelFormulae(self,fwhm,pixelSize,SNR):
        return 1.06 * math.sqrt(fwhm * math.fabs(pixelSize)) / SNR
    
    def convertVacuumToAir(self, lambda_vac):
        s = 1e4/lambda_vac
        n = 1 + 0.0000834254 + 0.02406147/(130 - s*s) + 0.00015998/(38.9 - s*s)
        lambda_air = lambda_vac/n
        return lambda_air
    
    def convertAirToVacuum(self, lambda_air):
        s = 1e4/lambda_air
        n = 1 + 0.00008336624212083 + 0.02408926869968/(130.1065924522 - s*s) + 0.0001599740894897/(38.92568793293 - s*s)
        lambda_vac = lambda_air*n
        return lambda_vac        

class blackBody(object):

    #Calculates de Integrated Flux by means of the Stefan-Boltzmann Law
    def StefanBoltzmann(self, teff):
        sigma = 5.67040 *1e-5 ##erg/(s cm2 K4) #Stefan-Boltzmann Constant
        sf = sigma * teff**4
        return sf 
    #######################################################################

    #Calculate Spectral Flux Density according to Planck's Law
    def calculateFluxSpectralDensity(self, inputLambda, teff):
        C1 = 3.742e+27 
        #C1 = 2*pi*h*c2 = 3.742e-5 (erg/s) cm2.
        # This need to be multiplied by an additional factor of 1e2 to take into account the unit conversion of inputLambda**5 
        C2 = 1.439e+08
        f = C1/ (inputLambda**5 * (math.exp(C2/inputLambda/teff)-1))
        return f
    #######################################################################

    
# lambda_air = 10830.2501
# print(lambda_vac, convertAirToVacuum(lambda_air))

