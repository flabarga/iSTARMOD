import numpy as np
import matplotlib.pyplot as plt
import copy as cpy
import math as math
import support as supp
from pathlib import Path
#from starot  import*
#from fshift import*
#from starmod_tools import * 
##from rotBroad import *
from iStarmod_tools import spcfParams as sm_parameters

from astropy.io import fits

                
class FITSObject(object):
    def __init__(self, inputFITS):
        self.specsLambdaValues = []
        self.workingLambdaValues = []
        self.dataValues = []
        self.workingDataValues1   = []
        self.workingDataValues2 = []
        self.workingDataValues3 = []        
        self.numberOfOrders = 0
        self.specs = []
        self.initialLambda = 0.0
        self.finalLambda   = 0.0
        self.deltaStep = 0.0
        try:
            self.fits_hdu_in = fits.open(inputFITS, lazy_load_hdus=False)
        except IOError:
            print ('File does not exist: {0!r}'.format(inputFITS))
            self.fits_hdu_in = None
            return
    def setPixelExclusionRanges(self, inputParams, wvl, skp, i, j, xlo, m):
        #The following seemed redundacy of code is intended to allow a 'disordered' specification 
        #of the exclusion zones
        if (wvl == True and i != -1):
            if((self.specsLambdaValues[i][j-1] <= inputParams.pixelExcl[0][0] and self.specsLambdaValues[i][j] > inputParams.pixelExcl[0][0]) or 
                (self.specsLambdaValues[i][j-1] <= inputParams.pixelExcl[0][1] and self.specsLambdaValues[i][j] > inputParams.pixelExcl[0][1])):
                skp[m] = j-xlo
                m += 1
            if((self.specsLambdaValues[i][j-1] <= inputParams.pixelExcl[1][0] and self.specsLambdaValues[i][j] > inputParams.pixelExcl[1][0]) or 
                (self.specsLambdaValues[i][j-1] <= inputParams.pixelExcl[1][1] and self.specsLambdaValues[i][j] > inputParams.pixelExcl[1][1])):
                skp[m] = j-xlo
                m += 1
            if((self.specsLambdaValues[i][j-1] <= inputParams.pixelExcl[2][0] and self.specsLambdaValues[i][j] > inputParams.pixelExcl[2][0]) or 
                (self.specsLambdaValues[i][j-1] <= inputParams.pixelExcl[2][1] and self.specsLambdaValues[i][j] > inputParams.pixelExcl[2][1])):
                skp[m] = j-xlo
                m += 1
            if((self.specsLambdaValues[i][j-1] <= inputParams.pixelExcl[3][0] and self.specsLambdaValues[i][j] > inputParams.pixelExcl[3][0]) or 
                (self.specsLambdaValues[i][j-1] <= inputParams.pixelExcl[3][1] and self.specsLambdaValues[i][j] > inputParams.pixelExcl[3][1])):
                skp[m] = j-xlo
                m += 1
            if((self.specsLambdaValues[i][j-1] <= inputParams.pixelExcl[4][0] and self.specsLambdaValues[i][j] > inputParams.pixelExcl[4][0]) or 
                (self.specsLambdaValues[i][j-1] <= inputParams.pixelExcl[4][1] and self.specsLambdaValues[i][j] > inputParams.pixelExcl[4][1])):
                skp[m] = j-xlo
                m += 1
        elif (wvl == True and i == -1):
            if((self.specsLambdaValues[j-1] <= inputParams.pixelExcl[0][0] and self.specsLambdaValues[j] > inputParams.pixelExcl[0][0]) or 
                (self.specsLambdaValues[j-1] <= inputParams.pixelExcl[0][1] and self.specsLambdaValues[j] > inputParams.pixelExcl[0][1])):
                skp[m] = j-xlo
                m += 1
            if((self.specsLambdaValues[j-1] <= inputParams.pixelExcl[1][0] and self.specsLambdaValues[j] > inputParams.pixelExcl[1][0]) or 
                (self.specsLambdaValues[j-1] <= inputParams.pixelExcl[1][1] and self.specsLambdaValues[j] > inputParams.pixelExcl[1][1])):
                skp[m] = j-xlo
                m += 1
            if((self.specsLambdaValues[j-1] <= inputParams.pixelExcl[2][0] and self.specsLambdaValues[j] > inputParams.pixelExcl[2][0]) or 
                (self.specsLambdaValues[j-1] <= inputParams.pixelExcl[2][1] and self.specsLambdaValues[j] > inputParams.pixelExcl[2][1])):
                skp[m] = j-xlo
                m += 1
            if((self.specsLambdaValues[j-1] <= inputParams.pixelExcl[3][0] and self.specsLambdaValues[j] > inputParams.pixelExcl[3][0]) or 
                (self.specsLambdaValues[j-1] <= inputParams.pixelExcl[3][1] and self.specsLambdaValues[j] > inputParams.pixelExcl[3][1])):
                skp[m] = j-xlo
                m += 1
            if((self.specsLambdaValues[j-1] <= inputParams.pixelExcl[4][0] and self.specsLambdaValues[j] > inputParams.pixelExcl[4][0]) or 
                (self.specsLambdaValues[j-1] <= inputParams.pixelExcl[4][1] and self.specsLambdaValues[j] > inputParams.pixelExcl[4][1])):
                skp[m] = j-xlo
                m += 1
        elif (wvl == False):
            if(j == inputParams.pixelExcl[0][0] or j == inputParams.pixelExcl[0][1]):
                skp[m] = j-xlo
                m += 1
            if(j == inputParams.pixelExcl[1][0] or j == inputParams.pixelExcl[1][1]):
                skp[m] = j-xlo
                m += 1
            if(j == inputParams.pixelExcl[2][0] or j == inputParams.pixelExcl[2][1]):
                skp[m] = j-xlo
                m += 1
            if(j == inputParams.pixelExcl[3][0] or j == inputParams.pixelExcl[3][1]):
                skp[m] = j-xlo
                m += 1
            if(j == inputParams.pixelExcl[4][0] or j == inputParams.pixelExcl[4][1]):
                skp[m] = j-xlo
                m += 1                  
        return m
    def selectPixelRange(self, inputParams, fig, subplot1, bPlot = True):
        self.workingLambdaValues = []
        self.workingDataValues1  = []
        print ("selectPixelRange///////////////////////////////////////////////////////////////")
        into = False
        xhi = 0
        xlo = 0
        pixelRange = []
        pixelRange.append(inputParams.pixelRange[0])
        pixelRange.append(inputParams.pixelRange[1])
        pixelRange.append(inputParams.pixelRange[2])
        skp = []
        for i in range(10):
            skp.append(0) 
        print (pixelRange[2])
        ipSpecFAp = inputParams.spectraFormatAperture
        if not pixelRange[2]:                        #Specified by Wavelenght in Angstroms
            m = 0
            if self.numberOfOrders > 1:
                for i in range(self.numberOfOrders):
                    print ("order: ",i, "//")
                    if i != ipSpecFAp:
                        continue
                    for j in range(len(self.specsLambdaValues[i])):
                        if (into == False): 
                            xlo = j
                        if ((self.specsLambdaValues[i][j] > pixelRange[0]) and (self.specsLambdaValues[i][j] < pixelRange[1])):
                                into = True
                                #print j,
                                self.workingLambdaValues.append(self.specsLambdaValues[i][j])
                                #if j<len(self.dataValues[i]):
                                self.workingDataValues1 .append(self.dataValues[i][j])
                                #else:
                                #    continue
                                ##FOR DEBUGGING PURPOSES
                                #print self.specsLambdaValues[i][j], self.dataValues[i][j]
                                ##------------------------------------------------------------
                        
                        if into == True :
                            m = self.setPixelExclusionRanges(inputParams, True, skp, i, j, xlo, m)
                            
                        if(self.specsLambdaValues[i][j] >= pixelRange[1]  and into == True):
                                self.initialLambda = self.workingLambdaValues[0]
                                self.finalLambda    = self.workingLambdaValues[len(self.workingLambdaValues)-1]
                                self.deltaStep = self.workingLambdaValues[1]-self.workingLambdaValues[0]
                                if len(self.workingDataValues1) < inputParams.numberOfWorkingDataValues:
                                    self.workingLambdaValues.append(self.specsLambdaValues[i][j])
                                    self.workingDataValues1 .append(self.dataValues[i][j])
                                xhi = j
                                break
                    #print ("//")
                    self.nanSanityCheck()
                if bPlot:
                    print ("/plot1/")
                    subplot1.plot(self.workingLambdaValues,self.workingDataValues1,'b-', lw = 1.5 )
            else:
                print (str(len(self.specsLambdaValues)))
                for j in range(len(self.specsLambdaValues)):
                            if (into == False): xlo = j
                            if ((self.specsLambdaValues[j] > pixelRange[0]) and (self.specsLambdaValues[j] < pixelRange[1])):
                                    #print j
                                    into = True
                                    self.workingLambdaValues.append(self.specsLambdaValues[j])
                                    self.workingDataValues1 .append(self.dataValues[j])
                                    ##FOR DEBUGGING PURPOSES
                                    #print self.specsLambdaValues[j], self.dataValues[j]
                                    ##----------------------------------------------------
                            if into == True :
                                m = self.setPixelExclusionRanges(inputParams, True, skp, -1, j, xlo, m) 
                            if(self.specsLambdaValues[j] >= pixelRange[1]):
                                xhi = j
                                break
                self.nanSanityCheck()
                if bPlot:
                    print ("/plot2/")
                    subplot1.plot(self.workingLambdaValues,self.workingDataValues1,'b-', lw = 1.5 )
                #xlo = xlo+1
                #xhi = xhi-1
        else:                                        #Specified by Pixels
            xlo = pixelRange[0]
            xhi = pixelRange[1]
            m = 0
            if self.numberOfOrders > 1:
                for j in range(xlo, xhi):
                        if (j < len(self.specsLambdaValues[ipSpecFAp])):
                            self.workingLambdaValues.append(self.specsLambdaValues[ipSpecFAp][j])
                            self.workingDataValues1 .append(self.dataValues[ipSpecFAp][j])
                            #print self.workingLambdaValues[j-xlo], '\t', self.workingDataValues1 [j-xlo]
                        m = self.setPixelExclusionRanges(inputParams, False, skp, 0, j, xlo, m)
                self.initialLambda = self.workingLambdaValues[0]
                self.deltaStep = self.workingLambdaValues[1]-self.workingLambdaValues[0]
            else:
                for j in range(xlo, xhi):
                        if (j < len(self.specsLambdaValues)):
                            #print j
                            self.workingLambdaValues.append(self.specsLambdaValues[j])
                            self.workingDataValues1 .append(self.dataValues[j])
                        m = self.setPixelExclusionRanges(inputParams, False, skp, 0, j, xlo, m)
                self.initialLambda = self.workingLambdaValues[0]
                self.deltaStep = self.workingLambdaValues[1]-self.workingLambdaValues[0]
            self.nanSanityCheck()
            if bPlot:
                print ("/plot3/")                
                plt.plot(self.workingLambdaValues,self.workingDataValues1 )
        
        return xlo, xhi, skp, fig, subplot1 
    
    def nanSanityCheck(self):
        num = len(self.workingLambdaValues)
        print (num)
        for iter in range(num):
            if math.isnan(self.workingDataValues1[iter]) and iter > 0 :
                self.workingDataValues1[iter] = self.workingDataValues1[iter-1]
        for iter in range(num):
            if (iter > 1 and iter < (num -1)) and (self.workingDataValues1[iter]/self.workingDataValues1[iter-2]) > 5.0 :
                # print "here???"
                self.workingDataValues1[iter] = (self.workingDataValues1[iter-2] + self.workingDataValues1[iter+2])/2.0
                
def readLambdaDataMultiSpec(header, numberOfCards,numberOfPointPerSpec,numberOfOrders, inlist):
    
    stringHeaders = ""
    strDebug = ""
    for iterh in range(numberOfCards):
        try:
            strDebug = str(header[iterh])
            stringHeaders += str(header[iterh])
            if (len(strDebug) == 67):
                stringHeaders += " "
        except fits.VerifyError:
            print ("error in format")

    retValueFind = stringHeaders.find("spec1")

    retValueFind1 = 0
    firstStep = 0
    secondStep = 0
    iterator = 1
    specs = []
    #find the characters between one 'spec' and another
    while retValueFind != -1 and retValueFind1 != -1 and iterator <= len(stringHeaders) :
        retValueFindDoubleQuotes = stringHeaders[secondStep+retValueFind:].find('\"')#secondStep in the first iteration is 0
        firstStep = secondStep+retValueFind+retValueFindDoubleQuotes #secondStep in the first iteration is 0
        searchString = "spec" + str(iterator+1)
        retValueFind1 = stringHeaders[firstStep:].find(searchString)
        if retValueFind1 != -1:
            secondStep = firstStep+retValueFind1
        else:
            secondStep = firstStep+68 #if it is not found the kwd 'spec' the use the whole line (we are at the end of the 'specs')
        stringFound = stringHeaders[firstStep+1:secondStep]
        valuesOfSpec = stringFound.split()
        ##FOR DEBUGGING PURPOSES##############################3
        #print iterator, valuesOfSpec
        ##------------------------------------
        if len(valuesOfSpec)>5:
            specs.append([float(valuesOfSpec[3]), float(valuesOfSpec[4])])#, int(valuesOfSpec[5])])
        retValueFind = stringHeaders[secondStep:].find("spec")
        iterator += 1
    
    specsLambdaValues = []
    N_orders = numberOfOrders
    
    for it1 in range(N_orders):
        specsLambdaValues.append([])
        for it2 in range(numberOfPointPerSpec):
            specsLambdaValues[it1].append(specs[it1][0]+it2*specs[it1][1])
    
    # remove borders
    inlist2 = inlist.tolist()
    for i in range(N_orders):
        del specsLambdaValues[i][0:40]
        del specsLambdaValues[i][-40:]
        del inlist2[i][0:40]
        del inlist2[i][-40:]
        
        
    return specsLambdaValues, inlist2, specs

def readFitsFile(inputFITS):
    
    objectFITS = FITSObject(inputFITS)
    if ".fits" in inputFITS.name:
        #Retrieving basic data of the fits format file, the headers and the proper data
        objectFITS.fits_hdu_in.info()
        header = fits.getheader(inputFITS)
        # headerKeys = header.keys()
        count = 0
        for iterh in header.keys():
            count += 1
        numberOfCards = count
        numberOfAxis = int(objectFITS.fits_hdu_in[0].header['NAXIS'])
        try:
            numberOfPointPerSpec = int(objectFITS.fits_hdu_in[0].header['NAXIS1'])
            print("NAXIS1: ", numberOfPointPerSpec)
        except KeyError:
            print ("KEY NOT FOUND")
            numberOfPointPerSpec = int(-1)
        if numberOfAxis >= 2:
            numberOfOrders = int(objectFITS.fits_hdu_in[0].header['NAXIS2'])        
        else:
            numberOfOrders = 1
        
        ##FOR DEBUGGING PURPOSES
        ##----------------------------------------------------------------------------------
        for iterh in header.keys():
            try:
                print("[ ", iterh, ",",header[iterh]," ]")
            except fits.VerifyError:
                print ("error in format")
        ##----------------------------------------------------------------------------------
        print ("No errors found")
        try:
            ctype1 = objectFITS.fits_hdu_in[0].header['CTYPE1']
        except KeyError:
            print ("KEY NOT FOUND - CTYPE1")
            ctype1 = -1
        inlist = objectFITS.fits_hdu_in[0].data
        print ("inlist: ", str(inlist))
        #Building the Array/s of spectrum's ordinate of wavelenghts
        if numberOfAxis == 2 and ctype1 == "MULTISPE":
            #In a multispec format, the wavelenghts are specified by means of
            #a series of paramters provided in the headers 
            specsLambdaValues, inlist2, objectFITS.specs = readLambdaDataMultiSpec(header,numberOfCards,numberOfPointPerSpec,numberOfOrders, inlist)
            #the way the spectrum is plotted depends on the format
            #for it in range(numberOfOrders):
            #    plt.plot(specsLambdaValues[it],inlist2[it])
            objectFITS.specsLambdaValues = specsLambdaValues
            objectFITS.dataValues  = inlist2 
            objectFITS.numberOfOrders = numberOfOrders
            return objectFITS
        elif numberOfAxis == 2 and (ctype1 == "LINEAR" or ctype1 == "Linear" or ctype1 == "WAVELENGTH"):
            #in a 'onedspec'-linear format, the wavelenghts are specified by means of 
            #only one set of parameters of the linear formulae
            objectFITS.initialLambda = float(objectFITS.fits_hdu_in[0].header['CRVAL1'])
            objectFITS.deltaStep = float(objectFITS.fits_hdu_in[0].header['CDELT1'])
            objectFITS.specsLambdaValues = objectFITS.fits_hdu_in[0].data[0]
            objectFITS.dataValues  = objectFITS.fits_hdu_in[0].data[1]
            objectFITS.numberOfOrders = 0
            
            return objectFITS
        elif numberOfAxis == 1 and (ctype1 == "LINEAR" or ctype1 == "Linear" or ctype1 == "WAVELENGTH"):
            #in a 'onedspec'-linear format, the wavelenghts are specified by means of 
            #only one set of parameters of the linear formulae
            objectFITS.initialLambda = float(objectFITS.fits_hdu_in[0].header['CRVAL1'])
            objectFITS.deltaStep = float(objectFITS.fits_hdu_in[0].header['CDELT1'])
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
            carmenesType = len(objectFITS.fits_hdu_in)
            print( "carmenesType = ", str(carmenesType))
            
            if carmenesType == 5: ########################ORIGINAL#######################
                objectFITS.specsLambdaValues = objectFITS.fits_hdu_in[4].data
                objectFITS.dataValues  = objectFITS.fits_hdu_in[1].data
            
            elif carmenesType == 4: #####################STACKED#########################
                objectFITS.specsLambdaValues = objectFITS.fits_hdu_in[3].data
                objectFITS.dataValues  = objectFITS.fits_hdu_in[1].data
                
                for i in range(len(objectFITS.specsLambdaValues)):
                    for j in range(len(objectFITS.specsLambdaValues[i])):
                        objectFITS.specsLambdaValues[i][j] = np.exp(objectFITS.specsLambdaValues[i][j])
                        
            elif carmenesType == 6: #####################TELLURIC CORRECTED ##############
                objectFITS.specsLambdaValues = objectFITS.fits_hdu_in[4].data
                objectFITS.dataValues  = objectFITS.fits_hdu_in[1].data
                print ("dataValues: ",len(objectFITS.dataValues[25]))
                
            objectFITS.numberOfOrders = objectFITS.specsLambdaValues.shape[0]
            for i in range(len(objectFITS.specsLambdaValues)):
                    print (i, "= [",objectFITS.specsLambdaValues[i][0],",",objectFITS.specsLambdaValues[i][len(objectFITS.specsLambdaValues[i])-1],"]")
            ################Normalization###############
            norm = 0.0
            count = 0
            for i in range(objectFITS.numberOfOrders):
                for j in range(len(objectFITS.dataValues[i])):
                    if not math.isnan(objectFITS.dataValues[i][j]) and objectFITS.dataValues[i][j] != 0.0:
                        norm += objectFITS.dataValues[i][j]
                        count += 1
            norm = norm/count
            if norm != 0.0:
                for i in range(objectFITS.numberOfOrders):
                    for j in range(len(objectFITS.dataValues[i])):
                        if not math.isnan(objectFITS.dataValues[i][j]):                        
                            objectFITS.dataValues[i][j] = objectFITS.dataValues[i][j]/norm
            ############################################
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
        print ("= [",objectFITS.specsLambdaValues[0],",",objectFITS.specsLambdaValues[len(objectFITS.specsLambdaValues)-1],"]")
        print ("Number of points of LambdaValues: ", len(objectFITS.specsLambdaValues))
        print ("Number of points of dataValues: ", len(objectFITS.dataValues))
        return objectFITS
    print("Bad Assignment of Headers. The fits object cannot be read")
    return


         
def readConfigFile(spcfFile):
    configParams = sm_parameters(spcfFile)
    configParams.readConfigFile(spcfFile)
    return configParams

def unittestreadFITS(spcfFile, fig, subplot1, debugging = False):

    ##################################################################################3
    #READING THE SPECIFICATION FILE
    ##################################################################################3
    inputParams = readConfigFile(spcfFile)
    
    p = Path(inputParams.workpath)
    print("Path = " + str(p))
    x = list(p.glob("*.fits"))
    print(x)
    
    inputObject = readFitsFile(inputParams.objectInputFits)
    # evaluate if 'pixel range' is given in pixel or in Angstroms
    # it is assumed that the start and end pixels specified belong to the same order 
    # the way the spectrum is plotted depends on the format of the fits or dat file
    xlo,xhi, skp, fig, subplot1 = inputObject.selectPixelRange(inputParams, fig, subplot1)
    if debugging: 
        print ("skp = ", str(skp))
        print (xlo, xhi)
    # inputObject.workingDataValues2 = cpy.copy(inputObject.workingDataValues1)    
    # inputObject.workingDataValues3 = cpy.copy(inputObject.workingDataValues1)
    # nptsInputObj = len(inputObject.workingDataValues1)
    # if debugging:   print ("###########33333########   ",nptsInputObj)
    # inputParams.finalLambda = inputObject.finalLambda ; inputParams.deltaStep   = inputObject.deltaStep
    # inputParams.numberOfWorkingDataValues = nptsInputObj
    # xObjmean = supp.gtmean(nptsInputObj,inputObject.workingDataValues1,0,nptsInputObj-1,skp)
    # print ("#########################################################")
    # print ("xObjmean =",xObjmean)
    # print ("#########################################################")
    # #for iter in range (len(inputObject.workingDataValues1)):
    # #    print "[",inputObject.workingLambdaValues[iter], ",", inputObject.workingDataValues1[iter],"]"
    # print ("#########################################################")
    # ### The initial lambda and the delta step are different in reading the FITS file and in the selection of the pixel range
    # #m  idLambda = inputObject.initialLambda + inputParams.deltaStep * len(inputObject.workingLambdaValues) / 2.0
    # ##  print "midLambda = ", midLambda
    # ##  cnst = (midLambda/299792.458)/inputParams.deltaStep     #Conversion Factor from Doppler Effect Formula
    # ##  print "const = ", cnst
    
    print ("///////////////////////////////////////////////////////////////")    
    subplot1.set_title(inputParams.objectInputFits0)
    plt.xlabel("Lambda")
    plt.ylabel("Normalised Flux")
    print ("///////////////////////////////////////////////////////////////")
        
    if inputObject.fits_hdu_in != None:
        inputObject.fits_hdu_in.close()
    return fig, subplot1
