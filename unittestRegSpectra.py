import math as math
import numpy as np
import scipy.interpolate as interpl
import matplotlib.pyplot as plt
import copy as cpy

from support import*
import support as supp
#from starot  import*
#from fshift import*
#from starmod_tools import * 
##from rotBroad import *
#

from astropy.io import fits

                
class FITSObject(object):
    def __init__(self, inputFITS):
        self.specsLambdaValues = []
        self.workingLambdaValues = []
        self.dataValues = []
        self.workingDataValues1   = []
        self.workingDataValues2 = []
        self.workingDataValues3 = []
        self.regularLambdaValues = [] 
        self.regularDataValues = []     
        self.numberOfOrders = 0
        self.specs = []
        self.initialLambda = 0.0
        self.finalLambda   = 0.0
        self.deltaStep = 0.0
        self.carmenesType = 0
        self.bCarmenes = False
        self.object = ""
        self.date = 0.0
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
    def selectPixelRange(self, inputParams, bPlot):
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
        print (pixelRange)
        ipSpecFAp = inputParams.spectraFormatAperture
        if not pixelRange[2]:                        #Specified by Wavelenght in Angstroms
            m = 0
            if self.numberOfOrders > 1:
                for i in range(self.numberOfOrders):
                    # print ("order: ",i, "//")
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
                    plt.plot(self.workingLambdaValues,self.workingDataValues1 )
            else:
                print (len(self.specsLambdaValues))
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
                    plt.plot(self.workingLambdaValues,self.workingDataValues1 )
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
        
        return xlo, xhi, skp
    
    def plot(self, inputParams):
        print ("/plotInplot/")        
        print (self.numberOfOrders)
        if self.numberOfOrders > 1:
            for it in range(self.numberOfOrders):
                plt.plot(self.specsLambdaValues[it],self.dataValues[it])
        else:
            plt.plot(self.specsLambdaValues,self.dataValues)
    def checkDeltaStep(self):
        delta  = 1.0
        deltaAvg = 0.0
        bRegularDeltaStep = True
        for iter in range(1,len(self.workingLambdaValues)):
            deltaTemp  = self.workingLambdaValues[iter] - self.workingLambdaValues[iter-1]
            if delta - deltaTemp  > 0.0000001:
                #print (str(iter)+"/ "+str(deltaTemp))
                if iter > 1:
                    bRegularDeltaStep = False
                delta = deltaTemp
            deltaAvg += deltaTemp
        print ("len = " +str(len(self.workingLambdaValues)))
        deltaAvg /= len(self.workingLambdaValues)
        print ("deltaAvg = "+ str(deltaAvg))
        
        return not bRegularDeltaStep, deltaAvg
    def regularizeDeltaStep(self,deltaAvg):
        initialLambda = self.workingLambdaValues[0]
        # regularLambdaValues = []
        for iter in range(len(self.workingLambdaValues)):
            self.regularLambdaValues.append(iter*deltaAvg+initialLambda)
        # self.workingLambdaValues = cpy.copy(regularLambdaValues)
        return
        
    def resampleDataValues(self, workingLambdaValues):
        ###### Needed to build an interpolator function
        # The interpolator is build from the original Lambda Values
        interpolate = interpl.interp1d(self.workingLambdaValues,self.workingDataValues1, "linear")
        ########################################
        # If not specified, use the own Lambda values array, read from FITS
        # When specified it is used the lambda values array of another spectra
        if workingLambdaValues == None:
            workingLambdaValues = cpy.copy(self.regularLambdaValues)
        for iter in range(len(workingLambdaValues)):
            try :
                self.regularDataValues.append(interpolate(workingLambdaValues[iter]))
            except ValueError:
                self.regularDataValues.append(self.workingDataValues1[iter])
            #print( str(iter) + " " + str(workingLambdaValues[iter]) +  " " + str(regularDataValues[iter]))
        print (str(len(workingLambdaValues)), str(len(self.regularDataValues)))
        print ("")
        print ("interpolation ended")
        self.workingDataValues1  = cpy.copy(self.regularDataValues)
        self.workingLambdaValues = cpy.copy(self.regularLambdaValues)
        return
        
    def nanSanityCheck(self, nValArray = 1):
        num = len(self.workingLambdaValues)
        
        #### Selecting the array from the input parameter
        if nValArray == 1:
            wDV = self.workingDataValues1
        elif nValArray == 2:
            wDV = self.workingDataValues2
        else:
            wDV = self.workingDataValues3
        
        #### If it is found a NaN (a defect in a pixel of the sensor)
        #### it is repeated the previous value of the spectrum
        for iter in range(num):
            if math.isnan(wDV[iter]) and iter > 0 :
                wDV[iter] = wDV[iter-1]
        #### In order to detect inputs from cosmic rays (nan from star spectrum!)
        #### In the middle we can average the values
        for iter in range(num):
            if (iter > 10 and iter < (num -10)) and (math.fabs(wDV[iter]/wDV[iter-1])) > 5.0 :
                wDV[iter] = (wDV[iter-3] + wDV[iter+3])/2.0
        #### In order to delete algorithm artifacts (nan from star spectrum!)
        #### At the borders the values cannot be averaged
        for iter in range(num):
            if iter < 10 and (math.fabs(wDV[iter]/wDV[iter+15]) > 5.0):
                wDV[iter] = wDV[iter+15]
            if iter > num-11 and (math.fabs(wDV[iter]/wDV[iter-15]) > 5.0):
                wDV[iter] = wDV[iter-15]
                
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
    if ".fits" in inputFITS:
        #Retrieving basic data of the fits format file, the headers and the proper data
        objectFITS.fits_hdu_in.info()
        header = fits.getheader(inputFITS)
        # headerKeys = header.keys()
        count = 0
        for iterh in header.keys():
            count += 1
        numberOfCards = count #len(headerKeys)
        numberOfAxis = int(objectFITS.fits_hdu_in[0].header['NAXIS'])
        try:
            numberOfPointPerSpec = int(objectFITS.fits_hdu_in[0].header['NAXIS1'])
        except KeyError:
            print ("KEY NOT FOUND")
            numberOfPointPerSpec = int(-1)
        if numberOfAxis >= 2:
            numberOfOrders = int(objectFITS.fits_hdu_in[0].header['NAXIS2'])        
        else:
            numberOfOrders = 1
        ##FOR DEBUGGING PURPOSES
        ##----------------------------------------------------------------------------------
        #for iterh in range(numberOfCards):
        #    try:
        #        print "[ ", headerKeys[iterh], ",",header[iterh]," ]"
        #    except fits.VerifyError:
        #        print "error in format"
        ##----------------------------------------------------------------------------------
        
        try:
            ctype1 = objectFITS.fits_hdu_in[0].header['CTYPE1']
        except KeyError:
            print ("KEY NOT FOUND - CTYPE1")
            ctype1 = -1
        inlist = objectFITS.fits_hdu_in[0].data
        print ("inlist: ", inlist)
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
        elif numberOfAxis == 1 and ctype1 == "LINEAR":
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
            objectFITS.bCarmenes = True
            carmenesType = len(objectFITS.fits_hdu_in)
            print ("carmenesType = ", carmenesType)
            
            if carmenesType == 5 or carmenesType == 6: ########################ORIGINAL#######################
                objectFITS.specsLambdaValues = objectFITS.fits_hdu_in[4].data
                objectFITS.dataValues  = objectFITS.fits_hdu_in[1].data
            
            elif carmenesType == 4: #####################STACKED#########################
                objectFITS.specsLambdaValues = objectFITS.fits_hdu_in[3].data
                objectFITS.dataValues  = objectFITS.fits_hdu_in[1].data
                
                for i in range(len(objectFITS.specsLambdaValues)):
                    for j in range(len(objectFITS.specsLambdaValues[i])):
                        objectFITS.specsLambdaValues[i][j] = np.exp(objectFITS.specsLambdaValues[i][j])

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
        self.spectraFormatAperture = -1.0
        self.spectraFormatBand = -1.0
        self.cnt = 1
        self.initialLambda = 0.0
        self.finalLambda   = 0.0
        self.deltaStep     = 0.0
        self.numberOfWorkingDataValues = 0
    def readConfigFile(self):

        #self.workpath = workpath
        pixelExcl = []
        for i in range(5):
            pixelExcl.append([0,0])
        #print pixelExcl
        
        i = 0
        with open(self.spcfFile, 'r') as inFile:
            for line in iter(inFile.readline,''):
                if 'IM_PATH' in line:
                    start = line.find('=') + 1
                    end   = line.find('#') - 1
                    impathStr = line[start:end].split()
                    self.workpath = impathStr[0]
                    print (impathStr )
                if 'OBJ_NAME' in line:
                    start = line.find('=') + 1
                    end   = line.find('#') - 1
                    objInStr = line[start:end].split()
                    objectInputFits0 = objInStr[0]
                    objectInputFits = self.workpath + objectInputFits0
                    self.objectInputFits0 = objectInputFits0
                    self.objectInputFits = objectInputFits 
                    #print objectInputFits
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
                    print ("PIX_ZONE = ", pixelRange)
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
                    print (i, pixelExcl[i])
                    if i<4: i+=1
                if 'APERTURE' in line:
                    start = line.find('=') + 1
                    end   = line.find(' #')  - 1   
                    spectraFormatAperture = int(line[start:end])
                    self.spectraFormatAperture = spectraFormatAperture
                    #print spectraFormatAperture
        self.pixelExcl = pixelExcl
        print (pixelExcl)
     
def readConfigFile(spcfFile):
    configParams = spcfParams(spcfFile)
    configParams.readConfigFile()
    return configParams


def plot(inObj, inputParams):
    if inObj.numberOfOrders > 1:
        for it in range(inObj.numberOfOrders):
            plt.plot(inObj.specsLambdaValues[inputParams.spectraFormatAperture],inObj.dataValues[inputParams.spectraFormatAperture])
            plt.plot(inObj.workingLambdaValues,inObj.workingDataValues1)
    else:
        plt.plot(inObj.specsLambdaValues,inObj.dataValues)
    #if inputObject.numberOfOrders > 1:
    #    for it in range(inputObject.numberOfOrders):
    #        plt.plot(inputObject.specsLambdaValues[inputParams.spectraFormatAperture],inputObject.dataValues[inputParams.spectraFormatAperture])
    #        plt.plot(inputObject.workingLambdaValues,inputObject.workingDataValues1)
    #else:
    #    plt.plot(inputObject.specsLambdaValues,inputObject.dataValues)
    ##print "npts in inputObject.workingDataValues1", len(inputObject.workingDataValues1)




def unittestRegSpectra(spcfFile):

    ##################################################################################3
    #READING THE SPECIFICATION FILE
    ##################################################################################3
    inputParams = readConfigFile(spcfFile)
    
    ##################################################################################3
    #READING and plotting THE OBJECT FITS FILE
    ##################################################################################3
        
    inputObject = readFitsFile(inputParams.objectInputFits)
    #evaluate if 'pixel range' is given in pixel or in Angstroms
    #it is assumed that the start and end pixels specified belong to the same order 
    xlo,xhi, skp = inputObject.selectPixelRange(inputParams, True)
    print ("skp = ", skp)
    print(len (inputObject.workingLambdaValues))
    print (xlo,xhi)
    if inputObject.bCarmenes == True:
        bDeltaStepChecked, deltaStepAvg = inputObject.checkDeltaStep()
        if bDeltaStepChecked:
            inputObject.regularizeDeltaStep(deltaStepAvg)
            inputObject.resampleDataValues(None)
    #print ("skp = " + str(skp))
    print(len (inputObject.workingLambdaValues))
    print("####################################################################")
    print("input object initialLambda =  " + str(inputObject.initialLambda))
    print("input object final Lambda =  " + str(inputObject.workingLambdaValues[len(inputObject.workingLambdaValues)-1]))
    print("Nr of points = " + str(len(inputObject.workingLambdaValues)))
    print("####################################################################")
    # inputObject.workingDataValues2 = cpy.copy(inputObject.workingDataValues1)    
    # inputObject.workingDataValues3 = cpy.copy(inputObject.workingDataValues1)
    nptsInputObj = len(inputObject.workingDataValues1)
    print ("###########33333########   ",nptsInputObj)
    inputParams.finalLambda = inputObject.finalLambda ; inputParams.deltaStep   = inputObject.deltaStep
    inputParams.numberOfWorkingDataValues = nptsInputObj
    #the way the spectrum is plotted depends on the format
    #inputObject.plot(inputParams)
    xObjmean = supp.gtmean(nptsInputObj,inputObject.workingDataValues1,0,nptsInputObj-1,skp)
    print ("#########################################################")
    print ("xObjmean =",xObjmean)
    print ("#########################################################")
    #for iter in range (len(inputObject.workingDataValues1)):
    #    print "[",inputObject.workingLambdaValues[iter], ",", inputObject.workingDataValues1[iter],"]"
    print ("#########################################################")
    ###The initial lambda and the delta step are different in reading the FITS file and in the selection of the pixel range
    #midLambda = inputObject.initialLambda + inputParams.deltaStep * len(inputObject.workingLambdaValues) / 2.0
    ##print "midLambda = ", midLambda
    #cnst = (midLambda/299792.458)/inputParams.deltaStep     #Conversion Factor from Doppler Effect Formula
    #print "const = ", cnst
    plt.plot(inputObject.workingLambdaValues,inputObject.workingDataValues1 )
    print ("///////////////////////////////////////////////////////////////")
    
    plt.title(inputParams.objectInputFits0)
    plt.xlabel("Lambda")
    plt.ylabel("Normalised Flux")
    plt.legend()
    plt.show()
    print ("///////////////////////////////////////////////////////////////")
    
    #print len(inlist), len(inlist[0])
    #print inlist[0]
    if inputObject.fits_hdu_in != None:
        inputObject.fits_hdu_in.close()

unittestRegSpectra("testnormalisation.sm")
