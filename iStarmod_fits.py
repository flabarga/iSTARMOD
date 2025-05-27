# -*- coding: utf-8 -*-
#=======================================================================
## Name:     iSTARMOD-FITS
# Filename: iStarmod_fits.py
# Type:     Creating and Processing FITS Objects in iSTARMOD
################################################################
# Language: Python 3.5
# Purpose:  Class definition of FITSObject:
#           Building FITS Objects from the read files
#           Processing the 1D fits data (spectra)
################################################################

from astropy.io import fits
import math
import numpy as np
import scipy.interpolate as interpl
import matplotlib.pyplot as plt
import copy as cpy


                                
class FITSObject(object):
    def __init__(self, inputFITS):
        self.specsLambdaValues     = []
        self.workingLambdaValues   = []
        self.regularLambdaValues   = [] 
        self.dataValues            = []
        self.workingDataValues1    = []
        self.workingDataValues2    = []
        self.workingDataValues3    = []
        self.regularDataValues     = []
        self.errWorkingDataValues  = []
        self.errorValues           = [] 
        self.numberOfOrders        = 0
        self.specs                 = []
        self.initialLambda         = 0.0
        self.finalLambda           = 0.0
        self.deltaStep             = 0.0
        self.carmenesType          = 0
        self.bCarmenes             = False
        self.object                = ""
        self.date                  = 0.0
        self.RVfromFITS            = 0.0
        self.BaryCorr              = 0.0
        self.SNR                   = 1.0   
        try:
            self.fits_hdu_in = fits.open(inputFITS, lazy_load_hdus=False)
        except IOError:
            print ('File does not exist: {0!r}'.format(inputFITS))
            return
        self.aperture_primary = 0
                    
    def readLambdaDataMultiSpec(self, header, numberOfCards,numberOfPointPerSpec,numberOfOrders, inlist):
    
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
                secondStep = firstStep+68 #if the kwd 'spec' is not found then use the whole line (we are at the end of the 'specs')
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
            
        #for i in range(N_orders):
        #    print len(specsLambdaValues[i]), len (inlist2[i])
        ##trying to erase spurious signals or artifacts
        #for i in range(N_orders):
        #    for it2 in range(4):
        #        meanValueFlux4it = np.fabs(inlist2[i][it2+1])
        #        ##FOR DEBUGGING PURPOSES
        #        #if i == 108:
        #            #print meanValueFlux4it , inlist2[i][it2] 
        #        ##------------------------------------------------------
        #        if (np.fabs(inlist2[i][it2]) > 5*meanValueFlux4it):
        #            inlist2[i][it2] = meanValueFlux4it
        #        if (meanValueFlux4it > 5*np.fabs(inlist2[i][it2])):
        #            inlist2[i][it2+1] = inlist2[i][it2]  
        #    for it2 in range(len(inlist2[i])-4):
        #        meanValueFlux4it = (np.fabs(inlist2[i][it2-4])+np.fabs(inlist2[i][it2+4]))/2 #seeking for 1 pixel variance
        #        ##FOR DEBUGGING PURPOSES
        #        #if i == 108:
        #        #    print meanValueFlux4it , inlist2[i][it2]
        #        #if inlist2[i][it2] > 100:
        #        #    print i, it2, inlist2[i][it2]
        #        ##------------------------------------------------------
        #        if (np.fabs(inlist2[i][it2]) > 5*meanValueFlux4it):
        #            inlist2[i][it2] = meanValueFlux4it
        #
        # previously, in versions prior to v7.x the following lines were outside the function, and this function were outside the class.
        # Now hter is no need for returning any value
        self.specsLambdaValues = specsLambdaValues
        self.dataValues  = inlist2 
        self.numberOfOrders = numberOfOrders
                
        return # specsLambdaValues, inlist2, specs
    
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
        
    def selectPixelRange(self, inputParams, primary, subPlot, bPlot):
        
        self.workingLambdaValues = []
        self.workingDataValues1  = []
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
        # print (str(pixelRange[2]))
        ipSpecFAp = inputParams.spectraFormatAperture
        if primary == True and inputParams.diffFormatPrimary == True:
           ipSpecFAp =  inputParams.spectraFormatAperturePrimary
        print(ipSpecFAp)
        if not pixelRange[2]:                        #Specified by Wavelenght in Angstroms
            m = 0
            if self.numberOfOrders > 1:
                for i in range(self.numberOfOrders):
                    #print ("order: ",i, "//")
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
                                if self.carmenesType == 5 or self.carmenesType == 6:
                                    self.errWorkingDataValues.append(self.errorValues[i][j])
                        if into == True :
                            m = self.setPixelExclusionRanges(inputParams, True, skp, i, j, xlo, m)
                            
                        if(self.specsLambdaValues[i][j] >= pixelRange[1]  and into == True):
                                self.initialLambda = self.workingLambdaValues[0]
                                self.finalLambda    = self.workingLambdaValues[len(self.workingLambdaValues)-1]
                                self.deltaStep = self.workingLambdaValues[1]-self.workingLambdaValues[0]
                                if len(self.workingDataValues1) < inputParams.numberOfWorkingDataValues:
                                    self.workingLambdaValues.append(self.specsLambdaValues[i][j])
                                    self.workingDataValues1 .append(self.dataValues[i][j])
                                    if self.carmenesType == 5 or self.carmenesType == 6:
                                        self.errWorkingDataValues.append(self.errorValues[i][j])
                                xhi = j
                                break
                self.nanSanityCheck()
                if bPlot:
                    subPlot.plot(self.workingLambdaValues,self.workingDataValues1 )
            else:
                #print (len(self.specsLambdaValues))
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
                    subPlot.plot(self.workingLambdaValues,self.workingDataValues1 )
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
                subPlot.plot(self.workingLambdaValues,self.workingDataValues1 )
        if (self.deltaStep == 0.0):
            self.deltaStep = self.workingLambdaValues[1]-self.workingLambdaValues[0]
        return xlo, xhi, skp
    
    def plot(self, inputParams):
        # print (str(self.numberOfOrders))
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
        deltaAvg /= len(self.workingLambdaValues)
        # print ("deltaAvg = "+ str(deltaAvg))
        # print ("len = " +str(len(self.workingLambdaValues)))
        return not bRegularDeltaStep, deltaAvg
    def regularizeDeltaStep(self,deltaAvg):
        initialLambda = self.workingLambdaValues[0]
        # regularLambdaValues = []
        for iter in range(len(self.workingLambdaValues)):
            self.regularLambdaValues.append(iter*deltaAvg+initialLambda)
        # self.workingLambdaValues = cpy.copy(regularLambdaValues)
        return
    def regularizeDeltaStepAndTrim(self,deltaAvg, nrpoints):
        if(nrpoints <= len(self.workingLambdaValues)):
            initialIndex = int((len(self.workingLambdaValues) - nrpoints)/2)
            initialLambda = self.workingLambdaValues[initialIndex]
            # regularLambdaValues = []
            for iter in range(nrpoints):
                self.regularLambdaValues.append(iter*deltaAvg+initialLambda)
            # del(self.workingLambdaValues)
            # self.workingLambdaValues = cpy.copy(self.regularLambdaValues)
        return
    def regularizeDeltaStepWvl(self,deltaStep):
        initialLambda = self.workingLambdaValues[0]
        finalLambda   = self.workingLambdaValues[len(self.workingLambdaValues)-1]
        _iter_ = 0
        valueLambda = 0.0
        while valueLambda <= finalLambda:
             self.regularLambdaValues.append( _iter_*deltaStep + initialLambda)
             _iter_ = _iter_ + 1
        return _iter_
        
    def resampleDataValues(self, workingLambdaValues):
        ###### Needed to build an interpolator function
        # The interpolator is build from the original Lambda Values
        interpolate = interpl.interp1d(self.workingLambdaValues,self.workingDataValues1, "linear")
        ########################################
        # If not specified, use the own Lambda values array, read from FITS
        # When specified it is used the lambda values array of another spectra
        if workingLambdaValues == None:
            workingLambdaValues = cpy.copy(self.regularLambdaValues)
        for _iter_ in range(len(workingLambdaValues)):
            try :
                self.regularDataValues.append(interpolate(workingLambdaValues[_iter_]))
            except ValueError:
                self.regularDataValues.append(self.workingDataValues1[len(self.workingDataValues1)-1])
            #print( str(iter) + " " + str(workingLambdaValues[iter]) +  " " + str(regularDataValues[iter]))
        # print (str(len(workingLambdaValues)), str(len(self.regularDataValues)))
        # print ("")
        # print ("interpolation ended")
        del self.workingDataValues1[:]
        self.workingDataValues1  = cpy.copy(self.regularDataValues)
        del self.workingLambdaValues[:]
        self.workingLambdaValues = cpy.copy(self.regularLambdaValues)
        return
     
    def nanSanityCheck(self, nValArray = 1):
        
        
        #### Selecting the array from the input parameter
        if nValArray == 1:
            wDV = self.workingDataValues1
        elif nValArray == 2:
            wDV = self.workingDataValues2
        else:
            wDV = self.workingDataValues3
        num = len(wDV)
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
            
    def writeToDat(self, specName, nValArray=3):

        numLambdas = len(self.workingLambdaValues)
        # end   = specName.find('fits')
        resOutfile = open(specName + ".dat" , 'w')
        
        if nValArray == 2:
            dataValues = self.workingDataValues2
        elif nValArray == 3:
            dataValues = self.workingDataValues3
        else:
            dataValues = self.workingDataValues1
        
        for it in range(numLambdas):
            resOutfile.write(str(self.workingLambdaValues[it]) + "    "+ str(dataValues[it]) + '\n')
        resOutfile.close()
        return
    
    def isInPixelExcl(self, fValue, inputParams):
        if ((fValue >= inputParams.pixelExcl[0][0] and fValue <= inputParams.pixelExcl[0][1]) or 
            (fValue >= inputParams.pixelExcl[1][0] and fValue <= inputParams.pixelExcl[1][1]) or
            (fValue >= inputParams.pixelExcl[2][0] and fValue <= inputParams.pixelExcl[2][1]) or 
            (fValue >= inputParams.pixelExcl[3][0] and fValue <= inputParams.pixelExcl[3][1]) or 
            (fValue >= inputParams.pixelExcl[4][0] and fValue <= inputParams.pixelExcl[4][1])):
            return True
        else:
            return False

    def normalizeSpectrum(self, inputParams, bWhole):
        norm = 0.0
        count = 1
        if(bWhole):
            for i in range(self.numberOfOrders):
                for j in range(len(self.dataValues[i])):
                    if (not math.isnan(self.dataValues[i][j]) and self.dataValues[i][j] != 0.0) and ():
                        norm += self.dataValues[i][j]
                        count += 1
            norm = norm/count
            if norm != 0.0:
                for i in range(self.numberOfOrders):
                    for j in range(len(self.dataValues[i])):
                        if not math.isnan(self.dataValues[i][j]):                        
                            self.dataValues[i][j] = self.dataValues[i][j]/norm
        else:
            order  = inputParams.spectraFormatAperture
            for j in range(len(self.dataValues[order])):
                if ((not math.isnan(self.dataValues[order][j]) and self.dataValues[order][j] != 0.0) and not self.isInPixelExcl(self.dataValues[order][j], inputParams)):
                        norm += self.dataValues[order][j]
                        count += 1
            norm = norm/count
            if norm != 0.0:
                for j in range(len(self.dataValues[order])):
                    if not math.isnan(self.dataValues[order][j]):                        
                        self.dataValues[order][j] = self.dataValues[order][j]/norm
        return 
        
    def evalReject (self, valueSpectrum, valueAdjustment, above,  reject_param):
        
        # In order to evaluate the rejection of points it is used a different 
        # expression if the point to be discarded is above the supposed continuum or below it
        if (not above):
            return (valueSpectrum > valueAdjustment * reject_param)
        else:
            return (valueSpectrum < valueAdjustment * (2.05 - reject_param))      
            
          
    def continuum_det_and_normalization (self, rejt, strain, inputParams):
        
        # This fragment of code has been inspired from ARESv2 of Sousa after being translated to python from C
        # The polynomial adjustment of the continuum  is made, in the original code, of order 2
        adjustValues1   = []
        adjustValues2   = []
        adjustValues3   = []
        adjustValues4   = []
        
        workingLambdas  = []
        workingLambdas2 = []        
        workingValues   = []
        workingValues2  = []
        nx_elements  = len(self.workingLambdaValues)
        print("nx_elements = ", nx_elements)
        for i in range(nx_elements):
            if ((self.workingLambdaValues[i] < inputParams.pixelExcl[1][0] or self.workingLambdaValues[i] > inputParams.pixelExcl[1][1]) and
                (self.workingLambdaValues[i] < inputParams.pixelExcl[2][0] or self.workingLambdaValues[i] > inputParams.pixelExcl[2][1]) and 
                (self.workingLambdaValues[i] < inputParams.pixelExcl[3][0] or self.workingLambdaValues[i] > inputParams.pixelExcl[3][1]) and 
                (self.workingLambdaValues[i] < inputParams.pixelExcl[4][0] or self.workingLambdaValues[i] > inputParams.pixelExcl[4][1]) ):
                workingLambdas.append(self.workingLambdaValues[i])
                workingValues.append(self.workingDataValues1[i])
        adj_coefs    = np.polyfit(workingLambdas,workingValues, 2)
        adjustValues1 = np.polyval(adj_coefs, self.workingLambdaValues)
        # splinefit1    = interpl.UnivariateSpline(workingLambdas,workingValues)
        # adjustValues1 = splinefit1(self.workingLambdaValues)
        #The array adjustValues1 has received the values of the first adjustment
        #In this new array it has to be done the exclusion of the emission zone of lambdas, gathered in adjustValues3
        for i in range(nx_elements):
            if ((self.workingLambdaValues[i] < inputParams.pixelExcl[1][0] or self.workingLambdaValues[i] > inputParams.pixelExcl[1][1]) and
                (self.workingLambdaValues[i] < inputParams.pixelExcl[2][0] or self.workingLambdaValues[i] > inputParams.pixelExcl[2][1]) and 
                (self.workingLambdaValues[i] < inputParams.pixelExcl[3][0] or self.workingLambdaValues[i] > inputParams.pixelExcl[3][1]) and 
                (self.workingLambdaValues[i] < inputParams.pixelExcl[4][0] or self.workingLambdaValues[i] > inputParams.pixelExcl[4][1]) ):
                adjustValues2.append(cpy.copy(adjustValues1[i]))
        
        # RUNS OF THE POINTS' REJECTION 
        # First Run. Discarding Points from below the continuum
        vecx = []
        vecy = []
        nvec = 0
        nx_elements2 = len(adjustValues2)
        for i in range(nx_elements2):
            #Thin out points belonging to spectral lines or noise. The remaining points supposedly belong to the continuum
            if (self.evalReject(workingValues[i], adjustValues2[i], False, rejt) and (i+1 < nx_elements2 and math.fabs(workingValues[i]-workingValues[i+1]) < strain * workingValues[i])):
                vecx.append(cpy.copy(workingLambdas[i]))
                vecy.append(cpy.copy(workingValues[i]))
                nvec+=1
        print("nvec = ", nvec)
        i = 0
        
        ## Re-arrange the arrays. In order to prepare for the second run
        for j in range(len(workingValues)):
            if(workingValues[j] == vecy[i] and i < nvec-1):
                workingValues2.append(cpy.copy(vecy[i]))
                workingLambdas2.append(cpy.copy(vecx[i]))
                i += 1
        adj_coefs    = np.polyfit(workingLambdas2,workingValues2, 2)
        adjustValues3 = np.polyval(adj_coefs, workingLambdas2)      
        # splinefit2    = interpl.UnivariateSpline(workingLambdas2,workingValues2)
        # adjustValues3 = splinefit2(workingLambdas2)        
        
        # Second Run. Discarding Points from above the continuum
        vecx2 = []
        vecy2 = []
        nvec = 0
        nx_elements2 = len(adjustValues3)
        for i in range(nx_elements2):
            #Thin out points belonging to spectral lines or noise. The remaining points supposedly belong to the continuum
            if (self.evalReject(workingValues2[i], adjustValues3[i], True, rejt) and (i+1 < nx_elements2 and math.fabs(workingValues2[i]-workingValues2[i+1]) < strain * workingValues2[i])):
                vecx2.append(cpy.copy(workingLambdas2[i]))
                vecy2.append(cpy.copy(workingValues2[i]))
                nvec+=1
        print("nvec = ", nvec)
        
        # DETERMINATION OF THE CONTINUUM
        
        #Once we have determined the relevant points to define the continuum, we must determine it making a second adjustment
        poly_fitn_coefs = np.polyfit(vecx2,vecy2, 2)
        adjustValues4    = np.polyval(poly_fitn_coefs,self.workingLambdaValues)
        # polyfitn = interpl.UnivariateSpline(vecx2,vecy2)
        # adjustValues4 = polyfitn(self.workingLambdaValues)
        j = 0
        weight = []
        for i in range(len(vecx2)):
            while(self.workingLambdaValues[j] != vecx2[i]):
                j += 1
            if(adjustValues4[j] > vecy2[i]):
                weight.append(0.35)
            else:
                weight.append(0.65)
        poly_fitn_coefs2 = np.polyfit(vecx2,vecy2, 2, rcond=None, full=False, w=weight)
        #adjustValues5    = np.polyval(poly_fitn_coefs2,self.workingLambdaValues)

        # NORMALISATION 
        # #Once we have determined the relevant points to define the continuum, we must determine the continuum making a third adjustment
        # #As in the prior run, the points of the continuum must be calculated for the WHOLE spectrum (including the exclusion zone
        # normalization by the just found continuum
        for i in range(nx_elements):
            if(i < nx_elements-1): 
                self.workingDataValues1[i] = self.workingDataValues1[i] / adjustValues4[i]
            else:
                self.workingDataValues1[i] = 1
            
        return #vecx2, vecy2, adjustValues4, adjustValues3

    def continuum_det_and_normalization_of_synth (self, rejt, strain, inputParams):
        
        # This fragment of code has been inspired from ARESv2 of Sousa after being translated to python from C
        # The polynomial adjustment of the continuum  is made, in the original code, of order 2
        adjustValues1   = []
        adjustValues2   = []
        adjustValues3   = []
        adjustValues4   = []
        adjustValues5   = []        
        workingLambdas  = []
        workingLambdas2 = []        
        workingValues   = []
        workingValues2  = []
        self.workingDataValues3 = cpy.copy(self.workingDataValues1)
        nx_elements  = len(self.workingLambdaValues)
        print("nx_elements = ", nx_elements)
        workingLambdas = cpy.copy(self.workingLambdaValues)
        workingValues  = cpy.copy(self.workingDataValues1)
        adj_coefs    = np.polyfit(workingLambdas,workingValues, 2)
        adjustValues1 = np.polyval(adj_coefs, self.workingLambdaValues)
        #The array adjustValues1 has received the values of the first adjustment
        
        # RUNS OF THE POINTS' REJECTION 
        # Discarding Points from below the continuum
        vecx = []
        vecy = []
        nvec = 0
        nx_elements2 = len(adjustValues1)
        for i in range(nx_elements2):
            #Thin out points belonging to spectral lines or noise. The remaining points supposedly belong to the continuum
            if (self.evalReject(workingValues[i], adjustValues1[i], False, rejt) and 
                (i+1 < nx_elements2 and math.fabs(workingValues[i]-workingValues[i+1]) < strain * workingValues[i])):
                vecx.append(cpy.copy(workingLambdas[i]))
                vecy.append(cpy.copy(workingValues[i]))
                nvec+=1
        print("nvec = ", nvec)
        i = 0
        
        ## Re-arrange the arrays. In order to prepare for the second run
        for j in range(len(workingValues)):
            if(workingValues[j] == vecy[i] and i < nvec-1):
                workingValues2.append(cpy.copy(vecy[i]))
                workingLambdas2.append(cpy.copy(vecx[i]))
                i += 1
        adj_coefs    = np.polyfit(workingLambdas2,workingValues2, 2)
        adjustValues3 = np.polyval(adj_coefs, workingLambdas2)      
        
        # DETERMINATION OF THE CONTINUUM
        
        #Once we have determined the relevant points to define the continuum, we must determine it making a second adjustment
        poly_fitn_coefs = np.polyfit(vecx,vecy, 2)
        adjustValues4    = np.polyval(poly_fitn_coefs,self.workingLambdaValues)
        j = 0
        weight = []
        for i in range(len(vecx)):
            while(self.workingLambdaValues[j] != vecx[i]):
                j += 1
            if(adjustValues4[j] > vecy[i]):
                weight.append(0.35)
            else:
                weight.append(0.65)
        poly_fitn_coefs2 = np.polyfit(vecx,vecy, 2, rcond=None, full=False, w=weight)
        adjustValues5    = np.polyval(poly_fitn_coefs2,self.workingLambdaValues)

        # NORMALISATION 
        # #Once we have determined the relevant points to define the continuum, we must determine the continuum making a third adjustment
        # #As in the prior run, the points of the continuum must be calculated for the WHOLE spectrum (including the exclusion zone
        # normalization by the just found continuum
        for i in range(nx_elements):
            if(i < nx_elements-1): 
                self.workingDataValues1[i] = self.workingDataValues1[i] / adjustValues4[i]
            else:
                self.workingDataValues1[i] = 1
            
        return #vecx2, vecy2, adjustValues4, adjustValues3
    
    def getMaxValue(self, workingValues):
        max_value = 0.0
        for it in range(len(workingValues)-1):
            if max_value < workingValues[it]:
                max_value = workingValues[it]
                lambLoc = it
        return max_value, lambLoc
    
    def getMinValue(self, workingValues):
        min_value = 9e99
        for it in range(len(workingValues)-1):
            if min_value > workingValues[it]:
                min_value = workingValues[it]
                lambLoc = it
        return min_value, lambLoc
    #Calculate Fbol in unprocessed data file by means of Trapezoidadl Rule
    def getRawFbol(self): 
        lendv = len(self.dataValues)-2
        dv1    = self.dataValues[0]
        dvlast = self.dataValues[lendv+1]
        Fbol = (dv1 + dvlast)/2
        for i in range(lendv):
            Fbol += self.dataValues[i+1]
        return Fbol

    ###########################################################################
    # Calculates de Mean Value of the flux in a wavelength interval
    # lambda_start: initial value of lambda interval
    # lambda_end  : final value of the lambda interval
    # inputObject : FITSObject containing the flux values
    ###########################################################################
    def CalculateMeanF_lambda(self, lambda_start, lambda_end):
        f_lambda_mean = 0
        count = 0
        for iterL in range (len(self.workingLambdaValues)):
            if (self.workingLambdaValues[iterL] <= lambda_start):
                continue
            elif (self.workingLambdaValues[iterL] > lambda_start and self.workingLambdaValues[iterL] <= lambda_end):
                f_lambda_mean += self.workingDataValues3[iterL]
                count += 1
            else:
                break
        return f_lambda_mean/count
    ###########################################################################
    
