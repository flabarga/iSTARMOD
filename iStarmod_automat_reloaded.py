#=======================================================================
## Name:     iSTARMOD_AUTOMAT_RELOADED
# Filename: 
# Type:     program
################################################################
# Language: Python 3.4
# Purpose:  Synthesize a composite spectrum from input templates and 
#           object spectrum, outputs composite and residual spectrum.
################################################################
# Call line:     none
# Subroutines:   starot (ndata, datain, vsini, lambda, scale, datout)
#                fshift (npts, datain, shift, datout)
#                readConfigFile(configfile)
#                plot(inObj, inputParams)
################################################################
# Author:           Orig.: Sam Barden 1984 / Final: Fernando Labarga 2017
# Modifications:    dph = David Huenemoerder     adw = Alan Welty
# Migration to Python  Fernando Labarga 2017
################################################################
#   15-feb-17  ffla Migration from FORTRAN code to python
#              Correction of the features of the kernel/broadening function in starot
################################################################################################
################################################################################################

import os 
import matplotlib.pyplot as plt
from iStarmod_fits import * 
import iStarmod_tools as tools
import support as supp
import starot as starot
import fshift as fshift
import writeFitsSpec as WFS
import copy as cpy
import math
from pathlib import Path
import iStarmod_rvvalues as rv


def starmod(spcfFile, debugging = False):

    ##################################################################################3
    #READING THE SPECIFICATION FILE
    ##################################################################################3
    try:
        inputParams = tools.readConfigFile(spcfFile)
    except FileNotFoundError:
        print("INPUT PARAMETERS FILE NOT FOUND. CHECK FILENAME OR PATH")
        return
        
    #ensure that the weights are correctly specified
    inputParams.checkWeights()
    tst = []
    #############################
    #Implementing the specification var/fix.
    #############################
    #If the parameter is 'fix', then tst[n] must be zero
    #When parameter is 'var', is assigned to tst[n] the input parameter value
    tst.append([])
    if inputParams.PRMRad[1]:     tst[0].append(0)
    else:                         tst[0].append(inputParams.PRMRad[0])
    if inputParams.PRMRot[1]:     tst[0].append(0)
    else:                         tst[0].append(inputParams.PRMRot[0])
    if inputParams.PRMWeight[1]:  tst[0].append(0)
    else:                         tst[0].append(inputParams.PRMWeight[0])
    print ("tst[0]: " + str(tst[0]))
    tst.append([])
    if inputParams.SECRad[1]:     tst[1].append(0)
    else:                         tst[1].append(inputParams.SECRad[0])
    if inputParams.SECRot[1]:     tst[1].append(0)
    else:                         tst[1].append(inputParams.SECRot[0])
    if inputParams.SECWeight[1]:  tst[1].append(0)
    else:                         tst[1].append(inputParams.SECWeight[0])
    print ("tst[1]: " + str(tst[1]))
    tst.append([])
    if inputParams.TERRad[1]:     tst[2].append(0)
    else:                         tst[2].append(inputParams.TERRad[0])
    if inputParams.TERRot[1]:     tst[2].append(0)
    else:                         tst[2].append(inputParams.TERRot[0])
    if inputParams.TERWeight[1]:  tst[2].append(0)
    else:                         tst[2].append(inputParams.TERWeight[0])
    print ("tst[2]: " + str(tst[2]))
    nstar = 3
    if ('NONE' in inputParams.tertiaryStarfits) or ('NONE' in inputParams.secondaryStarfits):
        nstar = 2
        tst[2][0] = 0
        tst[2][1] = 0
        tst[2][2] = 0
        inputParams.TERWeight[0] = 0.0
        inputParams.SECWeight[0] = float('{:5.3f}'.format(1.0 - inputParams.PRMWeight[0])) #To avoid floating point rounding errors
        
    if ('NONE' in inputParams.secondaryStarfits):
         nstar = 1
         tst[1][0] = 0
         tst[1][1] = 0
         tst[1][2] = 0
         inputParams.SECWeight[0] = 0.0
         inputParams.PRMWeight[0] = 1.0
         
    nra = tst[0][0] + tst[1][0] + tst[2][0]
    nro = tst[0][1] + tst[1][1] + tst[2][1]
    nwt = tst[0][2] + tst[1][2] + tst[2][2]
    pwt = inputParams.PRMWeight[0] 
    swt = inputParams.SECWeight[0] 
    twt = inputParams.TERWeight[0] 
    #ntest = nro + nra + nwt
    skp2 = []
    for i in range(10):
        skp2.append(0)
    #for i in range(5):
    #    skp.append(inputParams.pixelExcl[i][0])
    #    skp.append(inputParams.pixelExcl[i][1])#inputParams.pixelExcl[2]])
        
    #   Initialize best param array
    bstprm = []
    bstprm.append([])
    bstprm[0].append(inputParams.PRMRad[0])
    bstprm[0].append(inputParams.PRMRot[0])
    bstprm[0].append(inputParams.PRMWeight[0])
    bstprm.append([])
    bstprm[1].append(inputParams.SECRad[0])
    bstprm[1].append(inputParams.SECRot[0])
    bstprm[1].append(inputParams.SECWeight[0])
    bstprm.append([])
    bstprm[2].append(inputParams.TERRad[0])
    bstprm[2].append(inputParams.TERRot[0])
    bstprm[2].append(inputParams.TERWeight[0])
    
    xPrimStarmean = 0.0 
    xSecStarmean = 0.0
    xTerStarmean = 0.0
    
    deltaStepAvg = 0.0
    
        
    ######################## To check if RV Values have been calculated before #################################
    if (nstar == 1):
        rv_values_available = False
        rv_values = rv.RVValues()
        if ("nir-tac" not in inputParams.workpath) and ("NIR-TAC" not in inputParams.workpath) :
            filenameRV1 = (inputParams.workpath[:-4] + "RES\\" + "rvvalues.dat")
            filenameRV2 = (inputParams.workpath[:-4] + "RES\\")
        else:
            filenameRV1 = (inputParams.workpath[:-8] + "RES\\" + "rvvalues.dat")
            filenameRV2 = (inputParams.workpath[:-8] + "RES\\")
        if os.access(filenameRV1, os.F_OK):
            rv_values.readFile(filenameRV1)
            rv_values_available = True
            print ("filenameRV:" , filenameRV1)
        else:
            if not os.access(filenameRV2, os.F_OK):
                os.mkdir(filenameRV2) 
                print("Here!!!")       
            # rv_outfile = open(filenameRV2, 'a')
            print ("outfilenameRV:" , filenameRV2)
    ############################################################################################################
       
    ##################################################################################3
    #READING and plotting THE OBJECT FITS FILE
    ##################################################################################3
    dirPath = Path(inputParams.workpath)
    print("Path = " + str(dirPath))
    listFITSFiles = list(dirPath.glob("car-*.fits"))
    # To cope with specifying not CARMENES files
    if len(listFITSFiles) == 0 and '*' not in inputParams.objectInputFits:
        listFITSFiles.append(str(inputParams.objectInputFits))
    elif len(listFITSFiles) == 0 :
        listFITSFiles = list(dirPath.glob("*.fits"))
    
    if len(listFITSFiles) == 0:
        listFITSFiles = list(dirPath.glob("*.dat"))

    print(listFITSFiles)
    for objectInputFITSfromList in listFITSFiles:
        
        ##################  Newly added lines 18/04/2025  ###############################################
        inputParams.checkWeights()
        nwt = tst[0][2] + tst[1][2] + tst[2][2]
        pwt = inputParams.PRMWeight[0] 
        swt = inputParams.SECWeight[0] 
        twt = inputParams.TERWeight[0] 
        ##################################################################################################
        
        f, (subPlot1,subPlot2) = plt.subplots(2, figsize = (10,7), dpi = 150, sharex = True)#, sharey = True)
        subPlot1.clear()
        subPlot2.clear()    
        
        objectInputFITS = str(objectInputFITSfromList)
        inputObject = tools.readFitsFile(".\\" + objectInputFITS, inputParams.spectraFormatAperture)
        
        print("objectinputfits : ", objectInputFITS, ", MJD-OBS = ", inputObject.date)
        

        #Reading and preprocessing Input file
        #evaluate if 'pixel range' is given in pixel or in Angstroms
        #it is assumed that the start and end pixels specified belong to the same order 
        xlo,xhi, skp = inputObject.selectPixelRange(inputParams, False, subPlot2, False)
        inputObject.nanSanityCheck()
        print("SNR = ", inputObject.SNR)
        # The next lines, following (at least partially) the algorithms set out in Sousa et al., A&A 577, A67 (2015)
        # It is needed for CARMENES, and from other sources, spectra
        ##################################### NORMALIZATION ############
        rejt = 1. - (1/inputObject.SNR)
        if (inputObject.SNR < 80):
            strain = .1
        else:
            strain = .15
        inputObject.continuum_det_and_normalization(rejt,strain, inputParams)
        ##################################################################3
        ##################################### REBINNING     ############
        if inputObject.bCarmenes == True:
            if inputParams.deltaStep != -1.0:
                inputObject.regularizeDeltaStepAndTrim(inputParams.deltaStep, inputParams.nrpoints)
                inputObject.resampleDataValues(None)
            else:
                bDeltaStepChecked, deltaStepAvg = inputObject.checkDeltaStep()
                if bDeltaStepChecked:
                    inputObject.regularizeDeltaStep(deltaStepAvg)
                    inputObject.resampleDataValues(None)
        else:
        #WHAT IF IT IS THE CASE THAT THE INPUT OBJECT IS NOT A CARMENES SPECTRUM? IT IS NEEDED TO ASSIGN A VALUE TO deltaStepAvg
        #WE ARE ASSUMING (maybe it is too much to assume) THAT THE SPECTRUM HAS THE LAMBDA VALUES REGULARLY SPACED
            bDeltaStepChecked, deltaStepAvg = inputObject.checkDeltaStep()
            if bDeltaStepChecked:
                inputObject.regularizeDeltaStep(deltaStepAvg)
                inputObject.resampleDataValues(None)
            else:
                inputParams.deltaStep = deltaStepAvg
        inputObject.workingDataValues2 = cpy.copy(inputObject.workingDataValues1)    
        inputObject.workingDataValues3 = cpy.copy(inputObject.workingDataValues1)
        nptsInputObj = len(inputObject.workingDataValues1)
        inputParams.finalLambda = inputObject.finalLambda #inputParams.deltaStep   = inputObject.deltaStep
        inputParams.numberOfWorkingDataValues = nptsInputObj
        ##################################################################3

        #the way the spectrum is plotted depends on the format
        #inputObject.plot(inputParams)
        xObjmean = supp.gtmean(nptsInputObj,inputObject.workingDataValues1,0,nptsInputObj-1,skp)
        midLambda = inputObject.initialLambda + inputParams.deltaStep * len(inputObject.workingLambdaValues) / 2.0
        print(nptsInputObj, inputObject.workingLambdaValues[len(inputObject.workingLambdaValues)-1], midLambda)
        cnst = (midLambda/299792.458)/inputParams.deltaStep     #Conversion Factor from Doppler Effect Formula

        #Reading Prymary (Reference) Star spectrum
        primaryStar = tools.readFitsFile(inputParams.primaryStarfits, inputParams.spectraFormatAperture)
        # primaryStar.normalizeSpectrum(inputParams, False)
        if primaryStar != None:     
            primaryStar.selectPixelRange(inputParams, True, subPlot2, False)
            primaryStar.nanSanityCheck()
            print("SNR = ", primaryStar.SNR)
            ##################################### NORMALIZATION ############
            rejt = 1. - (1/primaryStar.SNR)
            if (primaryStar.SNR < 80):
                strain = .1
            else:
                strain = .15
            primaryStar.continuum_det_and_normalization(rejt,strain, inputParams)
            ###################################################################3
            #################################### REBINNING #################
            if primaryStar.bCarmenes == True: #or inputParams.diffFormatPrimary == True:
                if inputParams.deltaStep != -1.0:
                    primaryStar.regularizeDeltaStepAndTrim(inputParams.deltaStep, inputParams.nrpoints)
                    primaryStar.resampleDataValues(None)
                else:
                    bDeltaStepChecked1, deltaStepAvg1 = primaryStar.checkDeltaStep()
                    if bDeltaStepChecked or (deltaStepAvg1 != deltaStepAvg):
                        primaryStar.regularizeDeltaStep(deltaStepAvg)
                        primaryStar.resampleDataValues(inputObject.workingLambdaValues)
            else:
                #WHAT IF IT IS THE CASE THAT THE PRIMARY STAR OBJECT IS NOT A CARMENES SPECTRUM? IT IS NEEDED TO ASSIGN A VALUE TO deltaStepAvg
                #WE ARE ASSUMING (maybe it is too much to assume) THAT THE SPECTRUM HAS THE LAMBDA VALUES REGULARLY SPACED
                bDeltaStepChecked, deltaStepAvg = primaryStar.checkDeltaStep()
                if bDeltaStepChecked:
                    primaryStar.regularizeDeltaStep(deltaStepAvg)
                    primaryStar.resampleDataValues(None)

            nptsPrimSt = len(primaryStar.workingDataValues1)
            primaryStar.workingDataValues2 = cpy.copy(primaryStar.workingDataValues1)    
            primaryStar.workingDataValues3 = cpy.copy(primaryStar.workingDataValues1)
            xPrimStarmean = supp.gtmean(nptsPrimSt,primaryStar.workingDataValues1,0,nptsPrimSt-1,skp)
            ###################################################################3


        if 'NONE' not in inputParams.secondaryStarfits:
            secondaryStar = tools.readFitsFile(inputParams.secondaryStarfits, inputParams.spectraFormatAperture)
            secondaryStar.selectPixelRange(inputParams, True, subPlot2, False)
            ##################################### NORMALIZATION ############
            rejt = 1. - (1/secondaryStar.SNR)
            if (secondaryStar.SNR < 80):
                strain = .1
            else:
                strain = .15
            secondaryStar.continuum_det_and_normalization(rejt,strain, inputParams)
            ###################################################################3
            if secondaryStar.bCarmenes == True or inputParams.diffFormatPrimary == True:
                bDeltaStepChecked1, deltaStepAvg1 = secondaryStar.checkDeltaStep()
                if bDeltaStepChecked or (deltaStepAvg1 != deltaStepAvg):
                    secondaryStar.regularizeDeltaStep(deltaStepAvg1)
                    secondaryStar.resampleDataValues(inputObject.workingLambdaValues)
            else:
                #WHAT IF IT IS THE CASE THAT THE PRIMARY STAR OBJECT IS NOT A CARMENES SPECTRUM? IT IS NEEDED TO ASSIGN A VALUE TO deltaStepAvg
                #WE ARE ASSUMING (maybe it is too much to assume) THAT THE SPECTRUM HAS THE LAMBDA VALUES REGULARLY SPACED
                bDeltaStepChecked, deltaStepAvg = primaryStar.checkDeltaStep()
                if bDeltaStepChecked:
                    secondaryStar.regularizeDeltaStep(deltaStepAvg)
                    secondaryStar.resampleDataValues(None)
            secondaryStar.workingDataValues2 = cpy.copy(secondaryStar.workingDataValues1)    
            secondaryStar.workingDataValues3 = cpy.copy(secondaryStar.workingDataValues1)
            nptsSecSt = len(secondaryStar.workingDataValues1)
            xSecStarmean = supp.gtmean(nptsSecSt,secondaryStar.workingDataValues1,0,nptsSecSt-1,skp)
            #the way the spectrum is plotted depends on the format
            #plot(secondaryStar, inputParams)
        
        if 'NONE' not in inputParams.tertiaryStarfits:
            tertiaryStar = tools.readFitsFile(inputParams.tertiaryStarfits, inputParams.spectraFormatAperture)
            # tertiaryStar = tools.normalizeSpectrum(inputparams, False)
            tertiaryStar.selectPixelRange(inputParams, True, subPlot2, False)
            ##################################### NORMALIZATION ############
            rejt = 1. - (1/tertiaryStar.SNR)
            if (tertiaryStar.SNR < 80):
                strain = .1
            else:
                strain = .15
            tertiaryStar.continuum_det_and_normalization(rejt,strain, inputParams)
            ###################################################################3
            if tertiaryStar.bCarmenes == True: #or inputParams.diffFormatPrimary == True:
                bDeltaStepChecked1, deltaStepAvg1 = tertiaryStar.checkDeltaStep()
                if bDeltaStepChecked or (deltaStepAvg1 != deltaStepAvg):
                    tertiaryStar.regularizeDeltaStep(deltaStepAvg1)
                    tertiaryStar.resampleDataValues(inputObject.workingLambdaValues)
            else:
            #WHAT IF IT IS THE CASE THAT THE PRIMARY STAR OBJECT IS NOT A CARMENES SPECTRUM? IT IS NEEDED TO ASSIGN A VALUE TO deltaStepAvg
            #WE ARE ASSUMING (maybe it is too much to assume) THAT THE SPECTRUM HAS THE LAMBDA VALUES REGULARLY SPACED
                bDeltaStepChecked, deltaStepAvg = primaryStar.checkDeltaStep()
                if bDeltaStepChecked:
                    tertiaryStar.regularizeDeltaStep(deltaStepAvg)
                    tertiaryStar.resampleDataValues(None)
            tertiaryStar.workingDataValues2 = cpy.copy(tertiaryStar.workingDataValues1)    
            tertiaryStar.workingDataValues3 = cpy.copy(tertiaryStar.workingDataValues1)
            nptsTerSt = len(tertiaryStar.workingDataValues1)
            xTerStarmean = supp.gtmean(nptsTerSt,tertiaryStar.workingDataValues1,0,nptsTerSt-1,skp)
            #the way the spectrum is plotted depends on the format
            #plot(tertiaryStar, inputParams)
        
        ######################## NEW LINES 16/07 MODIFY THE WAY THE SPECTRA ARE SHIFTED ####
        ###################### New lines 10/05/2020
        if inputParams.deltaStep != -1.0:
            deltaStepProc = inputParams.deltaStep
        elif deltaStepAvg != 0 :
            deltaStepProc = deltaStepAvg
        else:
            deltaStepProc = deltaStepAvg1
        cnst = (midLambda/299792.458)/deltaStepProc
        ##########################################   10/07/2020  
        ##################################################################################3

        ##################################################################################3
        #PROCESSING THE SPECTRA
        ##################################################################################3
        ndata = 64
        npts  = len(inputObject.workingDataValues1)
        while (npts > ndata):
            ndata = ndata * 2
        zeroes = []
        for i in range(nptsInputObj):
            zeroes.append(0.0)
        
        ################################-----------------------------------------------------
        #Let the modeling begin!
        #Iteration loop begins here
        ###############################----------------------------------------------------3
        num = 123456
        print (' ')
        print ('####################################################')
        print (' The input parameters are:  ')
        print ("\t\t" + "PRI" + "\t" + "SEC" + "\t" + "TER")
        print (' V_rad   = ' + "\t" + str(inputParams.PRMRad[0]) + "\t" + str(inputParams.SECRad[0]) + "\t" + str(inputParams.TERRad[0]))
        print (' V_rot   = ' + "\t" + str(inputParams.PRMRot[0]) + "\t" + str(inputParams.SECRot[0]) + "\t" + str(inputParams.TERRot[0]))
        print (' Weight  = ' + "\t" + str('{:4.3f}'.format(inputParams.PRMWeight[0])) + "\t" + str('{:4.3f}'.format(inputParams.SECWeight[0])) + "\t" + str('{:4.3f}'.format(inputParams.TERWeight[0])))
        print ('-----------------------------------------------')
        print (' ')
        print (' Object    = ' + objectInputFITS)
        print (' Primary   = ' + str(inputParams.primaryStarfits))
        print (' Secondary = ' + str(inputParams.secondaryStarfits))
        print (' Tertiary  = ' + str(inputParams.tertiaryStarfits))
        print (' Order     = ' + str(inputParams.spectraFormatAperture))
        print 
        
        for niter in range(inputParams.numIter):
            ######################################################################################
            #                Set-up of the initial guess                                          
            ######################################################################################
            print ("ITERATION #" + str(niter+1))
            # print ("Weight parameter [primary star]: " + str(inputParams.PRMWeight[0]))
            #print primaryStar.workingDataValues3
            #print "//////pwt = ", pwt, "midLambda = ", midLambda, "cnst: ", cnst, "deltaStep = ", inputParams.deltaStep
            supp.scale(nptsPrimSt,primaryStar.workingDataValues1,pwt,primaryStar.workingDataValues3)
            
            ###############################################Rotationally broaden it
            # (workingDataValues2 holds scaled & broad prim)
            # print ("Rotation parameter [primary star]: " + str(inputParams.PRMRot[0]))
            primaryStar.workingDataValues2 = starot.starot(nptsPrimSt, primaryStar.workingDataValues3, inputParams.PRMRot[0],
                                                            midLambda, deltaStepProc)# ,primaryStar.workingDataValues2)
            #primaryStar.workingDataValues2 = rotBroad(primaryStar.workingLambdaValues, primaryStar.workingDataValues3, 0.6, inputParams.PRMRot[0])
            ###############################################Shift it in radial velocity
            # print ("Radial Velocity parameter [primary star]: " + str(inputParams.PRMRad[0]))
            pdch = inputParams.PRMRad[0] * cnst
            primaryStar.workingDataValues3 = fshift.fshift(nptsPrimSt,primaryStar.workingDataValues2,pdch)
            
            ############################################### Same for secondary:
            if (nstar > 1):
                
                supp.scale(npts,secondaryStar.workingDataValues1,swt,secondaryStar.workingDataValues3)
                # print ("Rotation parameter [secondary star]: " + str(inputParams.SECRot[0]))
                secondaryStar.workingDataValues2 = starot.starot(nptsSecSt, secondaryStar.workingDataValues3, inputParams.SECRot[0], 
                                                                midLambda, deltaStepProc)#, secondaryStar.workingDataValues2 )
                #secondaryStar.workingDataValues2 = rotBroad(secondaryStar.workingLambdaValues, secondaryStar.workingDataValues3, 0.6, inputParams.SECRot[0])
                # print ("Radial Velocity parameter [secondary star]: " + str(inputParams.SECRad[0]))
                sdch = inputParams.SECRad[0] * cnst
                secondaryStar.workingDataValues3 = fshift.fshift(nptsSecSt,secondaryStar.workingDataValues2,sdch)
    
            ############################################### And for the tertiary:
            if (nstar == 3):
                supp.scale(nptsTerSt,tertiaryStar.workingDataValues1,twt,tertiaryStar.workingDataValues3)
                # print ("Rotation parameter [tertiary star]: " + str(inputParams.TERRot[0]))
                tertiaryStar.workingDataValues2 = starot.starot(nptsTerSt, tertiaryStar.workingDataValues3,inputParams.TERRot[0], 
                                                                midLambda, deltaStepProc)#,tertiaryStar.workingDataValues2 )
                #tertiaryStar.workingDataValues2 = rotBroad(tertiaryStar.workingLambdaValues, tertiaryStar.workingDataValues3, 0.6, inputParams.TERRot[0])
                # print ("Radial Velocity parameter [tertiary star]: " + str(inputParams.TERRad[0]))
                tdch = inputParams.TERRad[0] * cnst
                tertiaryStar.workingDataValues3 = fshift.fshift(nptsTerSt,tertiaryStar.workingDataValues2,sdch)
    
    
            #////////////////////////////////////////////////////////////////////////////////////
            #-----------------------------------------------------------------------------------
            #Add to get first model,calculate the sum of the squared residuals (ssr) and
            #set initial ssr as minimum ssr
            #-----------------------------------------------------------------------------------
            #///////////////////////////////////////////////////////////////////////////////////
            xmean = xPrimStarmean*inputParams.PRMWeight[0] + xSecStarmean*inputParams.SECWeight[0] + xTerStarmean*inputParams.TERWeight[0]
            #minssr = 10.0
            if (niter == 0):
                if nstar == 1: #IF nstar = 1 AND pwt = 1 THERE IS NO NEED TO MAKE THE CALL TO sumthm
                    inputObject.workingDataValues2 = cpy.copy(primaryStar.workingDataValues3)
                elif nstar == 2:
                    inputObject.workingDataValues2 = supp.sumthm(npts,primaryStar.workingDataValues3,secondaryStar.workingDataValues3,zeroes)    
                else: #nstar == 3:
                    inputObject.workingDataValues2 = supp.sumthm(npts,primaryStar.workingDataValues3,secondaryStar.workingDataValues3,tertiaryStar.workingDataValues3)
                inputObject.workingDataValues2 = supp.scale2 (npts,inputObject.workingDataValues2,xObjmean/xmean)
                ssr,inputParams.cnt = supp.sumres(npts,inputObject.workingDataValues1,inputObject.workingDataValues2,
                                                0,npts-1,skp)
                minssr_init = ssr
            #end if
                    
            #####################################################
            #Start searching parameter space
            #####################################################
            #*** Find best radial velocities
            #####################################################
                        
            #----------------- If there exists a RV value available, use it ----------------------------------------------
            if nstar == 1 and rv_values_available :
                print("****************************************************    searching...")
                rv_value_iter = rv_values.find_mjd(inputObject.date)
                if rv_value_iter == None:
                    tst[0][0] = 1
                    rv_values_available = False
                elif nstar == 1:
                    print("************************************************    found...", rv_value_iter[1])
                    tst[0][0] = 0
                    inputParams.PRMRad[0] = rv_value_iter[1]
                    pdch = inputParams.PRMRad[0] * cnst
                    primaryStar.workingDataValues3 = fshift.fshift(nptsPrimSt, primaryStar.workingDataValues2, pdch)
            #------------------------------------------------------------------------------------------------------------
            #Shift the Primary
            if (tst[0][0] != 0):
                # print ("    Shifting the primary" )
                #Shift p2 and store in s3
                minssr = minssr_init
                #print ("minssr = ", str(minssr))
                parm = inputParams.PRMRad[0] - 50.0
                deltaparam = 5.0
                for j in range(8):
                    for i in range(20):    
                        pdch = parm * cnst
                        primaryStar.workingDataValues3 = fshift.fshift(nptsPrimSt,primaryStar.workingDataValues2,pdch)
                        if nstar == 2:
                            inputObject.workingDataValues2 = supp.sumthm(nptsPrimSt,primaryStar.workingDataValues3,
                                                                        secondaryStar.workingDataValues3,zeroes)
                        elif nstar == 3:
                            inputObject.workingDataValues2 = supp.sumthm(nptsPrimSt,primaryStar.workingDataValues3,secondaryStar.workingDataValues3,tertiaryStar.workingDataValues3)
                        else:
                            inputObject.workingDataValues2 = cpy.copy(primaryStar.workingDataValues3)
                        
                        supp.scale2 (nptsInputObj,inputObject.workingDataValues2,xObjmean/xmean)
                        
                        ###### Squared residuals and selection of bestparm (minimum)
                        ssr,inputParams.cnt = supp.sumres(nptsInputObj,inputObject.workingDataValues1,inputObject.workingDataValues2,
                                            0,nptsInputObj-1,skp)
                        if (ssr <= minssr):

                            minssr = ssr
                            bstprm[0][0] = parm
                        #end if
                        
                        parm = parm + deltaparam
                    #continue
                    deltaparam = deltaparam/10
                    parm = bstprm[0][0] - (10 * deltaparam)
                #continue
                
                # Once the final (best) parameter value has been obtained, the final spectrum working values are saved to primaryStar.workingDataValues3
                inputParams.PRMRad[0] = bstprm[0][0]
                pdch = inputParams.PRMRad[0] * cnst
                primaryStar.workingDataValues3 = fshift.fshift(nptsPrimSt, primaryStar.workingDataValues2, pdch)
            #end if

            #Shift the Secondary
            if (tst[1][0] != 0):
                # print ("    Shifting the Secondary")
                minssr = minssr_init
                parm = inputParams.SECRad[0] - 50.0
                deltaparam = 5.0
                for j in range(8):
                    for i in range(21):
                        sdch = parm * cnst
                        secondaryStar.workingDataValues3 = fshift.fshift(nptsSecSt,secondaryStar.workingDataValues2,sdch)
                        
                        if nstar == 2:
                            inputObject.workingDataValues2 = supp.sumthm(npts,primaryStar.workingDataValues3, secondaryStar.workingDataValues3,zeroes)
                        elif nstar == 3:
                            inputObject.workingDataValues2 = supp.sumthm(npts,primaryStar.workingDataValues3, secondaryStar.workingDataValues3,tertiaryStar.workingDataValues3)
                        
                        supp.scale2 (nptsSecSt,inputObject.workingDataValues2,xObjmean/xmean)
                        
                        ###### Squared residuals and selection of bestparm (minimum)
                        ssr,inputParams.cnt = supp.sumres(nptsInputObj,inputObject.workingDataValues1,inputObject.workingDataValues2,0,nptsInputObj-1,skp)
                        if (ssr < minssr):
                            minssr = ssr
                            bstprm[1][0]  = parm
                        #end if
                        parm = parm + deltaparam
                        #continue
                    deltaparam = 0.1 * deltaparam
                    parm = bstprm[1][0] - (10.0 * deltaparam)
                    #continue
                
                # Once the final (best) parameter value has been obtained, the final spectrum working values are saved to secondaryStar.workingDataValues3
                inputParams.SECRad[0] = bstprm[1][0]
                sdch = inputParams.SECRad[0] * cnst
                secondaryStar.workingDataValues3 = fshift.fshift(len(secondaryStar.workingDataValues2), secondaryStar.workingDataValues2, pdch)
            #end if  
                  
            #Shift the Tertiary
            if (tst[2][0] != 0):
                # print ("Shifting the tertiary")
                minssr = minssr_init
                parm = inputParams.TERRad[0] - 50.0
                deltaparam = 5.0
                for j in range(7):
                    for i in range(21):
                        tdch = parm * cnst
                        tertiaryStar.workingDataValues3 = fshift.fshift(nptsTerSt,tertiaryStar.workingDataValues2,tdch)
                        
                        inputObject.workingDataValues2 = supp.sumthm(npts,primaryStar.workingDataValues3, secondaryStar.workingDataValues3,tertiaryStar.workingDataValues3)
                        
                        supp.scale2 (npts,inputObject.workingDataValues2,xObjmean/xmean)
                        sr,ssr,inputParams.cnt = supp.sumres(nptsInputObj,inputObject.workingDataValues1,inputObject.workingDataValues2, 0,nptsInputObj-1,skp)
                        if (ssr < minssr):
                            minssr = ssr
                            bstprm[2][0]  = parm
                        #end if
                        parm = parm + deltaparam
                        #continue
                    deltaparam = 0.1 * deltaparam
                    parm = bstprm[2][0] - (10.0 * deltaparam)
                    #continue
                inputParams.TERRad[0] = bstprm[2][0]
                tdch = inputParams.TERRad[0] * cnst
                tertiaryStar.workingDataValues3 = fshift.fshift(nptsTerSt, tertiaryStar.workingDataValues2, tdch)
            #end if 


            #####################################################
            #*** Find the best vsini value (rotational velocity including the effect of the unknown inclination of the axis)
            #####################################################
            # print ("Finding the best rotational velocities")
            minssr = minssr_init

            #Rotate the primary
            if (tst[0][1] != 0):
                # print ("    Rotating the primary")
                #print ("minssr = ", str(minssr))
                supp.scale(nptsPrimSt,primaryStar.workingDataValues1,inputParams.PRMWeight[0], primaryStar.workingDataValues3)
                pdch = inputParams.PRMRad[0] * cnst
                primaryStar.workingDataValues2 = fshift.fshift(nptsPrimSt,primaryStar.workingDataValues3,pdch)
                parm = inputParams.PRMRot[0] - 50.0
                deltaparam = 5.0
                for j in range(7):
                    for i in range(20):
                        #print i
                        if (parm >= 0.0):#if <0 go to 509
                            primaryStar.workingDataValues3 = starot.starot(nptsPrimSt, primaryStar.workingDataValues2, parm, midLambda, deltaStepProc)#, primaryStar.workingDataValues3)
                            if nstar == 1:
                                inputObject.workingDataValues2 = cpy.copy(primaryStar.workingDataValues3)
                            elif nstar == 2:
                                inputObject.workingDataValues2 = supp.sumthm(nptsPrimSt,primaryStar.workingDataValues3, secondaryStar.workingDataValues3, zeroes)
                            else: #nstar == 3:
                                inputObject.workingDataValues2 = supp.sumthm(nptsPrimSt,primaryStar.workingDataValues3, secondaryStar.workingDataValues3,tertiaryStar.workingDataValues3)
                            
                            supp.scale2 (npts,inputObject.workingDataValues2,xObjmean/xmean)
                            ssr,inputParams.cnt = supp.sumres(npts,inputObject.workingDataValues1,inputObject.workingDataValues2, 0,npts-1,skp)
                            if (ssr <= minssr):
                                minssr = ssr
                                bstprm[0][1] = parm
                        parm = parm + deltaparam
                    #continue #510
                    deltaparam = deltaparam/10
                    parm = bstprm[0][1] - (10 * deltaparam)
                #continue #515
                inputParams.PRMRot[0] = bstprm[0][1]
                primaryStar.workingDataValues3 = starot.starot(nptsPrimSt, primaryStar.workingDataValues2,inputParams.PRMRot[0],
                                                            midLambda, deltaStepProc)#, primaryStar.workingDataValues3)
                #primaryStar.workingDataValues3 = rotBroad(primaryStar.workingLambdaValues, primaryStar.workingDataValues2,
                                                            #0.6, inputParams.PRMRot[0])
            #end if

            #Rotate the secundary 
            minssr = minssr_init
            if (tst[1][1] != 0):
                # print ("    Rotating the secondary")
                supp.scale(nptsSecSt,secondaryStar.workingDataValues1,inputParams.SECWeight[0],secondaryStar.workingDataValues3)
                sdch = inputParams.SECRad[0] * cnst
                secondaryStar.workingDataValues2 = fshift.fshift(nptsPrimSt,secondaryStar.workingDataValues3,sdch)
                parm = inputParams.SECRot[0] - 50.0
                deltaparam = 5.0
                for j in range(3):
                    for i in range(20):
                        if (parm >= 0.0): 
                            secondaryStar.workingDataValues3 = starot.starot(nptsSecSt, secondaryStar.workingDataValues2, parm, midLambda, deltaStepProc)#, secondaryStar.workingDataValues3)
                            
                            if nstar == 2:
                                    inputObject.workingDataValues2 = supp.sumthm(nptsSecSt,primaryStar.workingDataValues3, secondaryStar.workingDataValues3,zeroes)
                            elif nstar == 3:
                                    inputObject.workingDataValues2 = supp.sumthm(nptsSecSt,primaryStar.workingDataValues3, secondaryStar.workingDataValues3,tertiaryStar.workingDataValues3)
                            supp.scale2 (nptsPrimSt,inputObject.workingDataValues2,xObjmean/xmean)
                            ssr,inputParams.cnt = supp.sumres(npts,inputObject.workingDataValues1,inputObject.workingDataValues2,0,npts-1,skp)
                            if (ssr <= minssr):
                                minssr = ssr
                                bstprm[1][1] = parm
                        parm = parm + deltaparam #519
                        #continue #520
                    deltaparam = deltaparam/10
                    parm = bstprm [1][1] - (10 * deltaparam)
                #continue #525
                inputParams.SECRot[0] = bstprm[1][1]
                secondaryStar.workingDataValues3 = starot.starot(nptsSecSt, secondaryStar.workingDataValues2,inputParams.SECRot[0],
                                                            midLambda, deltaStepProc)#, secondaryStar.workingDataValues3)
            #end if
            
            #Rotate the tertiary
            minssr = minssr_init
            if (tst[2][1] != 0):
                # print ("    Rotating the tertiary")
                supp.scale(nptsSecSt,tertiaryStar.workingDataValues1,inputParams.TERWeight[0],tertiaryStar.workingDataValues3)
                sdch = inputParams.TERRad[0] * cnst
                tertiaryStar.workingDataValues2 = fshift.fshift(nptsPrimSt,tertiaryStar.workingDataValues3,sdch)
                parm = inputParams.TERRot[0] - 50.0
                deltaparam = 5.0
                for j in range(3): #do 525
                    for i in range(22): #do 520 i = 1,21
                        if (parm >= 0.0): #if <0 go to 519
                            tertiaryStar.workingDataValues3 = starot.starot(nptsTerSt, tertiaryStar.workingDataValues2, parm, midLambda, deltaStepProc)#, tertiaryStar.workingDataValues3)
                            inputObject.workingDataValues2 = supp.sumthm(nptsTerSt,primaryStar.workingDataValues3, secondaryStar.workingDataValues3,tertiaryStar.workingDataValues3)
                            supp.scale2 (nptsPrimSt,inputObject.workingDataValues2,xObjmean/xmean)
                            sr,ssr,inputParams.cnt = supp.sumres(nptsInputObj,inputObject.workingDataValues1, inputObject.workingDataValues2, 0, nptsPrimSt-1, skp)
                            if (ssr < minssr):
                                minssr = ssr
                                bstprm[1][1] = parm
                        parm = parm + deltaparam #519
                        #continue #520
                    deltaparam = deltaparam/10
                    parm = bstprm [2][1] - (10.0 * deltaparam)
                #continue #525
                inputParams.TERRot[0] = bstprm[1][1]
                tertiaryStar.workingDataValues3 = starot.starot(nptsTerSt, tertiaryStar.workingDataValues2, inputParams.TERRot[0], 
                                                        midLambda, deltaStepProc)#, tertiaryStar.workingDataValues3)
            #end if
            
            #####################################################                       
            #***     Find the best weights 
            #####################################################
            if ((nstar > 1) and (nwt >= 1)):
                # minssr = 100
                #Set up arrays to start
                primaryStar.workingDataValues3 = starot.starot(npts, primaryStar.workingDataValues1, inputParams.PRMRot[0],
                                                            midLambda, deltaStepProc)# former use of ndata instead of npts
                #primaryStar.workingDataValues3 = rotBroad(primaryStar.workingLambdaValues, primaryStar.workingDataValues1, 0.6, inputParams.PRMRot[0])
                pdch = inputParams.PRMRad[0] * cnst
                primaryStar.workingDataValues2 = fshift.fshift(npts,primaryStar.workingDataValues3,pdch)
                if (swt != 0.):
                    secondaryStar.workingDataValues3 = starot.starot(npts, secondaryStar.workingDataValues1, inputParams.SECRot[0],
                                                                midLambda, deltaStepProc)# former use of ndata instead of npts
                    pdch = inputParams.SECRad[0] * cnst
                    secondaryStar.workingDataValues2 = fshift.fshift(npts,secondaryStar.workingDataValues3,pdch)
                #end if
                if (twt != 0.):
                    tertiaryStar.workingDataValues3 = starot.starot(npts, tertiaryStar.workingDataValues1, inputParams.TERRot[0], 
                                                            midLambda, deltaStepProc)# former use of ndata instead of npts)
                    tdch = inputParams.TERRad[0] * cnst
                    tertiaryStar.workingDataValues2 = fshift.fshift(npts,tertiaryStar.workingDataValues3,pdch)
                #end if
                #Hold tertiary star weight fixed
                if (tst[1][2] != 0.0  and  tst[0][2] != 0.0):
                    twt = inputParams.TERWeight[0]
                    if (twt > 0.0):
                        supp.scale(npts,tertiaryStar.workingDataValues2,twt,tertiaryStar.workingDataValues3)
                    parm = 0.0
                    deltaparam = 0.10
                    swt = 1.0 - twt - parm
                    for i in range(3):                   # to 615
                        for j in range(11):              # to 610
                            if not ((parm < 0.0) or (parm > (1.0-twt))):
                                supp.scale(npts,primaryStar.workingDataValues2,parm-0.01,primaryStar.workingDataValues3)
                                supp.scale(npts,secondaryStar.workingDataValues2,swt-0.01,secondaryStar.workingDataValues3)
                                if (nstar < 3):
                                    inputObject.workingDataValues2 = supp.sumthm(npts,primaryStar.workingDataValues3,secondaryStar.workingDataValues3, zeroes)
                                else:
                                    inputObject.workingDataValues2 = supp.sumthm(npts,primaryStar.workingDataValues3,secondaryStar.workingDataValues3, tertiaryStar.workingDataValues3)
                                xmean = xPrimStarmean*parm + xSecStarmean*swt + xTerStarmean*twt
                                supp.scale2 (npts,inputObject.workingDataValues2,xObjmean/xmean)
                                ssr,inputParams.cnt = supp.sumres(npts,inputObject.workingDataValues1,
                                                                    inputObject.workingDataValues2,0,nptsPrimSt-1,skp)
                                if (ssr < minssr):
                                    minssr = ssr
                                    bstprm[0][2] = parm
                                    bstprm[1][2] = swt
                                #end if
                            parm = parm + deltaparam
                            swt = 1.0 - twt - parm
                        #610 continue 
                        deltaparam = 0.1 * deltaparam
                        parm = bstprm[0][2] - (5.0 * deltaparam)
                        swt = 1.0 - twt - parm
                    #615 continue
                    inputParams.PRMWeight[0] = bstprm[0][2]
                    inputParams.SECWeight[0] = bstprm[1][2]
                #end if
                #Hold primary weight fixed
                if ((nstar == 3) and (tst[2][2] == 0) and (tst[1][2] == 0)):
                    minssr = minssr_init
                    pwt = inputParams.PRMWeight[0]
                    supp.scale (npts,primaryStar.workingDataValues2,pwt,primaryStar.workingDataValues3)
                    parm = 0.0
                    deltaparam = 0.1
                    twt = 1.0 - pwt - parm
                    for i in range(3):                    # to 625
                        for j in range(11):               # to 620
                            if not ((parm < 0.0) or (parm > (1.0-pwt))):        #if not not go to 619
                                supp.scale(npts,secondaryStar.workingDataValues2,parm,secondaryStar.workingDataValues3)
                                supp.scale(npts,tertiaryStar.workingDataValues2,twt,tertiaryStar.workingDataValues3)
                                inputObject.workingDataValues2 = supp.sumthm(npts,primaryStar.workingDataValues3,
                                                                    secondaryStar.workingDataValues3,tertiaryStar.workingDataValues3)
                                xmean = xPrimStarmean*pwt + xSecStarmean*parm + xTerStarmean*twt
                                supp.scale2 (npts,inputObject.workingDataValues2,xObjmean/xmean)
                                ssr,sr,inputParams.cnt = supp.sumres(npts,inputObject.workingDataValues2, inputObject.workingDataValues2 ,
                                                            0,nptsPrimSt-1,skp)
                                if (ssr <= minssr):
                                    minssr = ssr
                                    bstprm[1][2] = parm
                                    bstprm[2][2] = twt
                            #end if
                            parm = parm + deltaparam    #  619
                            twt = 1.0 - pwt - parm
                        #continue   #620
                        deltaparam  = 0.1 * deltaparam
                        parm = bstprm[1][2] - (5.0 * deltaparam)
                        twt = 1.0 - pwt - parm
                    #continue   #625          
                    inputParams.SECWeight[0] = bstprm[1][2]
                    inputParams.TERWeight[0] = bstprm[2][2]
                #end if
                # Hold secondary weight fixed
                if ((nstar == 3) and (tst[0][1] != 0.) and (tst[2][2] != 0.)):
                    minssr = minssr_init
                    swt = inputParams[1][2]
                    supp.scale (npts,secondaryStar.workingDataValues2,swt,secondaryStar.workingDataValues3)
                    parm = 0.0
                    deltaparam = 0.1
                    twt = 1.0 - swt - parm
                    for i in range(3):
                        for j in range(11):
                            if not ((parm < 0.0) or (parm > (1.0-swt))): #if not not go to 629
                                supp.scale(npts,primaryStar.workingDataValues2,parm,primaryStar.workingDataValues3)
                                supp.scale(npts,tertiaryStar.workingDataValues2,twt,tertiaryStar.workingDataValues3)
                                inputObject.workingDataValues2 = supp.sumthm(npts,primaryStar.workingDataValues3,
                                                                    secondaryStar.workingDataValues3,tertiaryStar.workingDataValues3)
                                xmean = xPrimStarmean *parm + xSecStarmean*swt + xTerStarmean*twt
                                supp.scale2 (npts,inputObject.workingDataValues2,xObjmean/xmean)
                                ssr,sr,inputParams.cnt = supp.sumres(npts,inputObject.workingDataValues1,inputObject.workingDataValues2,
                                                                0,nptsPrimSt-1,skp)
                                if (ssr <= minssr):
                                    minssr = ssr
                                    bstprm[0][2] = parm
                                    bstprm[2][2] = twt
                                #end if
                            parm = parm + deltaparam #  629
                            twt = 1.0 - swt - parm
                        #end for
                        deltaparam = 0.1 * deltaparam
                        parm = bstprm[0][2] - (5.0 * deltaparam)
                        twt = 1.0 - swt - parm
                    #end for
                    inputParams.PRMWeight[0] = bstprm[0][2]
                    inputParams.TERWeight[0] = bstprm[2][2] 
                #end if
            #end if
        
            ######################################################
            #             Write iteration results 
            #             to the standard output
            #####################################################
            V_rot = int(inputParams.PRMRot[0]*1000)/1000       
            # print
            # print (' Results of iteration #' + str(niter+1))
            print (' Minimum RMS of the residuals = ' + str(math.sqrt (math.fabs(minssr) / inputParams.cnt)))
            print ('-----------------------------------------------')
            print ('             PRI      SEC      TER')
            print (' V_rad   = ' + str('{:6.3f}'.format(inputParams.PRMRad[0]))    + "    " + str('{:6.3f}'.format(inputParams.SECRad[0]))    + "    " + str('{:6.3f}'.format(inputParams.TERRad[0])) )
            print (' V_rot   = ' + str('{:6.3f}'.format(inputParams.PRMRot[0]))    + "    " + str('{:6.3f}'.format(inputParams.SECRot[0]))    + "    " + str('{:6.3f}'.format(inputParams.TERRot[0])) )
            print (' Weight  = ' + str('{:6.3f}'.format(inputParams.PRMWeight[0])) + "    " + str('{:6.3f}'.format(inputParams.SECWeight[0])) + "    " + str('{:6.3f}'.format(inputParams.TERWeight[0])) )
            print ('-----------------------------------------------')
            # print ('')
            # print ('')
            ################################################################
            #                   End iteration loop
            ################################################################
    
        ################################################################
        #                 Calculate final model                         
        ################################################################
        supp.scale(nptsPrimSt,primaryStar.workingDataValues1,inputParams.PRMWeight[0],primaryStar.workingDataValues3)
        primaryStar.workingDataValues2 = starot.starot(nptsPrimSt, primaryStar.workingDataValues3, inputParams.PRMRot[0], 
                                                    midLambda, deltaStepProc)# , primaryStar.workingDataValues2)
        #primaryStar.workingDataValues2 = rotBroad(primaryStar.workingLambdaValues, primaryStar.workingDataValues3, 0.6, inputParams.PRMRot[0])
        pdch = inputParams.PRMRad[0] * cnst
        primaryStar.workingDataValues3 = fshift.fshift(nptsPrimSt,primaryStar.workingDataValues2,pdch)
        
        if (swt != 0.0):
            supp.scale(nptsSecSt,secondaryStar.workingDataValues1,inputParams.SECWeight[0],secondaryStar.workingDataValues3)
            secondaryStar.workingDataValues2 = starot.starot(nptsSecSt, secondaryStar.workingDataValues3, inputParams.SECRot[0], midLambda, deltaStepProc)# , secondaryStar.workingDataValues2)
            sdch = inputParams.SECRad[0] * cnst
            secondaryStar.workingDataValues3 = fshift.fshift(nptsSecSt,secondaryStar.workingDataValues2,sdch)
        
        if (twt != 0.0):
            supp.scale(npts,tertiaryStar.workingDataValues1,inputParams.TERWeight[0],tertiaryStar.workingDataValues3)
            tertiaryStar.workingDataValues2 = starot.starot(nptsSecSt, tertiaryStar.workingDataValues3, inputParams.TERRot[0], 
                                                        midLambda, deltaStepProc)# , tertiaryStar.workingDataValues2)
            tdch = inputParams.TERRad[0] * cnst
            tertiaryStar.workingDataValues3 = fshift.fshift(nptsSecSt,tertiaryStar.workingDataValues2,tdch)
            
        if nstar == 2:
            inputObject.workingDataValues2 = supp.sumthm(nptsPrimSt,primaryStar.workingDataValues3,secondaryStar.workingDataValues3,zeroes)
        elif nstar == 1:
            inputObject.workingDataValues2 = cpy.copy(primaryStar.workingDataValues3) #sumthm(nptsPrimSt,primaryStar.workingDataValues3,zeroes,zeroes)        
        else:
            inputObject.workingDataValues2 = supp.sumthm(nptsPrimSt,primaryStar.workingDataValues3,
                                                secondaryStar.workingDataValues3,tertiaryStar.workingDataValues3)
                                                
        xmean = xPrimStarmean*inputParams.PRMWeight[0] + xSecStarmean*inputParams.SECWeight[0]+ xTerStarmean*inputParams.TERWeight[0]
        supp.scale2 (nptsInputObj,inputObject.workingDataValues2,xObjmean/xmean)
        inputObject.nanSanityCheck(2)
        #####################################################################
        #            Calculate the sum of the
        #             squared residuals (ssr)
        #####################################################################
            
        ssr,inputParams.cnt = supp.sumres(npts,inputObject.workingDataValues1,inputObject.workingDataValues2,0,nptsInputObj-1,skp)
        
        #####################################################################
        #     Write output specs if desired
        #####################################################################
        wp = {"NAXIS":1,"NAXIS1": nptsInputObj,"CRVAL1":str(inputObject.initialLambda), "CDELT1":str(inputObject.deltaStep), "CRPIX1":1, "CTYPE1":"LINEAR"}
        outputFilename = inputObject.object + "@" + str(inputObject.date) +"_" + inputParams.linespcf
        # if inputParams.synthSpecName != None:
        #     WFS.write1dFitsSpec(inputParams.synthSpecName + "_" + outputFilename + '.fits' , inputObject.workingDataValues2, waveParams=wp, clobber=True)
        #      inputObject.writeToDat(inputParams.synthSpecName + "_" + outputFilename,2)
        ####################################################################
        #         Substract model from object
        ####################################################################
        for i in range(10,nptsInputObj-10):
            #Substracted Spectrum = Initial Read Spectrum - Synthetic Spectrum
            inputObject.workingDataValues3[i] = (inputObject.workingDataValues1[i] - inputObject.workingDataValues2[i])
        #print inputObject.workingDataValues3[i], " = " , inputObject.workingDataValues1[i], " - ", inputObject.workingDataValues2[i]
        inputObject.nanSanityCheck(3)        
    
        ###################################################################
        #         Write subtracted spectrum
        ###################################################################
        if inputParams.substractedSpecName != None:
            #  WFS.write1dFitsSpec(inputParams.substractedSpecName + "_" + outputFilename + '.fits', inputObject.workingDataValues3, 
            #                      waveParams=wp, clobber=True)
             inputObject.writeToDat(inputParams.substractedSpecName + "_" + outputFilename, 3)
        ###################################################################
        #                That's all folks!
        ###################################################################
        print (' ')
        
        ##------------------------RV Out File-------------------------------------------##
        if nstar == 1:
            if (rv_values_available == False):
                rv_outfile = open(filenameRV1, 'a')
                rv_outfile.write(str(inputObject.date) + " " + str(inputParams.PRMRad[0]) + '\n')
                print ("#######################################################    outfilenameRV:" , filenameRV1)
                rv_outfile.close()
        ##------------------------------------------------------------------------------##

        ###############################################################################
        ########- Calculate the EW of the line-########################################
        if(inputParams.linespcf != 'NONE'):
            lambdas = tools.lambdaData("lambdas.dat")
            
            
            ldoLineInput = lambdas.setLineInputbyKey(inputParams.linespcf)
            print(inputParams.linespcf)
            
            if(nstar == 1):
                if inputObject.RVfromFITS != 0.0: ## The CARMENES provided spectrum contains these parameters in the fits headings
                    inputBaryCorr = inputObject.BaryCorr
                    inputValue = inputObject.RVfromFITS + inputBaryCorr
                    print ("inputValue = ",inputValue)
                else:                            ## General spectrum. No RV or barycentric correction available in fits headings
                    inputBaryCorr = inputObject.BaryCorr = 0
                    inputValue = inputParams.PRMRad[0]
                    print ("inputValue (no RVfromFITS) = ",inputValue)
            else:
                inputValue = inputParams.PRMRad[0]
                inputValue2 = inputParams.SECRad[0]

            if nstar == 1:
                ldoObs, _ = lambdas.calculateDopplerDisplacement(inputValue)
                print (ldoObs)
                ldoObs, param2, param5, param6, param3, param4 = lambdas.getMaximumValue(inputObject, inputValue, inputParams.FindMax_Tolerance) # ldoObs is modified within getMaximumValues
                ## FOR DEBUGGING PURPOSES. To draw the boundaries of the zone where to search for the maximum.
                if (debugging == True):
                    subPlot1.plot(ldoObs, param2, 'b+')
                    subPlot1.axvline(x=param5, c = 'r', linestyle='dashed')
                    subPlot1.axvline(x=param6, c = 'r', linestyle='dashed')
            elif (nstar == 2):
                ldoObs,ldoObs2 = lambdas.calculateDopplerDisplacement(inputValue,True, inputValue2)
                print (lambdas.ldoObs, lambdas.ldoObs2)
                if lambdas.ldoObs > lambdas.ldoObs2: #preparing for getMaximumValue
                    lambdas.ldoObs, lambdas.ldoObs2 = lambdas.ldoObs2 , lambdas.ldoObs
                print (lambdas.ldoObs, lambdas.ldoObs2)
                print (ldoObs, ldoObs2)
                ldoObs, param2, param5, param6, ldoObs2 , param4 = lambdas.getMaximumValue(inputObject, lambdas.ldoObs, inputParams.FindMax_Tolerance, 2, lambdas.ldoObs2) # ldoObs is modified within getMaximumValues
                if (debugging == True):
                    subPlot1.axvline(x=param5, c = 'r', linestyle='dashed')
                    subPlot1.axvline(x=param6, c = 'r', linestyle='dashed')
                    subPlot1.axvspan(param5, param6, color='green', alpha=0.1)

            if (inputObject.RVfromFITS != 0.0):
                print("ldoObs = ", ldoObs, " !!! ")
            else:
                print(ldoObs, param2)
                vRad = lambdas.calculateRadialVelocityfromDopplerDisplacement()
                ldoLineInput = lambdas.setLineInputbyValue(ldoObs)
            if nstar == 1:
                ldo =  ldoObs    
            elif nstar == 2:
                ldo =  (ldoObs2+ ldoObs)/2
                # ldo = (ldoObs + ldoObs2)/2
            
            ldoLineStart = lambdas.setLineStart(ldo - ldoLineInput[1])
            ldoLineEnd   = lambdas.setLineEnd  (ldo + ldoLineInput[1])
            if (inputObject.RVfromFITS == 0.0):
                print ("ldoLineInput= ", ldoLineInput[0], "\nldoLineStart= ", ldoLineStart, "ldoLineEnd= ", ldoLineEnd)

            print ("deltaStep = ", deltaStepProc)
            ## The following allows to select the subset of points of the subtracted spectrum on which we calculate the EW
            lambdas.selectLambdaRanges(inputObject.workingLambdaValues, deltaStepProc)
            if nstar == 1:
                EWpopt = lambdas.CalculateEW(inputObject)
            elif nstar == 2:
                EWpopt = lambdas.CalculateEW(inputObject,None, 
                        initial_weight=inputParams.PRMWeight[0], initial_weight2=inputParams.SECWeight[0],nstar=2)
            
            EW1, EW2, star1_spectrum, star2_spectrum, flambdavalues, summ_up, abs_error, rel_error = EWpopt 
            FWHM = lambdas.getFWHM(inputObject)
            SNR = max(inputObject.SNR, primaryStar.SNR)
            errorEW = lambdas.calculateErrorCayrelFormulae(FWHM,deltaStepProc,SNR)
        
            print("#################################################")
            print("     FWHM = ",FWHM)
            if nstar == 1:
                print("     EW   = " + str(EW1) + " +- " + str(errorEW))
            elif nstar == 2:
                print("     EW1   = " + str(EW1) + r" \\pm " + str(abs_error) + ' (' + str('{:3.1f}'.format(rel_error)) + '%)')
                print("     EW2   = " + str(EW2) + r" \\pm " + str(abs_error) + ' (' + str('{:3.1f}'.format(rel_error))+ '%)')
            print("#################################################")
            
            ##################### Addition 10/11/2020
            if(inputParams.deltaStep == -1.0):
                print("DELTA_STEP = ", deltaStepAvg) 
                print("CONST_ = ", cnst, '\n\n')
            #########################################

            # found = len(inputParams.workpath) - inputParams.workpath.find('\\')
            if (debugging == True): print(inputParams.workpath)
            if ("nir-tac" not in inputParams.workpath) and ("NIR-TAC" not in inputParams.workpath) :
                if not os.access(inputParams.workpath[:-4] + "\\RES", os.F_OK): os.mkdir(inputParams.workpath[:-4] + "\\RES")
                ewOutfile = open(inputParams.workpath[:-4] + "\\RES\\" + inputObject.object + "_" + inputParams.linespcf, 'a')
                print (inputParams.workpath[:-4] + "\\RES\\" + inputObject.object + "_" + inputParams.linespcf)
            else:
                if not os.access(inputParams.workpath[:-8] + "\\RES", os.F_OK): os.mkdir(inputParams.workpath[:-8] + "\\RES")
                ewOutfile = open(inputParams.workpath[:-8] + "\\RES\\" + inputObject.object + "_" + inputParams.linespcf, 'a')
            if nstar == 1:
                ewOutfile.write(str(inputObject.date)+ "  " + str(EW1) + "  " + str(errorEW) + "  " + str(inputParams.PRMRad[0]) + "  " + str(inputObject.BaryCorr) + "  " + str(inputParams.PRMRot[0]) + '\n')
            elif nstar == 2:
                ewOutfile.write(str(inputObject.date)+ "  " + str(EW1) + "  " + str(errorEW) + "  " + str(inputParams.PRMRad[0]) + "  " + str(inputObject.BaryCorr) + "  " + str(inputParams.PRMRot[0]) + str(inputParams.SECRot[0]) + str(inputParams.SECRad[0]) + '\n')
            else: ## nstar == 3
                ewOutfile.write(str(inputObject.date)+ "  " + str(EW1) + "  " + str(errorEW) + "  "+ str(inputParams.PRMRad[0]) + "  " + str(inputObject.BaryCorr) + "  " + str(inputParams.PRMRot[0]) + str(inputParams.SECRot[0]) + str(inputParams.SECRad[0]) + str(inputParams.TERRot[0]) + str(inputParams.TERRad[0]) + '\n')
            ###############################################################################
            ########## To draw the integration limits
            subPlot1.axvline(x=inputObject.workingLambdaValues[lambdas.start], linestyle='dashed')
            subPlot1.axvline(x=inputObject.workingLambdaValues[lambdas.end],   linestyle='dashed')
            if debugging == True:
                subPlot1.axvspan(inputObject.workingLambdaValues[lambdas.start], 
                                 inputObject.workingLambdaValues[lambdas.end], color='blue', alpha=0.1)
            subPlot1.axhline(y=0.0, linestyle='dashed')
            subPlot2.axvline(x=inputObject.workingLambdaValues[lambdas.start], linestyle='dashed')
            subPlot2.axvline(x=inputObject.workingLambdaValues[lambdas.end],   linestyle='dashed')

            if nstar ==2:
                    subPlot1.plot(flambdavalues, star1_spectrum, 'r-', lw = 1.0)
                    subPlot1.plot(flambdavalues, star2_spectrum, 'k-', lw = 1.0)
                    if (debugging == True):## FOR DEBUGGING PURPOSES
                        subPlot1.plot(ldoObs2, param4, 'r*') 
                        subPlot1.plot(ldoObs, param2, 'b+')
                        if summ_up != None:
                            subPlot1.plot(flambdavalues, summ_up, c = 'b', linestyle = 'dashed', lw = 1.25)
        
        if len(inputObject.workingLambdaValues) == len(inputObject.workingDataValues2):
                subPlot2.plot(inputObject.workingLambdaValues, inputObject.workingDataValues1, 'b-', lw = 1.0)
                subPlot2.plot(inputObject.workingLambdaValues, inputObject.workingDataValues2, 'r-', lw = 1.0)
                subPlot1.plot(inputObject.workingLambdaValues, inputObject.workingDataValues3, 'g-', lw = 1.0)
        elif len(inputObject.workingLambdaValues) >= len(inputObject.workingDataValues2):
            subPlot2.plot(inputObject.workingLambdaValues, inputObject.workingDataValues1, 'b-', lw = 1.0)
            subPlot2.plot(inputObject.workingLambdaValues[0:len(inputObject.workingDataValues2)],
                            inputObject.workingDataValues2, 'r-', lw = 1.0)
            subPlot1.plot(inputObject.workingLambdaValues, inputObject.workingDataValues3, 'g-', lw = 1.0)
        elif len(inputObject.workingLambdaValues) <= len(inputObject.workingDataValues2):
            subPlot2.plot(inputObject.workingLambdaValues, inputObject.workingDataValues1, 'b-', lw = 1.0)
            subPlot2.plot(inputObject.workingLambdaValues,
                            inputObject.workingDataValues2[0:len(inputObject.workingLambdaValues)], 'r-', lw = 1.0)
            subPlot1.plot(inputObject.workingLambdaValues, inputObject.workingDataValues3, 'g-', lw = 1.0)

        f.subplots_adjust(hspace = 0.3)

        dateTitle = int(inputObject.date*1000000)/1000000
        
        strTitle = str(inputObject.object) + "@" + str(dateTitle) + "_" + inputParams.linespcf 
        
        # # CaIIH&K
        if(inputParams.linespcf == "CaIIH" or inputParams.linespcf == "CaIIK"):
            subPlot1.set_title("Subtracted Spectrum for "+ inputParams.linespcf + r" / $v_{\rm rot}sini$ = " + str(V_rot) + " km/s")
            subPlot1.set_ylim(-1.0, inputParams.MaxFluxDisplayed_Sub)
            subPlot2.set_title(inputObject.object + "@" + str(int(dateTitle)) +"  &  Synthetic Spectra for " + inputParams.linespcf)
            subPlot2.set_ylim(-0.5, inputParams.MaxFluxDisplayed_Obj)
            # subPlot1.set_xlim(3945,3985)
            # subPlot2.set_xlim(3945,3985)
        # # Hbeta
        if(inputParams.linespcf == "Hbeta"):
            subPlot1.set_title("Subtracted Spectrum for "+ "$H\beta$" + r" / $v_{\rm rot}sini$ = " + str(V_rot) + " km/s")
            subPlot1.set_ylim(-1.0, inputParams.MaxFluxDisplayed_Sub)
            subPlot2.set_title(inputObject.object + "@" + str(int(dateTitle)) +"  &  Synthetic Spectra for " + r"$H\beta$")
            subPlot2.set_ylim(-0.5, inputParams.MaxFluxDisplayed_Obj)
            subPlot1.set_xlim(4847.5,4872.5)
            subPlot2.set_xlim(4847.5,4872.5)
            # print ("Here in Hbeta!!")
        # # HeD3 NaD2/NaD1
        if(inputParams.linespcf == 'HeD3' or inputParams.linespcf == 'NaD2' or inputParams.linespcf == 'NaD1'):
            subPlot1.set_title(r"Subtracted Spectrum for "+ inputParams.linespcf + r" / $v_{\rm rot}sini$ = " + str(V_rot) + " km/s")
            subPlot1.set_ylim(-0.25, inputParams.MaxFluxDisplayed_Sub)
            subPlot2.set_title(inputObject.object + "@" + str(float('{:6.5f}'.format(dateTitle))) +"  &  Synthetic Spectra for " + inputParams.linespcf)
            subPlot2.set_ylim(-0.2, inputParams.MaxFluxDisplayed_Obj)
            subPlot1.set_xlim(5860,5910)
            subPlot2.set_xlim(5860,5910)
        # # Halpha
        elif(inputParams.linespcf == 'Halpha'):
            subPlot1.set_title(r"Subtracted Spectrum for "+ r"$H\alpha$" + r" / $v_{\rm rot}sini$ = " + str('{:7.4f}'.format(inputParams.PRMRot[0])) + " km/s")
            subPlot2.set_title(inputObject.object + "@" + str(int(dateTitle)) +"  &  Synthetic Spectra for " + r"$H\alpha$")
            subPlot2.set_ylim(0.0, inputParams.MaxFluxDisplayed_Obj)
            subPlot1.set_ylim(-1.0, inputParams.MaxFluxDisplayed_Sub)
            subPlot1.set_xlim(6515,6605)
            subPlot2.set_xlim(6515,6605)
        #  CaIRT-a
        elif(inputParams.linespcf == 'CaIRT-a'):
            subPlot1.set_title(r"Subtracted Spectrum for "+ inputParams.linespcf + r" / $v_{\rm rot}sini$ = " + str(V_rot) + " km/s")
            subPlot2.set_title(inputObject.object + "@" + str(int(dateTitle)) +"  &  Synthetic Spectra for " + inputParams.linespcf)
            subPlot1.set_ylim(-0.5, inputParams.MaxFluxDisplayed_Sub)
            subPlot2.set_ylim(0.,inputParams.MaxFluxDisplayed_Obj)
            subPlot1.set_xlim(8460,8530)
            subPlot2.set_xlim(8460,8530)
        # # #################### Paschen D
        elif(inputParams.linespcf == 'PaschenD'):
            subPlot2.set_title(inputObject.object + "@" + str(int(dateTitle)) +"  &  Synthetic Spectra for " + inputParams.linespcf)
            subPlot2.set_ylim(0.25, inputParams.MaxFluxDisplayed_Obj)
            subPlot1.set_title(r"Subtracted Spectrum for "+ inputParams.linespcf + r" / $v_{\rm rot}sini$ = " + str(V_rot) + " km/s")
            subPlot1.set_ylim(-0.25, inputParams.MaxFluxDisplayed_Sub)
            subPlot1.set_xlim(10040,10065)
            subPlot2.set_xlim(10040,10065)
        #  ##################### PaschenB
        elif(inputParams.linespcf == 'PaschenB'):
            subPlot2.set_title(inputObject.object + "@" + str(int(dateTitle)) +"  &  Synthetic Spectra for " + inputParams.linespcf)
            subPlot2.set_ylim(0.5, inputParams.MaxFluxDisplayed_Obj)
            subPlot1.set_title(r"Subtracted Spectrum for "+ inputParams.linespcf + r" / $v_{\rm rot}sini$ = " + str(V_rot) + " km/s")
            subPlot1.set_ylim(-0.5, inputParams.MaxFluxDisplayed_Sub)
            subPlot1.set_xlim(12800, 12840)
            subPlot2.set_xlim(12800, 12840)
        else:
        # Any
            subPlot2.set_title(inputObject.object + "@" + str(int(dateTitle)) +"  &  Synthetic Spectra for " + inputParams.linespcf)
            subPlot2.set_ylim(0.0, inputParams.MaxFluxDisplayed_Obj)
            subPlot1.set_title(r"Subtracted Spectrum for "+ inputParams.linespcf + r" / $v_{\rm rot}sini$ = " + str(V_rot) + " km/s")
            subPlot1.set_ylim(-0.5, inputParams.MaxFluxDisplayed_Sub)

        plt.xlabel(r"$\lambda$  [${\rm \AA}$]")
        plt.ylabel("Normalised Flux")
        # plt.legend()
        
        if ("nir-tac" not in inputParams.workpath) and ("NIR-TAC" not in inputParams.workpath) :
            figureName = inputParams.workpath[:-4]  + "\\RES\\" + inputParams.linespcf + "\\" + strTitle
            if nstar == 1:
                figureName = figureName + ".png"
            elif nstar == 2:
                figureName = figureName + "_asSB2.png"
            if not os.access(inputParams.workpath[:-4]  + "\\RES\\" + inputParams.linespcf, os.F_OK):
                os.mkdir(inputParams.workpath[:-4] + "\\RES\\" + inputParams.linespcf)
        else:
            figureName = inputParams.workpath[:-8]  + "\\RES\\" + inputParams.linespcf + "\\" + strTitle
            if nstar == 1:  figureName = figureName + ".png"
            elif nstar == 2:    figureName = figureName + "_asSB2.png"
            if not os.access(inputParams.workpath[:-8] + "\\RES\\" + inputParams.linespcf, os.F_OK):
                os.mkdir(inputParams.workpath[:-8] + "\\RES\\" + inputParams.linespcf)        

        plt.savefig(figureName, dpi = 300)
        if (debugging == True): plt.show()
        
        
        if (".fits" in str(objectInputFITSfromList)):   inputObject.fits_hdu_in.close()

        primaryStar.fits_hdu_in.close()
        if (nstar == 2):    secondaryStar.fits_hdu_in.close()
        if (nstar == 3):    tertiaryStar.fits_hdu_in.close()
        
        if (inputParams.linespcf != 'NONE'):     ewOutfile.close()
        plt.close()
        inputParams = tools.readConfigFile(spcfFile)
    ## NEXT SPECTRUM FILE
    
## PROGRAM END 
