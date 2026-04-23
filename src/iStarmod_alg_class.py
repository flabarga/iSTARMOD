from iStarmod_fits import * 
import matplotlib.pyplot as plt
import os
from pathlib import Path
from queue import Queue, Empty
from threading import Timer
from time import sleep

class iStarmodGUI():

    
    def starmod(self, spcfFile, f, subPlot1, subPlot2, frm, debugging = False):

        import iStarmod_tools as tools
        import iStarmod_rvvalues as rv
        import support as supp
        import starot as starot
        import fshift as fshift
        import copy as cpy
        import math
        
        
        ##################################################################################3
        #READING THE SPECIFICATION FILE
        ##################################################################################3
        try:
            inputParams = tools.readConfigFile(spcfFile)
        except FileNotFoundError:
            print("INPUT PARAMETERS FILE NOT FOUND. CHECK FILENAME OR PATH")
            return
        except PermissionError as err:
            print(err)
            return
            
        #ensure that the weights are correctly specified
        inputParams.checkWeights()
        tst = []
        #############################
        #Implementing the specification var/fix.
        #############################
        #If the parameter is 'fix', then tst[n] must be zero
        #When parameter is 'var', is assigned to tst[n] the input parameter value
        print()
        tst.append([])
        if inputParams.PRMRad[1]:     tst[0].append(0)
        else:                         tst[0].append(inputParams.PRMRad[0])
        if inputParams.PRMRot[1]:     tst[0].append(0)
        else:                         tst[0].append(inputParams.PRMRot[0])
        if inputParams.PRMWeight[1]:  tst[0].append(0)
        else:                         tst[0].append(inputParams.PRMWeight[0])
        if debugging:   print ("tst[0]: " + str(tst[0]))
        tst.append([])
        if inputParams.SECRad[1]:     tst[1].append(0)
        else:                         tst[1].append(inputParams.SECRad[0])
        if inputParams.SECRot[1]:     tst[1].append(0)
        else:                         tst[1].append(inputParams.SECRot[0])
        if inputParams.SECWeight[1]:  tst[1].append(0)
        else:                         tst[1].append(inputParams.SECWeight[0])
        if debugging:   print ("tst[1]: " + str(tst[1]))
        tst.append([])
        if inputParams.TERRad[1]:     tst[2].append(0)
        else:                         tst[2].append(inputParams.TERRad[0])
        if inputParams.TERRot[1]:     tst[2].append(0)
        else:                         tst[2].append(inputParams.TERRot[0])
        if inputParams.TERWeight[1]:  tst[2].append(0)
        else:                         tst[2].append(inputParams.TERWeight[0])
        if debugging:   print ("tst[2]: " + str(tst[2]))
        nstar = 3
        if (inputParams.tertiaryStarfits != None) or ('NONE' in inputParams.secondaryStarfits0):
            nstar = 2
            tst[2][0] = 0
            tst[2][1] = 0
            tst[2][2] = 0
            inputParams.TERWeight[0] = 0.0
            inputParams.SECWeight[0] = float('{:5.3f}'.format(1.0 - inputParams.PRMWeight[0])) #To avoid floating point rounding errors
            
        if ('NONE' in inputParams.secondaryStarfits0):
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
        if (inputParams != None and inputParams.use_rv_values):
            if (nstar == 1):
                rv_values_available = False
                rv_values = rv.RVValues()
                if not os.access(Path(inputParams.resultspath), os.F_OK): 
                    os.mkdir(Path(inputParams.resultspath))
                if not os.access (Path(inputParams.resultspath) / Path('RES'), os.F_OK): 
                    os.mkdir(Path(inputParams.resultspath) / Path("RES"))

                if Path(inputParams.resultspath) == Path("\\"): inputParams.resultspath = "."
                filenameRV1 = inputParams.resultspath / Path("RES") / Path("rvvalues.dat")
                filenameRV2 = inputParams.resultspath / Path("RES")

                if os.access(filenameRV1, os.F_OK):
                    rv_values.readFile(filenameRV1)
                    rv_values_available = True
                    if (debugging == True): print ("filenameRV:" , filenameRV1)
                else:
                    if not os.access(filenameRV2, os.F_OK):
                        os.mkdir(filenameRV2) 
                    if (debugging == True): 
                        print ("outfilenameRV:" , filenameRV2)
        else:
            rv_values_available = False
            
        ############################################################################################################
        
        
        ##################################################################################3
        #READING and plotting THE OBJECT FITS FILE
        ##################################################################################3
        dirPath = inputParams.workpath        
        # print("Path =", dirPath)
        # print("Exists?", dirPath.exists())
        # print("Is dir?", dirPath.is_dir())
        print("Path = ", dirPath)
        listFITSFiles = list(dirPath.glob("car-*.fits"))
        # To cope with specifying not CARMENES files
        if len(listFITSFiles) == 0 and '*' not in inputParams.objectInputFits0:
            listFITSFiles.append(inputParams.objectInputFits)
        elif len(listFITSFiles) == 0 :
            listFITSFiles = list(dirPath.glob("*.fits"))
        
        if len(listFITSFiles) == 0:
            listFITSFiles = list(dirPath.glob("*.dat"))

        print("list FITS files: ", listFITSFiles)

        try:
            for objectInputFITSfromList in listFITSFiles:

                if frm.stop_istarmod_event.is_set():
                    print("Here5!!!. Stopped execution")
                    return
               
                ##################  Newly added lines 18/04/2025  ###############################################
                inputParams.checkWeights()
                nwt = tst[0][2] + tst[1][2] + tst[2][2]
                pwt = inputParams.PRMWeight[0] 
                swt = inputParams.SECWeight[0] 
                twt = inputParams.TERWeight[0] 
                ##################################################################################################
                
                # f, (subPlot1,subPlot2) = plt.subplots(2, figsize = (10,7), dpi = 150, sharex = True)#, sharey = True)
                # subPlot1.clear()
                # subPlot2.clear()    
                
                objectInputFITS = objectInputFITSfromList
                if debugging:  print("objectinputfits : ", objectInputFITS)
                inputObject = tools.readFitsFile(objectInputFITS, inputParams.spectraFormatAperture)
                # print("MJD-OBS = ", inputObject.date)
                #Reading and preprocessing Input file
                #evaluate if 'pixel range' is given in pixels or in Angstroms
                #it is assumed that the start and end pixels specified belong to the same order 
                xlo,xhi, skp = inputObject.selectPixelRange(inputParams, inputParams.spectraFormatAperture, debugging)#, subPlot2, False)
                inputObject.nanSanityCheck()
                if debugging: print("SNR = ", inputObject.SNR)
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

                xObjmean = supp.gtmean(nptsInputObj,inputObject.workingDataValues1,0,nptsInputObj-1,skp)
                midLambda = inputObject.initialLambda + inputParams.deltaStep * len(inputObject.workingLambdaValues) / 2.0
                if (debugging == True): print(nptsInputObj, inputObject.workingLambdaValues[len(inputObject.workingLambdaValues)-1], midLambda)
                cnst = (midLambda/299792.458)/inputParams.deltaStep     #Conversion Factor from Doppler Effect Formula

                #Reading Prymary (Reference) Star spectrum
                if (inputParams.spectraFormatAperture != -1 and inputParams.spectraFormatAperturePrimary == -1):
                    primaryStar = tools.readFitsFile(inputParams.primaryStarfits, inputParams.spectraFormatAperture)
                    if primaryStar != None:
                        primaryStar.selectPixelRange(inputParams, inputParams.spectraFormatAperture)#, subPlot2, False)
                elif inputParams.spectraFormatAperturePrimary != -1:
                    primaryStar = tools.readFitsFile(inputParams.primaryStarfits, inputParams.spectraFormatAperturePrimary)
                    if primaryStar != None:
                        primaryStar.selectPixelRange(inputParams, inputParams.spectraFormatAperturePrimary)#, subPlot2, False)
                if primaryStar != None:     
                    primaryStar.nanSanityCheck()
                    if debugging:   print("SNR = ", primaryStar.SNR)
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
                        if inputParams.deltaStep != -1.0: #if inputObject bCarmenes != True then inputParams.deltaStep has been assigned but inpuParams.nrpoints == 0... 
                            inputParams.deltaStep = -1.0
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


                if 'NONE' not in inputParams.secondaryStarfits.name:
                    if (inputParams.spectraFormatApertureSecondary != -1):
                        secondaryStar = tools.readFitsFile(inputParams.secondaryStarfits, inputParams.spectraFormatApertureSecondary)
                        if secondaryStar != None:
                            secondaryStar.selectPixelRange(inputParams, inputParams.spectraFormatApertureSecondary)#, subPlot2, False)
                    elif (inputParams.spectraFormatApertureSecondary == -1 and inputParams.spectraFormatAperture != -1):
                        secondaryStar = tools.readFitsFile(inputParams.secondaryStarfits, inputParams.spectraFormatAperture)
                        if secondaryStar != None:
                            secondaryStar.selectPixelRange(inputParams, inputParams.spectraFormatAperture) #, subPlot2, False)
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
                    
                if 'NONE' not in inputParams.tertiaryStarfits.name:
                    if (inputParams.spectraFormatApertureTertiary != -1):
                        tertiaryStar = tools.readFitsFile(inputParams.tertiaryStarfits, inputParams.spectraFormatApertureTertiary)
                        if tertiaryStar != None:
                            tertiaryStar.selectPixelRange(inputParams, inputParams.spectraFormatApertureTertiary)#, subPlot2, False)
                    elif (inputParams.spectraFormatAperture != -1 and inputParams.spectraFormatApertureTertiary == -1):
                        tertiaryStar = tools.readFitsFile(inputParams.tertiaryStarfits, inputParams.spectraFormatAperture)
                        if tertiaryStar != None:
                            tertiaryStar.selectPixelRange(inputParams, inputParams.spectraFormatAperture)#, subPlot2, False)
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
                print (' Object    = ', objectInputFITS)
                print (' Primary   = ', inputParams.primaryStarfits)
                print (' Secondary = ', inputParams.secondaryStarfits)
                print (' Tertiary  = ', inputParams.tertiaryStarfits)
                print (' Order     = ' + str(inputParams.spectraFormatAperture))
                # print ("")
                print ('####################################################')
                print ("")
                print ("")
                print ("-----------------------------------------------")
                
                for niter in range(inputParams.numIter):
                    if frm.stop_istarmod_event.is_set():
                        print("Here3!!!. Stopped execution")
                        return
                    ######################################################################################
                    #                Set-up of the initial guess                                          
                    ######################################################################################
                    print ("ITERATION #" + str(niter+1))
                    #print primaryStar.workingDataValues3
                    supp.scale(nptsPrimSt,primaryStar.workingDataValues1,pwt,primaryStar.workingDataValues3)
                    
                    ###############################################Rotationally broaden it
                    # (workingDataValues2 holds scaled & broad prim)
                    primaryStar.workingDataValues2 = starot.starot(nptsPrimSt, primaryStar.workingDataValues3, inputParams.PRMRot[0],
                                                                    midLambda, deltaStepProc)# ,primaryStar.workingDataValues2)
                    ###############################################Shift it in radial velocity
                    pdch = inputParams.PRMRad[0] * cnst
                    primaryStar.workingDataValues3 = fshift.fshift(nptsPrimSt,primaryStar.workingDataValues2,pdch)
                    
                    ############################################### Same for secondary:
                    if (nstar > 1):
                        
                        supp.scale(nptsSecSt,secondaryStar.workingDataValues1,swt,secondaryStar.workingDataValues3)
                        # print ("Rotation parameter [secondary star]: " + str(inputParams.SECRot[0]))
                        secondaryStar.workingDataValues2 = starot.starot(nptsSecSt, secondaryStar.workingDataValues3, inputParams.SECRot[0], 
                                                                        midLambda, deltaStepProc)#, secondaryStar.workingDataValues2 )
                        # print ("Radial Velocity parameter [secondary star]: " + str(inputParams.SECRad[0]))
                        sdch = inputParams.SECRad[0] * cnst
                        secondaryStar.workingDataValues3 = fshift.fshift(nptsSecSt,secondaryStar.workingDataValues2,sdch)
            
                    ############################################### And for the tertiary:
                    if (nstar == 3):
                        supp.scale(nptsTerSt,tertiaryStar.workingDataValues1,twt,tertiaryStar.workingDataValues3)
                        # print ("Rotation parameter [tertiary star]: " + str(inputParams.TERRot[0]))
                        tertiaryStar.workingDataValues2 = starot.starot(nptsTerSt, tertiaryStar.workingDataValues3,inputParams.TERRot[0], 
                                                                        midLambda, deltaStepProc)#,tertiaryStar.workingDataValues2 )
                        # print ("Radial Velocity parameter [tertiary star]: " + str(inputParams.TERRad[0]))
                        tdch = inputParams.TERRad[0] * cnst
                        tertiaryStar.workingDataValues3 = fshift.fshift(nptsTerSt,tertiaryStar.workingDataValues2,sdch)
            
            
                    #////////////////////////////////////////////////////////////////////////////////////
                    #-----------------------------------------------------------------------------------
                    # Add to get first model,calculate the sum of the squared residuals (ssr) and
                    # set initial ssr as minimum ssr
                    #-----------------------------------------------------------------------------------
                    #///////////////////////////////////////////////////////////////////////////////////
                    xmean = xPrimStarmean*inputParams.PRMWeight[0] + xSecStarmean*inputParams.SECWeight[0] + xTerStarmean*inputParams.TERWeight[0]
                    # minssr = 10.0
                    # It should be that: npts == nptsInputObj == nptsPrimSt == nptsSecSt == nptsTerSt 
                    # But just in case, is included the safeguard of the while loop in the case nstar == 1
                    if (niter == 0):
                        if nstar == 1: #IF nstar = 1 AND pwt = 1 THERE IS NO NEED TO MAKE THE CALL TO sumthm
                            while len(primaryStar.workingDataValues2) < len(inputObject.workingDataValues2):
                                primaryStar.workingDataValues2.append(0.0)
                            inputObject.workingDataValues2 = cpy.copy(primaryStar.workingDataValues3)
                        elif nstar == 2:
                            inputObject.workingDataValues2 = supp.sumthm(npts,primaryStar.workingDataValues3,secondaryStar.workingDataValues3,zeroes)    
                        else: #nstar == 3:
                            inputObject.workingDataValues2 = supp.sumthm(npts,primaryStar.workingDataValues3,secondaryStar.workingDataValues3,tertiaryStar.workingDataValues3)
                        
                        ndebugLen1 = len(inputObject.workingDataValues1)
                        ndebugLen2 = len(inputObject.workingDataValues2)
                        ndebugLen3 = len(inputObject.workingDataValues3)
                        
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
                        print("--------------------------------------------------->    searching...")
                        rv_value_iter = rv_values.find_mjd(inputObject.date)
                        if rv_value_iter == None:
                            tst[0][0] = 1
                            rv_values_available = False
                        elif nstar == 1:
                            print("----------------------------------------------->    found...", rv_value_iter[1])
                            tst[0][0] = 0
                            inputParams.PRMRad[0] = rv_value_iter[1]
                            pdch = inputParams.PRMRad[0] * cnst
                            primaryStar.workingDataValues3 = fshift.fshift(nptsPrimSt, primaryStar.workingDataValues2, pdch)
                    #------------------------------------------------------------------------------------------------------------
                    #Shift the Primary
                    if (tst[0][0] != 0):
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
                            deltaparam = deltaparam/10
                            parm = bstprm[0][1] - (10 * deltaparam)
                        inputParams.PRMRot[0] = bstprm[0][1]
                        primaryStar.workingDataValues3 = starot.starot(nptsPrimSt, primaryStar.workingDataValues2,inputParams.PRMRot[0],
                                                                    midLambda, deltaStepProc)#, primaryStar.workingDataValues3)
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
                            deltaparam = deltaparam/10
                            parm = bstprm [1][1] - (10 * deltaparam)
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
                            for i in range(22): 
                                if (parm >= 0.0):
                                    tertiaryStar.workingDataValues3 = starot.starot(nptsTerSt, tertiaryStar.workingDataValues2, parm, midLambda, deltaStepProc)#, tertiaryStar.workingDataValues3)
                                    inputObject.workingDataValues2 = supp.sumthm(nptsTerSt,primaryStar.workingDataValues3, secondaryStar.workingDataValues3,tertiaryStar.workingDataValues3)
                                    supp.scale2 (nptsPrimSt,inputObject.workingDataValues2,xObjmean/xmean)
                                    sr,ssr,inputParams.cnt = supp.sumres(nptsInputObj,inputObject.workingDataValues1, inputObject.workingDataValues2, 0, nptsPrimSt-1, skp)
                                    if (ssr < minssr):
                                        minssr = ssr
                                        bstprm[1][1] = parm
                                parm = parm + deltaparam 
                            deltaparam = deltaparam/10
                            parm = bstprm [2][1] - (10.0 * deltaparam)
                        inputParams.TERRot[0] = bstprm[1][1]
                        tertiaryStar.workingDataValues3 = starot.starot(nptsTerSt, tertiaryStar.workingDataValues2, inputParams.TERRot[0], 
                                                                midLambda, deltaStepProc)
                    #end if
                    
                    #####################################################                       
                    #***     Find the best weights 
                    #####################################################
                    if ((nstar > 1) and (nwt >= 1)):
                        # minssr = 100
                        #Set up arrays to start
                        primaryStar.workingDataValues3 = starot.starot(npts, primaryStar.workingDataValues1, inputParams.PRMRot[0],
                                                                    midLambda, deltaStepProc)# former use of ndata instead of npts
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
                    # print (' Results of iteration #' + str(niter+1))
                    print (' Minimum RMS of the residuals = ' + str('{:8.5f}'.format(math.sqrt (math.fabs(minssr) / inputParams.cnt))))
                    print ('-----------------------------------------------')
                    print ('             PRI      SEC      TER')
                    print (' V_rad   = ' + str('{:6.3f}'.format(inputParams.PRMRad[0]))    + "    " + str('{:6.3f}'.format(inputParams.SECRad[0]))    + "    " + str('{:6.3f}'.format(inputParams.TERRad[0])) )
                    print (' V_rot   = ' + str('{:6.3f}'.format(inputParams.PRMRot[0]))    + "    " + str('{:6.3f}'.format(inputParams.SECRot[0]))    + "    " + str('{:6.3f}'.format(inputParams.TERRot[0])) )
                    print (' Weight  = ' + str('{:6.3f}'.format(inputParams.PRMWeight[0])) + "    " + str('{:6.3f}'.format(inputParams.SECWeight[0])) + "    " + str('{:6.3f}'.format(inputParams.TERWeight[0])) )
                    print ('-----------------------------------------------')
                    partial_results_var = []
                    partial_results_var.append([niter+1, str('{:8.5f}'.format(math.sqrt (math.fabs(minssr) / inputParams.cnt))), 0])
                    partial_results_var.append([str('{:6.3f}'.format(inputParams.PRMRad[0])), str('{:6.3f}'.format(inputParams.PRMRot[0])), str('{:6.3f}'.format(inputParams.PRMWeight[0]))])
                    partial_results_var.append([str('{:6.3f}'.format(inputParams.SECRad[0])), str('{:6.3f}'.format(inputParams.SECRot[0])), str('{:6.3f}'.format(inputParams.SECWeight[0]))])
                    partial_results_var.append([str('{:6.3f}'.format(inputParams.TERRad[0])), str('{:6.3f}'.format(inputParams.TERRot[0])), str('{:6.3f}'.format(inputParams.TERWeight[0]))])
                    
                    frm.update_partial_results_frm(partial_results_var)
                    frm.idle()
                    ################################################################
                    #                   End iteration loop
                    ################################################################
            
                if (not frm.stop_istarmod_event.is_set()):
                        
                    ################################################################
                    #                 Calculate final model                         
                    ################################################################
                    supp.scale(nptsPrimSt,primaryStar.workingDataValues1,inputParams.PRMWeight[0],primaryStar.workingDataValues3)
                    primaryStar.workingDataValues2 = starot.starot(nptsPrimSt, primaryStar.workingDataValues3, inputParams.PRMRot[0], 
                                                                midLambda, deltaStepProc)# , primaryStar.workingDataValues2)
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
                        inputObject.workingDataValues2 = cpy.copy(primaryStar.workingDataValues3) 
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
                    if(inputParams.writeOutputSpec):
                        wp = {'OBJECT': str("SYNTHMODIF_" + primaryStar.object) ,'MJD':inputObject.date, "NAXIS":1,"NAXIS1": nptsInputObj,"CRVAL1":str(inputObject.initialLambda), "CDELT1":str(inputObject.deltaStep), "CRPIX1":1, "CTYPE1":"LINEAR"}
                        outputFilename = "@" + str(inputObject.date) +"_" + inputParams.linespcf
                        if inputParams.synthSpecName != None:
                            if (".fits" in inputParams.synthSpecName.name):
                                synthSpecName0 = inputParams.synthSpecName.name[:-5]
                            elif (".dat" in inputParams.synthSpecName.name):
                                synthSpecName0 = inputParams.synthSpecName.name[:-4]
                            else:
                                synthSpecName0 = inputParams.synthSpecName.name
                            synthSpecNamePath = inputParams.workpath / "OUTFITS"
                            if not os.access(synthSpecNamePath, os.F_OK):
                                os.mkdir(Path(synthSpecNamePath))
                            synthSpecName_fits = os.path.join(synthSpecNamePath, synthSpecName0 + outputFilename + '.fits')
                            inputObject.writeToFits(Path(synthSpecName_fits), wp, 2)
                            
                            synthSpecNamePath = os.path.join(inputParams.workpath, "OUTDAT")
                            if not os.access(synthSpecNamePath, os.F_OK):
                                os.mkdir(Path(synthSpecNamePath))
                            synthSpecName_dat = os.path.join(synthSpecNamePath, synthSpecName0 + outputFilename + '.dat')
                            inputObject.writeToDat (Path(synthSpecName_dat), 2)

                    ####################################################################
                    #         Substract model from object
                    ####################################################################
                    for i in range(nptsInputObj):
                        #Substracted Spectrum = Initial Read Spectrum - Synthetic Spectrum
                        inputObject.workingDataValues3[i] = (inputObject.workingDataValues1[i] - inputObject.workingDataValues2[i])
                    inputObject.nanSanityCheck(3)        
                
                    ###################################################################
                    #         Write subtracted spectrum if desired
                    ###################################################################
                    if(inputParams.writeOutputSpec):
                        wp = {'OBJECT': str("SUBTRACTED_" + inputObject.object) ,'MJD':inputObject.date, "NAXIS":1,"NAXIS1": nptsInputObj,"CRVAL1":str(inputObject.initialLambda), "CDELT1":str(inputObject.deltaStep), "CRPIX1":1, "CTYPE1":"LINEAR"}
                        if inputParams.subtractedSpecName != None:
                            
                            if (".fits" in inputParams.subtractedSpecName.name):
                                subtractedSpecName0 = inputParams.subtractedSpecName.name[:-5]
                            elif (".dat" in inputParams.subtractedSpecName.name):
                                subtractedSpecName0 = inputParams.subtractedSpecName.name[:-4]
                            else:
                                subtractedSpecName0 = inputParams.subtractedSpecName.name

                            subtractedSpecPath = inputParams.workpath / "OUTFITS"
                            if not os.access(subtractedSpecPath, os.F_OK):
                                os.mkdir(Path(subtractedSpecPath))
                            subtractedSpecName_fits = os.path.join(subtractedSpecPath, subtractedSpecName0 + outputFilename + '.fits')
                            inputObject.writeToFits(Path(subtractedSpecName_fits), wp, 3)
                            
                            subtractedSpecPath = os.path.join(inputParams.workpath, "OUTDAT")
                            if not os.access(subtractedSpecPath, os.F_OK):
                                os.mkdir(Path(subtractedSpecPath))
                            subtractedSpecName_dat  = os.path.join(subtractedSpecPath, subtractedSpecName0 + outputFilename + '.dat')
                            inputObject.writeToDat (Path(subtractedSpecName_dat) , 3)

                    ###################################################################
                    #                That's all folks!
                    ###################################################################
                    print (' ')
                    
                    ##------------------------RV Out File-------------------------------------------##
                    if nstar == 1 and (inputParams!= None and inputParams.use_rv_values):
                        if (rv_values_available == False):
                            rv_outfile = open(filenameRV1, 'a')
                            if inputParams.linespcf != 'NONE':
                                rv_outfile.write(str(inputObject.date) + " " + str(inputParams.PRMRad[0]) + '\n')
                            else:
                                rv_outfile.write(str(inputObject.date) + " " + str(inputParams.PRMRad[0]) + " " + str(inputObject.BaryCorr) + '\n')
                            print ("##########################    outfilenameRV:" , filenameRV1)
                            rv_outfile.close()
                    ##------------------------------------------------------------------------------##

                    ###############################################################################
                    ########--- Calculate the EW of the line (or lines) ---########################
                
                    lambdas = tools.lambdaData(Path("data/lambdas.dat").name, inputParams.ldo1_value, inputParams.ldo2_value)
                    if(inputParams.linespcf != 'NONE'):                             
                        ldoLineInput = lambdas.setLineInputbyKey(inputParams.linespcf)
                        if debugging: print("Line: ", inputParams.linespcf)
                        try:
                            if(nstar == 1):
                                if inputObject.RVfromFITS != 0.0: ## The CARMENES provided spectrum contains these parameters in the fits headings
                                    inputBaryCorr = inputObject.BaryCorr
                                    inputValue = inputObject.RVfromFITS + inputBaryCorr
                                    if debugging: print ("inputValue = ",inputValue)
                                else:                            ## General spectrum. No RV or barycentric correction available in fits headings
                                    inputBaryCorr = inputObject.BaryCorr = 0
                                    inputValue = inputParams.PRMRad[0]
                                    if debugging: print ("inputValue (no RVfromFITS) = ",inputValue)
                                ldoObs = lambdas.calculateDopplerDisplacement(inputValue)
                                if debugging: print ("LDO from Doppler displacement:", ldoObs)
                                ldoObs, param2, param5, param6, param3, param4 = lambdas.getMaximumValue(inputObject, inputValue, inputParams.FindMax_Tolerance) # ldoObs is modified within getMaximumValues
                            elif (nstar == 2 and inputParams.ldo1_value ==-1 and inputParams.ldo2_value == -1):
                                inputValue     = inputParams.PRMRad[0]
                                inputValue2    = inputParams.SECRad[0]
                                ldoObs,ldoObs2 = lambdas.calculateDopplerDisplacement(inputValue,inputValue2, True)
                                if debugging: print (lambdas.ldoObs, lambdas.ldoObs2)
                                if lambdas.ldoObs > lambdas.ldoObs2: #preparing for getMaximumValue
                                    lambdas.ldoObs, lambdas.ldoObs2 = lambdas.ldoObs2 , lambdas.ldoObs
                                    ldoObs , ldoObs2 = ldoObs2 , ldoObs 
                                if debugging: 
                                    print (lambdas.ldoObs, lambdas.ldoObs2)
                                    print (ldoObs, ldoObs2)
                                ldoObs, param2, param5, param6, ldoObs2, param4 = lambdas.getMaximumValue(inputObject, lambdas.ldoObs, inputParams.FindMax_Tolerance, 2, lambdas.ldoObs2) # ldoObs is modified within getMaximumValues
                            elif (nstar == 2 and inputParams.ldo1_value !=-1 and inputParams.ldo2_value != -1):
                                # Values "manually" selected. IT MUST BE SELECTED BOTH ldoX_value's 
                                lambdas.ldoObs  = ldoObs  = inputParams.ldo1_value
                                lambdas.ldoObs2 = ldoObs2 = inputParams.ldo2_value
                                param2 = 1
                                param4 = 1
                                factorStart = 1-inputParams.FindMax_Tolerance*10
                                factorEnd   = 1+inputParams.FindMax_Tolerance*10
                                param5 = factorStart*ldoObs2
                                param6 = factorEnd*ldoObs

                            if (inputObject.RVfromFITS != 0.0):
                                if debugging: print("ldoObs = ", ldoObs, " !!! ")
                            else:
                                if debugging: print(ldoObs, param2)
                                vRad = lambdas.calculateRadialVelocityfromDopplerDisplacement()
                                ldoLineInput = lambdas.setLineInputbyValue(ldoObs)
                            if nstar == 1:
                                ldo =  ldoObs    
                            elif nstar == 2:
                                if ldoObs2 != None:
                                    ldo =  (ldoObs2 + ldoObs)/2
                                else:
                                    ldo =  (lambdas.ldoObs2 + lambdas.ldoObs)/2
                                # ldo = (ldoObs + ldoObs2)/2
                            
                            ldoLineStart = lambdas.setLineStart(ldo - ldoLineInput[1])
                            ldoLineEnd   = lambdas.setLineEnd  (ldo + ldoLineInput[1])
                            if (inputObject.RVfromFITS == 0.0 and debugging):
                                print ("ldoLineInput= ", ldoLineInput[0], "\nldoLineStart= ", ldoLineStart, "ldoLineEnd= ", ldoLineEnd)

                            if debugging: print ("deltaStep = ", deltaStepProc)
                            ## The following allows to select the subset of points of the subtracted spectrum on which we calculate the EW
                            lambdas.selectLambdaRanges(inputObject.workingLambdaValues, deltaStepProc)
                            
                            EWpopt = lambdas.CalculateEW(inputObject, nstar, debugging, initial_weight=inputParams.PRMWeight[0], initial_weight2=inputParams.SECWeight[0])
                            
                            EW1, EW2, star1_spectrum, star2_spectrum, flambdavalues, summ_up, abs_error, rel_error = EWpopt
                            
                            FWHM = lambdas.getFWHM(inputObject, debugging)
                            SNR = max(inputObject.SNR, primaryStar.SNR)
                            errorEW = lambdas.calculateErrorCayrelFormulae(FWHM,deltaStepProc,SNR)
                            equivalent_widths_var = []
                            print("#################################################")
                            print("     FWHM = ",FWHM)
                            if nstar == 1:
                                print("     EW   = " + str(EW1) + " +- " + str(errorEW))
                                equivalent_widths_var.append(str('{:6.3f}'.format(EW1)))
                                equivalent_widths_var.append(str('--'))
                                equivalent_widths_var.append(str('--'))
                            elif nstar == 2: # and EW1 != None and EW2 != None:
                                if EW1 == None: EW1 = 0.0
                                if EW2 == None: EW2 = 0.0
                                if abs_error == None: abs_error = 0.0
                                if rel_error == None: rel_error = -1.0
                                print("     EW1   = " + str(EW1) + r" \pm " + str(abs_error) + ' (' + str('{:3.1f}'.format(rel_error)) + '%)')
                                print("     EW2   = " + str(EW2) + r" \pm " + str(abs_error) + ' (' + str('{:3.1f}'.format(rel_error))+ '%)')
                                equivalent_widths_var.append(str('{:6.3f}'.format(EW1)))
                                equivalent_widths_var.append(str('{:6.3f}'.format(EW2)))
                                equivalent_widths_var.append(str('--'))
                            print("#################################################")
                            frm.update_final_results(equivalent_widths_var)

                            ##################### Addition 10/11/2020
                            if(inputParams.deltaStep == -1.0):
                                print("DELTA_STEP = ", deltaStepAvg) 
                                print("CONST_ = ", cnst, '\n\n')
                            #########################################

                            # found = len(inputParams.workpath) - inputParams.workpath.find('\\')
                            if (debugging == True): print("IM_PATH: ", inputParams.workpath)
                            if (inputParams.resultspath == None or inputParams.resultspath.name == '') :
                                inputParams.resultspath = Path("./")
                            else:
                                if not os.access(Path(inputParams.resultspath), os.F_OK): 
                                    os.mkdir(Path(inputParams.resultspath))
                            if (debugging == True): print(Path(inputParams.resultspath))
                            respath = os.path.join(Path(inputParams.resultspath, "RES"))
                            if (debugging == True): print("RES_PATH: ", Path(respath))    

                            if not os.access(Path(respath), os.F_OK):
                                os.mkdir(Path(respath))
                            ewOutfile = open(os.path.join(respath, inputObject.object + "_" + inputParams.linespcf), 'a')
                            if (debugging == True): print (Path(ewOutfile.name))
                            if nstar == 1:
                                ewOutfile.write(str(inputObject.date)+ "  " + str(EW1) + "  " + str(errorEW) + "  " + str(inputParams.PRMRad[0]) + "  " + str(inputObject.BaryCorr) + "  " + str(inputParams.PRMRot[0]) + '\n')
                            elif nstar == 2 and EW1 != None and EW2 != None:
                                ewOutfile.write(str(inputObject.date)+ "  " + str(EW1) + "  " + str(abs_error) + "  " + str(EW2) + "  " + str(abs_error) + "  " + str(inputParams.PRMRad[0]) + "  " + str(inputObject.BaryCorr) + "  " + str(inputParams.PRMRot[0]) + "  " + str(inputParams.SECRad[0]) + "  " + str(inputParams.SECRot[0]) + '\n')
                            else: ## nstar == 3
                                ewOutfile.write(str(inputObject.date)+ "  " + str(EW1) + "  " + str(errorEW) + "  "+ str(inputParams.PRMRad[0]) + "  " + str(inputObject.BaryCorr) + "  " + str(inputParams.PRMRot[0]) + "  " + str(inputParams.SECRot[0]) + "  " + str(inputParams.SECRad[0]) + "  " + str(inputParams.TERRot[0]) + str(inputParams.TERRad[0]) + '\n')
                            
                        except ValueError:
                            PRINT ("##########################################################")
                            print ("\n NOT POSSIBLE TO DISCERN COMPONENTS OF THE BINARY/TRIPLE")
                            Print ("EWs NOT CALCULATED")
                            PRINT ("##########################################################")

                    if(inputParams.linespcf == 'NONE'):
                        ldoObs  = None
                        lambdas.ldoObs2 = None
                        param2 = None
                        param4 = None
                        param5 = None
                        param6 = None
                        summ_up = None
                        star1_spectrum = None
                        star2_spectrum = None
                        flambdavalues = None

                    payload = {
                        "inputParams": inputParams,
                        "nstar": nstar,
                        "date": inputObject.date, 
                        "object": inputObject.object,
                        "lambdas": lambdas,
                        "ldoObs": ldoObs,
                        "ldoObs2": lambdas.ldoObs2,
                        "param2": param2,
                        "param4": param4,
                        "param5": param5, 
                        "param6": param6,
                        "workingLambdaValues": inputObject.workingLambdaValues, 
                        "workingDataValues3":  inputObject.workingDataValues3,
                        "workingDataValues1":  inputObject.workingDataValues1, 
                        "workingDataValues2":  inputObject.workingDataValues2,
                        "summ_up":             summ_up, 
                        "flambdavalues":       flambdavalues, 
                        "star1_spectrum":      star1_spectrum, 
                        "star2_spectrum":      star2_spectrum,
                        "debugging": debugging
                    }


                    frm.plot_queue.put(payload)
                  
                    #### Closng all I/O files
                    if (".fits" in str(objectInputFITSfromList)):   inputObject.fits_hdu_in.close()

                    primaryStar.fits_hdu_in.close()
                    if (nstar == 2):    secondaryStar.fits_hdu_in.close()
                    if (nstar == 3):    tertiaryStar.fits_hdu_in.close()
                    
                    if (inputParams.linespcf != 'NONE'):     ewOutfile.close()
                    # plt.close()
                    
                    def _next_spectrum_():
                        sleep(.5)
                        print("WAITING FOR THE NEXT SPECTRUM TO BE PROCESSED" )
                        

                    # t = Timer(0.5,_next_spectrum_)
                    # t.start()
                    # sleep(0.5)
                    # inputParams = tools.readConfigFile(spcfFile)
                else: 
                    print("Here4!!!. Stopped execution")
                    return
            ## NEXT SPECTRUM FILE
        
        except FileNotFoundError as e:
            print("File: ", e.filename, " Not Found.")
            print("Please, Check the Input Parameters in the .sm file")
            # payload["error"] = str("File: ", e.filename, " Not Found.")
            return
        except PermissionError as e:
            if e.errno == EACCESS:
                print("Cannot access ", e.filename, "file")
                return
            payload["error"] = str("Cannot access ", e.filename, "file")
        
        return

    def _build_plot_(self, payload):

        ################ THIS FUNCTION ISOLATES PLOTTING FROM THE ALGORTIHM ##########
        ################ the payload conatains all paramters and objects needed to plot the different figures #######

        inputParams = payload["inputParams"]
        nstar       = payload["nstar"]
        _date       = payload["date"]
        _object     = payload["object"]
        lambdas     = payload["lambdas"]
        ldoObs      = payload["ldoObs"]
        ldoObs2     = payload["ldoObs2"]
        param2      = payload["param2"]
        param4      = payload["param4"]
        param5      = payload["param5"] 
        param6      = payload["param6"]
        debugging   = payload["debugging"]
        workingLambdaValues     = payload["workingLambdaValues"]
        workingDataValues13     = payload["workingDataValues3"]
        workingDataValues21     = payload["workingDataValues1"]
        workingDataValues22     = payload["workingDataValues2"]
        summ_up                 = payload["summ_up"]
        flambdavalues           = payload["flambdavalues"] 
        star1_spectrum          = payload["star1_spectrum"]
        star2_spectrum          = payload["star2_spectrum"]
        

        # These lines must/should be commented in the worker thread in order not to Starting a Matplotlib GUI outside of the main thread 
        # (And then fail)
        f, (subPlot1,subPlot2) = plt.subplots(2, figsize = (10,7), dpi = 150, sharex = True)#, sharey = True)
        # subPlot1.clear()
        # subPlot2.clear()
        if(inputParams.linespcf != 'NONE'):
            if nstar == 1:
            ## FOR DEBUGGING PURPOSES. To draw the boundaries of the zone where to search for the maximum.
                if (debugging == True):
                    subPlot1.plot(ldoObs, param2, 'b+')
                    subPlot1.axvline(x=param5, c = 'r', linestyle='dashed')
                    subPlot1.axvline(x=param6, c = 'r', linestyle='dashed')
            elif (nstar == 2):
                if (debugging == True):
                    subPlot1.axvline(x=param5, c = 'r', linestyle='dashed')
                    subPlot1.axvline(x=param6, c = 'r', linestyle='dashed')
                    subPlot1.axvspan(param5, param6, color='green', alpha=0.1)
            ###############################################################################
            ########## To draw the integration limits #####################################
            subPlot1.axvline(x=workingLambdaValues[lambdas.start], linestyle='dashed')
            subPlot1.axvline(x=workingLambdaValues[lambdas.end],   linestyle='dashed')
            if debugging == True:
                subPlot1.axvspan(workingLambdaValues[lambdas.start], 
                                    workingLambdaValues[lambdas.end]  , color='blue', alpha=0.1)
            subPlot1.axhline(y=0.0, linestyle='dashed')
            subPlot2.axvline(x=workingLambdaValues[lambdas.start], linestyle='dashed')
            subPlot2.axvline(x=workingLambdaValues[lambdas.end],   linestyle='dashed')
            ###############################################################################
            if nstar ==2:
                    subPlot1.plot(flambdavalues, star1_spectrum, 'c-', lw = 1.5)
                    subPlot1.plot(flambdavalues, star2_spectrum, 'y-', lw = 1.5)
                    if (debugging == True):## FOR DEBUGGING PURPOSES
                        if (ldoObs2 != None and param4 != None):
                            subPlot1.plot(ldoObs2, param4, 'r*') 
                        if(ldoObs != None and param2 != None):
                            subPlot1.plot(ldoObs, param2, 'b+')
                        if summ_up != None:
                            subPlot1.plot(flambdavalues, summ_up, c = 'm', linestyle = 'dashed', lw = 1.25)


        if len(workingLambdaValues) == len(workingDataValues22):
            subPlot2.plot(workingLambdaValues, workingDataValues21, 'b-', lw = 1.0)
            subPlot2.plot(workingLambdaValues, workingDataValues22, 'r-', lw = 1.0)
            subPlot1.plot(workingLambdaValues, workingDataValues13, 'g-', lw = 1.0)
        elif len(workingLambdaValues) >= len(workingDataValues22):
            subPlot2.plot(workingLambdaValues, workingDataValues21, 'b-', lw = 1.0)
            subPlot2.plot(workingLambdaValues[0:len(workingDataValues22)], workingDataValues22, 'r-', lw = 1.0)
            subPlot1.plot(workingLambdaValues, workingDataValues13, 'g-', lw = 1.0)
        elif len(workingLambdaValues) <= len(workingDataValues22):
            subPlot2.plot(workingLambdaValues, workingDataValues21, 'b-', lw = 1.0)
            subPlot2.plot(workingLambdaValues, workingDataValues22[0:len(workingLambdaValues)], 'r-', lw = 1.0)
            subPlot1.plot(workingLambdaValues, workingDataValues13, 'g-', lw = 1.0)

        f.subplots_adjust(hspace = 0.3)

        dateTitle = str('{:6.4f}'.format(_date))
        
        strTitle = str(_object) + "@" + str(dateTitle) + "_" + inputParams.linespcf 
        if (inputParams.wvl_display_range[0] != -1 and inputParams.wvl_display_range[1] != -1): 
            subPlot1.set_xlim(inputParams.wvl_display_range[0],inputParams.wvl_display_range[1])
        # # CaIIH&K
        if(inputParams.linespcf == "CaIIH" or inputParams.linespcf == "CaIIK"):
            subPlot1.set_title("Subtracted Spectrum for "+ inputParams.linespcf + r" / $v_{\rm rot}sini$ = " + str('{:7.4f}'.format(inputParams.PRMRot[0])) + " km/s", font = 'serif', size= 18)
            subPlot1.set_ylim(-1.0, inputParams.MaxFluxDisplayed_Sub)
            subPlot2.set_title(_object + "@" + dateTitle +"  &  Synthetic Spectra for " + inputParams.linespcf, font = 'serif', size= 18)
            subPlot2.set_ylim(-0.5, inputParams.MaxFluxDisplayed_Obj)
            # subPlot1.set_xlim(3945,3985)
            # subPlot2.set_xlim(3945,3985)
        # # Hbeta
        elif(inputParams.linespcf == "Hbeta"):
            subPlot1.set_title("Subtracted Spectrum for "+ r"$H\beta$" + r" / $v_{\rm rot}sini$ = " + str('{:7.4f}'.format(inputParams.PRMRot[0])) + " km/s", font = 'serif', size= 18)
            subPlot1.set_ylim(-1.0, inputParams.MaxFluxDisplayed_Sub)
            subPlot2.set_title(_object + "@" + dateTitle + "  &  Synthetic Spectra for " + r"$H\beta$", font = 'serif', size= 18)
            subPlot2.set_ylim(-0.5, inputParams.MaxFluxDisplayed_Obj)
            if (inputParams.wvl_display_range[0] == -1 and inputParams.wvl_display_range[1] == -1):
                subPlot1.set_xlim(4850,4890)
                
        # # HeD3 NaD2/NaD1
        elif(inputParams.linespcf == 'HeD3' or inputParams.linespcf == 'NaD2' or inputParams.linespcf == 'NaD1'):
            subPlot1.set_title(r"Subtracted Spectrum for "+ inputParams.linespcf + r" / $v_{\rm rot}sini$ = " + str('{:7.4f}'.format(inputParams.PRMRot[0])) + " km/s", font = 'serif', size= 18)
            subPlot1.set_ylim(-1.0, inputParams.MaxFluxDisplayed_Sub)
            subPlot2.set_title(_object + "@" + dateTitle + "  &  Synthetic Spectra for " + inputParams.linespcf, font = 'serif', size= 18)
            subPlot2.set_ylim(-0.5, inputParams.MaxFluxDisplayed_Obj)
            if (inputParams.wvl_display_range[0] == -1 and inputParams.wvl_display_range[1] == -1):
                subPlot1.set_xlim(5860,5910)

        # # Halpha
        elif(inputParams.linespcf == 'Halpha'):
            subPlot1.set_title(r"Subtracted Spectrum for "+ r"$H\alpha$" + r" / $v_{\rm rot}sini$ = " + str('{:7.4f}'.format(inputParams.PRMRot[0])) + " km/s", font = 'serif', size= 18)
            subPlot2.set_title(_object + "@" + dateTitle + "  &  Synthetic Spectra for " + r"$H\alpha$", font = 'serif', size= 18)
            subPlot2.set_ylim(0.0, inputParams.MaxFluxDisplayed_Obj)
            subPlot1.set_ylim(-1.0, inputParams.MaxFluxDisplayed_Sub)
            if (inputParams.wvl_display_range[0] == -1 and inputParams.wvl_display_range[1] == -1):
                subPlot1.set_xlim(6550,6580)

        #  CaIRT-a
        elif(inputParams.linespcf == 'CaIRT-a'):
            subPlot1.set_title(r"Subtracted Spectrum for "+ inputParams.linespcf + r" / $v_{\rm rot}sini$ = " + str('{:7.4f}'.format(inputParams.PRMRot[0])) + " km/s", font = 'serif', size= 18)
            subPlot2.set_title(_object + "@" + dateTitle +"  &  Synthetic Spectra for " + inputParams.linespcf, font = 'serif', size= 18)
            subPlot1.set_ylim(-0.75, inputParams.MaxFluxDisplayed_Sub)
            subPlot2.set_ylim(0.25,inputParams.MaxFluxDisplayed_Obj)
            if (inputParams.wvl_display_range[0] == -1 and inputParams.wvl_display_range[1] == -1):
                subPlot1.set_xlim(8465,8535)

        elif(inputParams.linespcf == 'CaIRT-b'):
            subPlot1.set_title(r"Subtracted Spectrum for "+ inputParams.linespcf + r" / $v_{\rm rot}sini$ = " + str('{:7.4f}'.format(inputParams.PRMRot[0])) + " km/s", font = 'serif', size= 18)
            subPlot2.set_title(_object + "@" + dateTitle +"  &  Synthetic Spectra for " + inputParams.linespcf, font = 'serif', size= 18)
            subPlot1.set_ylim(-0.75, inputParams.MaxFluxDisplayed_Sub)
            subPlot2.set_ylim(0.0,inputParams.MaxFluxDisplayed_Obj)
            if (inputParams.wvl_display_range[0] == -1 and inputParams.wvl_display_range[1] == -1):
                subPlot1.set_xlim(8535,8550)

        elif(inputParams.linespcf == 'CaIRT-c'):
            subPlot1.set_title(r"Subtracted Spectrum for "+ inputParams.linespcf + r" / $v_{\rm rot}sini$ = " + str('{:7.4f}'.format(inputParams.PRMRot[0])) + " km/s", font = 'serif', size= 18)
            subPlot2.set_title(_object + "@" + dateTitle +"  &  Synthetic Spectra for " + inputParams.linespcf, font = 'serif', size= 18)
            subPlot1.set_ylim(-0.5, inputParams.MaxFluxDisplayed_Sub)
            subPlot2.set_ylim(0.,inputParams.MaxFluxDisplayed_Obj)
            if (inputParams.wvl_display_range[0] == -1 and inputParams.wvl_display_range[1] == -1):
                subPlot1.set_xlim(8600,8680)

        # # #################### Paschen D
        elif(inputParams.linespcf == 'PaschenD'):
            subPlot2.set_title(_object + "@" + dateTitle +"  &  Synthetic Spectra for " + inputParams.linespcf, font = 'serif', size= 18)
            subPlot2.set_ylim(0.25, inputParams.MaxFluxDisplayed_Obj)
            subPlot1.set_title(r"Subtracted Spectrum for "+ inputParams.linespcf + r" / $v_{\rm rot}sini$ = " + str('{:7.4f}'.format(inputParams.PRMRot[0])) + " km/s", font = 'serif', size= 18)
            subPlot1.set_ylim(-0.25, inputParams.MaxFluxDisplayed_Sub)
            if (inputParams.wvl_display_range[0] == -1 and inputParams.wvl_display_range[1] == -1):
                subPlot1.set_xlim(10040,10065)

        #  ##################### HeI10830
        elif(inputParams.linespcf == 'HeI10833'):
            subPlot2.set_title(_object + "@" + dateTitle +"  &  Synthetic Spectra for " + inputParams.linespcf, font = 'serif', size= 18)
            subPlot2.set_ylim(0.25, inputParams.MaxFluxDisplayed_Obj)
            subPlot1.set_title(r"Subtracted Spectrum for "+ inputParams.linespcf + r" / $v_{\rm rot}sini$ = " + str('{:7.4f}'.format(inputParams.PRMRot[0])) + " km/s", font = 'serif', size= 18)
            subPlot1.set_ylim(-0.25, inputParams.MaxFluxDisplayed_Sub)
            if (inputParams.wvl_display_range[0] == -1 and inputParams.wvl_display_range[1] == -1):
                subPlot1.set_xlim(10820,10860)

        #  ##################### PaschenG
        elif(inputParams.linespcf == 'PaschenG'):
            subPlot2.set_title(_object + "@" + dateTitle +"  &  Synthetic Spectra for " + inputParams.linespcf, font = 'serif', size= 18)
            subPlot2.set_ylim(0.25, inputParams.MaxFluxDisplayed_Obj)
            subPlot1.set_title(r"Subtracted Spectrum for "+ inputParams.linespcf + r" / $v_{\rm rot}sini$ = " + str('{:7.4f}'.format(inputParams.PRMRot[0])) + " km/s", font = 'serif', size= 18)
            subPlot1.set_ylim(-0.25, inputParams.MaxFluxDisplayed_Sub)
            if (inputParams.wvl_display_range[0] == -1 and inputParams.wvl_display_range[1] == -1):
                subPlot1.set_xlim(10920,10955)

        #  ##################### PaschenB
        elif(inputParams.linespcf == 'PaschenB'):
            subPlot2.set_title(_object + "@" + dateTitle +"  &  Synthetic Spectra for " + inputParams.linespcf, font = 'serif', size= 18)
            subPlot2.set_ylim(0.5, inputParams.MaxFluxDisplayed_Obj)
            subPlot1.set_title(r"Subtracted Spectrum for "+ inputParams.linespcf + r" / $v_{\rm rot}sini$ = " + str('{:7.4f}'.format(inputParams.PRMRot[0])) + " km/s", font = 'serif', size= 18)
            subPlot1.set_ylim(-0.5, inputParams.MaxFluxDisplayed_Sub)
            if (inputParams.wvl_display_range[0] == -1 and inputParams.wvl_display_range[1] == -1):
                subPlot1.set_xlim(12800, 12840)
                
        else:
        # Any
            subPlot2.set_title(_object + "@" + dateTitle +"  &  Synthetic Spectra for " + inputParams.linespcf, font = 'serif', size= 18)
            subPlot2.set_ylim(0.0, inputParams.MaxFluxDisplayed_Obj)
            subPlot1.set_title(r"Subtracted Spectrum for "+ inputParams.linespcf + r" / $v_{\rm rot}sini$ = " + str('{:7.4f}'.format(inputParams.PRMRot[0])) + " km/s", font = 'serif', size= 18)
            subPlot1.set_ylim(-0.5, inputParams.MaxFluxDisplayed_Sub)
            # subPlot1.set_xlim(9326 , 9485)
            # subPlot2.set_xlim(9326 , 9485)
            
        plt.xlabel(r"$\lambda$  [${\rm \AA}$]", font = 'serif', size= 14)
        plt.ylabel("Normalised Flux", font= 'serif', size= 14)
        # plt.legend()
        # # At this point is supposed that a valid value of inputParams.resultspath has been assigned
        figure_repo = Path(inputParams.resultspath) / Path("RES") / Path(inputParams.linespcf)
        if not os.access(Path(figure_repo), os.F_OK): 
            os.mkdir(Path(figure_repo))
        figureName = os.path.join(figure_repo, strTitle)
        if nstar == 1:  
                figureName = figureName + ".png"
        elif nstar == 2:
                figureName = figureName + "_asSB2.png"
        
        figureName_p = Path(figureName)
        f.savefig(figureName_p, dpi = 300)
        
        return f, subPlot1, subPlot2
