# -*- coding: utf-8 -*-
######################################################################################################################
##################################### iStarmod_ChiFactor.py ##########################################################
#################################    Labarga, F // Montes, D.   ######################################################
######################################################################################################################

import math


def calculateLogXi(  Teff):
    # Calculate the value of Xi based on an aproximation by a exponential Schechter-type function(Teff)
    # Definition of Xi factor as in Walcowicz et al, 2008
    # log Xi = C + alpha* log(Teff) + beta*(Teff)
    
    polyCoeff = [-46.29, 12.81, -1.006e-3]
    LogXi = polyCoeff[0] + polyCoeff[1] * math.log10(Teff) + polyCoeff[2] * Teff 
    
    return LogXi 

def calculateXiHK(Teff):
    # The same calculation as the one for Halpha, but for CaII H&K
    # log Xi = C + alpha*log(Teff) + beta1*(Teff) + beta2*(Teff^2)
    # The function returns Xi
    polyCoeff = [-1.747e+02, 5.434e+01, -7.759e-03, 3.356e-07]
    LogXiHK = polyCoeff[0] + polyCoeff[1] * math.log10(Teff) + polyCoeff[2] * Teff + polyCoeff[3] * Teff*Teff
    return math.pow(10,LogXiHK)

def calculateXiHeD3(  Teff):
    # The same calculation as the one for Halpha, but for HeD3
    # log Xi = C + alpha*log(Teff) + beta1*(Teff) + beta2*(Teff^2) + beta3*(Teff^3)
    # The function returns Xi
    polyCoeff = [-447.10, 149.10, -3.328e-2, 2.888e-6, -1.125e-10]
    LogXi = polyCoeff[0] + polyCoeff[1]*math.log10(Teff) + polyCoeff[2]*Teff + polyCoeff[3]*math.pow(Teff, 2) + polyCoeff[4]*math.pow(Teff,3)
    return math.pow(10,LogXi)

def vis (Teff):
    Xi     = math.pow(10, calculateLogXi(Teff))
    return Xi   

def calculateLogXiCaIRT(  Teff):
    # The same calculation as the one for Halpha, but for CaIRT
    # log Xi = C + alpha*log(Teff) + beta1*(Teff) + beta2*(Teff^2) + beta3*(Teff^3)
    polyCoeff = [-85.93, 27.91, -0.006747, 0.0000006246 ,-0.00000000002619]
    LogXi = polyCoeff[0] + polyCoeff[1]*math.log10(Teff) + polyCoeff[2]*Teff + polyCoeff[3]*math.pow(Teff,2) + polyCoeff[4]*math.pow(Teff,3)
    return math.pow(10,  LogXi)

def calculateLogXiPaschenD(Teff):
    # The same calculation as the one for Halpha, but for PaschenDelta
    # log Xi = C + (alpha=0)*log(Teff) + beta1*(Teff) + beta2*(Teff^2) + beta3*(Teff^3) + beta4*(Teff^4) + beta5*(Teff^5)
    polyCoeff = [-1.261e+01, 9.693e-03, -4.249e-06, 9.020e-10, -9.344e-14, 3.777e-18]
    LogXi = polyCoeff[0] + polyCoeff[1]*Teff + polyCoeff[2]*Teff*Teff + polyCoeff[3]*math.pow(Teff,3) + polyCoeff[4]*math.pow(Teff,4) + polyCoeff[5]*math.pow(Teff,5)
    return math.pow(10,  LogXi)

def nir(Teff):
    # The same calculation as the one for Halpha, but for HeI10830/PaschenGamma
    # log Xi = C + alpha*log(Teff) + beta1*(Teff) + beta2*(Teff^2) + beta3*(Teff^3)
    polyCoeff = [4.525, -2.733, 3.790e-04, -2.360e-08] ## HeI 10830
    # polyCoeff = [4.306, -2.693, 4.175e-04, -2.852e-08] ## PaschenG 
    LogXi = polyCoeff[0] + polyCoeff[1]*math.log10(Teff) + polyCoeff[2]*Teff + polyCoeff[3]*math.pow(Teff,2)
    return math.pow(10,  LogXi)

def calculateLogXiPaschenB(Teff):
    # The same calculation as the one for Halpha, but for PaschenBeta
    # log Xi = C + (alpha=0)*log(Teff) + beta1*(Teff) + beta2*(Teff^2) + beta3*(Teff^3) + beta4*(Teff^4) + beta5*(Teff^5)
    # polyCoeff = [-3.732, 3.791e-04, -4.191e-07, 1.252e-10, -1.621e-14, 7.671e-19]
    polyCoeff = [6.372, -3.364, 5.096e-04, 3.533e-08, 0.0, 0.0, 0.0]
    LogXi = polyCoeff[0] + polyCoeff[1]*math.log10(Teff) + polyCoeff[2]*Teff + polyCoeff[3]*math.pow(Teff,2) + polyCoeff[4]*math.pow(Teff,3) + polyCoeff[5]*math.pow(Teff,4) + polyCoeff[6]*math.pow(Teff,5)
    return math.pow(10,  LogXi)

def calculateXi(Teff, spLine):
    switcher = {
        "H&K":  calculateXiHK,
        "HeD3":     calculateXiHeD3,
        "NaD1":     calculateXiHeD3,
        "NaD2":     calculateXiHeD3,
        "Halpha":   vis,
        "CaIRT-a":  calculateLogXiCaIRT,
        "CaIRT-b":  calculateLogXiCaIRT,
        "CaIRT-c":  calculateLogXiCaIRT,
        "PaschenD": calculateLogXiPaschenD,
        "HeI10833": nir,
        "PaschenG": nir,
        "PaschenB": calculateLogXiPaschenB
    }
    function = switcher.get(spLine, lambda: "Invalid spLine")
    Xi = function(Teff)
    return Xi
    
def calculateLogLHalpha_Lbol(  EWHalpha, Teff):
    if (EWHalpha != None):
        logEW  = math.log10(math.fabs(EWHalpha))
        logXi  =  calculateLogXi(Teff)
        return (logXi + logEW)
    else:
        return None

def calculateLogLLine_Lbol( Teff, spLine, EWLine):
    if (EWLine != None):
        logEW  = math.log10(math.fabs(EWLine))
        logXi =  math.log10(calculateXi(Teff, spLine))
        return (logXi + logEW)
    else:
        return None

def calculateLogLHalpha_Lbol_Err(  TeffErr, EW, EWErr):
    if (EW != None and EWErr != None):
        relativeEWErr =  EWErr / EW
        return (math.pow(TeffErr, 5)*3.90255*1e-16 + relativeEWErr)
    return 
    
def calculateLogLLine_Lbol_Err( TeffErr, EWLine, EWLineErr):
    if (EWLine != None and EWLineErr != None):
        LogEWErr  =  EWLineErr / EWLine
        return (math.pow(TeffErr, 5.0)*3.90255*1e-16 + LogEWErr)




#######################################################################
##################  EXAMPLE CODE  #####################################
#######################################################################
def calculadorChi_LogL(Teff, TeffErr, spLine, EW= 0.0, EWErr= 0.0):
    
    Xi = calculateXi(Teff, spLine)
    print("Xi = ", (Xi))

    LogLLine_Lbol = 0.0
    if(EW != 0.0 and EWErr != 0.0):
        LogLLine_Lbol      =  calculateLogLLine_Lbol(Teff, spLine,  EW )
        LogLLine_LbolErr   =  calculateLogLLine_Lbol_Err(TeffErr, EW, EWErr)
        print("LogLLine_LBol    =  ", '{:9.6f}'.format(LogLLine_Lbol))
        print("LogLLine_LBolErr =  ", '{:8.5f}'.format(LogLLine_LbolErr))
    return

# calculadorChi_LogL(3231, 82, "Halpha")
# calculadorChi_LogL(2764, 58, "H&K", -3.993, 0.003)

# print("\nHIP79796:")
# calculadorChi_LogL(3687.9, 50, "Halpha", -2.29800, 0.172) ## HIP79796 
# print("\nEV Lac:")
# calculadorChi_LogL(3306.0, 50, "Halpha", -4.48667, 0.003) ## HIP112460 EV Lac
# print("\nHD152751:")
# calculadorChi_LogL(3250.0, 50, "Halpha", -2.11000, 0.120) ## HD152751
# print("\nHIP84794:")
# calculadorChi_LogL(3100.0, 50, "Halpha", -2.06000, 0.115) ## HIP84794
# print("\nGJ388:")
# calculadorChi_LogL(3455.0, 50, "Halpha", -3.42200, 0.122) ## GJ388


# print("\n\nHIP118212:")
# calculadorChi_LogL(3733.0, 50, "Halpha", 0.1, 0.02) 
# print("\nHIP104383:")
# calculadorChi_LogL(3600, 50, "Halpha", 0.235, 0.0525) 
# print("\nHD216899:")
# calculadorChi_LogL(3798.0, 50, "Halpha", -0.11, 0.040) 
# print("\nHIP108752:")
# calculadorChi_LogL(3434.0, 50, "Halpha", -0.18, 0.030) 
# print("\nGJ734B:")
# calculadorChi_LogL(3175.0, 50, "Halpha", -0.19, 0.070) 

# print("\n\n\nHIP62686 original value:")
# calculadorChi_LogL(4784.0, 50, "Halpha", -0.47, 0.11) ## HIP 62686
# print("\n\n\nHIP62686 corrected value:")
# calculadorChi_LogL(4526, 66, "Halpha", -0.47, 0.11) ## HIP 62686

# print("\n\n\nHD21845 original value:")
# calculadorChi_LogL(5587.8, 50, "Halpha", -0.34, 0.04) ## HD 21845
# print("\n\n\nHD21845 corrected value:")
# calculadorChi_LogL(5084.0, 5.9, "Halpha", -0.34, 0.04) ## HD21845

# print("\n\n\nHIP45383 original value:")
# calculadorChi_LogL(4809.7, 50, "Halpha", -0.245, 0.023) ## HIP 45383
# print("\n\n\nHIP45383 corrected value:")
# calculadorChi_LogL(5181.0, 50, "Halpha", -0.245, 0.023) ## HIP 45383

# print("\n\n\nHD238224 original value:")
# calculadorChi_LogL(4372.5, 50, "Halpha", -0.38, 0.04) ## HD 238224
# print("\n\n\nHD238224 corrected value:")
# calculadorChi_LogL(4526, 66, "Halpha", -0.38, 0.04) ## HD 238224



print("\nHIP79796:")
calculadorChi_LogL(3687.9, 50, "CaIRT-a", -0.4000, 0.036) ## HIP79796 M1.0
print("\nEV Lac:")
calculadorChi_LogL(3306.0, 50, "CaIRT-a", -0.3900, 0.060) ## HIP112460 EV Lac
print("\nHD152751:")
calculadorChi_LogL(3250.0, 50, "CaIRT-a", -0.1000, 0.003) ## HD152751
print("\nHIP84794:")
calculadorChi_LogL(3100.0, 50, "CaIRT-a", -0.1200, 0.020) ## HIP84794
print("\nGJ388:")
calculadorChi_LogL(3455.0, 50, "CaIRT-a", -0.2900, 0.016) ## GJ388 AD Leo


# print("\n\nHIP118212:")
# calculadorChi_LogL(3733.0, 50, "CaIRT-a", 0.1, 0.02) 
print("\nHIP104383:")
calculadorChi_LogL(3600.0, 50, "CaIRT-a", -0.1125, 0.035) 
print("\nHD216899:")
calculadorChi_LogL(3798.0, 50, "CaIRT-a", -0.1300, 0.030) 
print("\nHIP108752:")
calculadorChi_LogL(3434.0, 50, "CaIRT-a", -0.0700, 0.020) 
print("\nGJ734B:")
calculadorChi_LogL(3175.0, 50, "CaIRT-a", -0.0600, 0.070) 

print("\n\n\nHIP62686 original value:")
calculadorChi_LogL(4784.0, 50, "CaIRT-a", -0.1600, 0.050) ## HIP 62686
print("\n\n\nHIP62686 corrected value:")
calculadorChi_LogL(4526.0, 66, "CaIRT-a", -0.1600, 0.050) ## HIP 62686

print("\n\n\nHD21845 original value:")
calculadorChi_LogL(5587.8, 50., "CaIRT-a", -0.2700, 0.040) ## HD 21845
print("\n\n\nHD21845 corrected value:")
calculadorChi_LogL(5084.0, 5.9, "CaIRT-a", -0.2700, 0.040) ## HD21845

print("\n\n\nHIP45383 original value:")
calculadorChi_LogL(4809.7, 50, "CaIRT-a", -0.1410, 0.0130) ## HIP 45383
print("\n\n\nHIP45383 corrected value:")
calculadorChi_LogL(5181.0, 50, "CaIRT-a", -0.1410, 0.0130) ## HIP 45383

print("\n\n\nHD238224 original value:")
calculadorChi_LogL(4372.5, 50, "CaIRT-a", -0.1400, 0.050) ## HD 238224
print("\n\n\nHD238224 corrected value:")
calculadorChi_LogL(4526.0, 66, "CaIRT-a", -0.1400, 0.050) ## HD 238224