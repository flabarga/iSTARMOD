# -*- coding: utf-8 -*-
import math
#import astropy
import numpy as np
#NEW###########
from fshift import*
#from scipy import signal 
#########################################################################
##-----------------------------------------------------------------------
##  This routine broadens a spectrum as a rotating star.
##  It assumes limb darkening from Gray (linear in mu) with eps = 0.6.
##  This may be made more general if something like Wade & Rucinski's
##  tables can be made available to the program.  This would require
##  only a new line in the starin file to specify the gravity of the
##  star.  (STARMOD already calculates rough wavelengths.)
##
##  Parameters:
##  ndata:        number of points of the data array
##  ydata:        data array to be rotationally broadened 
##  vsini:        rotational parameter (projected velocity)
##  midLambda:    midpoint of the wavelenght array (Lambda)
##  refinc:       Increment in wavelenght array per pixel
##  order:        resulting array from the perfomed operations (also returned)
##
##  Alan Welty              15-oct-89  new (for "disk")
##                          21-jun-93  modify for starmod
##                          03-jan-94  FFT version
##  Fernando Labarga        07-dec-16  Migration to Python
##                          17-apr-17  Modification of the setting of the broadening kernel
 
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
def starot (ndata, ydata, vsini, midLambda, refinc):
    
    #Initialising lists in order to allow their use as arrays à la FORTRAN
    order = []
    for i in range (ndata):
        order.append(0.0)

    phi  = []
    for i in range(ndata+1):
        phi.append(0.0)
    
    ans  = []
    for i in range(ndata):
        ans.append(0.0)
    
    c = 299792.458                     #speed of light in km. per sec.
    i=0
    dlm = midLambda * vsini / c        #Conversion Factor
    pix = dlm / refinc
                                       # if pix < 0.5 the value of the rotational velocity is too low 
                                       # to produce any effect of rotational bradening on the spectral lines
    if (pix < 0.5):
        for i in range (ndata):
            order [i] = ydata [i]
        return order                    
    
    npix =int(pix + 0.5)
    frac = (pix + 0.5) - float (npix)
    
    ##  Broadening function stuff
    eps = 0.6                                   #Gray parameter for Linear Limb Darkening Law
    cd = math.pi * dlm * (1.0 - (eps / 3.0))
    c1 = 2.0 * (1.0 - eps) / cd
    c2 = 0.5 * math.pi * eps / cd

    ##  CALCULATE GRAY PROFILE.
    
    ##  Central pixel
    #dl1 = 0.0
    #tmp1 = 1.0
    phi1 = c1 + c2 
    phi [0] = 0.0
    for j in range(1,11):
        dl2    = float(j) / 20.0 
        tmp2   = 1.0 - (dl2 / pix) ** 2
        phi2   = (c1 * math.sqrt(tmp2)) + (c2 * tmp2)
        phi[0] = phi[0] + phi1 + phi2
        dl1    = dl2
        tmp1   = tmp2
        phi1 = phi2

    phi[0] = phi[0] / 20.0
    norm = phi[0]

 
  ##  Whole pixels
    
    for i in range(1, npix):
        phi[i] = 0.0
        for j in range(-9,11,1):
            dl2 = float(i) + float(j) / 20.0
            tmp2 = 1.0 - (dl2 / pix) ** 2
            phi2 = (c1 * math.sqrt(tmp2)) + (c2 * tmp2)
            phi[i] = phi[i] + (phi1 + phi2) / 40.0
            #print "        phi[",i,"] =", phi[i], "j = ", j
            #dl1 = dl2
            #tmp1 = tmp2
            phi1 = phi2
        norm = norm + phi[i]
    
    ##  Partial pixel (remaining)
    phi [npix] = 0.0
    for j in range(-9,11,1):
        dl2 = float (npix) + float (j) / 20.0
        if (dl2 <= pix):
            tmp2 = 1.0 - (dl2 / pix) ** 2
            phi2 = (c1 * math.sqrt(tmp2)) + (c2 * tmp2)
            phi[npix] = phi[npix] + (phi1 + phi2) / 40.0
            #dl1 = dl2
            #tmp1 = tmp2
            phi1 = phi2
    norm = norm + phi[npix]
    m = npix * 2 + 1
    ## Wrap-around for FFT

    for i in range(npix+1, m):
        phi[i]  = phi[m-i]
        norm = norm + phi[i]
    
    #######################
    ## The Gray profile is halved so, we need to provide the other half of the line
    q = m+npix+1
    for i in range(m, q):
        phi[i]   = phi[i-m]
        phi[i-m] = 0.0
    for i in range(q):
        phi[i] = phi[i+npix+1]
        phi[i+npix+1] = 0.0
        #print "phi[",i,"] = ", phi[i]
    ################################
    
    ## Normalization 
    for i in range(m):
        phi[i] = phi[i] / norm
    
    ## Do the FFT convolution
    ans = np.convolve (ydata, phi, mode = 'full')
    
    #######################
    # Resulting from the execution of the algorithm 
    # (applying np.convolve) 
    # is an spectrum[array] broadened, but shifted. 
    #It is needed to rollback this shift
    cnst = (midLambda/299792.458)/refinc
    dchrot = -vsini * cnst
    ans = fshift(ndata,ans,dchrot)
    ################################
    
    for i in range(ndata):
        order[i] = ans[i]

    return order
