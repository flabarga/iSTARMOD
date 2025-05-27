# -*- coding: utf-8 -*-
import math
#import astropy
import numpy as np
#NEW###########
from fshift import*
#from scipy import signal 
#from rotBroad import *
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
def starot (ndata, ydata, vsini, midLambda, refinc): #, order):
    
    #print ndata
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
    #print "pix, npix = ", pix, npix

    ##  Broadening function stuff
    eps = 0.6                                   #Gray parameter for Linear Limb Darkening Law
    cd = math.pi * dlm * (1.0 - (eps / 3.0))
    c1 = 2.0 * (1.0 - eps) / cd
    c2 = 0.5 * math.pi * eps / cd

    ##  CALCULATE GRAY PROFILE.
    ##
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
    #print "phi[",0,"] = ", phi[0]
    #print
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
        #print "phi[",i,"] = ", phi[i]
    #print

    ##  Partial pixel (remaining)
    #npix += 1
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
    #print "phi[",npix,"] = ", phi[npix]
    #print
    m = npix * 2 + 1
    ## Wrap-around for FFT

    for i in range(npix+1, m):
        phi[i]  = phi[m-i]
        norm = norm + phi[i]
        #print "phi[",i,"] = ", phi[i]
    #print
    
    ############NEW LINES###########
    ## The Gray profile is halved so, we need to provide the other half of the line
    q = m+npix+1
    for i in range(m, q):
        phi[i]   = phi[i-m]
        phi[i-m] = 0.0
        #print "phi[",i,"] = ", phi[i]
    #print
    for i in range(q):
        phi[i] = phi[i+npix+1]
        phi[i+npix+1] = 0.0
        #print "phi[",i,"] = ", phi[i]
    ################################
    
    ## Normalization 
    for i in range(m):
        phi[i] = phi[i] / norm
    
    #midp = len(phi)/2
    #for i in range(m):
    #    phi[i],phi[midp-m+i] = swap(phi[i],phi[midp-m+i])
        
    #print phi
    ## Do the FFT convolution
    ##
    #ans = convlv (ydata, ndata, phi, m, 1)
    ans = np.convolve (ydata, phi, mode = 'full')
    #ans = signal.convolve(ydata,phi, mode= "full")
    #print "/////////////////////////////////////////////"
    #print "phi = ", phi
    #print len(order), len(ans), ndata, len(phi)
    #print "/////////////////////////////////////////////"
    
    ############NEW LINES###########
    # Resulting from the execution of the algorithm 
    # (applying np.convolve) 
    # is an spectrum[array] broadened, but shifted. 
    #It is needed to rollback this shift
    cnst = (midLambda/299792.458)/refinc
    dchrot = -vsini * cnst
    ans = fshift(ndata,ans,dchrot)
    ################################
    
    for i in range(ndata):
        #print phi[i]
        order[i] = ans[i]

    return order # ,phi

def swap(a,b):
    temp = a
    a = b; b = temp
    return a,b
    
###############################################################################
###############################################################################
## ****************************************************************************
## ** The following routines are from Numerical Recipes and python translated *
## ****************************************************************************
##              NOT USED
##############################################################################
def four1(data,nn,isign):
    
    n = 2*nn
    j = 0
    #print "n = ", n        
    for i in range (0,n,2):
        #print "j= ",j
        if j > i:
            tempr     = data[j]
            tempi     = data[j+1]
            data[j]   = data[i]
            data[j+1] = data[i+1]
            data[i]   = tempr
            data[i+1] = tempi
        m=n/2
        while (m >= 2 and j > m):
            j = j - m
            m = m/2
        j = j + m -1
    mmax = 2
    while (n > mmax):
        istep = 2*mmax
        theta = (2*math.pi)/(isign*mmax)
        wpr   = -2.*math.sin(0.50*theta)**2
        wpi   = math.sin(theta)
        wr    = 1.0000000000
        wi    = 0.0000000000
        for m in range(0,mmax,2):
            for i in range(m,n,istep):
                j=i+mmax
                #print "j=",j
                if j<len(data)-1:
                    tempr = float(wr)*data[j] - float(wi)*data[j+1]                                                                                    
                    tempi = float(wr)*data[j+1] + float(wi)*data[j]                                                                          
                    data[j]   = data[i]   - tempr
                    data[j+1] = data[i+1] - tempi
                    data[i]   = data[i]   + tempr
                    data[i+1] = data[i+1] + tempi
            wtemp = wr
            wr    = wr*wpr - wi*wpi + wr
            wi    = wi*wpr + wtemp*wpi + wi
        mmax = istep
    return

##############################################################################
##############################################################################

def twofft(data1,data2,n):
    
    #Initialising lists in order to allow their use as arrays à la FORTRAN
    fft1 = []
    fft2 = []

    c1 = complex(0.5, 0.0)
    c2 = complex(0.0, -0.5)
    
    ##Accessing dataX[j] in the usual way
    ##Initialising the arrays fft1 & fft2
    for j in range(n-1):#        do j=1,n
        fft1.append(complex(data1[j],data2[j]))
        fft2.append(complex(0.,0.))
        
    four1(fft1,n/2,1)
    fft2[0] = complex(fft1[0].imag,0.0)
    fft1[0] = complex(fft1[0].real,0.0)
    n2 = n/2+2
    for j in range(1,(n/2), 1):
        h1 = c1*(fft1[j] + fft1[n2-j].conjugate())
        h2 = c2*(fft1[j] - fft1[n2-j].conjugate())
        fft1[j]      = h1
        fft1[n2-j-2] = h1.conjugate()
        fft2[j]      = h2
        fft2[n2-j-2] = h2.conjugate()
        
    return fft1, fft2
#

##############################################################################
##############################################################################

def realft(data,n,isign):
    
    theta = math.pi/float(n)
    c1    = 0.5
    if (isign == 1):
        c2 = -0.5
        four1(data,n,+1)
    else:
        c2 = 0.5
        theta = -theta

    wpr  = -2.0*math.sin(0.5*theta)**2
    wpi  = math.sin(theta)
    wr   = 1.0 + wpr
    wi   = wpi
    n2p3 = 2*n+1
    
    for i in range(0,n/2):
        i1  = 2*i
        i2  = i1 + 1
        i3  = n2p3 - i2
        i4  = i3 + 1
        wrs = float(wr)
        wis = float(wi)
        h1r =  c1*(data[i1]  + data[i3])
        h1i =  c1*(data[i2]  - data[i4])
        h2r = -c2*(data[i2]  + data[i4])
        h2i =  c2*(data[i1]  - data[i3])
        data[i1] =  h1r + wrs*h2r - wis*h2i
        data[i2] =  h1i + wrs*h2i + wis*h2r
        data[i3] =  h1r - wrs*h2r + wis*h2i
        data[i4] = -h1i + wrs*h2i + wis*h2r
        wtemp = wr
        wr    = wr*wpr - wi*wpi+wr
        wi    = wi*wpr + wtemp*wpi + wi

    if (isign == 1):
        h1r     = data[1]
        data[1] = data[2] + h1r
        data[2] = data[2] - h1r
    else:
        h1r     = data[1]
        data[1] = c1*(h1r + data[2])
        data[2] = c1*(h1r - data[2])
        four1(data,n,-1)
#   endif
    return data
    
##############################################################################
##############################################################################

def convlv(data,n,respns,m,isign):
    
    NMAX = 4096            
    #Initialising the list allow us to use them as arrays à la FORTRAN
    fft = []
    for i in range (NMAX):
        fft.append(complex(0.0))
    ans = []
    for i in range (n):
        ans.append(complex(0.0))

    for i in range(0, int(math.floor((m-1)/2))-1):
        #print "i in convolve: ", i, "and n-i:", n-i, "also, m-i: ", m-i
        respns[n-i-1] = respns[m-i-1]

    for i in range((m+3)/2, n-(m-1)/2):
        respns[i] = 0.0

    #print "input parameters passed to convolve are (n):", n, "(m): ", m 
    fft, ans = twofft(data,respns,n) #twofft(data,respns,fft,ans,n)
    no2 = n/2

    for i in range(0,no2):
        if (isign == 1):
            ans[i] = fft[i]*ans[i]/no2
        elif (isign == -1):
            if (math.fabs(ans[i]) == 0.0):
                #print 'deconvolving at resp=0'
                break
            ans[i] = fft[i]/ans[i]/no2
        else:
            #print 'no meaning for isign'
            break

    ans[0] = complex(ans[0].real,ans[no2+1].real)
    ans    = realft(ans,no2,-1)
    
    return ans
#        END
##
#

##############################################################################
######################test code###############################################
#a = complex(1.,2.)
#print a.real
#print a.imag 
#print a
##st = ""
##st += str(a)
##print st
#b = a.conjugate()
#print b
#c = 2*a
#d = b**2
#print c
#print d

