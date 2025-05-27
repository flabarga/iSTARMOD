#########################################################################
#  
#  support.f 
#
#  This file contains several subroutines called by the main program.
#  They are:
#            scale   multiply an array by a constant
#            scale2  same as scale, but replaces the original array
#            sumthm  add three arrays point by point
#            sumres  calculate sum of squares of an array
#            gtmean  calculate average value of an array
#
#  The last two perform their calculations on selected regions of the
#  input arrays.  These regions are specified in the input file.
#
########################################################################
import math

def scale(npts,x,wt,y):
    for iter in range (npts):
        y[iter] = wt * x[iter]
    return
#######################################################################


def scale2(npts,x,wt):
    #print "scale2 npts: ",npts
    for i in range(npts):
        x[i] = wt * x[i]
    return x
#######################################################################


def sumthm(npts,a,b,c):
    y = []
    for i in range(npts):
        y.append(a[i] + b[i] + c[i])
    return y
#######################################################################


def sumres(npts,a,b,xlo,xhi,skp):
    
    resid = []
    ihi = -137
    ilo = -137
    cnt = 0
    j = 0
    
    while(ihi != xhi):
        if 2*j < len(skp) and skp[2*j] > ilo:
            ihi = skp[2*j] 
        else:
            ihi = xhi

        if (ihi > xhi):
            ihi = xhi
        for i in range(ilo, ihi):
            cnt += 1
            resid.append(b[i] - a[i])
        j+=1
        if 2*j-1 < len(skp):
            ilo = skp[2*j-1] +1
    
    ssr = 0.0
    for i in range(cnt+1,npts,1):
        resid.append(0.0)
        
    for i in range(cnt):
        ssr = ssr + resid[i]**2
 
    return ssr, cnt
#######################################################################
#
#      
def gtmean (npts,a,xlo,xhi,skp):
    #print
    #print "xlo = ", xlo, '    ', "xhi = ", xhi
    # print (skp)
    ilo = 0
    ihi = 0
    suma = 0.0
    cnt = 1
    j = 0
    xmean = 0.0
    ilo = xlo

    while(ihi != xhi):
        if 2*j < len(skp) and skp[2*j] > ilo:
            #print "j = ", j , "    skp[2*j] = " ,skp[2*j] 
            ihi = skp[2*j] 
            #print "new ihi = ", ihi
        else:
            #print "j = ", j 
            ihi = xhi

        if (ihi > xhi):
            ihi = xhi
        #print "ilo = ", ilo, "  ihi = ", ihi, "xhi = ", xhi
        for i in range(ilo,ihi):
            if math.isnan(a[i]) != True :
                cnt += 1
                suma += a[i]
           #print "+ ", a[i],
        #print
        j+=1
        if 2*j-1 < len(skp):
            ilo = skp[2*j-1] +1
            #print ilo , "#############"
        #print 'skp = ',skp[2*j-1]
        #if cnt != 0:
    xmean = suma / cnt 
    return xmean
#######################################################################

# x = []
# 
# #x[0] = 0
# y = []
# #y.append(137)
# npts = 25
# for itr in range(npts):
#     x.append(itr+1)
#     y.append(137)
#     print (itr, x[itr], y[itr]) 
# scale(npts,x,3,y)
# print (x)
# print (y)
# skp = []
# for itr in range(6):
#     skp.append(2*itr+3)
# print (skp)
# print ("gtmean = ", gtmean(npts,x,0,25,skp))
# for itr in range(npts):
#     y[itr]= 7
# print ('scale2 = ', scale2(npts,x,3))
# print (x)
# print (y)
# print ("gtmean = ", gtmean(npts,x,5,20,skp))
# print ("gtmean = ", gtmean(npts,y,0,25,skp))
# print (sumthm(npts,y,x,y))
# print (sumres(npts,x,y,5,20,skp))
# print (2**2)