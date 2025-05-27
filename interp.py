#############################################################################
# 
#  Name:     interp
#  Filename: interp.f
#
#  Language: Python
#  Purpose:  Interpolate between data points to evaluate a function, 
#            using a Lagrange Interpolation Polynomial 
# 
#
#  Call line: interp(x,y,npts,nterms,xin,yout) 
# 
#  Input Parameters:
#
#            float   x      := array of independent variable.
#            float   y      := array for dependent variable.
#            integer npts   := number of pairs of data points.
#            integer nterms := number of terms in fitting polynomial.
#            float   xin    := input x value.
#  Output Parameters:
#            float   yout   := interpolated value of y.
#            
#  Include files: none.
#  Subroutines/libraries: none.
# 
#  Author: Adapted by F. Labarga in Python
#          from       S. Barden in Fortran77
#             
#  Date:      Nov 84 / Dec 2016
#  Modified:
# 
#  Comments:.
# 
#############################################################################
# 
#import math
def interp (x,y,npts,nterms,xin):
    #print " npts = ", npts
    delta = [-1.0]*npts
    ##print delta
    #print "nterms = ", nterms
    yout = 0.0
#   Search for appropriate value of x(1);
#   Calculation of the closed interval [i1,i2]
#   Trying that xin be in the middle of the interval
    i1 = 0
    for i in range(npts):
        #print i, nterms/2, nterms
        if (xin-x[i] < 0): #130, 170, 190
            i1 = i - int(nterms/2) #130
            if (i1<=0):   #150, 150, 210
                i1 = 0   #150
                #print "< i1=",i1
            else: 
                #print "< i1=",i1, " and break"
                break
                #go to 210
        elif (xin-x[i] == 0):
            #print "== i1= ",i1
            yout = y[i] #170 
            return yout #go to 999
        else: #if (xin-x[i] >= 0):
            i1 = npts - nterms - 1  #190
            #print "> i1=",i1
    
    i2 = i1 + nterms   #210 
    # print ("i1 = ", i1)
    #print "i2 = ", i2
    
    #To avoid the effect of defining the interpolation at the ends of the sample    
    if (npts-i2 < 0): #230, 310, 310
        i2 = npts #230 
        #print "i2 prev = ",i2, "nterms = ", nterms
        i1 = i2 - nterms -1
        #print "i1 prev = ", i1
        if (i1 <= 0):  #260, 260, 310
            i1 = 0 #260 
        nterms = i2 - i1 
        # print ("i1 = ", i1)
        #print "i2 = ", i2
        
    ###Calculate the number of terms of the interpolating polynomial 
    ##print "nterms = ", nterms
    ## 310
    #print "***************************"
    # print(nterms)
    for j in range(nterms):
        #print "j =" ,j
        prod = 1.0
        for i in range(nterms):
            if (( i != j) and ((i+i1 < npts and j+i1 < npts))):# and ((x[j+i1] - x[i+i1])!= 0.0))):
            #The second condition is somewhat unneccesary
            #if ( i != j):
                prod = prod * ((xin - x[i+i1]) / (x[j+i1] - x[i+i1]))
        delta[j] = prod
        ##print delta
    ##print ""
    ##print delta
    Lagrange = 0.0
    for i in range(nterms):
        Lagrange += delta[i] * y[i+i1]
        ##print "delta[", i,"]= ", repr(delta[i]).rjust(22),"  y[", i,"]= ", repr(float(y[i])).rjust(17), "    Lagrange = ", repr(Lagrange).ljust(22)  
    ##print ""
    yout = Lagrange
    return yout

#test code
# x = [0.0, 0.785398163397448, 1.570796327, 2.35619449, 3.141592654, 3.926990817, 4.71238898, 5.497787144, 6.283185307]
# y = [0.0, 0.707106781186547 ,1.0, 0.707106781, 0.0, -0.707106781, -1.0, -0.707106781, 0] #1.22515E-16
# xin = 1.49
# print (x)
# print (y)
# yout = interp(x,y,9,12,xin)
# print ("")
# print ("f(", xin, ")= ",yout)