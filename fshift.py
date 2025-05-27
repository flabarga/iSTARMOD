# -*- coding: utf-8 -*-
########################################################################
#
#  Name:     fshift
#  Filename: fshift.py
#  Type:     subroutine.
#
#  Language: Python
#  Purpose:  Shift a data array by a non-integral amount.
#
#  #all line: fshift(npts,xin,shft,yout)
#                  npts (input) = # data points in array x.
#                  xin  (input) = input data array.
#                  shft (input) = shift value.
#                  yout (output) = shifted data array.
#
#  In#lude files: none.
#  Subroutines/libraries: ishift (npts,x,ishft,y) 
#                         interp (x,y,npts,nterms,xin,yout)  
#
#  Author: S. Barden   
#  Date:      nov 84
#  Modified: 
#
########################################################################

import math
from support import*
from interp import*
import scipy.interpolate as interpl
#

def ishift (npts,x,ishft):

    y = []
    for i in range(npts):
        y.append(0.0)
    for i in range(npts):
        j = i - ishft
        if (j < 0): j = 0
        if (j >= npts): j = npts-1
        y[i] = x[j]
    return y

      
def fshift (npts,xin,shft):
    
#-----------------------------------------------------------------------
    #In order to assign the list elements à la FORTRAN, we must initialise the lists
    yout = []
    for i in range(npts):
        yout.append(0.0)
    t  = []
    ty = []
    for i in range(npts):
        if i<4:
            ty.append(0.0)
        t.append(0.0)
    p  = []
    px = []
    for i in range(npts):
        if i<4:
            px.append(0.0)
        p.append(0.0)
#-----------------------------------------------------------------------
#                                       Break shft into integer 
#                                       and fractional parts;
    ishft = int(shft)
    fshft = shft % math.copysign(1.,shft)
#
#                                       Shift by integral amount;
    t = ishift(npts,xin,ishft)
#                                       Shift by fraction if non-zero.
    if (fshft == 0.0 ):
        for i in range(npts):
            yout[i] = t[i]

    else:
        for i in range(npts):
            p[i] = float(i)
        #interpolate = interpl.interp1d(p,t, "linear")
        for i in range(1,npts-1):
        #for i in range(npts-1):
            xinpt = float(i) - fshft
            #yout[i] = interpolate(xinpt)
            if (i-2 <= 0):
                ty[0] = 1000.0
                px[0] = -1.0
            else:
               ty[0] = t[i-2]
               px[0] = p[i-2]

            if (i-1 <= 0):
               ty[1] = 1000.0
               px[1] = 0.0
            else:
               ty[1] = t[i-1]
               px[1] = p[i-1]

            ty[2] = t[i]
            px[2] = p[i]
            if (i+1 > npts):
               ty[3] = 1000.0
               px[3] = npts + 1
            else:
               ty[3] = t[i+1]
               px[3] = p[i+1]

            inpts = 4
            intrms = 2
            # print(intrms)
            youtpt = interp(px,ty,inpts,intrms,xinpt)
            yout[i] = youtpt
        yout[0]  = yout[1]
        yout[-1] = yout[len(yout)-2]
    return yout
#-----------------------------------------------------------------------
#end

