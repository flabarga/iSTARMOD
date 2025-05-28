import math
import scipy.signal as signal 
import matplotlib.pyplot as plt

class RVValues():

    def __init__(self, filename = ""):
        self.mjd      = []
        self.rv_value = []
        self.values   = []
        self.filename = filename
        return 
        
    def readFile(self,filename):
        self.filename = filename
        with open(self.filename, 'r') as inFile:
            for line in iter(inFile.readline,''):
                itemStr       = line.split()
                self.mjd      = float(itemStr[0])
                self.rv_value = float(itemStr[1])
                self.values.append([self.mjd, self.rv_value])
        return 
        
    def find_mjd(self,mjd):
        for iter in self.values:
            if (math.fabs(mjd-iter[0]) < 0.01):
                print(mjd, iter[0])
                print ("####################################    Iter: ", iter[0])
                return iter
        return None
        
    def findPlanet(self, initial,step,nosteps, noread = False):
        
        if (not noread):
            self.readFile(self.filename)
        time     = []
        RV_kmday = []
        for item in self.values:
            RV_kmday.append(item[1]*86164)
        for item in self.values:
            if item[0] > 2400000:
                mjd = item[0] - 2400000
            else:
                mjd = item[0] 
            time.append(mjd)

        for it in range(len(self.values)):
            print( RV_kmday[it], time[it])   
        freq = []
        #Define the array of frequencies for which to compute the periodogram:
        # initial = 0.00005 
        # the input parameters say how we explore the frequencies dom 
        for i in range(nosteps):
            freq.append(initial+i*step)
            # print( freq[i]) 
        periodogram = signal.lombscargle(time[4:], RV_kmday[4:], freq, normalize = True, precenter = True)

        f, (subPlot1,subPlot2) = plt.subplots(2, figsize = (15,8), dpi = 90)
        
        subPlot1.plot(time[4:], RV_kmday[4:], 'r*')
        
        subPlot2.plot(freq, periodogram)
        print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^  ", len(self.values))        
        plt.show()
        
        
        
        
        
