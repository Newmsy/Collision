import matplotlib.pyplot as plt
import numpy as np

#CapillCollLastWall takes into account the last wall hit
#CapillCollY does not and just looks at the Y's
XPlotWidth=100

XVals = list(range(-XPlotWidth,XPlotWidth))

def Gauss(mean,std,xs,factor):
    ys=[]
    for i in xs:
        YVal = (1/(std*(2*np.pi)**0.5))*np.e**(-0.5*((i-mean)/std)**2)
        ys.append(YVal*factor)

    return ys

def SumGauss(mean,std,xs):
    YVals=Gauss(mean[0],std[0],xs)
    for n,m in enumerate(mean[1:]):
        YPlus=Gauss(m,std[n+1],xs)
        for i in range(len(xs)):
            YVals[i]+=YPlus[i]
    return YVals

def Lorentzian(mean,fwhm,xs,factor):
    ys=[]
    for i in xs:
        YVal=1/((np.pi*fwhm)*(1+((i-mean)/(fwhm))**2))
        ys.append(YVal*factor)
    return ys

means=[15,-15]
stds=[12.18,12.18]
plt.plot(XVals,Gauss(0,12.18,XVals,75000))
plt.plot(XVals,Lorentzian(0,9.5,XVals,115000))
#plt.show()

def SumLorentzian(means,stds,xs,factors):
    YVals=Lorentzian(means[0],stds[0],xs,factors[0])
    for n,m in enumerate(means[1:]):
        YPlus=Lorentzian(m,stds[n+1],xs,factors[n+1])
        for i in range(len(xs)):
            YVals[i]-=YPlus[i]
    return YVals

def MBoltz(mean,mkt,xs):
    ys=[]
    for i in xs:
        YVal=(mkt/np.pi)**(3/2)*(i-mean)**2 *np.e**(-mkt*(i-mean)**2)
        ys.append(YVal*500000)
    return ys

meansL=[0,0]
stdsL=[8.5,4]
factors=[161000,57500]
plt.plot(XVals,SumLorentzian(meansL,stdsL,XVals,factors))


with open(r'CollisionDataStore\FarCapillCollY.txt',mode='r') as filesaveY:
    HistHold=[]
    TotHold=[]
    AngHold=[]
    for n,line in enumerate(filesaveY):
        b=float(line[:-1])
        TotHold.append(b)
        if abs(b)<XPlotWidth:
            HistHold.append(b)
            AngHold.append(np.arctan(b/60)*180/np.pi)
    plt.hist(HistHold,100)
    #plt.axvline(3,color='r')
    #plt.axvline(-3,color='r')
    plt.show()
    print(np.std(HistHold))

    #plt.axvline(np.arctan(3/60)*180/np.pi,color='r')
    #plt.axvline(np.arctan(-3/60)*180/np.pi,color='r')
    plt.hist(AngHold,120)
    plt.plot(XVals,Gauss(0,12.18,XVals,50000))
    plt.plot(XVals,Lorentzian(0,10,XVals,60000))
    plt.show()
    print(n)
