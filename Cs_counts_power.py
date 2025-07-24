# -*- coding: utf-8 -*-


import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
from scipy.optimize import curve_fit

plt.rcParams.update({'axes.titlesize': 'large'})

def lorentz(x,h,c,w,o):
    return h/(1+(2*(x-c)/w)**2)+o

def rundir(dn,nfiles,maxcounts,magic):
    c1=[]
    c2=[]
    ts=[]
    plt.figure(dn)
    for i in range(nfiles):
        fn=dn+'/Data_%i_0.fit' % i
        #print(fn)
        a=fits.open(fn)
        b=a[1].data
        filter=b>0
        c=b[filter]
        #print(c.shape)
        d1=np.sort(c)
        d2=d1[:maxcounts]
        e=d2 % magic
        thetime=d2[-1]*5/1000
        #print(thetime)
        ts.append(thetime)
        y1,x1=np.histogram(e,20)
        counts1=np.sum(y1[2:8])
        c1.append(counts1/thetime)
        counts2=np.sum(y1[12:18])
        c2.append(counts2/thetime)
        if (i==11) or (i==0):
            #plt.figure()
            #fig,ax=plt.subplots()
            if (i==0):
                ax=plt.subplot(1,2,1)
            else:
                ax=plt.subplot(1,2,2)
            plt.bar(x1[:-1]*5,y1,width=20)#,'x')
            plt.setp(ax.spines.values(),linewidth=3)
            ax.xaxis.set_tick_params(width=3)
            ax.yaxis.set_tick_params(width=3)
            plt.xlabel('Time (ns)')
            if (i==0):
                plt.ylabel('Count rate')
                plt.text(300,300,r'$\Delta=30$ MHz')
            else:
                plt.text(300,300,r'$\Delta=0$ MHz')
        a.close()
    return(c1,c2)

def satcurve(x,h,s):
    return h*x/s/(1+x/s)


if __name__=='__main__':
    x = np.array([1.5, 2.5, 4, 6, 8, 12, 16, 22, 29, 38, 48, 61, 77, 94, 135, 158, 184, 210, 241, 285, 320,
         351, 385, 425, 461, 505, 545, 588, 625, 664, 702, 744, 775, 815])
    dirs = ['202301201612', '202301221029', '202301221325', '202301302102']
    c1=[]
    c2=[]
    for mydir in dirs:
        a1,a2=rundir(mydir,34,2000,120)
        c1.append(a1)
        c2.append(a2)
    mc1=np.mean(np.array(c1),0)
    sc1=np.std(np.array(c1),0)
    mc2=np.mean(np.array(c2),0)
    sc2=np.std(np.array(c2),0)
    mc1=mc1-mc1[0]
    mc2=mc2-mc2[0]
    print(mc1[0],mc2[0])
    p1,cv1=curve_fit(satcurve,x,mc1,p0=[10,200])
    print(p1)
    p2,cv2=curve_fit(satcurve,x,mc2,p0=[10,200])
    print(p2)
    x=x/p1[1]
    #plt.figure('Saturation')
    fig,ax1=plt.subplots()
    ax2=fig.add_axes([0.55,0.18,0.3,0.3])
    ax1.plot(x,mc1,'o',x,mc2,'x',x,satcurve(p1[1]*x,*p1),x,satcurve(p1[1]*x,*p2))
    plt.setp(ax1.spines.values(),linewidth=3)
    ax1.xaxis.set_tick_params(width=3)
    ax1.yaxis.set_tick_params(width=3)
    ax1.set_xlabel('$s_0$')
    ax1.set_ylabel('count rate (a.u.)')
    #.figure('Ratio_s')
    #ax=plt.subplot(336)
    plt.setp(ax2.spines.values(),linewidth=3)
    ax2.xaxis.set_tick_params(width=3)
    ax2.yaxis.set_tick_params(width=3)
    ax2.plot(x[1:],mc2[1:]/mc1[1:],'o')
    #plt.xlabel('$s_0$')
    plt.ylabel('Count ratio')
    plt.show()
