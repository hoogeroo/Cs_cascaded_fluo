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
        counts1=np.sum(y1[0:10])
        c1.append(counts1/thetime)
        counts2=np.sum(y1[10:20])
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

dns=sys.argv[1:]
nTrigs=120
#dn='202302031740'
fn=dns[0]+'/Data_0_0.fit'
data=fits.open(fn)
script=data[0].header['SCRIPT'].split(' ')
myxs=script[1]
#print(myxs)
myx=(eval('np.'+myxs)-80)*2
data.close()
nfiles=myx.shape[0]

#myx=myx[:nfiles]
c1=[]
c2=[]
shifts=[]
for dn in dns:
    t1,t2=rundir(dn,nfiles,2500,120)
    c1.append(t1)
    c2.append(t2)
    #plt.figure('Rates per mus')
    #plt.plot(myx,np.array(t2),'x',myx,np.array(t1),'o')#,myx,c2,'o')
    pguess=[10,1,10,1]
    p1,pc=curve_fit(lorentz,myx,t1,p0=pguess)
    #plt.plot(myx,lorentz(myx,*p1))
    print(p1)
    p2,pc=curve_fit(lorentz,myx,t2,p0=pguess)
    #plt.plot(myx,lorentz(myx,*p2))
    #plt.xlabel('Detuning (MHz)')
    #plt.ylabel('Counts per $\mu$s')
    print(p2)
    shift=p2[1]-p1[1]
    print('Shift is ',shift)
    shifts.append(shift)
plt.figure('Rates per mus',linewidth=2)
#plt.plot(myx,np.mean(np.array(c2),0),'x')
plt.errorbar(myx,np.mean(np.array(c2),0),np.std(np.array(c2),0),fmt='x',elinewidth=2)
plt.errorbar(myx,np.mean(np.array(c1),0),np.std(np.array(c1),0),fmt='o',elinewidth=2)
c1m=np.mean(np.array(c1),0)
c1s=np.std(np.array(c1),0)
p1,pc=curve_fit(lorentz,myx,c1m,p0=pguess)#,sigma=c1s)
plt.plot(myx,lorentz(myx,*p1))
print(p1)
c2m=np.mean(np.array(c2),0)
c2s=np.std(np.array(c2),0)
p2,pc=curve_fit(lorentz,myx,c2m,p0=pguess)#,sigma=c2s)
plt.plot(myx,lorentz(myx,*p2))
plt.xlabel('Detuning (MHz)')
plt.ylabel('Counts per $\mu$s')
print(p2)
#plt.plot(myx,np.mean(np.array(c1),0),'o')#,myx,c2,'o')
print("Mean shift is ",np.mean(shifts),'std', np.std(shifts))
cat=np.array(c1)
cab=np.array(c2)
#print(ct.shape)
rat=cab/cat
avrat=np.mean(rat,0)
strat=np.std(rat,0)/(np.sqrt(cat.shape[0]))

plt.figure('Ratio')
plt.errorbar(myx,avrat,strat,ls='none')
plt.plot(myx,avrat,'o')
p,pc=curve_fit(lorentz,myx,avrat,p0=[-1,0,10,1])
plt.plot(myx,lorentz(myx,*p),lw='2')
plt.xlabel('Detuning (MHz)')
plt.ylabel('Count ratio')
print('Mean parameters',p)
#plt.plot(myx,avrat,'o')
#plt.figure('Times (mus)')
#plt.plot(myx,ts)
plt.show()
