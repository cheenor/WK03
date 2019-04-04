#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 13 21:22:34 2015

@author: jhchen
"""
import string
import matplotlib.pyplot as plt
import numpy as np
import datetime
from pylab import *
import scipy.stats as scista
mpl.rcParams['ytick.labelsize'] = 20
#
def readAscii(fpath,iskp,*nl):
    #iskp  the total line skipped of the file
    # fpath   the full path of the file
    # usage: onedim=readAscii(fpaht,iskp)
    onedim=[]
    linesplit=[]
    f=open(fpath)
    if nl:
        nrl=nl[0]
        ff=f.readlines()[iskp:nrl]  ## first line in obs file is legend 
    else:
        ff=f.readlines()[iskp:]
    for line in ff:
        line=string.lstrip(line)
        linesplit.append(line[:-1].split(' '))
    for lnstrs in linesplit:
        for strs in lnstrs:
#            if strs!='':
#                onedim.append(string.atof(strs))
            try:
                onedim.append(string.atof(strs))
            except:
                strs=""
#            else:
#                onedim.append(string.atof(strs))
    del linesplit,ff
    f.close()
    print len(onedim)
    return onedim
def getmonmean(a,days):
    mon=np.zeros(shape=(12),dtype=float)
    monstd=np.zeros(shape=(12),dtype=float)
    k=0
    for im in range(0,12):
        nd=days[im]
        tmp=[]
        for idy in range(0,nd):
            mon[im]=mon[im]+a[k]/nd
            tmp.append(a[k])
            k=k+1
        monstd[im]=np.std(tmp)
    return mon,monstd
#
def getpct(a,pct):
    n=len(a)
    n1=int(n*(1-pct*0.01))
    b=sorted(a)
    bmin=b[0:n1]
    n2=n-n1-1
    bmax=b[n2:n]
    c=[]
    for i in range(0,n):
        tmp=a[i]
        if tmp not in bmin and tmp not in bmax:
            c.append(a[i])
    return c
def getday(cm,days):
    if cm==0:
        k0=0
        return k0
    else:
        k0=0
        for im in range(0,cm-1):
            k0=k0+days[im]
        return k0
def getseason(a,days):
    sen=np.zeros(shape=(4),dtype=float)
    std=np.zeros(shape=(4),dtype=float)
    for isn in range(0,4):
        if isn==0:
            tm=[2,3,4]
        elif isn==1:
            tm=[5,6,7]
        elif isn==2:
            tm=[8,9,10]
        elif isn==3:
            tm=[11,0,1]
        tmp=[]
        for imm in range(0,len(tm)):
            cm=tm[imm]
            nd=days[cm]
            k0=getday(cm,days)
            for idd in range(0,nd):
                k=k0+idd
                tmp.append(a[k])            
        sen[isn]=np.mean(tmp)
        std[isn]=np.std(tmp)
    return sen,std
def Mgetseason(a):
    sen=np.zeros(shape=(4),dtype=float)
    std=np.zeros(shape=(4),dtype=float)
    for isn in range(0,4):
        if isn==0:
            tm=[2,3,4]
        elif isn==1:
            tm=[5,6,7]
        elif isn==2:
            tm=[8,9,10]
        elif isn==3:
            tm=[11,0,1]
        tmp=[]
        for im in tm:
            tmp.append(a[im])            
        sen[isn]=np.mean(tmp)
        std[isn]=np.std(tmp)
    return sen,std
###########################################################
regions=['PRD','MLYR','ETP']
case='2D'
years=[2010]
months=[1]
days=[1]
nt=96*365  # 2880
nx=202
nga=len(regions)
dirsim='/Volumes/DATA1/Models/CRM/YEAR/'
direra='/Volumes/Disk_D/Work/MyPaper/2DForcing/Forcing/'
dirtrmm='/Volumes/MacHD/Work/Papers/WK03/DATA/'
dirout='/Volumes/MacHD/Work/Papers/WK03/pics/'
monstr=['Jan','Feb','Mar','Apr','May','June','July','Aug','Sept','Oct','Nov','Dec']
senstr=['Spring','Summer','Autumn','Winter']
ndays=np.zeros(shape=(12),dtype=int)
for i in range(0,12):
    ndays[i]=31
ndays[1]=28
ndays[3]=30
ndays[5]=30
ndays[8]=30
ndays[10]=30
color_cycle=['deeppink', 'lime', 'b', 'y','indigo', 'cyan']
wd=[2,2,2,2,2]
for iga in range(0,nga):
    iy=years[iga]
    ndays[1]=28
    if iy%400==0 or (iy%4==0 and iy%100!=0):
        ndays[1]=29
    else:
        ndays[1]=28   
    im=months[iga]
    jd=days[iga]
    area=regions[iga]
    yearstr="%4d"%iy
    datestr="%4d"%iy+"%2.2d"%im+"%2.2d"%jd+'_365d'
    dirin=dirsim+area+'/postdata/'
    dirin='/Volumes/MacHD/Work/Papers/WK03/DATA/'+regions[0]+'/'
#------------------------------------------------------------------------------
    fpath=dirin+'preci_'+area+case+'.txt'
    iskp=0
    onedim=readAscii(fpath,iskp)
    preci=np.ndarray(shape=(nt),dtype=float)
    pbl=np.ndarray(shape=(nt),dtype=float)
    fsh=np.ndarray(shape=(nt),dtype=float)
    flh=np.ndarray(shape=(nt),dtype=float)
    tmp=np.ndarray(shape=(4),dtype=float)
    for it in range(0,nt):
        tmp[0]=0.0
        tmp[1]=0.0
        tmp[2]=0.0
        tmp[3]=0.0
        for ix in range(0,nx): # for 2-201 200 grids
            if ix!=0 or ix!=nx-1 :
                k=it*(nx*4)+ix+nx*0
                tmp[0]=tmp[0]+onedim[k]*1e3*3600. # convert m/s to mm/hr
                k=it*(nx*4)+ix+nx*1
                tmp[1]=tmp[1]+onedim[k]
                k=it*(nx*4)+ix+nx*2
                tmp[2]=tmp[2]+onedim[k]
                k=it*(nx*4)+ix+nx*3
                tmp[3]=tmp[3]+onedim[k]
        preci[it]=tmp[0]/(nx*1.0)
        pbl[it]=tmp[1]/(nx*1.)
        fsh[it]=tmp[2]/(nx*1.)
        flh[it]=tmp[3]/(nx*1.)
    del onedim
    ntt=nt/(4*24)
    rainsim=np.ndarray(shape=(ntt),dtype=float)
    for it in range(0,nt/(4*24)):
        a=0.
        for i in range(0,96):
            k=it*96+i
            a=a+preci[k]
        rainsim[it]=a*24./96.  #rainfall_dailymean  mm/day
#------------------------------------------------------------------------------
    fpath=dirtrmm+area+"3B42RT_"+yearstr+"_rainfall_dailymean.txt"
    iskp=1
    onedim=readAscii(fpath,iskp)
    rainobs=[]
    for i in range(0,int(len(onedim)/2)):
        rainobs.append(onedim[i*2+1])
    nl=len(rainobs)
    fpath=direra+area+'/'+area+'_'+datestr+'_SHLH_ERA.43'
    iskp=0
    onedim=readAscii(fpath,iskp)
    rainect=np.ndarray(shape=(nl),dtype=float)
    rainecc=np.ndarray(shape=(nl),dtype=float)
    tmperact=[]
    tmperacc=[]
    for it in range(0,len(onedim)/4):
        k=it*4+2
        tmperact.append(onedim[k])
        tmperacc.append(onedim[k+1])
    for it in range(0,len(tmperact)/8):
        a=0.
        for i in range(0,8):
            k=it*8+i
            a=a+tmperact[k]
        rainect[it]=a*24./8.  ### mm/day
    for it in range(0,len(tmperacc)/8):
        a=0.
        for i in range(0,8):
            k=it*8+i
            a=a+tmperacc[k]
        rainecc[it]=a*24./8.  #### mm/day
    datestart=datetime.datetime(iy,im,jd,0,0,0)
    det=datetime.timedelta(hours=24)            
    dateiso=[]            
    for dt in range(0,nl):
        dateiso.append(datestart+dt*det)
    xdate=[]    
    xxx=range(0,nl)            
    for tm in dateiso:
        xdate.append(datetime.datetime.strftime(tm,"%d/%b")) 
#
    fig,ax=plt.subplots(nrows=1,ncols=1,figsize=(14,4))
    ns=min(nl,ntt)
#    ax=plt.subplot(1,1,1)    
    ax.plot(xxx[0:ns],rainsim[0:ns],color=color_cycle[0],lw=wd[0],label="CRM") 
    ax.plot(xxx[0:ns],rainobs[0:ns],color=color_cycle[1],lw=wd[1],label="TRMM")
    ax.plot(xxx[0:ns],rainect[0:ns],color=color_cycle[2],lw=wd[2],label="ERA")
    tilstr=area
    ax.set_title(tilstr, fontsize=18)
    xx1=[]
    yy1=[]
    yy2=[]
    for irn in range(0,len(rainsim)):
        if rainsim[irn]>=0.0:
            xx1.append(rainsim[irn])
            yy1.append(rainobs[irn])
            yy2.append(rainect[irn])
    r1,p1=scista.pearsonr(xx1,yy1)
    r1str="r1= "+"%.2f"%r1
    r2,p2=scista.pearsonr(xx1,yy2)
    r2str="r2= "+"%.2f"%r2
#    plt.ylabel(r'Precipitation ($mm$ $hr^{-1}$)', size=18)
    ax.set_yticks(range(0,25))
    ymajorLocator   = MultipleLocator(4)
    ax.yaxis.set_major_locator(ymajorLocator)
    ax.text(4,20,r1str,fontsize=18)
    ax.text(4,16,r2str,fontsize=18)
    ax.set_xticks(range(0,ns,35))
    xticklabels = [xdate[nn] for nn in range(0,ns,35)] 
    ax.set_xticklabels(xticklabels, size=18)
    fig.text(0.03, 0.85, r'Precipitation ($mm$ $day^{-1}$)', ha = 'left',fontsize=18,rotation=90)
    fig.subplots_adjust(left=0.1,bottom=0.12,right=1-0.05,top=1-0.1,hspace=0.4)
    plt.legend()
    #plt.show()                     
    plt.savefig(dirout+area+'_rainfall_Color_423_2.png',dpi=300)          
    #plt.show()
    plt.close()
    ############---- plot scatter --------------------------------------------------
    fig,ax=plt.subplots(nrows=1,ncols=2,figsize=(12,6))
    plt.subplot(1,2,1)
    ax[0]=plt.subplot(1,2,1)
    ax[0].set_ylim(0,25)
    ax[0].set_xlim(0,25)  
    plt.scatter(rainsim[0:ns],rainobs[0:ns],c='g',marker="+",
            label='TRMM')
    y=rainobs[0:ns]
    x=rainsim[0:ns]
    cf=np.polyfit(x,y,1)
    inter=cf[1]
    slope=cf[0]
    R0=np.corrcoef(x,y)
    RL=R0[0,1]
    R,p1=scista.pearsonr(x,y)
    xline=[0,25]
    yline=[]
    for x in xline:
        yline.append(x*slope+inter)
    ax[0].plot(xline,yline,color='g',lw=2,ls=':')
    term='slope='+'%.1f'%slope+' R='+' %.2f'%R
    ax[0].text(10,22,term,fontsize=18,color='g')  
    plt.scatter(rainsim[0:ns],rainect[0:ns],c=color_cycle[4],marker="o",
            linewidths=0.0,label='ERA',alpha=0.75)
    x=rainsim[0:ns]
    y=rainect[0:ns]
    cf=np.polyfit(x,y,1)
    inter=cf[1]
    slope=cf[0]
    R0=np.corrcoef(x,y)
    RL=R0[0,1]
    R,p1=scista.pearsonr(x,y)
    xline=[0,25]
    yline=[]
    for x in xline:
        yline.append(x*slope+inter)
    ax[0].plot(xline,yline,color=color_cycle[4],lw=2,ls=':')
    term='slope='+'%.1f'%slope+' R='+' %.2f'%R
    ax[0].text(10,20,term,fontsize=18,color=color_cycle[4])
    ymajorLocator   = MultipleLocator(4)
    ax[0].yaxis.set_major_locator(ymajorLocator)
    xmajorLocator   = MultipleLocator(4)
    ax[0].xaxis.set_major_locator(xmajorLocator)
    ######### get 90% of the rainfall ****** minus 5% smalleast and 5% greatest
    plt.subplot(1,2,2)
    ax[1]=plt.subplot(1,2,2)
    rainobs90=getpct(rainobs,90)
    rainsim90=getpct(rainsim,90)
    rainect90=getpct(rainect,90)
    nss=min(len(rainect90),len(rainobs90),len(rainsim90))
    ax[1].set_ylim(0,10)
    ax[1].set_xlim(0,10)  
    plt.scatter(rainsim90[0:nss],rainobs90[0:nss],c='g',marker="+",
            label='TRMM')
    y=rainobs90[0:nss]
    x=rainsim90[0:nss]
    cf=np.polyfit(x,y,1)
    inter=cf[1]
    slope=cf[0]
    R0=np.corrcoef(x,y)
    RL=R0[0,1]
    R,p1=scista.pearsonr(x,y)
    xline=[0,10]
    yline=[]
    for x in xline:
        yline.append(x*slope+inter)
    ax[1].plot(xline,yline,color='g',lw=2,ls='-')
    term='slope='+'%.1f'%slope+' R='+' %.2f'%R
    plt.ylabel(r'TRMM(3B42) and ERA', size=18)
    plt.xlabel(r'Simulation', size=18)
    ###############################################################################
    ax[1].text(3.5,8.5,term,fontsize=18,color='g')
    plt.scatter(rainsim90[0:nss],rainect90[0:nss],c=color_cycle[4],marker="o",
            linewidths=0.0,label='ERA',alpha=0.75)
    x=rainsim90[0:nss]
    y=rainect90[0:nss]
    cf=np.polyfit(x,y,1)
    inter=cf[1]
    slope=cf[0]
    R0=np.corrcoef(x,y)
    RL=R0[0,1]
    R,p1=scista.pearsonr(x,y)
    xline=[0,10]
    yline=[]
    for x in xline:
        yline.append(x*slope+inter)
    ax[1].plot(xline,yline,color=color_cycle[4],lw=2,ls=':')
    term='slope='+'%.1f'%slope+' R='+' %.2f'%R
    ax[1].text(3.5,7.5,term,fontsize=18,color=color_cycle[4])
    ymajorLocator   = MultipleLocator(4)
    ax[1].yaxis.set_major_locator(ymajorLocator)
    xmajorLocator   = MultipleLocator(4)
    ax[1].xaxis.set_major_locator(xmajorLocator)
    plt.legend()
    #plt.show()                     
    plt.savefig(dirout+area+'_rainfall_scatter_color.png',dpi=300)          
    #plt.show()
    plt.close()
    #################################################################
    rainobsmon,stdobs=getmonmean(rainobs,ndays)
    rainsimmon,stdsim=getmonmean(rainsim,ndays)
    rainectmon,stdect=getmonmean(rainect,ndays)
    fig,ax=plt.subplots(nrows=1,ncols=1,figsize=(12,5))
    ns=len(rainobsmon)
#    ax=plt.subplot(1,1,1)    
    ax.errorbar(xxx[0:ns],rainsimmon[0:ns],yerr=stdsim,fmt='-o',color=color_cycle[0],lw=wd[0],label="CRM") 
    ax.errorbar(xxx[0:ns],rainobsmon[0:ns],yerr=stdobs,fmt='-D',color=color_cycle[1],lw=wd[1],label="TRMM")
    ax.errorbar(xxx[0:ns],rainectmon[0:ns],yerr=stdect,fmt='-v',color=color_cycle[2],lw=wd[2],label="ERA")
    tilstr=area
    ax.set_title(tilstr, fontsize=18)
    xx1=[]
    yy1=[]
    yy2=[]
    for irn in range(0,len(rainsimmon)):
        if rainsim[irn]>=0.0:
            xx1.append(rainsimmon[irn])
            yy1.append(rainobsmon[irn])
            yy2.append(rainectmon[irn])
    r1,p1=scista.pearsonr(xx1,yy1)
    r1str="r1= "+"%.2f"%r1
    r2,p2=scista.pearsonr(xx1,yy2)
    r2str="r2= "+"%.2f"%r2
#    plt.ylabel(r'Precipitation ($mm$ $hr^{-1}$)', size=18)
    ax.set_yticks(range(0,13))
    ymajorLocator   = MultipleLocator(3)
    ax.yaxis.set_major_locator(ymajorLocator)
    #ax.text(2,12,r1str,fontsize=18)
    #ax.text(2,11,r2str,fontsize=18)
    ax.set_xticks(range(0,ns,1))
    xticklabels = [monstr[nn] for nn in range(0,ns,1)] 
    ax.set_xticklabels(xticklabels, size=18)
    fig.text(0.03, 0.85, r'Precipitation ($mm$ $day^{-1}$)', ha = 'left',fontsize=18,rotation=90)
    fig.subplots_adjust(left=0.1,bottom=0.12,right=1-0.05,top=1-0.1,hspace=0.4)
    plt.legend()
    #plt.show()                     
    plt.savefig(dirout+area+'monthmean_rainfall_Color_423_2.png',dpi=300)          
    #plt.show()
    plt.close()
    #####################################################################
    rainobssea,std1=getseason(rainobs,ndays)
    rainsimsea,std2=getseason(rainsim,ndays)
    rainectsea,std3=getseason(rainect,ndays)
    ####    
    rainobssea,std1=Mgetseason(rainobsmon)
    rainsimsea,std2=Mgetseason(rainsimmon)
    rainectsea,std3=Mgetseason(rainectmon)
    ind = np.arange(len(rainectsea))  # the x locations for the groups
    width = 0.3  # the width of the bars
    fig, ax = plt.subplots(nrows=1,ncols=1,figsize=(10,8))
    rects1 = ax.bar(ind - 2*width, rainobssea, width, yerr=std1,
                color='lime', label='TRMM')
    rects2 = ax.bar(ind - width, rainsimsea, width, yerr=std2,
                color='deeppink', label='SIM')
    rects3 = ax.bar(ind, rainectsea, width, yerr=std3,
                color='blue', label='ERA')
    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel(r'Precipitation ($mm$ $day^{-1}$)',fontsize=18)
    ax.set_title('Season Rainfall',fontsize=18)
    ax.set_xticks(ind-width)
    ax.set_xticklabels(senstr)
    ax.legend()
    plt.savefig(dirout+area+'Mseasonmean_rainfall_Color_hist.png',dpi=300)          
    #plt.show()
    plt.close()
#### ==========================================================================


