#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 02 20:02:42 2015

@author: chenjh
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import string
from pylab import *
from matplotlib.font_manager import FontProperties
mpl.rcParams['ytick.labelsize'] = 20
mpl.rcParams['xtick.labelsize'] = 20
#====================================================
def readAscii(fpath,iskp):
    #iskp  the total line skipped of the file
    # fpath   the full path of the file
    # usage: onedim=readAscii(fpaht,iskp)
    onedim=[]
    linesplit=[]
    f=open(fpath)
    ff=f.readlines()[iskp:]  ## first line in obs file is legend 
    for line in ff:
        line=string.lstrip(line)
        linesplit.append(line[:-1].split(' '))
    for lnstrs in linesplit:
        for strs in lnstrs:
            if strs!='':
                onedim.append(string.atof(strs))
    del linesplit
    f.close()
    return onedim
#====================================================
nz=34
ng=5
CASENMSTR=['PRD2D0','MLYR2D0','ETP2D0']
CASENM='PRD'
orderstr=[r'($a$)',r'($b$)',r'($c$)',r'($d$)',r'($e$)',r'($f$)']
DATESTR  =['20100101' ]
nga=len(CASENMSTR)
dirin0="/Volumes/MacHD/Work/Papers/WK03/DATA/"
dirpic="/Volumes/MacHD/Work/Papers/WK03/pics/"
monstr=['Jan.','Feb.','March','April','May','June','July','Aug.',
        'Sept.','Oct.','Nov.','Dec.']
nm=len(monstr)
mondays=[31,28,31,30,31,30,31,31,30,31,30,31]
if mondays[1]==28:
    ndays=365
elif mondays[1]==29:
    ndays=366
#-----------------------------------------------------------------------
zdat=[              0.0500000, 0.1643000, 0.3071000, 0.4786000
    , 0.6786000, 0.9071000, 1.1640000, 1.4500000, 1.7640001
    , 2.1070001, 2.4790001, 2.8789999, 3.3069999, 3.7639999
    , 4.2500000, 4.7639999, 5.3070002, 5.8790002, 6.4790001
    , 7.1069999, 7.7639999, 8.4499998, 9.1639996, 9.9069996
    ,10.6800003,11.4799995,12.3100004,13.1599998,14.0500002
    ,14.9600000,15.9099998,16.8799992,17.8799992,18.9099998]
height=[ -50.000 ,    50.000 ,   164.286,    307.143,    478.571  ,  678.571 ,
      907.143 ,  1164.286,   1450.000,   1764.286 ,  2107.143,   2478.572 ,
      2878.572,   3307.143,  3764.286,  4250.000,   4764.286,   5307.143, 
      5878.571,   6478.571,   7107.143,  7764.286,  8450.000,  9164.285,  
      9907.143,  10678.570,  11478.570,  12307.143,  13164.285,  14050.000,
      14964.285,  15907.143,  16878.572,  17878.572,  18907.145,  19964.285,
      21050.000,  22164.285,  23307.145,  24478.572,  25678.572,  26907.145,
      28164.285,  29450.000,  30764.285,  32107.145,  33478.570,  34878.570,
      36307.141,  37764.285,  39250.000,  40750.000]
nz=len(height)
#
#--------------------read data from cloudcell.f ---------------------------
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++            
font = {'family' : 'serif',
        'color'  : 'k',
        'weight' : 'normal',
        'size'   : 20,
        }  
cloudclors=['w','lightgray','plum','darkorchid','darkviolet','b','dodgerblue','skyblue','aqua',
            'greenyellow','lime','yellow','darkorange','chocolate','tomato','r']
daydata=np.ndarray(shape=(nga,4,ndays,nz),dtype=float)
mondata=np.ndarray(shape=(nga,4,nm,nz),dtype=float)
for iga in range(0,nga):    
    casenm=CASENMSTR[iga]
    if casenm[0:3]=='MLY' :
        areastr=casenm[0:4]
    else:
        areastr=casenm[0:3]
    dirin=dirin0+areastr+'/'
    fpath=dirin+casenm+'_daily_qcqaqbqr.txt'
    iskp=0
    onedim=readAscii(fpath,iskp)
    for ke in range(0,ndays):
        for iz in range(0,nz):
            for iv in range(0,4):
                k=ke*(nz*(4+1))+iz*(4+1)+iv+1
                daydata[iga,iv,ke,iz]=onedim[k] #/1000. #qctmp,qatmp,qbtmp,qrtmp
    for iv in range(0,4):        
        for iz in range(0,nz): 
            itt=0           
            for im in range(0,nm):
                nd=mondays[im]
                tmp=0.0
                for it in range(0,nd):
                    tmp=tmp+daydata[iga,iv,itt,iz]
                    itt=itt+1
                mondata[iga,iv,im,iz]=tmp/(nd*1.0)
#============================================================
#====================Ploting=================================
fig,ax=plt.subplots(nrows=3,ncols=2,figsize=(20,12))
xdat=range(0,ndays)
levss=[0,0.5,1.0,1.5,2,2.5,3,4,5,6,8,10]
levs=[]
for a in levss:
    levs.append(a*0.05)
jr=0
jc=0
ij=1
for iga in range(0,nga):
    casenm=CASENMSTR[iga]
    if casenm[0:3]=='MLY' :
        areastr=casenm[0:4]
    else:
        areastr=casenm[0:3]
    if jc==2:
        jc=0
        jr=jr+1
    plt.subplot(3,2,ij)
    ax[jr,jc]=plt.contourf(height,xdat,daydata[iga,0,:,:]+daydata[iga,3,:,:],
        cmap=plt.get_cmap('BuGn'),levels=levs, extend='both')
    marknm=orderstr[ij-1]+areastr+' Liquid Water'
    plt.title(marknm,fontsize=20)    
    plt.axis([0, ndays+1, 0, 16])
    axx=plt.subplot(3,2,ij)
    xmajorLocator   = MultipleLocator(30)
    axx.xaxis.set_major_locator(xmajorLocator) 
    ymajorLocator   = MultipleLocator(4)
    axx.yaxis.set_major_locator(ymajorLocator)
    if jr==2:
        plt.xlabel(r'Day of Year', fontdict=font)
    if jc==0:
        plt.ylabel(r'Height ($km$)', fontdict=font)
    jc=jc+1
    ij=ij+1
    plt.colorbar()
    #*******************************************************
    plt.subplot(3,2,ij)
    ax[jr,jc]=plt.contourf(height,xdat,daydata[iga,1,:,:]+daydata[iga,2,:,:],
        cmap=plt.get_cmap('BuGn'), levels=levs,extend='both')
    marknm=orderstr[ij-1]+areastr+' Ice Water'
    plt.title(marknm,fontsize=20)    
    plt.axis([0, ndays+1, 0, 16])
    axx=plt.subplot(3,2,ij)
    xmajorLocator   = MultipleLocator(30)
    axx.xaxis.set_major_locator(xmajorLocator) 
    ymajorLocator   = MultipleLocator(4)
    axx.yaxis.set_major_locator(ymajorLocator)
    if jr==2:
        plt.xlabel(r'Day of Year', fontdict=font)
    if jc==0:
        plt.ylabel(r'Height ($km$)', fontdict=font)
    jc=jc+1
    ij=ij+1
    plt.colorbar()
    #figtitle = r"Frequency of all cloud cells ($10^{-2}$ %)"
    #fig.text(0.5, 0.95, figtitle,
    #    horizontalalignment='center',
    #    fontproperties=FontProperties(size=20))
#cax = fig.add_axes([0.2, 0.025, 0.6, 0.02])
plt.subplots_adjust(left = 0.08, wspace = 0.2, hspace = 0.25, \
        bottom = 0.1, right=0.88, top = 0.90)
#cax = fig.add_axes([0.90, 0.2, 0.03, 0.6])
#fig.colorbar(ax[0,0],cax,orientation='vertical',extend='both', # 'horizontal'
#    extendfrac='auto',  spacing='uniform')
plt.savefig(dirpic+"AllRNG_Warm&Ice_Profiles_dayly.png",dpi=300)          
plt.close()
############################################################################### 
levss=[0,0.5,1.0,1.5,2,2.5,3,4,5,6,8,10]
levs=[]
for a in levss:
    levs.append(a*0.01)
fig,ax=plt.subplots(nrows=3,ncols=2,figsize=(20,12))
xdat=range(0,nm)
jr=0
jc=0
ij=1
for iga in range(0,nga):
    casenm=CASENMSTR[iga]
    if casenm[0:3]=='MLY' :
        areastr=casenm[0:4]
    else:
        areastr=casenm[0:3]
    if jc==2:
        jc=0
        jr=jr+1
    plt.subplot(3,2,ij)
    ax[jr,jc]=plt.contourf(height,xdat,mondata[iga,0,:,:]+mondata[iga,3,:,:],
        cmap=plt.get_cmap('BuGn'),levels=levs, extend='both')
    marknm=orderstr[ij-1]+areastr+' Liquid Water'
    plt.title(marknm,fontsize=20)    
    plt.axis([0, nm, 0, 16])
    axx=plt.subplot(3,2,ij)
    axx.set_xticks(range(0,nm,1))  
    xticklabels = [monstr[nn] for nn in range(0,nm,1)] 
    axx.set_xticklabels(xticklabels, size=20,rotation=90)
    #axx.xaxis.set_major_locator(xmajorLocator) 
    ymajorLocator   = MultipleLocator(4)
    axx.yaxis.set_major_locator(ymajorLocator)
    if jr==2:
        plt.xlabel(r'Month', fontdict=font)
    if jc==0:
        plt.ylabel(r'Height ($km$)', fontdict=font)
    jc=jc+1
    ij=ij+1
    plt.colorbar()
    #*******************************************************
    plt.subplot(3,2,ij)
    ax[jr,jc]=plt.contourf(height,xdat,mondata[iga,1,:,:]+mondata[iga,2,:,:],
        cmap=plt.get_cmap('BuGn'), levels=levs,extend='both')
    marknm=orderstr[ij-1]+areastr+' Ice Water'
    plt.title(marknm,fontsize=20)    
    plt.axis([0, nm, 0, 16])
    axx=plt.subplot(3,2,ij)
    axx.set_xticks(range(0,nm,1))  
    xticklabels = [monstr[nn] for nn in range(0,nm,1)] 
    axx.set_xticklabels(xticklabels, size=20,rotation=90)
    #axx.xaxis.set_major_locator(xmajorLocator) 
    ymajorLocator   = MultipleLocator(4)
    axx.yaxis.set_major_locator(ymajorLocator)
    if jr==2:
        plt.xlabel(r'Month', fontdict=font)
    if jc==0:
        plt.ylabel(r'Height ($km$)', fontdict=font)
    jc=jc+1
    ij=ij+1
    plt.colorbar()
    #figtitle = r"Frequency of all cloud cells ($10^{-2}$ %)"
    #fig.text(0.5, 0.95, figtitle,
    #    horizontalalignment='center',
    #    fontproperties=FontProperties(size=20))
#cax = fig.add_axes([0.2, 0.025, 0.6, 0.02])
plt.subplots_adjust(left = 0.08, wspace = 0.2, hspace = 0.25, \
        bottom = 0.1, right=0.88, top = 0.90)
#cax = fig.add_axes([0.90, 0.2, 0.03, 0.6])
#fig.colorbar(ax[0,0],cax,orientation='vertical',extend='both', # 'horizontal'
#    extendfrac='auto',  spacing='uniform')
plt.savefig(dirpic+"AllRNG_Warm&Ice_Profiles_monthly.png",dpi=300)          
plt.close()