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
nz=34
ng=5
CASENMSTR=['PRD2D0']
CASENM='PRD'
orderstr=[r'($a$)',r'($b$)',r'($c$)',r'($d$)',r'($e$)',r'($f$)']
DATESTR  =['20100101' ]
nga=len(CASENMSTR)
dirin="/Volumes/MacHD/Work/Papers/WK03/DATA/"+CASENM+'/'
dirpic="/Volumes/MacHD/Work/Papers/WK03/pics/"
season=['Spring','Summer','Autumn','Winter']
nsea=len(season)
#-----------------------------------------------------------------------
zdat=[              0.0500000, 0.1643000, 0.3071000, 0.4786000
    , 0.6786000, 0.9071000, 1.1640000, 1.4500000, 1.7640001
    , 2.1070001, 2.4790001, 2.8789999, 3.3069999, 3.7639999
    , 4.2500000, 4.7639999, 5.3070002, 5.8790002, 6.4790001
    , 7.1069999, 7.7639999, 8.4499998, 9.1639996, 9.9069996
    ,10.6800003,11.4799995,12.3100004,13.1599998,14.0500002
    ,14.9600000,15.9099998,16.8799992,17.8799992,18.9099998]
#
def readAscii(fpath,iskp,nrl):
    #iskp  the total line skipped of the file
    # fpath   the full path of the file
    # usage: onedim=readAscii(fpaht,iskp)
    onedim=[]
    linesplit=[]
    f=open(fpath)
    print iskp,nrl
    ff=f.readlines()[iskp:nrl]  ## first line in obs file is legend 
    for line in ff:
        line=string.lstrip(line)
        linesplit.append(line[:-1].split(' '))
    for lnstrs in linesplit:
        for strs in lnstrs:
            if strs!='':
                onedim.append(string.atof(strs))
    del linesplit,ff
    f.close()
    print len(onedim)
    return onedim
#--------------------read data from cloudcell.f ---------------------------
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            
#plot  ----------  ft  # xdat,ydat,zdat
font = {'family' : 'serif',
        'color'  : 'k',
        'weight' : 'normal',
        'size'   : 20,
        }  
cloudlevs=[5,10,15,20,30,40,50,60,70,80,90,100,110,120,130]
cloudclors=['w','lightgray','plum','darkorchid','darkviolet','b','dodgerblue','skyblue','aqua',
            'greenyellow','lime','yellow','darkorange','chocolate','tomato','r']
for iga in range(0,nga):
    casenm=CASENMSTR[iga]
    if casenm[0:3]=='MLY' :
        areastr=casenm[0:4]
    else:
        areastr=casenm[0:3]
    fig,ax=plt.subplots(nrows=2,ncols=2,figsize=(12,12))
    jr=0
    jc=0
    ij=1
    for ise in range(0,nsea):
        if jc==2:
            jc=0
            jr=jr+1
        ft=np.ndarray(shape=(nz,nz,ng),dtype=float)
        fpath=dirin+casenm+'_ALLCLOUDCELSS_FREQUENCY_F90'+season[ise]+'.TXT'
        for i in range(0,ng):
            iskp=i*(nz+2)+1
            nrl=nz+iskp
            onedim=readAscii(fpath,iskp,nrl)
            for ke in range(0,nz):
                for kb in range(0,nz):
                    k=ke*(nz+1)+kb+1
                    ft[kb,ke,i]=onedim[k]
        plt.subplot(2,2,ij)
#zdat[0,:]=0.0   ## the first level is below surface ground
        ft0=np.ndarray(shape=(nz,nz), dtype=float) #(km,km)  For exchange the dims
        for i1 in range(0,nz):
            i10=nz-i1-1
            for i2 in range(0,nz):
                ft0[i10,i2]=ft[i2,i1,4]
        ax[jr,jc]=plt.contourf(zdat,zdat,ft0,colors=cloudclors, levels=cloudlevs,extend='both')
        marknm=orderstr[ise]+' '+season[ise]
        plt.title(marknm,fontsize=20)    
        plt.axis([0, 20, 0, 20])
        axx=plt.subplot(2,2,ij)
        xmajorLocator   = MultipleLocator(4)
        axx.xaxis.set_major_locator(xmajorLocator) 
        ymajorLocator   = MultipleLocator(4)
        axx.yaxis.set_major_locator(ymajorLocator)
        if jr==1:
            plt.xlabel(r'Cloud Base Height ($km$)', fontdict=font)
        if jc==0:
            plt.ylabel(r'Cloud Top Height ($km$)', fontdict=font)
        jc=jc+1
        ij=ij+1
    figtitle = r"Frequency of all cloud cells ($10^{-2}$ %)"
    fig.text(0.5, 0.95, figtitle,
        horizontalalignment='center',
        fontproperties=FontProperties(size=20))
#cax = fig.add_axes([0.2, 0.025, 0.6, 0.02])
    plt.subplots_adjust(left = 0.08, wspace = 0.3, hspace = 0.25, \
        bottom = 0.1, right=0.88, top = 0.90)
    cax = fig.add_axes([0.90, 0.2, 0.03, 0.6])
    fig.colorbar(ax[0,0],cax,orientation='vertical',extend='both', # 'horizontal'
        extendfrac='auto',  spacing='uniform')
    #plt.show()
    plt.savefig(dirpic+areastr+"_CloudCellsTopBase_fortran_color.png",dpi=300)          
    #plt.show()
    plt.close()
############################################################################### 