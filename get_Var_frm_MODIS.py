#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  5 13:36:26 2018



The L3 grid cell indexed (1,1) in the SDS is located at the upper left corner of the map
and corresponds to a grid box with boundaries of 89° to 90°N latitude and 179° to 180°W
longitude
@author: chenjh
"""
from pyhdf.SD import SD,SDC
import string
import calendar
import os
import numpy as np
#------------------------------------------------------------------------------
def getfilenames(file_dir,suff): 
    L=[]   
    for root, dirs, files in os.walk(file_dir):  
        for file in files:  
            if os.path.splitext(file)[1] == suff:  
                L.append(os.path.join(root, file))  
    return L 
def isin(ig,lon,lat,slon,elon,slat,elat):
    pnt=0.
    if slon[ig]<=lon and elon[ig]>lon:
        if slat[ig]<=lat and elat[ig]>lat:
            pnt=1.
    return pnt
def getscal(c,cc):
    ax=1.
    bx=0.
    mi=-999.
    ma=-999.
    de=-999.
    n=len(c)
    for i in range(0,n):
        if cc[i]=='add_offset':
            bx=c[i]
        if cc[i]=='scale_factor':
            ax=c[i]
        if cc[i]=='valid_range':
            mi=c[i][0]
            ma=c[i][1]
        if cc[i]=='_FillValue':
            de=c[i]
    return ax,bx,mi,ma,de
###############################################################################       
dirin='/Volumes/DATA02/MODIS/Daily/MOD08_D3/'
dirout='/Volumes/SamWork/MODIS/MOD08_D3/'
regions=["ETP","WTP","PRD","MLYR","NPC","NEC"]
slon=[  90,   80,   110,  110,  112,  120 ]
elon=[  100,  90,   118,  122,  120,  130 ]
slat=[  27.5,  27.5,  27.5, 27, 34,   43  ]
elat=[  37.5,  37.5,  35,  33,  42,   49  ]
ng=len(regions)
default='-999'
vara=['Cloud_Fraction_Mean','Cloud_Fraction_Day_Mean','Cloud_Fraction_Night_Mean',
      'Cloud_Top_Temperature_Mean','Cloud_Top_Temperature_Day_Mean',
      'Cloud_Top_Temperature_Night_Mean', 
      'Cloud_Water_Path_Liquid_Mean','Cloud_Water_Path_Ice_Mean',
      'Cloud_Optical_Thickness_Liquid_Mean','Cloud_Optical_Thickness_Ice_Mean',
      'Cloud_Optical_Thickness_Combined_Mean']
scal=[0.0001,0.0001,0.0001,
       0.01,0.01,0.01,
       1.,1.,
       0.01,0.01,0.01 ]
addoff=[0,0,0,
        -15000,-15000,-15000,
        0.,0.,
        0.,0.,0. ]
nv=len(vara)
for iy in range(2010,2011):
    yearstr='%i'%iy
    if calendar.monthrange(iy,2)[1]==28:
        nday=365
    else:
        nday=366
    afile=[]
    for ig in range(0,ng):
        filename=regions[ig]+'_MODIS_dayliy_mean.txt'
        fpath=dirout+filename
        f=open(fpath,'w')
        item='Date'
        for va in vara:
            item=item+' '+va
        f.write(item)
        f.write('\n')
        afile.append(f)
        del f
    for idy in range(0,nday):
        daystr='%.3i'%(idy+1)
        datestr=yearstr+daystr
        tmdir=dirin+yearstr+'/'+daystr+'/'
        files=getfilenames(tmdir,'.hdf')
        cont=np.zeros(shape=(ng,nv),dtype=float)
        sums=np.zeros(shape=(ng,nv),dtype=float)
        for fpath in files:
            hdf=SD(fpath, SDC.READ)
            for iv in range(0,nv):
                tmvar=vara[iv]
                rawdata=hdf.select(tmvar)
                nlon=len(rawdata[0,:])
                nlat=len(rawdata[:,0])
                for ix in range(0,nlon):
                    lon=ix*1.0
                    for iy in range(0,nlat):
                        lat=90-iy*1.0
                        for ig in range(0,ng):
                            pnt=isin(ig,lon,lat,slon,elon,slat,elon)
                            if pnt>0:
                                c=rawdata.attributes().values()
                                cc=rawdata.attributes().keys()
                                ax,bx,mi,ma,de=getscal(c,cc)
                                #print ax,bx
                                if mi!=-999. :   
                                    if rawdata[iy,ix]>mi and rawdata[iy,ix]<ma:
                                        sums[ig,iv]=sums[ig,iv]+ax*(rawdata[iy,ix]-bx)
                                        cont[ig,iv]=cont[ig,iv]+1.0
                                else:
                                    if rawdata[iy,ix]!=de :
                                        sums[ig,iv]=sums[ig,iv]+ax*(rawdata[iy,ix]-bx)
                                        cont[ig,iv]=cont[ig,iv]+1.0                                    
                del rawdata
            del hdf
        for ig in range(0,ng):
            item=datestr
            for iv in range(0,nv):
                if cont[ig,iv]>0:
                    item=item+' %f'%(sums[ig,iv]/cont[ig,iv])
                else:
                    item=item+' '+default
            afile[ig].write(item)
            afile[ig].write('\n')
    for ig in range(0,ng):
        afile[ig].close()
                        
                        
                        
                
        
        
        
        
        
        
        
        
        
        