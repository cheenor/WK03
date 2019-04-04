#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  5 09:01:29 2018

@author: chenjh
"""
from pyhdf.SD import SD,SDC,SDS
import os
import pprint
def getfilenames(file_dir,suff): 
    L=[]   
    for root, dirs, files in os.walk(file_dir):  
        for file in files:  
            if os.path.splitext(file)[1] == suff:  
                L.append(os.path.join(root, file))  
    return L 
dirin='/Volumes/DATA02/MODIS/Daily/MOD08_D3/'
dirout='/Volumes/SamWork/MODIS/MOD08_D3/'
fpath=dirin+'2010/001/'+'MOD08_D3.A2010001.061.2017308151839.hdf'
hdf = SD(fpath, SDC.READ)
hdf2 = SDS(fpath, SDC.READ)
#fpath=dirout+'Vars_in_file.txt'
#f=open(fpath,'w')
#hdf.datasets()
'''
for n in range(0,len(hdf.datasets())):
    b=hdf.datasets().values()[n]
    a=hdf.datasets().keys()[n]
    btr=''
    for elm in b:
        #print elm
        try:
            nx=len(elm)
            for nn in range (0,nx):
                try:
                    btr=btr+' %i'%elm[nn] 
                except:
                    btr=btr+' '+elm[nn]
        except:
            try:
                btr=btr+' %i'%elm 
            except:
                btr=btr+' '+elm          
    f.write('%i '%n+a+' '+btr)
    f.write('\n')
'''
data3D = hdf.select('Cloud_Fraction_Mean')
pprint.pprint(data3D.attributes())
#data = data3D[11,:,:]
files=getfilenames(dirin,'.hdf')