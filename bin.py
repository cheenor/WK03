#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 06 16:35:37 2018

@author: jhchen
"""
import matplotlib.pyplot as plt
import matplotlib as mpl
import string
import numpy as np
import datetime
import struct
casenm='ETPCTR_EC'
nt=365*4  # 365 days
nz=33
starid=4  # this parameter can discard the first day
dirbin='/Work/MyPaper/2DForcing/Forcing/ETP/'
dircnrain='D:/MyPaper/PhD04/Data/RainCN05/'
pic_out='D:/Work/MyPaper/WK03/Pics/'
filename='ETP2D0_Raw_qcqaqbqr.bin'
fpath=dirs+filename
f=open(fpath,'rb')
data=struct.unpack('f',f.read(4))
print data
data1=f.read()[1:2]
print ord(data1)