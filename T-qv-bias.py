#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 31 20:17:12 2018

@author: jhchen
"""
import matplotlib.pyplot as plt
import matplotlib as mpl
import string
import numpy as np
import datetime
import struct
mpl.rcParams['ytick.labelsize'] = 20
mpl.rcParams['xtick.labelsize'] = 20
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
def readBin(fpath,N):
    # fpath   the full path of the file
    # usage: onedim=readAscii(fpaht,iskp)
    onedim=[]
    linesplit=[]
    f=open(fpath,'rb')
    data=struct.unpack('f',f.read(N:N+4))
    f.close()
    return onedim
def heigh2pressure(yy,tmp4prs,ydat,prelevel):
    nt=len(tmp4prs[0,:])
    nz=len(tmp4prs[:,0])
    menatmp=np.zeros(shape=(nz),dtype=float)
    for iz in range(0,nz):
        for it in range(0,nt):
            menatmp[iz]=menatmp[iz]+tmp4prs[iz,it]/nt    
    yyr=[]
    for h in yy:    
        if h==-50.:
            p=prelevel[0]
            yyr.append(p)
            #print p
            #return p
        else:
            for iz in range(1,nz):
                if h>ydat[iz-1] and h<=ydat[iz]:
                    dz=h-ydat[iz-1]
                    p0=prelevel[iz-1]
                    tmp=menatmp[iz]+(dz/(ydat[iz]-ydat[iz-1]))*(menatmp[iz]-menatmp[iz-1])
                    tmp2=dz/(18400*(1+(tmp-273.15)/273))
                    p=p0/(10**tmp2)
                    #print p
                    #return p
                    yyr.append(p)    
    return yyr
def pressure2heigh(pres,tmp4prs,ydat,prelevel):
    nt=len(tmp4prs[0,:])
    nz=len(tmp4prs[:,0])
    menatmp=np.zeros(shape=(nz),dtype=float)
    for iz in range(0,nz):
        for it in range(0,nt):
            menatmp[iz]=menatmp[iz]+tmp4prs[iz,it]/nt    
    for iz in range(1,nz):
        if pres<prelevel[iz-1] and pres>=prelevel[iz]:
            z1=ydat[iz-1]*1000. # km to m
            p1=prelevel[iz-1]
            at=((menatmp[iz]+menatmp[iz-1])/2.-273.15)*1./273.
            z2=z1+18400*(1+at)*math.log10(p1/pres)
            return z2
def convert_ax(ax_f,tmp4prs,ydat,prelevel):
    """
    Update second axis according with first axis.
    """
    y1, y2 = ax_f.get_ylim()
    yy=ax_f.get_yticks()
    #print yy
    #print y1,y2
    yyr=heigh2pressure(yy,tmp4prs,ydat,prelevel)
    #print yyr
    yyrstr=[]
    for y0 in yyr:
        yyrstr.append('%d'%y0)
    #n=len(yyrstr)
    #yyrstr[n-1]=yyrstr[n-1]+r'$(hPa)$'
    ax_c.set_yticklabels(yyrstr) #,heigh2pressure(y2,tmp4prs,ydat,prelevel))
    ax_c.figure.canvas.draw()            
#########################################################################
casenm='ETPCTR_EC'
nt=355*4  # 355 days
nz=33
starid=4  # this parameter can discard the first day
dirs='D:/Work/MyPaper/WK03/Data/2010ETP/'
diro='D:/Work/MyPaper/2DForcing/Forcing/ETP/'
dircnrain='D:/MyPaper/PhD04/Data/RainCN05/'
pic_out='D:/Work/MyPaper/WK03/Pics/'
dirpic=pic_out
dis_pressure_ea=[750.,550.,400.,250.,150.]
dis_pressure_tp=[450.,300.,200.,100.]
#
ydat_r=[ -50.000 ,    50.000 ,   164.286,    307.143,    478.571  ,  678.571 ,
      907.143 ,  1164.286,   1450.000,   1764.286 ,  2107.143,   2478.572 ,
      2878.572,   3307.143,  3764.286,  4250.000,   4764.286,   5307.143, 
      5878.571,   6478.571,   7107.143,  7764.286,  8450.000,  9164.285,  
      9907.143,  10678.570,  11478.570,  12307.143,  13164.285,  14050.000,
      14964.285,  15907.143,  16878.572,  17878.572,  18907.145,  19964.285,
      21050.000,  22164.285,  23307.145,  24478.572,  25678.572,  26907.145,
      28164.285,  29450.000,  30764.285,  32107.145,  33478.570,  34878.570,
      36307.141,  37764.285,  39250.000,  40750.000]
ydat=[]
for yd in ydat_r:
    ydat.append(yd*0.001)
############################################################################
for iga in range(0,nga):
    casenm=CASE[iga]
    if casenm[0:3]=='ETP' :
        iy,im,jd=2010,6,3
#        nt,nday=125,31
        namestr=casenm[0:3]
        marktr=r"($f$)"
        yearstr="%d"%iy
    elif casenm[0:3]=='WTP' :
        iy,im,jd=2010,7,3
#        nt,nday=125,31
        namestr=casenm[0:3]
        marktr=r"($e$)"
        yearstr="%d"%iy
    elif casenm[0:3]=='PRD' :
        iy,im,jd=2012,4,1
#        nt,nday=125,31
        namestr=casenm[0:3]
        marktr=r"($a$)"
        yearstr="%d"%iy
    elif casenm[0:3]=='MLY' :
        iy,im,jd=2010,6,2
#        nt,nday=125,31
        namestr=casenm[0:4]
        marktr=r"($b$)"
        yearstr="%d"%iy
    elif casenm[0:3]=='NPC' :
        iy,im,jd=2010,8,2
#        nt,nday=125,31
        namestr=casenm[0:3]
        marktr=r"($c$)"
        yearstr="%d"%iy
    elif casenm[0:3]=='NEC' :
        iy,im,jd=2012,7,6
#        nt,nday=125,31
        namestr=casenm[0:3]
        marktr=r"($d$)"
        yearstr="%d"%iy
    topstr='' #'_250'
    datestr="%4d"%iy+"%2.2d"%im+"%2.2d"%jd+"_365d"
    strnm=namestr+'_'+datestr
    dirin=dirin1+namestr+'/'
    folds="/CTREC"+"%4d"%iy+"%2.2d"%im+"%2.2d"%jd+"/Simulation/"
    datestr="%4d"%iy+"%2.2d"%im+"%2.2d"%jd+"_031d"
    dirin=dirs+namestr+folds
    dirobs=diro+namestr+'/'
    f43=namestr+'_'+datestr+"_ERA.43"
    f52=namestr+'_'+datestr+"_ERA_52pressure.52"    
    nameforcing=namestr+'_'+datestr+"_LSFORCING_ERA.37"
    titlestr.append(marktr+' '+namestr+' ('+yearstr+')')
#fpath=dirin+'micro_202_ETP2D3'
#fff=np.fromfile(fpath,dtype=float)
#fc=open(fpath,"rb")
#fff=fc.read()
#print fff[0],fff[1],fff[12],fff[78]
#print len(fff)
    datestart=datetime.datetime(iy,im,jd,0,0,0)
    det=datetime.timedelta(hours=6)            
    dateiso=[]            
    for dt in range(0,nt):
        dateiso.append(datestart+dt*det)
    xdate=[]    
    xdat=range(0,nt)            
    for tm in dateiso:
        xdate.append(datetime.datetime.strftime(tm,"%d/%b")) 
    xdatet.append(xdate)




