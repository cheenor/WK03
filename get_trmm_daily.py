#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 10 08:51:16 2018

@author: chenjh
"""
from pytrmm import TRMM3B42RTFile as TRM3B42
import os
import datetime
import calendar
import numpy as np
os.system("cls")
def getcont(a,crilt,crihv):
    n1=len(crilt)
    n2=len(crihv)
    contmp=np.zeros(shape=(n1+n2),dtype=float)
    for i in range(0,n1+n2):
        if i<n1:
            if i==0:
                if a<=crilt[i]:
                    contmp[i]=1
            elif a<crilt[i] and a>=crilt[i-1]:
                contmp[i]=1
        elif i==n1+n2-1 and a>=crihv[i-n1]:
            contmp[i]=1
        else:
            if a>=crihv[i-1-n1] and a < crihv[i-n1]:
                contmp[i]=1
    contmp[0]=0
    contmp[n1+n2-1]=0
    if a<=0.1:
        contmp[0]=1
    if a>=10:
        contmp[n1+n2-1]=1        
    return contmp
def getsas(tmprain):
    n=len(tmprain)
    h=int(n*0.20)
    if h>=1:
        mean=0
        for rain in tmprain:
            mean=mean+rain
        mean=mean/(n*1.0)
        print tmprain
        print 'AAAAA'
        tmprain.sort()
        print tmprain
        mi=0
        for i in range(0,h):
            mi=mi+tmprain[i]/(h*1.0)
        mx=0
        for i in range(n-h-1,n):
            mx=mx+tmprain[i]/(h*1.0)
    else:
        mean=-9
        mi=-9
        mx=-9
    return mean, mi,mx        
#################################################
dirin='/Volumes/DATA1/TRMM/'
dirout="/Volumes/MacHD/Work/Papers/WK03/DATA/"
#regions=['ETP','MLYR','PRD']
#slon =[90]
#elon =[100]
#slat =[27.5]
#elat =[37.5]
regions=["ETP","WTP","PRD","MLYR","NPC","NEC"]
slon=[  90,   80,   110,  110,  112,  120 ]
elon=[  100,  90,   118,  122,  120,  130 ]
slat=[  27.5,  27.5,  27.5, 27, 34,   43  ]
elat=[  37.5,  37.5,  35,  33,  42,   49  ]
# the start date of ever region and how many days for every region.
ayear=[2010]
amon =[1] 
aday =[1]
ndays=[365]
nag=len(regions)
#get the basic information of the data TRMM 3B40
fpath="/Volumes/DATA1/TRMM/2006/01/3B42RT.2006010100.7R2.bin.gz"
trmm_file = TRM3B42(fpath)
fpath=dirout+'3b42_FileHead.txt'
fout=open(fpath,'w')
#print(trmm_file.header())
ax=trmm_file.header()
headnm=ax.keys()
headva=ax.values()
na=len(headnm)
for a in range(0,na) :
    item="%s: "%headnm[a]
    fout.write(item)
    item="%s "%headva[a]
    fout.write(item)
    fout.write('\n')
precip = trmm_file.precip()
print('Array dimensions:', precip.shape)
fout.write('Array dimensions:')
nxy=[]
for a in precip.shape:
    nxy.append(a)
    item="%d "%a
    fout.write(item)
del trmm_file,precip
fout.close()
ny=nxy[0]
nx=nxy[1]
lon=[]
lat=[]
for ix in range(0,nx):
    lon.append(0.125+ix*0.25)
for iy in range(0,ny):
    lat.append(59.875-iy*0.25)
#
def getdatestr(yy,mm,dd,nd):
    datestart=datetime.datetime(yy,mm,dd,0,0,0)
    det=datetime.timedelta(hours=3)            
    dateiso=[]
    nt=nd*8          
    for dt in range(0,nt):
        dateiso.append(datestart+dt*det)
    xdate=[]               
    for tm in dateiso:
        xdate.append(datetime.datetime.strftime(tm,"%b/%d %H:%M"))
    return xdate
ndays=np.zeros(shape=(12),dtype=int)
for i in range(0,12):
    ndays[i]=31
ndays[1]=28
ndays[3]=30
ndays[5]=30
ndays[8]=30
ndays[10]=30
ig=3
#######
littlerain_L1=9.9
littlerain_L2=4.9
littlerain_L3=2.5
#######################
heavyrain_L1=25.
heavyrain_L2=50.
heavyrain_L3=175.
crilitt=[0.5,1,2]
criheav=[10,15,20]
crilevs=[0.01,0.1,0.5,3,6,9]
nl=len(crilitt)+len(criheav)
#####
for iyr in range(2010,2011):    
    #level 0 1  
    yearstr="%d"%iyr  #ayear[ig]
    #################
    fileout=regions[ig]+"3B42RT_"+yearstr+"_rainfall_dailymean.txt"
    fpath=dirout+fileout
    fout1=open(fpath,'w')
    item='Date mean'
    fout1.write(item)
    fout1.write('\n')
    ndays[1]=28
    if iyr%400==0 or (iyr%4==0 and iyr%100!=0):
        ndays[1]=29
    else:
        ndays[1]=28
    for im in range(0,12):
        monstr="%2.2d"%(im+1)
        nd=ndays[im]
        #tmprain=[]
        for day in range(0,nd):
            daystr="%2.2d"%(day+1)
            frq=np.zeros(shape=(nl),dtype=float)
            frqcont=np.zeros(shape=(nl),dtype=float)
            tmprain=[]
            cont=0.
            tmp=0. 
            for stp in (0,8):
                stpstr='%2.2d'%(stp*3)
                filename1="3B42RT."+yearstr+monstr+daystr+stpstr+".7R2.bin.gz"
                filename2="3B42RT."+yearstr+monstr+daystr+stpstr+".7.bin.gz"
                if iyr <2014:
                    fpath1=dirin+yearstr+'/'+monstr+'/'+filename1
                    fpath2=dirin+yearstr+'/'+monstr+'/'+filename2
                else:
                    fpath1=dirin+yearstr+'/'+filename1
                    fpath2=dirin+yearstr+'/'+filename2
                fpath=''
                if os.path.isfile(fpath1):
                    fpath=fpath1
                if os.path.isfile(fpath2):
                    fpath=fpath2
                if len(fpath)>0:
                    trmm_file=TRM3B42(fpath)
                    precip = trmm_file.precip()                   
                    for ix in range(0,nx):
                        if lon[ix]>(slon[ig]) and lon[ix]<(elon[ig]):
                            for iy in range(0,ny):
                                if lat[iy]>slat[ig] and lat[iy]<elat[ig]:
                                    if precip[iy,ix]>=0.:
                                        #print precip[iy,ix]
                                        tmp=tmp+precip[iy,ix]*24.
                                        #tmprain.append(precip[iy,ix])
                                        cont=cont+1.
                                        #contmp=getcont(precip[iy,ix],crilitt,criheav)
                                        #for ifq in range(0,nl):
                                        #    frq[ifq]=frq[ifq]+contmp[ifq]
            if cont>0:  
                item=yearstr+monstr+daystr+' '                
                item=item+'%.4f '%(tmp/cont)
            else:
                item=yearstr+monstr+daystr+' '
                item=item+'-9 '
            #fout.write(item)
            del trmm_file,precip
            #else:
            #    fout.write("-9")  # there is no data file
            #fout.write("\n")
            #if len(tmprain)>0:
            #    mean,mi,mx=getsas(tmprain)
            #else:
            #    mean=-9,
            #    mi=-9
            #    mx=-9
            #item=yearstr+monstr+daystr+' '+'%.2f '%mean+'%.2f '%mi+'%.2f '%mx
            fout1.write(item)
            fout1.write('\n')
    fout1.close()
'''   
    ndaystr="%3.3d"%ndays[ig]
    fpath=dirout+regions[ig]+"-"+yearstr+monstr+daystr+"_"+ndaystr+"d_TRMM3B40.txt"
    fout=open(fpath,"w")
    fout.write("Date ")
    fout.write("Precipitation")
    fout.write("\n")
    ndd=ndays[ig]
    xdate=getdatestr(ayear[ig],amon[ig],aday[ig],ndays[ig])
    iyy=ayear[ig]    
    imm=amon[ig]
    jd=aday[ig]
    monday = calendar.monthrange(iyy,imm)
    idt=0
    for idd in range(0,ndd):
        jdd=jd+idd
        if jdd>monday[1]:
            imm=imm+1
            monday = calendar.monthrange(iyy,imm)
            jdd=1;jd=1-idd
        yearstr="%d"%iyy
        monstr="%2.2d"%imm
        daystr="%2.2d"%jdd
        for ih in range(0,24,3):
            fout.write(xdate[idt]+" ")
            idt=idt+1
            hourstr="%2.2d"%ih
            filename="3B40RT."+yearstr+monstr+daystr+hourstr+".7R2.bin.gz"
            fold=yearstr+"/"+monstr+"/"
            fpath=dirin+fold+filename
            print fpath
            if os.path.isfile(fpath):
                trmm_file=TRM3B42(fpath)
                precip = trmm_file.precip()
                cont=0.
                tmp=0.
                for ix in range(0,nx):
                    if lon[ix]>(slon[ig]) and lon[ix]<(elon[ig]):
                        for iy in range(0,ny):
                            if lat[iy]>slat[ig] and lat[iy]<elat[ig]:
                                if precip[iy,ix]>=0.0 :
                                    tmp=tmp+precip[iy,ix]
                                    cont=cont+1.
                if cont>0:
                    item="%f "%(tmp/cont)
                else:
                    item='-9'
                fout.write(item)
                del trmm_file,precip
            else:
                fout.write("-9")  # there is no data file
            fout.write("\n")
    fout.close()
'''            
