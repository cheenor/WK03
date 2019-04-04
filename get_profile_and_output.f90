program digest_and_output
implicit none
integer nz,nx
integer i,j,k,iz,itt,ix,iday
real ntmp(52)
character*150 dirin,dirout,fpathin,fpathout
character*50 filename,regin(3)
real qc(200,52),qa(200,52),qb(200,52),qr(200,52)
real qctmp(52),qatmp(52),qbtmp(52),qrtmp(52)
nz=52
nx=200
dirin='/Volumes/DATA1/Models/CRM/YEAR/'
dirout='/Volumes/MacHD/Work/Papers/WK03/DATA/'
regin(1)='ETP'
regin(2)='PRD'
regin(3)='MLYR'
do i =1,3
    filename=trim(regin(i))//'2D0_Raw_qcqaqbqr.txt'
    fpathin=trim(dirin)//trim(regin(i))//'/'//trim(filename)
    filename=trim(regin(i))//'2D0_daily_qcqaqbqr.txt'
    fpathout=trim(dirout)//trim(regin(i))//'/'//trim(filename)
    open(10,file=trim(fpathin))
    open(20,file=trim(fpathout))
    itt=0           
    qctmp=0.0
    qatmp=0.0
    qbtmp=0.0
    qrtmp=0.0
    iday=1
    ntmp=0.0
    do j=1,73*480 ! one year,every 15 mins
        itt=itt+1
        do ix=1,nx
            read(10,99) qc(ix,:), qa(ix,:) &
        &         , qb(ix,:), qr(ix,:)
        enddo
        do iz =1,nz
            do ix=1,nx
                if((qc(ix,iz)+qc(ix,iz)+qc(ix,iz)+qc(ix,iz))>=0.001)then
                    qctmp(iz)=qctmp(iz)+qc(ix,iz)
                    qatmp(iz)=qatmp(iz)+qa(ix,iz)
                    qbtmp(iz)=qbtmp(iz)+qb(ix,iz)
                    qrtmp(iz)=qrtmp(iz)+qr(ix,iz)
                    ntmp(iz)=ntmp(iz)+1.0
                endif
            enddo
        enddo
        if (itt==96)then
            do iz=1,nz
                if(ntmp(iz)>0)then
                qctmp(iz)=qctmp(iz)/ntmp(iz)
                qatmp(iz)=qatmp(iz)/ntmp(iz)
                qbtmp(iz)=qbtmp(iz)/ntmp(iz)
                qrtmp(iz)=qrtmp(iz)/ntmp(iz)
                else
                qctmp(iz)=0.0
                qatmp(iz)=0.0
                qbtmp(iz)=0.0
                qrtmp(iz)=0.0
                endif
                write(20,*)iday,qctmp(iz),qatmp(iz),qbtmp(iz),qrtmp(iz)         
                qctmp(iz)=0.0
                qatmp(iz)=0.0
                qbtmp(iz)=0.0
                qrtmp(iz)=0.0
                ntmp(iz)=0.0
            enddo
            iday=iday+1
            itt=0 
        endif
    enddo
enddo
99    format(8e12.4)
end program