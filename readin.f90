!Read sphere information then generate grids
subroutine readin
use molecule
use comdata
use treecode
implicit double precision(a-h,o-z)
real*8 pos(3),vector(3)
integer nind(5)
! temp. var. for read in
CHARACTER(100) :: FHEAD
character(10) :: c1,c3,c4,c5,c6,c8
integer i2,i7,i5, max_ires, max_iatm, idx(3)
integer i,j,nremark,MEOF, nremark1, nremark2
real*8 xyzqr(5)

!------------------------
!added for multiple domain scheme 05/22/16
real*8 tvec(3),angle(3), zero_vec(3), center1(3), center2(3), diff(3), dist
real*8,allocatable:: sptpos_copy(:,:),sptnrm_copy(:,:),sptpos_copy2(:,:),sptnrm_copy2(:,:)
real*8,allocatable:: atmchr_copy(:),atmpos_copy(:,:),atmrad_copy(:),chrpos_copy(:,:)
integer,allocatable:: nvert_copy(:,:), nvert_copy2(:,:)

!added by gan!: for icosahedral grid
real*8 rad, center(3), xxx,yyy,zzz,rrr,aa,bb,cc
integer init_type,level,iwind,nedge,err,ierr
integer, allocatable:: nefv(:,:)

!gan: first define the level of deviding

read(den,'(i1)') level
!write(*,*) level

!imd=0 !0: monomer; 1: dimer by rot. & trans.; 2:dimer by reading cood.


!Obtain path
pathname='2sphere_data/'
lenpath = len(pathname)
do while (pathname(lenpath:lenpath) .eq. ' ')
    lenpath = lenpath - 1
enddo  

!Obtain filename
lenfname = len(fname)
do while (fname(lenfname:lenfname) .eq. ' ')
    lenfname = lenfname - 1
enddo 



if (imd==0) then
    nremark=0
    open(102,file=pathname(1:lenpath)//fname(1:lenfname)//".pqr")
    do
        READ(102,*) fhead
        if (fhead(1:6)=='REMARK') then
            nremark=nremark+1
        else
            exit
        endif
    enddo
    print *,'lines of remarks = ', nremark
    close(102)

    open(102,file=pathname(1:lenpath)//fname(1:lenfname)//".pqr")
    open(103,file=pathname(1:lenpath)//fname(1:lenfname)//".xyzr")

    do i=1,nremark
        read(102,*) fhead
    enddo

    natm=0
    do
        read(102,*,IOSTAT = MEOF) c1,i2,c3,c4,i5,xyzqr 

        if ((c1(1:3) .ne. 'END') .and. (MEOF .eq. 0)) then
            write(103,*) xyzqr(1:3),xyzqr(5)
            natm=natm+1
        endif
        IF(MEOF .LT. 0) EXIT
    enddo

    close(102)
    close(103)
    print *,'number of atoms = ', natm

    !Read atom coordinates and partial charges
    nchr=natm
    allocate(atmpos(3,natm),atmrad(natm),atmchr(nchr),chrpos(3,nchr),STAT=ierr)
    IF (ierr .NE. 0) THEN
        WRITE(6,*) 'Error allocating atmpos, atmrad, atmchr, chrpos!'
        STOP
    END IF

    open(102,file=pathname(1:lenpath)//fname(1:lenfname)//".pqr")

    do i=1,nremark
        read(102,*) fhead
    enddo

    do i=1,natm
        read(102,*,IOSTAT = MEOF) c1,i2,c3,c4,i5,xyzqr 
        atmpos(1:3,i)=xyzqr(1:3)
        atmrad(i)=xyzqr(5)
        chrpos(1:3,i)=xyzqr(1:3)
        atmchr(i)=xyzqr(4)
    enddo
    close(102)

    rslt=system('msms -if '//pathname(1:lenpath)//fname(1:lenfname)//".xyzr"//' -prob 1.4 -de ' &
    //den(1:5)//' -of '//pathname(1:lenpath)//fname(1:lenfname)//' > msms.out')    
      ! read the surface points
    OPEN(102,FILE=pathname(1:lenpath)//FNAME(1:lenfname)//".vert")

    READ(102,*) FHEAD
    READ(102,*) FHEAD
    READ(102,*) NSPT, ppp, qqq, rrr

    ALLOCATE(SPTPOS(3,NSPT), SPTNRM(3,NSPT), NATMAFF(NSPT), NSFTYPE(NSPT), STAT= ierr)
    IF (ierr .NE. 0) THEN
        WRITE(6,*) 'Error allocating SPTPOS, SPTNRM, NATMAFF, NSFTYPE!'
    STOP
    END IF

    SPTPOS=0.D0; SPTNRM=0.D0; NATMAFF=0; NSFTYPE=0;	
      
    DO I=1,NSPT
        READ(102,*) POS(1:3), VECTOR(1:3), KK, NAFF, NAFFT 
         
        SPTPOS(:,I) = POS;   SPTNRM(:,I) = VECTOR
        NATMAFF(I)  = NAFF;  NSFTYPE(I)  = NAFFT
    END DO
    CLOSE(102)

! read the surface triangulization

    OPEN(103,FILE=pathname(1:lenpath)//FNAME(1:lenfname)//".face")

    READ(103,*) FHEAD
    READ(103,*) FHEAD
    READ(103,*) NFACE, PPP, QQQ, RRR

    ALLOCATE(NVERT(3,NFACE), MFACE(NFACE), STAT=ierr)
    IF (ierr .NE. 0) THEN
        WRITE(6,*) 'Error allocating NVERT, MFACE'
    STOP
    END IF

    NVERT=0; MFACE=0
      
    DO I=1,NFACE 
        READ(103,*) NIND(1:5) 
        NVERT(1:3,I) = NIND(1:3);  MFACE(I) = NIND(4)
    END DO
    CLOSE(103)
    call surface_area(s_area) ! the post-MSMS code
    print *,'surface area=', real(s_area)

    rslt=system('rm '//pathname(1:lenpath)//fname(1:lenfname)//".xyzr")
    rslt=system('rm '//pathname(1:lenpath)//fname(1:lenfname)//".vert")
    rslt=system('rm '//pathname(1:lenpath)//fname(1:lenfname)//".face")


elseif (imd==1) then
    pi=acos(-1.d0)
    zero_vec=(/0.d0, 0.d0, 0.d0/)

    tvec=(/0.d0, 0.d0, 30.d0/)              !translation of the original domain
    angle=(/pi/2.d0,pi/4.d0,pi/3.d0/)       !Euler angles for the rotation
    write(*,*) 'traslation vector: ',real(tvec)
    write(*,*) 'rotation angles x,y,z: ',real(angle)

    nremark=0
    open(102,file=pathname(1:lenpath)//fname(1:lenfname)//".pqr")
    do
        READ(102,*) fhead
        if (fhead(1:6)=='REMARK') then
            nremark=nremark+1
        else
            exit
        endif
    enddo
    !print *,'lines of remarks = ', nremark
    close(102)

    open(102,file=pathname(1:lenpath)//fname(1:lenfname)//".pqr")
    open(103,file=pathname(1:lenpath)//fname(1:lenfname)//".xyzr")

    do i=1,nremark
        read(102,*) fhead
    enddo

    natm=0
    do
        read(102,*,IOSTAT = MEOF) c1,i2,c3,c4,i5,xyzqr 

        if ((c1(1:3) .ne. 'END') .and. (MEOF .eq. 0)) then
            write(103,*) xyzqr(1:3),xyzqr(5)
            natm=natm+1
        endif
        IF(MEOF .LT. 0) EXIT
    enddo

    close(102)
    close(103)


    !Read atom coordinates and partial charges
    nchr=natm
    allocate(atmpos(3,natm),atmrad(natm),atmchr(nchr),chrpos(3,nchr),STAT=ierr)
    IF (ierr .NE. 0) THEN
        WRITE(6,*) 'Error allocating atmpos, atmrad, atmchr, chrpos!'
        STOP
    END IF

    open(102,file=pathname(1:lenpath)//fname(1:lenfname)//".pqr")

    do i=1,nremark
        read(102,*) fhead
    enddo

    do i=1,natm
        read(102,*,IOSTAT = MEOF) c1,i2,c3,c4,i5,xyzqr 
        atmpos(1:3,i)=xyzqr(1:3)
        atmrad(i)=xyzqr(5)
        chrpos(1:3,i)=xyzqr(1:3)
        atmchr(i)=xyzqr(4)
    enddo

    close(102)

    allocate(atmpos_copy(3,natm),STAT=ierr)
    IF (ierr .NE. 0) THEN
        WRITE(6,*) 'Error allocating atmpos!'
        STOP
    END IF
        
    ! find coodinates of the 2nd monomer by given tranlation and rotation
    do i=1,natm
        call TransRot(atmpos_copy(:,i),atmpos(:,i),tvec,angle)
    enddo

    ! check if steric clashes exist
    iflag_steric_clashes=0
    do i=1,natm
        do j=1,natm
            center1=atmpos(:,i)
            center2=atmpos_copy(:,j)
            diff=center1-center2
            dist=sqrt(dot_product(diff,diff))
            if (dist<(atmrad(i)+atmrad(j)+0.0)) then
                write(*,*) 'WARNING: steric clashes, regenerate molecular surface using dimer'
                write(*,*) 'dist(i,j)<rad(i)+rad(j)+2*R_wt, if sig. drop of elements #, reduce R_wt'
                write(*,*) i,j,real(dist),real(atmrad(i)),real(atmrad(j)),real(atmrad(i)+atmrad(j)+2.2)
                iflag_steric_clashes=1
                goto 101
            endif 
        enddo
    enddo
    101 continue
        
    if (iflag_steric_clashes==1) then
        ! generate .xyzr file using dimer otherwise using monomer and its copy
        open(103,file=pathname(1:lenpath)//fname(1:lenfname)//".xyzr")
        do i=1,natm
            write(103,*) atmpos(1:3,i),atmrad(i)
        enddo
        do i=1,natm
            write(103,*) atmpos_copy(1:3,i),atmrad(i)
        enddo
        close(103)
            
        ! adjust geometry information by doubling the molecule size
        ALLOCATE(atmchr_copy(natm),atmrad_copy(natm),chrpos_copy(3,nchr),STAT=ierr)
        IF (ierr .NE. 0) THEN
            WRITE(6,*) 'Error allocating copies of atmchr, atmpos, atmrad, chrpos, for multiple domain'
            STOP
        END IF
        atmchr_copy=atmchr; atmpos_copy=atmpos; atmrad_copy=atmrad; chrpos_copy=chrpos

        deallocate(atmchr,atmpos,atmrad,chrpos,STAT=ierr)
        IF (ierr .NE. 0) THEN
            WRITE(6,*) 'Error deallocating atmchr,atmpos,atmrad,chrpos for multipole domain'
            STOP
        EndIF

        ALLOCATE(atmchr(2*natm),atmpos(3,2*natm),atmrad(2*natm),chrpos(3,2*nchr),STAT=ierr)
        IF (ierr .NE. 0) THEN
            WRITE(6,*) 'Error allocating copies of atmchr, atmpos, atmrad, chrpos, for multiple domain'
            STOP
        END IF

        atmchr(1:natm)=atmchr_copy;
        atmrad(1:natm)=atmrad_copy;
        atmpos(:,1:natm)=atmpos_copy;
        chrpos(:,1:nchr)=chrpos_copy;

        atmchr(natm+1:2*natm)=atmchr_copy;
        atmrad(natm+1:2*natm)=atmrad_copy;
        
        do i=1,natm
            call TransRot(atmpos(:,natm+i),atmpos_copy(:,i),tvec,angle)
        enddo
        do i=1,nchr
            call TransRot(chrpos(:,nchr+i),chrpos_copy(:,i),tvec,angle)
        enddo
    endif

    rslt=system('msms -if '//pathname(1:lenpath)//fname(1:lenfname)//".xyzr"//' -prob 1.4 -de ' &
    //den(1:5)//' -of '//pathname(1:lenpath)//fname(1:lenfname)//' > msms.out')    
      ! read the surface points
    OPEN(102,FILE=pathname(1:lenpath)//FNAME(1:lenfname)//".vert")

    READ(102,*) FHEAD
    READ(102,*) FHEAD
    READ(102,*) NSPT, ppp, qqq, rrr

    ALLOCATE(SPTPOS(3,NSPT), SPTNRM(3,NSPT), NATMAFF(NSPT), NSFTYPE(NSPT), STAT= ierr)
    IF (ierr .NE. 0) THEN
        WRITE(6,*) 'Error allocating SPTPOS, SPTNRM, NATMAFF, NSFTYPE!'
    STOP
    END IF

    SPTPOS=0.D0; SPTNRM=0.D0; NATMAFF=0; NSFTYPE=0;	
      
    DO I=1,NSPT
        READ(102,*) POS(1:3), VECTOR(1:3), KK, NAFF, NAFFT 
         
        SPTPOS(:,I) = POS;   SPTNRM(:,I) = VECTOR
        NATMAFF(I)  = NAFF;  NSFTYPE(I)  = NAFFT
    END DO
    CLOSE(102)

! read the surface triangulization

    OPEN(103,FILE=pathname(1:lenpath)//FNAME(1:lenfname)//".face")

    READ(103,*) FHEAD
    READ(103,*) FHEAD
    READ(103,*) NFACE, PPP, QQQ, RRR

    ALLOCATE(NVERT(3,NFACE), MFACE(NFACE), STAT=ierr)
    IF (ierr .NE. 0) THEN
        WRITE(6,*) 'Error allocating NVERT, MFACE'
    STOP
    END IF

    NVERT=0; MFACE=0
      
    DO I=1,NFACE 
        READ(103,*) NIND(1:5) 
        NVERT(1:3,I) = NIND(1:3);  MFACE(I) = NIND(4)
    END DO
    CLOSE(103)
    call surface_area(s_area) ! the post-MSMS code
    print *,'surface area=', real(s_area)

    rslt=system('rm '//pathname(1:lenpath)//fname(1:lenfname)//".xyzr")
    rslt=system('rm '//pathname(1:lenpath)//fname(1:lenfname)//".vert")
    rslt=system('rm '//pathname(1:lenpath)//fname(1:lenfname)//".face")

    write(*,*) 'cube dimensions for monomer (w/o steric clashes) or dimmer (w/ steric clahes): '
    write(*,*) 'x: [ ',real(minval(sptpos(1,:))),real(maxval(sptpos(1,:))),']' 
    write(*,*) 'y: [ ',real(minval(sptpos(2,:))),real(maxval(sptpos(2,:))),']' 
    write(*,*) 'z: [ ',real(minval(sptpos(3,:))),real(maxval(sptpos(3,:))),']' 

    if (iflag_steric_clashes==0) then 
        ALLOCATE(SPTPOS_copy(3,NSPT), SPTNRM_copy(3,NSPT), nvert_copy(3,nface), STAT= ierr)
        IF (ierr .NE. 0) THEN
            WRITE(6,*) 'Error allocating SPTPOS_copy, SPTNRM_copy, nvert_copy for multiple domain'
            STOP
        END IF
        SPTPOS_copy=sptpos; SPTNRM_copy=sptnrm; nvert_copy=nvert;
    
        ALLOCATE(atmchr_copy(natm),atmrad_copy(natm),chrpos_copy(3,nchr),STAT=ierr)
        IF (ierr .NE. 0) THEN
            WRITE(6,*) 'Error allocating copies of atmchr, atmpos, atmrad, chrpos, for multiple domain'
            STOP
        END IF
        atmchr_copy=atmchr; atmpos_copy=atmpos; atmrad_copy=atmrad; chrpos_copy=chrpos
        
        deallocate(nvert,sptpos,sptnrm,atmchr,atmpos,atmrad,chrpos,STAT=ierr)
        IF (ierr .NE. 0) THEN
            WRITE(6,*) 'Error deallocating nvert,sptpos,sptnrm,atmchr,atmpos,atmrad,chrpos for multipole domain'
            STOP
        EndIF

        ALLOCATE(SPTPOS(3,2*NSPT), SPTNRM(3,2*NSPT), nvert(3,2*nface), STAT= ierr)
        IF (ierr .NE. 0) THEN
            WRITE(6,*) 'Error allocating SPTPOS, SPTNRM, nvert for multiple domain'
            STOP
        END IF

        ALLOCATE(atmchr(2*natm),atmpos(3,2*natm),atmrad(2*natm),chrpos(3,2*nchr),STAT=ierr)
        IF (ierr .NE. 0) THEN
            WRITE(6,*) 'Error allocating copies of atmchr, atmpos, atmrad, chrpos, for multiple domain'
            STOP
        END IF

        
        sptpos(:,1:nspt)=sptpos_copy; 
        sptnrm(:,1:nspt)=sptnrm_copy;
        nvert(:,1:nface)=nvert_copy;
        atmchr(1:natm)=atmchr_copy;
        atmrad(1:natm)=atmrad_copy;
        atmpos(:,1:natm)=atmpos_copy;
        chrpos(:,1:nchr)=chrpos_copy;

        nvert(:,nface+1:2*nface)=nvert_copy+nspt
        atmchr(natm+1:2*natm)=atmchr_copy;
        atmrad(natm+1:2*natm)=atmrad_copy;
        
        ! translation and rotation
        do i=1,nspt
            call TransRot(sptpos(:,nspt+i),sptpos_copy(:,i),tvec,angle)
            call TransRot(sptnrm(:,nspt+i),sptnrm_copy(:,i),zero_vec,angle)
        enddo
        do i=1,natm
            call TransRot(atmpos(:,natm+i),atmpos_copy(:,i),tvec,angle)
        enddo
        do i=1,nchr
            call TransRot(chrpos(:,nchr+i),chrpos_copy(:,i),tvec,angle)
        enddo

       
        ! translation only
        !do i=1,3 
        !    sptpos(i,nspt+1:2*nspt)=sptpos_copy(i,:)+tvec(i); 
        !    sptnrm(i,nspt+1:2*nspt)=sptnrm_copy(i,:); 
        !    atmpos(i,natm+1:2*natm)=atmpos_copy(i,:)+tvec(i);
        !    chrpos(i,nchr+1:2*nchr)=chrpos_copy(i,:)+tvec(i);
        !enddo
        
        nface=2*nface
        nspt=2*nspt
    endif
        
    open(102,file=pathname(1:lenpath)//fname(1:lenfname)//".pqr")
    open(103,file=fname(1:lenfname)//"_monomer1.pqr")
    open(104,file=fname(1:lenfname)//"_monomer2.pqr")

    do i=1,nremark
        read(102,*) fhead
    enddo

    do i=1,natm
        read(102,*,IOSTAT = MEOF) c1,i2,c3,c4,i5

        if ((c1(1:3) .ne. 'END') .and. (MEOF .eq. 0)) then
            !read(102,'(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2)',IOSTAT = MEOF) & 
            write(103,'(A6,I5,1X,A5,A4,1X,I6,3X,5F10.5)') &
            c1,i2,c3,c4,i5,atmpos(:,i),atmchr(i),atmrad(i)
            write(104,'(A6,I5,1X,A5,A4,1X,I6,3X,5F10.5)') &
            c1,i2+max_iatm,c3,c4,i5+max_ires,atmpos(:,i+natm),atmchr(i+natm),atmrad(i+natm)
        endif
        IF(MEOF .LT. 0) EXIT
    enddo

    close(102)
    close(103)
    close(104)
    rslt=system('cat '//fname(1:lenfname)//"_monomer1.pqr " & 
    //fname(1:lenfname)//"_monomer2.pqr > "//fname(1:lenfname)//"_dimer.pqr")

    deALLOCATE(atmchr_copy,atmpos_copy,atmrad_copy,chrpos_copy,STAT=ierr)
    IF (ierr .NE. 0) THEN
    WRITE(6,*) 'Error deallocating copies of atmchr, atmpos, atmrad, chrpos, for multiple domain'
        STOP
    END IF

    natm=2*natm
    nchr=2*nchr
    write(*,*) 'cube dimensions for the updated dimer (w/o steric slashes): '
    write(*,*) 'x: [ ',real(minval(sptpos(1,:))),real(maxval(sptpos(1,:))),']' 
    write(*,*) 'y: [ ',real(minval(sptpos(2,:))),real(maxval(sptpos(2,:))),']' 
    write(*,*) 'z: [ ',real(minval(sptpos(3,:))),real(maxval(sptpos(3,:))),']' 


! imd=2: dimmer formed by reading two monomer information  !gan: changed into using icosahedral grid for spheres.
else            
    nremark1=0
    nremark2=0
    open(102,file=pathname(1:lenpath)//fname(1:lenfname)//"_sph1.pqr")
    do
        READ(102,*) fhead
        if (fhead(1:6)=='REMARK') then
            nremark1=nremark1+1
        else
            exit
        endif
    enddo
    close(102)
    open(103,file=pathname(1:lenpath)//fname(1:lenfname)//"_sph2.pqr")
    do
        READ(103,*) fhead
        if (fhead(1:6)=='REMARK') then
            nremark2=nremark2+1
        else
            exit
        endif
    enddo
    close(103)

    open(102,file=pathname(1:lenpath)//fname(1:lenfname)//"_sph1.pqr")
    do i=1,nremark1
        read(102,*) fhead
    enddo
    natm1=0
    do
        read(102,*,IOSTAT = MEOF) xyzqr

        if ((c1(1:3) .ne. 'END') .and. (MEOF .eq. 0)) then
            natm1=natm1+1
            !write(*,*) natm1,real(xyzqr)
        endif
        !print *,natm1,MEOF
        IF(MEOF .LT. 0) EXIT
    enddo
    close(102)

    open(103,file=pathname(1:lenpath)//fname(1:lenfname)//"_sph2.pqr")
    do i=1,nremark2
        read(103,*) fhead
    enddo
    natm2=0
    do
        read(103,*,IOSTAT = MEOF) xyzqr

        if ((c1(1:3) .ne. 'END') .and. (MEOF .eq. 0)) then
            natm2=natm2+1
        endif
        IF(MEOF .LT. 0) EXIT
    enddo
    close(103)

    !Read atom coordinates and partial charges
    nchr1=natm1
    nchr2=natm2
    natm=natm1+natm2
    nchr=natm
    allocate(atmpos(3,natm),atmrad(natm),atmchr(nchr),chrpos(3,nchr),STAT=ierr)
    IF (ierr .NE. 0) THEN
        WRITE(6,*) 'Error allocating atmpos, atmrad, atmchr, chrpos!'
        STOP
    END IF

    open(101,file=pathname(1:lenpath)//fname(1:lenfname)//"_sph1.pqr")
    do i=1,nremark1
        read(101,*) fhead
    enddo
    do i=1,natm1
        read(101,*,IOSTAT = MEOF) xyzqr
        atmpos(1:3,i)=xyzqr(1:3)
        atmrad(i)=xyzqr(5)
        chrpos(1:3,i)=xyzqr(1:3)
        atmchr(i)=xyzqr(4)
    enddo
    close(101)

    open(102,file=pathname(1:lenpath)//fname(1:lenfname)//"_sph2.pqr")
    do i=1,nremark2
        read(102,*) fhead
    enddo
    do i=1,natm2
        read(102,*,IOSTAT = MEOF) xyzqr
        atmpos(1:3,i+natm1)=xyzqr(1:3)
        atmrad(i+natm1)=xyzqr(5)
        chrpos(1:3,i+natm1)=xyzqr(1:3)
        atmchr(i+natm1)=xyzqr(4)
    enddo
    close(102)

    ! check if steric clashes exist
    iflag_steric_clashes=0
    do i=1,1 !gan
        do j=1,1 !gan
            center1=atmpos(:,i)
            center2=atmpos(:,j+natm1)
            diff=center1-center2
            dist=sqrt(dot_product(diff,diff))
            if (dist<(atmrad(i)+atmrad(j+natm1)+0.0)) then !changed into 0.0 to avoid merge --Gan
                write(*,*) 'WARNING: steric clashes, regenerate molecular surface using dimer'
                write(*,*) 'dist(i,j)<rad(i)+rad(j)+2*R_wt, if sig. drop of elements #, reduce R_wt'
                write(*,*) i,j,real(dist),real(atmrad(i)),real(atmrad(j+natm1)),real(atmrad(i)+atmrad(j+natm1)+2.2)
                iflag_steric_clashes=1
                goto 102
            endif 
        enddo
    enddo
    write(*,*) 'No steric clashes!'
    102 continue

    if (iflag_steric_clashes==1) then
        open(103,file=pathname(1:lenpath)//fname(1:lenfname)//".xyzr")
        do i=1,natm1
             write(103,*) atmpos(1:3,i),atmrad(i)
        enddo
        do j=1,natm2
             write(103,*) atmpos(1:3,j+natm1),atmrad(j+natm1)
        enddo
        close(103)

        rslt=system('msms -if '//pathname(1:lenpath)//fname(1:lenfname)//".xyzr"//' -prob 1.4 -de ' &
        //den(1:5)//' -of '//pathname(1:lenpath)//fname(1:lenfname)//' > msms.out')    
        ! read the surface points
        OPEN(102,FILE=pathname(1:lenpath)//FNAME(1:lenfname)//".vert")

        READ(102,*) FHEAD
        READ(102,*) FHEAD
        READ(102,*) NSPT, ppp, qqq, rrr

        ALLOCATE(SPTPOS(3,NSPT), SPTNRM(3,NSPT), NATMAFF(NSPT), NSFTYPE(NSPT), STAT= ierr)
        IF (ierr .NE. 0) THEN
            WRITE(6,*) 'Error allocating SPTPOS, SPTNRM, NATMAFF, NSFTYPE!'
        STOP
        END IF

        SPTPOS=0.D0; SPTNRM=0.D0; NATMAFF=0; NSFTYPE=0;	
      
        DO I=1,NSPT
            READ(102,*) POS(1:3), VECTOR(1:3), KK, NAFF, NAFFT 
         
            SPTPOS(:,I) = POS;   SPTNRM(:,I) = VECTOR
            NATMAFF(I)  = NAFF;  NSFTYPE(I)  = NAFFT
        END DO
        CLOSE(102)

        ! read the surface triangulization

        OPEN(103,FILE=pathname(1:lenpath)//FNAME(1:lenfname)//".face")

        READ(103,*) FHEAD
        READ(103,*) FHEAD
        READ(103,*) NFACE, PPP, QQQ, RRR

        ALLOCATE(NVERT(3,NFACE), MFACE(NFACE), STAT=ierr)
        IF (ierr .NE. 0) THEN
            WRITE(6,*) 'Error allocating NVERT, MFACE'
        STOP
        END IF

        NVERT=0; MFACE=0
      
        DO I=1,NFACE 
            READ(103,*) NIND(1:5) 
            NVERT(1:3,I) = NIND(1:3);  MFACE(I) = NIND(4)
        END DO
        CLOSE(103)
        call surface_area(s_area) ! the post-MSMS code
        print *,'surface area=', real(s_area)

        rslt=system('rm '//pathname(1:lenpath)//fname(1:lenfname)//".xyzr")
        rslt=system('rm '//pathname(1:lenpath)//fname(1:lenfname)//".vert")
        rslt=system('rm '//pathname(1:lenpath)//fname(1:lenfname)//".face")


    else !no steric clashes thus two monomers are generated by calling MSMS
        open(103,file=pathname(1:lenpath)//fname(1:lenfname)//".xyzr")
        do i=1,natm1
             write(103,*) atmpos(1:3,i),atmrad(i)
        enddo
close(103)

!Read atom coordinates and partial charges
OPEN(1,FILE=pathname(1:lenpath)//FNAME(1:lenfname)//".xyzr")
DO
READ(1,*,IOSTAT = MEOF) xxx, yyy, zzz, rrr
IF(MEOF .LT. 0) EXIT
END DO
CLOSE(1)
! These two are for future usage
init_type=1;	iwind=1

! For the icosahedron, initially 20 faces, 12 vertex and 32 edges
! With the relation V+F=E+2 and 12+20=30+2
nface=20;	nspt=12;	nedge=nface+nspt-2

! nefv contains the number of edges,faces and vertexes on each level
allocate(nefv(3,0:level),STAT=ierr)
IF (ierr .NE. 0) THEN
WRITE(6,*) 'Error allocating number of edges,vertex,and face'
STOP
END IF
nefv=0;	nefv(:,0)=(/nedge,nface,nspt/)

do i=1,level
nedge=nedge*2+3*nface
nface=nface*4
nspt=nedge+2-nface
nefv(:,i)=(/nedge,nface,nspt/)
enddo

! allocate varibles for triangulation
allocate(sptpos(3,nspt),sptnrm(3,nspt),nvert(3,nface), STAT=ierr)
IF (ierr .NE. 0) THEN
WRITE(6,*) 'Error allocating vertex and faces'
STOP
END IF
sptpos=0.d0;	nvert=0

ALLOCATE(NATMAFF(NSPT), NSFTYPE(NSPT), MFACE(NFACE), STAT= ierr)
IF (ierr .NE. 0) THEN
WRITE(6,*) 'Error allocating NATMAFF, NSFTYPE, MFACE!'
STOP
END IF

!rad=atmrad(1)
!center=center1
!rds=rad
center=(/0,0,0/)
rad=1.0
call sphere_tri(init_type,level,rad,center,iwind,nspt,nface,sptpos,nvert,nefv)
allocate(tr_area(2*nface),STAT=ierr)
tr_area=0.d0
do i=1,nface
idx=nvert(1:3,i) ! vertices index of the specific triangle
aa=sqrt(dot_product(sptpos(:,idx(1))-sptpos(:,idx(2)),sptpos(:,idx(1))-sptpos(:,idx(2))))
bb=sqrt(dot_product(sptpos(:,idx(1))-sptpos(:,idx(3)),sptpos(:,idx(1))-sptpos(:,idx(3))))
cc=sqrt(dot_product(sptpos(:,idx(2))-sptpos(:,idx(3)),sptpos(:,idx(2))-sptpos(:,idx(3))))
 tr_area(i)=triangle_area(aa,bb,cc) !gan: flat area
!tr_area(i)=SphereTriArea(sptpos(:,idx(1)),sptpos(:,idx(2)),sptpos(:,idx(3))) !gan: curved area
tr_area(nface+i)=tr_area(i)
end do

sptnrm=sptpos
rad=atmrad(1)
center1=atmpos(:,1)

sptpos(1,1:nspt)=sptpos(1,1:nspt)*rad+center1(1)
sptpos(2,1:nspt)=sptpos(2,1:nspt)*rad+center1(2)
sptpos(3,1:nspt)=sptpos(3,1:nspt)*rad+center1(3)

do i=1,nface
tr_area(i)=tr_area(i)*rad*rad
tr_area(nface+i)=tr_area(nface+i)*atmrad(natm1+1)*atmrad(natm1+1)
!gan: curved
end do






        ! readin data for monomer2
        ! first copy monomer1 information 
        NSPT1=NSPT; NFACE1=NFACE
        ALLOCATE(SPTPOS_copy(3,NSPT1), SPTNRM_copy(3,NSPT1), NVERT_copy(3,NFACE1), STAT= ierr)
        IF (ierr .NE. 0) THEN
            WRITE(6,*) 'Error allocating SPTPOS_copy, SPTNRM_copy, NVERT_copy!'
        STOP
        END IF
        SPTPOS_copy=SPTPOS
        SPTNRM_copy=SPTNRM
        NVERT_copy=NVERT

        deALLOCATE(SPTPOS, SPTNRM, NVERT, NATMAFF, NSFTYPE, MFACE, STAT= ierr)
        IF (ierr .NE. 0) THEN
            WRITE(6,*) 'Error deallocating SPTPOS, SPTNRM, NVERT, NATMAFF, NSFTYPE, MFACE!'
        STOP
        END IF

        open(103,file=pathname(1:lenpath)//fname(1:lenfname)//".xyzr")
        do j=1,natm2
             write(103,*) atmpos(1:3,j+natm1),atmrad(j+natm1)
        enddo
        close(103)

        !again, do that with icosahedral grid (gan)
!Read atom coordinates and partial charges
OPEN(1,FILE=pathname(1:lenpath)//FNAME(1:lenfname)//".xyzr")
DO
READ(1,*,IOSTAT = MEOF) xxx, yyy, zzz, rrr
IF(MEOF .LT. 0) EXIT
END DO
CLOSE(1)


! allocate varibles for triangulation
allocate(sptpos(3,nspt),sptnrm(3,nspt),nvert(3,nface), STAT=ierr)
IF (ierr .NE. 0) THEN
WRITE(6,*) 'Error allocating vertex and faces'
STOP
END IF
sptpos=0.d0;	nvert=0

ALLOCATE(NATMAFF(NSPT), NSFTYPE(NSPT), MFACE(NFACE), STAT= ierr)
IF (ierr .NE. 0) THEN
WRITE(6,*) 'Error allocating NATMAFF, NSFTYPE, MFACE!'
STOP
END IF

rad=atmrad(natm1+1)
center2=atmpos(:,natm1+1)

!call sphere_tri(init_type,level,rad,center,iwind,nspt,nface,sptpos,nvert,nefv)
sptpos(1,1:nspt)=(SPTPOS_copy(1,1:nspt)-center1(1))*rad/atmrad(1)+center2(1)
sptpos(2,1:nspt)=(SPTPOS_copy(2,1:nspt)-center1(2))*rad/atmrad(1)+center2(2)
sptpos(3,1:nspt)=(SPTPOS_copy(3,1:nspt)-center1(3))*rad/atmrad(1)+center2(3)
sptnrm=SPTNRM_copy
nvert=NVERT_copy


!finish twice for monomer two





        
        NSPT2=NSPT; NFACE2=NFACE
        ALLOCATE(SPTPOS_copy2(3,NSPT2), SPTNRM_copy2(3,NSPT2), NVERT_copy2(3,NFACE2), STAT= ierr)
        IF (ierr .NE. 0) THEN
            WRITE(6,*) 'Error allocating SPTPOS_copy, SPTNRM_copy, NVERT_copy!'
        STOP
        END IF
        SPTPOS_copy2=SPTPOS
        SPTNRM_copy2=SPTNRM
        NVERT_copy2=NVERT
        
        deALLOCATE(SPTPOS, SPTNRM, NVERT, STAT= ierr)
        IF (ierr .NE. 0) THEN
            WRITE(6,*) 'Error deallocating SPTPOS, SPTNRM, NVERT!'
        STOP
        END IF

        NSPT=NSPT1+NSPT2
        NFACE=NFACE1+NFACE2
        ALLOCATE(SPTPOS(3,NSPT), SPTNRM(3,NSPT), NVERT(3,NFACE), STAT= ierr)
        IF (ierr .NE. 0) THEN
            WRITE(6,*) 'Error allocating SPTPOS, SPTNRM, NVERT!'
        STOP
        END IF

        sptpos(:,1:nspt1)=sptpos_copy; 
        sptnrm(:,1:nspt1)=sptnrm_copy;
        nvert(:,1:nface1)=nvert_copy;


        sptpos(:,(nspt1+1):nspt)=sptpos_copy2; 
        sptnrm(:,(nspt1+1):nspt)=sptnrm_copy2;
        nvert(:,(nface1+1):nface)=nvert_copy2+nspt1;

        deALLOCATE(SPTPOS_copy, SPTNRM_copy, NVERT_copy, SPTPOS_copy2, SPTNRM_copy2, NVERT_copy2, STAT= ierr)
        IF (ierr .NE. 0) THEN
            WRITE(6,*) 'Error deallocating SPTPOS, SPTNRM, NVERT!'
        STOP
        END IF
    endif
endif 

!do i=1,80
!print *,'check location', nvert(:,i), nvert(:,i+80)
!end do

   End

!--------------------------------------------------------------------------------------
Subroutine sphere_tri(init_type,level,rad,center,iwind,nspt,nface,sptpos,nvert,nefv)
! init_type: 1: tetrahedron; 2:
implicit none
integer init_type,level,iwind,nvert(3,nface),nspt,nface,ico_f0(3,20),i,nefv(3,0:level)
real*8 rad,center(3),sptpos(3,nspt),t,tau,one,ico_v0(3,12)

!Twelve vertices of icosahedron on unit sphere

!tau = 0.8506508084; one = 0.5257311121
t=(1.d0+sqrt(5.d0))/2.d0
tau=t/sqrt(1.d0+t**2)
one=1/sqrt(1.d0+t**2)

! 12 vertex
ico_v0(:,1) = (/  tau,  one,	0.d0 /); ! ZA
ico_v0(:,2) = (/ -tau,  one,	0.d0 /); ! ZB
ico_v0(:,3) = (/ -tau, -one,	0.d0 /); ! ZC
ico_v0(:,4) = (/  tau, -one,	0.d0 /); ! ZD
ico_v0(:,5) = (/  one, 0.d0,	tau  /); ! YA
ico_v0(:,6) = (/  one, 0.d0,	-tau /); ! YB
ico_v0(:,7) = (/ -one, 0.d0,	-tau /); ! YC
ico_v0(:,8) = (/ -one, 0.d0,	tau  /); ! YD
ico_v0(:,9) = (/ 0.d0,  tau,	one  /); ! XA
ico_v0(:,10)= (/ 0.d0, -tau,	one  /); ! XB
ico_v0(:,11)= (/ 0.d0, -tau,	-one /); ! XC
ico_v0(:,12)= (/ 0.d0,  tau,	-one /); ! XD


! Structure for unit icosahedron
! 20 faces
ico_f0(:,1) =	(/5,  8,  9/)
ico_f0(:,2) =	(/5, 10,  8/)
ico_f0(:,3) =	(/6, 12,  7/)
ico_f0(:,4) =	(/6,  7, 11/)
ico_f0(:,5) =	(/1,  4,  5/)
ico_f0(:,6) =	(/1,  6,  4/)
ico_f0(:,7) =	(/3,  2,  8/)
ico_f0(:,8) =	(/3,  7,  2/)
ico_f0(:,9) =	(/9, 12,  1/)
ico_f0(:,10)=	(/9,  2, 12/)
ico_f0(:,11)=	(/10, 4, 11/)
ico_f0(:,12)=	(/10, 11, 3/)
ico_f0(:,13)=	(/9,  1,  5/)
ico_f0(:,14)=	(/12, 6,  1/)
ico_f0(:,15)=	(/5,  4, 10/)
ico_f0(:,16)=	(/6, 11,  4/)
ico_f0(:,17)=	(/8,  2,  9/)
ico_f0(:,18)=	(/7, 12,  2/)
ico_f0(:,19)=	(/8, 10,  3/)
ico_f0(:,20)=	(/7,  3, 11/)

sptpos(:,1:12)=ico_v0
nvert(:,1:20)=ico_f0


do i=1,level
call mesh_refine_tri4(i,nspt,nface,sptpos,nvert)
call sphere_project(nspt,sptpos,rad,nefv(3,i),center);
enddo

End

!---------------------------------------------------------------------------------------
Subroutine sphere_project(nspt,sptpos,rad,nspt_level,center)
implicit none
integer nspt,nspt_level,isptpos
real*8 sptpos(3,nspt),rad,center(3),x(3),xnorm

do isptpos=1,nspt_level
x=sptpos(:,isptpos)-center
xnorm=sqrt(dot_product(x,x))
x=x/xnorm
x=x*rad
sptpos(:,isptpos)=x+center
enddo


End

!---------------------------------------------------------------------------------------
Subroutine mesh_refine_tri4(ilevel,nspt,nface,sptpos,nvert)
implicit none
real*8, allocatable :: sptpos0(:,:)
integer, allocatable:: nvert0(:,:)

integer ilevel,i,nspt,nface,nface0,nspt0,nedge0,nvert(3,nface),err,NABC(3)
integer indx_sptpos,iface,NA,NB,NC,ifind
real*8 sptpos(3,nspt),ABC(3,3),A(3),B(3),C(3)

nface0=20
nspt0=12
nedge0=nface0+nspt0-2
do i=1,ilevel-1
nedge0=nedge0*2+3*nface0
nface0=nface0*4
nspt0=nedge0+2-nface0
enddo

allocate(sptpos0(3,nspt0),nvert0(3,nface0), STAT=err)
IF (err .NE. 0) THEN
WRITE(6,*) 'Error allocating local vertex and faces'
STOP
END IF
sptpos0=0.d0
nvert0=0

sptpos0(:,1:nspt0)=sptpos(:,1:nspt0)
nvert0(:,1:nface0)=nvert(:,1:nface0)

indx_sptpos=nspt0

do iface = 1,nface0
!Get the triangle vertex indices
NA = nvert0(1,iface)
NB = nvert0(2,iface)
NC = nvert0(3,iface)

!Get the triangle vertex coordinates
A = sptpos0(:,NA)
B = sptpos0(:,NB)
C = sptpos0(:,NC)

!Now find the midpoints between vertices
ABC(:,1) = (A + B) / 2;
ABC(:,2) = (B + C) / 2;
ABC(:,3) = (C + A) / 2;

!Store the midpoint vertices, while
!checking if midpoint vertex already exists

do i=1,3
call mesh_find_vertex(indx_sptpos,ifind,sptpos,nspt,ABC(:,i));
if (ifind==0) then
indx_sptpos=indx_sptpos+1
sptpos(:,indx_sptpos)=ABC(:,i)
NABC(i)=indx_sptpos
else
NABC(i)=ifind
endif
enddo

!Create new faces with orig vertices plus midpoints

nvert(:,iface*4-3) = (/ NA, NABC(1), NABC(3) /)
nvert(:,iface*4-2) = (/ NABC(1), NB, NABC(2) /)
nvert(:,iface*4-1) = (/ NABC(3), NABC(2), NC /)
nvert(:,iface*4-0) = (/ NABC(1), NABC(2), NABC(3)/)

enddo

deallocate(sptpos0,nvert0, STAT=err)
IF (err .NE. 0) THEN
WRITE(6,*) 'Error deallocating local vertex and faces'
STOP
END IF

End
!------------------------------------------------------------------------
subroutine surface_area(s_area)
use molecule
implicit double precision(a-h,o-z)
integer iface(3),jface(3),nfacenew,ialert
real*8 face(3,3),s_area,face_old(3,3),xx(3),yy(3),cpu1,cpu2
real*8,allocatable:: nvert_copy(:,:)

    print *,'# of surfaces=',nface,' # of surface points=',nspt
    call cpu_time(cpu1)
    s_area=0.d0
    nfacenew=nface
    allocate(nvert_copy(3,nface))
    nvert_copy=nvert
    do i=1,nface
        iface=nvert(:,i)
        xx=0.d0
        ialert=0;
        do j=1,3
            face(:,j)=sptpos(:,iface(j))
            xx=xx+1/3.d0*(face(:,j))
        enddo
        aa=sqrt(dot_product(face(:,1)-face(:,2),face(:,1)-face(:,2)))
        bb=sqrt(dot_product(face(:,1)-face(:,3),face(:,1)-face(:,3)))
        cc=sqrt(dot_product(face(:,2)-face(:,3),face(:,2)-face(:,3)))
        area_local=triangle_area(aa,bb,cc)
        !area_local=SphereTriArea(face(:,1),face(:,2),face(:,3)) !the function only works for 
!a unit sphere centered at origin

        do ii=max(1,i-10),i-1
            jface=nvert(:,ii)
            yy=0.d0
            do j=1,3
                face_old(:,j)=sptpos(:,jface(j))
                yy=yy+1/3.d0*(face_old(:,j))
            enddo
            dist_local=dot_product(xx-yy,xx-yy)
            if (dist_local<1.d-5) then
                ialert=1
                !print *,i,ii,'particles are too close',dist_local
            endif
        enddo
            if (area_local < 1.d-5 .or. ialert==1) then
                !print *,i,j,'small area=', area_local
                ichanged=nface-nfacenew
                nvert_copy(:,(i-ichanged):(nface-1))=nvert_copy(:,(i-ichanged+1):nface)
                nfacenew=nfacenew-1
            endif
            s_area=s_area+area_local
    enddo
    print *,nface-nfacenew,' ugly faces are deleted'
        
    nface=nfacenew
    deallocate(nvert,STAT=ierr)
    IF (ierr .NE. 0) THEN
        WRITE(6,*) 'Error deallocating nvert'
        STOP
    EndIF
       
    allocate(nvert(3,nface),STAT=ierr)
    IF (ierr .NE. 0) THEN
        WRITE(6,*) 'Error allocating nvert'
        STOP
    EndIF
    nvert=nvert_copy(:,1:nface)
    deallocate(nvert_copy, STAT=ierr)
    IF (ierr .NE. 0) THEN
        WRITE(6,*) 'Error deallocating nvert_copy'
        STOP
    EndIF
    call cpu_time(cpu2)
    print *, 'total MSMS post-processing time =',cpu2-cpu1
end

!------------------------------------------------------------------------
function triangle_area(aa,bb,cc)
implicit double precision(a-h,o-z)
s=0.5d0*(aa+bb+cc)
triangle_area=sqrt(s*(s-aa)*(s-bb)*(s-cc))
end

!------------------------------------------------------------------------
function SphereTriArea(xyzA, xyzB, xyzC) result(SphereTriAreaVector)
implicit double precision(a-h,o-z)
! Calculates the area of a spherical triangle on the unit sphere

!   NOTE : This function requires function SphereArcLength.

! Calling parameters

real*8 xyzA(3), xyzB(3), xyzC(3) !vertex coordinates

real*8 :: SphereTriAreaVector

! Local variables

real*8 :: side1, side2, side3, halfPerimeter, zz,rdssq



side1 = SphereArcLength(xyzA,xyzB)

side2 = SphereArcLength(xyzB,xyzC)

side3 = SphereArcLength(xyzC,xyzA)



halfPerimeter =(side1 + side2 + side3)/2.d0

rdssq=xyzA(1)*xyzA(1)+xyzA(2)*xyzA(2)+xyzA(3)*xyzA(3)

zz = tan(halfPerimeter/2.d0)*tan( (halfPerimeter-side1)/2.d0)*&

tan( (halfPerimeter - side2)/2.d0)*tan( (halfPerimeter - side3)/2.d0)


!print *,'check', side1, side2, side3, zz, 4.d0*atan2(sqrt(zz),1.d0)
SphereTriAreaVector = 4.d0*rdssq*atan2(sqrt(zz),1.d0)

end



function SphereArcLength (xyzA, xyzB) result(SphereArcLengthVector)
implicit double precision(a-h,o-z)
! returns the arc length between two points on the surface of a sphere.

! Calling parameters

real*8 xyzA(3), xyzB(3)

real*8 :: SphereArcLengthVector

! Local variables

real*8 crossProd(3)

real*8 :: dotProd , crossNorm



crossProd = [xyzA(2)*xyzB(3)-xyzB(2)*xyzA(3),xyzB(1)*xyzA(3)-xyzA(1)*xyzB(3),&

xyzA(1)*xyzB(2)-xyzB(1)*xyzA(2) ]



crossNorm = sqrt(sum(crossProd*crossProd))



dotProd = xyzA(1)*xyzB(1)+xyzA(2)*xyzB(2)+xyzA(3)*xyzB(3)



SphereArcLengthVector = atan2(crossNorm,dotProd)

end

!---------------------------------------------------------------------------------
! For a new vertex, compared with all the stored vertex
! ifind .ne. 0:	if the vertex is already stored
! ifind=0:		if the vertex is completely new
Subroutine mesh_find_vertex(indx_sptpos,ifind,sptpos,nspt,sptpos_new)
implicit none
real*8 sptpos(3,nspt),sptpos_new(3),diff(3)
integer indx_sptpos,ifind,nspt,isptpos

ifind=0

do isptpos=1,indx_sptpos
    diff=sptpos_new-sptpos(:,isptpos)
    if (sqrt(dot_product(diff,diff))<1.d-10) then
        ifind = isptpos
        return
    endif
enddo

End

!-----------------------------------------------------------------------
subroutine TransRot(tpos,pos,tvec,angle)
implicit none
real*8 pos(3),tvec(3),angle(3),tpos(3)

!local variables
real*8 COSA,SINA,COSB,SINB,COSG,SING,C(3,3) 
integer i 
!-----------------------------------------------------------------------
! purpose: compute euler rotation matrix 
!-----------------------------------------------------------------------
!
! C_ALPHA   : rotation about z-axis  
! C_BETA    : rotation about new y-axis  
! C_GAMMA   : rotation about new z-axis      
! C: rotation to model coordinates 
! C: rotation to geographic coordinates
!
!-----------------------------------------------------------------------

      !compute sin and cos of all angles in degrees
      COSA  = dcos(angle(1))
      SINA  = dsin(angle(1))
      COSB  = dcos(angle(2))
      SINB  = dsin(angle(2))
      COSG  = dcos(angle(3))
      SING  = dsin(angle(3))
 
      !compute rotation matrix into model coordinates
      !C(1,1) =   COSA * COSB * COSG  -  SINA * SING
      !C(1,2) =   SINA * COSB * COSG  +  COSA * SING
      !C(1,3) =        - SINB * COSG
      !C(2,1) = - COSA * COSB * SING  -  SINA * COSG
      !C(2,2) = - SINA * COSB * SING  +  COSA * COSG
      !C(2,3) =          SINB * SING
      !C(3,1) =   COSA * SINB
      !C(3,2) =   SINA * SINB
      !C(3,3) =          COSB
    
      !compute rotation matrix into geographical coordinates
      C(1,1) =   COSG * COSB * COSA  -  SING * SINA
      C(1,2) = - SING * COSB * COSA  -  COSG * SINA
      C(1,3) =          SINB * COSA
      C(2,1) =   COSG * COSB * SINA  +  SING * COSA
      C(2,2) = - SING * COSB * SINA  +  COSG * COSA
      C(2,3) =          SINB * SINA
      C(3,1) = - COSG * SINB
      C(3,2) =   SING * SINB
      C(3,3) =          COSB

      Do i=1,3
         tpos(i)=dot_product(C(i,:),pos)+tvec(i)
      EndDo
   
end subroutine TransRot


!#############################################################################################
!The face file contains three header lines followed by one triangle per line. 
!The first header line provides a comment and the file name of the sphere set. 
!The second header line holds comments about the content of the third line. 
!The third header line provides the number of triangles, the number of spheres in the set, 
!the triangulation density and the probe sphere radius. 

!The first three numbers are (1 based) vertex indices. 

!The next field can be: 
!1 for a triangle in a toric reen trant face, 
!2 for a triangle in a spheric reentrant face and 
!3 for a triangle in a contact face. 

!The last # on the line is the (1 based) face number in the analytical description of the solvent excluded surface. 
!These values are written in the following format ``%6d %6d %6d %2d %6d''.

!The vertex file contains three header lines (similar to the header in the .face file) 
!followed by one vertex per line and provides the coordinates (x,y,z) and the normals (nx,ny,nz) 
!followed by the number of the face (in the analytical description of the solvent excluded surface) 
!to which the vertex belongs. The vertices of the analytical surface have a value 0 in that field 
!and the vertices lying on edges of this surface have negative values. 
!The next field holds the (1 based) index of the closest sphere. 
!The next field is 
!1 for vertices which belong to toric reentrant faces (including ver tices of the analytical surface), 
!2 for vertices inside reentrant faces and 
!3 for vertices inside contact faces. 
!These values are written in the following format ``%9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %7d %7d %2d''.
