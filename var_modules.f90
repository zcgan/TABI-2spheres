MODULE MOLECULE

    integer :: natm, nspt, nface, nchr, nt
    integer :: imd, iflag_steric_clashes,nface1, nchr1 !/* for multiple domain*/
    integer :: natm1, natm2, nchr2, nface2, nspt1, nspt2

    real*8  :: eps, kappa, rds, eps0, eps1, para

    real*8,  dimension(:),     allocatable :: atmrad, atmchr
    real*8,  dimension(:,:),   allocatable :: atmpos , sptpos ,sptnrm, chrpos, chrpos_sph 
    real*8,  dimension(:,:,:), allocatable :: chgmnx
    integer, dimension(:),     allocatable :: natmaff, nsftype, mface 
    integer, dimension(:,:),   allocatable :: nvert

END MODULE MOLECULE

MODULE COMDATA   

   character(100) :: fname,pathname, den
   
   integer :: lenpath,lenfname

   real*8, dimension(:), allocatable :: bvct,xvct,F_exa, xtemp
   real*8, dimension(:,:), allocatable :: amtrx
   integer, dimension(:), allocatable :: indx
!######################################################
integer,allocatable,dimension(:):: iwork,IGWK
real*8,allocatable,dimension(:):: rwork,sb,sx,RGWK
!######################################################
           
END MODULE COMDATA 


module bicg
integer nmax
real*8, dimension(:), allocatable:: sa
integer, dimension(:), allocatable:: ijka
end module bicg


!------------------------------------------------------------------
module treecode

! r8 is 8-byte (double precision) real

      !INTEGER,PARAMETER :: r8=SELECTED_REAL_KIND(12)

! runtime parameters

      INTEGER :: numpars,order,maxparnode,iflag,forcedim
      
      !REAL(KIND=r8) :: theta 
      real*8 :: theta
! arrays for coordinates, charge, potential & force (tree and direct) 

      REAL*8,ALLOCATABLE,DIMENSION(:) :: x,y,z,q
      !REAL(KIND=r8),ALLOCATABLE,DIMENSION(:) :: tpoten,dpoten
      !REAL*8,ALLOCATABLE,DIMENSION(:,:) :: tforce,dforce
      INTEGER,ALLOCATABLE,DIMENSION(:) :: orderind

! timing variables

      REAL*8 :: timebeg,timeend

! local variables

      !INTEGER :: i,j,err
      REAL*8 :: xyzminmax(6)
      REAL*8 :: t1,abserr,relerr,absinf_err,relinf_err
      REAL*8,DIMENSION(3) :: f_inferr,f_relinferr,t


      real*8, dimension(:), allocatable:: tr_area
      real*8, dimension(:,:), allocatable:: tr_xyz, tr_q
      real*8, dimension(:,:,:), allocatable:: tchg,schg
      real*8, dimension(:,:,:,:), allocatable:: der_cof
      integer, dimension(:,:), allocatable:: kk
end module treecode
