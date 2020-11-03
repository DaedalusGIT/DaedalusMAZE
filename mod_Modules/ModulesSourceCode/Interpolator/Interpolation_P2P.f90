!=====================DAEDALUS_INTERPOLATION==================================!
! THESE SUBROUTINES ARE USED BY DAEDALUS INTERPOLATOR TO INTERPOLATE TIEGCM DATA
! TO DAEDALUS' ORBIT. THE PARALLEL VERSION INCORPORATES OpenMP FOR MULTITHREADING
! THE INTERPOLATING METHOD IS A 3D TRILINEAR SPATIAL INTERPOLATION. TO COMPILE
! WITH F2PY3 WRAPPER FOR PYTHON:
! Serial_Version:
! No Optimization: "f2py3 -c -lgomp --opt='-O2' Interpolation.f90 -m interpolation"
! Parallel_Version:
! No Optimization: "f2py3 -c --f90flags="-fopenmp " -lgomp  Interpolation_Parallel_Version.f95 -m interpolation"
! -O1 Level Optimization: "f2py3 -c --f90flags="-fopenmp " -lgomp --opt='-O1' Interpolation_Parallel_Version.f95 -m interpolation"
! -O2 Level Optimization: "f2py3 -c --f90flags="-fopenmp " -lgomp --opt='-O2' Interpolation_Parallel_Version.f95 -m interpolation"
! -O3 Level Optimization: "f2py3 -c --f90flags="-fopenmp " -lgomp --opt='-O3' Interpolation_Parallel_Version.f95 -m interpolation"
! -ffast Level Optimization: "f2py3 -c --f90flags="-fopenmp " -lgomp --opt='-ffast' Interpolation_Parallel_Version.f95 -m interpolation"
!KEP,2019
!============================================================================!
                  !  ____                 _       _
                 ! |  _ \  __ _  ___  __| | __ _| |_   _ ___
                 ! | | | |/ _` |/ _ \/ _` |/ _` | | | | / __|
                 ! | |_| | (_| |  __/ (_| | (_| | | |_| \__ \
                 ! |____/ \__,_|\___|\__,_|\__,_|_|\__,_|___/

subroutine interp_p2p(glat,glon,glev,dlat,dlon,dalt,ne,zg,&
                              dphi,dtheta,counter,res)
implicit none
! ======================Data Declarations=======================================
real*8,intent(in)::glat(:),glon(:),glev(:)
real*8,intent(in)::ne(:,:,:,:),zg(:,:,:,:)
real*8,intent(in)::dtheta,dphi,dlat,dlon,dalt
real*8:: dx,dy,dz,drho,w1,w2,w3,w4,w5,w6,w7,w8,m
integer,intent(in):: counter
integer::i,llat,llon,lalt
real*8,dimension(:),allocatable::alts
real*8,intent(inout)::res(:)

! ==============================================================================


! daedalus positions in LLA -->glat,glon,galt
! Find local grid positions and closeest neighbours in each dimesnion
! lat-->theta
! lon-->phi
! alt..>rho

! llat=(size(glat)/2) +1+ floor(dlat/dtheta)
! llon=(size(glon)/2) +1+ floor(dlon/dphi)

call local_positions(size(glat),dlat,glat,llat)
call local_positions(size(glon),dlon,glon,llon)





if (dlon .ge.177.5) then
  llon=size(glon)-1
endif
if (dlon .le.-177.5) then
  llon=1
endif
if (dlat .ge.87.5) then
  llat=size(glat)-1
endif
if (dlat .le.-87.5) then
  llat=1
endif

! Translate pressure levels to altitude from geographic height
allocate(alts(size(glev)))
alts(:)=zg(counter,:,llat,llon)/100000
call local_positions(size(alts),dalt,alts,lalt)
drho=alts(lalt+1)-alts(lalt)


! Find relative distances
dx= abs(dalt-alts(lalt))/drho
dy=abs(dlat-glat(llat))/dtheta
dz=abs(dlon-glon(llon))/dphi

!

!Find weights for interpolation
w1=(1-dx)*(1-dy)*(1-dz)
w2=(dx)*(1-dy)*(1-dz)
w3=(1-dx)*(dy)*(1-dz)
w4=(dx)*(dy)*(1-dz)
w5=(1-dx)*(1-dy)*(dz)
w6=(dx)*(1-dy)*(dz)
w7=(1-dx)*(dy)*(dz)
w8=(dx)*(dy)*(dz)

m=0.0
! Perform the Interpolation
m=       ne(counter,lalt,llat,llon)*w1
m=m+  ne(counter,lalt+1,llat,llon)*w2
m=m+  ne(counter,lalt,llat+1,llon)*w3
m=m+  ne(counter,lalt+1,llat+1,llon)*w4
m=m+  ne(counter,lalt,llat,llon+1)*w5
m=m+  ne(counter,lalt+1,llat,llon+1)*w6
m=m+  ne(counter,lalt,llat+1,llon+1)*w7
m=m+  ne(counter,lalt+1,llat+1,llon+1)*w8

!Store Result for Output
res(:)=m
return
end subroutine






subroutine local_positions(dim1,y,x,local_pos)
! Input:   dim1--> integer dimension of matrix x for allocation purposes
!          y--> real array of kind dp --> main grid pencil to find local positions in
!          x--> real kind dp -->global position
! Output:  localpos--> integer kind 4 local position in grid pencil
!KEP,2019

implicit none
integer,intent(in)::dim1
real*8,dimension(dim1),intent(in)::x
real*8,intent(in)::y
integer::ii,jj
integer,intent(out)::local_pos


do ii=1,size(x)-1
  if (y .ge. x(ii) .and. y.lt.x(ii+1)) then
    local_pos=ii
    exit
  endif

enddo

end subroutine
