!=====================DAEDALUS_INTERPOLATION====================================!
! THESE SUBROUTINES ARE USED BY DAEDALUS INTERPOLATOR TO INTERPOLATE TIEGCM DATA
! TO DAEDALUS' ORBIT. THE PARALLEL VERSION INCORPORATES OpenMP FOR MULTITHREADING
! THE INTERPOLATING METHOD IS A 3D TRILINEAR SPATIAL INTERPOLATION. TO COMPILE
! WITH F2PY3 WRAPPER FOR PYTHON:
! Serial_Version
! No Optimization: "f2py3 -c Interpolation_Serial.f90 -m interpolation"
!
!  Parallel_Version
!  No Optimization: "f2py3 -c --f90flags="-fopenmp " -lgomp  Interpolation_Parallel.f90 -m interpolation"
!  -O1 Level Optimization: "f2py3 -c --f90flags="-fopenmp " -lgomp --opt='-O1' Interpolation_Parallel.f90 -m interpolation_OpenMP"
!  -O2 Level Optimization: "f2py3 -c --f90flags="-fopenmp " -lgomp --opt='-O2' Interpolation_Parallel.f90 -m interpolation_OpenMP"
!  -O3 Level Optimization: "f2py3 -c --f90flags="-fopenmp " -lgomp --opt='-O3' Interpolation_Parallel.f90 -m interpolation_OpenMP"
!  -ffast Level Optimization: "f2py3 -c --f90flags="-fopenmp " -lgomp --opt='-ffast' Interpolation_Parallel.f90 -m interpolation_OpenMP"
! KEP 5/6/19
!===============================================================================!
                  !  ____                 _       _
                 ! |  _ \  __ _  ___  __| | __ _| |_   _ ___
                 ! | | | |/ _` |/ _ \/ _` |/ _` | | | | / __|
                 ! | |_| | (_| |  __/ (_| | (_| | | |_| \__ \
                 ! |____/ \__,_|\___|\__,_|\__,_|_|\__,_|___/


subroutine daed_interp(grid_lat,grid_lon,grid_lev,daed_lat,&
                         daed_lon,daed_alt,zg,model_dt,orbit_dt,ne,m)

implicit none
!=========================DATA DECLARATIONS=====================================
real*8,intent(in)::grid_lat(:),grid_lon(:),grid_lev(:),daed_lat(:),daed_lon(:),&
daed_alt(:)
real*8,intent(in)::ne(:,:,:,:),zg(:,:,:,:)
real*8,intent(in)::model_dt,orbit_dt
real*8::deltaphi,deltatheta,deltarho,dx,dy,dz,w1,w2,w3,w4,w5,w6,w7,w8,lat_temp
integer::i,j,k,jj,grid_update,dim1,dim2,counter,r_local,theta_local,phi_local,kk
real*8,dimension(:),allocatable::alts
real*8,intent(inout)::m(:)
real::t1,t2
! ==============================================================================


!Call Clock
call cpu_time(t1)


! **********PARAMETERS**********************
deltaphi= abs(grid_lon(3)-grid_lon(2))
deltatheta=abs(grid_lat(3)-grid_lat(2))
grid_update=floor(model_dt/orbit_dt)

j=1
counter=2
! *****************************************


! Some allocations
allocate(alts(size(grid_lev)))

! ~~~~~~~~~~~~~~~~~~~~~~~~~MAIN LOOP~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



! Go through daedalus' orbit
do i=1,size(daed_alt)


  ! Find Daedalus local positions r,theta
  call local_positions(size(grid_lon),daed_lon(i),grid_lon,phi_local)
  call local_positions(size(grid_lat),daed_lat(i),grid_lat,theta_local)



  ! If outside of grid modify neighbors
  if (daed_lon(i) .ge.177.5) then
    phi_local=size(grid_lon)-1
  endif
  if (daed_lon(i) .le.-177.5) then
    phi_local=1
  endif
  if (daed_lat(i) .ge.87.5) then
    theta_local=size(grid_lat)-1
  endif
  if (daed_lat(i) .le.-87.5) then
    theta_local=1
  endif



  ! Translate pressure levels to altitude from geographic height
  alts(:)=zg(counter,:,theta_local,phi_local)/100000
  call local_positions(size(alts),daed_alt(i),alts,r_local)

  ! Now find Delta-Rho
  deltarho=alts(r_local+1)-alts(r_local)


  dx=abs(((daed_alt(i)-alts(r_local))/deltarho))
  dy=abs(((daed_lat(i)-grid_lat(theta_local))/deltatheta))
  dz=abs(((daed_lon(i)-grid_lon(phi_local))/deltaphi))


  !Find weights for interpolation
  w1=(1-dx)*(1-dy)*(1-dz)
  w2=(dx)*(1-dy)*(1-dz)
  w3=(1-dx)*(dy)*(1-dz)
  w4=(dx)*(dy)*(1-dz)
  w5=(1-dx)*(1-dy)*(dz)
  w6=(dx)*(1-dy)*(dz)
  w7=(1-dx)*(dy)*(dz)
  w8=(dx)*(dy)*(dz)

  ! Perform the interpolation
  ! m(i)=0.0
  m(i)=       ne(counter,r_local,theta_local,phi_local)*w1
  m(i)=m(i)+  ne(counter,r_local+1,theta_local,phi_local)*w2
  m(i)=m(i)+  ne(counter,r_local,theta_local+1,phi_local)*w3
  m(i)=m(i)+  ne(counter,r_local+1,theta_local+1,phi_local)*w4
  m(i)=m(i)+  ne(counter,r_local,theta_local,phi_local+1)*w5
  m(i)=m(i)+  ne(counter,r_local+1,theta_local,phi_local+1)*w6
  m(i)=m(i)+  ne(counter,r_local,theta_local+1,phi_local+1)*w7
  m(i)=m(i)+  ne(counter,r_local+1,theta_local+1,phi_local+1)*w8


enddo

!Call Clock
call cpu_time(t2)
! print*, "Interpolation finished in:",t2-t1,"seconds"

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
