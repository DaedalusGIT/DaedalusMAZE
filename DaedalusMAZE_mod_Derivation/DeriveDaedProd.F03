/*
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

program DeriveDaedalusProducts
use apex
implicit none
integer, parameter :: ntimes =59, nlat = 72,nlon=144,nlev=57, nKP=24
! ntimes=1 ! for 074 file
! ntimes=11 ! for 148 file
! ntimes=59 ! for 221 file
! integer, parameter :: ntimes = 59, nlat = 72,nlon=144,nlev=57, nKP=24 !221 file
! integer, parameter :: ntimes = 59, nlat = 72,nlon=144,nlev=57, nKP=24 !221 file

integer:: required=MPI_THREAD_FUNNELED,provided,mpierror,numprocs,myid,root=0
integer, parameter :: dp = selected_real_kind(8)
integer, parameter :: dp1 = selected_real_kind(4)
integer, parameter :: dp2 = selected_real_kind(15)
character(len=54) :: LOG='logfile.txt',arg
character(len=6),parameter::outdir="Output"
character(len=256):: mkdirCmd,mkCpCmd,readdir
character(len=512)::rnCmd
character*200,dimension(:),allocatable::ncfile,localfiles
integer::timescnt,files,latcnt,levcnt,loncnt,filecnt,sndcnt,recvcnt,threads=1
integer::ncid,ncerror,varid,jj
real*8:: t0,t1,t2, t_init,t_finish

real(dp),parameter::deltalmd = 2.5  ! TIEGCM grid resolution in deg
real(dp),parameter::deltaphi = 2.5  ! TIEGCM grid resolution in deg
real(dp),parameter::Rearth = 6378.1370  ! Earths Radius (km)
real(dp),parameter::electron = 1.602176565 * 1e-19 !electron charge in C
real(dp),parameter::boltzmann=1.380645852 * 1e-16 !Boltzmann constant in cm^2*g*s^(-2)*K^(-1)
real(dp),parameter::me=9.10938356 * 1e-31 !electron mass in kg
integer, parameter ::mO=16 ! Oxygen atomic mass in g/mol
integer, parameter ::mN2=28  !N2 molecular mass in g/mol
integer, parameter ::mO2=32 !O2 molecular mass in g/mol
integer, parameter ::mNO=30 !NO molecular mass in g/mol
real(dp), parameter ::NA=6.02214076 * 1e23 !Avogadro's constant in mol^-1
real(dp), parameter ::fb=1.5 !Burnside facor (default=1.5)
real(dp), parameter ::q_e=1.60217657e-19  ! Coulomb
real(dp), parameter ::ArO = 15.9994  ! atomic mass O=16
real(dp), parameter ::ArN = 14.0067  ! atomic mass N=14
real(dp), parameter ::NAvog = 6.02214076e23
real(dp), parameter ::m_NOplus = (ArO + ArN) / (NAvog * 1000)!kg
real(dp), parameter ::m_O2plus = 2 * ArO / (NAvog * 1000)!kg
real(dp), parameter ::m_Oplus = 1 * ArO / (NAvog * 1000)!kg
real(dp), parameter ::m_e = 9.11e-31  ! kg
real(dp), parameter ::MOmass = mO / (NAvog * 1000)! Atomic oxygen mass in kg
real(dp), parameter ::M_ion = m_Oplus * 10 ** 3  ! in grams
integer(dp), parameter ::CubicCm2CubicM = 1e6
real(dp):: time_igrf = 2015.2! in python will be used this datetime.datetime(2015, 3, 15, 0, 0, 0), first time step of the TIEGCM run
integer:: initial_yr=2017
integer:: initial_mon=1
integer:: initial_day=1
integer:: initial_hour=0
integer:: initial_mins=0
real(dp):: initial_secs=0.0
integer,parameter:: secsOFyr=365*24*60*60
integer:: daysFRinit,restMINS,hrFRinit,minFRinit
integer::current_year,current_month,current_day,current_hour,current_mins
real(dp):: current_secs=0.0
logical :: dirExists

!apex variables
real(dp1):: bmag,si,qdlon,malat,vmp,W,D,Be3,sim,qdlat,F
real(dp1),dimension(3):: B,bhat,d1,d2,d3,e1,e2,e3,f1,f2,f3,g1,g2,g3
real(dp1):: sec,sbsllat,sbsllon,csza,sza,mlt
integer:: ier,i,iyr,iday,ihr,imn
real(dp1)::date
!netcdf variables
real(dp1),dimension(:),allocatable:: glat,glon,glev,gtime,KP
real(dp),dimension(:,:,:,:),allocatable::zg,zg2 !Geometric height in cm
real(dp),dimension(:,:,:,:),allocatable::Ne_all,Ne_all2 !electron density in cm^-3
real(dp),dimension(:,:,:,:),allocatable::PED_all,PED_all2 !Pedersen conductivity in Si/m
real(dp),dimension(:,:,:,:),allocatable::HALL_all,HALL_all2 !Hall conductivity in Si/m
real(dp),dimension(:,:,:,:),allocatable::Un_east,Un_east2 !neutral zonal wind (+East) in cm/s
real(dp),dimension(:,:,:,:),allocatable::Un_north,Un_north2 !neutral meridional wind (+North) in cm/s
real(dp),dimension(:,:,:,:),allocatable::Un_up,Un_up2 !neutral vertical wind (+Up) in cm/s
real(dp),dimension(:,:,:,:),allocatable::Ui_east,Ui_east2 !zonal ExB velocity in cm/s
real(dp),dimension(:,:,:,:),allocatable::Ui_north,Ui_north2 !meridional ExB velocity in cm/s
real(dp),dimension(:,:,:,:),allocatable::Ui_up,Ui_up2 !vertical ExB velocity in cm/s
real(dp),dimension(:,:,:,:),allocatable::Tn_all,Tn_all2 !Neutral Temperature in K
real(dp),dimension(:,:,:,:),allocatable::Te_all,Te_all2 !Electron Temperature in K
real(dp),dimension(:,:,:,:),allocatable::Ti_all,Ti_all2 !Ion Temperature in K
real(dp),dimension(:,:,:,:),allocatable::nO_all,nO_all2 !neutral atomic oxygen density in cm^-3
real(dp),dimension(:,:,:,:),allocatable::nO2_all,nO2_all2 !neutral molecular oxygen density in cm^-3
real(dp),dimension(:,:,:,:),allocatable::nN2_all,nN2_all2 !neutral molecular nitrogen density in cm^-3
real(dp),dimension(:,:,:,:),allocatable::NOp_all,NOp_all2 !O+ density in cm^-3
real(dp),dimension(:,:,:,:),allocatable::NO2p_all,NO2p_all2 !O2+ density in cm^-3
real(dp),dimension(:,:,:,:),allocatable::NNOp_all,NNOp_all2 !NO+ density in cm^-3
real(dp),dimension(:,:,:,:),allocatable::E_X,E_X2 !E field in V/cm
real(dp),dimension(:,:,:,:),allocatable::E_Y,E_Y2 !E field in V/cm
real(dp),dimension(:,:,:,:),allocatable::E_Z,E_Z2 !E field in V/cm

real(dp),dimension(:,:,:),allocatable::BX_BMAG,BX_BMAG2 !BX/BMAG ratio
real(dp),dimension(:,:,:),allocatable::BY_BMAG,BY_BMAG2 !BY/BMAG ratio
real(dp),dimension(:,:,:),allocatable::BZ_BMAG,BZ_BMAG2 !BZ/BMAG ratio
real(dp),dimension(:,:,:),allocatable::BMAG_G,BMAG2_G !magnetic field in Gauss

!other Variables
real(dp),dimension(:),allocatable::time_pf
real(dp),dimension(:,:,:,:),allocatable::zgf,zgout
real(dp),dimension(:,:,:,:),allocatable::mlat_qdf,mlat_qdfout
real(dp),dimension(:,:,:,:),allocatable::mlt_qdf,mlt_qdfout
real(dp),dimension(:,:,:,:),allocatable::BX_T ,BX_Tout!Tesla
real(dp),dimension(:,:,:,:),allocatable::BY_T,BY_Tout !Tesla
real(dp),dimension(:,:,:,:),allocatable::BZ_T,BZ_Tout !Tesla
real(dp),dimension(:,:,:,:),allocatable::BMAG_T,BMAG_Tout !Tesla
real(dp),dimension(:,:,:,:),allocatable::Un_east_si,Un_east_si_out !neutral zonal wind (+East) in m/s
real(dp),dimension(:,:,:,:),allocatable::Un_north_si,Un_north_si_out !neutral meridional wind (+North) in m/s
real(dp),dimension(:,:,:,:),allocatable::Un_up_si,Un_up_si_out !neutral vertical wind (+Up) in m/s
real(dp),dimension(:,:,:,:),allocatable::Ui_east_si,Ui_east_si_out !zonal ExB velocity in m/s
real(dp),dimension(:,:,:,:),allocatable::Ui_north_si,Ui_north_si_out !meridional ExB velocity in m/s
real(dp),dimension(:,:,:,:),allocatable::Ui_up_si,Ui_up_si_out !vertical ExB velocity in m/s
real(dp),dimension(:,:,:,:),allocatable::E_X_si,E_X_si_out !E field in V/m
real(dp),dimension(:,:,:,:),allocatable::E_Y_si,E_Y_si_out !E field in V/m
real(dp),dimension(:,:,:,:),allocatable::E_Z_si,E_Z_si_out !E field in V/m
!collision frequencies arrays
real(dp),dimension(:,:,:,:),allocatable::nu_O2p_O2,nu_O2p_O2_out
real(dp),dimension(:,:,:,:),allocatable::nu_O2p_O,nu_O2p_O_out
real(dp),dimension(:,:,:,:),allocatable::nu_O2p_N2,nu_O2p_N2_out
real(dp),dimension(:,:,:,:),allocatable::nu_O2p,nu_O2p_out
real(dp),dimension(:,:,:,:),allocatable::nu_Op_O,nu_Op_O_out
real(dp),dimension(:,:,:,:),allocatable::nu_Op_O2,nu_Op_O2_out
real(dp),dimension(:,:,:,:),allocatable::nu_Op_N2,nu_Op_N2_out
real(dp),dimension(:,:,:,:),allocatable::nu_Op,nu_Op_out
real(dp),dimension(:,:,:,:),allocatable::nu_NOp_O2,nu_NOp_O2_out
real(dp),dimension(:,:,:,:),allocatable::nu_NOp_O,nu_NOp_O_out
real(dp),dimension(:,:,:,:),allocatable::nu_NOp_N2,nu_NOp_N2_out
real(dp),dimension(:,:,:,:),allocatable::nu_NOp,nu_NOp_out
real(dp),dimension(:,:,:,:),allocatable::nu_e,nu_e_out
real(dp),dimension(:,:,:,:),allocatable::nu_in,nu_in_out
!conductivities calculated by Richmond
real(dp),dimension(:,:,:,:),allocatable::PED_calc,PED_calcout !Pedersen conductivity in Si/m
real(dp),dimension(:,:,:,:),allocatable::HALL_calc,HALL_calcout !Hall conductivity in Si/m
real(dp),dimension(:,:,:,:),allocatable::PAR_calc,PAR_calcout !Parallel conductivity in Si/m

!Heating Rates with TIEGCM conductivities
real(dp),dimension(:,:,:,:),allocatable::Ohmic,Ohmicout
real(dp),dimension(:,:,:,:),allocatable::Joule,Jouleout
real(dp),dimension(:,:,:,:),allocatable::Frictional,Frictionalout
real(dp),dimension(:,:,:,:),allocatable::Convenction_heating,Convenction_heating_out
real(dp2),dimension(:,:,:,:),allocatable::Wind_heating,Wind_heating_out
!Heating Rates with calculated conductivities
real(dp),dimension(:,:,:,:),allocatable::Ohmic_calc,Ohmic_calcout
real(dp),dimension(:,:,:,:),allocatable::Joule_calc,Joule_calcout
real(dp),dimension(:,:,:,:),allocatable::Frictional_calc,Frictional_calcout

!Cross Sections calculated
real(dp),dimension(:,:,:,:),allocatable::CS_Op,CS_Opout
real(dp),dimension(:,:,:,:),allocatable::CS,CSout
real(dp),dimension(:,:,:,:),allocatable::CS_O2p,CS_O2pout
real(dp),dimension(:,:,:,:),allocatable::CS_NOp,CS_NOpout
!Currents Calculations
real(dp),dimension(:,:,:,:),allocatable::Pedersen_current,Pedersen_current_out
real(dp),dimension(:,:,:,:),allocatable::Hall_current,Hall_current_out
real(dp),dimension(:,:,:,:),allocatable::Total_current,Total_current_out

!help variables
real(dp)::timeg,Tn,Te,Ti,k,mlat_qd,mlon_qd,mlt_qd,NO,NO2,NN2,Ne,NOp,NO2p,NNOp,temp,&
        &temp1,temp2,temp3,&
        & Une,Unn,Unu,Un_e,Un_n,Un_u,Uie,Uin,Uiu,Ui_e,Ui_n,Ui_u,BX_i,BY_i,BZ_i,bnorm,Be,Bn,Bu,Unvertx,Unverty,Unvertz,Unvertmag,&
        & Ri,Omega_e,Omega_Op,Omega_O2p,Omega_NOp,ratio_e,ratio_Op,ratio_O2p,ratio_NOP,&
        & alpha,beta,gamma,delta,alpha1,beta1,gamma1,delta1,N_neutral,M_ionSUM,MnN2,&
        & Evertx,Everty,Evertz,Estarnorm,Vi_starnorm,sigmaPED_file,sigmaHALL_file,UNxBnorm
real(dp1)::alt_p,lat_p,lon_p,geod_lat
real(dp),dimension(3)::UnENU,Un,Ui,UiENU,Benu,Bxyz,b_unit, Unvert,Unvertecef, Evert,Estar,Vi_star,CrossResult,&
        & Jpedersen,Jhall,UNxB,ExB,Efieldenu,Efieldecef
! --------------------------------------------------------------------------
!magnetic parameters variables
real(dp)::dipmag,sndip,decmag,sndec,sn2dec,csdec,tm1,tm2,tm3
!netcdf export help variables
integer, parameter :: NDIMS=4,NDIMS_B=3
integer::ncOUTid,ncOUTerror,varOUTid1,varOUTid2,varOUTid3,varOUTid4,varOUTid5,varOUTid6,varOUTid7,varOUTid8,varOUTid9,&
        &varOUTid10,varOUTid11,varOUTid12,varOUTid13,varOUTid14,varOUTid15,&
        &varOUTid17,varOUTid18,varOUTid19,varOUTid20,varOUTid21,varOUTid22,varOUTid23,varOUTid24,&
        &varOUTid25,varOUTid26,varOUTid27,varOUTid28,varOUTid29,varOUTid30,varOUTid31,varOUTid32,varOUTid33,&
        &varOUTid34,varOUTid35,varOUTid36,varOUTid37,varOUTid38,varOUTid39,&
        &varOUTid16,NCntimesID,NCnlevID,NCnlatID,NCnlonID
integer OUTDimids1(NDIMS),OUTDimids2(NDIMS_B)
integer::  times_dmind, lat_dmind,lon_dmind,lev_dmind
character*200::ncINfileNAMEsplited
character*250::ncOUTfileNAME
integer::name_index
! Init MPI
t0=MPI_Wtime()
call MPI_INIT_THREAD(required,provided,mpierror)
call check_MPI_status(mpierror)
call MPI_COMM_SIZE ( mpi_comm_world, numprocs, mpierror )
call check_MPI_status(mpierror)
call MPI_COMM_RANK ( mpi_comm_world, myid, mpierror )
call check_MPI_status(mpierror)
call omp_set_num_threads(threads)
if (provided .ne. required )then
   stop 'Abort:: MPI not initialized correclty'
endif
t1=MPI_Wtime()


! Root creates directory for Output Files
! if (myid == root  )then
!   inquire(file=trim(outdir)//'/.', exist=dirExists )
!   if (.not.dirExists) then
!     mkdirCmd = 'mkdir -p '//trim(outdir)
!     call system( mkdirCmd )
!   endif
! endif



! Root reads filenames
if (myid.eq.root) then
  ! to 61 to dialegoume tixaia?
  open(61,file=LOG,action='write')
  write(61,*) "========================= JouleHeatingDerivation==========================="
  write(61,*) 'MPI:: MPI Initialised with', numprocs ,'MPI Tasks and ', threads , 'number of threads in ', t1-t0, "seconds"
  write(61,*) "==========================================================================="
  close(61)

  i = 1
  do
    call get_command_argument(i, arg)
    if (len_trim(arg) == 0) exit
    i = i+1
  end do

  allocate (ncfile(i-1))
  filecnt=i-1
  i = 1

  do
    call get_command_argument(i, arg)
     if (len_trim(arg) == 0) exit
    ncfile(i)=trim(arg)
    i = i+1
  end do

  if (MOD(filecnt,numprocs) .eq. 0) then
    sndcnt=filecnt/numprocs
  else
    stop 'Use different number of MPI Tasks'
  endif
  recvcnt=sndcnt

  ! LOG
  open(61,file=LOG,action='write',position='append')
  write(61,*) "LOG::List Of Files to Run:"
    do jj=1,filecnt
    write(61,*) jj,"-->",ncfile(jj)
    enddo
  write(61,*) "LOG::Number of files per MPI Task :" ,recvcnt
  write(61,*) "==========================================================================="
  close(61)

endif

! LoadBalance
t0=MPI_Wtime()
call MPI_Bcast(sndcnt,  1, MPI_INTEGER, root, mpi_comm_world,mpierror)
call check_MPI_status(mpierror)
call MPI_Bcast(recvcnt, 1, MPI_INTEGER, root, mpi_comm_world,mpierror)
call check_MPI_status(mpierror)
call MPI_Bcast(filecnt, 1, MPI_INTEGER, root, mpi_comm_world,mpierror)
call check_MPI_status(mpierror)
if (myid .ne.root) then
   allocate(ncfile(filecnt))
endif
call MPI_Bcast(ncfile,filecnt*200, MPI_CHAR, root, mpi_comm_world,mpierror)
call check_MPI_status(mpierror)
allocate(localfiles(recvcnt))
call MPI_Scatter(ncfile,sndcnt*200,MPI_CHAR,localfiles,recvcnt*200,MPI_CHAR ,root,mpi_comm_world,mpierror)

!Memory Allocation for TIEGCM variables
!every processor job
allocate(glat(nlat))
allocate(glon(nlon))
allocate(glev(nlev))
allocate(gtime(ntimes))
! allocate(KP(nKP))
allocate(KP(ntimes))
allocate(zg(ntimes,nlev,nlat,nlon))
allocate(zg2(nlon,nlat,nlev,ntimes))
allocate(Ne_all(ntimes,nlev,nlat,nlon))
allocate(Ne_all2(nlon,nlat,nlev,ntimes))
allocate(PED_all(ntimes,nlev,nlat,nlon))
allocate(PED_all2(nlon,nlat,nlev,ntimes))
allocate(HALL_all(ntimes,nlev,nlat,nlon))
allocate(HALL_all2(nlon,nlat,nlev,ntimes))
allocate(Un_east(ntimes,nlev,nlat,nlon))
allocate(Un_east2(nlon,nlat,nlev,ntimes))
allocate(Un_north(ntimes,nlev,nlat,nlon))
allocate(Un_north2(nlon,nlat,nlev,ntimes))
allocate(Un_up(ntimes,nlev,nlat,nlon))
allocate(Un_up2(nlon,nlat,nlev,ntimes))
allocate(Ui_east(ntimes,nlev,nlat,nlon))
allocate(Ui_east2(nlon,nlat,nlev,ntimes))
allocate(Ui_north(ntimes,nlev,nlat,nlon))
allocate(Ui_north2(nlon,nlat,nlev,ntimes))
allocate(Ui_up(ntimes,nlev,nlat,nlon))
allocate(Ui_up2(nlon,nlat,nlev,ntimes))
allocate(Tn_all(ntimes,nlev,nlat,nlon))
allocate(Tn_all2(nlon,nlat,nlev,ntimes))
allocate(Te_all(ntimes,nlev,nlat,nlon))
allocate(Te_all2(nlon,nlat,nlev,ntimes))
allocate(Ti_all(ntimes,nlev,nlat,nlon))
allocate(Ti_all2(nlon,nlat,nlev,ntimes))
allocate(nO_all(ntimes,nlev,nlat,nlon))
allocate(nO_all2(nlon,nlat,nlev,ntimes))
allocate(nO2_all(ntimes,nlev,nlat,nlon))
allocate(nO2_all2(nlon,nlat,nlev,ntimes))
allocate(nN2_all(ntimes,nlev,nlat,nlon))
allocate(nN2_all2(nlon,nlat,nlev,ntimes))
allocate(NOp_all(ntimes,nlev,nlat,nlon))
allocate(NOp_all2(nlon,nlat,nlev,ntimes))
allocate(NO2p_all(ntimes,nlev,nlat,nlon))
allocate(NO2p_all2(nlon,nlat,nlev,ntimes))
allocate(NNOp_all(ntimes,nlev,nlat,nlon))
allocate(NNOp_all2(nlon,nlat,nlev,ntimes))
allocate(E_X(ntimes,nlev,nlat,nlon))
allocate(E_X2(nlon,nlat,nlev,ntimes))
allocate(E_Y(ntimes,nlev,nlat,nlon))
allocate(E_Y2(nlon,nlat,nlev,ntimes))
allocate(E_Z(ntimes,nlev,nlat,nlon))
allocate(E_Z2(nlon,nlat,nlev,ntimes))
allocate(BX_BMAG(ntimes,nlat,nlon))
allocate(BX_BMAG2(nlon,nlat,ntimes))
allocate(BY_BMAG(ntimes,nlat,nlon))
allocate(BY_BMAG2(nlon,nlat,ntimes))
allocate(BZ_BMAG(ntimes,nlat,nlon))
allocate(BZ_BMAG2(nlon,nlat,ntimes))
allocate(BMAG_G(ntimes,nlat,nlon))
allocate(BMAG2_G(nlon,nlat,ntimes))

!Memory Allocation for other variables
allocate(time_pf(ntimes))
allocate(zgf(ntimes,nlev,nlat,nlon))
allocate(zgout(nlon,nlat,nlev,ntimes))
allocate(mlat_qdf(ntimes,nlev,nlat,nlon))
allocate(mlt_qdf(ntimes,nlev,nlat,nlon))
allocate(BX_T(ntimes,nlev,nlat,nlon))
allocate(BY_T(ntimes,nlev,nlat,nlon))
allocate(BZ_T(ntimes,nlev,nlat,nlon))
allocate(BX_Tout(nlon,nlat,nlev,ntimes))
allocate(BY_Tout(nlon,nlat,nlev,ntimes))
allocate(BZ_Tout(nlon,nlat,nlev,ntimes))
allocate(BMAG_T(ntimes,nlev,nlat,nlon))
allocate(BMAG_Tout(ntimes,nlev,nlat,nlon))
allocate(E_X_si(ntimes,nlev,nlat,nlon))
allocate(E_X_si_out(nlon,nlat,nlev,ntimes))
allocate(E_Y_si(ntimes,nlev,nlat,nlon))
allocate(E_Y_si_out(nlon,nlat,nlev,ntimes))
allocate(E_Z_si(ntimes,nlev,nlat,nlon))
allocate(E_Z_si_out(nlon,nlat,nlev,ntimes))
allocate(Un_east_si(ntimes,nlev,nlat,nlon))
allocate(Un_east_si_out(nlon,nlat,nlev,ntimes))
allocate(Un_north_si(ntimes,nlev,nlat,nlon))
allocate(Un_north_si_out(nlon,nlat,nlev,ntimes))
allocate(Un_up_si(ntimes,nlev,nlat,nlon))
allocate(Un_up_si_out(nlon,nlat,nlev,ntimes))
allocate(Ui_east_si(ntimes,nlev,nlat,nlon))
allocate(Ui_east_si_out(nlon,nlat,nlev,ntimes))
allocate(Ui_north_si(ntimes,nlev,nlat,nlon))
allocate(Ui_north_si_out(nlon,nlat,nlev,ntimes))
allocate(Ui_up_si(ntimes,nlev,nlat,nlon))
allocate(Ui_up_si_out(nlon,nlat,nlev,ntimes))
allocate(nu_O2p_O2(ntimes,nlev,nlat,nlon))
allocate(nu_O2p_O2_out(nlon,nlat,nlev,ntimes))
allocate(nu_O2p_O(ntimes,nlev,nlat,nlon))
allocate(nu_O2p_O_out(nlon,nlat,nlev,ntimes))
allocate(nu_O2p_N2(ntimes,nlev,nlat,nlon))
allocate(nu_O2p_N2_out(nlon,nlat,nlev,ntimes))
allocate(nu_O2p(ntimes,nlev,nlat,nlon))
allocate(nu_O2p_out(nlon,nlat,nlev,ntimes))
allocate(nu_Op_O(ntimes,nlev,nlat,nlon))
allocate(nu_Op_O_out(nlon,nlat,nlev,ntimes))
allocate(nu_Op_O2(ntimes,nlev,nlat,nlon))
allocate(nu_Op_O2_out(nlon,nlat,nlev,ntimes))
allocate(nu_Op_N2(ntimes,nlev,nlat,nlon))
allocate(nu_Op_N2_out(nlon,nlat,nlev,ntimes))
allocate(nu_Op(ntimes,nlev,nlat,nlon))
allocate(nu_Op_out(nlon,nlat,nlev,ntimes))
allocate(nu_NOp_O2(ntimes,nlev,nlat,nlon))
allocate(nu_NOp_O2_out(nlon,nlat,nlev,ntimes))
allocate(nu_NOp_O(ntimes,nlev,nlat,nlon))
allocate(nu_NOp_O_out(nlon,nlat,nlev,ntimes))
allocate(nu_NOp_N2(ntimes,nlev,nlat,nlon))
allocate(nu_NOp_N2_out(nlon,nlat,nlev,ntimes))
allocate(nu_NOp(ntimes,nlev,nlat,nlon))
allocate(nu_NOp_out(nlon,nlat,nlev,ntimes))
allocate(nu_e(ntimes,nlev,nlat,nlon))
allocate(nu_e_out(nlon,nlat,nlev,ntimes))
allocate(nu_in(ntimes,nlev,nlat,nlon))
allocate(nu_in_out(nlon,nlat,nlev,ntimes))
allocate(Ohmic(ntimes,nlev,nlat,nlon))
allocate(Ohmicout(nlon,nlat,nlev,ntimes))
allocate(Joule(ntimes,nlev,nlat,nlon))
allocate(Jouleout(nlon,nlat,nlev,ntimes))
allocate(Frictional(ntimes,nlev,nlat,nlon))
allocate(Frictionalout(nlon,nlat,nlev,ntimes))
allocate(PED_calc(ntimes,nlev,nlat,nlon))
allocate(PED_calcout(nlon,nlat,nlev,ntimes))
allocate(HALL_calc(ntimes,nlev,nlat,nlon))
allocate(HALL_calcout(nlon,nlat,nlev,ntimes))
allocate(PAR_calc(ntimes,nlev,nlat,nlon))
allocate(PAR_calcout(nlon,nlat,nlev,ntimes))
allocate(CS_Op(ntimes,nlev,nlat,nlon))
allocate(CS_Opout(nlon,nlat,nlev,ntimes))
allocate(CS(ntimes,nlev,nlat,nlon))
allocate(CSout(nlon,nlat,nlev,ntimes))
allocate(CS_O2p(ntimes,nlev,nlat,nlon))
allocate(CS_O2pout(nlon,nlat,nlev,ntimes))
allocate(CS_NOp(ntimes,nlev,nlat,nlon))
allocate(CS_NOpout(nlon,nlat,nlev,ntimes))
allocate(Pedersen_current(ntimes,nlev,nlat,nlon))
allocate(Pedersen_current_out(nlon,nlat,nlev,ntimes))
allocate(Hall_current(ntimes,nlev,nlat,nlon))
allocate(Hall_current_out(nlon,nlat,nlev,ntimes))
allocate(Total_current(ntimes,nlev,nlat,nlon))
allocate(Total_current_out(nlon,nlat,nlev,ntimes))
allocate(Convenction_heating(ntimes,nlev,nlat,nlon))
allocate(Convenction_heating_out(nlon,nlat,nlev,ntimes))
allocate(Wind_heating(ntimes,nlev,nlat,nlon))
allocate(Wind_heating_out(nlon,nlat,nlev,ntimes))





if (myid==root) then
   ! LOG
   open(61,file=LOG,action='write',position='append')
   write(61,*) "LOG::LoadBalance and Memory Allocated"
   write(61,*) "==========================================================================="
   close(61)
   print*,"--------------------------------------------"
   print*, "Sanity Check:First file of each MPI Task:"
endif

print*, localfiles(1)
call MPI_Barrier(mpi_comm_world,mpierror)


!Main Loop
t0=MPI_Wtime()
do files=1,sndcnt

  !write(61,*) "LOG::MPI Task ",myid ,"Start with calculations of file",files,"out of",recvcnt
  ! print*, "LOG::MPI Task ",myid ,"Start with calculations of file:",localfiles(files),"out of",recvcnt


  ! rename output nc files
  name_index=len(trim(localfiles(files)))-3 ! count number of chars and put out the extension
  ncINfileNAMEsplited=localfiles(files)(1:name_index)! put out of name the extension
  ncOUTfileNAME=trim(ncINfileNAMEsplited)//"_JH_QD_AllVars"//".nc" !Output file name creation
  rnCmd='mv '//localfiles(files)//' '//ncOUTfileNAME !rename copied file
  ! print*, "rnCmd:", rnCmd
  print*, "LOG::MPI Task ",myid ,"Rename file",localfiles(files)
  call system (rnCmd)
  localfiles(files)=ncOUTfileNAME
  !Read File and Variables
  ! This is the name of the data file every processor will read each time
  ! ncerror= nf90_open(localfiles(files),NF90_NOWRITE,ncid)
  ncOUTerror= nf90_open(localfiles(files),nf90_write,ncOUTid)
  ncid=ncOUTid
  ncerror=nf90_inq_varid(ncid, "lat", varid)
  ncerror=nf90_get_var(ncid,varid,glat)
  call check_NetCDFvar_status(ncerror,"lat")

  ncerror=nf90_inq_varid(ncid, "lon", varid)
  ncerror=nf90_get_var(ncid,varid,glon)
  call check_NetCDFvar_status(ncerror,"lon")

  ncerror=nf90_inq_varid(ncid, "lev", varid)
  ncerror=nf90_get_var(ncid,varid,glev)
  call check_NetCDFvar_status(ncerror,"lev")

  ncerror=nf90_inq_varid(ncid, "time", varid)
  ncerror=nf90_get_var(ncid,varid,gtime)
  call check_NetCDFvar_status(ncerror,"time")

  ncerror=nf90_inq_varid(ncid, "Kp", varid)
  ncerror=nf90_get_var(ncid,varid,KP)
  call check_NetCDFvar_status(ncerror,"Kp")

  ncerror=nf90_inq_varid(ncid, "ZGMID", varid)
  ncerror=nf90_get_var(ncid,varid,zg2)
  call check_NetCDFvar_status(ncerror,"ZGMID")
  zg=reshape(zg2,(/ntimes,nlev,nlat,nlon/),ORDER=(/4,3,2,1/))


  ncerror=nf90_inq_varid(ncid, "ELECDEN", varid)
  ncerror=nf90_get_var(ncid,varid,Ne_all2)
  call check_NetCDFvar_status(ncerror,"ELECDEN")
  Ne_all=reshape(Ne_all2,(/ntimes,nlev,nlat,nlon/),ORDER=(/4,3,2,1/))


  ncerror=nf90_inq_varid(ncid, "SIGMA_PED", varid)
  ncerror=nf90_get_var(ncid,varid,PED_all2)
  call check_NetCDFvar_status(ncerror,"SIGMA_PED")
  PED_all=reshape(PED_all2,(/ntimes,nlev,nlat,nlon/),ORDER=(/4,3,2,1/))


  ncerror=nf90_inq_varid(ncid, "SIGMA_HAL", varid)
  ncerror=nf90_get_var(ncid,varid,HALL_all2)
  call check_NetCDFvar_status(ncerror,"SIGMA_HAL")
  HALL_all=reshape(HALL_all2,(/ntimes,nlev,nlat,nlon/),ORDER=(/4,3,2,1/))


  ncerror=nf90_inq_varid(ncid, "UN", varid)
  ncerror=nf90_get_var(ncid,varid,Un_east2)
  call check_NetCDFvar_status(ncerror,"UN")
  Un_east=reshape(Un_east2,(/ntimes,nlev,nlat,nlon/),ORDER=(/4,3,2,1/))


  ncerror=nf90_inq_varid(ncid, "VN", varid)
  ncerror=nf90_get_var(ncid,varid,Un_north2)
  call check_NetCDFvar_status(ncerror,"VN")
  Un_north=reshape(Un_north2,(/ntimes,nlev,nlat,nlon/),ORDER=(/4,3,2,1/))


  ncerror=nf90_inq_varid(ncid, "WN_lev", varid)
  ncerror=nf90_get_var(ncid,varid,Un_up2)
  call check_NetCDFvar_status(ncerror,"WN_lev")
  Un_up=reshape(Un_up2,(/ntimes,nlev,nlat,nlon/),ORDER=(/4,3,2,1/))


  ncerror=nf90_inq_varid(ncid, "Ui_lev", varid)
  ncerror=nf90_get_var(ncid,varid,Ui_east2)
  call check_NetCDFvar_status(ncerror,"Ui_lev")
  Ui_east=reshape(Ui_east2,(/ntimes,nlev,nlat,nlon/),ORDER=(/4,3,2,1/))

  ncerror=nf90_inq_varid(ncid, "Vi_lev", varid)
  ncerror=nf90_get_var(ncid,varid,Ui_north2)
  call check_NetCDFvar_status(ncerror,"Vi_lev")
  Ui_north=reshape(Ui_north2,(/ntimes,nlev,nlat,nlon/),ORDER=(/4,3,2,1/))

  ncerror=nf90_inq_varid(ncid, "Wi_lev", varid)
  ncerror=nf90_get_var(ncid,varid,Ui_up2)
  call check_NetCDFvar_status(ncerror,"Wi_lev")
  Ui_up=reshape(Ui_up2,(/ntimes,nlev,nlat,nlon/),ORDER=(/4,3,2,1/))

  ncerror=nf90_inq_varid(ncid, "TN", varid)
  ncerror=nf90_get_var(ncid,varid,Tn_all2)
  call check_NetCDFvar_status(ncerror,"TN")
  Tn_all=reshape(Tn_all2,(/ntimes,nlev,nlat,nlon/),ORDER=(/4,3,2,1/))

  ncerror=nf90_inq_varid(ncid, "TE", varid)
  ncerror=nf90_get_var(ncid,varid,Te_all2)
  call check_NetCDFvar_status(ncerror,"TE")
  Te_all=reshape(Te_all2,(/ntimes,nlev,nlat,nlon/),ORDER=(/4,3,2,1/))

  ncerror=nf90_inq_varid(ncid, "TI", varid)
  ncerror=nf90_get_var(ncid,varid,Ti_all2)
  call check_NetCDFvar_status(ncerror,"TI")
  Ti_all=reshape(Ti_all2,(/ntimes,nlev,nlat,nlon/),ORDER=(/4,3,2,1/))

  ncerror=nf90_inq_varid(ncid, "O_CM3", varid)
  ncerror=nf90_get_var(ncid,varid,nO_all2)
  call check_NetCDFvar_status(ncerror,"O_CM3")
  nO_all=reshape(nO_all2,(/ntimes,nlev,nlat,nlon/),ORDER=(/4,3,2,1/))

  ncerror=nf90_inq_varid(ncid, "O2_CM3", varid)
  ncerror=nf90_get_var(ncid,varid,nO2_all2)
  call check_NetCDFvar_status(ncerror,"O2_CM3")
  nO2_all=reshape(nO2_all2,(/ntimes,nlev,nlat,nlon/),ORDER=(/4,3,2,1/))

  ncerror=nf90_inq_varid(ncid, "N2_CM3", varid)
  ncerror=nf90_get_var(ncid,varid,nN2_all2)
  call check_NetCDFvar_status(ncerror,"N2_CM3")
  nN2_all=reshape(nN2_all2,(/ntimes,nlev,nlat,nlon/),ORDER=(/4,3,2,1/))

  ncerror=nf90_inq_varid(ncid, "OP", varid)
  ncerror=nf90_get_var(ncid,varid,NOp_all2)
  call check_NetCDFvar_status(ncerror,"Op")
  NOp_all=reshape(NOp_all2,(/ntimes,nlev,nlat,nlon/),ORDER=(/4,3,2,1/))


  ncerror=nf90_inq_varid(ncid, "O2P", varid)
  ncerror=nf90_get_var(ncid,varid,NO2p_all2)
  call check_NetCDFvar_status(ncerror,"O2P")
  NO2p_all=reshape(NO2p_all2,(/ntimes,nlev,nlat,nlon/),ORDER=(/4,3,2,1/))

  ncerror=nf90_inq_varid(ncid, "NOP_LAM", varid)
  ncerror=nf90_get_var(ncid,varid,NNOp_all2)
  call check_NetCDFvar_status(ncerror,"NOP_LAM")
  NNOp_all=reshape(NNOp_all2,(/ntimes,nlev,nlat,nlon/),ORDER=(/4,3,2,1/))

  ncerror=nf90_inq_varid(ncid, "EEX", varid)
  ncerror=nf90_get_var(ncid,varid,E_X2)
  call check_NetCDFvar_status(ncerror,"EEX")
  E_X=reshape(E_X2,(/ntimes,nlev,nlat,nlon/),ORDER=(/4,3,2,1/))


  ncerror=nf90_inq_varid(ncid, "EEY", varid)
  ncerror=nf90_get_var(ncid,varid,E_Y2)
  call check_NetCDFvar_status(ncerror,"EEY")
  E_Y=reshape(E_Y2,(/ntimes,nlev,nlat,nlon/),ORDER=(/4,3,2,1/))


  ncerror=nf90_inq_varid(ncid, "EEZ", varid)
  ncerror=nf90_get_var(ncid,varid,E_Z2)
  call check_NetCDFvar_status(ncerror,"EEX")
  E_Z=reshape(E_Z2,(/ntimes,nlev,nlat,nlon/),ORDER=(/4,3,2,1/))


  ncerror=nf90_inq_varid(ncid, "BX", varid)
  ncerror=nf90_get_var(ncid,varid,BX_BMAG2)
  call check_NetCDFvar_status(ncerror,"BX")
  BX_BMAG=reshape(BX_BMAG2,(/ntimes,nlat,nlon/),ORDER=(/3,2,1/))


  ncerror=nf90_inq_varid(ncid, "BY", varid)
  ncerror=nf90_get_var(ncid,varid,BY_BMAG2)
  call check_NetCDFvar_status(ncerror,"BY")
  BY_BMAG=reshape(BY_BMAG2,(/ntimes,nlat,nlon/),ORDER=(/3,2,1/))


  ncerror=nf90_inq_varid(ncid, "BZ", varid)
  ncerror=nf90_get_var(ncid,varid,BZ_BMAG2)
  call check_NetCDFvar_status(ncerror,"BZ")
  BZ_BMAG=reshape(BZ_BMAG2,(/ntimes,nlat,nlon/),ORDER=(/3,2,1/))


  ncerror=nf90_inq_varid(ncid, "BMAG", varid)
  ncerror=nf90_get_var(ncid,varid,BMAG2_G)
  call check_NetCDFvar_status(ncerror,"BMAG")
  BMAG_G=reshape(BMAG2_G,(/ntimes,nlat,nlon/),ORDER=(/3,2,1/))

  ! print*,"load complete"
  ! pause
  ! ncerror=nf90_close(ncid)
  ! print*, "file readed"
  !Timer
  t1=MPI_Wtime()
  ! for first file of every year
  ! do levcnt=1,nlev
  !    do latcnt=1,nlat
  !       do loncnt=1,nlon
  !         zgf(1,levcnt,latcnt,loncnt)=0
  !         mlat_qdf(1,levcnt,latcnt,loncnt)=0
  !         mlt_qdf(1,levcnt,latcnt,loncnt)=0
  !         Un_east_si(1,levcnt,latcnt,loncnt)=0
  !         Un_north_si(1,levcnt,latcnt,loncnt)=0
  !         Un_up_si(1,levcnt,latcnt,loncnt)=0
  !         Ui_east_si(1,levcnt,latcnt,loncnt)=0
  !         Ui_north_si(1,levcnt,latcnt,loncnt)=0
  !         Ui_up_si(1,levcnt,latcnt,loncnt)=0
  !         BX_T(1,levcnt,latcnt,loncnt)=0
  !         BY_T(1,levcnt,latcnt,loncnt)=0
  !         BZ_T(1,levcnt,latcnt,loncnt)=0
  !         BMAG_T(1,levcnt,latcnt,loncnt)=0
  !         E_X_si(1,levcnt,latcnt,loncnt)=0
  !         E_Y_si(1,levcnt,latcnt,loncnt)=0
  !         E_Z_si(1,levcnt,latcnt,loncnt)=0
  !         Ohmic(1,levcnt,latcnt,loncnt)=0
  !         Joule(1,levcnt,latcnt,loncnt)=0
  !         Frictional(1,levcnt,latcnt,loncnt)=0
  !         PED_calc(1,levcnt,latcnt,loncnt)=0
  !         HALL_calc(1,levcnt,latcnt,loncnt)=0
  !         PAR_calc(1,levcnt,latcnt,loncnt)=0
  !         CS_Op(1,levcnt,latcnt,loncnt)=0
  !         CS(1,levcnt,latcnt,loncnt)=0
  !         CS_O2p(1,levcnt,latcnt,loncnt)=0
  !         CS_NOp(1,levcnt,latcnt,loncnt)=0
  !         Pedersen_current(1,levcnt,latcnt,loncnt)=0
  !         Hall_current(1,levcnt,latcnt,loncnt)=0
  !         Total_current(1,levcnt,latcnt,loncnt)=0
  !         Convenction_heating(1,levcnt,latcnt,loncnt)=0
  !         Wind_heating(1,levcnt,latcnt,loncnt)=0
  !       end do
  !    end do
  ! end do
  ! Iterare Over Gridded Domain
  do timescnt=1,ntimes
  ! do timescnt=1,3
      ! print*,"timescnt:",timescnt
    ! Handle Time-MPI Independent
    timeg=gtime(timescnt)!gtime, minutes since 2015-03-15
    ! timeg=525600.0 !minutes since 2015-03-15
    !convert minutes to years, months, days,hours, mins,secs
    daysFRinit=int(timeg*60/86400)! days from initial date
    if (daysFRinit.gt.365) then !calculate if days are more than a year
      current_year=initial_yr+int(daysFRinit/365)
      daysFRinit=daysFRinit-int(daysFRinit/365)*365
      current_day=initial_day+daysFRinit
    else
      current_year=initial_yr
      current_day=initial_day+daysFRinit
    end if
    restMINS=mod(int(timeg)*60,86400)/60
    hrFRinit=int(restMINS/60)
    current_hour=initial_hour+hrFRinit
    minFRinit=mod(restMINS,60)
    current_mins=initial_mins+minFRinit
    date=current_year+(current_day*24*60*60+current_hour*60*60+current_mins*60+current_secs)/secsOFyr

    ! SetUp Apex
    call apex_setup(date,1000.0)
    do levcnt=1,nlev
       do latcnt=1,nlat
          do loncnt=1,nlon

            ! GEO coordinates of desired point
            alt_p=zg(timescnt,levcnt,latcnt,loncnt)/1.0e5 !cm to km
            zgf(timescnt,levcnt,latcnt,loncnt)=alt_p !Output file in km
            lat_p = glat(latcnt)   !deg
            ! print*, "lat::::::",lat_p
            lon_p = glon(loncnt)   !deg

            !Temperatures
            Tn=Tn_all(timescnt,levcnt,latcnt,loncnt)     !Kelvin
            Te=Te_all(timescnt,levcnt,latcnt,loncnt)     !Kelvin
            Ti=Ti_all(timescnt,levcnt,latcnt,loncnt)     !Kelvin

            !Geo2Geod
            call geo_lat2geod_lat(lat_p,geod_lat)
            !Apex
            call apex_mall(lat_p,lon_p,1000.0,alt_p, B,bhat,bmag,si,qdlon,malat,vmp,W,&
                           D,Be3,sim,d1,d2,d3,e1,e2,e3,qdlat,F,f1,f2,f3,g1,g2,g3,ier)
            mlat_qd=qdlat
            mlon_qd=qdlon
            call cofrm(date)
            iyr=current_year
            iday=current_day
            ihr=current_hour
            imn=current_mins
            sec=current_secs
            call subsol(iyr,iday,ihr,imn,sec,sbsllat,sbsllon)
            call magloctm(qdlon,sbsllat,sbsllon,colat,elon,mlt)
            mlt_qd=mlt

            mlat_qdf(timescnt,levcnt,latcnt,loncnt)=mlat_qd !Output file in deg
            mlt_qdf(timescnt,levcnt,latcnt,loncnt)=mlt_qd !Output file in deg

            !neutral Densities
            NO=nO_all(timescnt,levcnt,latcnt,loncnt)!cm^-3
            NO2=nO2_all(timescnt,levcnt,latcnt,loncnt)!cm^-3
            NN2=nN2_all(timescnt,levcnt,latcnt,loncnt)!cm^-3

            ! electron density
            Ne = Ne_all(timescnt,levcnt,latcnt,loncnt) * CubicCm2CubicM ! m^-3

            !ion Densities
            NOp=NOp_all(timescnt,levcnt,latcnt,loncnt) * CubicCm2CubicM  !m^-3
            NO2p=NO2p_all(timescnt,levcnt,latcnt,loncnt) * CubicCm2CubicM !m^-3
            NNOp=NNOp_all(timescnt,levcnt,latcnt,loncnt) * CubicCm2CubicM !m^-3

            !neutral wind
            Une=Un_east(timescnt,levcnt,latcnt,loncnt)
            Unn=Un_north(timescnt,levcnt,latcnt,loncnt)
            Unu=Un_up(timescnt,levcnt,latcnt,loncnt)

            Un_e = Une / 100  ! m/s
            Un_n = Unn / 100  ! m/s
            Un_u = Unu / 100  ! m/s

            Un_east_si(timescnt,levcnt,latcnt,loncnt)=Un_e
            Un_north_si(timescnt,levcnt,latcnt,loncnt)=Un_n
            Un_up_si(timescnt,levcnt,latcnt,loncnt)=Un_u

            !neutral wind in ENU in m/s
            UnENU=[Un_e,Un_n,Un_u]
            call enu_ecef(lat_p,lon_p,UnENU,Un) !Un in neutral velocity in ecef coordinates
            !1 ExB velocity
            Uie=Ui_east(timescnt,levcnt,latcnt,loncnt)!Output in si (m/s)
            Uin=Ui_north(timescnt,levcnt,latcnt,loncnt)!Output in si (m/s)
            Uiu=Ui_up(timescnt,levcnt,latcnt,loncnt)!Output in si (m/s)

            Ui_e = Uie / 100  ! m/sec
            Ui_n = Uin / 100  ! m/sec
            Ui_u = Uiu / 100  ! m/sec

            Ui_east_si(timescnt,levcnt,latcnt,loncnt)=Ui_e!Output in si (m/s)
            Ui_north_si(timescnt,levcnt,latcnt,loncnt)=Ui_n!Output in si (m/s)
            Ui_up_si(timescnt,levcnt,latcnt,loncnt)=Ui_u!Output in si (m/s)

            ! ExB drift in ENU in m/s
            UiENU=[Ui_e,Ui_n,Ui_u]
            call enu_ecef(lat_p,lon_p,UiENU,Ui) !Ui in neutral velocity in ecef coordinates
            ! Magnetic Field ENU
            Be=BX_BMAG(timescnt,latcnt,loncnt)*BMAG_G(timescnt,latcnt,loncnt)*1e-4! 1e-4 factor from Gouss to Tesla
            Bn=BY_BMAG(timescnt,latcnt,loncnt)*BMAG_G(timescnt,latcnt,loncnt)*1e-4! 1e-4 factor from Gouss to Tesla
            Bu=BZ_BMAG(timescnt,latcnt,loncnt)*BMAG_G(timescnt,latcnt,loncnt)*1e-4! 0.0001 factor from Gouss to Tesla
            Benu= [Be,Bn,Bu]
            call enu_ecef(lat_p,lon_p,Benu,Bxyz)
            ! Magnetic field unit vector
            bnorm=BMAG_G(timescnt,latcnt,loncnt)*1e-4 ! 0.0001 factor from Gauss  to Tesla
            b_unit(1)= Bxyz(1)/bnorm
            b_unit(2)= Bxyz(2)/bnorm
            b_unit(3)= Bxyz(3)/bnorm
            BX_i=Be / bnorm
            BY_i= Bn / bnorm
            BZ_i= Bu / bnorm

            BX_T(timescnt,levcnt,latcnt,loncnt)=Be !Output in si (Tesla)
            BY_T(timescnt,levcnt,latcnt,loncnt)=Bn !Output in si (Tesla)
            BZ_T(timescnt,levcnt,latcnt,loncnt)=Bu !Output in si (Tesla)

            BMAG_T(timescnt,levcnt,latcnt,loncnt)=bnorm !Output in si (Tesla)

            !neutral wind perpendicular to the magnetic field
            call CROSS(Un,b_unit,Unvert)
            Unvertx=Unvert(1)
            Unverty=Unvert(2)
            Unvertz=Unvert(3)
            Unvertmag=sqrt(Unvert(1) * Unvert(1) + Unvert(2) * Unvert(2) + Unvert(3) * Unvert(3))

            ! conductivities
            Ri = (Ti + Tn) / 2
            !nu_O2p=nu_O2p_O2+nu_O2p_O+nu_O2p_N2
            !1O2+-O2 collision frequency
            nu_O2p_O2(timescnt,levcnt,latcnt,loncnt) = 2.59 * 1e-11 * NO2 * sqrt(Ri) *((1 - 0.073 * log10(Ri)) ** 2)
            !O2+-O collision frequency
            nu_O2p_O(timescnt,levcnt,latcnt,loncnt) = 2.31 * 1e-10 * NO
            ! temp=2.31 * (10 ** (-10)) * NO
            !O2+-N2 collision frequency
            nu_O2p_N2(timescnt,levcnt,latcnt,loncnt) = 4.13 * 1e-10 * NN2
            ! print*, "NOP"
            ! print*,NO2,NO,NN2
            ! print*,nu_O2p_O(timescnt,levcnt,latcnt,loncnt),temp
            !O2+ collision frequency
            nu_O2p(timescnt,levcnt,latcnt,loncnt) = nu_O2p_O2(timescnt,levcnt,latcnt,loncnt) &
                                                    + nu_O2p_O(timescnt,levcnt,latcnt,loncnt)+ &
                                                    nu_O2p_N2(timescnt,levcnt,latcnt,loncnt)

            !nu_Op=nu_Op_O+nu_Op_O2+nu_Op_N
            !O+-O2 collision frequency
            nu_Op_O(timescnt,levcnt,latcnt,loncnt) = 3.67 * 1e-11 * NO * sqrt(Ri) *((1 - 0.064 * log10(Ri)) ** 2) * fb
            !O+-O collision frequency
            nu_Op_O2(timescnt,levcnt,latcnt,loncnt) = 6.64 * 1e-10 * NO2
            !O+-N2 collision frequency
            nu_Op_N2(timescnt,levcnt,latcnt,loncnt) = 6.82 * 1e-10 * NN2
            !O+ collision frequency
            nu_Op(timescnt,levcnt,latcnt,loncnt) = nu_Op_O(timescnt,levcnt,latcnt,loncnt) &
                                                  + nu_Op_O2(timescnt,levcnt,latcnt,loncnt)+ &
                                                  nu_Op_N2(timescnt,levcnt,latcnt,loncnt)


            !nu_NOp=nu_NOp_O2+nu_NOp_O+nu_NOp_N2
            !NO+-O2 collision frequency
            nu_NOp_O2(timescnt,levcnt,latcnt,loncnt)=(4.27 * 1e-10 * NO2)
            !NO+-O collision frequency
            nu_NOp_O(timescnt,levcnt,latcnt,loncnt)=(2.44 * 1e-10 * NO)
            !NO+-N2 collision frequency
            nu_NOp_N2(timescnt,levcnt,latcnt,loncnt)=(4.34 * 1e-10 * NN2)
            nu_NOp(timescnt,levcnt,latcnt,loncnt)=nu_NOp_O2(timescnt,levcnt,latcnt,loncnt) &
                                                  + nu_NOp_O(timescnt,levcnt,latcnt,loncnt) &
                                                  + nu_NOp_N2(timescnt,levcnt,latcnt,loncnt)
            nu_in(timescnt,levcnt,latcnt,loncnt)=nu_O2p(timescnt,levcnt,latcnt,loncnt)+&
                                                nu_Op(timescnt,levcnt,latcnt,loncnt)+&
                                                nu_NOp(timescnt,levcnt,latcnt,loncnt)
            !nu_e
            nu_e(timescnt,levcnt,latcnt,loncnt) = 2.33 * 1e-11 * NN2 * Te * (1 - 1.21 * 1e-4 * Te)  &
                                                  + 1.82 * 1e-10 * NO2 * sqrt(Te) * (1 + (3.6 * 1e-2 * sqrt(Te)))  &
                                                  + 8.9 * 1e-11 * NO *sqrt(Te) * (1 + (5.7 * (1e-4 * Te)))

            ! GyroFrequencies
            Omega_e = q_e * bnorm / m_e
            Omega_Op = q_e * bnorm / m_Oplus
            Omega_O2p = q_e * bnorm / m_O2plus
            Omega_NOp = q_e * bnorm / m_NOplus

            ! Ratios
            ratio_e = nu_e(timescnt,levcnt,latcnt,loncnt) / Omega_e
            ratio_Op = nu_Op(timescnt,levcnt,latcnt,loncnt) / Omega_Op
            ratio_O2p = nu_O2p(timescnt,levcnt,latcnt,loncnt) / Omega_O2p
            ratio_NOP=nu_NOp(timescnt,levcnt,latcnt,loncnt) / Omega_NOp

            ! Pedersen Conductivity
            Ne=NOp+NO2p+NNOp !Force Charge Neutrality
            alpha = (NOp) * (ratio_Op /(1 + ratio_Op ** 2))
            beta = (Ne) * (ratio_e /(1 + ratio_e ** 2))
            gamma = (NO2p) * (ratio_O2p /(1 + ratio_O2p ** 2))
            delta = (NNOp) * (ratio_NOP /(1 + ratio_NOP ** 2))
            ! #############################################################
            PED_calc(timescnt,levcnt,latcnt,loncnt) = (q_e / bnorm) * (alpha + beta + gamma +delta)

            ! Hall Conductivity
            alpha1 = (NOp) * (1 / (1 + ratio_Op ** 2))

            beta1 = (Ne) * (1 / (1 + ratio_e ** 2))
            gamma1 = (NO2p) * (1 / (1 + ratio_O2p ** 2))
            delta1 = (NNOp) * (1 / (1 + ratio_NOp ** 2))
            ! ############################################################
            HALL_calc(timescnt,levcnt,latcnt,loncnt) = (q_e / bnorm) * (beta1 - alpha1 - gamma1-delta1)

            ! Parallel Conductivity
            PAR_calc(timescnt,levcnt,latcnt,loncnt) = ((Ne) * q_e * q_e) /(m_e * nu_e(timescnt,levcnt,latcnt,loncnt))

            !Electric Field perpendicular to the magnetic field
            Evertx=E_X(timescnt,levcnt,latcnt,loncnt)*100  !Output file in V/m
            Everty=E_Y(timescnt,levcnt,latcnt,loncnt)*100    !Output file in V/m
            Evertz=E_Z(timescnt,levcnt,latcnt,loncnt)*100    !Output file in V/m
            E_X_si(timescnt,levcnt,latcnt,loncnt)=Evertx !Output file in V/m
            E_Y_si(timescnt,levcnt,latcnt,loncnt)=Everty !Output file in V/m
            E_Z_si(timescnt,levcnt,latcnt,loncnt)=Evertz !Output file in V/m
            Efieldenu(1)=Evertx
            Efieldenu(2)=Everty
            Efieldenu(3)=Evertz

            call enu_ecef(lat_p,lon_p,Efieldenu,Efieldecef)

            call CROSS(Efieldecef,b_unit,CrossResult)
            Evert=CrossResult
            call enu_ecef(lat_p,lon_p,Unvert,Unvertecef)
            call CROSS(Unvertecef,Bxyz,CrossResult)
            Estar=Evert+CrossResult
            Estarnorm=sqrt(Estar(1)*Estar(1)+Estar(2)*Estar(2)+Estar(3)*Estar(3))
            sigmaPED_file=PED_all(timescnt,levcnt,latcnt,loncnt)
            ! Ohmic(timescnt,levcnt,latcnt,loncnt)=sigmaPED_file*Estarnorm*Estarnorm
            Ohmic(timescnt,levcnt,latcnt,loncnt)=PED_calc(timescnt,levcnt,latcnt,loncnt)*&
                                                  dot_product(Estar,Estar)!*Estarnorm*Estarnorm


            ! Ion Velocity in neutral reference frame (formula derived through ion momentum equation)
            call CROSS(Estar,b_unit,CrossResult)
            Vi_star = (nu_Op(timescnt,levcnt,latcnt,loncnt) * Omega_Op * Estar + Omega_Op ** 2 * CrossResult)&
                        /(bnorm * (nu_Op(timescnt,levcnt,latcnt,loncnt) ** 2 +Omega_Op ** 2))


            Vi_starnorm = sqrt(Vi_star(1) * Vi_star(1) + Vi_star(2) * Vi_star(2) + Vi_star(3) * Vi_star(3))
            Joule(timescnt,levcnt,latcnt,loncnt)=q_e*Ne*dot_product(Vi_star,Estar)
            Frictional(timescnt,levcnt,latcnt,loncnt) = MOmass * Ne *nu_Op(timescnt,levcnt,latcnt,loncnt) &
                                                      * Vi_starnorm * Vi_starnorm

            ! ############################# Cross sections #############################################
            N_neutral = (NN2 + NO +NO2)*CubicCm2CubicM  ! in m^-3
            MnN2=(mN2)/(NA*1000)
            M_ionSUM=MO+MO2+MNO

            CS_Op(timescnt,levcnt,latcnt,loncnt) = (nu_Op(timescnt,levcnt,latcnt,loncnt) / N_neutral)/ &
                                                    sqrt((2 * boltzmann * Ti) / MO)

            CS(timescnt,levcnt,latcnt,loncnt)=(nu_in(timescnt,levcnt,latcnt,loncnt)/N_neutral)/&
                                              (sqrt((2*boltzmann*Ti)/M_ionSUM))


            CS_O2p(timescnt,levcnt,latcnt,loncnt)=(nu_O2p(timescnt,levcnt,latcnt,loncnt)/N_neutral)/(sqrt((2*boltzmann*Ti)/MO2))


            CS_NOp(timescnt,levcnt,latcnt,loncnt)=(nu_O2p(timescnt,levcnt,latcnt,loncnt)/N_neutral)/(sqrt((2*boltzmann*Ti)/MNO))



            ! ############################### Currents calculation ####################################
            ! Jpedersen = sigmaPED_file * Estar
            Jpedersen = PED_calc(timescnt,levcnt,latcnt,loncnt) * Estar
            Pedersen_current(timescnt,levcnt,latcnt,loncnt) = sqrt(Jpedersen(1) ** 2 + Jpedersen(2) ** 2 + Jpedersen(3) ** 2)

            call CROSS(b_unit,Estar,CrossResult)
            sigmaHALL_file=HALL_all(timescnt,levcnt,latcnt,loncnt)
            ! Jhall = sigmaHALL_file *CrossResult
            Jhall = HALL_calc(timescnt,levcnt,latcnt,loncnt) *CrossResult
            Hall_current(timescnt,levcnt,latcnt,loncnt) = sqrt(Jhall(1) ** 2 + Jhall(2) ** 2 + Jhall(3) ** 2)

            Total_current(timescnt,levcnt,latcnt,loncnt) = sqrt((Jpedersen(1) + Jhall(1)) ** 2 + (Jpedersen(2) + Jhall(2)) ** 2 &
                                          +(Jpedersen(3) + Jhall(3)) ** 2)

            !Magnetic Parameters
            dipmag=atan(Bxyz(3)/sqrt(Bxyz(2)**2+Bxyz(1)**2))
            sndip=sin(dipmag)
            decmag=atan2(Bxyz(2),Bxyz(1))
            sndec=sin(decmag)
            sn2dec=sndec**2
            csdec=cos(decmag)
            tm1=sn2dec+(1-sn2dec)*sndip**2
            tm2=(1-sn2dec)+sn2dec*sndip**2
            tm3=sndec*csdec*(1-sndip**2)


            ! Convenction_heating(timescnt,levcnt,latcnt,loncnt)=sigmaPED_file*dot_product(Evert,Evert)
            Convenction_heating(timescnt,levcnt,latcnt,loncnt)=PED_calc(timescnt,levcnt,latcnt,loncnt)*&
                                                                dot_product(Evert,Evert)

            call CROSS(Unvertecef,Bxyz,CrossResult)
            UNxB=CrossResult
            UNxBnorm=sqrt(UNxB(1)*UNxB(1)+UNxB(2)*UNxB(2)+UNxB(3)*UNxB(3))
            call CROSS(Evert,Bxyz,CrossResult)
            ExB=CrossResult
            ! Wind_heating(timescnt,levcnt,latcnt,loncnt)=sigmaPED_file*UNxBnorm*UNxBnorm-2*sigmaPED_file*dot_product(Un,ExB)
            ! Wind_heating(timescnt,levcnt,latcnt,loncnt)=PED_calc(timescnt,levcnt,latcnt,loncnt)*&
            !                                             UNxBnorm*UNxBnorm-2*sigmaPED_file*dot_product(Un,ExB)
            ! temp1=tm2*Un(1)**2+tm1*Un(2)**2
            ! temp2=2*tm3*Un(1)*Un(2)
            ! temp3=2*(Un(1)*Ui(1)+Un(2)*Ui(2))
            Wind_heating(timescnt,levcnt,latcnt,loncnt)=PED_calc(timescnt,levcnt,latcnt,loncnt)*(bnorm**2)*&
            (tm2*Un(1)**2+tm1*Un(2)**2-2*tm3*Un(1)*Un(2)-2*(Un(1)*Ui(1)+Un(2)*Ui(2)))

            ! if (lat_p.eq.63.75) then
            !   if (lon_p.eq.(-117.5)) then
            !     print*,levcnt-1
            !     print*,"Pedersen",PED_calc(timescnt,levcnt,latcnt,loncnt)
            !     print*,"bnorm",bnorm
            !     print*, tm1,tm2,tm3,Un(1),Un(2),Un(3)
            !     print*,Ui(1),Ui(2)
            !     print*,"temp1",temp1
            !     print*,"temp2",temp2
            !     print*,"temp3",temp3
            !     print*,Wind_heating(timescnt,levcnt,latcnt,loncnt)
            !     print*,"-----------------------"
            !   end if
            ! end if

         enddo
       enddo
     enddo
  enddo
  ! print*, "LOG::MPI Task ",myid ,"End with calculations of file:",localfiles(files)

  !Smooth out 2 last levels
  Ohmic(:,nlev,:,:)=Ohmic(:,nlev-2,:,:)
  Ohmic(:,nlev-1,:,:)=Ohmic(:,nlev-2,:,:)



  ncOUTerror = nf90_redef(ncOUTid) !define mode

  ncOUTerror=nf90_inq_dimid(ncOUTid, "time", NCntimesID)
  ncOUTerror=nf90_inq_dimid(ncOUTid, "lev", NCnlevID)
  ncOUTerror=nf90_inq_dimid(ncOUTid, "lat", NCnlatID)
  ncOUTerror=nf90_inq_dimid(ncOUTid, "lon", NCnlonID)
  ! ncOUTerror=nf90_inquire_dimension(ncOUTid, NCntimesID, len=times_dmind)
  ! ncOUTerror=nf90_inquire_dimension(ncOUTid, NCnlevID, len=lev_dmind)
  ! ncOUTerror=nf90_inquire_dimension(ncOUTid, NCnlatID, len=lat_dmind)
  ! ncOUTerror=nf90_inquire_dimension(ncOUTid, NCnlonID, len=lon_dmind)

  OUTDimids1=(/ NCnlonID,NCnlatID,NCnlevID,NCntimesID/)
  OUTDimids2=(/ NCnlonID,NCnlatID,NCntimesID/)

  ncOUTerror= nf90_def_var(ncOUTid, "Ohmic", NF90_FLOAT, OUTDimids1, varOUTid1)
  ncOUTerror= nf90_def_var(ncOUTid, "mlt_qdf", NF90_FLOAT, OUTDimids1, varOUTid2)
  ncOUTerror= nf90_def_var(ncOUTid, "mlat_qdf", NF90_FLOAT, OUTDimids1, varOUTid3)
  ncOUTerror= nf90_def_var(ncOUTid, "BX_si", NF90_FLOAT, OUTDimids1, varOUTid4)
  ncOUTerror= nf90_def_var(ncOUTid, "BY_si", NF90_FLOAT, OUTDimids1, varOUTid5)
  ncOUTerror= nf90_def_var(ncOUTid, "BZ_si", NF90_FLOAT, OUTDimids1, varOUTid6)
  ncOUTerror= nf90_def_var(ncOUTid, "EEX_si", NF90_FLOAT, OUTDimids1, varOUTid7)
  ncOUTerror= nf90_def_var(ncOUTid, "EEY_si", NF90_FLOAT, OUTDimids1, varOUTid8)
  ncOUTerror= nf90_def_var(ncOUTid, "EEZ_si", NF90_FLOAT, OUTDimids1, varOUTid9)
  ncOUTerror= nf90_def_var(ncOUTid, "UN_si", NF90_FLOAT, OUTDimids1, varOUTid10)
  ncOUTerror= nf90_def_var(ncOUTid, "VN_si", NF90_FLOAT, OUTDimids1, varOUTid11)
  ncOUTerror= nf90_def_var(ncOUTid, "WN_si", NF90_FLOAT, OUTDimids1, varOUTid12)
  ncOUTerror= nf90_def_var(ncOUTid, "Ui_lev_si", NF90_FLOAT, OUTDimids1, varOUTid13)
  ncOUTerror= nf90_def_var(ncOUTid, "Vi_lev_si", NF90_FLOAT, OUTDimids1, varOUTid14)
  ncOUTerror= nf90_def_var(ncOUTid, "Wi_lev_si", NF90_FLOAT, OUTDimids1, varOUTid15)
  ncOUTerror= nf90_def_var(ncOUTid, "BMAG_si", NF90_FLOAT, OUTDimids1, varOUTid16)
  ncOUTerror= nf90_def_var(ncOUTid, "Joule", NF90_FLOAT, OUTDimids1, varOUTid17)
  ncOUTerror= nf90_def_var(ncOUTid, "Frictional", NF90_FLOAT, OUTDimids1, varOUTid18)
  ncOUTerror= nf90_def_var(ncOUTid, "Pedersen_calc", NF90_FLOAT, OUTDimids1, varOUTid19)
  ncOUTerror= nf90_def_var(ncOUTid, "Hall_calc", NF90_FLOAT, OUTDimids1, varOUTid20)
  ncOUTerror= nf90_def_var(ncOUTid, "Parallel_calc", NF90_FLOAT, OUTDimids1, varOUTid21)
  ncOUTerror= nf90_def_var(ncOUTid, "Op_CrossSection", NF90_FLOAT, OUTDimids1, varOUTid22)
  ncOUTerror= nf90_def_var(ncOUTid, "O2p_CrossSection", NF90_FLOAT, OUTDimids1, varOUTid33)
  ncOUTerror= nf90_def_var(ncOUTid, "NOp_CrossSection", NF90_FLOAT, OUTDimids1, varOUTid34)
  ncOUTerror= nf90_def_var(ncOUTid, "Vin_CrossSection", NF90_FLOAT, OUTDimids1, varOUTid35)
  ncOUTerror= nf90_def_var(ncOUTid, "Pedersen_Current", NF90_FLOAT, OUTDimids1, varOUTid23)
  ncOUTerror= nf90_def_var(ncOUTid, "Hall_Current", NF90_FLOAT, OUTDimids1, varOUTid24)
  ncOUTerror= nf90_def_var(ncOUTid, "Total_Current", NF90_FLOAT, OUTDimids1, varOUTid25)
  ncOUTerror= nf90_def_var(ncOUTid, "VOp", NF90_FLOAT, OUTDimids1, varOUTid26)
  ncOUTerror= nf90_def_var(ncOUTid, "VO2p", NF90_FLOAT, OUTDimids1, varOUTid27)
  ncOUTerror= nf90_def_var(ncOUTid, "VNOp", NF90_FLOAT, OUTDimids1, varOUTid28)
  ncOUTerror= nf90_def_var(ncOUTid, "Ven", NF90_FLOAT, OUTDimids1, varOUTid29)
  ncOUTerror= nf90_def_var(ncOUTid, "Convection_heating", NF90_FLOAT, OUTDimids1, varOUTid30)
  ncOUTerror= nf90_def_var(ncOUTid, "Wind_heating", NF90_FLOAT, OUTDimids1, varOUTid31)
  ncOUTerror= nf90_def_var(ncOUTid, "Vin", NF90_FLOAT, OUTDimids1, varOUTid32)
  ncOUTerror=nf90_put_att(ncOUTid,varOUTid1,"units","W/m^3")
  ncOUTerror=nf90_put_att(ncOUTid,varOUTid2,"long_name" , "magnetic local time")
  ncOUTerror=nf90_put_att(ncOUTid,varOUTid2,"units" , "degrees")
  ncOUTerror=nf90_put_att(ncOUTid,varOUTid3,"long_name" , "quasi dipole magnetic latitude")
  ncOUTerror=nf90_put_att(ncOUTid,varOUTid3,"units" , "degrees")
  ncOUTerror=nf90_put_att(ncOUTid,varOUTid4,"long_name" , "Eastward component of magnetic field ")
  ncOUTerror=nf90_put_att(ncOUTid,varOUTid4,"units" , "Tesla")
  ncOUTerror=nf90_put_att(ncOUTid,varOUTid5,"long_name" , "Northward component of magnetic field ")
  ncOUTerror=nf90_put_att(ncOUTid,varOUTid5,"units" , "Tesla")
  ncOUTerror=nf90_put_att(ncOUTid,varOUTid6,"long_name" , "Upward component of magnetic field ")
  ncOUTerror=nf90_put_att(ncOUTid,varOUTid6,"units" , "Tesla")
  ncOUTerror=nf90_put_att(ncOUTid,varOUTid16,"long_name" , "magnetic field magnitude")
  ncOUTerror=nf90_put_att(ncOUTid,varOUTid16,"units" , "Tesla")
  ncOUTerror=nf90_put_att(ncOUTid,varOUTid7,"long_name" , "Ex ")
  ncOUTerror=nf90_put_att(ncOUTid,varOUTid7,"units" , "V/m")
  ncOUTerror=nf90_put_att(ncOUTid,varOUTid8,"long_name" , "Ey ")
  ncOUTerror=nf90_put_att(ncOUTid,varOUTid8,"units" , "V/m")
  ncOUTerror=nf90_put_att(ncOUTid,varOUTid9,"long_name" , "Ez ")
  ncOUTerror=nf90_put_att(ncOUTid,varOUTid9,"units" , "V/m")
  ncOUTerror=nf90_put_att(ncOUTid,varOUTid10,"long_name" , "NEUTRAL ZONAL WIND (+EAST) ")
  ncOUTerror=nf90_put_att(ncOUTid,varOUTid10,"units" , "m/sec")
  ncOUTerror=nf90_put_att(ncOUTid,varOUTid11,"long_name" , "NEUTRAL MERIDIONAL WIND (+NORTH) ")
  ncOUTerror=nf90_put_att(ncOUTid,varOUTid11,"units" , "m/sec")
  ncOUTerror=nf90_put_att(ncOUTid,varOUTid12,"long_name" , "Neutral Vertical Wind (plus up)")
  ncOUTerror=nf90_put_att(ncOUTid,varOUTid12,"units" , "m/sec")
  ncOUTerror=nf90_put_att(ncOUTid,varOUTid13,"units" , "m/sec")
  ncOUTerror=nf90_put_att(ncOUTid,varOUTid14,"units" , "m/sec")
  ncOUTerror=nf90_put_att(ncOUTid,varOUTid15,"units" , "m/sec")
  ncOUTerror=nf90_put_att(ncOUTid,varOUTid17,"units","W/m^3")
  ncOUTerror=nf90_put_att(ncOUTid,varOUTid18,"units","W/m^3")
  ncOUTerror=nf90_put_att(ncOUTid,varOUTid19,"units","S/m")
  ncOUTerror=nf90_put_att(ncOUTid,varOUTid19,"Description","Pedersen Conductivity Calculated from Richmond equations")
  ncOUTerror=nf90_put_att(ncOUTid,varOUTid20,"units","S/m")
  ncOUTerror=nf90_put_att(ncOUTid,varOUTid20,"Description","Hall Conductivity Calculated from Richmond equations")
  ncOUTerror=nf90_put_att(ncOUTid,varOUTid21,"units","S/m")
  ncOUTerror=nf90_put_att(ncOUTid,varOUTid21,"Description","Parallel Conductivity Calculated from Richmond equations")
  ! ncOUTerror=nf90_put_att(ncOUTid,varOUTid22,"units","cm^2")
  ncOUTerror=nf90_put_att(ncOUTid,varOUTid22,"units","m^2")
  ncOUTerror=nf90_put_att(ncOUTid,varOUTid33,"units","m^2")
  ncOUTerror=nf90_put_att(ncOUTid,varOUTid34,"units","m^2")
  ncOUTerror=nf90_put_att(ncOUTid,varOUTid35,"units","m^2")
  ncOUTerror=nf90_put_att(ncOUTid,varOUTid23,"units","A/m^2")
  ncOUTerror=nf90_put_att(ncOUTid,varOUTid24,"units","A/m^2")
  ncOUTerror=nf90_put_att(ncOUTid,varOUTid25,"units","A/m^2")
  ncOUTerror=nf90_put_att(ncOUTid,varOUTid26,"units","sec^-1")
  ncOUTerror=nf90_put_att(ncOUTid,varOUTid26,"Description","VOp Collision frequency")
  ncOUTerror=nf90_put_att(ncOUTid,varOUTid27,"units","sec^-1")
  ncOUTerror=nf90_put_att(ncOUTid,varOUTid27,"Description","VO2p Collision frequency")
  ncOUTerror=nf90_put_att(ncOUTid,varOUTid28,"units","sec^-1")
  ncOUTerror=nf90_put_att(ncOUTid,varOUTid28,"Description","VNOp Collision frequency")
  ncOUTerror=nf90_put_att(ncOUTid,varOUTid29,"units","sec^-1")
  ncOUTerror=nf90_put_att(ncOUTid,varOUTid29,"Description","Ve Collision frequency")
  ncOUTerror=nf90_put_att(ncOUTid,varOUTid30,"units","W/m^3")
  ncOUTerror=nf90_put_att(ncOUTid,varOUTid31,"units","W/m^3")
  ncOUTerror=nf90_put_att(ncOUTid,varOUTid32,"units","sec^-1")
  ncOUTerror=nf90_put_att(ncOUTid,varOUTid32,"Description","Ve Collision frequency")
  ncOUTerror=nf90_enddef(ncOUTid)
  Ohmicout=reshape(Ohmic,(/nlon,nlat,nlev,ntimes/),ORDER=(/4,3,2,1/))
  mlt_qdfout=reshape(mlt_qdf,(/nlon,nlat,nlev,ntimes/),ORDER=(/4,3,2,1/))
  mlat_qdfout=reshape(mlat_qdf,(/nlon,nlat,nlev,ntimes/),ORDER=(/4,3,2,1/))
  BX_Tout=reshape(BX_T,(/nlon,nlat,nlev,ntimes/),ORDER=(/4,3,2,1/))
  BY_Tout=reshape(BY_T,(/nlon,nlat,nlev,ntimes/),ORDER=(/4,3,2,1/))
  BZ_Tout=reshape(BZ_T,(/nlon,nlat,nlev,ntimes/),ORDER=(/4,3,2,1/))
  BMAG_Tout=reshape(BMAG_T,(/nlon,nlat,nlev,ntimes/),ORDER=(/4,3,2,1/))
  E_X_si_out=reshape(E_X_si,(/nlon,nlat,nlev,ntimes/),ORDER=(/4,3,2,1/))
  E_Y_si_out=reshape(E_Y_si,(/nlon,nlat,nlev,ntimes/),ORDER=(/4,3,2,1/))
  E_Z_si_out=reshape(E_Z_si,(/nlon,nlat,nlev,ntimes/),ORDER=(/4,3,2,1/))
  Un_east_si_out=reshape(Un_east_si,(/nlon,nlat,nlev,ntimes/),ORDER=(/4,3,2,1/))
  Un_north_si_out=reshape(Un_north_si,(/nlon,nlat,nlev,ntimes/),ORDER=(/4,3,2,1/))
  Un_up_si_out=reshape(Un_up_si,(/nlon,nlat,nlev,ntimes/),ORDER=(/4,3,2,1/))
  Ui_east_si_out=reshape(Ui_east_si,(/nlon,nlat,nlev,ntimes/),ORDER=(/4,3,2,1/))
  Ui_north_si_out=reshape(Ui_north_si,(/nlon,nlat,nlev,ntimes/),ORDER=(/4,3,2,1/))
  Ui_up_si_out=reshape(Ui_up_si,(/nlon,nlat,nlev,ntimes/),ORDER=(/4,3,2,1/))
  Jouleout=reshape(Joule,(/nlon,nlat,nlev,ntimes/),ORDER=(/4,3,2,1/))
  Frictionalout=reshape(Frictional,(/nlon,nlat,nlev,ntimes/),ORDER=(/4,3,2,1/))
  PED_calcout=reshape(PED_calc,(/nlon,nlat,nlev,ntimes/),ORDER=(/4,3,2,1/))
  HALL_calcout=reshape(HALL_calc,(/nlon,nlat,nlev,ntimes/),ORDER=(/4,3,2,1/))
  PAR_calcout=reshape(PAR_calc,(/nlon,nlat,nlev,ntimes/),ORDER=(/4,3,2,1/))
  CS_Opout=reshape(CS_Op,(/nlon,nlat,nlev,ntimes/),ORDER=(/4,3,2,1/))
  CSout=reshape(CS,(/nlon,nlat,nlev,ntimes/),ORDER=(/4,3,2,1/))
  CS_O2pout=reshape(CS_O2p,(/nlon,nlat,nlev,ntimes/),ORDER=(/4,3,2,1/))
  CS_NOpout=reshape(CS_NOp,(/nlon,nlat,nlev,ntimes/),ORDER=(/4,3,2,1/))
  Pedersen_current_out=reshape(Pedersen_current,(/nlon,nlat,nlev,ntimes/),ORDER=(/4,3,2,1/))
  Hall_current_out=reshape(Hall_current,(/nlon,nlat,nlev,ntimes/),ORDER=(/4,3,2,1/))
  Total_current_out=reshape(Total_current,(/nlon,nlat,nlev,ntimes/),ORDER=(/4,3,2,1/))
  nu_Op_out=reshape(nu_Op,(/nlon,nlat,nlev,ntimes/),ORDER=(/4,3,2,1/))
  nu_O2p_out=reshape(nu_O2p,(/nlon,nlat,nlev,ntimes/),ORDER=(/4,3,2,1/))
  nu_NOp_out=reshape(nu_NOp,(/nlon,nlat,nlev,ntimes/),ORDER=(/4,3,2,1/))
  nu_e_out=reshape(nu_e,(/nlon,nlat,nlev,ntimes/),ORDER=(/4,3,2,1/))
  nu_in_out=reshape(nu_in,(/nlon,nlat,nlev,ntimes/),ORDER=(/4,3,2,1/))
  Convenction_heating_out=reshape(Convenction_heating,(/nlon,nlat,nlev,ntimes/),ORDER=(/4,3,2,1/))
  Wind_heating_out=reshape(Wind_heating,(/nlon,nlat,nlev,ntimes/),ORDER=(/4,3,2,1/))

  ncOUTerror=nf90_put_var(ncOUTid, varOUTid1, Ohmicout)
  ncOUTerror=nf90_put_var(ncOUTid,  varOUTid2, mlt_qdfout)
  ncOUTerror=nf90_put_var(ncOUTid,  varOUTid3, mlat_qdfout)
  ncOUTerror=nf90_put_var(ncOUTid,  varOUTid4, BX_Tout)
  ncOUTerror=nf90_put_var(ncOUTid,  varOUTid5, BY_Tout)
  ncOUTerror=nf90_put_var(ncOUTid,  varOUTid6, BZ_Tout)
  ncOUTerror=nf90_put_var(ncOUTid,  varOUTid7, E_X_si_out)
  ncOUTerror=nf90_put_var(ncOUTid,  varOUTid8, E_Y_si_out)
  ncOUTerror=nf90_put_var(ncOUTid,  varOUTid9, E_Z_si_out)
  ncOUTerror=nf90_put_var(ncOUTid,  varOUTid10, Un_east_si_out)
  ncOUTerror=nf90_put_var(ncOUTid,  varOUTid11, Un_north_si_out)
  ncOUTerror=nf90_put_var(ncOUTid,  varOUTid12, Un_up_si_out)
  ncOUTerror=nf90_put_var(ncOUTid,  varOUTid13, Ui_east_si_out)
  ncOUTerror=nf90_put_var(ncOUTid,  varOUTid14, Ui_north_si_out)
  ncOUTerror=nf90_put_var(ncOUTid,  varOUTid15, Ui_up_si_out)
  ncOUTerror=nf90_put_var(ncOUTid,  varOUTid16, BMAG_Tout)
  ncOUTerror=nf90_put_var(ncOUTid,  varOUTid17, Jouleout)
  ncOUTerror=nf90_put_var(ncOUTid,  varOUTid18, Frictionalout)
  ncOUTerror=nf90_put_var(ncOUTid,  varOUTid19, PED_calcout)
  ncOUTerror=nf90_put_var(ncOUTid,  varOUTid20, HALL_calcout)
  ncOUTerror=nf90_put_var(ncOUTid,  varOUTid21, PAR_calcout)
  ncOUTerror=nf90_put_var(ncOUTid,  varOUTid22, CS_Opout)
  ncOUTerror=nf90_put_var(ncOUTid,  varOUTid33, CS_O2pout)
  ncOUTerror=nf90_put_var(ncOUTid,  varOUTid34, CS_NOpout)
  ncOUTerror=nf90_put_var(ncOUTid,  varOUTid35, CSout)
  ncOUTerror=nf90_put_var(ncOUTid,  varOUTid23, Pedersen_current_out)
  ncOUTerror=nf90_put_var(ncOUTid,  varOUTid24, Hall_current_out)
  ncOUTerror=nf90_put_var(ncOUTid,  varOUTid25, Total_current_out)
  ncOUTerror=nf90_put_var(ncOUTid,  varOUTid26, nu_Op_out)
  ncOUTerror=nf90_put_var(ncOUTid,  varOUTid27, nu_O2p_out)
  ncOUTerror=nf90_put_var(ncOUTid,  varOUTid28, nu_NOp_out)
  ncOUTerror=nf90_put_var(ncOUTid,  varOUTid29, nu_e_out)
  ncOUTerror=nf90_put_var(ncOUTid,  varOUTid30, Convenction_heating_out)
  ncOUTerror=nf90_put_var(ncOUTid,  varOUTid31, Wind_heating_out)
  ncOUTerror=nf90_put_var(ncOUTid,  varOUTid32, nu_in_out)

  ncOUTerror= nf90_close(ncOUTid)


  t2=MPI_Wtime()
  open(61,file=LOG,action='write',position='append')
  ! print*, "LOG::MPI Task ",myid ,"Done with file:",files,"out of",recvcnt,"in",t2-t1,"s"
  write(61,*) "LOG::MPI Task ",myid ,"Done with file:",files,"out of",recvcnt,"in",t2-t1,"s"
  close(61)
  ! print*, "LOG::MPI Task ",myid ,"Done with file:",files,"out of",recvcnt,"in",t2-t1,"s"

end do

call MPI_Barrier(mpi_comm_world,mpierror)
t_finish=MPI_Wtime()
if (myid==root) then
  print*, "DONE!"
  print*, "Housekeeping"
  open(61,file=LOG,action='write',position='append')
  write(61,*) "LOG:: All files calculated in " ,t_finish-t0 ,"seconds"
  close(61)
endif

!deallocate arrays
deallocate(glat)
deallocate(glon)
deallocate(glev)
deallocate(gtime)
deallocate(KP)
deallocate(zg)
deallocate(zg2)
deallocate(Ne_all)
deallocate(Ne_all2)
deallocate(PED_all)
deallocate(PED_all2)
deallocate(HALL_all)
deallocate(HALL_all2)
deallocate(Un_east)
deallocate(Un_east2)
deallocate(Un_north)
deallocate(Un_north2)
deallocate(Un_up)
deallocate(Un_up2)
deallocate(Ui_east)
deallocate(Ui_east2)
deallocate(Ui_north)
deallocate(Ui_north2)
deallocate(Ui_up)
deallocate(Ui_up2)
deallocate(Tn_all)
deallocate(Tn_all2)
deallocate(Te_all)
deallocate(Te_all2)
deallocate(Ti_all)
deallocate(Ti_all2)
deallocate(nO_all)
deallocate(nO_all2)
deallocate(nO2_all)
deallocate(nO2_all2)
deallocate(nN2_all)
deallocate(nN2_all2)
deallocate(NOp_all)
deallocate(NOp_all2)
deallocate(NO2p_all)
deallocate(NO2p_all2)
deallocate(NNOp_all)
deallocate(NNOp_all2)
deallocate(E_X)
deallocate(E_X2)
deallocate(E_Y)
deallocate(E_Y2)
deallocate(E_Z)
deallocate(E_Z2)
deallocate(BX_BMAG)
deallocate(BX_BMAG2)
deallocate(BY_BMAG)
deallocate(BY_BMAG2)
deallocate(BZ_BMAG)
deallocate(BZ_BMAG2)
deallocate(BMAG_G)
deallocate(BMAG2_G)

!Memory Allocation for other variables
deallocate(time_pf)
deallocate(zgf)
deallocate(zgout)
deallocate(mlat_qdf)
deallocate(mlt_qdf)
deallocate(BX_T)
deallocate(BY_T)
deallocate(BZ_T)
deallocate(BX_Tout)
deallocate(BY_Tout)
deallocate(BZ_Tout)
deallocate(BMAG_T)
deallocate(BMAG_Tout)
deallocate(E_X_si)
deallocate(E_X_si_out)
deallocate(E_Y_si)
deallocate(E_Y_si_out)
deallocate(E_Z_si)
deallocate(E_Z_si_out)
deallocate(Un_east_si)
deallocate(Un_east_si_out)
deallocate(Un_north_si)
deallocate(Un_north_si_out)
deallocate(Un_up_si)
deallocate(Un_up_si_out)
deallocate(Ui_east_si)
deallocate(Ui_east_si_out)
deallocate(Ui_north_si)
deallocate(Ui_north_si_out)
deallocate(Ui_up_si)
deallocate(Ui_up_si_out)
deallocate(nu_O2p_O2)
deallocate(nu_O2p_O2_out)
deallocate(nu_O2p_O)
deallocate(nu_O2p_O_out)
deallocate(nu_O2p_N2)
deallocate(nu_O2p_N2_out)
deallocate(nu_O2p)
deallocate(nu_O2p_out)
deallocate(nu_Op_O)
deallocate(nu_Op_O_out)
deallocate(nu_Op_O2)
deallocate(nu_Op_O2_out)
deallocate(nu_Op_N2)
deallocate(nu_Op_N2_out)
deallocate(nu_Op)
deallocate(nu_Op_out)
deallocate(nu_NOp_O2)
deallocate(nu_NOp_O2_out)
deallocate(nu_NOp_O)
deallocate(nu_NOp_O_out)
deallocate(nu_NOp_N2)
deallocate(nu_NOp_N2_out)
deallocate(nu_NOp)
deallocate(nu_NOp_out)
deallocate(nu_e)
deallocate(nu_e_out)
deallocate(nu_in)
deallocate(nu_in_out)
deallocate(Ohmic)
deallocate(Ohmicout)
deallocate(Joule)
deallocate(Jouleout)
deallocate(Frictional)
deallocate(Frictionalout)
deallocate(PED_calc)
deallocate(PED_calcout)
deallocate(HALL_calc)
deallocate(HALL_calcout)
deallocate(PAR_calc)
deallocate(PAR_calcout)
deallocate(CS_Op)
deallocate(CS_Opout)
deallocate(CS)
deallocate(CSout)
deallocate(CS_O2p)
deallocate(CS_O2pout)
deallocate(CS_NOp)
deallocate(CS_NOpout)
deallocate(Pedersen_current)
deallocate(Pedersen_current_out)
deallocate(Hall_current)
deallocate(Hall_current_out)
deallocate(Total_current)
deallocate(Total_current_out)
deallocate(Convenction_heating)
deallocate(Convenction_heating_out)
deallocate(Wind_heating)
deallocate(Wind_heating_out)



call MPI_Finalize(mpierror) !end mpi
end program




subroutine check_MPI_status(ierror)
  implicit none
  integer, intent(in):: ierror


  if (ierror .eq. 1) then
    stop 'COMM:MPI directive returned error status=1'
  endif

end subroutine

subroutine check_NetCDFvar_status(ierror,varname)
  implicit none
  integer, intent(in):: ierror
  character(len=*),intent(in):: varname
  if (ierror .ne. 0) then
    print*, 'COMM:NetCDF variable "', varname, '" returned error status=1'
    stop
  endif

end subroutine

subroutine geo_lat2geod_lat(phi_deg,geod_lat_deg)
  ! calculate geocentric latitude from geodetic latitude
  ! according to WGS 84
  implicit none
  integer, parameter :: dp = selected_real_kind(4)
  real(dp), intent(in)::phi_deg
  real(dp), intent(out):: geod_lat_deg
  real(dp):: geod_lat_rad
  real(dp)::a,f,b,e
  real(dp)::phi_rad
  real(dp),parameter::pi=3.14159
  a=6378137 !meter semi major axis of earth
  f=1/298.257 !flattening
  b=a-f*a !semi minor axis
  e=((a**2-b**2)**(1/2))/a
  phi_rad=phi_deg*(pi/180)
  geod_lat_rad=atan(tan(phi_rad)/(1-e**2))
  geod_lat_deg=geod_lat_rad*180/pi
end subroutine

subroutine enu_ecef(lat_phi,lon_lmd,Fenu,Fxyz)
  ! calculate geocentric latitude from geodetic latitude
  ! according to WGS 84
  implicit none
  integer, parameter :: dp = selected_real_kind(4)
  integer, parameter :: dp1 = selected_real_kind(8)
  real(dp1),dimension(3),intent(in)::Fenu
  real(dp), intent(in)::lat_phi, lon_lmd
  real(dp1),dimension(3),intent(out)::Fxyz
  real(dp),parameter::pi=3.14159
  real(dp)::fac,lat_ph,lon_lm
  real(dp1),dimension(3)::north_temp_unit,east_temp_unit,up_temp_unit
  real(dp1),dimension(3)::Fnorth_vector,Feast_vector,Fup_vector

  fac=pi/180
  lat_ph=lat_phi*fac
  lon_lm = lon_lmd * fac

  north_temp_unit = (/-cos(lon_lm) * sin(lat_ph), - sin(lon_lm) * sin(lat_ph), cos(lat_ph)/)
  east_temp_unit = (/-sin(lon_lm), cos(lon_lm), 0.0 /)
  up_temp_unit = (/cos(lon_lm) * cos(lat_ph), sin(lon_lm) * cos(lat_ph), sin(lat_ph)/)

  Fnorth_vector = (/Fenu(2) * north_temp_unit(1), Fenu(2) * north_temp_unit(2), Fenu(2) * north_temp_unit(3)/)
  Feast_vector = (/Fenu(1) * east_temp_unit(1), Fenu(1) * east_temp_unit(2), Fenu(1) * east_temp_unit(3)/)
  Fup_vector = (/Fenu(3) * up_temp_unit(1), Fenu(3) * up_temp_unit(2), Fenu(3) * up_temp_unit(3)/)

  Fxyz(1) = Fnorth_vector(1) + Feast_vector(1) + Fup_vector(1)
  Fxyz(2) = Fnorth_vector(2) + Feast_vector(2) + Fup_vector(2)
  Fxyz(3) = Fnorth_vector(3) + Feast_vector(3) + Fup_vector(3)

end subroutine

SUBROUTINE CROSS(A,B,D)
!**************************** DESCRIPTION*************************************
! CROSS PRODUCT ROUTINE
! INPUTS:
! A--> 3D MATRIX
! B--> 3D MATRIX
! OUTPUTS:
! C--> A X B CROSS PRODUCT 3D MATRIX
! ******************************************************************************
IMPLICIT NONE
INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(8)
REAL(dp),DIMENSION(3),INTENT(IN)::A,B
REAL(dp),DIMENSION(3),INTENT(OUT)::D
D(1)=A(2)*B(3)-A(3)*B(2)
D(2)=A(3)*B(1)-A(1)*B(3)
D(3)=A(1)*B(2)-A(2)*B(1)
END SUBROUTINE
