program main
  use, intrinsic :: iso_fortran_env, only:stdout => output_unit
  use nrtype
  use nr, only : vegas
  use comvar
  use kinematics
  implicit none
  real :: tstart,tstop
  integer :: u,stat

  ! vegas variables 
  integer(i4b) :: ncall1,itmax1
  integer(i4b) :: ncall2,itmax2
  integer(i4b) :: init,ndim,nprn
  real(wp) :: avgi,chi2a,sd
  real(wp), dimension(30) :: region

  ! declare external interface
  interface
    function func(dx,wgt)
    use nrtype
    implicit none
    real(wp), dimension(:), intent(in) :: dx
    real(wp), intent(in) :: wgt
    real(wp) :: func
    end function func
    function getu() result(u)
    integer :: u
    end function getu
    subroutine print_prog_info(u)
    integer :: u
    end subroutine print_prog_info
  end interface

  call cpu_time(tstart)

  ! some initial settings
  CME = 5020d0
  !bpmin = 0d0;  bpmax = 20d0; bpn = 40
  !bpmin = bpmin / gevfm
  !bpmax = bpmax / gevfm
  qtmin = 0.0d0;  qtmax = 0.1d0;  qtn = 50
  ylmin = -1d0; ylmax = +1d0
  pltmin = 4d0; pltmax = 20d0
  Mllmin = 4d0; Mllmax = 45d0
  ktmin = 0d0;  ktmax = 10d0
  atomA = 207;  atomZ = 82
  CME2 = CME*CME
  m_lep = M_mu
  !RA = 6.62d0 / gevfm
  RA = 6.62d0 / gevfm
  ncall1 = 100000;   itmax1 = 10
  ncall2 = 1000000;  itmax2 = 1

  ! initialize vegas variables
  nprn = -1
  avgi = 0d0; sd = 0d0; chi2a = 0d0
  ndim = 9;
  region(1) = ktmin;      region(ndim+1) = ktmax
  region(2) = 0d0;        region(ndim+2) = twoPI
  region(3) = ktmin;      region(ndim+3) = ktmax
  region(4) = 0d0;        region(ndim+4) = twoPI
  region(5) = ylmin;      region(ndim+5) = ylmax
  region(6) = ylmin;      region(ndim+6) = ylmax
  region(7) = pltmin;     region(ndim+7) = pltmax
  region(8) = 0d0;        region(ndim+8) = twoPI
  region(9) = 0d0;        region(ndim+9) = twoPI

  ! initialize output files
  u = getu()
  open(newunit=u,file=ofile,status='replace')
  close(u)

  ! loop definition
  nbin = qtn  !bpn
  hmin = qtmin  !bpmin
  hmax = qtmax  !bpmax
  bin = merge( (hmax-hmin)/nbin , log(hmax/hmin)/nbin , .true. )
  do ind = 1,nbin
    binL = merge( hmin+bin*(ind-1) , hmin*exp(bin*(ind-1)) , .true. )
    binR = merge( hmin+bin*ind     , hmin*exp(bin*ind)     , .true. )
    binM = (binL+binR)/2d0
    qt_imb = binM

    !accevnt=0; totevnt=0;
    init = -1
    !startcount = .false.
    call vegas(region(1:2*ndim),func,init,ncall1,itmax1,nprn,avgi,sd,chi2a)
    init = +1
    !startcount = .true.
    call vegas(region(1:2*ndim),func,init,ncall2,itmax2,nprn,avgi,sd,chi2a)

    !eff = real(accevnt,wp)/real(totevnt,wp)*100d0
    u = stdout
    write(u,*) qt_imb,avgi,sd!,real(eff),'%'
    u = getu()
    open(u,file=ofile,position='append')
    write(u,*) qt_imb,avgi,sd/avgi!,real(eff),'%'
    close(u)
  enddo

  call cpu_time(tstop)
  write(stdout,*) 'time elapsed: ',tstop-tstart,'seconds.'
end program main



function func(dx,wgt)
  use nrtype
  use comvar
  use kinematics
  implicit none
  real(wp), dimension(:), intent(in) :: dx
  real(wp), intent(in) :: wgt
  real(wp) :: func
  real(wp) :: costerm
  !interface
  !  subroutine print_debug_info()
  !  end subroutine print_debug_info
  !end interface

  if(startcount) totevnt=totevnt+1
  func = 0d0

  k1t     = dx(1)
  phik1   = dx(2)
  kat     = dx(3)
  phika   = dx(4)
  y1      = dx(5)
  y2      = dx(6)
  pl_avg  = dx(7)
  phi_pl  = dx(8)
  phi_qt  = dx(9)

  ! k1t = 2d0
  ! phik1 = 1d0
  ! kat = 2d0
  ! phika = 1d0
  ! y1 = 0d0
  ! y2 = 0d0
  ! pl_avg = 10d0
  ! phi_pl = 1d0
  ! phi_qt = 2d0
  ! qt_imb = 0.05d0

  bp = 10d0 / gevfm

  call setkinematics 

  if(x1.le.0d0 .or. x1.ge.1d0) return
  if(x2.le.0d0 .or. x2.ge.1d0) return
  if(Mll.le.Mllmin .or. Mll.ge.Mllmax) return
  !if(aj.ge.ajcut) return

  

  bessel = bp / twoPI * bessel_jn(0,bp*k1ka)
  !bessel = bp / twoPI * bessel_jn(2,bp*k1ka)
  !bessel = bp / twoPI * bessel_jn(4,bp*k1ka)
  !bessel = bp / twoPI * bessel_jn(6,bp*k1ka)
  !bessel = bp / twoPI * bessel_jn(8,bp*k1ka)

  costerm = 1d0
  !costerm = 2d0*cos(0d0*phi_pl + 0d0*phi_qt + 0d0*phi_dt)    ! 000
  !costerm = -2d0*cos(4d0*phi_pl - 4d0*phi_qt)

  func = pl_avg * qt_imb * k1t * kat * bessel
  func = func * x1f1 * x2f2 * sigma0
  !func = func * gevfm**2 * fm2barn / micro    ! GeV^-2 -> mb
  func = func * costerm

  ! write(*,*) x1f1, x2f2

  !if(isnan(func) .or. debug) then
    !write(*,*) 'something wrong'
    !call print_debug_info()
  !endif

  if(startcount) accevnt=accevnt+1
  return
end function func




function getu() result(u)
  implicit none
  integer :: u
  logical :: ex
  do u=10,100
    inquire(unit=u,opened=ex)
    if(.not. ex) return
  enddo
end function getu




subroutine print_prog_info(u)
  use, intrinsic :: iso_fortran_env, only:compiler_version
  use kinematics
  implicit none
  integer :: u
  write(u,*) '***********************************'
  write(u,*) '     dilepton photo-production     '
  write(u,*) 'Phys. Rev. Lett. 122, 132301 (2019)'
  write(u,*) 'Phys. Rev. D     102, 094013 (2020)'
  write(u,*) '***********************************'
  write(u,*) 'program compiled by:',compiler_version()
  !write(u,*) 'with compiler options:',compiler_options()
  write(u,*) '***********************************'
end subroutine print_prog_info
