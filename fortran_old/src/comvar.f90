module comvar

  use nrtype
  implicit none
  public
 
  ! conversion units
  real(wp), parameter :: planck = 6.62607015d-34    ! Planck constant (J.s)
  real(wp), parameter :: planckbar = planck/twoPI   ! reduced Planck (J.s)
  real(wp), parameter :: charge = 1.602176634d-19   ! chrg mag (J/eV)
  real(wp), parameter :: hbar = planckbar/charge    ! hbar (eV.s)
  real(wp), parameter :: light = 299792458d0        ! light speed (m/s)
  real(wp), parameter :: hbarc = hbar*light         ! conversion unit (eV.m)
  real(wp), parameter :: giga = 1d+9
  real(wp), parameter :: mili = 1d-3
  real(wp), parameter :: micro = 1d-6
  real(wp), parameter :: femto = 1d-15
  real(wp), parameter :: fm2barn = 1d-2             ! fm^2->barns (barn/fm^2)
  real(wp), parameter :: gevfm = hbarc/giga/femto   ! common conversion (GeV.fm)

  ! constant mass (in GeV)
  real(wp), parameter :: M_pro = 0.93827208816d0
  real(wp), parameter :: M_neu = 0.93956542052d0
  real(wp), parameter :: M_ele = 0.00051d0
  real(wp), parameter :: M_mu  = 0.10566d0
  real(wp), parameter :: M_tau = 1.77686d0

  ! io-file name
  character(len=*), parameter :: ifile='input.dat'
  character(len=*), parameter :: ofile='output.dat'

  ! loop variables
  integer(i4b) :: nbin,ind
  real(wp) :: hmin,hmax,bin
  real(wp) :: binL,binM,binR

  ! temporary test variables, counters
  real    :: tempr1,tempr2
  integer :: tempi1,tempi2

  ! event counters
  integer(8) :: totevnt,accevnt
  real(wp) :: eff

  ! flags
  logical :: debug
  logical :: startcount

  ! nucleus info
  integer(i2b) :: atomA,atomZ
  real(wp) :: RA
  real(wp), parameter :: a0 = 1d0/sqrt(2d0) /gevfm !~0.71 fm
  real(wp) :: bupc

  ! coupling 
  real(wp), parameter :: alphae = 1d0/137.036d0
  real(wp) :: alphas

  ! overall
  real(wp) :: sigma0,x1f1,x2f2,bessel

end module comvar
