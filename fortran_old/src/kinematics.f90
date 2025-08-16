module kinematics

  use nrtype
  implicit none
 
  ! center of mass energy
  real(wp) :: CME,CME2

  ! observable
  real(wp) :: almin,almax,alpha
  integer :: nalpha

  ! transverse momentum imbalances
  real(wp) :: pt_imb,ptmin,ptmax,phi_pt
  real(wp) :: pl_avg,plmin,plmax,phi_pl
  real(wp) :: qt_imb,qtmin,qtmax,phi_qt
  integer :: qtn

  ! lepton transverse momentum, azi_angle, rapidity, mass
  real(wp) :: p1t,p2t
  real(wp) :: phip1,phip2
  real(wp) :: y1,y2
  real(wp) :: pltmin,pltmax
  real(wp) :: ylmin,ylmax

  ! lepton transverse mass
  real(wp) :: mtsq,mt
  real(wp) :: m1tsq,m2tsq
  real(wp) :: m1t,m2t

  ! photon transverse momentum, angle
  real(wp) :: ktmin,ktmax
  real(wp) :: k1t,kat
  real(wp) :: phik1,phika
  real(wp) :: k2t,kbt
  real(wp) :: phik2,phikb
  real(wp) :: k1sq,kasq,k2sq,kbsq
  real(wp) :: k1ka,k2kb
  real(wp) :: dotterm,crossterm
  real(wp) :: d1a2b,d12ab,d1ba2
  real(wp) :: d1pap,d2pbp
  real(wp) :: Oisoterm, Omegaterm
  real(wp) :: phi_dt
  !real(wp) :: k1x,k1y,k2x,k2y
  !real(wp) :: kax,kay,kbx,kby

  ! lepton masses
  integer(i1b) :: Mli
  real(wp) :: m_lep
  real(wp) :: Mllmin,Mll,Mllmax

  ! dilepton rapidity gap
  real(wp) :: yllmin,yll,yllmax

  ! dilepton asymmetry cut
  real(wp) :: aj,ajcut

  ! momentum fraction
  real(wp) :: x1,x2

  ! Mandelstam
  real(wp) :: mans,mant,manu

  ! impact parameter
  real(wp) :: bpmin,bp,bpmax
  integer :: bpn

contains

  subroutine setkinematics
  use nrtype
  use comvar
  use photonff
  use vector2d
  implicit none
  real(wp) :: ml2,mp2,deltaphi
  type(vec2d) :: p1v,p2v,ptv,qtv
  type(vec2d) :: k1v,k2v,kav,kbv
  type(vec2d) :: k1av
  real(wp) :: pitermij,pitermkl
  real(wp) :: ptmax

  !deltaphi = (1d0-alpha)*PI
  !phip1 = 0d0
  !phip2 = phip1 + deltaphi
  !call p1v%setpol(p1t,phip1)
  !call p2v%setpol(p2t,phip2)
  call ptv%setpol(pl_avg,phi_pl)
  call qtv%setpol(qt_imb,phi_qt)
  p1v = 0.5d0*qtv + ptv
  p2v = 0.5d0*qtv - ptv

  ml2 = m_lep*m_lep
  mp2 = M_pro*M_pro
  !m1t = sqrt( ml2 + p1v%m2() )
  !m2t = sqrt( ml2 + p2v%m2() )
  p1t = p1v%m()
  p2t = p2v%m()

  !x1 = (m1t*exp(+y1) + m2t*exp(+y2))/CME
  !x2 = (m1t*exp(-y1) + m2t*exp(-y2))/CME
  ptmax = pl_avg !max(p1v%m(),p2v%m())
  x1 = ptmax * (exp(+y1) + exp(+y2))/CME
  x2 = ptmax * (exp(-y1) + exp(-y2))/CME

  mans = +x1*x2*CME2
  !mant = ml2 - x1*CME*m1t*exp(-y1)
  !manu = ml2 - x1*CME*m2t*exp(-y2)
  mant = - x1*CME*ptmax*exp(-y1)
  manu = - x1*CME*ptmax*exp(-y2)

  !yll = abs(y1-y2)
  Mll = sqrt(mans)
  !aj = abs(p1t-p2t)/(p1t+p2t)

  !ptv=p1v+p2v
  !qtv=ptv   

  call k1v%setpol(k1t,phik1)
  call kav%setpol(kat,phika)
  k2v=qtv-k1v
  kbv=qtv-kav
  k2t=k2v%m()
  kbt=kbv%m()
  phik2 = k2v%f()
  phikb = kbv%f()

  k1sq = x1*x1*mp2 + k1v%m2()
  kasq = x1*x1*mp2 + kav%m2()
  k2sq = x2*x2*mp2 + k2v%m2()
  kbsq = x2*x2*mp2 + kbv%m2()

  k1av=k1v-kav
  k1ka=abs(k1av)
  phi_dt = k1av%f()

  !d1a2b = (k1v*kav)*(k2v*kbv)
  !d12ab = (k1v*k2v)*(kav*kbv)
  !d1ba2 = (k1v*kbv)*(kav*k2v)

  !d1pap = (k1v*ptv)*(kav*ptv)
  !d2pbp = (k2v*ptv)*(kbv*ptv)

  !pitermij = 2d0*d1pap/(ptv*ptv)-(k1v*kav)
  !pitermkl = 2d0*d2pbp/(ptv*ptv)-(k2v*kbv)

  !Oisoterm = 2d0*(mant/manu + manu/mant)*(d1a2b-d12ab+d1ba2)
  !Omegaterm = 2d0*pitermij*pitermkl-(d1ba2+d12ab-d1a2b)

  !sigma0 = alphae**2/mans/mans * (Oisoterm - 4d0*Omegaterm)
  sigma0 = 2d0*alphae**2/mans/mans * k1t*kat*k2t*kbt &
          * ( (mant/manu+manu/mant)*cos(phik1-phika+phik2-phikb) &
            - 2d0*cos(phik1+phika+phik2+phikb-4d0*phi_pl) )

  x1f1 = xf(sqrt(k1sq) , sqrt(kasq))
  x2f2 = xf(sqrt(k2sq) , sqrt(kbsq))

  end subroutine setkinematics

end module kinematics

