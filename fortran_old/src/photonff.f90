module photonff
  
  use nrtype
  use comvar
  implicit none

contains

  ! photon distribution
  function xf(k,kp) result(res)
  implicit none
  real(wp) :: res
  real(wp), intent(in) :: k,kp
  real(wp) :: nff,fac
  res = 0d0
  if(k.le.0d0 .or. kp.le.0d0) return
  fac = atomZ*atomZ * alphae /PI/PI
  nff = ffnuc(k)/k/k * ffnuc(kp)/kp/kp
  res = fac * nff
  return
  end function xf

  ! nucleus form factor
  function ffnuc(k) result(res)
  implicit none
  real(wp) :: res
  real(wp), intent(in) :: k
  res = sin(k*RA) - k*RA * cos(k*RA)
  res = res * 3d0 / (k*RA)**3 / (1d0+a0**2*k**2)
  return
  end function ffnuc

end module photonff
