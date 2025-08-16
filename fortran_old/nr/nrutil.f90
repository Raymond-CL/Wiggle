MODULE nrutil
  USE nrtype
  IMPLICIT NONE
  INTEGER(I4B), PARAMETER :: NPAR_ARTH=16,NPAR2_ARTH=8
  INTERFACE reallocate
    MODULE PROCEDURE reallocate_iv,reallocate_im
  END INTERFACE
  INTERFACE arth
    MODULE PROCEDURE arth_i
  END INTERFACE
CONTAINS
!BL
  FUNCTION reallocate_iv(p,n)
  INTEGER(I4B), DIMENSION(:), POINTER :: p, reallocate_iv
  INTEGER(I4B), INTENT(IN) :: n
  INTEGER(I4B) :: nold,ierr
  allocate(reallocate_iv(n),stat=ierr)
  if (ierr /= 0) call &
    nrerror('reallocate_iv: problem in attempt to allocate memory')
  if (.not. associated(p)) RETURN
  nold=size(p)
  reallocate_iv(1:min(nold,n))=p(1:min(nold,n))
  deallocate(p)
  END FUNCTION reallocate_iv
!BL
  FUNCTION reallocate_im(p,n,m)
  INTEGER(I4B), DIMENSION(:,:), POINTER :: p, reallocate_im
  INTEGER(I4B), INTENT(IN) :: n,m
  INTEGER(I4B) :: nold,mold,ierr
  allocate(reallocate_im(n,m),stat=ierr)
  if (ierr /= 0) call &
    nrerror('reallocate_im: problem in attempt to allocate memory')
  if (.not. associated(p)) RETURN
  nold=size(p,1)
  mold=size(p,2)
  reallocate_im(1:min(nold,n),1:min(mold,m))= &
    p(1:min(nold,n),1:min(mold,m))
  deallocate(p)
  END FUNCTION reallocate_im
!BL
  SUBROUTINE nrerror(string)
  CHARACTER(LEN=*), INTENT(IN) :: string
  write (*,*) 'nrerror: ',string
  STOP 'program terminated by nrerror'
  END SUBROUTINE nrerror
!BL
  FUNCTION arth_i(first,increment,n)
  INTEGER(I4B), INTENT(IN) :: first,increment,n
  INTEGER(I4B), DIMENSION(n) :: arth_i
  INTEGER(I4B) :: k,k2,temp
  if (n > 0) arth_i(1)=first
  if (n <= NPAR_ARTH) then
    do k=2,n
      arth_i(k)=arth_i(k-1)+increment
    end do
  else
    do k=2,NPAR2_ARTH
      arth_i(k)=arth_i(k-1)+increment
    end do
    temp=increment*NPAR2_ARTH
    k=NPAR2_ARTH
    do
      if (k >= n) exit
      k2=k+k
      arth_i(k+1:min(k2,n))=temp+arth_i(1:min(k,n-k))
      temp=temp+temp
      k=k2
    end do
  end if
  END FUNCTION arth_i
END MODULE nrutil
