MODULE nr
  INTERFACE ran1
    SUBROUTINE ran1_s(harvest)
      USE nrtype
      REAL(WP), INTENT(OUT) :: harvest
    END SUBROUTINE ran1_s
  END INTERFACE
  INTERFACE
    SUBROUTINE vegas(region,func,init,ncall,itmx,nprn,tgral,sd,chi2a)
      USE nrtype
      REAL(WP), DIMENSION(:), INTENT(IN) :: region
      INTEGER(I4B), INTENT(IN) :: init,ncall,itmx,nprn
      REAL(WP), INTENT(OUT) :: tgral,sd,chi2a
      INTERFACE
        FUNCTION func(pt,wgt)
          USE nrtype
          REAL(WP), DIMENSION(:), INTENT(IN) :: pt
          REAL(WP), INTENT(IN) :: wgt
          REAL(WP) :: func
        END FUNCTION func
      END INTERFACE
    END SUBROUTINE vegas
  END INTERFACE
END MODULE nr
