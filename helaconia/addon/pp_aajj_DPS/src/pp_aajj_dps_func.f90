MODULE pp_aajj_dps_func
  USE pp_aajj_dps_global
  IMPLICIT NONE
CONTAINS
  SUBROUTINE setscale_dps(ipip,scale)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::ipip
    REAL(KIND(1d0)),INTENT(OUT)::scale
    scale=dps_pmom(ipip,3,1)**2+dps_pmom(ipip,3,2)**2
    scale=DSQRT(scale)
    RETURN
  END SUBROUTINE setscale_dps

  SUBROUTINE Combine_scale_dps(scale1,scale2,scale)
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::scale1,scale2
    REAL(KIND(1d0)),INTENT(OUT)::scale
    scale=DSQRT(scale1*scale2)
    RETURN
  END SUBROUTINE Combine_scale_dps
END MODULE pp_aajj_dps_func
