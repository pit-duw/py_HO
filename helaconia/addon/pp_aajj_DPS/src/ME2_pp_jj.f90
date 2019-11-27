MODULE ME2_pp_jj
  USE pp_aajj_dps_global
  USE Helac_Global
  IMPLICIT NONE
CONTAINS
  SUBROUTINE ME2_gg_gg(wme)
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(OUT)::wme
    REAL(KIND(1d0))::s,t
    s=s12_jj
    t=t13_jj
    ! it has been averaged over initial color and helicity
    ! and final state symmetry factor
    wme=(9d0*Gstrong**4*(s**2 + s*t + t**2)**3)/(4.*s**2*t**2*(s + t)**2)
    RETURN
  END SUBROUTINE ME2_gg_gg

  SUBROUTINE CF_gg_gg(icol,rcol)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::icol
    REAL(KIND(1d0)),INTENT(IN)::rcol
    REAL(KIND(1d0))::s,t,amp2col,amp2tot
    INTEGER::init=0,offset,i,j_col
    INTEGER,DIMENSION(6,4,2)::icol_gg_save
    REAL(KIND(1d0)),DIMENSION(6)::amp2
    SAVE init,icol_gg_save
    IF(init.EQ.0)THEN
       icol_gg_save(1:6,1,1)=1
       icol_gg_save(1:6,2,1)=2
       icol_gg_save(1:6,3,1)=3
       icol_gg_save(1:6,4,1)=4
       ! Tr(1,2,3,4)
       icol_gg_save(1,1,2)=4
       icol_gg_save(1,2,2)=1
       icol_gg_save(1,3,2)=2
       icol_gg_save(1,4,2)=3
       ! Tr(1,2,4,3)
       icol_gg_save(2,1,2)=3
       icol_gg_save(2,2,2)=1
       icol_gg_save(2,3,2)=4
       icol_gg_save(2,4,2)=2
       ! Tr(1,3,2,4)
       icol_gg_save(3,1,2)=4
       icol_gg_save(3,2,2)=3
       icol_gg_save(3,3,2)=1
       icol_gg_save(3,4,2)=2
       ! Tr(1,3,4,2)
       icol_gg_save(4,1,2)=2
       icol_gg_save(4,2,2)=4
       icol_gg_save(4,3,2)=1
       icol_gg_save(4,4,2)=3
       ! Tr(1,4,2,3)
       icol_gg_save(5,1,2)=3
       icol_gg_save(5,2,2)=4
       icol_gg_save(5,3,2)=2
       icol_gg_save(5,4,2)=1
       ! Tr(1,4,3,2)
       icol_gg_save(6,1,2)=2
       icol_gg_save(6,2,2)=3
       icol_gg_save(6,3,2)=4
       icol_gg_save(6,4,2)=1
       init=1
    ENDIF
    !IF(icol.EQ.1)THEN
    !   offset=0
    !ELSE
    !   offset=4
    !ENDIF
    s=s12_jj
    t=t13_jj
    ! no initial average and final symmetry
    ! no Gstrong**4
    ! Tr(1,2,3,4)
    amp2(1)=(64d0*(s**2 + s*t + t**2)**2)/(s**2*(s + t)**2)
    ! Tr(1,2,4,3)
    amp2(2)=(64d0*(s**2 + s*t + t**2)**2)/(s**2*t**2)
    ! Tr(1,3,2,4)
    amp2(3)=(64d0*(s**2 + s*t + t**2)**2)/(t**2*(s + t)**2)
    ! Tr(1,3,4,2)
    amp2(4)=(64d0*(s**2 + s*t + t**2)**2)/(s**2*t**2)
    ! Tr(1,4,2,3)
    amp2(5)=(64d0*(s**2 + s*t + t**2)**2)/(t**2*(s + t)**2)
    ! Tr(1,4,3,2)
    amp2(6)=(64d0*(s**2 + s*t + t**2)**2)/(s**2*(s + t)**2)
    amp2tot=0d0
    DO i=1,6
       amp2tot=amp2tot+amp2(i)
    ENDDO
    amp2col=0
    j_col=0
    IF(amp2tot.LE.0d0)THEN
       ! randomly choosing
       DO i=1,6
          amp2col=amp2col+1d0/6d0
          IF(rcol.LE.amp2col)THEN
             j_col=i
             EXIT
          ENDIF
       ENDDO
    ELSE
       ! choosing following amp2
       DO i=1,6
          amp2col=amp2col+amp2(i)/amp2tot
          IF(rcol.LE.amp2col)THEN
             j_col=i
             EXIT
          ENDIF
       ENDDO
    ENDIF
    DO i=1,4
       icolun_dps(icol,i,1)=icol_gg_save(j_col,i,2)
       icolun_dps(icol,i,2)=icol_gg_save(j_col,i,1)
    ENDDO
    RETURN
  END SUBROUTINE CF_gg_gg

  SUBROUTINE ME2_gg_qq(wme)
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(OUT)::wme
    REAL(KIND(1d0))::s,t
    s=s12_jj
    t=t13_jj
    wme=-(Gstrong**4*(s**2 + 2*s*t + 2*t**2)*&
         (4*s**2 + 9*s*t + 9*t**2))/(24.*s**2*t*(s + t))
    RETURN
  END SUBROUTINE ME2_gg_qq

  SUBROUTINE CF_gg_qq(icol,rcol)
    ! make sure it is g g > q q~
    IMPLICIT NONE
    INTEGER,INTENT(IN)::icol
    REAL(KIND(1d0)),INTENT(IN)::rcol
    INTEGER::init=0,j_col,offset,i,j
    SAVE init
    INTEGER,DIMENSION(2,4,2)::icol_gg_save
    SAVE icol_gg_save
    REAL(KIND(1d0)),DIMENSION(2)::amp2
    REAL(KIND(1d0))::amp2tot,amp2col,s,t
    IF(init.EQ.0)THEN
       ! T(1,2)
       icol_gg_save(1,1,1)=1
       icol_gg_save(1,2,1)=2
       icol_gg_save(1,3,1)=3
       icol_gg_save(1,4,1)=0
       icol_gg_save(1,1,2)=3
       icol_gg_save(1,2,2)=1
       icol_gg_save(1,3,2)=0
       icol_gg_save(1,4,2)=2
       ! T(2,1)
       icol_gg_save(2,1,1)=1
       icol_gg_save(2,2,1)=2
       icol_gg_save(2,3,1)=3
       icol_gg_save(2,4,1)=0
       icol_gg_save(2,1,2)=2
       icol_gg_save(2,2,2)=3
       icol_gg_save(2,3,2)=0
       icol_gg_save(2,4,2)=1
       init=1
    ENDIF
    !IF(icol.EQ.1)THEN
    !   offset=0
    !ELSE
    !   offset=4
    !ENDIF
    s=s12_jj
    t=t13_jj
    ! no initial avarage and final symmetry
    ! no Gstrong**4
    ! T(1,2)
    amp2(1)=(-8d0*(s + t)*(s**2 + 2*s*t + 2*t**2))/(s**2*t)
    ! T(2,1)
    amp2(2)=(-8d0*t*(s**2 + 2*s*t + 2*t**2))/(s**2*(s + t))
    amp2tot=0d0
    DO i=1,2
       amp2tot=amp2tot+amp2(i)
    ENDDO
    amp2col=0
    j_col=0
    IF(amp2tot.LE.0d0)THEN
       ! randomly choosing
       DO i=1,2
          amp2col=amp2col+1d0/6d0
          IF(rcol.LE.amp2col)THEN
             j_col=i
             EXIT
          ENDIF
       ENDDO
    ELSE
       ! choosing following amp2
       DO i=1,2
          amp2col=amp2col+amp2(i)/amp2tot
          IF(rcol.LE.amp2col)THEN
             j_col=i
             EXIT
          ENDIF
       ENDDO
    ENDIF
    DO i=1,4
       icolun_dps(icol,i,1)=icol_gg_save(j_col,i,2)
       icolun_dps(icol,i,2)=icol_gg_save(j_col,i,1)
    ENDDO
    RETURN
  END SUBROUTINE CF_gg_qq

  SUBROUTINE ME2_qq_gg(wme)
    ! q q~ > g g
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(OUT)::wme
    REAL(KIND(1d0))::s,t
    s=s12_jj
    t=t13_jj
    wme=(-4d0*Gstrong**4*(s**2 + 2*s*t + 2*t**2)*&
         (4*s**2 + 9*s*t + 9*t**2))/(27.*s**2*t*(s + t))
    RETURN
  END SUBROUTINE ME2_qq_gg

  SUBROUTINE CF_qq_gg(iqqx,icol,rcol)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::iqqx ! iqqx: 1 q qx, 2 qx q
    INTEGER,INTENT(IN)::icol
    REAL(KIND(1d0)),INTENT(IN)::rcol
    REAL(KIND(1d0))::s,t,amp2tot,amp2col
    INTEGER::init=0,offset,i,j_col
    SAVE init
    INTEGER,DIMENSION(2,2,4,2)::icol_qq_save
    SAVE icol_qq_save
    REAL(KIND(1d0)),DIMENSION(2)::amp2
    IF(init.EQ.0)THEN
       ! iqqx = 1
       ! q qx
       ! T(3,4)
       icol_qq_save(1,1,1,1)=0
       icol_qq_save(1,1,2,1)=2
       icol_qq_save(1,1,3,1)=3
       icol_qq_save(1,1,4,1)=1
       icol_qq_save(1,1,1,2)=1
       icol_qq_save(1,1,2,2)=0
       icol_qq_save(1,1,3,2)=2
       icol_qq_save(1,1,4,2)=3
       ! T(4,3)
       icol_qq_save(1,2,1,1)=0
       icol_qq_save(1,2,2,1)=2
       icol_qq_save(1,2,3,1)=3
       icol_qq_save(1,2,4,1)=1
       icol_qq_save(1,2,1,2)=3
       icol_qq_save(1,2,2,2)=0
       icol_qq_save(1,2,3,2)=1
       icol_qq_save(1,2,4,2)=2
       ! iqqx = 2
       ! qx q
       ! T(3,4)
       icol_qq_save(2,1,1,1:2)=icol_qq_save(1,1,2,1:2)
       icol_qq_save(2,1,2,1:2)=icol_qq_save(1,1,1,1:2)
       icol_qq_save(2,1,3:4,1:2)=icol_qq_save(1,1,3:4,1:2)
       ! T(4,3)
       icol_qq_save(2,2,1,1:2)=icol_qq_save(1,2,2,1:2)
       icol_qq_save(2,2,2,1:2)=icol_qq_save(1,2,1,1:2)
       icol_qq_save(2,2,3:4,1:2)=icol_qq_save(1,2,3:4,1:2)
       init=1
    ENDIF
    !IF(icol.EQ.1)THEN
    !   offset=0
    !ELSE
    !   offset=4
    !ENDIF
    s=s12_jj
    IF(iqqx.EQ.1)THEN
       t=t13_jj
    ELSE
       t=-t13_jj-s12_jj
    ENDIF
    ! no initial avarage and final symmetry
    ! no Gstrong**4
    ! T(3,4)
    amp2(1)=(-8d0*t*(s**2 + 2*s*t + 2*t**2))/(s**2*(s + t))
    ! T(4,3)
    amp2(2)=(-8d0*(s + t)*(s**2 + 2*s*t + 2*t**2))/(s**2*t)
        amp2tot=0d0
    DO i=1,2
       amp2tot=amp2tot+amp2(i)
    ENDDO
    amp2col=0
    j_col=0
    IF(amp2tot.LE.0d0)THEN
       ! randomly choosing
       DO i=1,2
          amp2col=amp2col+1d0/6d0
          IF(rcol.LE.amp2col)THEN
             j_col=i
             EXIT
          ENDIF
       ENDDO
    ELSE
       ! choosing following amp2
       DO i=1,2
          amp2col=amp2col+amp2(i)/amp2tot
          IF(rcol.LE.amp2col)THEN
             j_col=i
             EXIT
          ENDIF
       ENDDO
    ENDIF
    DO i=1,4
       icolun_dps(icol,i,1)=icol_qq_save(iqqx,j_col,i,2)
       icolun_dps(icol,i,2)=icol_qq_save(iqqx,j_col,i,1)
    ENDDO
    RETURN
  END SUBROUTINE CF_qq_gg

  SUBROUTINE ME2_gq_gq(iqqx,wme)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::iqqx ! iqqx: 1, q+g, 2, qx+g, 3, g+q, 4, g+qx
    REAL(KIND(1d0)),INTENT(OUT)::wme
    REAL(KIND(1d0))::s,t
    s=s12_jj
    IF(iqqx.EQ.3.OR.iqqx.EQ.4)THEN
       ! g+q or g+qx
       t=t13_jj
    ELSEIF(iqqx.EQ.1.OR.iqqx.EQ.2)THEN
       ! q+g or qx+g
       t=-t13_jj-s12_jj
    ELSE
       WRITE(*,*)"ERROR:Unknown iqqx=",iqqx
       STOP
    ENDIF
    wme=(Gstrong**4*(2*s**2 + 2*s*t + t**2)*(9*s**2 &
         + 9*s*t + 4*t**2))/(9.*s*t**2*(s + t))
    RETURN
  END SUBROUTINE ME2_gq_gq

  SUBROUTINE CF_gq_gq(iqqx,icol,rcol)
    ! g q > g q
    IMPLICIT NONE
    INTEGER,INTENT(IN)::iqqx ! iqqx: 1, q+g, 2, qx+g, 3, g+q, 4, g+qx
    INTEGER,INTENT(IN)::icol
    REAL(KIND(1d0)),INTENT(IN)::rcol
    INTEGER::init=0,offset,i,j_col
    SAVE init
    INTEGER,DIMENSION(4,2,4,2)::icol_gq_save
    SAVE icol_gq_save
    REAL(KIND(1d0)),DIMENSION(2)::amp2
    REAL(KIND(1d0))::s,t,amp2tot,amp2col
    IF(init.EQ.0)THEN
       ! iqqx = 3
       ! g+q
       ! T(1,3)
       icol_gq_save(3,1,1,1)=1
       icol_gq_save(3,1,2,1)=0
       icol_gq_save(3,1,3,1)=2
       icol_gq_save(3,1,4,1)=3
       icol_gq_save(3,1,1,2)=3
       icol_gq_save(3,1,2,2)=2
       icol_gq_save(3,1,3,2)=1
       icol_gq_save(3,1,4,2)=0
       ! T(3,1)
       icol_gq_save(3,2,1,1)=1
       icol_gq_save(3,2,2,1)=0
       icol_gq_save(3,2,3,1)=2
       icol_gq_save(3,2,4,1)=3
       icol_gq_save(3,2,1,2)=2
       icol_gq_save(3,2,2,2)=1
       icol_gq_save(3,2,3,2)=3
       icol_gq_save(3,2,4,2)=0
       ! iqqx = 1
       ! q+g
       ! T(2,3)
       icol_gq_save(1,1,1,1:2)=icol_gq_save(3,1,2,1:2)
       icol_gq_save(1,1,2,1:2)=icol_gq_save(3,1,1,1:2)
       icol_gq_save(1,1,3:4,1:2)=icol_gq_save(3,1,3:4,1:2)
       ! T(3,2)
       icol_gq_save(1,2,1,1:2)=icol_gq_save(3,2,2,1:2)
       icol_gq_save(1,2,2,1:2)=icol_gq_save(3,2,1,1:2)
       icol_gq_save(1,2,3:4,1:2)=icol_gq_save(3,2,3:4,1:2)
       ! iqqx = 4
       ! g+qx
       ! T(1,3)
       icol_gq_save(4,1,1,1)=1
       icol_gq_save(4,1,2,1)=2
       icol_gq_save(4,1,3,1)=3
       icol_gq_save(4,1,4,1)=0
       icol_gq_save(4,1,1,2)=2
       icol_gq_save(4,1,2,2)=0
       icol_gq_save(4,1,3,2)=1
       icol_gq_save(4,1,4,2)=3
       ! T(3,1)
       icol_gq_save(4,1,1,1)=1 
       icol_gq_save(4,1,2,1)=2
       icol_gq_save(4,1,3,1)=3
       icol_gq_save(4,1,4,1)=0
       icol_gq_save(4,1,1,2)=3
       icol_gq_save(4,1,2,2)=0
       icol_gq_save(4,1,3,2)=2
       icol_gq_save(4,1,4,2)=1
       ! iqqx = 2
       ! qx+g
       ! T(2,3)
       icol_gq_save(2,1,1,1:2)=icol_gq_save(4,1,2,1:2)
       icol_gq_save(2,1,2,1:2)=icol_gq_save(4,1,1,1:2)
       icol_gq_save(2,1,3:4,1:2)=icol_gq_save(4,1,3:4,1:2)
       ! T(3,2)
       icol_gq_save(2,2,1,1:2)=icol_gq_save(4,2,2,1:2)
       icol_gq_save(2,2,2,1:2)=icol_gq_save(4,2,1,1:2)
       icol_gq_save(2,2,3:4,1:2)=icol_gq_save(4,2,3:4,1:2)
       init=1
    ENDIF
    !IF(icol.EQ.1)THEN
    !   offset=0
    !ELSE
    !   offset=4
    !ENDIF
    s=s12_jj
    IF(iqqx.EQ.3.OR.iqqx.EQ.4)THEN
       t=t13_jj
    ELSEIF(iqqx.EQ.1.OR.iqqx.EQ.2)THEN
       t=-t13_jj-s12_jj
    ENDIF
    ! no initial avarage and final symmetry
    ! no Gstrong**4
    ! T(1,3)
    amp2(1)=(8d0*s*(2*s**2 + 2*s*t + t**2))/(t**2*(s + t))
    ! T(3,1)
    amp2(2)=(8d0*(s + t)*(2*s**2 + 2*s*t + t**2))/(s*t**2)
    IF(iqqx.EQ.4.OR.iqqx.EQ.2)THEN
       amp2tot=amp2(1)
       amp2(1)=amp2(2)
       amp2(2)=amp2tot
    ENDIF
    amp2tot=0d0
    DO i=1,2
       amp2tot=amp2tot+amp2(i)
    ENDDO
    amp2col=0
    j_col=0
    IF(amp2tot.LE.0d0)THEN
       ! randomly choosing
       DO i=1,2
          amp2col=amp2col+1d0/6d0
          IF(rcol.LE.amp2col)THEN
             j_col=i
             EXIT
          ENDIF
       ENDDO
    ELSE
       ! choosing following amp2
       DO i=1,2
          amp2col=amp2col+amp2(i)/amp2tot
          IF(rcol.LE.amp2col)THEN
             j_col=i
             EXIT
          ENDIF
       ENDDO
    ENDIF
    DO i=1,4
       icolun_dps(icol,i,1)=icol_gq_save(iqqx,j_col,i,2)
       icolun_dps(icol,i,2)=icol_gq_save(iqqx,j_col,i,1)
    ENDDO
    RETURN
  END SUBROUTINE CF_gq_gq

  SUBROUTINE ME2_q1q1_q2q2(wme)
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(OUT)::wme
    REAL(KIND(1d0))::s,t
    s=s12_jj
    t=t13_jj
    wme=(4d0*Gstrong**4*(s**2 + 2*s*t + 2*t**2))/(9.*s**2)
    RETURN
  END SUBROUTINE ME2_q1q1_q2q2

  SUBROUTINE CF_q1q1_q2q2(iqqx,icol)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::iqqx ! iqqx, 1 q+qx, 2 qx+q
    INTEGER,INTENT(IN)::icol
    INTEGER::init=0,i
    SAVE init
    INTEGER,DIMENSION(2,4,2)::icol_qq_save
    SAVE icol_qq_save
    IF(init.EQ.0)THEN
       ! iqqx=1
       ! q+qx
       icol_qq_save(1,1,1)=1
       icol_qq_save(1,2,1)=0
       icol_qq_save(1,3,1)=0
       icol_qq_save(1,4,1)=2
       icol_qq_save(1,1,2)=0
       icol_qq_save(1,2,2)=2
       icol_qq_save(1,3,2)=1
       icol_qq_save(1,4,2)=0
       ! iqqx=2
       ! qx+q
       icol_qq_save(2,1,1)=0
       icol_qq_save(2,2,1)=1
       icol_qq_save(2,3,1)=0
       icol_qq_save(2,4,1)=2
       icol_qq_save(2,1,2)=2
       icol_qq_save(2,2,2)=0
       icol_qq_save(2,3,2)=1
       icol_qq_save(2,4,2)=0
       init=1
    ENDIF
    !IF(icol.EQ.1)THEN
    !   offset=0
    !ELSE
    !   offset=4
    !ENDIF
    DO i=1,4
       icolun_dps(icol,i,1:2)=icol_qq_save(iqqx,i,1:2)
    ENDDO
    RETURN
  END SUBROUTINE CF_q1q1_q2q2

  SUBROUTINE ME2_q1q1x_q1q1x(iqqx,wme)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::iqqx ! iqqx, 1 q+qx, 2 qx+q
    REAL(KIND(1d0)),INTENT(OUT)::wme
    REAL(KIND(1d0))::s,t
    s=s12_jj
    IF(iqqx.EQ.1)THEN
       t=t13_jj
    ELSE
       t=-t13_jj-s12_jj
    ENDIF
    wme=(8d0*Gstrong**4*(3*s**4 + 2*s**3*t + s**2*t**2 + &
         2*s*t**3 + 3*t**4))/(27.*s**2*t**2)
    RETURN
  END SUBROUTINE ME2_q1q1x_q1q1x

  SUBROUTINE CF_q1q1x_q1q1x(iqqx,icol,rcol)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::iqqx ! iqqx, 1 q+qx, 2 qx+q 
    INTEGER,INTENT(IN)::icol
    REAL(KIND(1d0)),INTENT(IN)::rcol
    INTEGER::init=0,offset,i,j_col
    SAVE init
    INTEGER,DIMENSION(2,2,4,2)::icol_qq_save
    SAVE icol_qq_save
    REAL(KIND(1d0))::s,t,amp2tot,amp2col
    REAL(KIND(1d0)),DIMENSION(2)::amp2
    IF(init.EQ.0)THEN
       ! iqqx=1
       ! q+qx
       ! delta(1,3)delta(2,4)
       icol_qq_save(1,1,1,1)=1
       icol_qq_save(1,1,2,1)=0
       icol_qq_save(1,1,3,1)=0
       icol_qq_save(1,1,4,1)=2
       icol_qq_save(1,1,1,2)=0
       icol_qq_save(1,1,2,2)=2
       icol_qq_save(1,1,3,2)=1
       icol_qq_save(1,1,4,2)=0
       ! delta(1,2)delta(3,4)
       icol_qq_save(1,2,1,1)=1
       icol_qq_save(1,2,2,1)=0
       icol_qq_save(1,2,3,1)=0
       icol_qq_save(1,2,4,1)=2
       icol_qq_save(1,2,1,2)=0
       icol_qq_save(1,2,2,2)=1
       icol_qq_save(1,2,3,2)=2
       icol_qq_save(1,2,4,2)=0
       ! iqqx=2
       ! qx+q
       ! delta(1,4)delta(2,3)
       icol_qq_save(2,1,1,1)=0
       icol_qq_save(2,1,2,1)=1
       icol_qq_save(2,1,3,1)=0
       icol_qq_save(2,1,4,1)=2
       icol_qq_save(2,1,1,2)=2
       icol_qq_save(2,1,2,2)=0
       icol_qq_save(2,1,3,2)=1
       icol_qq_save(2,1,4,2)=0
       ! delta(1,2)delta(3,4)
       icol_qq_save(2,2,1,1)=0
       icol_qq_save(2,2,2,1)=1
       icol_qq_save(2,2,3,1)=0
       icol_qq_save(2,2,4,1)=2
       icol_qq_save(2,2,1,2)=1
       icol_qq_save(2,2,2,2)=0
       icol_qq_save(2,2,3,2)=2
       icol_qq_save(2,2,4,2)=0
       init=1
    ENDIF
    !IF(icol.EQ.1)THEN
    !   offset=0
    !ELSE
    !   offset=4
    !ENDIF
    s=s12_jj
    IF(iqqx.EQ.1)THEN
       t=t13_jj
    ELSE
       t=-t13_jj-s12_jj
    ENDIF
    ! no initial avarage and final symmetry
    ! no Gstrong**4
    ! delta(1,3)delta(2,4)
    amp2(1)=(2d0*(s**2 + 2*s*t + 2*t**2))/s**2
    ! delta(1,2)delta(3,4)
    amp2(2)=(2d0*(2*s**2 + 2*s*t + t**2))/t**2
    amp2tot=0d0
    DO i=1,2
       amp2tot=amp2tot+amp2(i)
    ENDDO
    amp2col=0
    j_col=0
    IF(amp2tot.LE.0d0)THEN
       ! randomly choosing
       DO i=1,2
          amp2col=amp2col+1d0/6d0
          IF(rcol.LE.amp2col)THEN
             j_col=i
             EXIT
          ENDIF
       ENDDO
    ELSE
       ! choosing following amp2
       DO i=1,2
          amp2col=amp2col+amp2(i)/amp2tot
          IF(rcol.LE.amp2col)THEN
             j_col=i
             EXIT
          ENDIF
       ENDDO
    ENDIF
    DO i=1,4
       icolun_dps(icol,i,1:2)=icol_qq_save(iqqx,j_col,i,1:2)
    ENDDO
    RETURN
  END SUBROUTINE CF_q1q1x_q1q1x

  SUBROUTINE ME2_q1q1_q1q1(wme)
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(OUT)::wme
    REAL(KIND(1d0))::s,t
    s=s12_jj
    t=t13_jj
    wme=(4d0*Gstrong**4*(3*s**4 + 10*s**3*t + 13*s**2*t**2 +& 
         6*s*t**3 + 3*t**4))/(27.*t**2*(s + t)**2)
    RETURN
  END SUBROUTINE ME2_q1q1_q1q1

  SUBROUTINE CF_q1q1_q1q1(iqqx,icol,rcol)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::iqqx ! iqqx, 1 q+q, 2 qx+qx
    INTEGER,INTENT(IN)::icol
    REAL(KIND(1d0)),INTENT(IN)::rcol
    INTEGER::init=0,offset,i,j_col
    SAVE init
    INTEGER,DIMENSION(2,2,4,2)::icol_qq_save
    SAVE icol_qq_save
    REAL(KIND(1d0))::s,t,amp2tot,amp2col
    REAL(KIND(1d0)),DIMENSION(2)::amp2
    IF(init.EQ.0)THEN
       ! iqqx=1
       ! q+q
       ! delta(1,4)delta(2,3)
       icol_qq_save(1,1,1,1)=1
       icol_qq_save(1,1,2,1)=2
       icol_qq_save(1,1,3,1)=0
       icol_qq_save(1,1,4,1)=0
       icol_qq_save(1,1,1,2)=0
       icol_qq_save(1,1,2,2)=0
       icol_qq_save(1,1,3,2)=2
       icol_qq_save(1,1,4,2)=1
       ! delta(1,3)delta(2,4)
       icol_qq_save(1,2,1,1)=1
       icol_qq_save(1,2,2,1)=2
       icol_qq_save(1,2,3,1)=0
       icol_qq_save(1,2,4,1)=0
       icol_qq_save(1,2,1,2)=0
       icol_qq_save(1,2,2,2)=0
       icol_qq_save(1,2,3,2)=1
       icol_qq_save(1,2,4,2)=2
       ! iqqx=2
       ! qx+qx
       ! delta(1,4)delta(2,3)
       icol_qq_save(2,1,1,1)=0
       icol_qq_save(2,1,2,1)=0
       icol_qq_save(2,1,3,1)=1
       icol_qq_save(2,1,4,1)=2
       icol_qq_save(2,1,1,2)=2
       icol_qq_save(2,1,2,2)=1
       icol_qq_save(2,1,3,2)=0
       icol_qq_save(2,1,4,2)=0
       ! delta(1,3)delta(2,4)
       icol_qq_save(2,2,1,1)=0
       icol_qq_save(2,2,2,1)=0
       icol_qq_save(2,2,3,1)=1
       icol_qq_save(2,2,4,1)=2
       icol_qq_save(2,2,1,2)=1
       icol_qq_save(2,2,2,2)=2
       icol_qq_save(2,2,3,2)=0
       icol_qq_save(2,2,4,2)=0
       init=1
    ENDIF
    !IF(icol.EQ.1)THEN
    !   offset=0
    !ELSE
    !   offset=4
    !ENDIF
    s=s12_jj
    t=t13_jj
    ! no initial avarage and final symmetry
    ! no Gstrong**4
    ! delta(1,4)delta(2,3)
    amp2(1)=(2d0*(2*s**2 + 2*s*t + t**2))/t**2
    ! delta(1,3)delta(2,4)
    amp2(2)=(2d0*(s**2 + t**2))/(s + t)**2
    amp2tot=0d0
    DO i=1,2
       amp2tot=amp2tot+amp2(i)
    ENDDO
    amp2col=0
    j_col=0
    IF(amp2tot.LE.0d0)THEN
       ! randomly choosing
       DO i=1,2
          amp2col=amp2col+1d0/6d0
          IF(rcol.LE.amp2col)THEN
             j_col=i
             EXIT
          ENDIF
       ENDDO
    ELSE
       ! choosing following amp2
       DO i=1,2
          amp2col=amp2col+amp2(i)/amp2tot
          IF(rcol.LE.amp2col)THEN
             j_col=i
             EXIT
          ENDIF
       ENDDO
    ENDIF
    DO i=1,4
       icolun_dps(icol,i,1:2)=icol_qq_save(iqqx,j_col,i,1:2)
    ENDDO
    RETURN
  END SUBROUTINE CF_q1q1_q1q1

  SUBROUTINE ME2_q1q2_q1q2(iqqx,wme)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::iqqx ! iqqx, 1 q1+q2, 2 q2+q1, 3 q1+q2x, 4 q2x+q1, 5 q1x+q2, 6 q2+q1x, 7 q1x+q2x, 8 q2x+q1x 
    REAL(KIND(1d0)),INTENT(OUT)::wme
    REAL(KIND(1d0))::s,t
    s=s12_jj
    IF(MOD(iqqx,2).EQ.0)THEN
       t=-t13_jj-s12_jj
    ELSE
       t=t13_jj
    ENDIF
    wme=(4d0*Gstrong**4*(2*s**2 + 2*s*t + t**2))/(9.*t**2)
    RETURN
  END SUBROUTINE ME2_q1q2_q1q2

  SUBROUTINE CF_q1q2_q1q2(iqqx,icol)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::iqqx ! iqqx, 1 q1+q2, 2 q2+q1, 3 q1+q2x, 4 q2x+q1, 5 q1x+q2, 6 q2+q1x, 7 q1x+q2x, 8 q2x+q1x
    INTEGER,INTENT(IN)::icol
    INTEGER::init=0,offset,i
    SAVE init
    INTEGER,DIMENSION(8,4,2)::icol_qq_save
    SAVE icol_qq_save
    IF(init.EQ.0)THEN
       ! iqqx=1
       ! q1+q2
       ! delta(1,4)delta(2,3)
       icol_qq_save(1,1,1)=1
       icol_qq_save(1,2,1)=2
       icol_qq_save(1,3,1)=0
       icol_qq_save(1,4,1)=0
       icol_qq_save(1,1,2)=0
       icol_qq_save(1,2,2)=0
       icol_qq_save(1,3,2)=2
       icol_qq_save(1,4,2)=1
       ! iqqx=2
       ! q2+q1
       ! delta(1,3)delta(2,4)
       icol_qq_save(2,1,1)=1
       icol_qq_save(2,2,1)=2
       icol_qq_save(2,3,1)=0
       icol_qq_save(2,4,1)=0
       icol_qq_save(2,1,2)=0
       icol_qq_save(2,2,2)=0
       icol_qq_save(2,3,2)=1
       icol_qq_save(2,4,2)=2
       ! iqqx=3
       ! q1+q2x
       ! delta(1,2)delta(3,4)
       icol_qq_save(3,1,1)=1
       icol_qq_save(3,2,1)=0
       icol_qq_save(3,3,1)=0
       icol_qq_save(3,4,1)=2
       icol_qq_save(3,1,2)=0
       icol_qq_save(3,2,2)=1
       icol_qq_save(3,3,2)=2
       icol_qq_save(3,4,2)=0
       ! iqqx=4
       ! q2x+q1
       ! delta(1,2)delta(3,4)
       icol_qq_save(4,1,1)=0
       icol_qq_save(4,2,1)=1
       icol_qq_save(4,3,1)=0
       icol_qq_save(4,4,1)=2
       icol_qq_save(4,1,2)=1
       icol_qq_save(4,2,2)=0
       icol_qq_save(4,3,2)=2
       icol_qq_save(4,4,2)=0
       ! iqqx=5
       ! q1x+q2
       ! delta(1,2)delta(3,4)
       icol_qq_save(5,1,1)=0
       icol_qq_save(5,2,1)=1
       icol_qq_save(5,3,1)=2
       icol_qq_save(5,4,1)=0
       icol_qq_save(5,1,2)=1
       icol_qq_save(5,2,2)=0
       icol_qq_save(5,3,2)=0
       icol_qq_save(5,4,2)=2
       ! iqqx=6
       ! q2+q1x
       ! delta(1,2)delta(3,4)
       icol_qq_save(6,1,1)=1
       icol_qq_save(6,2,1)=0
       icol_qq_save(6,3,1)=2
       icol_qq_save(6,4,1)=0
       icol_qq_save(6,1,2)=0
       icol_qq_save(6,2,2)=1
       icol_qq_save(6,3,2)=0
       icol_qq_save(6,4,2)=2
       ! iqqx=7
       ! q1x+q2x
       ! delta(1,4)delta(2,3)
       icol_qq_save(7,1,1)=0
       icol_qq_save(7,2,1)=0
       icol_qq_save(7,3,1)=1
       icol_qq_save(7,4,1)=2
       icol_qq_save(7,1,2)=2
       icol_qq_save(7,2,2)=1
       icol_qq_save(7,3,2)=0
       icol_qq_save(7,4,2)=0
       ! iqqx=8
       ! q2x+q1x
       ! delta(1,3)delta(2,4)
       icol_qq_save(8,1,1)=0
       icol_qq_save(8,2,1)=0
       icol_qq_save(8,3,1)=1
       icol_qq_save(8,4,1)=2
       icol_qq_save(8,1,2)=1
       icol_qq_save(8,2,2)=2
       icol_qq_save(8,3,2)=0
       icol_qq_save(8,4,2)=0
       init=1
    ENDIF
    !IF(icol.EQ.1)THEN
    !   offset=0
    !ELSE
    !   offset=4
    !ENDIF
    DO i=1,4
       icolun_dps(icol,i,1:2)=icol_qq_save(iqqx,i,1:2)
    ENDDO
    RETURN
  END SUBROUTINE CF_q1q2_q1q2
END MODULE ME2_pp_jj
