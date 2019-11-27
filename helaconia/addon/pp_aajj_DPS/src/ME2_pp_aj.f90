MODULE ME2_pp_aj
  USE pp_aajj_dps_global
  USE Helac_Global
  IMPLICIT NONE
CONTAINS
  SUBROUTINE ME2_qq_ag(idq,wme)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::idq
    REAL(KIND(1d0)),INTENT(OUT)::wme
    REAL(KIND(1d0))::s,t
    INTEGER::init=0,i
    REAL(KIND(1d0)),DIMENSION(5)::charge
    SAVE init,charge
    IF(init.EQ.0)THEN
       ! d
       charge(1)=-1d0/3d0
       ! u
       charge(2)=2d0/3d0
       ! s
       charge(3)=-1d0/3d0
       ! c
       charge(4)=2d0/3d0
       ! b
       charge(5)=-1d0/3d0
       init=1
    ENDIF
    s=s12_aj
    t=t13_aj
    ! it has been averaged over initial color and helicity
    ! and final state symmetry factor
    wme=(-32d0*EL**2*Gstrong**2*(s**2 + 2*s*t + 2*t**2))/(81.*t*(s + t))
    wme=wme*charge(abs(idq))**2/(4d0/9d0)
    RETURN
  END SUBROUTINE ME2_qq_ag

  SUBROUTINE CF_qq_ag(iqqx,icol)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::iqqx ! iqqx: 1 qqx, 2 qxq
    INTEGER,INTENT(IN)::icol
    REAL(KIND(1d0))::s,t
    INTEGER::init=0,offset,i
    INTEGER,DIMENSION(2,4,2)::icol_ag_save
    SAVE init,icol_ag_save
    IF(init.EQ.0)THEN
       ! iqqx = 1
       ! qqx
       icol_ag_save(1,1,1)=1
       icol_ag_save(1,1,2)=0
       icol_ag_save(1,2,1)=0
       icol_ag_save(1,2,2)=2
       icol_ag_save(1,3,1:2)=0
       icol_ag_save(1,4,2)=1
       icol_ag_save(1,4,1)=2
       ! iqqx = 2
       ! qxq
       icol_ag_save(2,1,1)=0
       icol_ag_save(2,1,2)=1
       icol_ag_save(2,2,1)=2
       icol_ag_save(2,2,2)=0
       icol_ag_save(2,3,1:2)=0
       icol_ag_save(2,4,2)=2
       icol_ag_save(2,4,1)=1
       init=1
    END IF
    !IF(icol.EQ.1)THEN
    !   offset=0
    !ELSE
    !   offset=4
    !ENDIF
    DO i=1,4
       icolun_dps(icol,i,1:2)=icol_ag_save(iqqx,i,1:2)
    ENDDO
    RETURN
  END SUBROUTINE CF_qq_ag

  SUBROUTINE ME2_qg_aq(idq,iqqx,wme)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::idq
    INTEGER,INTENT(IN)::iqqx ! iqqx: 1, q+g, 2, qx+g, 3, g+q, 4, g+qx
    REAL(KIND(1d0)),INTENT(OUT)::wme
    REAL(KIND(1d0))::s,t
    INTEGER::init=0,i
    REAL(KIND(1d0)),DIMENSION(5)::charge
    SAVE init,charge
    IF(init.EQ.0)THEN
       ! d
       charge(1)=-1d0/3d0
       ! u
       charge(2)=2d0/3d0
       ! s
       charge(3)=-1d0/3d0
       ! c
       charge(4)=2d0/3d0
       ! b
       charge(5)=-1d0/3d0
       init=1
    ENDIF
    s=s12_aj
    IF(iqqx.EQ.1.OR.iqqx.EQ.2)THEN
       t=t13_aj
    ELSEIF(iqqx.EQ.3.OR.iqqx.EQ.4)THEN
       t=-t13_aj-s12_aj
    ELSE
       WRITE(*,*)"ERROR:Unknown iqqx=",iqqx
       STOP
    ENDIF
    wme=(-4d0*EL**2*Gstrong**2*(s**2 + t**2))/(27.*s*t)
    wme=wme*charge(abs(idq))**2/(4d0/9d0)
    RETURN
  END SUBROUTINE ME2_qg_aq

  SUBROUTINE CF_qg_aq(iqqx,icol)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::iqqx ! iqqx: 1, q+g, 2, qx+g, 3, g+q, 4, g+qx
    INTEGER,INTENT(IN)::icol
    INTEGER::init=0,offset,i
    INTEGER,DIMENSION(4,4,2)::icol_aq_save
    SAVE init,icol_aq_save
    IF(init.EQ.0)THEN
       ! iqqx = 1
       ! q+g
       icol_aq_save(1,1,1)=1
       icol_aq_save(1,1,2)=0
       icol_aq_save(1,2,1)=2
       icol_aq_save(1,2,2)=1
       icol_aq_save(1,3,1:2)=0
       icol_aq_save(1,4,2)=0
       icol_aq_save(1,4,1)=2
       ! iqqx = 2
       ! qx+g
       icol_aq_save(2,1,1)=0
       icol_aq_save(2,1,2)=1
       icol_aq_save(2,2,1)=1
       icol_aq_save(2,2,2)=2
       icol_aq_save(2,3,1:2)=0
       icol_aq_save(2,4,2)=0
       icol_aq_save(2,4,1)=2
       ! iqqx = 3
       ! g+q
       icol_aq_save(3,1,1)=1
       icol_aq_save(3,1,2)=2
       icol_aq_save(3,2,1)=2
       icol_aq_save(3,2,2)=0
       icol_aq_save(3,3,1:2)=0
       icol_aq_save(3,4,2)=1
       icol_aq_save(3,4,1)=0
       ! iqqx = 4
       ! g+qx
       icol_aq_save(4,1,1)=1
       icol_aq_save(4,1,2)=2
       icol_aq_save(4,2,1)=0
       icol_aq_save(4,2,2)=1
       icol_aq_save(4,3,1:2)=0
       icol_aq_save(4,4,2)=0
       icol_aq_save(4,4,1)=2
       init=1
    ENDIF
    !IF(icol.EQ.1)THEN
    !   offset=0
    !ELSE
    !   offset=4
    !ENDIF
    DO i=1,4
       icolun_dps(icol,i,1:2)=icol_aq_save(iqqx,i,1:2)
    ENDDO
    RETURN
  END SUBROUTINE CF_qg_aq

END MODULE ME2_pp_aj
