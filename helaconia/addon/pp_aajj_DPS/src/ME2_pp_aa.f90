MODULE ME2_pp_aa
  USE pp_aajj_dps_global
  USE Helac_Global
  IMPLICIT NONE
CONTAINS
  SUBROUTINE ME2_qq_aa(idq,wme)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::idq
    REAL(KIND(1d0)),INTENT(OUT)::wme
    REAL(KIND(1d0))::s,t
    INTEGER::init=0
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
    s=s12_aa
    t=t13_aa
    ! it has been averaged over initial color and helicity
    ! and final state symmetry factor
    wme=(-16d0*EL**4*(s**2 + 2*s*t + 2*t**2))/(243.*t*(s + t))
    wme=wme*charge(abs(idq))**4/(16d0/81d0)
    RETURN
  END SUBROUTINE ME2_qq_aa

  SUBROUTINE CF_qq_aa(iqqx,icol)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::iqqx ! iqqx: 1 qqx, 2 qxq
    INTEGER,INTENT(IN)::icol
    INTEGER::init=0,offset,i
    INTEGER,DIMENSION(2,4,2)::icol_aa_save
    SAVE init,icol_aa_save
    IF(init.EQ.0)THEN
       ! iqqx = 1
       ! qqx
       icol_aa_save(1,1,1)=1
       icol_aa_save(1,1,2)=0
       icol_aa_save(1,2,1)=0
       icol_aa_save(1,2,2)=1
       icol_aa_save(1,3,1:2)=0
       icol_aa_save(1,4,1:2)=0
       ! iqqx = 2
       ! qxq
       icol_aa_save(2,1,1)=0
       icol_aa_save(2,1,2)=1
       icol_aa_save(2,2,1)=1
       icol_aa_save(2,2,2)=0
       icol_aa_save(2,3,1:2)=0
       icol_aa_save(2,4,1:2)=0
       init=1
    ENDIF
    !IF(icol.EQ.1)THEN
    !   offset=0
    !ELSE
    !   offset=4
    !ENDIF
    DO i=1,4
       icolun_dps(icol,i,1:2)=icol_aa_save(iqqx,i,1:2)
    ENDDO
    RETURN
  END SUBROUTINE CF_qq_aa

  SUBROUTINE ME2_gg_aa(wme)
    ! without top loop
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(OUT)::wme
    REAL(KIND(1d0))::s,t
    INTEGER::i,j
    REAL(KIND(1d0)),DIMENSION(14,14)::pairmatrix
    COMPLEX(KIND(1d0)),EXTERNAL::camp_gg_aa
    s=s12_aa
    t=t13_aa
    pairmatrix(1,1)=t**2*(s + t)**2
    pairmatrix(1,2)=t**2*(s + t)**2
    pairmatrix(1,3)=s*t*(s + t)
    pairmatrix(1,4)=s*t*(s + t)
    pairmatrix(1,5)=t**2*(s + t)
    pairmatrix(1,6)=t**2*(s + t)
    pairmatrix(1,7)=t*(s + t)**2
    pairmatrix(1,8)=t**2*(s + t)
    pairmatrix(1,9)=t**2*(s + t)
    pairmatrix(1,10)=t*(s + t)
    pairmatrix(1,11)=t*(s + t)**2
    pairmatrix(1,12)=t*(s + t)
    pairmatrix(1,13)=(t**2*(s + t)**2)/s
    pairmatrix(1,14)=t*(s + t)
    pairmatrix(2,2)=t**2*(s + t)**2
    pairmatrix(2,3)=s*t*(s + t)
    pairmatrix(2,4)=s*t*(s + t)
    pairmatrix(2,5)=t**2*(s + t)
    pairmatrix(2,6)=t**2*(s + t)
    pairmatrix(2,7)=t*(s + t)**2
    pairmatrix(2,8)=t**2*(s + t)
    pairmatrix(2,9)=t**2*(s + t)
    pairmatrix(2,10)=t*(s + t)
    pairmatrix(2,11)=t*(s + t)**2
    pairmatrix(2,12)=t*(s + t)
    pairmatrix(2,13)=(t**2*(s + t)**2)/s
    pairmatrix(2,14)=t*(s + t)
    pairmatrix(3,3)=2*s**2
    pairmatrix(3,4)=2*s**2
    pairmatrix(3,5)=s*t
    pairmatrix(3,6)=s*t
    pairmatrix(3,7)=s*(s + t)
    pairmatrix(3,8)=s*t
    pairmatrix(3,9)=s*t
    pairmatrix(3,10)=s
    pairmatrix(3,11)=s*(s + t)
    pairmatrix(3,12)=s
    pairmatrix(3,13)=t*(s + t)
    pairmatrix(3,14)=2*s
    pairmatrix(4,4)=2*s**2
    pairmatrix(4,5)=s*t
    pairmatrix(4,6)=s*t
    pairmatrix(4,7)=s*(s + t)
    pairmatrix(4,8)=s*t
    pairmatrix(4,9)=s*t
    pairmatrix(4,10)=s
    pairmatrix(4,11)=s*(s + t)
    pairmatrix(4,12)=s
    pairmatrix(4,13)=t*(s + t)
    pairmatrix(4,14)=2*s
    pairmatrix(5,5)=2*t**2
    pairmatrix(5,6)=2*t**2
    pairmatrix(5,7)=t*(s + t)
    pairmatrix(5,8)=t**2
    pairmatrix(5,9)=t**2
    pairmatrix(5,10)=t
    pairmatrix(5,11)=t*(s + t)
    pairmatrix(5,12)=2*t
    pairmatrix(5,13)=(t**2*(s + t))/s
    pairmatrix(5,14)=t
    pairmatrix(6,6)=2*t**2
    pairmatrix(6,7)=t*(s + t)
    pairmatrix(6,8)=t**2
    pairmatrix(6,9)=t**2
    pairmatrix(6,10)=t
    pairmatrix(6,11)=t*(s + t)
    pairmatrix(6,12)=2*t
    pairmatrix(6,13)=(t**2*(s + t))/s
    pairmatrix(6,14)=t
    pairmatrix(7,7)=2*(s + t)**2
    pairmatrix(7,8)=t*(s + t)
    pairmatrix(7,9)=t*(s + t)
    pairmatrix(7,10)=2*(s + t)
    pairmatrix(7,11)=(s + t)**2
    pairmatrix(7,12)=s + t
    pairmatrix(7,13)=(t*(s + t)**2)/s
    pairmatrix(7,14)=s + t
    pairmatrix(8,8)=2*t**2
    pairmatrix(8,9)=2*t**2
    pairmatrix(8,10)=2*t
    pairmatrix(8,11)=t*(s + t)
    pairmatrix(8,12)=t
    pairmatrix(8,13)=(t**2*(s + t))/s
    pairmatrix(8,14)=t
    pairmatrix(9,9)=2*t**2
    pairmatrix(9,10)=2*t
    pairmatrix(9,11)=t*(s + t)
    pairmatrix(9,12)=t
    pairmatrix(9,13)=(t**2*(s + t))/s
    pairmatrix(9,14)=t
    pairmatrix(10,10)=4d0
    pairmatrix(10,11)=s + t
    pairmatrix(10,12)=2d0
    pairmatrix(10,13)=(t*(s + t))/s
    pairmatrix(10,14)=2d0
    pairmatrix(11,11)=2*(s + t)**2
    pairmatrix(11,12)=2*(s + t)
    pairmatrix(11,13)=(t*(s + t)**2)/s
    pairmatrix(11,14)=s + t
    pairmatrix(12,12)=4d0
    pairmatrix(12,13)=(t*(s + t))/s
    pairmatrix(12,14)=2d0
    pairmatrix(13,13)=(2*t**2*(s + t)**2)/s**2
    pairmatrix(13,14)=(2*t*(s + t))/s
    pairmatrix(14,14)=4d0
    DO i=2,14
       DO j=1,i-1
          pairmatrix(i,j)=pairmatrix(j,i)
       ENDDO
    END DO
    wme=0d0
    DO i=1,14
       DO j=1,14
          wme=wme+camp_gg_aa(i,s,t)*pairmatrix(i,j)*CONJG(camp_gg_aa(j,s,t))
       ENDDO
    ENDDO
    wme=wme*8d0/2d0/4d0/64d0*EL**4*Gstrong**4
    RETURN
  END SUBROUTINE ME2_gg_aa

  SUBROUTINE CF_gg_aa(icol)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::icol
    INTEGER::init=0,offset,i,j
    INTEGER,DIMENSION(4,2)::icol_aa_save
    SAVE init,icol_aa_save
    IF(init.EQ.0)THEN
       icol_aa_save(1,1)=1
       icol_aa_save(1,2)=2
       icol_aa_save(2,1)=2
       icol_aa_save(2,2)=1
       icol_aa_save(3,1:2)=0
       icol_aa_save(4,1:2)=0
       init=1
    ENDIF
    !IF(icol.EQ.1)THEN
    !   offset=0
    !ELSE
    !   offset=4
    !ENDIF
    DO i=1,4
       icolun_dps(icol,i,1:2)=icol_aa_save(i,1:2)
    ENDDO
    RETURN
  END SUBROUTINE CF_gg_aa

END MODULE ME2_pp_aa
