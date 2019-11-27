MODULE Projectors
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                  !
!   This module provides projectors for heavy quarkoniums                          !
!                                                                                  !
!   Author: Hua-Sheng Shao                                                         !
!   Organization: Physics School of Peking University                              !
!   Creation Time: 20:20 22.Sep.2011                                               !
!   Last Modification: 20:20 22.Sep.2011                                           !
!                                                                                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
USE Helac_wavef
IMPLICIT NONE
CONTAINS
SUBROUTINE SpinProj(zq1,zq2,hel,prj) 
! The Spin Projector For Production
IMPLICIT NONE
COMPLEX(KIND(1d0)),DIMENSION(1:5),INTENT(IN)::zq1,zq2
REAL(KIND(1d0))::E1,E2,m1,m2
! hel(1),hel(2) represents the spin of Quark1 and Anti-Quark2,respectively
! hel(3),hel(4) are the total spin and the third components of the Hadron
INTEGER,DIMENSION(4),INTENT(IN)::hel
COMPLEX(KIND(1d0)),DIMENSION(1:4)::zy1,zy2,zy
COMPLEX(KIND(1d0)),DIMENSION(1:5)::zP
COMPLEX(KIND(1d0)),INTENT(OUT)::prj
COMPLEx(KIND(1d0))::apr,zzP0,zzPZ
E1=(DREAL(zq1(1)+zq1(2)))/2
E2=(DREAL(zq2(1)+zq2(2)))/2
m1=DREAL(zq1(5))
m2=DREAL(zq2(5))
zP(1:5)=zq1(1:5)+zq2(1:5)
zzP0=DCMPLX(E1+E2)
zzPZ=DCMPLX(DIMAG(zP(1)))
apr=DCMPLX(E1+E2)
IF(hel(1).GE.0)THEN
	CALL Helac_uplus(zy1,zq1)
ELSE
	CALL Helac_uminus(zy1,zq1)
ENDIF 
IF(hel(2).GE.0)THEN
	CALL Helac_vbplus(zy2,zq2)
ELSE
	CALL Helac_vbminus(zy2,zq2)
ENDIF
SELECT CASE(hel(3))
CASE(0)
	prj=apr*(zy2(1)*zy1(1)+zy2(2)*zy1(2)-zy2(3)*zy1(3)-zy2(4)*zy1(4))&
	-zzP0*(zy2(1)*zy1(3)+zy2(2)*zy1(4)-zy2(3)*zy1(1)-zy2(4)*zy1(2))&
	-zzPZ*(zy2(1)*zy1(3)-zy2(2)*zy1(4)+zy2(3)*zy1(1)-zy2(4)*zy1(2))&
	-zP(3)*(zy2(2)*zy1(3)+zy2(4)*zy1(1))-zP(4)*(zy2(1)*zy1(4)+zy2(3)*zy1(2))
CASE(1)
	SELECT CASE(hel(4))
	CASE(-1)
		CALL Helac_eminus(zy,zP)
	CASE(1)
		CALL Helac_eplus(zy,zP)
	CASE(0)
		CALL Helac_elong(zy,zP)
	END SELECT
	prj=apr*(zy(1)*(-zy2(1)*zy1(3)-zy2(4)*zy1(2))+zy(2)*(-zy2(2)*zy1(4)-zy2(3)*zy1(1))&
	+zy(3)*(zy2(4)*zy1(1)-zy2(2)*zy1(3))+zy(4)*(-zy2(1)*zy1(4)+zy2(3)*zy1(2)))&
	+DCMPLX(DREAL(zP(2)))*(zy(1)*(zy2(1)*zy1(1)+zy2(4)*zy1(4))+zy(3)*zy2(2)*zy1(1)&
	-zy(4)*zy2(3)*zy1(4))+DCMPLX(DREAL(zP(1)))*(zy(2)*(zy2(2)*zy1(2)+zy2(3)*zy1(3))&
	-zy(3)*zy2(4)*zy1(3)+zy(4)*zy2(1)*zy1(2))+zP(4)*(-zy(1)*zy2(1)*zy1(2)+zy(2)*zy2(3)*zy1(4)&
	-zy(3)*(zy2(2)*zy1(2)+zy2(4)*zy1(4)))+zP(3)*(zy(1)*zy2(4)*zy1(3)-zy(2)*zy2(2)*zy1(1)&
	-zy(4)*(zy2(1)*zy1(1)+zy2(3)*zy1(3)))
END SELECT 
prj=prj/apr/SQRT(8*(E1+m1)*(E2+m2))
END SUBROUTINE SpinProj

FUNCTION ClebschGordan(LL,SS,JJ)
IMPLICIT NONE
INTEGER,DIMENSION(2),INTENT(IN)::LL,SS,JJ
COMPLEX(KIND(1d0))::ClebschGordan
IF((LL(2)+SS(2).NE.JJ(2)).OR.(LL(1).LT.0).OR.(SS(1).LT.0).OR.(JJ(1).LT.0) &
.OR.(ABS(LL(2)).GT.LL(1)).OR.(ABS(SS(2)).GT.SS(1)).OR.(ABS(JJ(2)).GT.JJ(1))&
.OR.(LL(1)+SS(1).LT.JJ(1)).OR.(ABS(LL(1)-SS(1)).GT.JJ(1)))THEN
	ClebschGordan=DCMPLX(0)
	RETURN
END IF
IF(LL(1).EQ.0)THEN
	IF((SS(1).NE.JJ(1)).OR.(SS(2).NE.JJ(2)))THEN
		ClebschGordan=DCMPLX(0)
		RETURN
	ELSE
		ClebschGordan=DCMPLX(1)
		RETURN
	ENDIF
ENDIF
IF(SS(1).EQ.0)THEN
	IF((LL(1).NE.JJ(1)).OR.(LL(2).NE.JJ(2)))THEN
		ClebschGordan=DCMPLX(0)
		RETURN
	ELSE
		ClebschGordan=DCMPLX(1)
		RETURN
	ENDIF
ENDIF
IF((LL(1).EQ.1).AND.(SS(1).EQ.1))THEN
	SELECT CASE(JJ(1))
	CASE(1)
		SELECT CASE(LL(2))
		CASE(0)
			ClebschGordan=DCMPLX(-SS(2)*1d0/DSQRT(2d0))
			RETURN
		CASE(1,-1)
			ClebschGordan=DCMPLX(DBLE(LL(2))*1d0/DSQRT(2d0))
			RETURN
		END SELECT
	CASE(2)
		SELECT CASE(LL(2))
		CASE(1,-1)
			IF(SS(2).EQ.LL(2))THEN
				ClebschGordan=DCMPLX(1d0)
				RETURN
			ENDIF
			IF(SS(2).EQ.0)THEN
				ClebschGordan=DCMPLX(1d0/DSQRT(2d0))
				RETURN
			ENDIF
			IF(SS(2)+LL(2).EQ.0)THEN
				ClebschGordan=DCMPLX(1d0/DSQRT(6d0))
				RETURN
			ENDIF
		CASE(0)
			SELECT CASE(SS(2))
			CASE(1,-1)
				ClebschGordan=DCMPLX(1d0/DSQRT(2d0))
				RETURN
			CASE(0)
				ClebschGordan=DCMPLX(DSQRT(2d0/3d0))
				RETURN
			END SELECT
		END SELECT
	CASE(0)
		SELECT CASE(LL(2))
		CASE(1,-1)
			ClebschGordan=DCMPLX(1d0/DSQRT(3d0))
			RETURN
		CASE(0)
			ClebschGordan=DCMPLX(-1d0/DSQRT(3d0))
			RETURN
		END SELECT
	END SELECT
ENDIF
PRINT *,"Wrong in ClebschGordan ! STOP !"
STOP
END FUNCTION ClebschGordan
END MODULE Projectors
