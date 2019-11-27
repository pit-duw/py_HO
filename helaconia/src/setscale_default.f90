MODULE setscale
USE Helac_Global
IMPLICIT NONE
CONTAINS      
SUBROUTINE qcdscale(scale1)
REAL(KIND=DBL),INTENT(OUT)::scale1
REAL(KIND=DBL)::px,py,ptw2,mas
INTEGER::ij
SELECT CASE(nscheme)
CASE(0)
! fixed scheme
	scale1=fschemevalue
CASE(1)
! Sqrt(pt1**2+m1**2)
	IF(ABS(iflh(3)).GT.100)THEN
	! for quarkonium
		px=Phegas_pmom(3,1)+Phegas_pmom(4,1)
		py=Phegas_pmom(3,2)+Phegas_pmom(4,2)
		ptw2=px**2+py**2
		mas=parmas(ifl(3))+parmas(ifl(4))
	ELSE
		px=Phegas_pmom(3,1)
		py=Phegas_pmom(3,2)
		ptw2=px**2+py**2
		mas=parmas(ifl(3))
	ENDIF
	scale1=DSQRT(mas**2+ptw2)
CASE(2)
! Sqrt(pt1**2+sum(mi,i=3,n)**2) in Jpsi ccbar by F.Maltoni
	mas=0d0
	DO ij=3,n
		mas=mas+parmas(ifl(ij))
	ENDDO
	IF(ABS(iflh(3)).GT.100)THEN
		px=Phegas_pmom(3,1)+Phegas_pmom(4,1)
		py=Phegas_pmom(3,2)+Phegas_pmom(4,2)
	ELSE
		px=Phegas_pmom(3,1)
		py=Phegas_pmom(3,2)
	ENDIF
	ptw2=px**2+py**2
	scale1=DSQRT(mas**2+ptw2)
END SELECT
!
! default Z_mass
!
!scale1=parmas(32)
!Scale1=parmas(11)
!scale1=12.8d0

!
! example 1: pt of particles 3 and 4 (=W) + W_mass ^2
!
!      px=Phegas_pmom(3,1)+Phegas_pmom(4,1)
!      py=Phegas_pmom(3,2)+Pegas_pmom(4,2)
!      ptw2=px**2+py**2
!      scale2=parmas(33)**2+ptw2

!
! example 2: pt of 3 + pt of 4 + 2 times the top mass
!
!      pt1=Phegas_pmom(3,1)**2+Phegas_pmom(3,2)**2
!      pt2=Phegas_pmom(4,1)**2+Phegas_pmom(4,2)**2
!      scale2=2*parmas(11)**2+pt1+pt2
!      scale1=dsqrt(scale2)

scale1=scale1*scalefactor

END SUBROUTINE qcdscale

END MODULE Setscale
