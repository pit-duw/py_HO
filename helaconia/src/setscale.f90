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
CASE(3)
   IF(emep_ISR_shower.EQ.0)THEN
      STOP "Scale cannot be 3 (only 0,1,2) without e-e+ shower."
   ENDIF
   scale1=DSQRT(Q2OUT_QEDPS)
CASE(4)
   ! HT/2
   mas=0d0
   DO ij=3,n
      mas=mas+DSQRT(parmas(ifl(ij))**2+Phegas_pmom(ij,1)**2+&
           Phegas_pmom(ij,2)**2)
   ENDDO
   scale1=mas/2d0
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
scale1=MAX(scale1,1d0)
scale1=scale1*scalefactor

END SUBROUTINE qcdscale

END MODULE Setscale
