INTEGER,INTENT(IN)::ihelm
INTEGER,INTENT(OUT)::iheld1,iheld2
REAL(KIND(1d0))::rrr,rrr2,costh,sinth,phi
INTEGER::init=0,i,j,k
SAVE init
REAL(KIND(1d0)),PARAMETER::pi=3.1415926535897932385d0
REAL(KIND(1d0)),DIMENSION(4)::fac1
SAVE fac1
REAL(KIND(1d0)),DIMENSION(4)::fac2
SAVE fac2
REAL(KIND(1d0)),DIMENSION(4)::fac3
SAVE fac3
IF(init.EQ.0)THEN
	fac1(1)=(1 + 2*aJ10*aJ11)/4.
	fac1(2)=0.5
	fac1(3)=(3 - 2*aJ10*aJ11)/4.
	fac1(4)=1
	fac2(1)=(1 + 2*aJ10*aJ11)/4.
	fac2(2)=0.5
	fac2(3)=(3 - 2*aJ10*aJ11)/4.
	fac2(4)=1
	fac3(1)=(1 + 2*aJ10*aJ11)/4.
	fac3(2)=0.5
	fac3(3)=(3 - 2*aJ10*aJ11)/4.
	fac3(4)=1
	init=1
ENDIF
SELECT CASE(ihelm)
CASE(-1)
	rrr=Helac_rnmy(0)
	k=1
	DO i=1,4
		IF(rrr.LE.fac1(i))THEN
			k=i
			EXIT
		ENDIF
		IF(i.EQ.4)THEN
			k=4
			EXIT
		ENDIF
	ENDDO
	SELECT CASE(k)
	CASE(1)
	! decay via ((1 + 2*aJ10*aJ11)*Sin(Subscript(\[Theta],13))**2)/8.
		iheld1=-1
		iheld2=-1
		rrr2=2d0*(Helac_rnmy(0)-0.5d0)
		IF(rrr2.LT.-1d0)THEN
			rrr2=-1d0
		ELSEIF(rrr2.GT.1d0)THEN
			rrr2=1d0
		ENDIF
		costh=SOLVEVll00(rrr2)
	CASE(2)
	! decay via -((-1 + 2*aJ10*aJ11)*Sin(Subscript(\[Theta],13)/2.)**4)/4.
		iheld1=0
		iheld2=-1
		rrr2=8d0*Helac_rnmy(0)
		costh=1-rrr2**(1d0/3d0)
	CASE(3)
	! decay via -((-1 + 2*aJ10*aJ11)*Cos(Subscript(\[Theta],13)/2.)**4)/4.
		iheld1=0
		iheld2=1
		rrr2=8d0*Helac_rnmy(0)
		costh=rrr2**(1d0/3d0)-1
	CASE(4)
	! decay via ((1 + 2*aJ10*aJ11)*Sin(Subscript(\[Theta],13))**2)/8.
		iheld1=1
		iheld2=1
		rrr2=2d0*(Helac_rnmy(0)-0.5d0)
		IF(rrr2.LT.-1d0)THEN
			rrr2=-1d0
		ELSEIF(rrr2.GT.1d0)THEN
			rrr2=1d0
		ENDIF
		costh=SOLVEVll00(rrr2)
	END SELECT
CASE(0)
	rrr=Helac_rnmy(0)
	k=1
	DO i=1,4
		IF(rrr.LE.fac2(i))THEN
			k=i
			EXIT
		ENDIF
		IF(i.EQ.4)THEN
			k=4
			EXIT
		ENDIF
	ENDDO
	SELECT CASE(k)
	CASE(1)
	! decay via ((1 + 2*aJ10*aJ11)*Cos(Subscript(\[Theta],13))**2)/4.
		iheld1=-1
		iheld2=-1
		rrr2=2d0/3d0*Helac_rnmy(0)-1d0/3d0
		IF(rrr2.GE.0d0)THEN
		   costh=(rrr2*3d0)**(1d0/3d0)
		ELSE
                   costh=-(-rrr2*3d0)**(1d0/3d0)
                ENDIF
	CASE(2)
	! decay via -((-1 + 2*aJ10*aJ11)*Sin(Subscript(\[Theta],13))**2)/8.
		iheld1=0
		iheld2=-1
		rrr2=2d0*(Helac_rnmy(0)-0.5d0)
		IF(rrr2.LT.-1d0)THEN
			rrr2=-1d0
		ELSEIF(rrr2.GT.1d0)THEN
			rrr2=1d0
		ENDIF
		costh=SOLVEVll00(rrr2)
	CASE(3)
	! decay via -((-1 + 2*aJ10*aJ11)*Sin(Subscript(\[Theta],13))**2)/8.
		iheld1=0
		iheld2=1
		rrr2=2d0*(Helac_rnmy(0)-0.5d0)
		IF(rrr2.LT.-1d0)THEN
			rrr2=-1d0
		ELSEIF(rrr2.GT.1d0)THEN
			rrr2=1d0
		ENDIF
		costh=SOLVEVll00(rrr2)
	CASE(4)
	! decay via ((1 + 2*aJ10*aJ11)*Cos(Subscript(\[Theta],13))**2)/4.
		iheld1=1
		iheld2=1
		rrr2=2d0/3d0*Helac_rnmy(0)-1d0/3d0
                IF(rrr2.GE.0d0)THEN
		   costh=(rrr2*3d0)**(1d0/3d0)
                ELSE
                   costh=-(-rrr2*3d0)**(1d0/3d0)
                ENDIF
	END SELECT
CASE(1)
	rrr=Helac_rnmy(0)
	k=1
	DO i=1,4
		IF(rrr.LE.fac3(i))THEN
			k=i
			EXIT
		ENDIF
		IF(i.EQ.4)THEN
			k=4
			EXIT
		ENDIF
	ENDDO
	SELECT CASE(k)
	CASE(1)
	! decay via ((1 + 2*aJ10*aJ11)*Sin(Subscript(\[Theta],13))**2)/8.
		iheld1=-1
		iheld2=-1
		rrr2=2d0*(Helac_rnmy(0)-0.5d0)
		IF(rrr2.LT.-1d0)THEN
			rrr2=-1d0
		ELSEIF(rrr2.GT.1d0)THEN
			rrr2=1d0
		ENDIF
		costh=SOLVEVll00(rrr2)
	CASE(2)
	! decay via -((-1 + 2*aJ10*aJ11)*Cos(Subscript(\[Theta],13)/2.)**4)/4.
		iheld1=0
		iheld2=-1
		rrr2=8d0*Helac_rnmy(0)
		costh=rrr2**(1d0/3d0)-1
	CASE(3)
	! decay via -((-1 + 2*aJ10*aJ11)*Sin(Subscript(\[Theta],13)/2.)**4)/4.
		iheld1=0
		iheld2=1
		rrr2=8d0*Helac_rnmy(0)
		costh=1-rrr2**(1d0/3d0)
	CASE(4)
	! decay via ((1 + 2*aJ10*aJ11)*Sin(Subscript(\[Theta],13))**2)/8.
		iheld1=1
		iheld2=1
		rrr2=2d0*(Helac_rnmy(0)-0.5d0)
		IF(rrr2.LT.-1d0)THEN
			rrr2=-1d0
		ELSEIF(rrr2.GT.1d0)THEN
			rrr2=1d0
		ENDIF
		costh=SOLVEVll00(rrr2)
	END SELECT
END SELECT
sinth=DSQRT(ABS(1d0-costh**2))
phi=2*pi*Helac_rnmy(0)
