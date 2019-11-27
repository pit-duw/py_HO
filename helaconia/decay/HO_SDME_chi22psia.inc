INTEGER,INTENT(IN)::ihelm
INTEGER,INTENT(OUT)::iheld1,iheld2
REAL(KIND(1d0))::rrr,rrr2,costh,sinth,phi,prob,delta0,delta1,fmax
INTEGER::init=0,i,j,k
SAVE init
LOGICAL::lexit
REAL(KIND(1d0)),PARAMETER::pi=3.1415926535897932385d0
REAL(KIND(1d0)),DIMENSION(6)::fac1
SAVE fac1
REAL(KIND(1d0)),DIMENSION(6)::fac2
SAVE fac2
REAL(KIND(1d0)),DIMENSION(6)::fac3
SAVE fac3
REAL(KIND(1d0)),DIMENSION(6)::fac4
SAVE fac4
REAL(KIND(1d0)),DIMENSION(6)::fac5
SAVE fac5
IF(init.EQ.0)THEN
	delta0=1d0+2d0*aJ21*(DSQRT(5d0)*aJ22+2d0*aJ23)&
	+4d0*aJ22*(aJ22+DSQRT(5d0)*aJ23)+3d0*aJ23**2
	delta0=delta0/10d0
	delta1=9d0+6d0*aJ21*(DSQRT(5d0)*aJ22-4d0*aJ23)&
	-4d0*aJ22*(aJ22+2d0*DSQRT(5d0)*aJ23)+7d0*aJ23**2
	delta1=delta1/30d0
	fac1(1)=delta0/2.
	fac1(2)=(1 - delta1)/2.
	fac1(3)=0.5
	fac1(4)=(1 + delta1)/2.
	fac1(5)=1 - delta0/2.
	fac1(6)=1
	fac2(1)=delta0/2.
	fac2(2)=(1 - delta1)/2.
	fac2(3)=0.5
	fac2(4)=(1 + delta1)/2.
	fac2(5)=1 - delta0/2.
	fac2(6)=1
	fac3(1)=delta0/2.
	fac3(2)=(1 - delta1)/2.
	fac3(3)=0.5
	fac3(4)=(1 + delta1)/2.
	fac3(5)=1 - delta0/2.
	fac3(6)=1
	fac4(1)=delta0/2.
	fac4(2)=(1 - delta1)/2.
	fac4(3)=0.5
	fac4(4)=(1 + delta1)/2.
	fac4(5)=1 - delta0/2.
	fac4(6)=1
	fac5(1)=delta0/2.
	fac5(2)=(1 - delta1)/2.
	fac5(3)=0.5
	fac5(4)=(1 + delta1)/2.
	fac5(5)=1 - delta0/2.
	fac5(6)=1
	init=1
ENDIF
SELECT CASE(ihelm)
CASE(-2)
	rrr=Helac_rnmy(0)
	k=1
	DO i=1,6
		IF(rrr.LE.fac1(i))THEN
			k=i
			EXIT
		ENDIF
		IF(i.EQ.6)THEN
			k=6
			EXIT
		ENDIF
	ENDDO
	SELECT CASE(k)
	CASE(1)
	! decay via (3*delta0*Sin(Subscript(\[Theta],13))**4)/16.
		iheld1=-1
		iheld2=-1
		fmax=1d0
		lexit=.FALSE.
		DO WHILE(.NOT.lexit)
			costh=2d0*Helac_rnmy(0)-1d0
			prob=(1d0-costh**2)**2/fmax
			rrr2=Helac_rnmy(0)
			IF(prob.GT.rrr2)lexit=.TRUE.
		ENDDO
	CASE(2)
	! decay via -((-1 + delta0 + delta1)*Cos(Subscript(\[Theta],13)/2.)**8)/2.
		iheld1=-1
		iheld2=1
		fmax=1d0
		lexit=.FALSE.
		DO WHILE(.NOT.lexit)
			costh=2d0*Helac_rnmy(0)-1d0
			prob=(1d0+costh)**4/16d0/fmax
			rrr2=Helac_rnmy(0)
			IF(prob.GT.rrr2)lexit=.TRUE.
		ENDDO
	CASE(3)
	! decay via (delta1*Sin(Subscript(\[Theta],13)/2.)**4*Sin(Subscript(\[Theta],13))**2)/2.
		iheld1=0
		iheld2=-1
		fmax=0.43d0
		lexit=.FALSE.
		DO WHILE(.NOT.lexit)
			costh=2d0*Helac_rnmy(0)-1d0
			prob=(1d0-costh)**2*(1d0-costh**2)/4d0/fmax
			rrr2=Helac_rnmy(0)
			IF(prob.GT.rrr2)lexit=.TRUE.
		ENDDO
	CASE(4)
	! decay via (delta1*Csc(Subscript(\[Theta],13)/2.)**4*Sin(Subscript(\[Theta],13))**6)/32.
		iheld1=0
		iheld2=1
		fmax=6.75d0
		lexit=.FALSE.
		DO WHILE(.NOT.lexit)
			costh=2d0*Helac_rnmy(0)-1d0
			prob=4d0*(1d0-costh)*(1d0+costh)**3/fmax
			rrr2=Helac_rnmy(0)
			IF(prob.GT.rrr2)lexit=.TRUE.
		ENDDO
	CASE(5)
	! decay via -((-1 + delta0 + delta1)*Sin(Subscript(\[Theta],13)/2.)**8)/2.
		iheld1=1
		iheld2=-1
		fmax=1d0
		lexit=.FALSE.
		DO WHILE(.NOT.lexit)
			costh=2d0*Helac_rnmy(0)-1d0
			prob=(1d0-costh)**4/16d0/fmax
			rrr2=Helac_rnmy(0)
			IF(prob.GT.rrr2)lexit=.TRUE.
		ENDDO
	CASE(6)
	! decay via (3*delta0*Sin(Subscript(\[Theta],13))**4)/16.
		iheld1=1
		iheld2=1
		fmax=1d0
		lexit=.FALSE.
		DO WHILE(.NOT.lexit)
			costh=2d0*Helac_rnmy(0)-1d0
			prob=(1d0-costh**2)**2/fmax
			rrr2=Helac_rnmy(0)
			IF(prob.GT.rrr2)lexit=.TRUE.
		ENDDO
	END SELECT
CASE(-1)
	rrr=Helac_rnmy(0)
	k=1
	DO i=1,6
		IF(rrr.LE.fac2(i))THEN
			k=i
			EXIT
		ENDIF
		IF(i.EQ.6)THEN
			k=6
			EXIT
		ENDIF
	ENDDO
	SELECT CASE(k)
	CASE(1)
	! decay via (3*delta0*Sin(2*Subscript(\[Theta],13))**2)/16.
		iheld1=-1
		iheld2=-1
		fmax=1d0
		lexit=.FALSE.
		DO WHILE(.NOT.lexit)
			costh=2d0*Helac_rnmy(0)-1d0
			prob=4d0*(1d0-costh**2)*costh**2/fmax
			rrr2=Helac_rnmy(0)
			IF(prob.GT.rrr2)lexit=.TRUE.
		ENDDO
	CASE(2)
	! decay via -((-1 + delta0 + delta1)*Csc(Subscript(\[Theta],13)/2.)**4*Sin(Subscript(\[Theta],13))**6)/32.
		iheld1=-1
		iheld2=1
		fmax=6.75d0
		lexit=.FALSE.
		DO WHILE(.NOT.lexit)
			costh=2d0*Helac_rnmy(0)-1d0
			prob=4d0*(1d0-costh)*(1d0+costh)**3/fmax
			rrr2=Helac_rnmy(0)
			IF(prob.GT.rrr2)lexit=.TRUE.
		ENDDO
	CASE(3)
	! decay via (delta1*(1 + 2*Cos(Subscript(\[Theta],13)))**2*Sin(Subscript(\[Theta],13)/2.)**4)/2.
		iheld1=0
		iheld2=-1
		fmax=1d0
		lexit=.FALSE.
		DO WHILE(.NOT.lexit)
			costh=2d0*Helac_rnmy(0)-1d0
			prob=(1d0+2d0*costh)**2*(1d0-costh)**2/4d0/fmax
			rrr2=Helac_rnmy(0)
			IF(prob.GT.rrr2)lexit=.TRUE.
		ENDDO
	CASE(4)
	! decay via (delta1*Cos(Subscript(\[Theta],13)/2.)**4*(1 - 2*Cos(Subscript(\[Theta],13)))**2)/2.
		iheld1=0
		iheld2=1
		fmax=1d0
		lexit=.FALSE.
		DO WHILE(.NOT.lexit)
			costh=2d0*Helac_rnmy(0)-1d0
			prob=(1d0+costh)**2/4d0*(1d0-2d0*costh)**2/fmax
			rrr2=Helac_rnmy(0)
			IF(prob.GT.rrr2)lexit=.TRUE.
		ENDDO
	CASE(5)
	! decay via -((-1 + delta0 + delta1)*Sin(Subscript(\[Theta],13)/2.)**4*Sin(Subscript(\[Theta],13))**2)/2.
		iheld1=1
		iheld2=-1
		fmax=0.43d0
		lexit=.FALSE.
		DO WHILE(.NOT.lexit)
			costh=2d0*Helac_rnmy(0)-1d0
			prob=(1d0-costh)**2*(1d0-costh**2)/4d0/fmax
			rrr2=Helac_rnmy(0)
			IF(prob.GT.rrr2)lexit=.TRUE.
		ENDDO
	CASE(6)
	! decay via (3*delta0*Sin(2*Subscript(\[Theta],13))**2)/16.
		iheld1=1
		iheld2=1
		fmax=1d0
		lexit=.FALSE.
		DO WHILE(.NOT.lexit)
			costh=2d0*Helac_rnmy(0)-1d0
			prob=4d0*(1d0-costh**2)*costh**2/fmax
			rrr2=Helac_rnmy(0)
			IF(prob.GT.rrr2)lexit=.TRUE.
		ENDDO
	END SELECT
CASE(0)
	rrr=Helac_rnmy(0)
	k=1
	DO i=1,6
		IF(rrr.LE.fac3(i))THEN
			k=i
			EXIT
		ENDIF
		IF(i.EQ.6)THEN
			k=6
			EXIT
		ENDIF
	ENDDO
	SELECT CASE(k)
	CASE(1)
	! decay via (delta0*(1 + 3*Cos(2*Subscript(\[Theta],13)))**2)/32.
		iheld1=-1
		iheld2=-1
		fmax=16d0
		lexit=.FALSE.
		DO WHILE(.NOT.lexit)
			costh=2d0*Helac_rnmy(0)-1d0
			prob=(-2d0+6d0*costh**2)**2/fmax
			rrr2=Helac_rnmy(0)
			IF(prob.GT.rrr2)lexit=.TRUE.
		ENDDO
	CASE(2)
	! decay via (-3*(-1 + delta0 + delta1)*Sin(Subscript(\[Theta],13))**4)/16.
		iheld1=-1
		iheld2=1
		fmax=1d0
		lexit=.FALSE.
		DO WHILE(.NOT.lexit)
			costh=2d0*Helac_rnmy(0)-1d0
			prob=(1d0-costh**2)**2/fmax
			rrr2=Helac_rnmy(0)
			IF(prob.GT.rrr2)lexit=.TRUE.
		ENDDO
	CASE(3)
	! decay via (3*delta1*Sin(2*Subscript(\[Theta],13))**2)/16.
		iheld1=0
		iheld2=-1
		fmax=1d0
		lexit=.FALSE.
		DO WHILE(.NOT.lexit)
			costh=2d0*Helac_rnmy(0)-1d0
			prob=4d0*(1d0-costh**2)*costh**2/fmax
			rrr2=Helac_rnmy(0)
			IF(prob.GT.rrr2)lexit=.TRUE.
		ENDDO
	CASE(4)
	! decay via (3*delta1*Sin(2*Subscript(\[Theta],13))**2)/16.
		iheld1=0
		iheld2=1
		fmax=1d0
		lexit=.FALSE.
		DO WHILE(.NOT.lexit)
			costh=2d0*Helac_rnmy(0)-1d0
			prob=4d0*(1d0-costh**2)*costh**2/fmax
			rrr2=Helac_rnmy(0)
			IF(prob.GT.rrr2)lexit=.TRUE.
		ENDDO
	CASE(5)
	! decay via (-3*(-1 + delta0 + delta1)*Sin(Subscript(\[Theta],13))**4)/16.
		iheld1=1
		iheld2=-1
		fmax=1d0
		lexit=.FALSE.
		DO WHILE(.NOT.lexit)
			costh=2d0*Helac_rnmy(0)-1d0
			prob=(1d0-costh**2)**2/fmax
			rrr2=Helac_rnmy(0)
			IF(prob.GT.rrr2)lexit=.TRUE.
		ENDDO
	CASE(6)
	! decay via (delta0*(1 + 3*Cos(2*Subscript(\[Theta],13)))**2)/32.
		iheld1=1
		iheld2=1
		fmax=16d0
		lexit=.FALSE.
		DO WHILE(.NOT.lexit)
			costh=2d0*Helac_rnmy(0)-1d0
			prob=(-2d0+6d0*costh**2)**2/fmax
			rrr2=Helac_rnmy(0)
			IF(prob.GT.rrr2)lexit=.TRUE.
		ENDDO
	END SELECT
CASE(1)
	rrr=Helac_rnmy(0)
	k=1
	DO i=1,6
		IF(rrr.LE.fac4(i))THEN
			k=i
			EXIT
		ENDIF
		IF(i.EQ.6)THEN
			k=6
			EXIT
		ENDIF
	ENDDO
	SELECT CASE(k)
	CASE(1)
	! decay via (3*delta0*Sin(2*Subscript(\[Theta],13))**2)/16.
		iheld1=-1
		iheld2=-1
		fmax=1d0
		lexit=.FALSE.
		DO WHILE(.NOT.lexit)
			costh=2d0*Helac_rnmy(0)-1d0
			prob=4d0*(1d0-costh**2)*costh**2/fmax
			rrr2=Helac_rnmy(0)
			IF(prob.GT.rrr2)lexit=.TRUE.
		ENDDO
	CASE(2)
	! decay via -((-1 + delta0 + delta1)*Sin(Subscript(\[Theta],13)/2.)**4*Sin(Subscript(\[Theta],13))**2)/2.
		iheld1=-1
		iheld2=1
		fmax=0.43d0
		lexit=.FALSE.
		DO WHILE(.NOT.lexit)
			costh=2d0*Helac_rnmy(0)-1d0
			prob=(1d0-costh)**2*(1d0-costh**2)/4d0/fmax
			rrr2=Helac_rnmy(0)
			IF(prob.GT.rrr2)lexit=.TRUE.
		ENDDO
	CASE(3)
	! decay via (delta1*Cos(Subscript(\[Theta],13)/2.)**4*(1 - 2*Cos(Subscript(\[Theta],13)))**2)/2.
		iheld1=0
		iheld2=-1
		fmax=1d0
		lexit=.FALSE.
		DO WHILE(.NOT.lexit)
			costh=2d0*Helac_rnmy(0)-1d0
			prob=(1d0+costh)**2/4d0*(1d0-2d0*costh)**2/fmax
			rrr2=Helac_rnmy(0)
			IF(prob.GT.rrr2)lexit=.TRUE.
		ENDDO
	CASE(4)
	! decay via (delta1*(1 + 2*Cos(Subscript(\[Theta],13)))**2*Sin(Subscript(\[Theta],13)/2.)**4)/2.
		iheld1=0
		iheld2=1
		fmax=1d0
		lexit=.FALSE.
		DO WHILE(.NOT.lexit)
			costh=2d0*Helac_rnmy(0)-1d0
			prob=(1d0+2d0*costh)**2*(1d0-costh)**2/4d0/fmax
			rrr2=Helac_rnmy(0)
			IF(prob.GT.rrr2)lexit=.TRUE.
		ENDDO
	CASE(5)
	! decay via -((-1 + delta0 + delta1)*Csc(Subscript(\[Theta],13)/2.)**4*Sin(Subscript(\[Theta],13))**6)/32.
		iheld1=1
		iheld2=-1
		fmax=6.75d0
		lexit=.FALSE.
		DO WHILE(.NOT.lexit)
			costh=2d0*Helac_rnmy(0)-1d0
			prob=4d0*(1d0-costh)*(1d0+costh)**3/fmax
			rrr2=Helac_rnmy(0)
			IF(prob.GT.rrr2)lexit=.TRUE.
		ENDDO
	CASE(6)
	! decay via (3*delta0*Sin(2*Subscript(\[Theta],13))**2)/16.
		iheld1=1
		iheld2=1
		fmax=1d0
		lexit=.FALSE.
		DO WHILE(.NOT.lexit)
			costh=2d0*Helac_rnmy(0)-1d0
			prob=4d0*(1d0-costh**2)*costh**2/fmax
			rrr2=Helac_rnmy(0)
			IF(prob.GT.rrr2)lexit=.TRUE.
		ENDDO
	END SELECT
CASE(2)
	rrr=Helac_rnmy(0)
	k=1
	DO i=1,6
		IF(rrr.LE.fac5(i))THEN
			k=i
			EXIT
		ENDIF
		IF(i.EQ.6)THEN
			k=6
			EXIT
		ENDIF
	ENDDO
	SELECT CASE(k)
	CASE(1)
	! decay via (3*delta0*Sin(Subscript(\[Theta],13))**4)/16.
		iheld1=-1
		iheld2=-1
		fmax=1d0
		lexit=.FALSE.
		DO WHILE(.NOT.lexit)
			costh=2d0*Helac_rnmy(0)-1d0
			prob=(1d0-costh**2)**2/fmax
			rrr2=Helac_rnmy(0)
			IF(prob.GT.rrr2)lexit=.TRUE.
		ENDDO
	CASE(2)
	! decay via -((-1 + delta0 + delta1)*Sin(Subscript(\[Theta],13)/2.)**8)/2.
		iheld1=-1
		iheld2=1
		fmax=1d0
		lexit=.FALSE.
		DO WHILE(.NOT.lexit)
			costh=2d0*Helac_rnmy(0)-1d0
			prob=(1d0-costh)**4/16d0/fmax
			rrr2=Helac_rnmy(0)
			IF(prob.GT.rrr2)lexit=.TRUE.
		ENDDO
	CASE(3)
	! decay via (delta1*Csc(Subscript(\[Theta],13)/2.)**4*Sin(Subscript(\[Theta],13))**6)/32.
		iheld1=0
		iheld2=-1
		fmax=6.75d0
		lexit=.FALSE.
		DO WHILE(.NOT.lexit)
			costh=2d0*Helac_rnmy(0)-1d0
			prob=4d0*(1d0-costh)*(1d0+costh)**3/fmax
			rrr2=Helac_rnmy(0)
			IF(prob.GT.rrr2)lexit=.TRUE.
		ENDDO
	CASE(4)
	! decay via (delta1*Sin(Subscript(\[Theta],13)/2.)**4*Sin(Subscript(\[Theta],13))**2)/2.
		iheld1=0
		iheld2=1
		fmax=0.43d0
		lexit=.FALSE.
		DO WHILE(.NOT.lexit)
			costh=2d0*Helac_rnmy(0)-1d0
			prob=(1d0-costh)**2*(1d0-costh**2)/4d0/fmax
			rrr2=Helac_rnmy(0)
			IF(prob.GT.rrr2)lexit=.TRUE.
		ENDDO
	CASE(5)
	! decay via -((-1 + delta0 + delta1)*Cos(Subscript(\[Theta],13)/2.)**8)/2.
		iheld1=1
		iheld2=-1
		fmax=1d0
		lexit=.FALSE.
		DO WHILE(.NOT.lexit)
			costh=2d0*Helac_rnmy(0)-1d0
			prob=(1d0+costh)**4/16d0/fmax
			rrr2=Helac_rnmy(0)
			IF(prob.GT.rrr2)lexit=.TRUE.
		ENDDO
	CASE(6)
	! decay via (3*delta0*Sin(Subscript(\[Theta],13))**4)/16.
		iheld1=1
		iheld2=1
		fmax=1d0
		lexit=.FALSE.
		DO WHILE(.NOT.lexit)
			costh=2d0*Helac_rnmy(0)-1d0
			prob=(1d0-costh**2)**2/fmax
			rrr2=Helac_rnmy(0)
			IF(prob.GT.rrr2)lexit=.TRUE.
		ENDDO
	END SELECT
END SELECT
sinth=DSQRT(ABS(1d0-costh**2))
phi=2*pi*Helac_rnmy(0)
