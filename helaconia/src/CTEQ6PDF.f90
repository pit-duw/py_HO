MODULE CTEQ6PDF
!============================================================================
!                             April 10, 2002, v6.01
!                             February 23, 2003, v6.1
!
!   Ref[1]: "New Generation of Parton Distributions with Uncertainties from Global QCD Analysis"
!       By: J. Pumplin, D.R. Stump, J.Huston, H.L. Lai, P. Nadolsky, W.K. Tung
!       JHEP 0207:012(2002), hep-ph/0201195
!
!   Ref[2]: "Inclusive Jet Production, Parton Distributions, and the Search for New Physics"
!       By : D. Stump, J. Huston, J. Pumplin, W.K. Tung, H.L. Lai, S. Kuhlmann, J. Owens
!       hep-ph/0303013
!
!   This package contains
!   (1) 4 standard sets of CTEQ6 PDF's (CTEQ6M, CTEQ6D, CTEQ6L, CTEQ6L1) ;
!   (2) 40 up/down sets (with respect to CTEQ6M) for uncertainty studies from Ref[1];
!   (3) updated version of the above: CTEQ6.1M and its 40 up/down eigenvector sets from Ref[2].
!
!  The CTEQ6.1M set provides a global fit that is almost equivalent in every respect
!  to the published CTEQ6M, Ref[1], although some parton distributions (e.g., the gluon)
!  may deviate from CTEQ6M in some kinematic ranges by amounts that are well within the
!  specified uncertainties.
!  The more significant improvements of the new version are associated with some of the
!  40 eigenvector sets, which are made more symmetrical and reliable in (3), compared to (2).
!
!  Details about calling convention are:
! ---------------------------------------------------------------------------
!  Iset   PDF-set     Description       Alpha_s(Mz)**Lam4  Lam5   Table_File
! ===========================================================================
! Standard, "best-fit", sets:
! --------------------------
!   1    CTEQ6M   Standard MSbar scheme   0.118     326   226    cteq6m.tbl
!   2    CTEQ6D   Standard DIS scheme     0.118     326   226    cteq6d.tbl
!   3    CTEQ6L   Leading Order           0.118**   326** 226    cteq6l.tbl
!   4    CTEQ6L1  Leading Order           0.130**   215** 165    cteq6l1.tbl
! ============================================================================
! For uncertainty calculations using eigenvectors of the Hessian:
! ---------------------------------------------------------------
!     central + 40 up/down sets along 20 eigenvector directions
!                             -----------------------------
!                Original version, Ref[1]:  central fit: CTEQ6M (=CTEQ6M.00)
!                             -----------------------
!  1xx  CTEQ6M.xx  +/- sets               0.118     326   226    cteq6m1xx.tbl
!        where xx = 01-40: 01/02 corresponds to +/- for the 1st eigenvector, ... etc.
!        e.g. 100      is CTEQ6M.00 (=CTEQ6M),
!             101/102 are CTEQ6M.01/02, +/- sets of 1st eigenvector, ... etc.
!                              -----------------------
!                Updated version, Ref[2]:  central fit: CTEQ6.1M (=CTEQ61.00)
!                              -----------------------
!  2xx  CTEQ61.xx  +/- sets               0.118     326   226    ctq61.xx.tbl
!        where xx = 01-40: 01/02 corresponds to +/- for the 1st eigenvector, ... etc.
!        e.g. 200      is CTEQ61.00 (=CTEQ6.1M),
!             201/202 are CTEQ61.01/02, +/- sets of 1st eigenvector, ... etc.
! ===========================================================================
!   ** ALL fits are obtained by using the same coupling strength
!   \alpha_s(Mz)=0.118 and the NLO running \alpha_s formula, except CTEQ6L1
!   which uses the LO running \alpha_s and its value determined from the fit.
!   For the LO fits, the evolution of the PDF and the hard cross sections are
!   calculated at LO.  More detailed discussions are given in the references.
!
!   The table grids are generated for 10^-6 < x < 1 and 1.3 < Q < 10,000 (GeV).
!   PDF values outside of the above range are returned using extrapolation.
!   Lam5 (Lam4) represents Lambda value (in MeV) for 5 (4) flavors.
!   The matching alpha_s between 4 and 5 flavors takes place at Q=4.5 GeV,
!   which is defined as the bottom quark mass, whenever it can be applied.
!
!   The Table_Files are assumed to be in the working directory.
!
!   Before using the PDF, it is necessary to do the initialization by
!       Call SetCtq6(Iset)
!   where Iset is the desired PDF specified in the above table.
!
!   The function Ctq6Pdf (Iparton, X, Q)
!   returns the parton distribution inside the proton for parton [Iparton]
!   at [X] Bjorken_X and scale [Q] (GeV) in PDF set [Iset].
!   Iparton  is the parton label (5, 4, 3, 2, 1, 0, -1, ......, -5)
!                            for (b, c, s, d, u, g, u_bar, ..., b_bar),
!
!   For detailed information on the parameters used, e.q. quark masses,
!   QCD Lambda, ... etc.,  see info lines at the beginning of the
!   Table_Files.
!
!   These programs, as provided, are in double precision.  By removing the
!   "Implicit Double Precision" lines, they can also be run in single
!   precision.
!
!   If you have detailed questions concerning these CTEQ6 distributions,
!   or if you find problems/bugs using this package, direct inquires to
!   Pumplin@pa.msu.edu or Tung@pa.msu.edu.
!
!===========================================================================
IMPLICIT NONE
SAVE
!INTEGER,PARAMETER::SGL=SELECTED_REAL_KIND(p=6)
INTEGER,PARAMETER,PRIVATE::DBL=SELECTED_REAL_KIND(p=13)
REAL(KIND=DBL),PARAMETER,PRIVATE::pi=3.1415926535897932385d0
REAL(KIND=DBL),PRIVATE::Alambda,AL,Qini, Qmax, Xmin
INTEGER,PRIVATE::Nx, Nt, NfMx, Nfl, Iorder
INTEGER,PARAMETER,PRIVATE::MXX = 96, MXQ = 20, MXF = 5
INTEGER,PARAMETER,PRIVATE::MXPQX = (MXF + 3) * MXQ * MXX
INTEGER,PARAMETER,PRIVATE::MXQX= MXQ * MXX
REAL(KIND=DBL),DIMENSION(0:MXX),PRIVATE::XV
REAL(KIND=DBL),DIMENSION(0:MXQ),PRIVATE::TV
REAL(KIND=DBL),DIMENSION(MXPQX),PRIVATE::UPD
REAL(KIND=DBL),DIMENSION(6),PRIVATE::Amass
CONTAINS
FUNCTION Ctq6Pdf_f90 (Iparton, X, Q)
REAL(KIND=DBL)::Ctq6Pdf_f90
REAL(KIND=DBL),INTENT(IN)::X,Q
INTEGER,INTENT(IN)::Iparton
LOGICAL::Warn=.TRUE.
SAVE Warn

IF(X .LT. 0d0 .OR. X .GT. 1d0) THEN
   WRITE(*,*) 'X out of range in Ctq6Pdf_f90: ', X
   STOP
ENDIF
IF (Q.LT.Alambda)THEN
   WRITE(*,*) 'Q out of range in Ctq6Pdf_f90: ', Q
   STOP
ENDIF

IF((Iparton.LT. -NfMx .OR. Iparton .GT. NfMx)) THEN
   IF (Warn) THEN
!        put a warning for calling extra flavor.
	  Warn = .FALSE.
	  WRITE(*,*) 'Warning: Iparton out of range in Ctq6Pdf_f90: ', Iparton
   ENDIF
   CTQ6PDF_f90 = 0d0
   RETURN
ENDIF

Ctq6Pdf_f90 = PartonX6 (Iparton, X, Q)
IF(Ctq6Pdf_f90.LT.0.d0)  Ctq6Pdf_f90 = 0.d0

!                             ********************
END FUNCTION Ctq6Pdf_f90

SUBROUTINE SETCTQ6f90 (Iset)
IMPLICIT NONE
INTEGER,INTENT(IN)::Iset
INTEGER,PARAMETER::Isetmax0=5
CHARACTER(len=6),DIMENSION(Isetmax0)::Flnm=&
                   (/ 'cteq6m', 'cteq6d', 'cteq6l', 'cteq6l','ctq61.'/)
CHARACTER(40)::Tablefile='test.tbl'
CHARACTER(len=24),PARAMETER::cteq_data_dir="./pdf/cteq/"
CHARACTER(3)::nn
INTEGER::Isetold=-987,Isetmin0=1,Isetmin1=100,Isetmax1=140,Isetmin2=200,Isetmax2=240
INTEGER::IU,status
SAVE

!             If data file not initialized, do so.
IF(Iset.NE.Isetold) THEN
	IU= NextUn()
    IF (Iset.GE.Isetmin0 .AND. Iset.LE.3) THEN
		Tablefile=TRIM(cteq_data_dir)//Flnm(Iset)//'.tbl'
    ELSEIF (Iset.EQ.4) THEN
		Tablefile=TRIM(cteq_data_dir)//Flnm(Iset)//'1.tbl'
    ELSEIF (Iset.GE.Isetmin1 .AND. Iset.LE.Isetmax1) THEN
		WRITE(nn,'(I3)') Iset
        Tablefile=TRIM(cteq_data_dir)//Flnm(1)//nn//'.tbl'
    ELSEIF (Iset.GE.Isetmin2 .AND. Iset.LE.Isetmax2) THEN
		WRITE(nn,'(I3)') Iset
        Tablefile=TRIM(cteq_data_dir)//Flnm(5)//nn(2:3)//'.tbl'
    ELSE
        PRINT *, 'Invalid Iset number in SetCtq6 :', Iset
        STOP
    ENDIF
    OPEN(IU, File=Tablefile, Status='OLD', Err=100)
	CALL ReadTbl (IU)
	CLOSE (IU)
    Isetold=Iset
ENDIF
RETURN

100 PRINT *, ' Data file ', Tablefile, ' cannot be opened '//'in SetCtq6f90!!'
STOP
!                             ********************
END SUBROUTINE SETCTQ6f90 

SUBROUTINE ReadTbl (Nu)
IMPLICIT NONE
INTEGER,INTENT(IN)::Nu
CHARACTER(len=80)::Line
INTEGER::Nblk,Npts,I,Iq,IRET
REAL(KIND=DBL)::Dr,Fl

READ(Nu, '(A)') Line
READ(Nu, '(A)') Line
READ(Nu, *) Dr, Fl, Al, (Amass(I),I=1,6)
Iorder = NINT(Dr)
Nfl = NINT(Fl)
Alambda = Al

READ(Nu, '(A)') Line
READ(Nu, *) NX,  NT, NfMx

READ(Nu, '(A)') Line
READ(Nu, *) QINI, QMAX, (TV(I), I =0, NT)

READ(Nu, '(A)') Line
READ(Nu, *) XMIN, (XV(I), I =0, NX)

DO Iq = 0, NT
	TV(Iq) = LOG(LOG (TV(Iq) /Al))
ENDDO
!
! Since quark = anti-quark for nfl>2 at this stage,
! we Read  out only the non-redundent data points
! No of flavors = NfMx (sea) + 1 (gluon) + 2 (valence)

Nblk = (NX+1) * (NT+1)
Npts =  Nblk  * (NfMx+3)
READ  (Nu, '(A)') Line
READ  (Nu, *, IOSTAT=IRET) (UPD(I), I=1,Npts)

RETURN
!                        ****************************
END SUBROUTINE ReadTbl

FUNCTION NextUn()
! Returns an unallocated FORTRAN i/o unit.
INTEGER::NextUn
LOGICAL::EX
INTEGER::N

DO N = 10, 300
   INQUIRE (UNIT=N, OPENED=EX)
   IF (.NOT. EX)THEN
      NextUn = N
      RETURN
   ENDIF
ENDDO
WRITE(*,*) ' There is no available I/O unit. '
STOP
!               *************************
END FUNCTION NextUn

! interpolation with polynominals, correct!
SUBROUTINE POLINT (XA,YA,N,X,Y,DY)
IMPLICIT NONE
!   Adapted from "Numerical Recipes" 
INTEGER,PARAMETER::NMAX=10
REAL(KIND=DBL),INTENT(IN)::X
REAL(KIND=DBL),INTENT(OUT)::Y,DY
INTEGER,INTENT(IN)::N
REAL(KIND=DBL),DIMENSION(N),INTENT(IN)::XA,YA
REAL(KIND=DBL),DIMENSION(NMAX)::C,D
INTEGER::NS,I,M
REAL(KIND=DBL)::DIF,DIFT,HO,HP,W,DEN
NS=1
DIF=DABS(X-XA(1))
DO I=1,N
   DIFT=DABS(X-XA(I))
   IF (DIFT.LT.DIF) THEN
      NS=I
      DIF=DIFT
   ENDIF
   C(I)=YA(I)
   D(I)=YA(I)
ENDDO
Y=YA(NS)
NS=NS-1
DO M=1,N-1
   DO I=1,N-M
      HO=XA(I)-X
      HP=XA(I+M)-X
      W=C(I+1)-D(I)
      DEN=HO-HP
!      IF(DEN.EQ.0.d0)PAUSE
      DEN=W/DEN
      D(I)=HP*DEN
      C(I)=HO*DEN
   ENDDO
   IF (2*NS.LT.N-M)THEN
      DY=C(NS+1)
   ELSE
      DY=D(NS)
      NS=NS-1
   ENDIF
   Y=Y+DY
ENDDO
END SUBROUTINE POLINT

FUNCTION PartonX6 (IPRTN, XX, QQ)
!  Given the parton distribution function in the array U in
!  COMMON / PEVLDT / , this routine interpolates to find
!  the parton distribution at an arbitray point in x and q.
!
IMPLICIT NONE
INTEGER,INTENT(IN)::IPRTN
REAL(KIND=DBL),INTENT(IN)::XX,QQ
REAL(KIND=DBL)::PartonX6
REAL(KIND=DBL)::OneP=1.0001
REAL(KIND=DBL),DIMENSION(0:MXX)::xvpow
REAL(KIND=DBL),DIMENSION(4)::fvec,fij
INTEGER::ientry=0,nqvec=4
REAL(KIND=DBL)::xpow=0.3d0  !**** choice of interpolation variable
REAL(KIND=DBL)::X,Q,tt,ss,svec1,svec2,svec3,svec4,s12,s13,s23,s24,s34,sy2,sy3,&
				t12,tdet,t34,tmp1,tmp2,t24,t13,ty3,tvec1,tvec3,tvec4,tf2,tf3,const1,const2,&
				const3,const4,const5,const6,s1213,s2434,sdet,tmp,tvec2,t23,ty2,fx,Dfx,sf2,sf3,&
				Dfq,ff,g1,g4,h00
INTEGER::i,JLx,JU,JM,Jx,JLq,it,jtmp,Jq,Ip,J1
SAVE ientry,xvpow

! store the powers used for interpolation on first call...
IF(ientry .EQ. 0) THEN
	ientry = 1
    xvpow(0) = 0d0
    DO i = 1, Nx
		xvpow(i) = xv(i)**xpow
    ENDDO
ENDIF

X = XX
Q = QQ
tt = LOG(LOG(Q/Al))

!      -------------    find lower end of interval containing x, i.e.,
!                       get jx such that xv(jx) .le. x .le. xv(jx+1)...
JLx = -1
JU = Nx+1
DO     ! 11
	IF(JU-JLx .GT. 1) THEN
		JM = (JU+JLx) / 2
		IF (X .GE. XV(JM)) THEN
			JLx = JM
        ELSE
            JU = JM
        ENDIF
!         Goto 11
	ELSE
		EXIT
    ENDIF
ENDDO
!                     Ix    0   1   2      Jx  JLx         Nx-2     Nx
!                           |---|---|---|...|---|-x-|---|...|---|---|
!                     x     0  Xmin               x                 1
!
IF(JLx .LE. -1) THEN
	PRINT '(A,1pE12.4)', 'Severe error: x <= 0 in PartonX6! x = ', x
    STOP
ELSEIF (JLx .EQ. 0) THEN
	Jx = 0
ELSEIF (JLx .LE. Nx-2) THEN
!                For interrior points, keep x in the middle, as shown above
	Jx = JLx - 1
ELSEIF (JLx.EQ.Nx-1 .OR. x.LT.OneP) THEN
!                  We tolerate a slight over-shoot of one (OneP=1.00001),
!              perhaps due to roundoff or whatever, but not more than that.
!                                      Keep at least 4 points >= Jx
	Jx = JLx - 2
ELSE
    PRINT '(A,1pE12.4)', 'Severe error: x > 1 in PartonX6! x = ', x
    STOP
ENDIF
!          ---------- Note: JLx uniquely identifies the x-bin; Jx does not.
!                       This is the variable to be interpolated in
ss = x**xpow

IF (JLx.GE.2 .AND. JLx.LE.Nx-2) THEN
!     initiation work for "interior bins": store the lattice points in s...
	svec1 = xvpow(jx)
	svec2 = xvpow(jx+1)
	svec3 = xvpow(jx+2)
	svec4 = xvpow(jx+3)

	s12 = svec1 - svec2
	s13 = svec1 - svec3
	s23 = svec2 - svec3
	s24 = svec2 - svec4
	s34 = svec3 - svec4

	sy2 = ss - svec2
	sy3 = ss - svec3

! constants needed for interpolating in s at fixed t lattice points...
	const1 = s13/s23
	const2 = s12/s23
	const3 = s34/s23
	const4 = s24/s23
	s1213 = s12 + s13
	s2434 = s24 + s34
	sdet = s12*s34 - s1213*s2434
	tmp = sy2*sy3/sdet
	const5 = (s34*sy2-s2434*sy3)*tmp/s12
	const6 = (s1213*sy2-s12*sy3)*tmp/s34

ENDIF

!         --------------Now find lower end of interval containing Q, i.e.,
!                          get jq such that qv(jq) .le. q .le. qv(jq+1)...
JLq = -1
JU = NT+1
DO  ! 12
	IF(JU-JLq .GT. 1) THEN
		JM = (JU+JLq) / 2
		IF (tt .GE. TV(JM)) THEN
			JLq = JM
        ELSE
			JU = JM
        ENDIF
!         Goto 12
	ELSE
		EXIT
	ENDIF
ENDDO

IF(JLq .LE. 0) THEN
	Jq = 0
ELSEIF(JLq .LE. Nt-2) THEN
!keep q in the middle, as shown above
	Jq = JLq - 1
ELSE
!JLq .GE. Nt-1 case:  Keep at least 4 points >= Jq.
	Jq = Nt - 3
ENDIF
!This is the interpolation variable in Q

IF(JLq.GE.1 .AND. JLq.LE.Nt-2) THEN
! store the lattice points in t...
	tvec1 = Tv(jq)
	tvec2 = Tv(jq+1)
	tvec3 = Tv(jq+2)
	tvec4 = Tv(jq+3)

	t12 = tvec1 - tvec2
	t13 = tvec1 - tvec3
	t23 = tvec2 - tvec3
	t24 = tvec2 - tvec4
	t34 = tvec3 - tvec4

	ty2 = tt - tvec2
	ty3 = tt - tvec3

	tmp1 = t12 + t13
	tmp2 = t24 + t34

	tdet = t12*t34 - tmp1*tmp2

ENDIF

! get the pdf function values at the lattice points...

IF(Iprtn .GE. 3) THEN
	Ip = - Iprtn
ELSE
	Ip = Iprtn
ENDIF
jtmp = ((Ip + NfMx)*(NT+1)+(jq-1))*(NX+1)+jx+1

DO it = 1, nqvec
	J1  = jtmp + it*(NX+1)

    IF(Jx .EQ. 0) THEN
!  For the first 4 x points, interpolate x^2*f(x,Q)
!  This applies to the two lowest bins JLx = 0, 1
!  We can not put the JLx.eq.1 bin into the "interrior" section
!  (as we do for q), since Upd(J1) is undefined.
		fij(1) = 0
		fij(2) = Upd(J1+1) * XV(1)**2
		fij(3) = Upd(J1+2) * XV(2)**2
		fij(4) = Upd(J1+3) * XV(3)**2
!
!  Use Polint which allows x to be anywhere w.r.t. the grid

		CALL Polint (XVpow(0), Fij(1), 4, ss, Fx, Dfx)

        IF (x .GT. 0d0)  Fvec(it) =  Fx / x**2
!  Pdf is undefined for x.eq.0
	ELSEIF (JLx .EQ. Nx-1) THEN
!  This is the highest x bin:
		CALL Polint (XVpow(Nx-3), Upd(J1), 4, ss, Fx, Dfx)
        Fvec(it) = Fx
	ELSE
!  for all interior points, use Jon's in-line function
!  This applied to (JLx.Ge.2 .and. JLx.Le.Nx-2)
		sf2 = Upd(J1+1)
		sf3 = Upd(J1+2)

		g1 =  sf2*const1 - sf3*const2
		g4 = -sf2*const3 + sf3*const4

		Fvec(it) = (const5*(Upd(J1)-g1) &
                 + const6*(Upd(J1+3)-g4) &
                 + sf2*sy3 - sf3*sy2) / s23

	ENDIF

ENDDO
! We now have the four values Fvec(1:4)
! interpolate in t...

IF(JLq .LE. 0) THEN
! 1st Q-bin, as well as extrapolation to lower Q
	CALL Polint (TV(0), Fvec(1), 4, tt, ff, Dfq)
ELSEIF (JLq .GE. Nt-1) THEN
! Last Q-bin, as well as extrapolation to higher Q
    CALL Polint (TV(Nt-3), Fvec(1), 4, tt, ff, Dfq)
ELSE
! Interrior bins : (JLq.GE.1 .and. JLq.LE.Nt-2)
! which include JLq.Eq.1 and JLq.Eq.Nt-2, since Upd is defined for
! the full range QV(0:Nt)  (in contrast to XV)
	tf2 = fvec(2)
	tf3 = fvec(3)

	g1 = ( tf2*t13 - tf3*t12) / t23
	g4 = (-tf2*t34 + tf3*t24) / t23

	h00 = ((t34*ty2-tmp2*ty3)*(fvec(1)-g1)/t12 &
         +  (tmp1*ty2-t12*ty3)*(fvec(4)-g4)/t34)

	ff = (h00*ty2*ty3/tdet + tf2*ty3 - tf3*ty2) / t23
ENDIF

PartonX6 = ff

!                                       ********************
END FUNCTION PartonX6

END MODULE CTEQ6PDF
