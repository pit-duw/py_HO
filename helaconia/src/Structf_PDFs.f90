MODULE Structf_PDFs
USE Helac_Global
USE CTEQ6PDF
USE setscale
IMPLICIT NONE
CONTAINS      
SUBROUTINE strf_pdf(wsf)
!include 'declare.h'
!include 'common_int.h'
!include 'common_strf.h'
!include 'parameter_pp.h'
INTEGER::ipp=1  ! ipp=1 ppbar ; ipp=2 pp
INTEGER::ih=1   ! ih=1 no photon PDF,ih=2, photon from proton/anti-proton,ih=3 photon from electron/positron 
REAL(KIND=DBL),INTENT(OUT)::wsf          
CHARACTER(len=20),DIMENSION(20)::parm
REAL(KIND=DBL),DIMENSION(20)::val
REAL(KIND=DBL),DIMENSION(-7:7)::pdflist
INTEGER::init=0
LOGICAL::use_cteq6_f90=.TRUE.      
SAVE init,ipp,ih,use_cteq6_f90
REAL(KIND=DBL)::glu1_ct,glu2_ct,u1_ct,u2_ct,d1_ct,d2_ct,s1_ct,s2_ct,c1_ct,c2_ct,b1_ct,b2_ct,&
                ub1_ct,ub2_ct,db1_ct,db2_ct,sb1_ct,sb2_ct,cb1_ct,cb2_ct,bb1_ct,bb2_ct,sf_ct
INCLUDE "../input/lhapdf/call_strf_lhapdf"
IF(init.EQ.0)THEN
   SELECT CASE(iPDFSUP1)
   CASE(10000)
      CALL SetCtq6f90(1)
      pdlabel="cteq6_m"
      nloop=2
      alphaQCD2=0.118d0
      use_cteq6_f90=.TRUE.
   CASE(10041)
      CALL SetCtq6f90(3)
      pdlabel="cteq6_l"
      nloop=2
      alphaQCD2=0.118d0
      use_cteq6_f90=.TRUE.
   CASE(10042)
      CALL SetCtq6f90(4)
      pdlabel="cteq6l1"
      nloop=1
      alphaQCD2=0.130d0
      use_cteq6_f90=.TRUE.
   CASE DEFAULT
      CALL pdfset_internal
      use_cteq6_f90=.FALSE.
   END SELECT
   IF(Coll_Type.EQ.1)THEN
      ipp=2
   ELSE
      ipp=1
   ENDIF
   init=1
ENDIF

CALL qcdscale(scale)
IF(reweight_scale)THEN
   IF(reweight_Fscale_phase.EQ.1)THEN
      scale=scale*rw_Fscale_up
   ELSEIF(reweight_Fscale_phase.EQ.-1)THEN
      scale=scale*rw_Fscale_down
   ENDIF
ENDIF   

IF(use_cteq6_f90)THEN
   glu1_ct = xp1*Ctq6Pdf_f90(0,xp1,scale)
   glu2_ct = xp2*Ctq6Pdf_f90(0,xp2,scale)
   u1_ct   = xp1*Ctq6Pdf_f90(1,xp1,scale)
   u2_ct   = xp2*Ctq6Pdf_f90(1,xp2,scale)
   d1_ct   = xp1*Ctq6Pdf_f90(2,xp1,scale)
   d2_ct   = xp2*Ctq6Pdf_f90(2,xp2,scale)
   s1_ct   = xp1*Ctq6Pdf_f90(3,xp1,scale)
   s2_ct   = xp2*Ctq6Pdf_f90(3,xp2,scale)
   c1_ct   = xp1*Ctq6Pdf_f90(4,xp1,scale)
   c2_ct   = xp2*Ctq6Pdf_f90(4,xp2,scale)
   b1_ct   = xp1*Ctq6Pdf_f90(5,xp1,scale)
   b2_ct   = xp2*Ctq6Pdf_f90(5,xp2,scale)
   ub1_ct  = xp1*Ctq6Pdf_f90(-1,xp1,scale)
   ub2_ct  = xp2*Ctq6Pdf_f90(-1,xp2,scale)
   db1_ct  = xp1*Ctq6Pdf_f90(-2,xp1,scale)
   db2_ct  = xp2*Ctq6Pdf_f90(-2,xp2,scale)
   sb1_ct  = xp1*Ctq6Pdf_f90(-3,xp1,scale)
   sb2_ct  = xp2*Ctq6Pdf_f90(-3,xp2,scale)
   cb1_ct  = xp1*Ctq6Pdf_f90(-4,xp1,scale)
   cb2_ct  = xp2*Ctq6Pdf_f90(-4,xp2,scale)
   bb1_ct  = xp1*Ctq6Pdf_f90(-5,xp1,scale)
   bb2_ct  = xp2*Ctq6Pdf_f90(-5,xp2,scale)
ELSE
   CALL fdist(ih,xp1,scale,pdflist(-7:7))
   glu1_ct = xp1*pdflist(0)
   u1_ct   = xp1*pdflist(2)
   d1_ct   = xp1*pdflist(1)
   s1_ct   = xp1*pdflist(3)
   c1_ct   = xp1*pdflist(4)
   b1_ct   = xp1*pdflist(5)
   ub1_ct  = xp1*pdflist(-2)
   db1_ct  = xp1*pdflist(-1)
   sb1_ct  = xp1*pdflist(-3)
   cb1_ct  = xp1*pdflist(-4)
   bb1_ct  = xp1*pdflist(-5)
   CALL fdist(ih,xp2,scale,pdflist(-7:7))
   glu2_ct = xp2*pdflist(0)
   u2_ct   = xp2*pdflist(2)
   d2_ct   = xp2*pdflist(1)
   s2_ct   = xp2*pdflist(3)
   c2_ct   = xp2*pdflist(4)
   b2_ct   = xp2*pdflist(5)
   ub2_ct  = xp2*pdflist(-2)
   db2_ct  = xp2*pdflist(-1)
   sb2_ct  = xp2*pdflist(-3)
   cb2_ct  = xp2*pdflist(-4)
   bb2_ct  = xp2*pdflist(-5)
ENDIF
! recover it
IF(reweight_scale)THEN
   IF(reweight_Fscale_phase.EQ.1)THEN
      scale=scale/rw_Fscale_up
   ELSEIF(reweight_Fscale_phase.EQ.-1)THEN
      scale=scale/rw_Fscale_down
   ENDIF
ENDIF
sf_ct=1
IF(qsumQ)THEN
	IF(ifl(1).EQ.35)THEN
		sf_ct=sf_ct*glu1_ct
	ELSEIF(ABS(ifl(1)).LE.12.AND.(MOD(ABS(ifl(1)),4).EQ.0.OR.MOD(ABS(ifl(1)),4).EQ.3) &
	.AND.ifl(2).EQ.35)THEN
		IF(onlyqcd)THEN
			SELECT CASE(iqnum)
			CASE(1)
				sf_ct=sf_ct*(u1_ct+ub1_ct)
			CASE(2)
				sf_ct=sf_ct*(u1_ct+ub1_ct+d1_ct+db1_ct)
			CASE(3)
				sf_ct=sf_ct*(u1_ct+ub1_ct+d1_ct+db1_ct+s1_ct+sb1_ct)
			CASE(4)
				sf_ct=sf_ct*(c1_ct+cb1_ct+u1_ct+ub1_ct+d1_ct+db1_ct+s1_ct+sb1_ct)
			CASE DEFAULT
				sf_ct=sf_ct*(b1_ct+bb1_ct+c1_ct+cb1_ct+u1_ct+ub1_ct+d1_ct+db1_ct+s1_ct+sb1_ct)
			END SELECT
		ELSE
			IF(MOD(ABS(ifl(1)),4).EQ.3)THEN
				SELECT CASE(iqnum)
				CASE(1,2,3)
					sf_ct=sf_ct*(u1_ct+ub1_ct)
				CASE DEFAULT
					sf_ct=sf_ct*(u1_ct+ub1_ct+c1_ct+cb1_ct)
				END SELECT
			ELSEIF(MOD(ABS(ifl(1)),4).EQ.0)THEN
				SELECT CASE(iqnum)
				CASE(1)
					sf_ct=sf_ct
				CASE(2)
					sf_ct=sf_ct*(d1_ct+db1_ct)
				CASE(3,4)
					sf_ct=sf_ct*(d1_ct+db1_ct+s1_ct+sb1_ct)
				CASE DEFAULT
					sf_ct=sf_ct*(d1_ct+db1_ct+s1_ct+sb1_ct+b1_ct+bb1_ct)
				END SELECT
			ENDIF
		ENDIF
	ENDIF

	IF(ifl(2).EQ.35)THEN
		sf_ct=sf_ct*glu2_ct
	ELSEIF(ABS(ifl(2)).LE.12.AND.(MOD(ABS(ifl(2)),4).EQ.0.OR.MOD(ABS(ifl(2)),4).EQ.3) &
	.AND.ifl(1).EQ.35)THEN
		IF(onlyqcd)THEN
			SELECT CASE(iqnum)
			CASE(1)
				sf_ct=sf_ct*(u2_ct+ub2_ct)
			CASE(2)
				sf_ct=sf_ct*(u2_ct+ub2_ct+d2_ct+db2_ct)
			CASE(3)
				sf_ct=sf_ct*(u2_ct+ub2_ct+d2_ct+db2_ct+s2_ct+sb2_ct)
			CASE(4)
				sf_ct=sf_ct*(c2_ct+cb2_ct+u2_ct+ub2_ct+d2_ct+db2_ct+s2_ct+sb2_ct)
			CASE DEFAULT
				sf_ct=sf_ct*(b2_ct+bb2_ct+c2_ct+cb2_ct+u2_ct+ub2_ct+d2_ct+db2_ct+s2_ct+sb2_ct)
			END SELECT
		ELSE
			IF(MOD(ABS(ifl(2)),4).EQ.3)THEN
				SELECT CASE(iqnum)
				CASE(1,2,3)
					sf_ct=sf_ct*(u2_ct+ub2_ct)
				CASE DEFAULT
					sf_ct=sf_ct*(u2_ct+ub2_ct+c2_ct+cb2_ct)
				END SELECT
			ELSEIF(MOD(ABS(ifl(2)),4).EQ.0)THEN
				SELECT CASE(iqnum)
				CASE(1)
					sf_ct=sf_ct
				CASE(2)
					sf_ct=sf_ct*(d2_ct+db2_ct)
				CASE(3,4)
					sf_ct=sf_ct*(d2_ct+db2_ct+s2_ct+sb2_ct)
				CASE DEFAULT
					sf_ct=sf_ct*(d2_ct+db2_ct+s2_ct+sb2_ct+b2_ct+bb2_ct)
				END SELECT
			ENDIF
		ENDIF
	ENDIF

	IF(ifl(1)+ifl(2).EQ.0.AND.(MOD(ABS(ifl(1)),4).EQ.0.OR.MOD(ABS(ifl(1)),4).EQ.3)&
      .AND.ABS(ifl(1)).LE.12)THEN
                IF(sf_ct.NE.1d0)THEN
                   WRITE(*,*)"ERROR 1 in Structf_PDFs.f90"
                   STOP
                ENDIF
		IF(onlyqcd)THEN
			SELECT CASE(iqnum)
			CASE(1)
				sf_ct=sf_ct*(u1_ct*ub2_ct+ub1_ct*u2_ct)
			CASE(2)
				sf_ct=sf_ct*(u1_ct*ub2_ct+ub1_ct*u2_ct+d1_ct*db2_ct+d2_ct*db1_ct)
			CASE(3)
				sf_ct=sf_ct*(u1_ct*ub2_ct+ub1_ct*u2_ct+d1_ct*db2_ct+d2_ct*db1_ct&
				+s1_ct*sb2_ct+s2_ct*sb1_ct)
			CASE(4)
				sf_ct=sf_ct*(u1_ct*ub2_ct+ub1_ct*u2_ct+d1_ct*db2_ct+d2_ct*db1_ct&
				+s1_ct*sb2_ct+s2_ct*sb1_ct+c1_ct*cb2_ct+cb1_ct*c2_ct)
			CASE DEFAULT
				sf_ct=sf_ct*(u1_ct*ub2_ct+ub1_ct*u2_ct+d1_ct*db2_ct+d2_ct*db1_ct&
				+s1_ct*sb2_ct+s2_ct*sb1_ct+c1_ct*cb2_ct+cb1_ct*c2_ct+b1_ct*bb2_ct+bb1_ct*b2_ct)
			END SELECT
		ELSE
			IF(MOD(ABS(ifl(1)),4).EQ.3)THEN
				SELECT CASE(iqnum)
				CASE(1,2,3)
					sf_ct=sf_ct*(u1_ct*ub2_ct+ub1_ct*u2_ct)
				CASE DEFAULT
					sf_ct=sf_ct*(u1_ct*ub2_ct+ub1_ct*u2_ct+c1_ct*cb2_ct+c2_ct*cb1_ct)
				END SELECT
			ELSEIF(MOD(ABS(ifl(1)),4).EQ.0)THEN
				SELECT CASE(iqnum)
				CASE(1)
					sf_ct=sf_ct
				CASE(2)
					sf_ct=sf_ct*(d1_ct*db2_ct+d2_ct*db1_ct)
				CASE(3,4)
					sf_ct=sf_ct*(s1_ct*sb2_ct+sb1_ct*s2_ct+d1_ct*db2_ct+d2_ct*db1_ct)
				CASE DEFAULT
					sf_ct=sf_ct*(s1_ct*sb2_ct+sb1_ct*s2_ct+d1_ct*db2_ct+d2_ct*db1_ct&
					+b1_ct*bb2_ct+b2_ct*bb1_ct)
				END SELECT
			ENDIF
		ENDIF
        ELSEIF(ifl(1)-ifl(2).EQ.0.AND.(MOD(ABS(ifl(1)),4).EQ.0.OR.MOD(ABS(ifl(1)),4).EQ.3)&
             .AND.onlyqcd.AND.ABS(ifl(1)).LE.12)THEN
                IF(sf_ct.NE.1d0)THEN
                   WRITE(*,*)"ERROR 2 in Structf_PDFs.f90"
                   STOP
                ENDIF
           ! q q > q q in QCD
           SELECT CASE(iqnum)
           CASE(1)
              sf_ct=sf_ct*(u1_ct*u2_ct+ub1_ct*ub2_ct) ! u u + u~ u~
           CASE(2)
              sf_ct=sf_ct*(u1_ct*u2_ct+ub1_ct*ub2_ct+d1_ct*d2_ct+db2_ct*db1_ct) ! u u, d d,u~ u~, d~ d~
           CASE(3)
              sf_ct=sf_ct*(u1_ct*u2_ct+ub1_ct*ub2_ct+d1_ct*d2_ct+db2_ct*db1_ct&
                   +s1_ct*s2_ct+sb2_ct*sb1_ct) ! u u, d d, s s, u~ u~, d~ d~, s~ s~
           CASE(4)
              sf_ct=sf_ct*(u1_ct*u2_ct+ub1_ct*ub2_ct+d1_ct*d2_ct+db2_ct*db1_ct&
                   +s1_ct*s2_ct+sb2_ct*sb1_ct+c1_ct*c2_ct+cb1_ct*cb2_ct)
           CASE DEFAULT
              sf_ct=sf_ct*(u1_ct*u2_ct+ub1_ct*ub2_ct+d1_ct*d2_ct+db2_ct*db1_ct&
                   +s1_ct*s2_ct+sb2_ct*sb1_ct+c1_ct*c2_ct+cb1_ct*cb2_ct+b1_ct*b2_ct+bb1_ct*bb2_ct)
           END SELECT
        ELSEIF((MOD(ABS(ifl(1)),4).EQ.0.OR.MOD(ABS(ifl(1)),4).EQ.3).AND.onlyqcd.AND.iqnum.GT.1&
             .AND.ABS(ifl(1)).LE.12.AND.ABS(ifl(2)).LE.12&
             .AND.(MOD(ABS(ifl(2)),4).EQ.0.OR.MOD(ABS(ifl(2)),4).EQ.3))THEN
                IF(sf_ct.NE.1d0)THEN
                   WRITE(*,*)"ERROR 3 in Structf_PDFs.f90"
                   STOP
                ENDIF
           ! q q' > q q' in QCD
           SELECT CASE(iqnum)
           CASE(2)
              sf_ct=sf_ct*((u1_ct+ub1_ct)*(d2_ct+db2_ct)+(d1_ct+db1_ct)*(u2_ct+ub2_ct))
           CASE(3)
              sf_ct=sf_ct*((u1_ct+ub1_ct)*(d2_ct+db2_ct+s2_ct+sb2_ct)&
                   +(d1_ct+db1_ct)*(u2_ct+ub2_ct+s2_ct+sb2_ct)&
                   +(s1_ct+sb1_ct)*(u2_ct+ub2_ct+d2_ct+db2_ct))
           CASE(4)
              sf_ct=sf_ct*((u1_ct+ub1_ct)*(d2_ct+db2_ct+s2_ct+sb2_ct+c2_ct+cb2_ct)&
                   +(d1_ct+db1_ct)*(u2_ct+ub2_ct+s2_ct+sb2_ct+c2_ct+cb2_ct)&
                   +(s1_ct+sb1_ct)*(u2_ct+ub2_ct+d2_ct+db2_ct+c2_ct+cb2_ct)&
                   +(c1_ct+cb1_ct)*(u2_ct+ub2_ct+d2_ct+db2_ct+s2_ct+sb2_ct))
           CASE DEFAULT
              sf_ct=sf_ct*((u1_ct+ub1_ct)*(d2_ct+db2_ct+s2_ct+sb2_ct+c2_ct+cb2_ct+b2_ct+bb2_ct)&
                   +(d1_ct+db1_ct)*(u2_ct+ub2_ct+s2_ct+sb2_ct+c2_ct+cb2_ct+b2_ct+bb2_ct)&
                   +(s1_ct+sb1_ct)*(u2_ct+ub2_ct+d2_ct+db2_ct+c2_ct+cb2_ct+b2_ct+bb2_ct)&
                   +(c1_ct+cb1_ct)*(u2_ct+ub2_ct+d2_ct+db2_ct+s2_ct+sb2_ct+b2_ct+bb2_ct)&
                   +(b1_ct+bb1_ct)*(u2_ct+ub2_ct+d2_ct+db2_ct+s2_ct+sb2_ct+c2_ct+cb2_ct))
           END SELECT
	ELSEIF(ifl(1).NE.35.AND.ifl(2).NE.35)THEN
		PRINT *,"Please choose quarksumQ=.FALSE."
		STOP
	ENDIF

ELSE      
	IF(ipp.EQ.1)THEN
	! ppbar
		SELECT CASE(ifl(1))  
		CASE(3)
			sf_ct=sf_ct*u1_ct
		CASE(-3) 
			sf_ct=sf_ct*ub1_ct
		CASE(4)  
			sf_ct=sf_ct*d1_ct
		CASE(-4) 
			sf_ct=sf_ct*db1_ct
		CASE(7)  
			sf_ct=sf_ct*c1_ct
		CASE(-7) 
			sf_ct=sf_ct*cb1_ct
		CASE(8)  
			sf_ct=sf_ct*s1_ct
		CASE(-8) 
			sf_ct=sf_ct*sb1_ct
!   IF(ifl(1).EQ.11) sf_ct=sf_ct*t1_ct
!   IF(ifl(1).EQ.-11)sf_ct=sf_ct*tb1_ct
		CASE(12) 
			sf_ct=sf_ct*b1_ct
		CASE(-12)
			sf_ct=sf_ct*bb1_ct
		CASE(35) 
			sf_ct=sf_ct*glu1_ct
		END SELECT

		SELECT CASE(ifl(2))
		CASE(3)
			sf_ct=sf_ct*ub2_ct
		CASE(-3) 
			sf_ct=sf_ct*u2_ct
		CASE(4)  
			sf_ct=sf_ct*db2_ct
		CASE(-4) 
			sf_ct=sf_ct*d2_ct
		CASE(7)  
			sf_ct=sf_ct*cb2_ct
		CASE(-7) 
			sf_ct=sf_ct*c2_ct
		CASE(8)  
			sf_ct=sf_ct*sb2_ct
		CASE(-8) 
			sf_ct=sf_ct*s2_ct
!       if(ifl(2).eq.11) sf=sf*t2
!       if(ifl(2).eq.-11)sf=sf*tb2
		CASE(12) 
			sf_ct=sf_ct*bb2_ct
		CASE(-12)
			sf_ct=sf_ct*b2_ct
		CASE(35) 
			sf_ct=sf_ct*glu2_ct
		END SELECT
	ELSEIF(ipp.EQ.2)THEN
! pp
		SELECT CASE(ifl(1))
		CASE(3)  
			sf_ct=sf_ct*u1_ct
		CASE(-3) 
			sf_ct=sf_ct*ub1_ct
		CASE(4)  
			sf_ct=sf_ct*d1_ct
		CASE(-4) 
			sf_ct=sf_ct*db1_ct
		CASE(7)  
			sf_ct=sf_ct*c1_ct
		CASE(-7) 
			sf_ct=sf_ct*cb1_ct
		CASE(8)  
			sf_ct=sf_ct*s1_ct
		CASE(-8) 
			sf_ct=sf_ct*sb1_ct
!       if(ifl(1).eq.11) sf=sf*t1
!       if(ifl(1).eq.-11)sf=sf*tb1
		CASE(12) 
			sf_ct=sf_ct*b1_ct
		CASE(-12)
			sf_ct=sf_ct*bb1_ct
		CASE(35) 
			sf_ct=sf_ct*glu1_ct
		END SELECT

		SELECT CASE(ifl(2))
		CASE(3)  
			sf_ct=sf_ct*u2_ct
		CASE(-3) 
			sf_ct=sf_ct*ub2_ct
		CASE(4)  
			sf_ct=sf_ct*d2_ct
		CASE(-4) 
			sf_ct=sf_ct*db2_ct
		CASE(7)  
			sf_ct=sf_ct*c2_ct
		CASE(-7) 
			sf_ct=sf_ct*cb2_ct
		CASE(8)  
			sf_ct=sf_ct*s2_ct
		CASE(-8) 
			sf_ct=sf_ct*sb2_ct
!       if(ifl(2).eq.11) sf=sf*t2
!       if(ifl(2).eq.-11)sf=sf*tb2
		CASE(12) 
			sf_ct=sf_ct*b2_ct
		CASE(-12)
			sf_ct=sf_ct*bb2_ct
		CASE(35) 
			sf_ct=sf_ct*glu2_ct
		END SELECT
	ENDIF
ENDIF

!      if(ifl(1).eq. 3.and.ifl(2).eq.-3.or.
!    &    ifl(1).eq.-3.and.ifl(2).eq. 3)then
!      if(init.eq.0)print*,'qq',ifl(1),ifl(2)
!      sf =uv1*us2+us1*uv2+dv1*ds2+ds1*dv2
!    &        +2.d0*(us1*us2+ds1*ds2)            !q+q or qb+qb
!    &        +2.d0*(st1*st2+ch1*ch2+bo1*bo2)
!      endif
!      if(ifl(1).eq.35.and.ifl(2).eq.35)then 
!      if(init.eq.0)print*,'gg',ifl(1),ifl(2)
!      sf =gl1*gl2                           !g+g
!      endif                             

wsf=wjac*sf_ct/xp1/xp2
IF(init.EQ.0)init=1
       
END SUBROUTINE strf_pdf

!       include 'pdfset.h'
SUBROUTINE pdfset_internal
!       include 'declare.h'
! choose the pdfset
  IMPLICIT NONE
  SELECT CASE(iPDFSUP1)
  CASE(20200)
     pdlabel="mrs02nl"
  CASE(20270)
     pdlabel="mrs02nn"
  CASE(20000)
     pdlabel="mrs0119"
  CASE(20002)
     pdlabel="mrs0117"
  CASE(20003)
     pdlabel="mrs0121"
  CASE(20004)
     pdlabel="mrs01_j"
  CASE(28000)
     pdlabel="mrs99_1"
  CASE(28002)
     pdlabel="mrs99_2"
  CASE(28001)
     pdlabel="mrs99_3"
  CASE(28003)
     pdlabel="mrs99_4"
  CASE(28004)
     pdlabel="mrs99_5"
  CASE(28005)
     pdlabel="mrs99_6"
  CASE(28006)
     pdlabel="mrs99_7"
  CASE(28007)
     pdlabel="mrs99_8"
  CASE(28008)
     pdlabel="mrs99_9"
  CASE(28009)
     pdlabel="mrs9910"
  CASE(28010)
     pdlabel="mrs9911"
  CASE(28011)
     pdlabel="mrs9912"
  CASE(29000)
     pdlabel="mrs98z1"
  CASE(29002)
     pdlabel="mrs98z2"
  CASE(29001)
     pdlabel="mrs98z3"
  CASE(29003)
     pdlabel="mrs98z4"
  CASE(29004)
     pdlabel="mrs98z5"
  CASE(29071)
     pdlabel="mrs98ht"
  CASE(29041)
     pdlabel="mrs98l1"
  CASE(29043)
     pdlabel="mrs98l2"
  CASE(29042)
     pdlabel="mrs98l3"
  CASE(29044)
     pdlabel="mrs98l4"
  CASE(29045)
     pdlabel="mrs98l5"
  CASE(10000)
     pdlabel="cteq6_m"
  CASE(10043)
     pdlabel="cteq6_d"
  CASE(10041)
     pdlabel="cteq6_l"
  CASE(10042)
     pdlabel="cteq6l1"
  CASE(19050)
     pdlabel="cteq5_m"
  CASE(19060)
     pdlabel="cteq5_d"
  CASE(19070)
     pdlabel="cteq5_l"
  CASE(19071)
     pdlabel="cteq5l1"
  CASE(19072)
     pdlabel="cteq5hj"
  CASE(19073)
     pdlabel="cteq5hq"
  CASE(19053)
     pdlabel="cteq5f3"
  CASE(19054)
     pdlabel="cteq5f4"
  CASE(19051)
     pdlabel="cteq5m1"
  CASE(19075)
     pdlabel="ctq5hq1"
  CASE(19150)
     pdlabel="cteq4_m"
  CASE(19160)
     pdlabel="cteq4_d"
  CASE(19170)
     pdlabel="cteq4_l"
  CASE(19151)
     pdlabel="cteq4a1"
  CASE(19152)
     pdlabel="cteq4a2"
  CASE(19153)
     pdlabel="cteq4a3"
  CASE(19154)
     pdlabel="cteq4a4"
  CASE(19155)
     pdlabel="cteq4a5"
  CASE(19156)
     pdlabel="cteq4hj"
  CASE(19157)
     pdlabel="cteq4lq"
  CASE(19250)
     pdlabel="cteq3_m"
  CASE(19270)
     pdlabel="cteq3_l"
  CASE(19260)
     pdlabel="cteq3_d"
  END SELECT
  CALL pdfwrap
  RETURN
END SUBROUTINE pdfset_internal
!      subroutine structm
!SUBROUTINE structm(xp1,scale,uv1,dv1,us1,ds1,st1,ch1,bo1,to1,gl1) 
!       include 'declare.h'
!PRINT *,'YOU SHOULD TAKE CARE OF PDF SET LIBRARY'
!STOP
!END SUBROUTINE structm

FUNCTION alphas2(xxxx)
!       include 'declare.h'
REAL(KIND=DBL)::ALPHAS2
REAL(KIND=DBL),INTENT(IN)::xxxx
alphas2=0.130d0
END FUNCTION alphas2

SUBROUTINE pdfwrap
  IMPLICIT NONE
  nloop=2 ! NLO running unless set otherwise

!
!  MRST2002
!  1     NLO   0.1197    0.00949 
!  2     NNLO  0.1154    0.00685
!  
  SELECT CASE(pdlabel)
  CASE('mrs02nl')
     ! 'MRST2002nlo (Standard MSbar)'
     alphaQCD2=0.1197d0
  CASE('mrs02nn')
     ! 'MRST2002nnlo (NNLO fit)'
     alphaQCD2=0.1154d0
!
!  MRST2001
!  1     alf119  central gluon, a_s       323      0.119    0.00927  
!  2     alf117  lower a_s                290      0.117    0.00953  
!  3     alf121  higher a_s               362      0.121    0.00889  
!  4     j121    better fit to jet data   353      0.121    0.00826  
!
  CASE('mrs0119')
     ! 'MRST2001nlo (Standard MSbar)'
     alphaQCD2=0.119d0
  CASE('mrs0117')
     ! 'MRST2001nlo (lower alpha_S)'
     alphaQCD2=0.117d0
  CASE('mrs0121')
     ! 'MRST2001nlo (higher alpha_S)'
     alphaQCD2=0.121d0
  CASE('mrs01_j')
     ! 'MRST2001nlo (Jet Fit)'
     alphaQCD2=0.121d0
!
! MRS99NLO
!  1     COR01  central gluon, a_s    300      0.1175   0.00537  C
!  2     COR02  higher gluon          300      0.1175   0.00497  C
!  3     COR03  lower gluon           300      0.1175   0.00398  C
!  4     COR04  lower a_s             229      0.1125   0.00585  C
!  5     COR05  higher a_s            383      0.1225   0.00384  C
!  6     COR06  quarks up             303.3    0.1178   0.00497  C
!  7     COR07  quarks down           290.3    0.1171   0.00593  C
!  8     COR08  strange up            300      0.1175   0.00524  C
!  9     COR09  strange down          300      0.1175   0.00524  C
!  10    C0R10  charm up              300      0.1175   0.00525  C
!  11    COR11  charm down            300      0.1175   0.00524  C
!  12    COR12  larger d/u            300      0.1175   0.00515  C
!
  CASE('mrs99_1')
     alphaQCD2=0.1175d0
  CASE('mrs99_2')
     alphaQCD2=0.1175d0
  CASE('mrs99_3')
     alphaQCD2=0.1175d0
  CASE('mrs99_4')
     alphaQCD2=0.1125d0
  CASE('mrs99_5')
     alphaQCD2=0.1225d0
  CASE('mrs99_6')
     alphaQCD2=0.1178d0
  CASE('mrs99_7')
     alphaQCD2=0.1171d0
  CASE('mrs99_8')
     alphaQCD2=0.1175d0
  CASE('mrs99_9')
     alphaQCD2=0.1175d0
  CASE('mrs9910')
     alphaQCD2=0.1175d0
  CASE('mrs9911')
     alphaQCD2=0.1175d0
  CASE('mrs9912')
     alphaQCD2=0.1175d0
!
! MRS98NLO
!    ft08a  central gluon, a_s  300      0.1175   0.00561  
!    ft09a  higher gluon        300      0.1175   0.00510  
!    ft11a  lower gluon         300      0.1175   0.00408  
!    ft24a  lower a_s           229      0.1125   0.00586  
!    ft23a  higher a_s          383      0.1225   0.00410  
!
  CASE('mrs98z1')
     ! 'MRST98nlo (central gluon/alphas)'
     alphaQCD2=0.1175d0
  CASE('mrs98z2')
     ! 'MRST98nlo (higher gluon)'
     alphaQCD2=0.1175d0
  CASE('mrs98z3')
     ! 'MRST98nlo (lower gluon)'
     alphaQCD2=0.1175d0
  CASE('mrs98z4')
     ! 'MRST98nlo (lower alpha_S)'
     alphaQCD2=0.1125d0
  CASE('mrs98z5')
     ! 'MRST98nlo (higher alpha_S)'
     alphaQCD2=0.1225d0
  CASE('mrs98ht')
     ! 'MRST98nloht (HT)'
!-- real value
     alphaQCD2=0.1170d0
!-- modified - DEBUG
     alphaQCD2=0.1175d0
     WRITE(*,*) 'alpha_s(MZ) for mrs98ht has been modified from'
     WRITE(*,*) 'the inherent 0.1170 to a new value of 0.1175'    
!
!  MRS98LO
!    lo05a  central gluon, a_s  174      0.1250   0.01518  
!    lo09a  higher gluon        174      0.1250   0.01616  
!    lo10a  lower gluon         174      0.1250   0.01533  
!    lo01a  lower a_s           136      0.1200   0.01652  
!    lo07a  higher a_s          216      0.1300   0.01522  
!
  CASE('mrs98l1')
     ! 'MRST98lo (central gluon/alphas)'
     alphaQCD2=0.125d0
     nloop=1
  CASE('mrs98l2')
     ! 'MRST98lo (higher gluon)'
     alphaQCD2=0.125d0
     nloop=1
  CASE('mrs98l3')
     ! 'MRST98lo (lower gluon)'
     alphaQCD2=0.125d0
     nloop=1
  CASE('mrs98l4')
     ! 'MRST98lo (lower a_s)'
     alphaQCD2=0.120d0
     nloop=1
  CASE('mrs98l5')
     ! 'MRST98lo (higher a_s)'
     alphaQCD2=0.130d0
     nloop=1
!
! CTEQ4
!   1      CTEQ4M   Standard MSbar scheme   0.116        1.6      cteq4m.tbl
!   2      CTEQ4D   Standard DIS scheme     0.116        1.6      cteq4d.tbl
!   3      CTEQ4L   Leading Order           0.116        1.6      cteq4l.tbl
!   4      CTEQ4A1  Alpha_s series          0.110        1.6      cteq4a1.tbl
!   5      CTEQ4A2  Alpha_s series          0.113        1.6      cteq4a2.tbl
!   6      CTEQ4A3  same as CTEQ4M          0.116        1.6      cteq4m.tbl
!   7      CTEQ4A4  Alpha_s series          0.119        1.6      cteq4a4.tbl
!   8      CTEQ4A5  Alpha_s series          0.122        1.6      cteq4a5.tbl
!   9      CTEQ4HJ  High Jet                0.116        1.6      cteq4hj.tbl
!   10     CTEQ4LQ  Low Q0                  0.114        0.7      cteq4lq.tbl
!
  CASE('cteq3_m')
     ! 'CTEQ3m (Standard MSbar)'
     alphaQCD2=0.112d0
  CASE('cteq3_l')
     ! 'CTEQ3l (LO fit)'
     alphaQCD2=0.112d0
     nloop=1
  CASE('cteq3_d')
     ! 'CTEQ3d (Standard DIS)'
     alphaQCD2=0.112d0
  CASE('cteq4_m')
     ! 'CTEQ4m (Standard MSbar)'
     alphaQCD2=0.116d0
  CASE('cteq4_d')
     ! 'CTEQ4d (Standard DIS)'
     alphaQCD2=0.116d0
  CASE('cteq4_l')
     ! 'CTEQ4l (LO fit)'
     alphaQCD2=0.132d0
     nloop=1
  CASE('cteq4a1')
     ! 'CTEQ4A1'
     alphaQCD2=0.110d0
  CASE('cteq4a2')
     ! 'CTEQ4A2'
     alphaQCD2=0.113d0
  CASE('cteq4a3')
     ! 'CTEQ4A3'
     alphaQCD2=0.116d0
  CASE('cteq4a4')
     ! 'CTEQ4A4'
     alphaQCD2=0.119d0
  CASE('cteq4a5')
     ! 'CTEQ4A5'
     alphaQCD2=0.122d0
  CASE('cteq4hj')
     ! 'CTEQ4HJ (high jet)'
     alphaQCD2=0.116d0
  CASE('cteq4lq')
     ! 'CTEQ4LQ (low Q0)'
     alphaQCD2=0.114d0
!
! ---------------------------------------------------------------------------
!  Iset   PDF        Description       Alpha_s(Mz)  Lam4  Lam5   Table_File
! ---------------------------------------------------------------------------
!   1    CTEQ5M   Standard MSbar scheme   0.118     326   226    cteq5m.tbl
!   2    CTEQ5D   Standard DIS scheme     0.118     326   226    cteq5d.tbl
!   3    CTEQ5L   Leading Order           0.127     192   146    cteq5l.tbl
!   4    CTEQ5HJ  Large-x gluon enhanced  0.118     326   226    cteq5hj.tbl
!   5    CTEQ5HQ  Heavy Quark             0.118     326   226    cteq5hq.tbl
!   6    CTEQ5F3  Nf=3 FixedFlavorNumber  0.106     (Lam3=395)   cteq5f3.tbl
!   7    CTEQ5F4  Nf=4 FixedFlavorNumber  0.112     309   XXX    cteq5f4.tbl
!         --------------------------------------------------------
!   8    CTEQ5M1  Improved CTEQ5M         0.118     326   226    cteq5m1.tbl
!   9    CTEQ5HQ1 Improved CTEQ5HQ        0.118     326   226    ctq5hq1.tbl
! ---------------------------------------------------------------------------
!
  CASE('cteq5_m')
     ! 'CTEQ5m (Standard MSbar)'
     Call SetCtq5(1)
     alphaQCD2=0.118d0
  CASE('cteq5_d')
     ! 'CTEQ5d (Standard DIS)'
     Call SetCtq5(2)
     alphaQCD2=0.118d0
  CASE('cteq5_l')
     ! 'CTEQ5l (LO fit)'
     Call SetCtq5(3)
     alphaQCD2=0.127d0
     nloop=1
  CASE('cteq5l1')
     ! 'CTEQ5l1 (LO fit)'
     alphaQCD2=0.127d0
     nloop=1
  CASE('cteq5hj')
     ! 'CTEQ5hj (large-x gluon enhanced)'
     Call SetCtq5(4)
     alphaQCD2=0.118d0
  CASE('cteq5hq')
     ! 'CTEQ5hq (heavy quark)'
     Call SetCtq5(5)
     alphaQCD2=0.118d0
  CASE('cteq5f3')
     ! 'CTEQ5f3 (3-flav-DIS)'
     Call SetCtq5(6)
     alphaQCD2=0.106d0
  CASE('cteq5f4')
     ! 'CTEQ5f4 (4-flav-DIS)'
     Call SetCtq5(7)
     alphaQCD2=0.112d0
  CASE('cteq5m1')
     ! 'CTEQ5m1 (updated CTEQ5m)'
     Call SetCtq5(8)
     alphaQCD2=0.118d0
  CASE('ctq5hq1')
     Call SetCtq5(9)
     alphaQCD2=0.118d0
!
!   1    CTEQ6M   Standard MSbar scheme   0.118     326   226    cteq6m.tbl
!   2    CTEQ6D   Standard DIS scheme     0.118     326   226    cteq6d.tbl
!   3    CTEQ6L   Leading Order           0.118**   326** 226    cteq6l.tbl
!   4    CTEQ6L1  Leading Order           0.130**   215** 165    cteq6l1.tbl
!
! Note:CTEQ6L1 uses the LO running alpha_s 
!
  CASE('cteq6_m')
     ! 'CTEQ6m (Standard MSbar)'
     alphaQCD2=0.118d0
     Call SetCtq6(1)
  CASE('cteq6_d')
     ! 'CTEQ6d (Standard DIS)'
     alphaQCD2=0.118d0
     Call SetCtq6(2)
  CASE('cteq6_l')
     ! 'CTEQ6l (LO fit/NLO alphas)'
     alphaQCD2=0.118d0
     Call SetCtq6(3)
  CASE('cteq6l1')
     ! 'CTEQ6l1 (LO fit/LO alphas)'
     alphaQCD2=0.130d0
     nloop=1
     Call SetCtq6(4)
  CASE DEFAULT
     WRITE(*,*) 'Unimplemented distribution= ',pdlabel
     WRITE(*,*) 'Implemented are: ',&
          'mrs02nl,','mrs02nn,',&
          'mrs0119,','mrs0117,','mrs0121,','mrs01_j,',&
          'mrs99_1,','mrs99_2,','mrs99_3,','mrs99_4,','mrs99_5,','mrs99_6,',&
          'mrs99_7,','mrs99_8,','mrs99_9,','mrs9910,','mrs9911,','mrs9912,',&
          'mrs98z1,','mrs98z2,','mrs98z3,','mrs98z4,','mrs98z5,','mrs98ht,',&
          'mrs98l1,','mrs98l2,','mrs98l3,','mrs98l4,','mrs98l5,',&
          'cteq3_m,','cteq3_l,','cteq3_d,',&
          'cteq4_m,','cteq4_d,','cteq4_l,','cteq4a1,','cteq4a2,',&
          'cteq4a3,','cteq4a4,','cteq4a5,','cteq4hj,','cteq4lq,',&
          'cteq5_m,','cteq5_d,','cteq5_l,','cteq5hj,','cteq5hq,',&
          'cteq5f3,','cteq5f4,','cteq5m1,','ctq5hq1,','cteq5l1,',&
          'cteq6_m,','cteq6_d,','cteq6_l,','cteq6l1,'
     STOP
  END SELECT
  RETURN
END SUBROUTINE pdfwrap

FUNCTION pdg2pdf(ih,ipdg,x,xmu)
!***************************************************************************
!     Based on pdf.f, wrapper for calling the pdf of MCFM
!***************************************************************************
  IMPLICIT NONE
  REAL(KIND(1d0))::pdg2pdf
  REAL(KIND(1d0)),INTENT(IN)::x,xmu
  INTEGER,INTENT(IN)::ih,ipdg
  REAL(KIND(1d0)),EXTERNAL::Ctq3df,Ctq4Fn,Ctq5Pdf,Ctq6Pdf,Ctq5L
  INTEGER::mode,Irt,i,j
  REAL(KIND(1d0)),DIMENSION(2)::xlast,xmulast
  REAL(KIND(1d0)),DIMENSION(-7:7,2)::pdflast
  REAL(KIND(1d0))::q2max
  CHARACTER(len=7),DIMENSION(2)::pdlabellast
!  REAL(KIND(1d0))::epa_electron,epa_proton
  INTEGER::ipart,ireuse,iporg
  INTEGER,DIMENSION(2)::ihlast
  LOGICAL::init=.TRUE.
  SAVE xlast,xmulast,pdflast,pdlabellast,ihlast,init

  IF(init)THEN
     xlast(1:2)=-99d9
     xmulast(1:2)=-99d9
     pdflast(-7:7,2)=-99d9
     pdlabellast(1:2)='abcdefg'
     ihlast(1:2)=-99
     init=.FALSE.
  ENDIF
!     Make sure we have a reasonable Bjorken x. Note that even though
!     x=0 is not reasonable, we prefer to simply return pdg2pdf=0
!     instead of stopping the code, as this might accidentally happen.
  IF(x.eq.0d0)THEN
     pdg2pdf=0d0
     RETURN
  ELSEIF(x.LT.0d0.OR.x.GT.1d0)THEN
     WRITE(*,*) 'PDF not supported for Bjorken x ', x
     STOP
  ENDIF

  ipart=ipdg
  IF(IABS(ipart).EQ.21) ipart=0
  IF(IABS(ipart).EQ.22) ipart=7
  iporg=ipart

!     This will be called for any PDG code, but we only support up to 7
  IF(IABS(ipart).GT.7)THEN
     WRITE(*,*) 'PDF not supported for pdg ',ipdg
     WRITE(*,*) 'For lepton colliders, please set the colpar'//&
          'variables to 3 in the user.inp'  
     STOP
  ENDIF

  ireuse = 0
  DO i=1,2
     ! Check if result can be reused since any of last two calls
     IF(x.EQ.xlast(i).AND.xmu.EQ.xmulast(i).AND.&
          pdlabel.EQ.pdlabellast(i).AND.ih.EQ.ihlast(i))THEN
        ireuse = i
     ENDIF
  ENDDO

  ! Reuse previous result, if possible
  IF(ireuse.GT.0)THEN
     IF(pdflast(iporg,ireuse).NE.-99d9)THEN
        pdg2pdf=pdflast(iporg,ireuse)
        RETURN
     ENDIF
  ENDIF

!     Bjorken x and/or facrorization scale and/or PDF set are not
!     identical to the saved values: this means a new event and we
!     should reset everything to compute new PDF values. Also, determine
!     if we should fill ireuse=1 or ireuse=2.
  IF(ireuse.EQ.0.AND.xlast(1).NE.-99d9.AND.xlast(2).NE.-99d9)THEN
     DO i=1,2
        xlast(i)=-99d9
        xmulast(i)=-99d9
        DO j=-7,7
           pdflast(j,i)=-99d9
        ENDDO
        pdlabellast(i)='abcdefg'
        ihlast(i)=-99
     ENDDO
!     everything has been reset. Now set ireuse=1 to fill the first
!     arrays of saved values below
     ireuse=1
  ELSEIF(ireuse.EQ.0.AND.xlast(1).NE.-99d9)THEN
!     This is first call after everything has been reset, so the first
!     arrays are already filled with the saved values (hence
!     xlast(1).ne.-99d9). Fill the second arrays of saved values (done
!     below) by setting ireuse=2
     ireuse=2
  ELSEIF(ireuse.EQ.0)THEN
!     Special: only used for the very first call to this function:
!     xlast(i) are initialized as data statements to be equal to -99d9
     ireuse=1
  ENDIF

!     Give the current values to the arrays that should be
!     saved. 'pdflast' is filled below.
  xlast(ireuse)=x
  xmulast(ireuse)=xmu
  pdlabellast(ireuse)=pdlabel
  ihlast(ireuse)=ih

  IF(IABS(ipart).EQ.7.AND.ih.GT.1)THEN
     ! photon PDF
     q2max=xmu*xmu
     IF(ih.EQ.3)THEN       !from the electron
        pdg2pdf=epa_electron(x,q2max)
     ELSEIF(ih.EQ.2)THEN !from a proton without breaking
        pdg2pdf=epa_proton(x,q2max)
     ENDIF
     pdflast(iporg,ireuse)=pdg2pdf
     RETURN
  ENDIF
      
  IF(pdlabel(1:5).EQ.'cteq3')THEN
     IF(pdlabel.EQ.'cteq3_m')THEN
        mode=1
     ELSEIF(pdlabel.EQ.'cteq3_l')THEN
        mode=2
     ELSEIF(pdlabel.EQ.'cteq3_d')THEN
        mode=3
     ENDIF

         
     IF(IABS(ipart).GE.1.AND.IABS(ipart).LE.2)&
          ipart=SIGN(3-IABS(ipart),ipart)
     pdg2pdf=Ctq3df(mode,ipart,x,xmu,Irt)/x
     ! valence part and sea part are seperately
     ! for u/d quark there are valence quark part and sea quark part 
     ! add sea quark part
     IF(ipdg.GE.1.AND.ipdg.LE.2)&
          pdg2pdf=pdg2pdf+Ctq3df(mode,-ipart,x,xmu,Irt)/x

  ELSEIF(pdlabel(1:5).EQ.'cteq4')THEN
     IF(pdlabel.EQ.'cteq4_m')THEN
        mode=1
     ELSEIF(pdlabel.EQ.'cteq4_d')THEN
        mode=2
     ELSEIF(pdlabel.EQ.'cteq4_l')THEN
        mode=3
     ELSEIF(pdlabel.EQ.'cteq4a1')THEN
        mode=4
     ELSEIF(pdlabel.EQ.'cteq4a2')THEN
        mode=5
     ELSEIF(pdlabel.EQ.'cteq4a3')THEN
        mode=6
     ELSEIF(pdlabel.EQ.'cteq4a4')THEN
        mode=7
     ELSEIF(pdlabel.EQ.'cteq4a5')THEN
        mode=8
     ELSEIF(pdlabel.EQ.'cteq4hj')THEN
        mode=9
     ELSEIF(pdlabel.EQ.'cteq4lq')THEN
        mode=10
     ENDIF
     ! u,d in cteq are different with pdg         
     IF(IABS(ipart).GE.1.AND.IABS(ipart).LE.2)&
          ipart=SIGN(3-iabs(ipart),ipart)
     pdg2pdf=Ctq4Fn(mode,ipart,x,xmu)

  ELSEIF(pdlabel.eq.'cteq5l1')THEN
     ! u,d in cteq are different with pdg
     IF(IABS(ipart).GE.1.AND.IABS(ipart).LE.2)&
          ipart=SIGN(3-IABS(ipart),ipart)

     pdg2pdf=Ctq5L(ipart,x,xmu)         
  ELSEIF((pdlabel(1:5).EQ.'cteq5').OR.& 
       (pdlabel(1:4).EQ.'ctq5'))THEN
     !  u,d in cteq are different with pdg
     IF(IABS(ipart).GE.1.AND.IABS(ipart).LE.2)&
          ipart=SIGN(3-IABS(ipart),ipart)

     pdg2pdf=Ctq5Pdf(ipart,x,xmu)                  
  ELSEIF(pdlabel(1:5).EQ.'cteq6')THEN
     !  u,d in cteq are different with pdg
     IF(IABS(ipart).GE.1.AND.IABS(ipart).LE.2)&
          ipart=SIGN(3-IABS(ipart),ipart)

     pdg2pdf=Ctq6Pdf(ipart,x,xmu)
  ELSE
     CALL pftopdg(ih,x,xmu,pdflast(-7:7,ireuse))
     pdg2pdf=pdflast(iporg,ireuse)
  ENDIF

  pdflast(iporg,ireuse)=pdg2pdf
  RETURN
END FUNCTION pdg2pdf

SUBROUTINE pftopdg(ih,x,q,pdf)
!***************************************************************************
!     Wrapper for calling the pdf of MCFM
!***************************************************************************
  IMPLICIT NONE
  REAL(KIND(1d0)),INTENT(IN)::x,q
  INTEGER,INTENT(IN)::ih
  REAL(KIND(1d0)),DIMENSION(-7:7),INTENT(OUT)::pdf

  CALL fdist(ih,x,q,pdf)
      
  RETURN
END SUBROUTINE pftopdg

SUBROUTINE fdist(ih,x,xmu,fx)
!***********************************************************************
!     MCFM PDF CALLING ROUTINE
!***********************************************************************
  IMPLICIT NONE
  INTEGER,INTENT(IN)::ih
  REAL(KIND(1d0)),INTENT(IN)::x,xmu
  REAL(KIND(1d0)),DIMENSION(-7:7)::fx
  REAL(KIND(1d0))::u_val,d_val,u_sea,d_sea,s_sea,c_sea,b_sea,gluon
  REAL(KIND(1d0)),EXTERNAL::Ctq3df,Ctq4Fn,Ctq5Pdf,Ctq6Pdf,Ctq5L
  REAL(KIND(1d0))::q2max
!  REAL(KIND(1d0))::epa_electron,epa_proton
  INTEGER::mode,Iprtn,Irt

  DO Iprtn=-7,7
     fx(Iprtn)=0d0
  ENDDO
!---set to zero if x out of range
  IF(x.GE.1d0)THEN
     RETURN
  ENDIF
       
  IF((pdlabel(1:3).EQ.'mrs')&
       .OR.(pdlabel(2:4).EQ.'mrs'))THEN
     IF(pdlabel.EQ.'mrs02nl')THEN
        mode=1
        CALL mrst2002(x,xmu,mode,u_val,d_val,u_sea,d_sea,&
             s_sea,c_sea,b_sea,gluon)
     ELSEIF(pdlabel.EQ.'mrs02nn')THEN
        mode=2
        CALL mrst2002(x,xmu,mode,u_val,d_val,u_sea,d_sea,&
             s_sea,c_sea,b_sea,gluon)
     ELSEIF(pdlabel.EQ.'mrs0119')THEN
        mode=1
        CALL mrst2001(x,xmu,mode,u_val,d_val,u_sea,d_sea,&
             s_sea,c_sea,b_sea,gluon)
     ELSEIF(pdlabel.EQ.'mrs0117')THEN
        mode=2
        CALL mrst2001(x,xmu,mode,u_val,d_val,u_sea,d_sea,&
             s_sea,c_sea,b_sea,gluon)
     ELSEIF(pdlabel.EQ.'mrs0121')THEN
        mode=3
        CALL mrst2001(x,xmu,mode,u_val,d_val,u_sea,d_sea,&
             s_sea,c_sea,b_sea,gluon)
     ELSEIF(pdlabel.EQ.'mrs01_j')THEN
        mode=4
        CALL mrst2001(x,xmu,mode,u_val,d_val,u_sea,d_sea,&
             s_sea,c_sea,b_sea,gluon)
     ELSEIF(pdlabel.EQ.'mrs99_1')THEN
        mode=1
        CALL mrs99(x,xmu,mode,u_val,d_val,u_sea,d_sea,&
             s_sea,c_sea,b_sea,gluon)
     ELSEIF(pdlabel.EQ.'mrs99_2')THEN
        mode=2
        CALL mrs99(x,xmu,mode,u_val,d_val,u_sea,d_sea,&
             s_sea,c_sea,b_sea,gluon)
     ELSEIF(pdlabel.EQ.'mrs99_3')THEN
        mode=3
        CALL mrs99(x,xmu,mode,u_val,d_val,u_sea,d_sea,&
             s_sea,c_sea,b_sea,gluon)
     ELSEIF(pdlabel.EQ.'mrs99_4')THEN
        mode=4
        CALL mrs99(x,xmu,mode,u_val,d_val,u_sea,d_sea,&
             s_sea,c_sea,b_sea,gluon)
     ELSEIF(pdlabel.EQ.'mrs99_5')THEN
        mode=5
        CALL mrs99(x,xmu,mode,u_val,d_val,u_sea,d_sea,&
             s_sea,c_sea,b_sea,gluon)
     ELSEIF(pdlabel.EQ.'mrs99_6')THEN
        mode=6
        CALL mrs99(x,xmu,mode,u_val,d_val,u_sea,d_sea,&
             s_sea,c_sea,b_sea,gluon)
     ELSEIF(pdlabel.EQ.'mrs99_7')THEN
        mode=7
        CALL mrs99(x,xmu,mode,u_val,d_val,u_sea,d_sea,&
             s_sea,c_sea,b_sea,gluon)
     ELSEIF(pdlabel.EQ.'mrs99_8')THEN
        mode=8
        CALL mrs99(x,xmu,mode,u_val,d_val,u_sea,d_sea,&
             s_sea,c_sea,b_sea,gluon)
     ELSEIF(pdlabel.EQ.'mrs99_9')THEN
        mode=9
        CALL mrs99(x,xmu,mode,u_val,d_val,u_sea,d_sea,&
             s_sea,c_sea,b_sea,gluon)
     ELSEIF(pdlabel.EQ.'mrs9910')THEN
        mode=10
        CALL mrs99(x,xmu,mode,u_val,d_val,u_sea,d_sea,&
             s_sea,c_sea,b_sea,gluon)
     ELSEIF(pdlabel.EQ.'mrs9911')THEN
        mode=11
        CALL mrs99(x,xmu,mode,u_val,d_val,u_sea,d_sea,&
             s_sea,c_sea,b_sea,gluon)
     ELSEIF(pdlabel.EQ.'mrs9912')THEN
        mode=12
        CALL mrs99(x,xmu,mode,u_val,d_val,u_sea,d_sea,&
             s_sea,c_sea,b_sea,gluon)
     ELSEIF(pdlabel.EQ.'mrs98z1')THEN
        mode=1
        CALL mrs98(x,xmu,mode,u_val,d_val,u_sea,d_sea,&
             s_sea,c_sea,b_sea,gluon)
     ELSEIF(pdlabel.EQ.'mrs98z2')THEN
        mode=2 
        CALL mrs98(x,xmu,mode,u_val,d_val,u_sea,d_sea,&
             s_sea,c_sea,b_sea,gluon)
     ELSEIF(pdlabel.EQ.'mrs98z3')THEN
        mode=3
        CALL mrs98(x,xmu,mode,u_val,d_val,u_sea,d_sea,&
             s_sea,c_sea,b_sea,gluon)
     ELSEIF(pdlabel.EQ.'mrs98z4')THEN
        mode=4
        CALL mrs98(x,xmu,mode,u_val,d_val,u_sea,d_sea,&
             s_sea,c_sea,b_sea,gluon)
     ELSEIF(pdlabel.EQ.'mrs98z5')THEN
        mode=5
        CALL mrs98(x,xmu,mode,u_val,d_val,u_sea,d_sea,&
             s_sea,c_sea,b_sea,gluon)
     ELSEIF(pdlabel.EQ.'mrs98l1')THEN
        mode=1
        CALL mrs98lo(x,xmu,mode,u_val,d_val,u_sea,d_sea,&
             s_sea,c_sea,b_sea,gluon)
     ELSEIF(pdlabel.EQ.'mrs98l2')THEN
        mode=2 
        CALL mrs98lo(x,xmu,mode,u_val,d_val,u_sea,d_sea,&
             s_sea,c_sea,b_sea,gluon)
     ELSEIF(pdlabel.EQ.'mrs98l3')THEN
        mode=3
        CALL mrs98lo(x,xmu,mode,u_val,d_val,u_sea,d_sea,&
             s_sea,c_sea,b_sea,gluon)
     ELSEIF(pdlabel.EQ.'mrs98l4')THEN
        mode=4
        CALL mrs98lo(x,xmu,mode,u_val,d_val,u_sea,d_sea,&
             s_sea,c_sea,b_sea,gluon)
     ELSEIF(pdlabel.EQ.'mrs98l5')THEN
        mode=5
        CALL mrs98lo(x,xmu,mode,u_val,d_val,u_sea,d_sea,&
             s_sea,c_sea,b_sea,gluon)
     ELSEIF(pdlabel.EQ.'mrs98ht')THEN
        mode=1
        CALL mrs98ht(x,xmu,mode,u_val,d_val,u_sea,d_sea,&
             s_sea,c_sea,b_sea,gluon)
     ENDIF
!-----assign mrs to standard grid
     fx(-5)=b_sea/x
     fx(-4)=c_sea/x
     fx(-3)=s_sea/x
     fx( 0)=gluon/x
     fx(+3)=fx(-3)
     fx(+4)=fx(-4)
     fx(+5)=fx(-5)
     ! valence part and sea part are seperately
     ! for u/d quark there are valence quark part and sea quark part
     fx(1)=(d_val+d_sea)/x
     fx(2)=(u_val+u_sea)/x
     fx(-1)=d_sea/x
     fx(-2)=u_sea/x

  ELSEIF(pdlabel(1:5).EQ.'cteq3')THEN     
     IF(pdlabel.EQ.'cteq3_m')THEN
        mode=1
     ELSEIF(pdlabel.EQ.'cteq3_l')THEN
        mode=2
     ELSEIF(pdlabel.EQ.'cteq3_d')THEN
        mode=3
     ENDIF
     fx(-5)=Ctq3df(mode,-5,x,xmu,Irt)/x
     fx(-4)=Ctq3df(mode,-4,x,xmu,Irt)/x
     fx(-3)=Ctq3df(mode,-3,x,xmu,Irt)/x
         
     fx(0)=Ctq3df(mode,0,x,xmu,Irt)/x
         
     fx(+3)=Ctq3df(mode,+3,x,xmu,Irt)/x
     fx(+4)=Ctq3df(mode,+4,x,xmu,Irt)/x
     fx(+5)=Ctq3df(mode,+5,x,xmu,Irt)/x
     fx(-1)=Ctq3df(mode,-2,x,xmu,Irt)/x
     fx(-2)=Ctq3df(mode,-1,x,xmu,Irt)/x
     ! valence part and sea part are seperately
     ! for u/d quark there are valence quark part and sea quark part
     fx(1)=Ctq3df(mode,+2,x,xmu,Irt)/x+fx(-1)
     fx(2)=Ctq3df(mode,+1,x,xmu,Irt)/x+fx(-2)
     
  ELSEIF(pdlabel(1:5).EQ.'cteq4')THEN
     IF(pdlabel.EQ.'cteq4_m')THEN
        mode=1
     ELSEIF(pdlabel.EQ.'cteq4_d')THEN
        mode=2
     ELSEIF(pdlabel.EQ.'cteq4_l')THEN
        mode=3
     ELSEIF(pdlabel.EQ.'cteq4a1')THEN
        mode=4
     ELSEIF(pdlabel.EQ.'cteq4a2')THEN
        mode=5
     ELSEIF(pdlabel.EQ.'cteq4a3')THEN
        mode=6
     ELSEIF(pdlabel.EQ.'cteq4a4')THEN
        mode=7
     ELSEIF(pdlabel.EQ.'cteq4a5')THEN
        mode=8
     ELSEIF(pdlabel.EQ.'cteq4hj')THEN
        mode=9
     ELSEIF(pdlabel.EQ.'cteq4lq')THEN
        mode=10
     ENDIF
         
     fx(-5)=Ctq4Fn(mode,-5,x,xmu)
     fx(-4)=Ctq4Fn(mode,-4,x,xmu)
     fx(-3)=Ctq4Fn(mode,-3,x,xmu)
         
     fx(0)=Ctq4Fn(mode,0,x,xmu)
     
     fx(+3)=Ctq4Fn(mode,+3,x,xmu)
     fx(+4)=Ctq4Fn(mode,+4,x,xmu)
     fx(+5)=Ctq4Fn(mode,+5,x,xmu)
     fx(1)=Ctq4Fn(mode,+2,x,xmu)
     fx(2)=Ctq4Fn(mode,+1,x,xmu)
     fx(-1)=Ctq4Fn(mode,-2,x,xmu)
     fx(-2)=Ctq4Fn(mode,-1,x,xmu)

  ELSEIF(pdlabel.EQ.'cteq5l1')THEN
     fx(-5)=Ctq5L(-5,x,xmu)
     fx(-4)=Ctq5L(-4,x,xmu)
     fx(-3)=Ctq5L(-3,x,xmu)
         
     fx(0)=Ctq5L(0,x,xmu)
     
     fx(+3)=Ctq5L(+3,x,xmu)
     fx(+4)=Ctq5L(+4,x,xmu)
     fx(+5)=Ctq5L(+5,x,xmu)
         
     fx(1)=Ctq5L(+2,x,xmu)
     fx(2)=Ctq5L(+1,x,xmu)
     fx(-1)=Ctq5L(-2,x,xmu)
     fx(-2)=Ctq5L(-1,x,xmu)
         
  ELSEIF((pdlabel(1:5).EQ.'cteq5').OR.& 
       (pdlabel(1:4).EQ.'ctq5'))THEN
         
     fx(-5)=Ctq5Pdf(-5,x,xmu)
     fx(-4)=Ctq5Pdf(-4,x,xmu)
     fx(-3)=Ctq5Pdf(-3,x,xmu)
         
     fx(0)=Ctq5Pdf(0,x,xmu)
     
     fx(+3)=Ctq5Pdf(+3,x,xmu)
     fx(+4)=Ctq5Pdf(+4,x,xmu)
     fx(+5)=Ctq5Pdf(+5,x,xmu)
         
     fx(1)=Ctq5Pdf(+2,x,xmu)
     fx(2)=Ctq5Pdf(+1,x,xmu)
     fx(-1)=Ctq5Pdf(-2,x,xmu)
     fx(-2)=Ctq5Pdf(-1,x,xmu)
                  
  ELSEIF(pdlabel(1:5).EQ.'cteq6')THEN
         
     fx(-5)=Ctq6Pdf(-5,x,xmu)
     fx(-4)=Ctq6Pdf(-4,x,xmu)
     fx(-3)=Ctq6Pdf(-3,x,xmu)
         
     fx(0)=Ctq6Pdf(0,x,xmu)
         
     fx(+3)=Ctq6Pdf(+3,x,xmu)
     fx(+4)=Ctq6Pdf(+4,x,xmu)
     fx(+5)=Ctq6Pdf(+5,x,xmu)
         
     fx(1)=Ctq6Pdf(+2,x,xmu)
     fx(2)=Ctq6Pdf(+1,x,xmu)
     fx(-1)=Ctq6Pdf(-2,x,xmu)
     fx(-2)=Ctq6Pdf(-1,x,xmu)
  ENDIF
!
!  a "diffractive" photon
!      
  q2max=xmu*xmu
  IF(ih.EQ.3)THEN  !from the electron
     fx(7)=epa_electron(x,q2max)
  ELSEIF(ih.EQ.2)THEN  !from a proton without breaking
     fx(7)=epa_proton(x,q2max)
  ENDIF
      
  RETURN
END SUBROUTINE fdist

! Note: electron structure function see QEDPS,i.e. qedps.pdf
FUNCTION electron_strf(x,q2max)
  ! electron structure function, see Eq.(15) in qedps.pdf
  IMPLICIT NONE
  REAL(KIND(1d0)),INTENT(IN)::x,q2max ! q2max is the Q^2 of total energy of e-e+
  REAL(KIND(1d0))::electron_strf
  REAL(KIND(1d0))::beta
  REAL(KIND(1d0))::xin=0.511d-3 ! electron mass in GeV
  REAL(KIND(1d0)),PARAMETER::PI=3.14159265358979323846d0
  REAL(KIND(1d0))::alpha,al,pi2,x1
  IF(aqedup.GT.0d0)THEN
     alpha=aqedup
  ELSE
     alpha = 0.0072992701d0
  ENDIF
  al=DLOG(q2max/xin**2)
  beta = 2d0*alpha/PI*(al-1d0) ! beta ~ 2 beta
  pi2=pi**2
  x1=1d0-x
  electron_strf=1d0+0.75d0*beta+beta**2/4d0*(9d0/8d0-pi2/3d0)&
       *beta*x1**(beta-1d0)+(-beta*(1d0-x1/2d0)+beta*beta/8d0*&
       (-4d0*(2d0-x1)*DLOG(x1)-(1d0+3d0*(1d0-x1)**2)/x1&
       *DLOG(1d0-x1)-6d0+x1))
  RETURN
END FUNCTION electron_strf
!/* ********************************************************* */
!/*  Equivalent photon approximation structure function.   * */
!/*     Improved Weizsaecker-Williams formula              * */
!/*   V.M.Budnev et al., Phys.Rep. 15C (1975) 181          * */
!/* ********************************************************* */
!   provided by Tomasz Pierzchala - UCL

FUNCTION epa_electron(x,q2max)
  ! photon PDF from electron
  IMPLICIT NONE
  REAL(KIND(1d0))::epa_electron
  REAL(KIND(1d0)),INTENT(IN)::x,q2max
  INTEGER::i
  REAL(KIND(1d0))::xin=0.511d-3 ! electron mass in GeV
  REAL(KIND(1d0))::alpha
  REAL(KIND(1d0))::f,q2min
  REAL(KIND(1d0)),PARAMETER::PI=3.14159265358979323846d0

  IF(aqedup.GT.0d0)THEN
     alpha=aqedup
  ELSE
     alpha = 0.0072992701d0
  ENDIF
!     // x = omega/E = (E-E')/E
  IF(x.LT.1)THEN
     q2min= xin*xin*x*x/(1-x)
     IF(q2min.LT.q2max)THEN 
        f = alpha/2d0/PI*&
             (2d0*xin*xin*x*(-1/q2min+1/q2max)+&
             (2-2d0*x+x*x)/x*dlog(q2max/q2min))            
     ELSE
        f = 0d0 
     ENDIF
  ELSE
     f= 0d0
  ENDIF
!      write (*,*) x,dsqrt(q2min),dsqrt(q2max),f
  IF(f.LT.0d0)f = 0d0
  epa_electron= f

END FUNCTION epa_electron

FUNCTION epa_proton(x,q2max)
  IMPLICIT NONE
  REAL(KIND(1d0))::epa_proton
  INTEGER::i
  REAL(KIND(1d0)),INTENT(IN)::x,q2max
  REAL(KIND(1d0))::xin=0.938d0 ! proton mass in GeV
  REAL(KIND(1d0))::alpha,qz
  REAL(KIND(1d0))::f,qmi,qma
  REAL(KIND(1d0)),PARAMETER::PI=3.14159265358979323846d0

  IF(aqedup.GT.0d0)THEN
     alpha = aqedup
  ELSE
     alpha = 0.0072992701d0
  ENDIF
  qz = 0.71d0
    
!     // x = omega/E = (E-E')/E
  IF(x.LT.1)THEN
     qmi= xin*xin*x*x/(1-x)
     IF(qmi.LT.q2max)THEN          
        f = alpha/PI*(phi_f(x,q2max/qz)-phi_f(x,qmi/qz))*(1-x)/x
     ELSE
        f=0d0
     ENDIF
  ELSE
     f= 0d0
  ENDIF
  IF(f.LT.0d0) f = 0d0
  epa_proton= f
END FUNCTION epa_proton

FUNCTION phi_f(x,qq)
  IMPLICIT NONE
  REAL(KIND(1d0))::phi_f
  REAL(KIND(1d0)),INTENT(IN)::x,qq
  REAL(KIND(1d0))::y,qq1,f,a,b,c

  a = 7.16d0
  b = -3.96d0
  c = 0.028d0

  qq1=1+qq
  y= x*x/(1-x)
  f=(1+a*y)*(-DLOG(qq1/qq)+1/qq1+1/(2*qq1*qq1)+1/(3*qq1*qq1*qq1))
  f=f + (1-b)*y/(4*qq*qq1*qq1*qq1);
  f=f+ c*(1+y/4)*(DLOG((qq1-b)/qq1)+b/qq1+b*b/(2*qq1*qq1)+&
       b*b*b/(3*qq1*qq1*qq1))
  phi_f= f
END FUNCTION phi_f
END MODULE Structf_PDFs
