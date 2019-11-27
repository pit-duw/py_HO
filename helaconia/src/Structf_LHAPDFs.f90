SUBROUTINE strf_lhapdf(wsf)
USE Helac_Global
USE setscale
INTEGER::ipp=1  ! ipp=1 ppbar ; ipp=2 pp
INTEGER::ih=1   ! ih=1 no photon PDF,ih=2, photon from proton/anti-proton,ih=3 photon from electron/positron 
REAL(KIND=DBL),INTENT(OUT)::wsf          
CHARACTER(len=20),DIMENSION(20)::parm
REAL(KIND=DBL),DIMENSION(20)::val
REAL(KIND=DBL),DIMENSION(-7:7)::pdflist
INTEGER::init=0,ijij      
SAVE init,ipp,ih
REAL(KIND=DBL)::glu1_ct,glu2_ct,u1_ct,u2_ct,d1_ct,d2_ct,s1_ct,s2_ct,c1_ct,c2_ct,b1_ct,b2_ct,&
                ub1_ct,ub2_ct,db1_ct,db2_ct,sb1_ct,sb2_ct,cb1_ct,cb2_ct,bb1_ct,bb2_ct,sf_ct
!REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE::glu1_et,glu2_et,u1_et,u2_et,d1_et,d2_et,s1_et,s2_et,c1_et,c2_et,b1_et,b2_et,&
!     ub1_et,ub2_et,db1_et,db2_et,sb1_et,sb2_et,cb1_et,cb2_et,bb1_et,bb2_et
IF(init.EQ.0)THEN
   IF(.NOT.lhapdfwrapinit)THEN
      CALL lhapdfwrap
      lhapdfwrapinit=.TRUE.
   ENDIF
   IF(Coll_Type.EQ.1)THEN
      ipp=2
   ELSE
      ipp=1
   ENDIF
!   IF(reweight_pdf.AND.ho_npdf.GT.0)THEN
!      ALLOCATE(wgtxsecpdf(ho_npdf),glu1_et(ho_npdf),glu2_et(ho_npdf),&
!           u1_et(ho_npdf),u2_et(ho_npdf),d1_et(ho_npdf),d2_et(ho_npdf),&
!           s1_et(ho_npdf),s2_et(ho_npdf),c1_et(ho_npdf),c2_et(ho_npdf),&
!           b1_et(ho_npdf),b2_et(ho_npdf),ub1_et(ho_npdf),ub2_et(ho_npdf),&
!           db1_et(ho_npdf),db2_et(ho_npdf),sb1_et(ho_npdf),sb2_et(ho_npdf),&
!           cb1_et(ho_npdf),cb2_et(ho_npdf),bb1_et(ho_npdf),bb2_et(ho_npdf))
!   ENDIF
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

CALL pftopdglha(ih,xp1,scale,pdflist(-7:7))
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
CALL pftopdglha(ih,xp2,scale,pdflist(-7:7))
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
!IF(reweight_pdf.AND.ho_npdf.GT.0)THEN
!   DO ijij=1,ho_npdf
!      CALL pdfset(parm,value)
!   ENDDO
   
!ENDIF
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
                   WRITE(*,*)"ERROR 1 in Structf_LHAPDFs.f90"
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
                   WRITE(*,*)"ERROR 2 in Structf_LHAPDFs.f90"
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
                   WRITE(*,*)"ERROR 3 in Structf_LHAPDFs.f90"
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

wsf=wjac*sf_ct/xp1/xp2
IF(init.EQ.0)init=1
       
END SUBROUTINE strf_lhapdf

FUNCTION pdg2lhapdf(ih,ipdg,x,xmu)
!***************************************************************************                          
!     wrapper for calling the lhapdf                                             
!***************************************************************************                          
  IMPLICIT NONE
  REAL(KIND(1d0))::pdg2lhapdf
  INTEGER,INTENT(IN)::ih,ipdg
  REAL(KIND(1d0)),INTENT(IN)::x,xmu
  INTEGER::i,j,ipart,iporg,ireuse,iset,imem
  INTEGER,DIMENSION(2)::ihlast,imemlast
  REAL(KIND(1d0)),DIMENSION(2)::xlast,xmulast
  REAL(KIND(1d0)),DIMENSION(-7:7,2)::pdflast
  LOGICAL::init=.TRUE.
  SAVE ihlast,xlast,xmulast,pdflast,imemlast,init
  IF(init)THEN
     ihlast(1:2)=-99
     xlast(1:2)=-99d9
     xmulast(1:2)=-99d9
     pdflast(-7:7,1:2)=-99d9
     imemlast(1:2)=-99
  ENDIF
  !     Make sure we have a reasonable Bjorken x. Note that even though                                 
  !     x=0 is not reasonable, we prefer to simply return pdg2pdf=0                                     
  !     instead of stopping the code, as this might accidentally happen.                                
  IF(x.EQ.0d0)THEN
     pdg2lhapdf=0d0
     RETURN
  ELSEIF(x.LT.0d0.OR.x.GT.1d0)THEN
     WRITE(*,*) 'PDF not supported for Bjorken x ', x
     STOP
  ENDIF

  ipart=ipdg
  IF(ipart.EQ.21) ipart=0
  IF(IABS(ipart).EQ.22) ipart=7
  iporg=ipart

  ! This will be called for any PDG code, but we only support up to 7                               
  IF(IABS(ipart).GT.7)THEN
     WRITE(*,*) 'PDF not supported for pdg ',ipdg
     STOP
  ENDIF

  ! Determine the iset used in lhapdf                                                               
  CALL getnset(iset)
  IF(iset.NE.1)THEN
     WRITE(*,*) 'PDF not supported for Bjorken x ', x
     STOP
  ENDIF

  ! Determine the member of the set (function of lhapdf)                                            
  CALL getnmem(iset,imem)

  ireuse = 0
  DO i=1,2
     ! Check if result can be reused since any of last two calls                                       
     IF(x.EQ.xlast(i) .AND. xmu.EQ.xmulast(i) .AND.&
          imem.EQ.imemlast(i) .AND. ih.EQ.ihlast(i))THEN
        ireuse = i
     ENDIF
  ENDDO
  !  Reuse previous result, if possible                                                              
  IF(ireuse.GT.0)THEN
     IF(pdflast(iporg,ireuse).NE.-99d9)THEN
        pdg2lhapdf=pdflast(iporg,ireuse)
        RETURN
     ENDIF
  ENDIF

  ! Bjorken x and/or facrorization scale and/or PDF set are not                                     
  ! identical to the saved values: this means a new event and we                                    
  ! should reset everything to compute new PDF values. Also, determine                              
  ! if we should fill ireuse=1 or ireuse=2.                                                         
  IF(ireuse.EQ.0.AND.xlast(1).NE.-99d9.AND.xlast(2).NE.-99d9)THEN
     DO i=1,2
        xlast(i)=-99d9
        xmulast(i)=-99d9
        DO j=-7,7
           pdflast(j,i)=-99d9
        ENDDO
        imemlast(i)=-99
        ihlast(i)=-99
     ENDDO
     ! everything has been reset. Now set ireuse=1 to fill the first                                   
     ! arrays of saved values below                                                                    
     ireuse=1
  ELSEIF(ireuse.EQ.0.AND.xlast(1).NE.-99d9)THEN
     ! This is first call after everything has been reset, so the first                                
     ! arrays are already filled with the saved values (hence                                          
     ! xlast(1).ne.-99d9). Fill the second arrays of saved values (done                                
     ! below) by setting ireuse=2                                                                      
     ireuse=2
  ELSEIF(ireuse.EQ.0)THEN
     ! Special: only used for the very first call to this function:                                    
     ! xlast(i) are initialized as data statements to be equal to -99d9                                
     ireuse=1
  ENDIF

  ! Call lhapdf and give the current values to the arrays that should                               
  ! be saved                                                                                        
  CALL pftopdglha(ih,x,xmu,pdflast(-7:7,ireuse))
  xlast(ireuse)=x
  xmulast(ireuse)=xmu
  ihlast(ireuse)=ih
  imemlast(ireuse)=imem
  pdg2lhapdf=pdflast(ipart,ireuse)
  RETURN
END FUNCTION pdg2lhapdf

SUBROUTINE lhapdfwrap
  USE Helac_Global
  IMPLICIT NONE
  REAL(KIND(1d0)),PARAMETER::zmass=91.188d0
  CHARACTER(len=20),DIMENSION(20)::parm
  REAL(KIND(1d0)),DIMENSION(20)::value
  REAL(KIND(1d0)),EXTERNAL::alphasPDF ! the lhapdf internal function
  ! initialize the pdf set                                                                          
  CALL FindPDFPath(LHAPath)
  CALL SetPDFPath(LHAPath) ! Location of the LHAPDF library of PDFs
  ! they are defined in share/lhapdf/PDFsets.index
  value(1)=iPDFSUP1   ! The LHAPDF set/number is selected depending on ABS(INT(value(1))) in HERWIG and Standalone
  parm(1)='DEFAULT' ! PYTHIA: 'NPTYPE'; HERWIG: 'HWLHAPDF'; Standalone "DEFAULT'
  CALL pdfset(parm,value) ! PDFSET is a function by using lhaglue interface in lhapdf
  CALL GetOrderAs(nloop) ! the order of running alpha_S, 0:LO,1:NLO etc
  nloop=nloop+1 ! nloop=1:LO; nloop=2:NLO
  alphaQCD2=alphasPDF(zmass) ! internal function of LHAPDF
  RETURN
END SUBROUTINE lhapdfwrap

SUBROUTINE FindPDFPath(LHAPath)
!********************************************************************
! generic subroutine to open the table files in the right directories
!********************************************************************
  IMPLICIT NONE
  CHARACTER(len=150),INTENT(INOUT)::LHAPath
  CHARACTER(len=3)::up
  LOGICAL::exists
  INTEGER::i
  ! try in the LHAPath first
  INQUIRE(File=LHAPath,exist=exists)
  IF(exists)RETURN
  !first try in the current directory                                                              
  LHAPath='PDFsets'
  INQUIRE(File=LHAPath, exist=exists)
  IF(exists)RETURN
  ! then try one directory up                                                                       
  LHAPath=up//LHAPath
  INQUIRE(File=LHAPath, exist=exists)
  IF(exists)RETURN
  ! finally try in the lib directory                                                                
  LHAPath='lib/PDFsets'
  INQUIRE(File=LHAPath, exist=exists)
  IF(exists)RETURN
  DO i=1,6
     LHAPath=up//LHAPath
     INQUIRE(File=LHAPath, exist=exists)
     IF(exists)RETURN
  ENDDO
  PRINT *,'ERROR:Could not find PDFsets directory in LHAPDF, quitting'
  STOP

  RETURN
END SUBROUTINE FindPDFPath

SUBROUTINE pftopdglha(ih,x,q,pdf)
!***************************************************************************                          
!     Wrapper for calling the pdf from lhapdf                                                         
!***************************************************************************                          
  IMPLICIT NONE
  REAL(KIND(1d0)),INTENT(IN)::x,q
  REAL(KIND(1d0)),DIMENSION(-7:7),INTENT(OUT)::pdf
  REAL(KIND(1d0)),DIMENSION(-6:6)::f ! momentum density x*PDF
  INTEGER,INTENT(IN)::IH
  INTEGER::I
  REAL(KIND(1d0))::photon
  LOGICAL,EXTERNAL::has_photon ! the lhapdf internal function

  IF(ABS(ih).EQ.1)THEN
     pdf(-7)=0d0
     IF(has_photon())THEN
        CALL evolvePDFphoton(x, q, f, photon)
        pdf(7)= photon
     ELSE
        pdf(7) = 0d0
        CALL evolvePDF(x, q, f)
     ENDIF
     DO i=-6,6
        pdf(i)=f(i)/x
     ENDDO
  ELSE
     WRITE(*,*) 'beam type is not supported in lhadpf'
     DO i=-7,7
        pdf(i)=0d0
     ENDDO
  ENDIF

  RETURN
END SUBROUTINE pftopdglha

SUBROUTINE get_pdfup(pdfin,pdfgup,pdfsup,lhaid)
!-------------------------------------------------
!   Convert pdlabel name to LHAPDF number 
!------------------------------------------------- 
  IMPLICIT NONE
  CHARACTER(len=*),INTENT(IN)::pdfin
  INTEGER::mpdf
  INTEGER,DIMENSION(2),INTENT(OUT)::PDFGUP,PDFSUP
  INTEGER,INTENT(IN)::lhaid
  INTEGER,PARAMETER::npdfs=13
  INTEGER::i
  CHARACTER(len=7),DIMENSION(npdfs)::pdflabs=(/'none   ',&
       'mrs02nl',&
       'mrs02nn',&
       'cteq4_m',&
       'cteq4_l',&
       'cteq4_d',&
       'cteq5_m',&
       'cteq5_d',&
       'cteq5_l',&
       'cteq5m1',&
       'cteq6_m',&
       'cteq6_l',&
       'cteq6l1'/)
  INTEGER,DIMENSION(npdfs)::numspdf=(/00000,&
       20250,&
       20270,&
       19150,&
       19170,&
       19160,&
       19050,&
       19060,&
       19070,&
       19051,&
       10000,&
       10041,&
       10042/)
  
  IF(pdfin.EQ."lhapdf")THEN
     WRITE(*,*)'using LHAPDF'
     DO i=1,2
        pdfgup(i)=0
        pdfsup(i)=lhaid
     ENDDO
     RETURN
  ENDIF


  mpdf=-1
  DO i=1,npdfs
     IF(pdfin(1:LEN_TRIM(pdfin)) .EQ. pdflabs(i))THEN
        mpdf=numspdf(i)
     ENDIF
  ENDDO

  IF(mpdf.EQ.-1)THEN
     WRITE(*,*)'pdf ',pdfin,' not implemented in get_pdfup.'
     WRITE(*,*)'known pdfs are'
     WRITE(*,*) pdflabs
     WRITE(*,*)'using ',pdflabs(12)
     mpdf=numspdf(12)
  ENDIF

  DO i=1,2
     pdfgup(i)=0
     pdfsup(i)=mpdf
  ENDDO

  RETURN
END SUBROUTINE get_pdfup
