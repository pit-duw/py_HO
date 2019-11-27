MODULE Helac_wavef
USE Helac_Global
USE Helac_Func_1
IMPLICIT NONE
CONTAINS
!*****************************************************************
!  Be careful that dconjg(eminus) =/= eplus
!  because we are working in the rotated lightcone frame
!  where for the 3- and 4-  component there is an extra
!  imaginary unit in its definition
!  in lightcone frame p=(p0+pz,p0-pz,px+ipy,px-ipy)
!*****************************************************************
! external wave_functions for vector bosons
!SUBROUTINE Helac_wave_functions()
!IMPLICIT NONE
!COMPLEX(KIND=DBL),DIMENSION(4)::zy
!COMPLEX(KIND=DBL),DIMENSION(5)::zq
!return
! here quantization aixes are all its momentum,i.e. helicity and gamma=0
! entry
! Eq(I.2.41) in arXiv:0812.1594v2 [hep-ph]
! in normal frame
! eps(-1)=(0,cos(th)cos(phi)+i sin(phi),cos(th)sin(phi)-i cos(phi),-sin(th))/Sqrt(2)
! zy1 returns the lightcone frame form
SUBROUTINE Helac_eminus(zy1,zq1)
COMPLEX(KIND=DBL),DIMENSION(5),INTENT(IN)::zq1
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(OUT)::zy1
REAL(KIND=DBL)::pt,pz,p,p0       
pt=DIMAG(zq1(2))
pz=DIMAG(zq1(1))
p =DIMAG(zq1(5))
p0=DREAL(zq1(1))-pz
      
IF(p.GT.eps(1))THEN
   IF(pt.GT.eps(1))THEN
       zy1(1)=-DCMPLX(pt/p/DSQRT(dnou(2)))
       zy1(2)= DCMPLX(pt/p/DSQRT(dnou(2)))
       zy1(3)=zq1(3)*DCMPLX(( p+pz)/p/pt/DSQRT(dnou(2)))
       zy1(4)=zq1(4)*DCMPLX((-p+pz)/p/pt/DSQRT(dnou(2)))
   ELSE
   ! we choose phi=0 here!
       IF(pz.GT.dnou(0))THEN
          zy1(1)=DCMPLX(dnou(0),dnou(0))
          zy1(2)=DCMPLX(dnou(0),dnou(0))
          zy1(3)=DCMPLX(DSQRT(dnou(2)),dnou(0))
          zy1(4)=DCMPLX(dnou(0),dnou(0))
       ELSE
          zy1(1)=DCMPLX(dnou(0),dnou(0))
          zy1(2)=DCMPLX(dnou(0),dnou(0))
          zy1(3)=DCMPLX(dnou(0),dnou(0))
          zy1(4)=DCMPLX(DSQRT(dnou(2)),dnou(0))
       ENDIF
   ENDIF
ELSE
   IF(p0.GT.eps(1))THEN
   ! we choose z direction
       zy1(1)=DCMPLX(dnou(0),dnou(0))
       zy1(2)=DCMPLX(dnou(0),dnou(0))
       zy1(3)=DCMPLX(DSQRT(dnou(2)),dnou(0))
       zy1(4)=DCMPLX(dnou(0),dnou(0))
   ELSE
       zy1(1)=DCMPLX(dnou(0),dnou(0))
       zy1(2)=DCMPLX(dnou(0),dnou(0))
       zy1(3)=DCMPLX(dnou(0),dnou(0))
       zy1(4)=DCMPLX(dnou(0),dnou(0))
   ENDIF
ENDIF
END SUBROUTINE Helac_eminus
! Eq(I.2.41) in arXiv:0812.1594v2 [hep-ph]
! in normal frame
! eps(+1)=(0,-cos(th)cos(phi)+i sin(phi),-cos(th)sin(phi)-i cos(phi),+sin(th))/Sqrt(2)
! zy1 returns the lightcone frame form       
SUBROUTINE Helac_eplus(zy1,zq1)
COMPLEX(KIND=DBL),DIMENSION(5),INTENT(IN)::zq1
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(OUT)::zy1
REAL(KIND=DBL)::pt,pz,p,p0      
pt=DIMAG(zq1(2))
pz=DIMAG(zq1(1))
p =DIMAG(zq1(5))
p0=DREAL(zq1(1))-pz
IF(p.GT.eps(1))THEN
   IF(pt.GT.eps(1))THEN
      zy1(1)= DCMPLX(pt/p/DSQRT(dnou(2)))
      zy1(2)=-DCMPLX(pt/p/dsqrt(dnou(2)))
      zy1(3)=zq1(3)*DCMPLX(( p-pz)/p/pt/DSQRT(dnou(2)))
      zy1(4)=zq1(4)*DCMPLX((-p-pz)/p/pt/DSQRT(dnou(2)))
   ELSE
      IF(pz.GT.dnou(0))THEN
        zy1(1)=DCMPLX(dnou(0),dnou(0))
        zy1(2)=DCMPLX(dnou(0),dnou(0))
        zy1(3)=DCMPLX(dnou(0),dnou(0))
        zy1(4)=-DCMPLX(DSQRT(dnou(2)),dnou(0))
       ELSE
        zy1(1)=DCMPLX(dnou(0),dnou(0))
        zy1(2)=DCMPLX(dnou(0),dnou(0))
        zy1(3)=-DCMPLX(DSQRT(dnou(2)),dnou(0))
        zy1(4)=DCMPLX(dnou(0),dnou(0))
       ENDIF
    ENDIF
ELSE
    IF(p0.GT.eps(1))THEN
       zy1(1)=DCMPLX(dnou(0),dnou(0))
       zy1(2)=DCMPLX(dnou(0),dnou(0))
       zy1(3)=DCMPLX(dnou(0),dnou(0))
       zy1(4)=-DCMPLX(DSQRT(dnou(2)),dnou(0))
    ELSE
       zy1(1)=DCMPLX(dnou(0),dnou(0))
       zy1(2)=DCMPLX(dnou(0),dnou(0))
       zy1(3)=DCMPLX(dnou(0),dnou(0))
       zy1(4)=DCMPLX(dnou(0),dnou(0))
    ENDIF
ENDIF
END SUBROUTINE Helac_eplus
! Eq(I.2.43) in arXiv:0812.1594v2 [hep-ph]
! in normal frame
! eps(+1)=(p,E sin(th) cos(phi),E sin(th)sin(phi),E cos(th))/m
! zy1 returns the lightcone frame form 
SUBROUTINE Helac_elong(zy1,zq1)
COMPLEX(KIND=DBL),DIMENSION(5),INTENT(IN)::zq1
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(OUT)::zy1
REAL(KIND=DBL)::pt,pz,p,p0,p2      
pt=DIMAG(zq1(2))
pz=DIMAG(zq1(1))
p =DIMAG(zq1(5))
p0=DREAL(zq1(1))-pz       
p2=p0*p0-p*p
IF(p2.GT.eps(1))THEN
   p2=DSQRT(p2)
   IF(p.GT.eps(1))THEN
       zy1(1)=DCMPLX(p/p2+pz*p0/(p*p2))
       zy1(2)=DCMPLX(p/p2-pz*p0/(p*p2))
       zy1(3)=zq1(3)*DCMPLX(p0/(p*p2))
       zy1(4)=zq1(4)*DCMPLX(p0/(p*p2))
   ELSE
       zy1(1)=DCMPLX(dnou(1),dnou(0))
       zy1(2)=DCMPLX(-dnou(1),dnou(0))
       zy1(3)=DCMPLX(dnou(0),dnou(0))
       zy1(4)=DCMPLX(dnou(0),dnou(0))
   ENDIF
ELSE
   zy1(1:4)=DCMPLX(dnou(0),dnou(0))
ENDIF
END SUBROUTINE Helac_elong   
    
! external wave_functions for fermions
! Eqs(2.75)-(2.76),(3.1.19)-(3.1.22),(C.1.11),(G.4.11)-(G.4.12)
! in arXiv:0812.1594v2 [hep-ph] with quantization aix its momentum (i.e.helicity)  and gamma=0
! It is a spinor with  momentum (p0,p sin(th) cos(phi),p sin(th) sin(phi),p cos(th))
! and we choose gamma=-phi in  (C.111).
SUBROUTINE Helac_uplus(zy1,zq1)
COMPLEX(KIND=DBL),DIMENSION(5),INTENT(IN)::zq1
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(OUT)::zy1
REAL(KIND=DBL)::emass,a,b,c,r      
emass=DREAL(zq1(5))
IF(emass.NE.dnou(0))THEN
   !b=p+pz
   b=DIMAG(zq1(5))+DIMAG(zq1(1))
   IF(b.GT.eps(1))THEN
   !a=p+p0
      a=b+DREAL(zq1(2))
   !c=2*p
      c=a+b-DREAL(zq1(1))
      r=DSQRT(a*b*c)
	  zy1(1)= DCMPLX(r/c)
      zy1(2)= zq1(3)*DCMPLX(a/r)
      zy1(3)=-DCMPLX(emass*b/r)
      zy1(4)=-zq1(3)*DCMPLX(emass/r)
    ELSE
	  !a=p+p0
	  a=b+DREAL(zq1(2))
      zy1(1)=DCMPLX(dnou(0),dnou(0))
      zy1(2)=-DCMPLX(DSQRT(a))
      zy1(3)=DCMPLX(dnou(0),dnou(0))
      zy1(4)=DCMPLX(emass/DSQRT(a))
    ENDIF
ELSE
    CALL Helac_ru(zy1,zq1)
ENDIF
END SUBROUTINE Helac_uplus

SUBROUTINE Helac_ubplus(zy1,zq1)
COMPLEX(KIND=DBL),DIMENSION(5),INTENT(IN)::zq1
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(OUT)::zy1
REAL(KIND=DBL)::emass,a,b,c,r      
emass=DREAL(zq1(5))
IF(emass.NE.dnou(0))THEN
   b=DIMAG(zq1(5))+DIMAG(zq1(1))
   !b=p+pz
   IF(b.GT.eps(1))THEN
   !a=p+p0
     a=b+DREAL(zq1(2))
   !c=2*p
     c=a+b-DREAL(zq1(1))
     r=DSQRT(a*b*c)
     zy1(1)= DCMPLX(emass*b/r)
     zy1(2)= zq1(4)*DCMPLX(emass/r)
     zy1(3)=-DCMPLX(r/c)
     zy1(4)=-zq1(4)*DCMPLX(a/r)
   ELSE
     a=b+DREAL(zq1(2))
     zy1(1)=DCMPLX(dnou(0),dnou(0))
     zy1(2)=-DCMPLX(emass/DSQRT(a))
     zy1(3)=DCMPLX(dnou(0),dnou(0))
     zy1(4)= DCMPLX(DSQRT(a))
   ENDIF
ELSE
   CALL Helac_rub(zy1,zq1)
ENDIF
END SUBROUTINE Helac_ubplus
       
SUBROUTINE Helac_uminus(zy1,zq1)
COMPLEX(KIND=DBL),DIMENSION(5),INTENT(IN)::zq1
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(OUT)::zy1
REAL(KIND=DBL)::emass,a,b,c,r        
emass=DREAL(zq1(5))
IF(emass.NE.dnou(0))THEN
   b=DIMAG(zq1(5))+DIMAG(zq1(1))
   IF(b.GT.eps(1))THEN
      a=b+DREAL(zq1(2))
      c=a+b-DREAL(zq1(1))
      r=DSQRT(a*b*c)
      zy1(1)= zq1(4)*DCMPLX(emass/r)
      zy1(2)=-DCMPLX(emass*b/r)
      zy1(3)=-zq1(4)*DCMPLX(a/r)
      zy1(4)=DCMPLX(r/c)
   ELSE
      a=b+DREAL(zq1(2))
      zy1(1)=-DCMPLX(emass/DSQRT(a))
      zy1(2)=DCMPLX(dnou(0),dnou(0))
      zy1(3)= DCMPLX(DSQRT(a))
      zy1(4)=DCMPLX(dnou(0),dnou(0))
   ENDIF
ELSE
   CALL Helac_lu(zy1,zq1)
ENDIF
END SUBROUTINE Helac_uminus
       
SUBROUTINE Helac_ubminus(zy1,zq1)
COMPLEX(KIND=DBL),DIMENSION(5),INTENT(IN)::zq1
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(OUT)::zy1
REAL(KIND=DBL)::emass,a,b,c,r        
emass=DREAL(zq1(5))
IF(emass.NE.dnou(0))THEN
  b=DIMAG(zq1(5))+DIMAG(zq1(1))
  IF(b.GT.eps(1))THEN
     a=b+DREAL(zq1(2))
     c=a+b-DREAL(zq1(1))
     r=DSQRT(a*b*c)

     zy1(1)= zq1(3)*DCMPLX(a/r)
     zy1(2)=-DCMPLX(r/c)
     zy1(3)=-zq1(3)*DCMPLX(emass/r)
     zy1(4)= DCMPLX(emass*b/r)
  ELSE
     a=b+DREAL(zq1(2))
     zy1(1)=-DSQRT(a)
     zy1(2)=DCMPLX(dnou(0),dnou(0))
     zy1(3)= DCMPLX(emass/DSQRT(a))
     zy1(4)=DCMPLX(dnou(0),dnou(0))
  ENDIF
ELSE
  CALL Helac_lub(zy1,zq1)
ENDIF
END SUBROUTINE Helac_ubminus
       
SUBROUTINE Helac_vplus(zy1,zq1)
COMPLEX(KIND=DBL),DIMENSION(5),INTENT(IN)::zq1
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(OUT)::zy1
REAL(KIND=DBL)::emass,a,b,c,r         
emass=DREAL(zq1(5))
IF(emass.NE.dnou(0))THEN
   b=DIMAG(zq1(5))+DIMAG(zq1(1))   
   IF(b.GT.eps(1))THEN
     a=b+DREAL(zq1(2))
     c=a+b-DREAL(zq1(1))
     r=DSQRT(a*b*c)

     zy1(1)=-zq1(4)*DCMPLX(emass/r) 
     zy1(2)= DCMPLX(emass*b/r)
     zy1(3)=-zq1(4)*DCMPLX(a/r)
     zy1(4)= DCMPLX(r/c)
   ELSE
     a=b+DREAL(zq1(2))
     zy1(1)= DCMPLX(emass/DSQRT(a))
     zy1(2)=DCMPLX(dnou(0),dnou(0))
     zy1(3)=DCMPLX(DSQRT(a))
     zy1(4)=DCMPLX(dnou(0),dnou(0))
   ENDIF
ELSE
   CALL Helac_lu(zy1,zq1)
ENDIF
END SUBROUTINE Helac_vplus 

SUBROUTINE Helac_vbplus(zy1,zq1)
COMPLEX(KIND=DBL),DIMENSION(5),INTENT(IN)::zq1
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(OUT)::zy1
REAL(KIND=DBL)::emass,a,b,c,r      
emass=DREAL(zq1(5))
IF(emass.NE.dnou(0))THEN
   b=DIMAG(zq1(5))+DIMAG(zq1(1))
   IF(b.GT.eps(1))THEN
      a=b+DREAL(zq1(2))
      c=a+b-DREAL(zq1(1))
      r=DSQRT(a*b*c)
      zy1(1)= zq1(3)*DCMPLX(a/r) 
      zy1(2)=-DCMPLX(r/c)
      zy1(3)= zq1(3)*DCMPLX(emass/r)
      zy1(4)=-DCMPLX(emass*b/r)
   ELSE
      a=b+DREAL(zq1(2))
      zy1(1)=-DCMPLX(DSQRT(a))
      zy1(2)=DCMPLX(dnou(0),dnou(0))
      zy1(3)=-DCMPLX(emass/DSQRT(a))
      zy1(4)=DCMPLX(dnou(0),dnou(0))
   ENDIF
ELSE
   CALL Helac_lub(zy1,zq1)
ENDIF
END SUBROUTINE Helac_vbplus
       
SUBROUTINE Helac_vminus(zy1,zq1)
COMPLEX(KIND=DBL),DIMENSION(5),INTENT(IN)::zq1
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(OUT)::zy1
REAL(KIND=DBL)::emass,a,b,c,r        
emass=DREAL(zq1(5))
IF(emass.NE.dnou(0))THEN
   b=DIMAG(zq1(5))+DIMAG(zq1(1))
   IF(b.GT.eps(1))THEN
      a=b+DREAL(zq1(2))
      c=a+b-DREAL(zq1(1))
      r=DSQRT(a*b*c)
      zy1(1)= DCMPLX(r/c)
      zy1(2)= zq1(3)*DCMPLX(a/r)
      zy1(3)= emass*DCMPLX(b/r)
      zy1(4)= zq1(3)*DCMPLX(emass/r)
   ELSE
      a=b+DREAL(zq1(2))
      zy1(1)=DCMPLX(dnou(0),dnou(0))
      zy1(2)=-DCMPLX(DSQRT(a))
      zy1(3)=DCMPLX(dnou(0),dnou(0))
      zy1(4)=-DCMPLX(emass/DSQRT(a))
   ENDIF
ELSE
   CALL Helac_ru(zy1,zq1)
ENDIF
END SUBROUTINE Helac_vminus

SUBROUTINE Helac_vbminus(zy1,zq1)
COMPLEX(KIND=DBL),DIMENSION(5),INTENT(IN)::zq1
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(OUT)::zy1
REAL(KIND=DBL)::emass,a,b,c,r         
emass=DREAL(zq1(5))
IF(emass.NE.dnou(0))THEN
   b=DIMAG(zq1(5))+DIMAG(zq1(1))
   IF(b.GT.eps(1))THEN
     a=b+DREAL(zq1(2))
     c=a+b-DREAL(zq1(1))
     r=DSQRT(a*b*c)
     zy1(1)=-DCMPLX(emass*b/r)
     zy1(2)=-zq1(4)*DCMPLX(emass/r)
     zy1(3)=-DCMPLX(r/c)
     zy1(4)=-zq1(4)*DCMPLX(a/r)
   ELSE
     a=b+DREAL(zq1(2))
     zy1(1)=DCMPLX(dnou(0),dnou(0))
     zy1(2)= DCMPLX(emass/DSQRT(a))
     zy1(3)=DCMPLX(dnou(0),dnou(0))
     zy1(4)= DCMPLX(DSQRT(a))
   ENDIF
ELSE
   CALL Helac_rub(zy1,zq1)
ENDIF
END SUBROUTINE Helac_vbminus

! massless spinors       
SUBROUTINE Helac_ru(zy1,zq1)
COMPLEX(KIND=DBL),DIMENSION(5),INTENT(IN)::zq1
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(OUT)::zy1
REAL(KIND=DBL)::pp,a
zy1(3)=DCMPLX(dnou(0),dnou(0))
zy1(4)=DCMPLX(dnou(0),dnou(0))
pp=DREAL(zq1(1))
IF(pp.GT.eps(1))THEN
  a=DSQRT(pp)
  zy1(2)= zq1(3)/DCMPLX(a)
  zy1(1)= DCMPLX(a)
ELSE
  zy1(2)=DCMPLX(-DSQRT(DREAL(zq1(2))),dnou(0))
  zy1(1)=DCMPLX(dnou(0),dnou(0))
ENDIF
END SUBROUTINE Helac_ru       
       
SUBROUTINE Helac_lu(zy1,zq1)
COMPLEX(KIND=DBL),DIMENSION(5),INTENT(IN)::zq1
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(OUT)::zy1
REAL(KIND=DBL)::pp,a        
zy1(1)=DCMPLX(dnou(0),dnou(0))
zy1(2)=DCMPLX(dnou(0),dnou(0))
!p0+pz
pp=DREAL(zq1(1))
IF(pp.GT.eps(1))THEN
   a=DSQRT(pp)
   zy1(3)=-zq1(4)/DCMPLX(a)
   zy1(4)=DCMPLX(a)
ELSE  
   zy1(3)=DCMPLX(DSQRT(DREAL(zq1(2))),dnou(0))
   zy1(4)=DCMPLX(dnou(0),dnou(0))
ENDIF
END SUBROUTINE Helac_lu
     
SUBROUTINE Helac_lub(zy1,zq1)
COMPLEX(KIND=DBL),DIMENSION(5),INTENT(IN)::zq1
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(OUT)::zy1
REAL(KIND=DBL)::pp,a 
zy1(3)=DCMPLX(dnou(0),dnou(0))
zy1(4)=DCMPLX(dnou(0),dnou(0))
pp=DREAL(zq1(1))
IF(pp.GT.eps(1))THEN
  a=DSQRT(pp)
  zy1(1)= zq1(3)/DCMPLX(a)
  zy1(2)=-DCMPLX(a)
ELSE
  zy1(1)=DCMPLX(-DSQRT(DREAL(zq1(2))),dnou(0))
  zy1(2)=DCMPLX(dnou(0),dnou(0))
ENDIF
END SUBROUTINE Helac_lub
       
SUBROUTINE Helac_rub(zy1,zq1)
COMPLEX(KIND=DBL),DIMENSION(5),INTENT(IN)::zq1
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(OUT)::zy1
REAL(KIND=DBL)::pp,a 
zy1(1)=DCMPLX(dnou(0),dnou(0))
zy1(2)=DCMPLX(dnou(0),dnou(0))
pp=DREAL(zq1(1))
IF(pp.GT.eps(1))THEN
   a=DSQRT(pp)
   zy1(4)=-zq1(4)/DCMPLX(a)
   zy1(3)=-DCMPLX(a)   
ELSE      
   zy1(4)= DCMPLX(DSQRT(DREAL(zq1(2))),dnou(0))
   zy1(3)=DCMPLX(dnou(0),dnou(0))
ENDIF
END SUBROUTINE Helac_rub

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
! External wave functions for hadrons production                              !
! in nonrealavistic QCD  up to P-wave                                         !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE  Helac_Qvplus(zy1,zq1,zP,spin,nderi,llz,ssz)
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(OUT)::zy1
COMPLEX(KIND=DBL),DIMENSION(5),INTENT(IN)::zq1,zP
INTEGER,INTENT(IN)::spin,nderi,llz,ssz ! spin , rank of derivation, L_z,S_z
COMPLEX(KIND=DBL),DIMENSION(4)::zPy,zqq,zeps
COMPLEX(KIND=DBL),DIMENSION(5)::zP2
COMPLEX(KIND=DBL)::zmm
zmm=DCMPLX(DREAL(zq1(5)))
CALL Helac_uplus(zPy,zP)
! return from (p0,-p) to (p0,p) as the definition in 
! Eq.(A16) in Phys.Rev.D  57 (1998)  4258
!IF(PolarFrame.EQ.1)THEN
zP2(2)=DCMPLX(-DREAL(zP(2)),DIMAG(zP(2)))! p0-pz+i*pt  incoming -p0+pz+i*pt outgoing
zP2(1)=-zP(1)                  ! po+pz+i*pz incoming -po-pz-i*pz outgoing
zP2(3)=-zP(3)                  ! px+i*py incoming -px-i*py outgoing
zP2(4)=-zP(4)                  ! px-i*py incoming -px+i*py outgoing
zP2(5)=zP(5)
!ELSE
!	zP2(1:5)=zP(1:5)
!ENDIF       
IF(nderi.EQ.0)THEN
	zqq(1)=DCMPLX(DREAL(zq1(1)))
	zqq(2)=DCMPLX(DREAL(zq1(2)))
	zqq(3)=zq1(3)
	zqq(4)=zq1(4)
ELSEIF(nderi.EQ.1)THEN
	IF(imode.EQ.1)THEN
		CALL PolarHadron(zqq,zP2,llz)
	ELSE
		SELECT CASE(llz)
		CASE(3)
			CALL Helac_elong(zqq,zP2)
		CASE(2)
			CALL Helac_eplus(zqq,zP2)
		CASE(1)
			CALL Helac_eminus(zqq,zP2)
		END SELECT
	ENDIF
	zqq(1:4)=-zqq(1:4)
	zmm=0
ELSE
	PRINT *,"Warning (nderi) in Helac_Qvplus ! STOP !"
	STOP
ENDIF
IF(spin.EQ.0)THEN
	zy1(1)=zqq(1)*zPy(3)+zqq(4)*zPy(4)-zmm*zPy(1)
	zy1(2)=zqq(3)*zPy(3)+zqq(2)*zPy(4)-zmm*zPy(2)
	zy1(3)=-zqq(2)*zPy(1)+zqq(4)*zPy(2)+zmm*zPy(3)
	zy1(4)=zqq(3)*zPy(1)-zqq(1)*zPy(2)+zmm*zPy(4)
	RETURN
ELSEIF(spin.EQ.1)THEN
	CALL PolarHadron(zeps,zP2,ssz)
	zy1(1)=(zeps(2)*zqq(1)-zeps(3)*zqq(4))*zPy(1)+(-zeps(4)*zqq(1)+zeps(1)*zqq(4))*zPy(2)&
	+zeps(1)*zmm*zPy(3)+zeps(4)*zmm*zPy(4)
	zy1(2)=(-zeps(3)*zqq(2)+zeps(2)*zqq(3))*zPy(1)+(zeps(1)*zqq(2)-zeps(4)*zqq(3))*zPy(2)&
	+zeps(3)*zmm*zPy(3)+zeps(2)*zmm*zPy(4)
	zy1(3)=zeps(2)*zmm*zPy(1)-zeps(4)*zmm*zPy(2)+(zeps(1)*zqq(2)-zeps(3)*zqq(4))*zPy(3)&
	+(zeps(4)*zqq(2)-zeps(2)*zqq(4))*zPy(4)
	zy1(4)=-zeps(3)*zmm*zPy(1)+zeps(1)*zmm*zPy(2)+(zeps(3)*zqq(1)-zeps(1)*zqq(3))*zPy(3)&
	+(zeps(2)*zqq(1)-zeps(4)*zqq(3))*zPy(4)
	RETURN
ENDIF
PRINT *,"Wrong in Helac_Qvplus ! STOP !"
STOP
END SUBROUTINE Helac_Qvplus 

SUBROUTINE  Helac_Qvminus(zy1,zq1,zP,spin,nderi,llz,ssz)
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(OUT)::zy1
COMPLEX(KIND=DBL),DIMENSION(5),INTENT(IN)::zq1,zP
INTEGER,INTENT(IN)::spin,nderi,llz,ssz  ! spin , rank of derivation, L_z,S_z
COMPLEX(KIND=DBL),DIMENSION(4)::zPy,zqq,zeps
COMPLEX(KIND=DBL),DIMENSION(5)::zP2
COMPLEX(KIND=DBL)::zmm
zmm=DCMPLX(DREAL(zq1(5)))
CALL Helac_uminus(zPy,zP)
! return from (p0,-p) to (p0,p) as the definition in 
! Eq.(A16) in Phys.Rev.D  57 (1998)  4258
!IF(PolarFrame.EQ.1)THEN
zP2(2)=DCMPLX(-DREAL(zP(2)),DIMAG(zP(2)))! p0-pz+i*pt  incoming -p0+pz+i*pt outgoing
zP2(1)=-zP(1)                  ! po+pz+i*pz incoming -po-pz-i*pz outgoing
zP2(3)=-zP(3)                  ! px+i*py incoming -px-i*py outgoing
zP2(4)=-zP(4)                  ! px-i*py incoming -px+i*py outgoing
zP2(5)=zP(5)
!ELSE
!	zP2(1:5)=zP(1:5)
!ENDIF          
IF(nderi.EQ.0)THEN
	zqq(1)=DCMPLX(DREAL(zq1(1)))
	zqq(2)=DCMPLX(DREAL(zq1(2)))
	zqq(3)=zq1(3)
	zqq(4)=zq1(4)
ELSEIF(nderi.EQ.1)THEN
	IF(imode.EQ.1)THEN        
		CALL PolarHadron(zqq,zP2,llz)
	ELSE
		SELECT CASE(llz)
		CASE(3)
			CALL Helac_elong(zqq,zP2)
		CASE(2)
			CALL Helac_eplus(zqq,zP2)
		CASE(1)
			CALL Helac_eminus(zqq,zP2)
		END SELECT
	ENDIF
	zqq(1:4)=-zqq(1:4)
	zmm=0
ELSE
	PRINT *,"Warning (nderi) in Helac_Qvminus ! STOP !"
	STOP
ENDIF
IF(spin.EQ.0)THEN
	zy1(1)=zqq(1)*zPy(3)+zqq(4)*zPy(4)-zmm*zPy(1)
	zy1(2)=zqq(3)*zPy(3)+zqq(2)*zPy(4)-zmm*zPy(2)
	zy1(3)=-zqq(2)*zPy(1)+zqq(4)*zPy(2)+zmm*zPy(3)
	zy1(4)=zqq(3)*zPy(1)-zqq(1)*zPy(2)+zmm*zPy(4)
	RETURN
ELSEIF(spin.EQ.1)THEN
	CALL PolarHadron(zeps,zP2,ssz)

	zy1(1)=(zeps(2)*zqq(1)-zeps(3)*zqq(4))*zPy(1)+(-zeps(4)*zqq(1)+zeps(1)*zqq(4))*zPy(2)&
	+zeps(1)*zmm*zPy(3)+zeps(4)*zmm*zPy(4)
	zy1(2)=(-zeps(3)*zqq(2)+zeps(2)*zqq(3))*zPy(1)+(zeps(1)*zqq(2)-zeps(4)*zqq(3))*zPy(2)&
	+zeps(3)*zmm*zPy(3)+zeps(2)*zmm*zPy(4)
	zy1(3)=zeps(2)*zmm*zPy(1)-zeps(4)*zmm*zPy(2)+(zeps(1)*zqq(2)-zeps(3)*zqq(4))*zPy(3)&
	+(zeps(4)*zqq(2)-zeps(2)*zqq(4))*zPy(4)
	zy1(4)=-zeps(3)*zmm*zPy(1)+zeps(1)*zmm*zPy(2)+(zeps(3)*zqq(1)-zeps(1)*zqq(3))*zPy(3)&
	+(zeps(2)*zqq(1)-zeps(4)*zqq(3))*zPy(4)
	RETURN
ENDIF
PRINT *,"Wrong in Helac_Qvminus ! STOP !"
STOP
END SUBROUTINE Helac_Qvminus   

SUBROUTINE  Helac_Qubplus(zy1,zq1,zP,nderi,llz)
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(OUT)::zy1
COMPLEX(KIND=DBL),DIMENSION(5),INTENT(IN)::zq1,zP
INTEGER,INTENT(IN)::nderi,llz  ! rank of derivation, L_z
COMPLEX(KIND=DBL),DIMENSION(4)::zPy,zqq,zeps
COMPLEX(KIND=DBL),DIMENSION(5)::zP2
COMPLEX(KIND=DBL)::zmm,zMM2,zmm3
zmm=DCMPLX(DREAL(zq1(5)))
zMM2=DCMPLX(DREAL(zP(5)))-zmm
zmm3=zmm
CALL Helac_ubplus(zPy,zP)
IF(nderi.EQ.0)THEN
	zqq(1)=DCMPLX(DREAL(zq1(1)))
	zqq(2)=DCMPLX(DREAL(zq1(2)))
	zqq(3)=zq1(3)
	zqq(4)=zq1(4)
ELSEIF(nderi.EQ.1)THEN
! return from (p0,-p) to (p0,p) as the definition in 
! Eq.(A16) in Phys.Rev.D  57 (1998)  4258
!	IF(PolarFrame.EQ.1)THEN
	zP2(2)=DCMPLX(-DREAL(zP(2)),DIMAG(zP(2)))! p0-pz+i*pt  incoming -p0+pz+i*pt outgoing
	zP2(1)=-zP(1)                  ! po+pz+i*pz incoming -po-pz-i*pz outgoing
	zP2(3)=-zP(3)                  ! px+i*py incoming -px-i*py outgoing
	zP2(4)=-zP(4)                  ! px-i*py incoming -px+i*py outgoing
	zP2(5)=zP(5)
!	ELSE
!		zP2(1:5)=zP(1:5)
!	ENDIF
	IF(imode.EQ.1)THEN                 
		CALL PolarHadron(zqq,zP2,llz)
	ELSE
		SELECT CASE(llz)
		CASE(3)
			CALL Helac_elong(zqq,zP2)
		CASE(2)
			CALL Helac_eplus(zqq,zP2)
		CASE(1)
			CALL Helac_eminus(zqq,zP2)
		END SELECT
	ENDIF
	zmm=0
ELSE
	PRINT *,"Warning (nderi) in Helac_Qubplus ! STOP !"
	STOP
ENDIF
zy1(1)=-zqq(2)*zPy(3)+zqq(3)*zPy(4)+zmm*zPy(1)
zy1(2)=zqq(4)*zPy(3)-zqq(1)*zPy(4)+zmm*zPy(2)
zy1(3)=-zqq(1)*zPy(1)-zqq(3)*zPy(2)+zmm*zPy(3)
zy1(4)=-zqq(4)*zPy(1)-zqq(2)*zPy(2)+zmm*zPy(4)
zy1(1:4)=zy1(1:4)/CDSQRT(32*zMM2*zmm3)/(zMM2+zmm3)
RETURN	
PRINT *,"Wrong in Helac_Qubplus ! STOP !"
STOP
END SUBROUTINE Helac_Qubplus 

SUBROUTINE  Helac_Qubminus(zy1,zq1,zP,nderi,llz)
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(OUT)::zy1
COMPLEX(KIND=DBL),DIMENSION(5),INTENT(IN)::zq1,zP
INTEGER,INTENT(IN)::nderi,llz  ! rank of derivation, L_z
COMPLEX(KIND=DBL),DIMENSION(4)::zPy,zqq,zeps
COMPLEX(KIND=DBL),DIMENSION(5)::zP2
COMPLEX(KIND=DBL)::zmm,zMM2,zmm3
zmm=DCMPLX(DREAL(zq1(5)))
zMM2=DCMPLX(DREAL(zP(5)))-zmm
zmm3=zmm
CALL Helac_ubminus(zPy,zP)
IF(nderi.EQ.0)THEN
	zqq(1)=DCMPLX(DREAL(zq1(1)))
	zqq(2)=DCMPLX(DREAL(zq1(2)))
	zqq(3)=zq1(3)
	zqq(4)=zq1(4)
ELSEIF(nderi.EQ.1)THEN
! return from (p0,-p) to (p0,p) as the definition in 
! Eq.(A16) in Phys.Rev.D  57 (1998)  4258
!	IF(PolarFrame.EQ.1)THEN
	zP2(2)=DCMPLX(-DREAL(zP(2)),DIMAG(zP(2)))! p0-pz+i*pt  incoming -p0+pz+i*pt outgoing
	zP2(1)=-zP(1)                  ! p0+pz+i*pz incoming -po-pz-i*pz outgoing
	zP2(3)=-zP(3)                  ! px+i*py incoming -px-i*py outgoing
	zP2(4)=-zP(4)                  ! px-i*py incoming -px+i*py outgoing
	zP2(5)=zP(5)
!	ELSE
!		zP2(1:5)=zP(1:5)
!	ENDIF
	IF(imode.EQ.1)THEN                  
		CALL PolarHadron(zqq,zP2,llz)
	ELSE
		SELECT CASE(llz)
		CASE(3)
			CALL Helac_elong(zqq,zP2)
		CASE(2)
			CALL Helac_eplus(zqq,zP2)
		CASE(1)
			CALL Helac_eminus(zqq,zP2)
		END SELECT
	ENDIF
	zmm=0
ELSE
	PRINT *,"Warning (nderi) in Helac_Qubminus ! STOP !"
	STOP
ENDIF
zy1(1)=-zqq(2)*zPy(3)+zqq(3)*zPy(4)+zmm*zPy(1)
zy1(2)=zqq(4)*zPy(3)-zqq(1)*zPy(4)+zmm*zPy(2)
zy1(3)=-zqq(1)*zPy(1)-zqq(3)*zPy(2)+zmm*zPy(3)
zy1(4)=-zqq(4)*zPy(1)-zqq(2)*zPy(2)+zmm*zPy(4)
zy1(1:4)=zy1(1:4)/CDSQRT(32*zMM2*zmm3)/(zMM2+zmm3)
RETURN	
PRINT *,"Wrong in Helac_Qubminus ! STOP !"
STOP
END SUBROUTINE Helac_Qubminus
! Eq.(A16) in Phys.Rev.D  57 (1998)  4258
SUBROUTINE PolarHadron(zy,zP,lamda)
IMPLICIT NONE
INTEGER,INTENT(IN)::lamda
COMPLEX(KIND=DBL),DIMENSION(5),INTENT(IN)::zP
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(OUT)::zy
COMPLEX(KIND=DBL)::alphaz,betaz,alphax,betax
COMPLEX(KIND=DBL),DIMENSION(4)::ZAA,ZBB,ZPP
COMPLEX(KIND=DBL),DIMENSION(4)::ZXX,ZYY,ZZZ,ztemp
COMPLEX(KIND=DBL)::zmass,zAP,zBP,Iuni,ziss
COMPLEX(KIND=DBL),DIMENSION(5)::zPP1,zPP2
! Helicity frame
IF(PolarFrame.EQ.1)THEN
	SELECT CASE(lamda)
	CASE(3)
		CALL Helac_elong(zy,zP)
	CASE(2)
		CALL Helac_eplus(zy,zP)
	CASE(1)
		CALL Helac_eminus(zy,zP)
	CASE DEFAULT
		PRINT *,"Wrong (lamda, helicity) in PolarHadron ! STOP !"
		STOP
	END SELECT
	RETURN
ELSE
	zPP1(1:5)=zq(1,1:5)/xp1
	zPP2(1:5)=zq(2,1:5)/xp2
	Iuni=DCMPLX(0,1d0)
	ZAA(1)=(zPP1(3)+zPP2(3)+zPP1(4)+zPP2(4))/DCMPLX(2d0)
	ZAA(2)=(zPP1(3)+zPP2(3)-zPP1(4)-zPP2(4))/DCMPLX(2d0)/Iuni
	ZAA(3)=DCMPLX(DIMAG(zPP1(1))+DIMAG(zPP2(1)))
	ZAA(4)=DCMPLX(DREAL(zPP1(1))+DREAL(zPP2(1))+DREAL(zPP1(2))+DREAL(zPP2(2)))/DCMPLX(2d0)
	ZBB(1)=(zPP1(3)-zPP2(3)+zPP1(4)-zPP2(4))/DCMPLX(2d0)
	ZBB(2)=(zPP1(3)-zPP2(3)-zPP1(4)+zPP2(4))/DCMPLX(2d0)/Iuni
	ZBB(3)=DCMPLX(DIMAG(zPP1(1))-DIMAG(zPP2(1)))
	ZBB(4)=DCMPLX(DREAL(zPP1(1))-DREAL(zPP2(1))+DREAL(zPP1(2))-DREAL(zPP2(2)))/DCMPLX(2d0)
	ZPP(1)=(zP(3)+zP(4))/DCMPLX(2d0)
	ZPP(2)=(zP(3)-zP(4))/DCMPLX(2d0)/Iuni
	ZPP(3)=DCMPLX(DIMAG(zP(1)))
	ZPP(4)=DCMPLX(DREAL(zP(1))+DREAL(zP(2)))/DCMPLX(2d0)
	zmass=DCMPLX(DREAL(zP(5)))
	zAP=ZPP(4)*ZAA(4)-ZPP(1)*ZAA(1)-ZPP(2)*ZAA(2)-ZPP(3)*ZAA(3)
	zBP=ZPP(4)*ZBB(4)-ZPP(1)*ZBB(1)-ZPP(2)*ZBB(2)-ZPP(3)*ZBB(3)
	ziss=ZAA(4)**2-ZAA(1)**2-ZAA(2)**2-ZAA(3)**2
	SELECT CASE(PolarFrame)
!	CASE(1)   	! Helicity frame
!		alphaz=-zmass/CDSQRT(zAP**2-zmass*2*ziss)
!		betaz=DCMPLX(0d0)
!		alphax=(zAP*zBP)/CDSQRT(ziss*(zAP**2-zmass**2*ziss)*(zAP**2-zBP**2-zmass**2*ziss))
!		betax=-CDSQRT(zAP**2-zmass**2*ziss)/CDSQRT(ziss*(zAP**2-zBP**2-zmass**2*ziss))
	CASE(2)     ! Collins-Soper frame
		alphaz=-zBP/CDSQRT(ziss*(zAP**2-zBP**2))
		betaz=zAP/CDSQRT(ziss*(zAP**2-zBP**2))
		alphax=-zmass*zAP/CDSQRT((zAP**2-zBP**2)*(zAP**2-zBP**2-zmass**2*ziss))
		betax=zmass*zBP/CDSQRT((zAP**2-zBP**2)*(zAP**2-zBP**2-zmass**2*ziss))
	CASE(3)     ! Gottfried-Jackson frame
		alphaz=zmass/(zAP+zBP)
		betaz=alphaz
		alphax=-(zBP**2+zAP*zBP+zmass**2*ziss)/(zAP+zBP)/&
		CDSQRT(ziss*(zAP**2-zBP**2-zmass**2*ziss))
		betax=(zAP**2+zAP*zBP-zmass**2*ziss)/(zAP+zBP)/&
		CDSQRT(ziss*(zAP**2-zBP**2-zmass**2*ziss))
	CASE(4)     ! Target frame
		alphaz=-zmass/(zAP-zBP)
		betaz=-alphaz
		alphax=-(zBP**2-zAP*zBP+zmass**2*ziss)/(zAP-zBP)/&
		CDSQRT(ziss*(zAP**2-zBP**2-zmass**2*ziss))
		betax=-(zAP**2-zAP*zBP-zmass**2*ziss)/(zAP-zBP)/&
		CDSQRT(ziss*(zAP**2-zBP**2-zmass**2*ziss))
	CASE DEFAULT
		PRINT *,"Wrong (PolarFrame) in PolarHadron ! STOP !"
		STOP
	END SELECT

	ZZZ(1:4)=alphaz*(zAA(1:4)-zAP/zmass**2*(zPP(1:4)))&
	+betaz*(zBB(1:4)-zBP/zmass**2*(zPP(1:4)))

	ZXX(1:4)=alphax*(zAA(1:4)-zAP/zmass**2*(zPP(1:4)))&
	+betax*(zBB(1:4)-zBP/zmass**2*(zPP(1:4)))
	
	ZYY(4)=(-zPP(1)*ZXX(2)*ZZZ(3)+zPP(1)*ZXX(3)*ZZZ(2)&
	-zPP(3)*ZXX(1)*ZZZ(2)+zPP(3)*ZXX(2)*ZZZ(1)&
	-zPP(2)*ZXX(3)*ZZZ(1)+zPP(2)*ZXX(1)*ZZZ(3))/zmass
	ZYY(1)=(-zPP(4)*zXX(2)*ZZZ(3)+zPP(4)*ZXX(3)*ZZZ(2)&
	-zPP(3)*ZXX(4)*ZZZ(2)+zPP(3)*ZXX(2)*ZZZ(4)&
	-zPP(2)*ZXX(3)*ZZZ(4)+zPP(2)*ZXX(4)*ZZZ(3))/zmass
	ZYY(2)=(zPP(4)*ZXX(1)*ZZZ(3)-zPP(4)*ZXX(3)*ZZZ(1)&
	+zPP(3)*ZXX(4)*ZZZ(1)-zPP(3)*ZXX(1)*ZZZ(4)&
	+zPP(1)*ZXX(3)*ZZZ(4)-zPP(1)*ZXX(4)*ZZZ(3))/zmass
	ZYY(3)=(-zPP(4)*ZXX(1)*ZZZ(2)+zPP(4)*ZXX(2)*ZZZ(1)&
	-zPP(2)*ZXX(4)*ZZZ(1)+zPP(2)*ZXX(1)*ZZZ(4)&
	-zPP(1)*ZXX(2)*ZZZ(4)+zPP(1)*ZXX(4)*ZZZ(2))/zmass

	SELECT CASE(lamda)
	CASE(3)
		zy(1:4)=ZZZ(1:4)
	CASE(2)
		zy(1:4)=(-ZXX(1:4)-Iuni*ZYY(1:4))/DCMPLX(DSQRT(2d0))
	CASE(1)
		zy(1:4)=(ZXX(1:4)-Iuni*ZYY(1:4))/DCMPLX(DSQRT(2d0))
	CASE DEFAULT
		PRINT *,"Wrong(lamda) in PolarHadron ! STOP !"
		STOP
	END SELECT

	CALL transf_lor_lc(ztemp,zy)
	zy(1:4)=ztemp(1:4)
	RETURN
ENDIF
END SUBROUTINE PolarHadron
! transfer lorentz representation p=(px,py,pz,p0) 
! to light cone representation pA=(p0+pz,p0-pz,px+i*py,px-i*py)
SUBROUTINE transf_lor_lc(zpA1,zp1)
IMPLICIT NONE
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(IN)::zp1
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(OUT)::zpA1
COMPLEX(KIND=DBL),PARAMETER::Iuni=DCMPLX(0,1d0)
zpA1(1)=zp1(4)+zp1(3)
zpA1(2)=zp1(4)-zp1(3)
zpA1(3)=zp1(1)+Iuni*zp1(2)
zpA1(4)=zp1(1)-Iuni*zp1(2)
END SUBROUTINE transf_lor_lc
! return whether two color i,j are singlet
FUNCTION ColorSinglet(col1,col2)
IMPLICIT NONE
INTEGER,INTENT(IN)::col1,col2
LOGICAL::ColorSinglet
INTEGER::i,j,init=0
INTEGER,DIMENSION(1:20)::colorlist
INTEGER,DIMENSION(10)::octetlist
SAVE init,colorlist
IF(init.EQ.0)THEN
	colorlist(1:20)=0
	j=1
	DO i=1,nhad
		IF(i.LE.2)THEN
			IF(iflh(i).EQ.-3.OR.iflh(i).EQ.-4.OR.iflh(i).EQ.-7& 
			.OR.iflh(i).EQ.-8.OR.iflh(i).EQ.-11.OR.iflh(i).EQ.-12.OR.iflh(i).EQ.35)j=j+1
		ELSE
			IF(iflh(i).EQ.3.OR.iflh(i).EQ.4.OR.iflh(i).EQ.7& 
			.OR.iflh(i).EQ.8.OR.iflh(i).EQ.11.OR.iflh(i).EQ.12.OR.iflh(i).EQ.35)THEN
				j=j+1
			ELSEIF(ABS(iflh(i)).GT.100)THEN
				colorlist(j)=i
				j=j+1
			ENDIF
		ENDIF
	ENDDO
	init=1
ENDIF
IF(col1.NE.col2)THEN
	ColorSinglet=.FALSE.
	RETURN
ENDIF
IF(colorlist(col1).EQ.0)THEN
	ColorSinglet=.FALSE.
	RETURN
ENDIF
octetlist(1:10)=0
CALL Helac_bin(nhad,octetmsinglet,octetlist(1:nhad))
IF(SubInteger(iflh(colorlist(col1)),1,1).EQ.1.OR.octetlist(nhad-colorlist(col1)+1).EQ.1)THEN
	ColorSinglet=.TRUE.
	RETURN
ENDIF
!j=1
!DO i=1,nhad
!	IF(i.LE.2)THEN
!		IF(iflh(i).EQ.-3.OR.iflh(i).EQ.-4.OR.iflh(i).EQ.-7& 
!        .OR.iflh(i).EQ.-8.OR.iflh(i).EQ.-11.OR.iflh(i).EQ.-12.OR.iflh(i).EQ.35)j=j+1
!	ELSE
!		IF(iflh(i).EQ.3.OR.iflh(i).EQ.4.OR.iflh(i).EQ.7& 
!        .OR.iflh(i).EQ.8.OR.iflh(i).EQ.11.OR.iflh(i).EQ.12.OR.iflh(i).EQ.35)THEN
!			j=j+1
!		ELSEIF(ABS(iflh(i)).GT.100)THEN
!			IF(SubInteger(iflh(i),1,1).EQ.1.AND.col1.EQ.j)THEN
!				ColorSinglet=.TRUE.
!				RETURN
!			ELSEIF(SubInteger(iflh(i),1,1).NE.1.AND.col1.EQ.j)THEN
!				ColorSinglet=.FALSE.
!				RETURN
!			ENDIF
!			j=j+1
!		ENDIF
!	ENDIF
!	IF(col1.LT.j)THEN
!		ColorSinglet=.FALSE.
!		RETURN
!	ENDIF
!ENDDO
ColorSinglet=.FALSE.
RETURN
END FUNCTION ColorSinglet

FUNCTION ColorSinglet2(col1,col2)
IMPLICIT NONE
INTEGER,INTENT(IN)::col1,col2
LOGICAL::ColorSinglet2
INTEGER::i,j,init=0
INTEGER,DIMENSION(1:20)::colorlist
SAVE init,colorlist
IF(init.EQ.0)THEN
	colorlist(1:20)=0
	j=1
	DO i=1,nhad
		IF(i.LE.2)THEN
			IF(iflh(i).EQ.-3.OR.iflh(i).EQ.-4.OR.iflh(i).EQ.-7& 
			.OR.iflh(i).EQ.-8.OR.iflh(i).EQ.-11.OR.iflh(i).EQ.-12.OR.iflh(i).EQ.35)j=j+1
		ELSE
			IF(iflh(i).EQ.3.OR.iflh(i).EQ.4.OR.iflh(i).EQ.7& 
			.OR.iflh(i).EQ.8.OR.iflh(i).EQ.11.OR.iflh(i).EQ.12.OR.iflh(i).EQ.35)THEN
				j=j+1
			ELSEIF(ABS(iflh(i)).GT.100)THEN
				colorlist(j)=i
				j=j+1
			ENDIF
		ENDIF
	ENDDO
	init=1
ENDIF
IF(col1.LE.0.OR.col2.LE.0)THEN
	ColorSinglet2=.FALSE.
	RETURN
ENDIF
IF(colorlist(col2).EQ.0)THEN
	ColorSinglet2=.FALSE.
	RETURN
ENDIF
i=Quarkonium2(colorlist(col2))
i=Quarkonium(i,3)
IF(SubInteger(iflh(colorlist(col2)),1,1).EQ.1.AND.icol(i,1).EQ.col1)THEN
	ColorSinglet2=.TRUE.
	RETURN
ENDIF
ColorSinglet2=.FALSE.
RETURN
END FUNCTION ColorSinglet2
! return the quark position in color array of the related antiquark positionin color array
FUNCTION QuarkColorPos(k)
IMPLICIT NONE
INTEGER,INTENT(IN)::k
INTEGER::QuarkColorPos
INTEGER::i,j1,j2,init=0
INTEGER,DIMENSION(1:20)::colorlist
SAVE init,colorlist
IF(init.EQ.0)THEN
	j1=0
	j2=0
	colorlist(1:20)=0
	DO i=1,nhad
		IF(i.LE.2)THEN
			IF(iflh(i).EQ.-3.OR.iflh(i).EQ.-4.OR.iflh(i).EQ.-7& 
			.OR.iflh(i).EQ.-8.OR.iflh(i).EQ.-11.OR.iflh(i).EQ.-12.OR.iflh(i).EQ.35)j2=j2+1
			IF(iflh(i).EQ.3.OR.iflh(i).EQ.4.OR.iflh(i).EQ.7& 
			.OR.iflh(i).EQ.8.OR.iflh(i).EQ.11.OR.iflh(i).EQ.12.OR.iflh(i).EQ.35)j1=j1+1
		ELSE
			IF(iflh(i).EQ.-3.OR.iflh(i).EQ.-4.OR.iflh(i).EQ.-7& 
			.OR.iflh(i).EQ.-8.OR.iflh(i).EQ.-11.OR.iflh(i).EQ.-12.OR.iflh(i).EQ.35)j1=j1+1
			IF(iflh(i).EQ.3.OR.iflh(i).EQ.4.OR.iflh(i).EQ.7& 
			.OR.iflh(i).EQ.8.OR.iflh(i).EQ.11.OR.iflh(i).EQ.12.OR.iflh(i).EQ.35)j2=j2+1
			IF(AbS(iflh(i)).GT.100)THEN
				j1=j1+1
				j2=j2+1
				colorlist(j2)=j1
			ENDIF		
		ENDIF
	ENDDO
	init=1
ENDIF
QuarkColorPos=colorlist(k)
RETURN
END FUNCTION QuarkColorPos
END MODULE Helac_wavef
