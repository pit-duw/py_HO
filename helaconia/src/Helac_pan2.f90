MODULE Helac_pan2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  This module includes  all vertex and propagator functions !
!  Multiply a factor -i to the feynman rules in arXiv:0709.1075[hep-ph] 
!  and Phys.Rev.D 67,014026(2003)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
USE Helac_Global
USE Helac_Func_1
USE Helac_wavef
IMPLICIT NONE
!COMPLEX(KIND=DBL)::fmas
!COMPLEX(KIND=DBL),DIMENSION(4)::zb0,zb1,zb2,zb3,zpb0,zpb1,zpb2,zps0,zps2,zp0,zp1,zp2
CONTAINS
!SUBROUTINE Helac_vertices
!       include 'declare.h'       
!       include 'common_flags.h'
!       include 'common_masses.h'
!END SUBROUTINE Helac_vertices
SUBROUTINE Helac_xv3(zb0,zp0,zb1,zp1,zb2,zp2,mm,zg,isg)
IMPLICIT NONE
INTEGER,INTENT(IN)::isg,mm
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(OUT)::zb0
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(IN)::zp0,zb1,zp1,zp2,zb2
COMPLEX(KIND=DBL),INTENT(IN)::zg
COMPLEX(KIND=DBL)::z,z1,z2,z3,z4,z5 
INTEGER::k 

z=zg*isg
z1=Helac_zprod(zp1,zb2)*z
z2=Helac_zprod(zp0,zb2)*z
z3=Helac_zprod(zp2,zb1)*z
z4=Helac_zprod(zp0,zb1)*z
z5=Helac_zprod(zb1,zb2)*z
       
DO k=1,4
   zb0(k)=(z1+z2)*zb1(k)-(z3+z4)*zb2(k)+z5*( zp2(k) - zp1(k) ) 
ENDDO

END SUBROUTINE Helac_xv3

SUBROUTINE Helac_xv3d(zb0,zp0,zb1,zp1,zb2,zp2,mm,zg,isg,deriv,ip1array,ip2array)
IMPLICIT NONE
INTEGER,INTENT(IN)::isg,mm
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(OUT)::zb0
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(IN)::zp0,zb1,zp1,zp2,zb2
COMPLEX(KIND=DBL),INTENT(IN)::zg
INTEGER,INTENT(IN)::deriv   ! 2**(i-1)+2**(j-1) represents derive to i-th and j-th p-waves
INTEGER,DIMENSION(10),INTENT(IN)::ip1array,ip2array
COMPLEX(KIND=DBL),DIMENSION(4)::zzp0,zzp1,zzp2,zeps
INTEGER::kk1,k

IF(deriv.LE.0)THEN
	CALL Helac_xv3(zb0,zp0,zb1,zp1,zb2,zp2,mm,zg,isg)
	RETURN
ENDIF
IF(Helac_level(Pwavenum,deriv).GT.1)THEN
	zb0(1:4)=0
	RETURN
ENDIF
CALL Deriv_eps(zeps,k,deriv)
DO kk1=1,4
	zzp0(kk1)=(ip1array(k)+ip2array(k))*zeps(kk1)
	zzp1(kk1)=ip1array(k)*zeps(kk1)
	zzp2(kk1)=ip2array(k)*zeps(kk1)
ENDDO
CALL Helac_xv3(zb0,zzp0,zb1,zzp1,zb2,zzp2,mm,zg,isg)

RETURN

END SUBROUTINE Helac_xv3d

SUBROUTINE Helac_xvff(zb0,zp0,zpb1,zps2,mm,zgl,zgr,isg)
IMPLICIT NONE
INTEGER,INTENT(IN)::isg,mm
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(OUT)::zb0
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(IN)::zp0,zpb1,zps2
COMPLEX(KIND=DBL),INTENT(IN)::zgl,zgr
COMPLEX(KIND=DBL)::z,zl,zr
z=dnou(2)*isg
zl=zgl*z
zr=zgr*z
! in light-cone expression see Eq(16) in hep-ph/0002082
zb0(1)=-zps2(1)*zpb1(3)*zr-zps2(4)*zpb1(2)*zl
     
zb0(2)=-zps2(2)*zpb1(4)*zr-zps2(3)*zpb1(1)*zl
     
zb0(3)=-zps2(2)*zpb1(3)*zr+zps2(4)*zpb1(1)*zl
     
zb0(4)=-zps2(1)*zpb1(4)*zr+zps2(3)*zpb1(2)*zl
      
END SUBROUTINE Helac_xvff
           
SUBROUTINE Helac_xffvr(zps0,zp0,zb1,zps2,mm,zgr,isg,fmas)
IMPLICIT NONE
INTEGER,INTENT(IN)::isg,mm
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(OUT)::zps0
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(IN)::zp0,zb1,zps2
COMPLEX(KIND=DBL),INTENT(IN)::zgr,fmas
COMPLEX(KIND=DBL)::z,zz
zps0(3)=DCMPLX(dnou(0),dnou(0))
zps0(4)=DCMPLX(dnou(0),dnou(0))
z=isg*zgr
! in light-cone expression see Eq(18) in hep-ph/0002082         
zps0(1)=((-zb1(2)*zp0(1)+zb1(3)*zp0(4))*zps2(1)&
 +( zb1(4)*zp0(1)-zb1(1)*zp0(4))*zps2(2))*z

zps0(2)=(( zb1(3)*zp0(2)-zb1(2)*zp0(3))*zps2(1)&
+(-zb1(1)*zp0(2)+zb1(4)*zp0(3))*zps2(2))*z
           
IF(fmas.NE.dnou(0))THEN
   zz=z*fmas
   zps0(3)=( zb1(2)*zps2(1)-zb1(4)*zps2(2))*zz
   zps0(4)=(-zb1(3)*zps2(1)+zb1(1)*zps2(2))*zz
ENDIF
        
END SUBROUTINE Helac_xffvr

SUBROUTINE Helac_xffvrd(zps0,zp0,zb1,zps2,mm,zgr,isg,fmas,deriv,ip0array)
IMPLICIT NONE
INTEGER,INTENT(IN)::isg,mm
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(OUT)::zps0
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(IN)::zp0,zb1,zps2
COMPLEX(KIND=DBL),INTENT(IN)::zgr,fmas
INTEGER,DIMENSION(10),INTENT(IN)::ip0array
COMPLEX(KIND=DBL),DIMENSION(4)::zeps
INTEGER,INTENT(IN)::deriv
INTEGER::k

IF(deriv.LE.0)THEN
	CALL Helac_xffvr(zps0,zp0,zb1,zps2,mm,zgr,isg,fmas)
	RETURN
ENDIF
IF(Helac_level(Pwavenum,deriv).GT.1)THEN
	zps0(1:4)=0
	RETURN
ENDIF
CALL Deriv_eps(zeps,k,deriv)
IF(ip0array(k).EQ.0)THEN
	zps0(1:4)=0
	RETURN
ENDIF
zeps(1:4)=ip0array(k)*zeps(1:4)
CALL Helac_xffvr(zps0,zeps,zb1,zps2,mm,zgr,isg,DCMPLX(0d0))
        
END SUBROUTINE Helac_xffvrd

SUBROUTINE Helac_xffvl(zps0,zp0,zb1,zps2,mm,zgl,isg,fmas)
IMPLICIT NONE
INTEGER,INTENT(IN)::isg,mm
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(OUT)::zps0
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(IN)::zp0,zb1,zps2
COMPLEX(KIND=DBL),INTENT(IN)::zgl,fmas
COMPLEX(KIND=DBL)::z,zz
zps0(1)=DCMPLX(dnou(0),dnou(0))
zps0(2)=DCMPLX(dnou(0),dnou(0))
z=isg*zgl

zps0(3)=((-zb1(1)*zp0(2)+zb1(3)*zp0(4))*zps2(3)&
  +(-zb1(4)*zp0(2)+zb1(2)*zp0(4))*zps2(4))*z

zps0(4)=((-zb1(3)*zp0(1)+zb1(1)*zp0(3))*zps2(3)&
  +(-zb1(2)*zp0(1)+zb1(4)*zp0(3))*zps2(4))*z
           
IF(fmas.NE.dnou(0))THEN
   zz=z*fmas
   zps0(1)=( zb1(1)*zps2(3)+zb1(4)*zps2(4))*zz
   zps0(2)=( zb1(3)*zps2(3)+zb1(2)*zps2(4))*zz
ENDIF
END SUBROUTINE Helac_xffvl

SUBROUTINE Helac_xffvld(zps0,zp0,zb1,zps2,mm,zgl,isg,fmas,deriv,ip0array)
IMPLICIT NONE
INTEGER,INTENT(IN)::isg,mm
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(OUT)::zps0
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(IN)::zp0,zb1,zps2
COMPLEX(KIND=DBL),INTENT(IN)::zgl,fmas
INTEGER,DIMENSION(10),INTENT(IN)::ip0array
COMPLEX(KIND=DBL),DIMENSION(4)::zeps
INTEGER,INTENT(IN)::deriv
INTEGER::k

IF(deriv.LE.0)THEN
	CALL Helac_xffvl(zps0,zp0,zb1,zps2,mm,zgl,isg,fmas)
	RETURN
ENDIF
IF(Helac_level(Pwavenum,deriv).GT.1)THEN
	zps0(1:4)=0
	RETURN
ENDIF
CALL Deriv_eps(zeps,k,deriv)
IF(ip0array(k).EQ.0)THEN
	zps0(1:4)=0
	RETURN
ENDIF
zeps(1:4)=ip0array(k)*zeps(1:4)
CALL Helac_xffvl(zps0,zeps,zb1,zps2,mm,zgl,isg,DCMPLX(0d0))
        
END SUBROUTINE Helac_xffvld

SUBROUTINE Helac_xffvl0(zps0,zb1,zps2,zgl,isg)
IMPLICIT NONE
INTEGER,INTENT(IN)::isg
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(OUT)::zps0
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(IN)::zb1,zps2
COMPLEX(KIND=DBL),INTENT(IN)::zgl
COMPLEX(KIND=DBL)::z       
zps0(3)=DCMPLX(dnou(0),dnou(0))
zps0(4)=DCMPLX(dnou(0),dnou(0))
       
z=isg*zgl

zps0(1)=( zb1(1)*zps2(3)+zb1(4)*zps2(4))*z

zps0(2)=( zb1(3)*zps2(3)+zb1(2)*zps2(4))*z
        
END SUBROUTINE Helac_xffvl0
 
SUBROUTINE Helac_xffvr0(zps0,zb1,zps2,zgr,isg)
IMPLICIT NONE
INTEGER,INTENT(IN)::isg
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(OUT)::zps0
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(IN)::zb1,zps2
COMPLEX(KIND=DBL),INTENT(IN)::zgr
COMPLEX(KIND=DBL)::z           
zps0(1)=DCMPLX(dnou(0),dnou(0))
zps0(2)=DCMPLX(dnou(0),dnou(0))
z=isg*zgr
zps0(3)=(zb1(2)*zps2(1)-zb1(4)*zps2(2))*z
zps0(4)=(-zb1(3)*zps2(1)+zb1(1)*zps2(2))*z
END SUBROUTINE Helac_xffvr0

SUBROUTINE Helac_xaffvl(zpb0,zp0,zb1,zpb2,mm,zgl,isg,fmas)
IMPLICIT NONE
INTEGER,INTENT(IN)::isg,mm
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(OUT)::zpb0
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(IN)::zp0,zb1,zpb2
COMPLEX(KIND=DBL),INTENT(IN)::zgl,fmas
COMPLEX(KIND=DBL)::z,zz
zpb0(3)=DCMPLX(dnou(0),dnou(0))
zpb0(4)=DCMPLX(dnou(0),dnou(0))
z=isg*zgl
zpb0(1)= (( zb1(1)*zp0(2)-zb1(4)*zp0(3))*zpb2(1)&
+( zb1(3)*zp0(2)-zb1(2)*zp0(3))*zpb2(2))*z
zpb0(2)= (( zb1(4)*zp0(1)-zb1(1)*zp0(4))*zpb2(1)&
+( zb1(2)*zp0(1)-zb1(3)*zp0(4))*zpb2(2))*z
IF(fmas.NE.dnou(0))THEN
   zz=z*fmas
   zpb0(3)=( zb1(1)*zpb2(1)+zb1(3)*zpb2(2))*zz
   zpb0(4)=( zb1(4)*zpb2(1)+zb1(2)*zpb2(2))*zz
ENDIF
END SUBROUTINE Helac_xaffvl

SUBROUTINE Helac_xaffvld(zpb0,zp0,zb1,zpb2,mm,zgl,isg,fmas,deriv,ip0array)
IMPLICIT NONE
INTEGER,INTENT(IN)::isg,mm
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(OUT)::zpb0
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(IN)::zp0,zb1,zpb2
COMPLEX(KIND=DBL),INTENT(IN)::zgl,fmas
INTEGER,DIMENSION(10),INTENT(IN)::ip0array
COMPLEX(KIND=DBL),DIMENSION(4)::zeps
INTEGER,INTENT(IN)::deriv
INTEGER::k

IF(deriv.LE.0)THEN
	CALL Helac_xaffvl(zpb0,zp0,zb1,zpb2,mm,zgl,isg,fmas)
	RETURN
ENDIF
IF(Helac_level(Pwavenum,deriv).GT.1)THEN
	zpb0(1:4)=0
	RETURN
ENDIF
CALL Deriv_eps(zeps,k,deriv)
IF(ip0array(k).EQ.0)THEN
	zpb0(1:4)=0
	RETURN
ENDIF
zeps(1:4)=ip0array(k)*zeps(1:4)
CALL Helac_xaffvl(zpb0,zeps,zb1,zpb2,mm,zgl,isg,DCMPLX(0d0))
        
END SUBROUTINE Helac_xaffvld
      
SUBROUTINE Helac_xaffvr(zpb0,zp0,zb1,zpb2,mm,zgr,isg,fmas)
IMPLICIT NONE
INTEGER,INTENT(IN)::isg,mm
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(OUT)::zpb0
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(IN)::zp0,zb1,zpb2
COMPLEX(KIND=DBL),INTENT(IN)::zgr,fmas
COMPLEX(KIND=DBL)::z,zz
zpb0(1)=DCMPLX(dnou(0),dnou(0))
zpb0(2)=DCMPLX(dnou(0),dnou(0))
z=isg*zgr
zpb0(3)= (( zb1(2)*zp0(1)-zb1(4)*zp0(3))*zpb2(3)&
+( zb1(1)*zp0(3)-zb1(3)*zp0(1))*zpb2(4))*z
zpb0(4)= (( zb1(2)*zp0(4)-zb1(4)*zp0(2))*zpb2(3)&
+( zb1(1)*zp0(2)-zb1(3)*zp0(4))*zpb2(4))*z
         
IF(fmas.NE.dnou(0))THEN
   zz=z*fmas       
   zpb0(1)=( zb1(2)*zpb2(3)-zb1(3)*zpb2(4))*zz
   zpb0(2)=(-zb1(4)*zpb2(3)+zb1(1)*zpb2(4))*zz
ENDIF
END SUBROUTINE Helac_xaffvr

SUBROUTINE Helac_xaffvrd(zpb0,zp0,zb1,zpb2,mm,zgr,isg,fmas,deriv,ip0array)
IMPLICIT NONE
INTEGER,INTENT(IN)::isg,mm
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(OUT)::zpb0
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(IN)::zp0,zb1,zpb2
COMPLEX(KIND=DBL),INTENT(IN)::zgr,fmas
INTEGER,DIMENSION(10),INTENT(IN)::ip0array
COMPLEX(KIND=DBL),DIMENSION(4)::zeps
INTEGER,INTENT(IN)::deriv
INTEGER::k

IF(deriv.LE.0)THEN
	CALL Helac_xaffvr(zpb0,zp0,zb1,zpb2,mm,zgr,isg,fmas)
	RETURN
ENDIF
IF(Helac_level(Pwavenum,deriv).GT.1)THEN
	zpb0(1:4)=0
	RETURN
ENDIF
CALL Deriv_eps(zeps,k,deriv)
IF(ip0array(k).EQ.0)THEN
	zpb0(1:4)=0
	RETURN
ENDIF
zeps(1:4)=ip0array(k)*zeps(1:4)
CALL Helac_xaffvr(zpb0,zeps,zb1,zpb2,mm,zgr,isg,DCMPLX(0d0))
        
END SUBROUTINE Helac_xaffvrd
        
SUBROUTINE Helac_xaffvr0(zpb0,zb1,zpb2,zgr,isg)
IMPLICIT NONE
INTEGER,INTENT(IN)::isg
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(OUT)::zpb0
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(IN)::zb1,zpb2
COMPLEX(KIND=DBL),INTENT(IN)::zgr
COMPLEX(KIND=DBL)::z
zpb0(3)=DCMPLX(dnou(0),dnou(0))
zpb0(4)=DCMPLX(dnou(0),dnou(0))
z=isg*zgr     
zpb0(1)=( zb1(2)*zpb2(3)-zb1(3)*zpb2(4))*z  
zpb0(2)=(-zb1(4)*zpb2(3)+zb1(1)*zpb2(4))*z
END SUBROUTINE Helac_xaffvr0
        
SUBROUTINE Helac_xaffvl0(zpb0,zb1,zpb2,zgl,isg)
IMPLICIT NONE
INTEGER,INTENT(IN)::isg
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(OUT)::zpb0
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(IN)::zb1,zpb2
COMPLEX(KIND=DBL),INTENT(IN)::zgl
COMPLEX(KIND=DBL)::z
zpb0(1)=DCMPLX(dnou(0),dnou(0))
zpb0(2)=DCMPLX(dnou(0),dnou(0))
z=isg*zgl   
zpb0(3)=( zb1(1)*zpb2(1)+zb1(3)*zpb2(2))*z
zpb0(4)=( zb1(4)*zpb2(1)+zb1(2)*zpb2(2))*z
END SUBROUTINE Helac_xaffvl0
         
SUBROUTINE Helac_xv4(zb0,zp0,zb1,zb2,zb3,mm,zg,isg)
IMPLICIT NONE
INTEGER,INTENT(IN)::isg,mm
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(OUT)::zb0
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(IN)::zp0,zb1,zb2,zb3
COMPLEX(KIND=DBL),INTENT(IN)::zg
COMPLEX(KIND=DBL)::z
INTEGER::k     
z=zg*isg
DO k=1,4
   zb0(k)=(dnou(2)*Helac_zprod(zb2,zb3)*zb1(k)-Helac_zprod(zb1,zb3)*zb2(k)&
          -Helac_zprod(zb1,zb2)*zb3(k))*z
ENDDO
END SUBROUTINE Helac_xv4
      
SUBROUTINE Helac_xffsl(zps0,zp0,zs1,zps2,mm,zgl,isg,fmas)
IMPLICIT NONE
INTEGER,INTENT(IN)::isg,mm
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(OUT)::zps0
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(IN)::zp0,zps2
COMPLEX(KIND=DBL),INTENT(IN)::zgl,fmas,zs1
COMPLEX(KIND=DBL)::z               
zps0(3)=DCMPLX(dnou(0),dnou(0))
zps0(4)=DCMPLX(dnou(0),dnou(0))
z=isg*zgl
zps0(1)=(zp0(1)*zps2(3)+zp0(4)*zps2(4))*zs1*z
zps0(2)=(zp0(3)*zps2(3)+zp0(2)*zps2(4))*zs1*z
IF(fmas.NE.dnou(0))THEN
    zps0(3)=-zs1*zps2(3)*fmas*z            
    zps0(4)=-zs1*zps2(4)*fmas*z
ENDIF
END SUBROUTINE Helac_xffsl

SUBROUTINE Helac_xffsld(zps0,zp0,zs1,zps2,mm,zgl,isg,fmas,deriv,ip0array)
IMPLICIT NONE
INTEGER,INTENT(IN)::isg,mm
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(OUT)::zps0
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(IN)::zp0,zps2
COMPLEX(KIND=DBL),INTENT(IN)::zgl,fmas,zs1
INTEGER,DIMENSION(10),INTENT(IN)::ip0array
COMPLEX(KIND=DBL),DIMENSION(4)::zeps
INTEGER,INTENT(IN)::deriv
INTEGER::k

IF(deriv.LE.0)THEN
	CALL Helac_xffsl(zps0,zp0,zs1,zps2,mm,zgl,isg,fmas)
	RETURN
ENDIF
IF(Helac_level(Pwavenum,deriv).GT.1)THEN
	zps0(1:4)=0
	RETURN
ENDIF
CALL Deriv_eps(zeps,k,deriv)
IF(ip0array(k).EQ.0)THEN
	zps0(1:4)=0
	RETURN
ENDIF
zeps(1:4)=ip0array(k)*zeps(1:4)
CALL Helac_xffsl(zps0,zeps,zs1,zps2,mm,zgl,isg,DCMPLX(0d0))
        
END SUBROUTINE Helac_xffsld
      
SUBROUTINE Helac_xffsr(zps0,zp0,zs1,zps2,mm,zgr,isg,fmas)
IMPLICIT NONE
INTEGER,INTENT(IN)::isg,mm
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(OUT)::zps0
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(IN)::zp0,zps2
COMPLEX(KIND=DBL),INTENT(IN)::zgr,fmas,zs1
COMPLEX(KIND=DBL)::z  
zps0(1)=DCMPLX(dnou(0),dnou(0))
zps0(2)=DCMPLX(dnou(0),dnou(0))
z=isg*zgr
zps0(3)=(zp0(2)*zps2(1)-zp0(4)*zps2(2))*zs1*z
zps0(4)=(-zp0(3)*zps2(1)+zp0(1)*zps2(2))*zs1*z
IF(fmas.NE.dnou(0))THEN
   zps0(1)=-zs1*zps2(1)*fmas*z
   zps0(2)=-zs1*zps2(2)*fmas*z
ENDIF
END SUBROUTINE Helac_xffsr

SUBROUTINE Helac_xffsrd(zps0,zp0,zs1,zps2,mm,zgr,isg,fmas,deriv,ip0array)
IMPLICIT NONE
INTEGER,INTENT(IN)::isg,mm
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(OUT)::zps0
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(IN)::zp0,zps2
COMPLEX(KIND=DBL),INTENT(IN)::zgr,fmas,zs1
INTEGER,DIMENSION(10),INTENT(IN)::ip0array
COMPLEX(KIND=DBL),DIMENSION(4)::zeps
INTEGER,INTENT(IN)::deriv
INTEGER::k

IF(deriv.LE.0)THEN
	CALL Helac_xffsr(zps0,zp0,zs1,zps2,mm,zgr,isg,fmas)
	RETURN
ENDIF
IF(Helac_level(Pwavenum,deriv).GT.1)THEN
	zps0(1:4)=0
	RETURN
ENDIF
CALL Deriv_eps(zeps,k,deriv)
IF(ip0array(k).EQ.0)THEN
	zps0(1:4)=0
	RETURN
ENDIF
zeps(1:4)=ip0array(k)*zeps(1:4)
CALL Helac_xffsr(zps0,zeps,zs1,zps2,mm,zgr,isg,DCMPLX(0d0))
        
END SUBROUTINE Helac_xffsrd
      
SUBROUTINE Helac_xffsr0(zps0,zs1,zps2,zgr,isg)
IMPLICIT NONE
INTEGER,INTENT(IN)::isg
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(OUT)::zps0
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(IN)::zps2
COMPLEX(KIND=DBL),INTENT(IN)::zgr,zs1
COMPLEX(KIND=DBL)::z
zps0(3)=DCMPLX(dnou(0),dnou(0))
zps0(4)=DCMPLX(dnou(0),dnou(0))
z=zgr*isg
zps0(1)=-zs1*zps2(1)*z
zps0(2)=-zs1*zps2(2)*z
END SUBROUTINE Helac_xffsr0
             
SUBROUTINE Helac_xffsl0(zps0,zs1,zps2,zgl,isg)
IMPLICIT NONE
INTEGER,INTENT(IN)::isg
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(OUT)::zps0
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(IN)::zps2
COMPLEX(KIND=DBL),INTENT(IN)::zgl,zs1
COMPLEX(KIND=DBL)::z       
zps0(1)=DCMPLX(dnou(0),dnou(0))
zps0(2)=DCMPLX(dnou(0),dnou(0))
z=zgl*isg
zps0(3)=-zs1*zps2(3)*z
zps0(4)=-zs1*zps2(4)*z
END SUBROUTINE Helac_xffsl0
                
SUBROUTINE Helac_xaffsl(zpb0,zp0,zs1,zpb2,mm,zgl,isg,fmas)
IMPLICIT NONE
INTEGER,INTENT(IN)::isg,mm
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(OUT)::zpb0
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(IN)::zp0,zpb2
COMPLEX(KIND=DBL),INTENT(IN)::zgl,fmas,zs1
COMPLEX(KIND=DBL)::z    
zpb0(3)=DCMPLX(dnou(0),dnou(0))
zpb0(4)=DCMPLX(dnou(0),dnou(0))
z=isg*zgl
zpb0(1)=(-zp0(2)*zpb2(3)+zp0(3)*zpb2(4))*zs1*z
zpb0(2)=( zp0(4)*zpb2(3)-zp0(1)*zpb2(4))*zs1*z
IF(fmas.NE.dnou(0))THEN
    zpb0(3)=-zs1*zpb2(3)*fmas*z          
    zpb0(4)=-zs1*zpb2(4)*fmas*z
ENDIF
END SUBROUTINE Helac_xaffsl

SUBROUTINE Helac_xaffsld(zpb0,zp0,zs1,zpb2,mm,zgl,isg,fmas,deriv,ip0array)
IMPLICIT NONE
INTEGER,INTENT(IN)::isg,mm
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(OUT)::zpb0
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(IN)::zp0,zpb2
COMPLEX(KIND=DBL),INTENT(IN)::zgl,fmas,zs1
INTEGER,DIMENSION(10),INTENT(IN)::ip0array
COMPLEX(KIND=DBL),DIMENSION(4)::zeps
INTEGER,INTENT(IN)::deriv
INTEGER::k

IF(deriv.LE.0)THEN
	CALL Helac_xaffsl(zpb0,zp0,zs1,zpb2,mm,zgl,isg,fmas)
	RETURN
ENDIF
IF(Helac_level(Pwavenum,deriv).GT.1)THEN
	zpb0(1:4)=0
	RETURN
ENDIF
CALL Deriv_eps(zeps,k,deriv)
IF(ip0array(k).EQ.0)THEN
	zpb0(1:4)=0
	RETURN
ENDIF
zeps(1:4)=ip0array(k)*zeps(1:4)
CALL Helac_xaffsl(zpb0,zeps,zs1,zpb2,mm,zgl,isg,DCMPLX(0d0))
        
END SUBROUTINE Helac_xaffsld
        
SUBROUTINE Helac_xaffsr(zpb0,zp0,zs1,zpb2,mm,zgr,isg,fmas)
IMPLICIT NONE
INTEGER,INTENT(IN)::isg,mm
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(OUT)::zpb0
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(IN)::zp0,zpb2
COMPLEX(KIND=DBL),INTENT(IN)::zgr,fmas,zs1
COMPLEX(KIND=DBL)::z  
zpb0(1)=DCMPLX(dnou(0),dnou(0))
zpb0(2)=DCMPLX(dnou(0),dnou(0))
z=isg*zgr
zpb0(3)=(-zp0(1)*zpb2(1)-zp0(3)*zpb2(2))*zs1*z
zpb0(4)=(-zp0(4)*zpb2(1)-zp0(2)*zpb2(2))*zs1*z
IF(fmas.NE.dnou(0))THEN
    zpb0(1)=-zs1*zpb2(1)*fmas*z 
    zpb0(2)=-zs1*zpb2(2)*fmas*z
ENDIF
END SUBROUTINE Helac_xaffsr

SUBROUTINE Helac_xaffsrd(zpb0,zp0,zs1,zpb2,mm,zgr,isg,fmas,deriv,ip0array)
IMPLICIT NONE
INTEGER,INTENT(IN)::isg,mm
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(OUT)::zpb0
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(IN)::zp0,zpb2
COMPLEX(KIND=DBL),INTENT(IN)::zgr,fmas,zs1
INTEGER,DIMENSION(10),INTENT(IN)::ip0array
COMPLEX(KIND=DBL),DIMENSION(4)::zeps
INTEGER,INTENT(IN)::deriv
INTEGER::k

IF(deriv.LE.0)THEN
	CALL Helac_xaffsr(zpb0,zp0,zs1,zpb2,mm,zgr,isg,fmas)
	RETURN
ENDIF
IF(Helac_level(Pwavenum,deriv).GT.1)THEN
	zpb0(1:4)=0
	RETURN
ENDIF
CALL Deriv_eps(zeps,k,deriv)
IF(ip0array(k).EQ.0)THEN
	zpb0(1:4)=0
	RETURN
ENDIF
zeps(1:4)=ip0array(k)*zeps(1:4)
CALL Helac_xaffsr(zpb0,zeps,zs1,zpb2,mm,zgr,isg,DCMPLX(0d0))
        
END SUBROUTINE Helac_xaffsrd
      
SUBROUTINE Helac_xaffsr0(zpb0,zs1,zpb2,zgr,isg)
IMPLICIT NONE
INTEGER,INTENT(IN)::isg
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(OUT)::zpb0
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(IN)::zpb2
COMPLEX(KIND=DBL),INTENT(IN)::zgr,zs1
COMPLEX(KIND=DBL)::z    
zpb0(3)=DCMPLX(dnou(0),dnou(0))
zpb0(4)=DCMPLX(dnou(0),dnou(0))
z=isg*zgr
zpb0(1)=-zs1*zpb2(1)*z
zpb0(2)=-zs1*zpb2(2)*z
END SUBROUTINE Helac_xaffsr0
      
SUBROUTINE Helac_xaffsl0(zpb0,zs1,zpb2,zgl,isg)
IMPLICIT NONE
INTEGER,INTENT(IN)::isg
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(OUT)::zpb0
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(IN)::zpb2
COMPLEX(KIND=DBL),INTENT(IN)::zgl,zs1
COMPLEX(KIND=DBL)::z
zpb0(1)=DCMPLX(dnou(0),dnou(0))
zpb0(2)=DCMPLX(dnou(0),dnou(0))
z=isg*zgl
zpb0(3)=-zs1*zpb2(3)*z
zpb0(4)=-zs1*zpb2(4)*z
END SUBROUTINE Helac_xaffsl0
      
SUBROUTINE Helac_xsff(zs0,zp0,zpb1,zps2,mm,zgl,zgr,isg)
IMPLICIT NONE
INTEGER,INTENT(IN)::isg,mm
COMPLEX(KIND=DBL),INTENT(OUT)::zs0
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(IN)::zp0,zpb1,zps2 
COMPLEX(KIND=DBL),INTENT(IN)::zgl,zgr
COMPLEX(KIND=DBL)::z      
z=isg
zs0=(-zpb1(1)*zps2(1)-zpb1(2)*zps2(2))*z*zgr+(-zpb1(3)*zps2(3)-zpb1(4)*zps2(4))*z*zgl
END SUBROUTINE Helac_xsff

SUBROUTINE Helac_xvvs(zb0,zp0,zb1,zs2,mm,zg,isg)
IMPLICIT NONE
INTEGER,INTENT(IN)::isg,mm
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(OUT)::zb0
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(IN)::zp0,zb1 
COMPLEX(KIND=DBL),INTENT(IN)::zg,zs2
COMPLEX(KIND=DBL)::z
INTEGER::k        
z=zg*isg
DO k=1,4
  zb0(k)=zb1(k)*zs2*z
ENDDO
END SUBROUTINE Helac_xvvs
 
SUBROUTINE Helac_xsvv(zs0,zp0,zb1,zb2,mm,zg,isg)
IMPLICIT NONE
INTEGER,INTENT(IN)::isg,mm
COMPLEX(KIND=DBL),INTENT(OUT)::zs0
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(IN)::zp0,zb1,zb2 
COMPLEX(KIND=DBL),INTENT(IN)::zg
COMPLEX(KIND=DBL)::z       
z=zg*isg
zs0=-Helac_zprod(zb1,zb2)*z
END SUBROUTINE Helac_xsvv
       
SUBROUTINE Helac_xvvss(zb0,zp0,zb1,zs2,zs3,mm,zg,isg)
IMPLICIT NONE
INTEGER,INTENT(IN)::isg,mm
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(OUT)::zb0
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(IN)::zp0,zb1
COMPLEX(KIND=DBL),INTENT(IN)::zg,zs2,zs3 
COMPLEX(KIND=DBL)::z
INTEGER::k              
z=zg*isg
DO k=1,4 
   zb0(k)=zb1(k)*zs2*zs3*z 
ENDDO       
END SUBROUTINE Helac_xvvss

SUBROUTINE Helac_xssvv(zs0,zp0,zs1,zb2,zb3,mm,zg,isg)
IMPLICIT NONE
INTEGER,INTENT(IN)::isg,mm
COMPLEX(KIND=DBL),INTENT(OUT)::zs0
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(IN)::zp0,zb2,zb3
COMPLEX(KIND=DBL),INTENT(IN)::zg,zs1 
COMPLEX(KIND=DBL)::z       
z=zg*isg
zs0=-Helac_zprod(zb2,zb3)*zs1*z
END SUBROUTINE Helac_xssvv
       
SUBROUTINE Helac_xs3(zs0,zp0,zs1,zs2,mm,zg,isg)
IMPLICIT NONE
INTEGER,INTENT(IN)::isg,mm
COMPLEX(KIND=DBL),INTENT(OUT)::zs0
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(IN)::zp0
COMPLEX(KIND=DBL),INTENT(IN)::zg,zs1,zs2
COMPLEX(KIND=DBL)::z
z=zg*isg
zs0=-zs1*zs2*z
END SUBROUTINE Helac_xs3
     
SUBROUTINE Helac_xs4(zs0,zp0,zs1,zs2,zs3,mm,zg,isg)
IMPLICIT NONE
INTEGER,INTENT(IN)::isg,mm
COMPLEX(KIND=DBL),INTENT(OUT)::zs0
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(IN)::zp0
COMPLEX(KIND=DBL),INTENT(IN)::zg,zs1,zs2,zs3
COMPLEX(KIND=DBL)::z
z=zg*isg
zs0=-zs1*zs2*zs3*z
END SUBROUTINE Helac_xs4
       
SUBROUTINE Helac_xssv(zs0,zp0,zs1,zp1,zb2,mm,zg,isg)
IMPLICIT NONE
INTEGER,INTENT(IN)::isg,mm
COMPLEX(KIND=DBL),INTENT(OUT)::zs0
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(IN)::zp0,zp1,zb2
COMPLEX(KIND=DBL),INTENT(IN)::zg,zs1
COMPLEX(KIND=DBL)::z
z=zg*isg
zs0=(Helac_zprod(zp0,zb2)+Helac_zprod(zp1,zb2))*zs1*z
END SUBROUTINE Helac_xssv

SUBROUTINE Helac_xssvd(zs0,zp0,zs1,zp1,zb2,mm,zg,isg,deriv,ip0array,ip1array)
IMPLICIT NONE
INTEGER,INTENT(IN)::isg,mm
COMPLEX(KIND=DBL),INTENT(OUT)::zs0
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(IN)::zp0,zp1,zb2
COMPLEX(KIND=DBL),INTENT(IN)::zg,zs1
INTEGER,DIMENSION(10),INTENT(IN)::ip0array,ip1array
COMPLEX(KIND=DBL),DIMENSION(4)::zeps,zeps1
INTEGER,INTENT(IN)::deriv
INTEGER::k

IF(deriv.LE.0)THEN
	CALL Helac_xssv(zs0,zp0,zs1,zp1,zb2,mm,zg,isg)
	RETURN
ENDIF
IF(Helac_level(Pwavenum,deriv).GT.1)THEN
	zs0=0
	RETURN
ENDIF
CALL Deriv_eps(zeps,k,deriv)
zeps1(1:4)=ip1array(k)*zeps(1:4)
zeps(1:4)=ip0array(k)*zeps(1:4)
CALL Helac_xssv(zs0,zeps,zs1,zeps1,zb2,mm,zg,isg)
        
END SUBROUTINE Helac_xssvd
       
SUBROUTINE Helac_xvss(zb0,zp0,zs1,zp1,zs2,zp2,mm,zg,isg)
IMPLICIT NONE
INTEGER,INTENT(IN)::isg,mm
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(OUT)::zb0
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(IN)::zp0,zp1,zp2
COMPLEX(KIND=DBL),INTENT(IN)::zg,zs1,zs2
COMPLEX(KIND=DBL)::z
INTEGER::k
z=zg*isg
DO k=1,4
   zb0(k)=zs1*zs2*(zp1(k)-zp2(k))*z
ENDDO
END SUBROUTINE Helac_xvss

SUBROUTINE Helac_xvssd(zb0,zp0,zs1,zp1,zs2,zp2,mm,zg,isg,deriv,ip1array,ip2array)
IMPLICIT NONE
INTEGER,INTENT(IN)::isg,mm
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(OUT)::zb0
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(IN)::zp0,zp1,zp2
COMPLEX(KIND=DBL),INTENT(IN)::zg,zs1,zs2
INTEGER,DIMENSION(10),INTENT(IN)::ip1array,ip2array
COMPLEX(KIND=DBL),DIMENSION(4)::zeps,zeps1
INTEGER,INTENT(IN)::deriv
INTEGER::k

IF(deriv.LE.0)THEN
	CALL Helac_xvss(zb0,zp0,zs1,zp1,zs2,zp2,mm,zg,isg)
	RETURN
ENDIF
IF(Helac_level(Pwavenum,deriv).GT.1)THEN
	zb0(1:4)=0
	RETURN
ENDIF
CALL Deriv_eps(zeps,k,deriv)
zeps1(1:4)=ip1array(k)*zeps(1:4)
zeps(1:4)=ip2array(k)*zeps(1:4)
CALL Helac_xvss(zb0,zp0,zs1,zeps1,zs2,zeps,mm,zg,isg)
        
END SUBROUTINE Helac_xvssd

SUBROUTINE Helac_propag(zy,m0,zp0)       
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(INOUT)::zy      
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(IN)::zp0
INTEGER,INTENT(IN)::m0
INTEGER::mm
COMPLEX(KIND=DBL)::zp,z,zmas,zex
INTEGER::i
zp =Helac_zprod(zp0,zp0)
!  PHOTON,GLUON 
IF(m0.EQ.31.OR.m0.EQ.35)THEN
   z=zp
   IF(iunitary.EQ.0)THEN
        DO i=1,4
           zy(i)=zy(i)/z
        ENDDO
   ELSEIF(iunitary.EQ.1)THEN
        DO i=1,4
           zy(i)=zy(i)/z
        ENDDO
   ENDIF
   RETURN
ENDIF
!  Z
IF(m0.EQ.32)THEN
!  iwidth=1
   IF(iwidth.EQ.1)THEN
        mm=m0
	! include'compl_mass.h'
		IF(parmas(mm).GT.dnou(0).AND.iwidth.EQ.1)THEN
            zmas=parmas(mm)*DCMPLX(DSQRT( dnou(1)/dnou(2)&
			+dnou(1)/dnou(2)*DSQRT(dnou(1)+parwid(mm)**2/parmas(mm)**2)), &
			-DSQRT(-dnou(1)/dnou(2)+&
			dnou(1)/dnou(2)*DSQRT(dnou(1)+parwid(mm)**2/parmas(mm)**2)))
        ELSE
            zmas=DCMPLX(parmas(mm),dnou(0))
        ENDIF
        z=zp-zmas*zmas
        IF(iunitary.EQ.0)THEN
            DO i=1,4
               zy(i)=zy(i)/z
            ENDDO
        ELSEIF(iunitary.EQ.1)THEN
            zex=Helac_zprod(zy,zp0)
            DO i=1,4
               zy(i)=(zy(i)-zp0(i)*zex/zmas/zmas)/z
            ENDDO
        ENDIF
!  iwidth=0
    ELSEIF(iwidth.EQ.0)THEN
        z=zp-parmas(m0)**2+DCMPLX(dnou(0),dnou(1))*parmas(m0)*parwid(m0)
        zmas=DCMPLX(parmas(m0),dnou(0))
        IF(iunitary.EQ.0)THEN
            DO i=1,4
               zy(i)=zy(i)/z
            ENDDO
        ELSEIF(iunitary.EQ.1)THEN
            zex=Helac_zprod(zy,zp0)
            DO i=1,4
               zy(i)=(zy(i)-zp0(i)*zex/zmas/zmas)/z
            ENDDO
        ENDIF
    ENDIF
    RETURN
ENDIF
!  W
IF(m0.EQ.33.OR.m0.EQ.34)THEN
!  iwidth=1
    IF(iwidth.EQ.1)THEN
        mm=m0
!        include 'compl_mass.h'
        IF(parmas(mm).GT.dnou(0).AND.iwidth.EQ.1)THEN
            zmas=parmas(mm)*DCMPLX(DSQRT( dnou(1)/dnou(2)&
			+dnou(1)/dnou(2)*DSQRT(dnou(1)+parwid(mm)**2/parmas(mm)**2)), &
			-DSQRT(-dnou(1)/dnou(2)+&
			dnou(1)/dnou(2)*DSQRT(dnou(1)+parwid(mm)**2/parmas(mm)**2)))
        ELSE
            zmas=DCMPLX(parmas(mm),dnou(0))
        ENDIF
        z=zp-zmas*zmas

        IF(iunitary.EQ.0)THEN
            DO i=1,4
                zy(i)=zy(i)/z
            ENDDO
        ELSEIF(iunitary.EQ.1)THEN
            zex=Helac_zprod(zy,zp0)
            DO i=1,4
                zy(i)=(zy(i)-zp0(i)*zex/zmas/zmas)/z
            ENDDO
        ENDIF
!  iwidth=0
     ELSEIF(iwidth.EQ.0)THEN
        z=zp-parmas(m0)**2+DCMPLX(dnou(0),dnou(1))*parmas(m0)*parwid(m0)
        zmas=DCMPLX(parmas(m0),dnou(0))

        IF(iunitary.EQ.0)THEN
            DO i=1,4
                zy(i)=zy(i)/z
            ENDDO
        ELSEIF(iunitary.EQ.1)THEN
            zex=Helac_zprod(zy,zp0)
            DO i=1,4
                zy(i)=(zy(i)-zp0(i)*zex/zmas/zmas)/z
            ENDDO
        ENDIF
     ENDIF
     RETURN
ENDIF
! S
IF(m0.GE.41)THEN
   IF(iwidth.eq.1)THEN
        mm=m0
!        include 'compl_mass.h'
        IF(parmas(mm).GT.dnou(0).AND.iwidth.EQ.1)THEN
            zmas=parmas(mm)*DCMPLX(SQRT( dnou(1)/dnou(2)&
			+dnou(1)/dnou(2)*SQRT(dnou(1)+parwid(mm)**2/parmas(mm)**2)), &
			-SQRT(-dnou(1)/dnou(2)+&
			dnou(1)/dnou(2)*SQRT(dnou(1)+parwid(mm)**2/parmas(mm)**2)))
        ELSE
            zmas=DCMPLX(parmas(mm),dnou(0))
        ENDIF
        z=zp-zmas*zmas
   ELSEIF(iwidth.EQ.0)THEN
        z=zp-parmas(m0)**2+DCMPLX(dnou(0),dnou(1))*parmas(m0)*parwid(m0)
   ENDIF
   DO i=1,1
        zy(i)=zy(i)/z
   ENDDO
   RETURN
ENDIF
! F
z=zp-parmas(m0)**2+DCMPLX(dnou(0),dnou(1))*parmas(m0)*parwid(m0)
DO i=1,4
   zy(i)=zy(i)/z
ENDDO
END SUBROUTINE Helac_propag

SUBROUTINE Helac_propagd(zy,m0,zp0,deriv,ip0array)
IMPLICIT NONE
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(INOUT)::zy      
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(IN)::zp0
INTEGER,INTENT(IN)::m0,deriv
INTEGER,DIMENSION(10),INTENT(IN)::ip0array
INTEGER::mm,i,k,num,k1,k2
COMPLEX(KIND=DBL)::zp,z,z1,zmas,zex,zex1,z2,zex2,zex3,zex4
COMPLEX(KIND=DBL),DIMENSION(4)::zeps,zeps1,zeps2

IF(deriv.LE.0)THEN
	CALL Helac_propag(zy,m0,zp0)
	RETURN
ENDIF
num=Helac_level(Pwavenum,deriv)
SELECT CASE(num)
CASE(1)
	CALL Deriv_eps(zeps,k,deriv)
	zeps(1:4)=ip0array(k)*zeps(1:4)
	zp =Helac_zprod(zp0,zp0)
!  PHOTON,GLUON 
	IF(m0.EQ.31.OR.m0.EQ.35)THEN
		z=zp**2
		z1=DCMPLX(2d0)*Helac_zprod(zp0,zeps)
		IF(iunitary.EQ.0)THEN
			DO i=1,4
				zy(i)=-z1*zy(i)/z
			ENDDO
		ELSEIF(iunitary.EQ.1)THEN
			DO i=1,4
				zy(i)=-z1*zy(i)/z
			ENDDO
		ENDIF
		RETURN
	ENDIF
!  Z
	IF(m0.EQ.32)THEN
!  iwidth=1
		IF(iwidth.EQ.1)THEN
			mm=m0
	! include'compl_mass.h'
			IF(parmas(mm).GT.dnou(0).AND.iwidth.EQ.1)THEN
				zmas=parmas(mm)*DCMPLX(DSQRT( dnou(1)/dnou(2)&
				+dnou(1)/dnou(2)*DSQRT(dnou(1)+parwid(mm)**2/parmas(mm)**2)), &
				-DSQRT(-dnou(1)/dnou(2)+&
				dnou(1)/dnou(2)*DSQRT(dnou(1)+parwid(mm)**2/parmas(mm)**2)))
			ELSE
				zmas=DCMPLX(parmas(mm),dnou(0))
			ENDIF
			z=zp-zmas*zmas
			z1=DCMPLX(2d0)*Helac_zprod(zp0,zeps)
			IF(iunitary.EQ.0)THEN
				DO i=1,4
					zy(i)=-z1*zy(i)/z**2
				ENDDO
			ELSEIF(iunitary.EQ.1)THEN
				zex=Helac_zprod(zy,zp0)
				zex1=Helac_zprod(zy,zeps)
				DO i=1,4
					zy(i)=-z1*(zy(i)-zp0(i)*zex/zmas/zmas)/z**2&
					-zeps(i)*zex/zmas/zmas/z-zp0(i)*zex1/zmas/zmas/z
				ENDDO
			ENDIF
!  iwidth=0
		ELSEIF(iwidth.EQ.0)THEN
			z=zp-parmas(m0)**2+DCMPLX(dnou(0),dnou(1))*parmas(m0)*parwid(m0)
			zmas=DCMPLX(parmas(m0),dnou(0))
			z1=DCMPLX(2d0)*Helac_zprod(zp0,zeps)
			IF(iunitary.EQ.0)THEN
				DO i=1,4
					zy(i)=-z1*zy(i)/z**2
				ENDDO
			ELSEIF(iunitary.EQ.1)THEN
				zex=Helac_zprod(zy,zp0)
				zex1=Helac_zprod(zy,zeps)
				DO i=1,4
					zy(i)=-z1*(zy(i)-zp0(i)*zex/zmas/zmas)/z**2&
					-zeps(i)*zex/zmas/zmas/z-zp0(i)*zex1/zmas/zmas/z
				ENDDO
			ENDIF
		ENDIF
		RETURN
	ENDIF
!  W
	IF(m0.EQ.33.OR.m0.EQ.34)THEN
!  iwidth=1
		IF(iwidth.EQ.1)THEN
			mm=m0
!        include 'compl_mass.h'
			IF(parmas(mm).GT.dnou(0).AND.iwidth.EQ.1)THEN
				zmas=parmas(mm)*DCMPLX(DSQRT( dnou(1)/dnou(2)&
				+dnou(1)/dnou(2)*DSQRT(dnou(1)+parwid(mm)**2/parmas(mm)**2)), &
				-DSQRT(-dnou(1)/dnou(2)+&
				dnou(1)/dnou(2)*DSQRT(dnou(1)+parwid(mm)**2/parmas(mm)**2)))
			ELSE
				zmas=DCMPLX(parmas(mm),dnou(0))
			ENDIF
			z=zp-zmas*zmas
			z1=DCMPLX(2d0)*Helac_zprod(zp0,zeps)

			IF(iunitary.EQ.0)THEN
				DO i=1,4
					zy(i)=-z1*zy(i)/z**2
				ENDDO
			ELSEIF(iunitary.EQ.1)THEN
				zex=Helac_zprod(zy,zp0)
				zex1=Helac_zprod(zy,zeps)
				DO i=1,4
					zy(i)=-z1*(zy(i)-zp0(i)*zex/zmas/zmas)/z**2&
					-zeps(i)*zex/zmas/zmas/z-zp0(i)*zex1/zmas/zmas/z
				ENDDO
			ENDIF
!  iwidth=0
		ELSEIF(iwidth.EQ.0)THEN
			z=zp-parmas(m0)**2+DCMPLX(dnou(0),dnou(1))*parmas(m0)*parwid(m0)
			zmas=DCMPLX(parmas(m0),dnou(0))
			z1=DCMPLX(2d0)*Helac_zprod(zp0,zeps)

			IF(iunitary.EQ.0)THEN
				DO i=1,4
					zy(i)=-z1*zy(i)/z**2
				ENDDO
			ELSEIF(iunitary.EQ.1)THEN
				zex=Helac_zprod(zy,zp0)
				zex1=Helac_zprod(zy,zeps)
				DO i=1,4
					zy(i)=-z1*(zy(i)-zp0(i)*zex/zmas/zmas)/z**2&
					-zeps(i)*zex/zmas/zmas/z-zp0(i)*zex1/zmas/zmas/z
				ENDDO
			ENDIF
		ENDIF
		RETURN
	ENDIF
! S
	IF(m0.GE.41)THEN
		IF(iwidth.eq.1)THEN
			mm=m0
!        include 'compl_mass.h'
			IF(parmas(mm).GT.dnou(0).AND.iwidth.EQ.1)THEN
				zmas=parmas(mm)*DCMPLX(SQRT( dnou(1)/dnou(2)&
				+dnou(1)/dnou(2)*SQRT(dnou(1)+parwid(mm)**2/parmas(mm)**2)), &
				-SQRT(-dnou(1)/dnou(2)+&
				dnou(1)/dnou(2)*SQRT(dnou(1)+parwid(mm)**2/parmas(mm)**2)))
			ELSE
				zmas=DCMPLX(parmas(mm),dnou(0))
			ENDIF
			z=zp-zmas*zmas
		ELSEIF(iwidth.EQ.0)THEN
			z=zp-parmas(m0)**2+DCMPLX(dnou(0),dnou(1))*parmas(m0)*parwid(m0)
		ENDIF
		z=z**2
		z1=DCMPLX(2d0)*Helac_zprod(zp0,zeps)
		DO i=1,1
			zy(i)=-z1*zy(i)/z
		ENDDO
		RETURN
	ENDIF
! F
	z=zp-parmas(m0)**2+DCMPLX(dnou(0),dnou(1))*parmas(m0)*parwid(m0)
	z=z**2
	z1=DCMPLX(2d0)*Helac_zprod(zp0,zeps)
	DO i=1,4
		zy(i)=-z1*zy(i)/z
	ENDDO
	RETURN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE(2)
	CALL Deriv_eps(zeps1,k1,deriv)
	zeps1(1:4)=ip0array(k1)*zeps1(1:4)
	CALL Deriv_eps(zeps2,k2,deriv-2**(k1-1))
	zeps2(1:4)=ip0array(k2)*zeps2(1:4)
	zp =Helac_zprod(zp0,zp0)
!  PHOTON,GLUON 
	IF(m0.EQ.31.OR.m0.EQ.35)THEN
		z=zp
		z1=DCMPLX(8d0)*Helac_zprod(zp0,zeps1)*Helac_zprod(zp0,zeps2)
		z2=DCMPLX(2d0)*Helac_zprod(zeps2,zeps1)
		IF(iunitary.EQ.0)THEN
			DO i=1,4
				zy(i)=-z2*zy(i)/z**2+z1*zy(i)/z**3
			ENDDO
		ELSEIF(iunitary.EQ.1)THEN
			DO i=1,4
				zy(i)=-z2*zy(i)/z**2+z1*zy(i)/z**3
			ENDDO
		ENDIF
		RETURN
	ENDIF
!  Z
	IF(m0.EQ.32)THEN
!  iwidth=1
		IF(iwidth.EQ.1)THEN
			mm=m0
	! include'compl_mass.h'
			IF(parmas(mm).GT.dnou(0).AND.iwidth.EQ.1)THEN
				zmas=parmas(mm)*DCMPLX(DSQRT( dnou(1)/dnou(2)&
				+dnou(1)/dnou(2)*DSQRT(dnou(1)+parwid(mm)**2/parmas(mm)**2)), &
				-DSQRT(-dnou(1)/dnou(2)+&
				dnou(1)/dnou(2)*DSQRT(dnou(1)+parwid(mm)**2/parmas(mm)**2)))
			ELSE
				zmas=DCMPLX(parmas(mm),dnou(0))
			ENDIF
			z=zp-zmas*zmas
			z1=DCMPLX(8d0)*Helac_zprod(zp0,zeps1)*Helac_zprod(zp0,zeps2)
			z2=DCMPLX(2d0)*Helac_zprod(zeps2,zeps1)
			IF(iunitary.EQ.0)THEN
				DO i=1,4
					zy(i)=-z2*zy(i)/z**2+z1*zy(i)/z**3
				ENDDO
			ELSEIF(iunitary.EQ.1)THEN
				zex=Helac_zprod(zy,zp0)
				zex1=Helac_zprod(zy,zeps1)
				zex2=Helac_zprod(zy,zeps2)
				zex3=DCMPLX(2d0)*Helac_zprod(zp0,zeps2)
				zex4=DCMPLX(2d0)*Helac_zprod(zp0,zeps1)
				DO i=1,4
					zy(i)=-z2*(zy(i)-zp0(i)*zex/zmas/zmas)/z**2&
					+z1*(zy(i)-zp0(i)*zex/zmas/zmas)/z**3&
					+zex3*(zeps1(i)*zex+zp0(i)*zex1)/zmas/zmas/z**2&
					+zex4*(zeps2(i)*zex+zp0(i)*zex2)/zmas/zmas/z**2&
					-(zeps1(i)*zex2+zeps2(i)*zex1)/zmas/zmas/z
				ENDDO
			ENDIF
!  iwidth=0
		ELSEIF(iwidth.EQ.0)THEN
			z=zp-parmas(m0)**2+DCMPLX(dnou(0),dnou(1))*parmas(m0)*parwid(m0)
			zmas=DCMPLX(parmas(m0),dnou(0))
			z1=DCMPLX(8d0)*Helac_zprod(zp0,zeps1)*Helac_zprod(zp0,zeps2)
			z2=DCMPLX(2d0)*Helac_zprod(zeps2,zeps1)
			IF(iunitary.EQ.0)THEN
				DO i=1,4
					zy(i)=-z2*zy(i)/z**2+z1*zy(i)/z**3
				ENDDO
			ELSEIF(iunitary.EQ.1)THEN
				zex=Helac_zprod(zy,zp0)
				zex1=Helac_zprod(zy,zeps1)
				zex2=Helac_zprod(zy,zeps2)
				zex3=DCMPLX(2d0)*Helac_zprod(zp0,zeps2)
				zex4=DCMPLX(2d0)*Helac_zprod(zp0,zeps1)
				DO i=1,4
					zy(i)=-z2*(zy(i)-zp0(i)*zex/zmas/zmas)/z**2&
					+z1*(zy(i)-zp0(i)*zex/zmas/zmas)/z**3&
					+zex3*(zeps1(i)*zex+zp0(i)*zex1)/zmas/zmas/z**2&
					+zex4*(zeps2(i)*zex+zp0(i)*zex2)/zmas/zmas/z**2&
					-(zeps1(i)*zex2+zeps2(i)*zex1)/zmas/zmas/z
				ENDDO
			ENDIF
		ENDIF
		RETURN
	ENDIF
!  W
	IF(m0.EQ.33.OR.m0.EQ.34)THEN
!  iwidth=1
		IF(iwidth.EQ.1)THEN
			mm=m0
!        include 'compl_mass.h'
			IF(parmas(mm).GT.dnou(0).AND.iwidth.EQ.1)THEN
				zmas=parmas(mm)*DCMPLX(DSQRT( dnou(1)/dnou(2)&
				+dnou(1)/dnou(2)*DSQRT(dnou(1)+parwid(mm)**2/parmas(mm)**2)), &
				-DSQRT(-dnou(1)/dnou(2)+&
				dnou(1)/dnou(2)*DSQRT(dnou(1)+parwid(mm)**2/parmas(mm)**2)))
			ELSE
				zmas=DCMPLX(parmas(mm),dnou(0))
			ENDIF
			z=zp-zmas*zmas
			z1=DCMPLX(8d0)*Helac_zprod(zp0,zeps1)*Helac_zprod(zp0,zeps2)
			z2=DCMPLX(2d0)*Helac_zprod(zeps2,zeps1)
			IF(iunitary.EQ.0)THEN
				DO i=1,4
					zy(i)=-z2*zy(i)/z**2+z1*zy(i)/z**3
				ENDDO
			ELSEIF(iunitary.EQ.1)THEN
				zex=Helac_zprod(zy,zp0)
				zex1=Helac_zprod(zy,zeps1)
				zex2=Helac_zprod(zy,zeps2)
				zex3=DCMPLX(2d0)*Helac_zprod(zp0,zeps2)
				zex4=DCMPLX(2d0)*Helac_zprod(zp0,zeps1)
				DO i=1,4
					zy(i)=-z2*(zy(i)-zp0(i)*zex/zmas/zmas)/z**2&
					+z1*(zy(i)-zp0(i)*zex/zmas/zmas)/z**3&
					+zex3*(zeps1(i)*zex+zp0(i)*zex1)/zmas/zmas/z**2&
					+zex4*(zeps2(i)*zex+zp0(i)*zex2)/zmas/zmas/z**2&
					-(zeps1(i)*zex2+zeps2(i)*zex1)/zmas/zmas/z
				ENDDO
			ENDIF
!  iwidth=0
		ELSEIF(iwidth.EQ.0)THEN
			z=zp-parmas(m0)**2+DCMPLX(dnou(0),dnou(1))*parmas(m0)*parwid(m0)
			zmas=DCMPLX(parmas(m0),dnou(0))
			z1=DCMPLX(8d0)*Helac_zprod(zp0,zeps1)*Helac_zprod(zp0,zeps2)
			z2=DCMPLX(2d0)*Helac_zprod(zeps2,zeps1)
			IF(iunitary.EQ.0)THEN
				DO i=1,4
					zy(i)=-z2*zy(i)/z**2+z1*zy(i)/z**3
				ENDDO
			ELSEIF(iunitary.EQ.1)THEN
				zex=Helac_zprod(zy,zp0)
				zex1=Helac_zprod(zy,zeps1)
				zex2=Helac_zprod(zy,zeps2)
				zex3=DCMPLX(2d0)*Helac_zprod(zp0,zeps2)
				zex4=DCMPLX(2d0)*Helac_zprod(zp0,zeps1)
				DO i=1,4
					zy(i)=-z2*(zy(i)-zp0(i)*zex/zmas/zmas)/z**2&
					+z1*(zy(i)-zp0(i)*zex/zmas/zmas)/z**3&
					+zex3*(zeps1(i)*zex+zp0(i)*zex1)/zmas/zmas/z**2&
					+zex4*(zeps2(i)*zex+zp0(i)*zex2)/zmas/zmas/z**2&
					-(zeps1(i)*zex2+zeps2(i)*zex1)/zmas/zmas/z
				ENDDO
			ENDIF
		ENDIF
		RETURN
	ENDIF
! S
	IF(m0.GE.41)THEN
		IF(iwidth.eq.1)THEN
			mm=m0
!        include 'compl_mass.h'
			IF(parmas(mm).GT.dnou(0).AND.iwidth.EQ.1)THEN
				zmas=parmas(mm)*DCMPLX(SQRT( dnou(1)/dnou(2)&
				+dnou(1)/dnou(2)*SQRT(dnou(1)+parwid(mm)**2/parmas(mm)**2)), &
				-SQRT(-dnou(1)/dnou(2)+&
				dnou(1)/dnou(2)*SQRT(dnou(1)+parwid(mm)**2/parmas(mm)**2)))
			ELSE
				zmas=DCMPLX(parmas(mm),dnou(0))
			ENDIF
			z=zp-zmas*zmas
		ELSEIF(iwidth.EQ.0)THEN
			z=zp-parmas(m0)**2+DCMPLX(dnou(0),dnou(1))*parmas(m0)*parwid(m0)
		ENDIF
		z1=DCMPLX(8d0)*Helac_zprod(zp0,zeps1)*Helac_zprod(zp0,zeps2)
		z2=DCMPLX(2d0)*Helac_zprod(zeps2,zeps1)
		DO i=1,1
			zy(i)=-z2*zy(i)/z**2+z1*zy(i)/z**3
		ENDDO
		RETURN
	ENDIF
! F
	z=zp-parmas(m0)**2+DCMPLX(dnou(0),dnou(1))*parmas(m0)*parwid(m0)
	z1=DCMPLX(8d0)*Helac_zprod(zp0,zeps1)*Helac_zprod(zp0,zeps2)
	z2=DCMPLX(2d0)*Helac_zprod(zeps2,zeps1)
	DO i=1,4
		zy(i)=-z2*zy(i)/z**2+z1*zy(i)/z**3
	ENDDO
	RETURN
CASE DEFAULT
	PRINT *,"Wrong (deriv>2) of Helac_propagd in Helac_pan2.f90 ! STOP !"
	STOP
END SELECT
END SUBROUTINE Helac_propagd

FUNCTION Helac_zprod(za,zb)
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(IN)::za,zb
COMPLEX(KIND=DBL)::Helac_zprod
Helac_zprod=(za(1)*zb(2)+za(2)*zb(1)-za(3)*zb(4)-za(4)*zb(3) )/dnou(2)
END FUNCTION Helac_zprod
! from the derivation to eps(L_z)
!SUBROUTINE Deriv_eps(zeps,k,deriv)
!IMPLICIT NONE
!INTEGER,INTENT(IN)::deriv
!COMPLEX(KIND=DBL),DIMENSION(4),INTENT(OUT)::zeps
!INTEGER,DIMENSION(Pwavenum)::yyii
!INTEGER,INTENT(OUT)::k
!INTEGER::kk1,kk2
!COMPLEX(KIND=DBL),PARAMETER::zi=DCMPLX(0d0,1d0)
!REAL(KIND=DBL)::rr
!COMPLEX(KIND=DBL),DIMENSION(5)::zPPQQ
!INTEGER,DIMENSION(3)::index2
!COMPLEX(KIND=DBL),DIMENSION(3,4)::zeps1
!CALL Helac_bin(Pwavenum,deriv,yyii)
!DO kk1=1,Pwavenum
!	IF(yyii(Pwavenum-kk1+1).EQ.1)THEN
!		k=kk1
!		EXIT
!	ENDIF
!ENDDO
!kk1=Pwave(k)
!kk2=Quarkonium2(kk1)
!kk2=Quarkonium(kk2,2)  ! the number of first quark in parton array
!CALL PlusZq(zPPQQ(1:5),zq(kk2,1:5),zq(kk2+1,1:5))
!zPPQQ(2)=DCMPLX(DREAL(zPPQQ(2))*io(kk2),DIMAG(zPPQQ(2)))! p0-pz+i*pt  incoming -p0+pz+i*pt outgoing
!zPPQQ(1)=zPPQQ(1)*io(kk2)          ! po+pz+i*pz incoming -po-pz-i*pz outgoing
!zPPQQ(3)=zPPQQ(3)*io(kk2)         ! px+i*py incoming -px-i*py outgoing
!zPPQQ(4)=zPPQQ(4)*io(kk2)         ! px-i*py incoming -px+i*py outgoing       
!IF(iranhel.GE.2)THEN	
!	rr=Qhelran(kk1)
!	DO kk2=1,3
!		CALL PolarHadron(zeps1(kk2,1:4),zPPQQ,kk2)
!	ENDDO
!	CALL Helac_elong(zeps1(3,1:4),zPPQQ)
!	CALL Helac_eplus(zeps1(2,1:4),zPPQQ)
!	CALL Helac_eminus(zeps1(1,1:4),zPPQQ)
!	zeps(1:4)=(DCOS(rr)+zi*DSIN(rr))*zeps1(1,1:4)+(DCOS(rr)-zi*DSIN(rr))*zeps1(2,1:4) &
!	+zeps1(3,1:4)
!ELSE
!	CALL Qeveryhel(index2,ipol(kk1),kk1)
!	IF(index2(2).EQ.3)CALL Helac_elong(zeps,zPPQQ)
!	IF(index2(2).EQ.2)CALL Helac_eplus(zeps,zPPQQ)
!	IF(index2(2).EQ.1)CALL Helac_eminus(zeps,zPPQQ)
!ENDIF
!END SUBROUTINE Deriv_eps

SUBROUTINE Deriv_eps(zeps,k,deriv)
IMPLICIT NONE
INTEGER,INTENT(IN)::deriv
COMPLEX(KIND=DBL),DIMENSION(4),INTENT(OUT)::zeps
INTEGER,DIMENSION(Pwavenum)::yyii
INTEGER,INTENT(OUT)::k
INTEGER::kk1,kk2
COMPLEX(KIND=DBL),PARAMETER::zi=DCMPLX(0d0,1d0)
REAL(KIND=DBL)::rr,rr2
COMPLEX(KIND=DBL),DIMENSION(5)::zPPQQ
INTEGER,DIMENSION(3)::index2
COMPLEX(KIND=DBL),DIMENSION(3,4)::zeps1
COMPLEX(KIND=DBL),DIMENSION(1:3,1:7,4)::zepslist
INTEGER::Lz
SAVE zepslist
IF(GLOBALINIT_mader.EQ.0)THEN
	zepslist(1:3,1:7,1:4)=0
	GLOBALINIT_mader=1
ENDIF
CALL Helac_bin(Pwavenum,deriv,yyii)
DO kk1=1,Pwavenum
	IF(yyii(Pwavenum-kk1+1).EQ.1)THEN
		k=kk1
		EXIT
	ENDIF
ENDDO
kk1=Pwave(k)
IF(iranhel.GE.2)THEN
	rr2=ABS(zepslist(1,deriv,1))+ABS(zepslist(1,deriv,2))&
	+ABS(zepslist(1,deriv,3))+ABS(zepslist(1,deriv,4))
	IF(rr2.GT.0)THEN
		zeps(1:4)=zepslist(1,deriv,1:4)
		RETURN
	ENDIF
ELSE
	CALL Qeveryhel(index2,ipol(kk1),kk1)
	Lz=index2(2)
	rr2=ABS(zepslist(Lz,deriv,1))+ABS(zepslist(Lz,deriv,2))&
	+ABS(zepslist(Lz,deriv,3))+ABS(zepslist(Lz,deriv,4))
	IF(rr2.GT.0)THEN
		zeps(1:4)=zepslist(Lz,deriv,1:4)
		RETURN
	ENDIF
ENDIF
kk2=Quarkonium2(kk1)
kk2=Quarkonium(kk2,2)  ! the number of first quark in parton array
CALL PlusZq(zPPQQ(1:5),zq(kk2,1:5),zq(kk2+1,1:5))
!IF(PolarFrame.EQ.1)THEN
!STOP
zPPQQ(2)=DCMPLX(DREAL(zPPQQ(2))*io(kk2),DIMAG(zPPQQ(2)))! p0-pz+i*pt  incoming -p0+pz+i*pt outgoing
zPPQQ(1)=zPPQQ(1)*io(kk2)          ! po+pz+i*pz incoming -po-pz-i*pz outgoing
zPPQQ(3)=zPPQQ(3)*io(kk2)         ! px+i*py incoming -px-i*py outgoing
zPPQQ(4)=zPPQQ(4)*io(kk2)         ! px-i*py incoming -px+i*py outgoing
!ENDIF     
IF(iranhel.GE.2)THEN	
	rr=Qhelran(kk1)
	IF(imode.EQ.1)THEN
		DO kk2=1,3
			CALL PolarHadron(zeps1(kk2,1:4),zPPQQ,kk2)
		ENDDO
	ELSE
		CALL Helac_elong(zeps1(3,1:4),zPPQQ)
		CALL Helac_eplus(zeps1(2,1:4),zPPQQ)
		CALL Helac_eminus(zeps1(1,1:4),zPPQQ)
	ENDIF
	zeps(1:4)=(DCOS(rr)+zi*DSIN(rr))*zeps1(1,1:4)+(DCOS(rr)-zi*DSIN(rr))*zeps1(2,1:4) &
	+zeps1(3,1:4)
	zepslist(1,deriv,1:4)=zeps(1:4)
ELSE
	IF(imode.EQ.1)THEN
		CALL PolarHadron(zeps,zPPQQ,lz)
	ELSE
		IF(Lz.EQ.3)CALL Helac_elong(zeps,zPPQQ)
		IF(Lz.EQ.2)CALL Helac_eplus(zeps,zPPQQ)
		IF(Lz.EQ.1)CALL Helac_eminus(zeps,zPPQQ)
	ENDIF
	zepslist(Lz,deriv,1:4)=zeps(1:4)
ENDIF
RETURN
END SUBROUTINE Deriv_eps
END MODULE Helac_pan2