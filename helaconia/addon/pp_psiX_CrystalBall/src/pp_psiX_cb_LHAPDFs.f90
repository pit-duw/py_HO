SUBROUTINE strf_lhapdf_pp_psiX_cb(wsf)
USE Helac_Global
USE pp_psiX_cb_global
!INTEGER,INTENT(IN)::ipip
INTEGER::ipp=1  ! ipp=1 ppbar ; ipp=2 pp
INTEGER::ih=1   ! ih=1 no photon PDF,ih=2, photon from proton/anti-proton,ih=3 photon from electron/positron 
REAL(KIND=DBL),INTENT(OUT)::wsf          
CHARACTER(len=20),DIMENSION(20)::parm
REAL(KIND=DBL),DIMENSION(20)::val
REAL(KIND=DBL),DIMENSION(-7:7)::pdflist
INTEGER::init=0      
SAVE init,ipp,ih
REAL(KIND=DBL)::glu1_ct,glu2_ct,u1_ct,u2_ct,d1_ct,d2_ct,s1_ct,s2_ct,c1_ct,c2_ct,b1_ct,b2_ct,&
                ub1_ct,ub2_ct,db1_ct,db2_ct,sb1_ct,sb2_ct,cb1_ct,cb2_ct,bb1_ct,bb2_ct,sf_ct
       
IF(init.EQ.0)THEN
   CALL lhapdfwrap
   IF(Coll_Type.EQ.1)THEN
      ipp=2
   ELSE
      ipp=1
   ENDIF
   init=1
ENDIF

scale=cb_pmom(3,1)**2+cb_pmom(3,2)**2+mpsi**2
scale=DSQRT(scale)

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

wsf=glu1_ct*glu2_ct*wjac/xp1/xp2
IF(cb_includeqq)THEN
   wsf=wsf+(u1_ct*ub2_ct+d1_ct*db2_ct+s1_ct*sb2_ct+c1_ct*cb2_ct+&
        ub1_ct*u2_ct+db1_ct*d2_ct+sb1_ct*s2_ct+cb1_ct*c2_ct)*wjac/xp1/xp2
ENDIF

IF(init.EQ.0)init=1
       
END SUBROUTINE strf_lhapdf_pp_psiX_cb
