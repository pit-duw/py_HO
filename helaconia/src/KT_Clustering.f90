MODULE KT_Clustering
! MLM-scheme matching, refer to arXiv:0706.2569
USE Helac_Global
USE Structf_PDFs
USE Helac_Func_1
IMPLICIT NONE
INTEGER,DIMENSION(20)::ifl_cl
INTEGER,DIMENSION(20,2)::icol_cl
REAL(KIND=DBL),DIMENSION(20,4)::pcl
CONTAINS       
SUBROUTINE ktclus(kt_pmom)
!INTEGER,INTENT(IN)::n1
REAL(KIND=DBL),DIMENSION(1:n,1:4),INTENT(IN)::kt_pmom
INTEGER::i,j,np,ifound,i0,j0
REAL(KIND=DBL)::dmin,d
!     common/clus/io(1:n),icol_cl(20,2),pcl(20,4)
DO i=1,nhad
   IF(ABS(iflh(i)).GT.100)THEN
		pcl(i,1:4)=kt_pmom(hadron2parton(i),1:4)+kt_pmom(hadron2parton(i)+1,1:4)
   ELSE
		pcl(i,1:4)=kt_pmom(hadron2parton(i),1:4)
   ENDIF
ENDDO
np=nhad
! 1    continue
DO
  ifound=0
  i0=0
  j0=0
  dmin=1d13
  IF(np.EQ.0)RETURN

  DO i=1,np-1
     DO j=i+1,np
        d=kt_dist(i,j)
        IF(d.LT.dmin)THEN
          ifound=1
          i0=i
          j0=j
          dmin=d
        ENDIF
     ENDDO
  ENDDO

  IF(ifound.EQ.1)THEN
     CALL kt_combine(i0,j0)
     np=np-1
! call weight       
     IF(dmin.LT.0)WRITE(*,*)dmin,&   !w,
	                        np,i0,j0
!       goto 1
  ELSE
     EXIT
  ENDIF
ENDDO

END SUBROUTINE ktclus

FUNCTION kt_dist(i,j)
IMPLICIT NONE
INTEGER,INTENT(IN)::i,j
REAL(KIND=DBL)::kt_dist
REAL(KIND=DBL)::pt1,pt2,y1,y2,dphi

pt1=pcl(i,1)**2+pcl(i,2)**2
pt2=pcl(j,1)**2+pcl(j,2)**2
! Eq.(7) in Nucl.Phys.B 406 (1993) 187-224
! Eqs.(9,10) in 0706.2569 
! initial-final
IF(ioh(i).EQ. 1.AND.ioh(j).EQ.-1)THEN
    kt_dist=pt2
ENDIF
! final-initial
IF(ioh(i).EQ.-1.AND.ioh(j).EQ. 1)THEN
    kt_dist=pt1
ENDIF
! final-final
IF(ioh(i).EQ.-1.AND.ioh(j).EQ.-1)THEN
    y1=-0.5d0*DLOG((pcl(i,4)-pcl(i,3))/(pcl(i,4)+pcl(i,3))) ! rapidity
    y2=-0.5d0*dlog((pcl(j,4)-pcl(j,3))/(pcl(j,4)+pcl(j,3)))
    dphi=(pcl(i,1)*pcl(j,1)+pcl(i,2)*pcl(j,2))/DSQRT(pt1*pt2)
    IF(dphi.LT.-1d0.AND.dphi.GT.-1d0-1d-6)dphi=-1.d0
    dphi=DACOS(dphi)
    kt_dist=MIN(pt1,pt2)*((y1-y2)**2+dphi**2)
    IF(kt_dist.LT.0)WRITE(*,*)i,j,y1,y2,dphi,pt1,pt2
    IF(kt_dist.LT.0)WRITE(*,*)i,pcl(i,1:4)
    IF(kt_dist.LT.0)WRITE(*,*)j,pcl(j,1:4)
ENDIF
! initial-initial
IF(ioh(i).EQ. 1.AND.ioh(j).EQ. 1)THEN
    kt_dist=1d14
ENDIF
END FUNCTION kt_dist

! Covariant E-scheme
SUBROUTINE kt_combine(i,j)
IMPLICIT NONE
REAL(KIND=DBL),DIMENSION(4)::p
INTEGER,INTENT(IN)::i,j
INTEGER::i1,j1

i1=MIN(i,j)
j1=MAX(i,j)

p(1:4)=ioh(i)*pcl(i,1:4)+ioh(j)*pcl(j,1:4)
IF(p(4).LT.0)p(1:4)=-p(1:4)

pcl(i1,1:4)=p(1:4)
pcl(j1:nhad-1,1:4)=pcl(j1+1:nhad,1:4)

END SUBROUTINE kt_combine

!-----------------------------------------------------------       
SUBROUTINE ktreweight(w)
REAL(KIND=DBL),INTENT(OUT)::w
INTEGER::i,j,nino,np,npo,ifound,i0,j0,if1,if2
LOGICAL,DIMENSION(20)::ino
REAL(KIND=DBL)::dmin,d,Q,Q1,Q2

ino(1:20)=.TRUE.
DO i=1,nhad
   IF(ABS(iflh(i)).GT.100)THEN
		pcl(i,1:4)=phegas_pmom(hadron2parton(i),1:4)+phegas_pmom(hadron2parton(i)+1,1:4)
		ifl_cl(i)=iflh(i)
   ELSE
		pcl(i,1:4)=phegas_pmom(hadron2parton(i),1:4)
		ifl_cl(i)=MOD(iflh(i),35)
   ENDIF
   icol_cl(i,1:2)=icol_un(i,1:2)
   ! only the QCD partons
   IF(ABS(iflh(i)).EQ.3)ino(i)=.FALSE.
   IF(ABS(iflh(i)).EQ.4)ino(i)=.FALSE.
   IF(ABS(iflh(i)).EQ.7)ino(i)=.FALSE.
   IF(ABS(iflh(i)).EQ.8)ino(i)=.FALSE.
   IF(ABS(iflh(i)).EQ.11)ino(i)=.FALSE.
   IF(ABS(iflh(i)).EQ.12)ino(i)=.FALSE.
   IF(ABS(iflh(i)).EQ.35)ino(i)=.FALSE.
ENDDO

nino=0
DO i=1,nhad
   IF(ino(i))nino=nino+1
ENDDO

np=nhad
w=1.0
Q1=1d14
Q2=1d14
DO      
! 1    continue
  ifound=0
  i0=0
  j0=0
  dmin=1d13
  IF(np.EQ.0)RETURN

  DO i=1,np-1
     IF(ino(i))CYCLE
     if1=ifl_cl(i)*ioh(i)
     DO j=i+1,np
        IF(ino(j))CYCLE
        if2=ifl_cl(j)*ioh(j)
!  combination only relating to gluon to quarks or incoming/outgoing quark 
!  to incoming/outgoing antiquark ( outgoing/incoming quark) 
        IF((if1+if2)*if1*if2.NE.0)CYCLE
! check color clusterability
        IF(icol_cl(i,1).NE.icol_cl(j,2).AND.icol_cl(i,2).NE.icol_cl(j,1))CYCLE
        d=kt_dist_2(i,j)
        IF(d.LT.dmin)THEN
           ifound=1
           i0=i
           j0=j
           dmin=d
        ENDIF
    ENDDO
  ENDDO
! combine the smallest one first
! np=4 is the core process
  IF(ifound.EQ.1.AND.np.GT.4)THEN
    CALL kt_combine_2(i0,j0)
    np=np-1

    Q=DSQRT(dmin)
    ! alpha_S reweighting
    w=w*alphas_clu(Q)/aqcdup   !0.130d0 !alphas(scale)
    ! PDF reweighting
    IF(istruc.EQ.1.AND.(i0.LE.2.AND.j0.GE.3).OR.(i0.GE.3.AND.j0.LE.2))THEN
       IF(MIN(i0,j0).EQ.1.AND.Q1.GT.Q)Q1=Q
       IF(MIN(i0,j0).EQ.2.AND.Q2.GT.Q)Q2=Q
    ENDIF

    IF(dmin.LT.0)WRITE(*,*)dmin,w,np,i0,j0
!       goto 1
  ELSE
     EXIT
  ENDIF
ENDDO
! no PDF reweighting for MLM
! PDF reweighting
!IF(istruc.EQ.1)CALL PDFreweight(Q1,Q2,w)
! no reweigting again since the remaing vertices uses hard scale
! which is the default one aqcdup
!CALL qcdscale(Q)

! the number of final-state jets
!npo=np-nino-2
!IF(npo.GT.0)w=w*(alphas_clu(Q)/aqcdup)**(npo) !0.130d0 !alphas(scale)

END SUBROUTINE ktreweight
 
FUNCTION kt_dist_2(i,j)
IMPLICIT NONE
INTEGER,INTENT(IN)::i,j
REAL(KIND=DBL)::kt_dist_2,pt1,pt2,y1,y2,p1,p2,dphi,pttemp

pt1=pcl(i,1)**2+pcl(i,2)**2
pt2=pcl(j,1)**2+pcl(j,2)**2
   
! initial-final
IF(ioh(i).EQ. 1.AND.ioh(j).EQ.-1)THEN
   kt_dist_2=pt2
   pttemp=ptc(hadron2parton(j))**2
   IF(ABS(iflh(j)).GT.100)THEN
      pttemp=(ptc(hadron2parton(j))+ptc(hadron2parton(j)+1))**2
   ENDIF
   IF(kt_dist_2.LT.pttemp)kt_dist_2=1d14 !ptc(j)**2
ENDIF
! final-initial
IF(ioh(i).EQ.-1.AND.ioh(j).EQ. 1)THEN
   kt_dist_2=pt1
   pttemp=ptc(hadron2parton(i))**2
   IF(ABS(iflh(i)).GT.100)THEN
      pttemp=(ptc(hadron2parton(i))+ptc(hadron2parton(i)+1))**2
   ENDIF
   IF(kt_dist_2.LT.pttemp)kt_dist_2=1d14 !ptc(i)**2
ENDIF
! final-final
IF(ioh(i).EQ.-1.AND.ioh(j).EQ.-1)THEN

   p1=DSQRT(pcl(i,1)**2+pcl(i,2)**2+pcl(i,3)**2)
   y1=-0.5d0*DLOG((p1-pcl(i,3))/(p1+pcl(i,3)))  ! pesudorapdity
   p2=DSQRT(pcl(j,1)**2+pcl(j,2)**2+pcl(j,3)**2)
   y2=-0.5d0*DLOG((p2-pcl(j,3))/(p2+pcl(j,3)))

   dphi=(pcl(i,1)*pcl(j,1)+pcl(i,2)*pcl(j,2))/DSQRT(pt1*pt2)
   IF(dphi.LT.-1d0.AND.dphi.GT.-1d0-1d-6)dphi=-1.d0
   dphi=DACOS(dphi)


   IF(ktmeasure.EQ.1)kt_dist_2=MIN(pt1,pt2)*((y1-y2)**2+dphi**2)
   IF(ktmeasure.EQ.2)kt_dist_2=2*MIN(pt1,pt2)*(DCOSH(y1-y2)-DCOS(dphi))
   kt_dist_2=kt_dist_2/RR**2
   pttemp=ptc(hadron2parton(i))**2
   IF(ABS(iflh(i)).GT.100)THEN
      pttemp=(ptc(hadron2parton(i))+ptc(hadron2parton(i)+1))**2
   ENDIF        
   IF(kt_dist_2.LT.drc(hadron2parton(i),hadron2parton(j))**2*pttemp)kt_dist_2=1d14 !drc(i,j)**2*ptc(i)**2
   IF(kt_dist_2.LT.0)WRITE(*,*)i,j,y1,y2,dphi,pt1,pt2
   IF(kt_dist_2.LT.0)WRITE(*,*)i,pcl(i,1:4)
   IF(kt_dist_2.LT.0)WRITE(*,*)j,pcl(j,1:4)
   IF(kt_dist_2.LT.0)WRITE(*,*)ioh(1:nhad)
ENDIF
! initial-initial
IF(ioh(i).EQ. 1.AND.ioh(j).EQ. 1)THEN
    kt_dist_2=1d14
ENDIF
END FUNCTION kt_dist_2

SUBROUTINE kt_combine_2(i,j)
IMPLICIT NONE
INTEGER,INTENT(IN)::i,j
REAL(KIND=DBL),DIMENSION(4)::p
INTEGER::i1,j1,if1,if2,if0

i1=MIN(i,j)
j1=MAX(i,j)

p(1:4)=ioh(i)*pcl(i,1:4)+ioh(j)*pcl(j,1:4)
IF(p(4).LT.0)p(1:4)=-p(1:4)

pcl(i1,1:4)=p(1:4)
pcl(j1:nhad-1,1:4)=pcl(j1+1:nhad,1:4)
! ajust the color flow
IF(icol_cl(i1,1).EQ.icol_cl(j1,2))THEN
   icol_cl(i1,1)=icol_cl(j1,1)
ELSEIF(icol_cl(i1,2).EQ.icol_cl(j1,1))THEN
   icol_cl(i1,2)=icol_cl(j1,2)
ENDIF
IF(icol_cl(i1,1).EQ.icol_cl(i1,2))icol_cl(i1,1:2)=0
icol_cl(j1:nhad-1,1:2)=icol_cl(j1+1:nhad,1:2)

if1=ifl_cl(i1)*ioh(i1)
if2=ifl_cl(j1)*ioh(j1)

IF(if1+if2.EQ.0)if0=0
IF(if1.EQ.0)if0=if2
IF(if2.EQ.0)if0=if1

ifl_cl(i1)=if0
ifl_cl(j1:nhad-1)=ifl_cl(j1+1:nhad)

END SUBROUTINE kt_combine_2

!-----------------------------------------------------------------------      
! running alphas in the clustering procedure      
FUNCTION ALPHAS_clu(Q)
IMPLICIT NONE
!      include 'alphas.h' 
REAL(KIND=DBL),INTENT(IN)::Q
REAL(KIND=DBL)::Q2
REAL(KIND=DBL)::ALPHAS_clu
alphas_clu=ALPHAS(Q)
RETURN
Q2=Q**2
alphas_clu=alfas(Q2,0.165d0,1,5)

END FUNCTION ALPHAS_clu

!-------------------------------------------------------------------
!------- ALPHA QCD -------------------------------------
! Program to calculate alfa strong with nf flavours,
! as a function of lambda with 5 flavors.
! The value of alfa is matched at the thresholds q = mq.
! When invoked with nf < 0 it chooses nf as the number of
! flavors with mass less than q.
!
FUNCTION alfas(q2,xlam,nloop,inf)
IMPLICIT NONE
REAL(KIND=DBL)::olam=0.d0,pi=3.14159d0,xmb=4.5d0,xmc=1.5d0
INTEGER,INTENT(IN)::nloop,inf
REAL(KIND=DBL),INTENT(IN)::xlam,q2
REAL(KIND=DBL)::alfas
REAL(KIND=DBL)::b3,b4,b5,bp5,bp4,bp3,xlc,xlb,xllc,xllb,c45,c35,q,xlq,xllq
INTEGER::nf
SAVE
IF(xlam.NE.olam)THEN
    olam = xlam
    b5  = (33-2*5)/pi/12   ! beta0/4/pi
    b4  = (33-2*4)/pi/12
    b3  = (33-2*3)/pi/12
    IF(nloop.EQ.1)THEN
       bp5=0
       bp4=0
       bp3=0
    ELSEIF(nloop.EQ.2)THEN
       bp5 = (153 - 19*5) / pi / 2 / (33 - 2*5)   ! beta1/beta0/4/pi
       bp4 = (153 - 19*4) / pi / 2 / (33 - 2*4)
       bp3 = (153 - 19*3) / pi / 2 / (33 - 2*3)
    ENDIF
    xlc = 2 * LOG(xmc/xlam)
    xlb = 2 * LOG(xmb/xlam)
    xllc = LOG(xlc)
    xllb = LOG(xlb)
    c45  =  1/( 1/(b5 * xlb) - xllb*bp5/(b5 * xlb)**2 ) &
            - 1/( 1/(b4 * xlb) - xllb*bp4/(b4 * xlb)**2 )
    c35  =  1/( 1/(b4 * xlc) - xllc*bp4/(b4 * xlc)**2 ) &
           - 1/( 1/(b3 * xlc) - xllc*bp3/(b3 * xlc)**2 ) + c45
ENDIF
q   = SQRT(q2)
xlq = 2 * LOG( q/xlam )
xllq = LOG( xlq )
nf = inf
IF( nf .LT. 0) THEN
    IF( q .GT. xmb ) THEN
       nf = 5
    ELSEIF( q .GT. xmc ) THEN
       nf = 4
    ELSE
       nf = 3
    ENDIF
ENDIF
IF( nf .EQ. 5 ) THEN
    alfas = 1/(b5 * xlq) -  bp5/(b5 * xlq)**2 * xllq
ELSEIF( nf .EQ. 4 ) THEN
    alfas = 1/( 1/(1/(b4 * xlq) - bp4/(b4 * xlq)**2 * xllq) + c45 )
ELSEIF( nf .EQ. 3 ) THEN
    alfas = 1/( 1/(1/(b3 * xlq) - bp3/(b3 * xlq)**2 * xllq) + c35 )
ELSE
    WRITE(*,*)'error in alfa: unimplemented # of light flavours',nf
    STOP
ENDIF
END FUNCTION alfas
!-------------------------------------------
! Program to calculate as with nf flavours
! as a function of lambda with nf flavours
!
FUNCTION alfa(q,xlam,nloop,nf)
REAL(KIND=DBL),INTENT(IN)::q,xlam
INTEGER,INTENT(IN)::nloop,nf
REAL(KIND=DBL)::xlp,b,bp,t,xlt
REAL(KIND=DBL)::pi=3.14159d0
REAL(KIND=DBL)::alfa
xlp = FLOAT(nloop-1)
b  = (33-2*nf)/pi/12
bp = (153 - 19*nf) / pi / 2 / (33 - 2*nf)  * xlp
t = 2 * LOG( q/xlam )
xlt = LOG( t )
alfa = 1/(b * t) -  bp/(b * t)**2 * xlt
END FUNCTION alfa

!----------------------------------------------------------
! Program to get lambda_nf from as_nf at the scale q
!
FUNCTION xlambd(as,q,nloop,nf)
REAL(KIND=DBL),INTENT(IN)::as,q
INTEGER,INTENT(IN)::nloop,nf
REAL(KIND=DBL)::xlambd
REAL(KIND=DBL)::pi=3.14159d0
REAL(KIND=DBL)::xlp,b,bp,t,xlt,ot,as0,as1
xlp = FLOAT(nloop-1)
b  = (33-2*nf)/pi/12
bp = (153 - 19*nf) / pi / 2 / (33 - 2*nf) * xlp
t  = 1/b/as
!-----------------------------------------------------------
! Solve the equation
!
DO
  xlt = LOG(t)
  ot = t
!-----------------------------------------------------------
! Solve the equation
! Value and Derivative of alfa with respect to t
!
  as0  = 1/b/t - bp*xlt/(b*t)**2
  as1  = - 1/b/t**2 -bp/b**2*(1-2*xlt)/t**3
  t  = (as-as0)/as1 + t
  IF(ABS(ot-t)/ot.LE..00001d0)EXIT
ENDDO
xlambd = q/EXP(t/2)
END FUNCTION xlambd

SUBROUTINE PDFreweight(Q1,Q2,w)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                          !
!     14.Aug.2012                                          !
!     Hua-Sheng Shao                                       !
!                                                          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL(KIND=DBL),INTENT(IN)::Q1,Q2
REAL(KIND=DBL),INTENT(INOUT)::w
INTEGER::ipp=1  ! ipp=1 ppbar ; ipp=2 pp
INTEGER::ih=1   ! ih=1 no photon PDF, ih=2 photon from proton/antiproton,ih=3 photon from electron/positron         
CHARACTER(len=20),DIMENSION(20)::parm
REAL(KIND=DBL),DIMENSION(20)::val
REAL(KIND=DBL),DIMENSION(-7:7)::pdflist
INTEGER::init=0
LOGICAL::use_cteq6_f90=.TRUE.      
SAVE init,ipp,ih,use_cteq6_f90
REAL(KIND=DBL)::glu11_ct,glu12_ct,u11_ct,u12_ct,d11_ct,d12_ct,s11_ct,s12_ct,c11_ct,c12_ct,&
				b11_ct,b12_ct,ub11_ct,ub12_ct,db11_ct,db12_ct,sb11_ct,sb12_ct,&
				cb11_ct,cb12_ct,bb11_ct,bb12_ct,sf_ct,glu21_ct,glu22_ct,u21_ct,&
				u22_ct,d21_ct,d22_ct,s21_ct,s22_ct,c21_ct,c22_ct,&
				b21_ct,b22_ct,ub21_ct,ub22_ct,db21_ct,db22_ct,sb21_ct,sb22_ct,&
				cb21_ct,cb22_ct,bb21_ct,bb22_ct,QQ1,QQ2       
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
QQ1=MIN(Q1,scale)
QQ2=MIN(Q2,scale)

IF(use_cteq6_f90)THEN
   glu11_ct = Ctq6Pdf_f90(0,xp1,QQ1)
   glu12_ct = Ctq6Pdf_f90(0,xp1,scale)
   glu21_ct = Ctq6Pdf_f90(0,xp2,QQ2)
   glu22_ct = Ctq6Pdf_f90(0,xp2,scale)
   u11_ct   = Ctq6Pdf_f90(1,xp1,QQ1)
   u12_ct   = Ctq6Pdf_f90(1,xp1,scale)
   u21_ct   = Ctq6Pdf_f90(1,xp2,QQ2)
   u22_ct   = Ctq6Pdf_f90(1,xp2,scale)
   d11_ct   = Ctq6Pdf_f90(2,xp1,QQ1)
   d12_ct   = Ctq6Pdf_f90(2,xp1,scale)
   d21_ct   = Ctq6Pdf_f90(2,xp2,QQ2)
   d22_ct   = Ctq6Pdf_f90(2,xp2,scale)
   s11_ct   = Ctq6Pdf_f90(3,xp1,QQ1)
   s12_ct   = Ctq6Pdf_f90(3,xp1,scale)
   s21_ct   = Ctq6Pdf_f90(3,xp2,QQ2)
   s22_ct   = Ctq6Pdf_f90(3,xp2,scale)
   c11_ct   = Ctq6Pdf_f90(4,xp1,QQ1)
   c12_ct   = Ctq6Pdf_f90(4,xp1,scale)
   c21_ct   = Ctq6Pdf_f90(4,xp2,QQ2)
   c22_ct   = Ctq6Pdf_f90(4,xp2,scale)
   b11_ct   = Ctq6Pdf_f90(5,xp1,QQ1)
   b12_ct   = Ctq6Pdf_f90(5,xp1,scale)
   b21_ct   = Ctq6Pdf_f90(5,xp2,QQ2)
   b22_ct   = Ctq6Pdf_f90(5,xp2,scale)
   ub11_ct  = Ctq6Pdf_f90(-1,xp1,QQ1)
   ub12_ct  = Ctq6Pdf_f90(-1,xp1,scale)
   ub21_ct  = Ctq6Pdf_f90(-1,xp2,QQ2)
   ub22_ct  = Ctq6Pdf_f90(-1,xp2,scale)
   db11_ct  = Ctq6Pdf_f90(-2,xp1,QQ1)
   db12_ct  = Ctq6Pdf_f90(-2,xp1,scale)
   db21_ct  = Ctq6Pdf_f90(-2,xp2,QQ2)
   db22_ct  = Ctq6Pdf_f90(-2,xp2,scale)
   sb11_ct  = Ctq6Pdf_f90(-3,xp1,QQ1)
   sb12_ct  = Ctq6Pdf_f90(-3,xp1,scale)
   sb21_ct  = Ctq6Pdf_f90(-3,xp2,QQ2)
   sb22_ct  = Ctq6Pdf_f90(-3,xp2,scale)
   cb11_ct  = Ctq6Pdf_f90(-4,xp1,QQ1)
   cb12_ct  = Ctq6Pdf_f90(-4,xp1,scale)
   cb21_ct  = Ctq6Pdf_f90(-4,xp2,QQ2)
   cb22_ct  = Ctq6Pdf_f90(-4,xp2,scale)
   bb11_ct  = Ctq6Pdf_f90(-5,xp1,QQ1)
   bb12_ct  = Ctq6Pdf_f90(-5,xp1,scale)
   bb21_ct  = Ctq6Pdf_f90(-5,xp2,QQ2)
   bb22_ct  = Ctq6Pdf_f90(-5,xp2,scale)
ELSE
   CALL fdist(ih,xp1,scale,pdflist(-7:7))
   glu12_ct = pdflist(0)
   u12_ct   = pdflist(2)
   d12_ct   = pdflist(1)
   s12_ct   = pdflist(3)
   c12_ct   = pdflist(4)
   b12_ct   = pdflist(5)
   ub12_ct  = pdflist(-2)
   db12_ct  = pdflist(-1)
   sb12_ct  = pdflist(-3)
   cb12_ct  = pdflist(-4)
   bb12_ct  = pdflist(-5)
   IF(QQ1.NE.scale)THEN
      CALL fdist(ih,xp1,QQ1,pdflist(-7:7))
   ENDIF
   glu11_ct = pdflist(0)
   u11_ct   = pdflist(2)
   d11_ct   = pdflist(1)
   s11_ct   = pdflist(3)
   c11_ct   = pdflist(4)
   b11_ct   = pdflist(5)
   ub11_ct  = pdflist(-2)
   db11_ct  = pdflist(-1)
   sb11_ct  = pdflist(-3)
   cb11_ct  = pdflist(-4)
   bb11_ct  = pdflist(-5)
   CALL fdist(ih,xp2,scale,pdflist(-7:7))
   glu22_ct = pdflist(0)
   u22_ct   = pdflist(2)
   d22_ct   = pdflist(1)
   s22_ct   = pdflist(3)
   c22_ct   = pdflist(4)
   b22_ct   = pdflist(5)
   ub22_ct  = pdflist(-2)
   db22_ct  = pdflist(-1)
   sb22_ct  = pdflist(-3)
   cb22_ct  = pdflist(-4)
   bb22_ct  = pdflist(-5)
   IF(QQ2.NE.scale)THEN
      CALL fdist(ih,xp2,QQ2,pdflist(-7:7))
   ENDIF
   glu21_ct = pdflist(0)
   u21_ct   = pdflist(2)
   d21_ct   = pdflist(1)
   s21_ct   = pdflist(3)
   c21_ct   = pdflist(4)
   b21_ct   = pdflist(5)
   ub21_ct  = pdflist(-2)
   db21_ct  = pdflist(-1)
   sb21_ct  = pdflist(-3)
   cb21_ct  = pdflist(-4)
   bb12_ct  = pdflist(-5)
ENDIF

sf_ct=1
IF(qsumQ)THEN
	IF(ifl(1).EQ.35)THEN
		sf_ct=sf_ct*glu11_ct/glu12_ct
	ELSEIF(ABS(ifl(1)).LE.12.AND.(MOD(ABS(ifl(1)),4).EQ.0.OR.MOD(ABS(ifl(1)),4).EQ.3) &
	.AND.ifl(2).EQ.35)THEN
		IF(onlyqcd)THEN
			SELECT CASE(iqnum)
			CASE(1)
				sf_ct=sf_ct*(u11_ct+ub11_ct)/(u12_ct+ub12_ct)
			CASE(2)
				sf_ct=sf_ct*(u11_ct+ub11_ct+d11_ct+db11_ct)/(u12_ct+ub12_ct+d12_ct+db12_ct)
			CASE(3)
				sf_ct=sf_ct*(u11_ct+ub11_ct+d11_ct+db11_ct+s11_ct+sb11_ct)&
				/(u12_ct+ub12_ct+d12_ct+db12_ct+s12_ct+sb12_ct)
			CASE(4)
				sf_ct=sf_ct*(c11_ct+cb11_ct+u11_ct+ub11_ct+d11_ct+db11_ct+s11_ct+sb11_ct)&
				/(c12_ct+cb12_ct+u12_ct+ub12_ct+d12_ct+db12_ct+s12_ct+sb12_ct)
			CASE DEFAULT
				sf_ct=sf_ct*(b11_ct+bb11_ct+c11_ct+cb11_ct+u11_ct+ub11_ct+d11_ct+db11_ct+s11_ct+sb11_ct)&
				/(b12_ct+bb12_ct+c12_ct+cb12_ct+u12_ct+ub12_ct+d12_ct+db12_ct+s12_ct+sb12_ct)
			END SELECT
		ELSE
			IF(MOD(ABS(ifl(1)),4).EQ.3)THEN
				SELECT CASE(iqnum)
				CASE(1,2,3)
					sf_ct=sf_ct*(u11_ct+ub11_ct)/(u12_ct+ub12_ct)
				CASE DEFAULT
					sf_ct=sf_ct*(u11_ct+ub11_ct+c11_ct+cb11_ct)/(u12_ct+ub12_ct+c12_ct+cb12_ct)
				END SELECT
			ELSEIF(MOD(ABS(ifl(1)),4).EQ.0)THEN
				SELECT CASE(iqnum)
				CASE(1)
					sf_ct=sf_ct
				CASE(2)
					sf_ct=sf_ct*(d11_ct+db11_ct)/(d12_ct+db12_ct)
				CASE(3,4)
					sf_ct=sf_ct*(d11_ct+db11_ct+s11_ct+sb11_ct)/(d12_ct+db12_ct+s12_ct+sb12_ct)
				CASE DEFAULT
					sf_ct=sf_ct*(d11_ct+db11_ct+s11_ct+sb11_ct+b11_ct+bb11_ct)&
					/(d12_ct+db12_ct+s12_ct+sb12_ct+b12_ct+bb12_ct)
				END SELECT
			ENDIF
		ENDIF
	ENDIF

	IF(ifl(2).EQ.35)THEN
		sf_ct=sf_ct*glu21_ct/glu22_ct
	ELSEIF(ABS(ifl(2)).LE.12.AND.(MOD(ABS(ifl(2)),4).EQ.0.OR.MOD(ABS(ifl(2)),4).EQ.3) &
	.AND.ifl(1).EQ.35)THEN
		IF(onlyqcd)THEN
			SELECT CASE(iqnum)
			CASE(1)
				sf_ct=sf_ct*(u21_ct+ub21_ct)/(u22_ct+ub22_ct)
			CASE(2)
				sf_ct=sf_ct*(u21_ct+ub21_ct+d21_ct+db21_ct)/(u22_ct+ub22_ct+d22_ct+db22_ct)
			CASE(3)
				sf_ct=sf_ct*(u21_ct+ub21_ct+d21_ct+db21_ct+s21_ct+sb21_ct)&
				/(u22_ct+ub22_ct+d22_ct+db22_ct+s22_ct+sb22_ct)
			CASE(4)
				sf_ct=sf_ct*(c21_ct+cb21_ct+u21_ct+ub21_ct+d21_ct+db21_ct+s21_ct+sb21_ct)&
				/(c22_ct+cb22_ct+u22_ct+ub22_ct+d22_ct+db22_ct+s22_ct+sb22_ct)
			CASE DEFAULT
				sf_ct=sf_ct*(b21_ct+bb21_ct+c21_ct+cb21_ct+u21_ct+ub21_ct+d21_ct+db21_ct+s21_ct+sb21_ct)&
				/(b22_ct+bb22_ct+c22_ct+cb22_ct+u22_ct+ub22_ct+d22_ct+db22_ct+s22_ct+sb22_ct)
			END SELECT
		ELSE
			IF(MOD(ABS(ifl(2)),4).EQ.3)THEN
				SELECT CASE(iqnum)
				CASE(1,2,3)
					sf_ct=sf_ct*(u21_ct+ub21_ct)/(u22_ct+ub22_ct)
				CASE DEFAULT
					sf_ct=sf_ct*(u21_ct+ub21_ct+c21_ct+cb21_ct)/(u22_ct+ub22_ct+c22_ct+cb22_ct)
				END SELECT
			ELSEIF(MOD(ABS(ifl(2)),4).EQ.0)THEN
				SELECT CASE(iqnum)
				CASE(1)
					sf_ct=sf_ct
				CASE(2)
					sf_ct=sf_ct*(d21_ct+db21_ct)/(d22_ct+db22_ct)
				CASE(3,4)
					sf_ct=sf_ct*(d21_ct+db21_ct+s21_ct+sb21_ct)/(d22_ct+db22_ct+s22_ct+sb22_ct)
				CASE DEFAULT
					sf_ct=sf_ct*(d21_ct+db21_ct+s21_ct+sb21_ct+b21_ct+bb21_ct)&
					/(d22_ct+db22_ct+s22_ct+sb22_ct+b22_ct+bb22_ct)
				END SELECT
			ENDIF
		ENDIF
	ENDIF


	IF(ifl(1)+ifl(2).EQ.0.AND.(MOD(ABS(ifl(1)),4).EQ.0.OR.MOD(ABS(ifl(1)),4).EQ.3))THEN
		IF(onlyqcd)THEN
			SELECT CASE(iqnum)
			CASE(1)
				sf_ct=sf_ct*(u11_ct*ub21_ct+ub11_ct*u21_ct)/(u12_ct*ub22_ct+ub12_ct*u22_ct)
			CASE(2)
				sf_ct=sf_ct*(u11_ct*ub21_ct+ub11_ct*u21_ct+d11_ct*db21_ct+d21_ct*db11_ct)&
				/(u12_ct*ub22_ct+ub12_ct*u22_ct+d12_ct*db22_ct+d22_ct*db12_ct)
			CASE(3)
				sf_ct=sf_ct*(u11_ct*ub21_ct+ub11_ct*u21_ct+d11_ct*db21_ct+d21_ct*db11_ct&
				+s11_ct*sb21_ct+s21_ct*sb11_ct)&
				/(u12_ct*ub22_ct+ub12_ct*u22_ct+d12_ct*db22_ct+d22_ct*db12_ct+s12_ct*sb22_ct+s22_ct*sb12_ct)
			CASE(4)
				sf_ct=sf_ct*(u11_ct*ub21_ct+ub11_ct*u21_ct+d11_ct*db21_ct+d21_ct*db11_ct&
				+s11_ct*sb21_ct+s21_ct*sb11_ct+c11_ct*cb21_ct+cb11_ct*c21_ct)&
				/(u12_ct*ub22_ct+ub12_ct*u22_ct+d12_ct*db22_ct+d22_ct*db12_ct&
				+s12_ct*sb22_ct+s22_ct*sb12_ct+c12_ct*cb22_ct+cb12_ct*c22_ct)
			CASE DEFAULT
				sf_ct=sf_ct*(u11_ct*ub21_ct+ub11_ct*u21_ct+d11_ct*db21_ct+d21_ct*db11_ct&
				+s11_ct*sb21_ct+s21_ct*sb11_ct+c11_ct*cb21_ct+cb11_ct*c21_ct+b11_ct*bb21_ct+bb11_ct*b21_ct)&
				/(u12_ct*ub22_ct+ub12_ct*u22_ct+d12_ct*db22_ct+d22_ct*db12_ct&
				+s12_ct*sb22_ct+s22_ct*sb12_ct+c12_ct*cb22_ct+cb12_ct*c22_ct+b12_ct*bb22_ct+bb12_ct*b22_ct)
			END SELECT
		ELSE
			IF(MOD(ABS(ifl(1)),4).EQ.3)THEN
				SELECT CASE(iqnum)
				CASE(1,2,3)
					sf_ct=sf_ct*(u11_ct*ub21_ct+ub11_ct*u21_ct)/(u12_ct*ub22_ct+ub12_ct*u22_ct)
				CASE DEFAULT
					sf_ct=sf_ct*(u11_ct*ub21_ct+ub11_ct*u21_ct+c11_ct*cb21_ct+c21_ct*cb11_ct)&
					/(u12_ct*ub22_ct+ub12_ct*u22_ct+c12_ct*cb22_ct+c22_ct*cb12_ct)
				END SELECT
			ELSEIF(MOD(ABS(ifl(1)),4).EQ.0)THEN
				SELECT CASE(iqnum)
				CASE(1)
					sf_ct=sf_ct
				CASE(2)
					sf_ct=sf_ct*(d11_ct*db21_ct+d21_ct*db11_ct)/(d12_ct*db22_ct+d22_ct*db12_ct)
				CASE(3,4)
					sf_ct=sf_ct*(s11_ct*sb21_ct+sb11_ct*s21_ct+d11_ct*db21_ct+d21_ct*db11_ct)&
					/(s12_ct*sb22_ct+sb12_ct*s22_ct+d12_ct*db22_ct+d22_ct*db12_ct)
				CASE DEFAULT
					sf_ct=sf_ct*(s11_ct*sb21_ct+sb11_ct*s21_ct+d11_ct*db21_ct+d21_ct*db11_ct&
					+b11_ct*bb21_ct+b21_ct*bb11_ct)/(s12_ct*sb22_ct+sb12_ct*s22_ct+&
					d12_ct*db22_ct+d22_ct*db12_ct+b12_ct*bb22_ct+b22_ct*bb12_ct)
				END SELECT
			ENDIF
		ENDIF
	ELSEIF(ifl(1).NE.35.AND.ifl(2).NE.35)THEN
		PRINT *,"Please choose quarksumQ=.FALSE."
		STOP
	ENDIF

ELSE      
	IF(ipp.EQ.1)THEN
	! ppbar
		SELECT CASE(ifl(1))  
		CASE(3)
			sf_ct=sf_ct*u11_ct/u12_ct
		CASE(-3) 
			sf_ct=sf_ct*ub11_ct/ub12_ct
		CASE(4)  
			sf_ct=sf_ct*d11_ct/d12_ct
		CASE(-4) 
			sf_ct=sf_ct*db11_ct/db12_ct
		CASE(7)  
			sf_ct=sf_ct*c11_ct/c12_ct
		CASE(-7) 
			sf_ct=sf_ct*cb11_ct/cb12_ct
		CASE(8)  
			sf_ct=sf_ct*s11_ct/s12_ct
		CASE(-8) 
			sf_ct=sf_ct*sb11_ct/sb12_ct
!   IF(ifl(1).EQ.11) sf_ct=sf_ct*t1_ct
!   IF(ifl(1).EQ.-11)sf_ct=sf_ct*tb1_ct
		CASE(12) 
			sf_ct=sf_ct*b11_ct/b12_ct
		CASE(-12)
			sf_ct=sf_ct*bb11_ct/bb12_ct
		CASE(35) 
			sf_ct=sf_ct*glu11_ct/glu12_ct
		END SELECT

		SELECT CASE(ifl(2))
		CASE(3)
			sf_ct=sf_ct*ub21_ct/ub22_ct
		CASE(-3) 
			sf_ct=sf_ct*u21_ct/u22_ct
		CASE(4)  
			sf_ct=sf_ct*db21_ct/db22_ct
		CASE(-4) 
			sf_ct=sf_ct*d21_ct/d22_ct
		CASE(7)  
			sf_ct=sf_ct*cb21_ct/cb22_ct
		CASE(-7) 
			sf_ct=sf_ct*c21_ct/c22_ct
		CASE(8)  
			sf_ct=sf_ct*sb21_ct/sb22_ct
		CASE(-8) 
			sf_ct=sf_ct*s21_ct/s22_ct
!       if(ifl(2).eq.11) sf=sf*t2
!       if(ifl(2).eq.-11)sf=sf*tb2
		CASE(12) 
			sf_ct=sf_ct*bb21_ct/bb22_ct
		CASE(-12)
			sf_ct=sf_ct*b21_ct/b22_ct
		CASE(35) 
			sf_ct=sf_ct*glu21_ct/glu22_ct
		END SELECT
	ELSEIF(ipp.EQ.2)THEN
! pp
		SELECT CASE(ifl(1))
		CASE(3)  
			sf_ct=sf_ct*u11_ct/u12_ct
		CASE(-3) 
			sf_ct=sf_ct*ub11_ct/ub12_ct
		CASE(4)  
			sf_ct=sf_ct*d11_ct/d12_ct
		CASE(-4) 
			sf_ct=sf_ct*db11_ct/db12_ct
		CASE(7)  
			sf_ct=sf_ct*c11_ct/c12_ct
		CASE(-7) 
			sf_ct=sf_ct*cb11_ct/cb12_ct
		CASE(8)  
			sf_ct=sf_ct*s11_ct/s12_ct
		CASE(-8) 
			sf_ct=sf_ct*sb11_ct/sb12_ct
!       if(ifl(1).eq.11) sf=sf*t1
!       if(ifl(1).eq.-11)sf=sf*tb1
		CASE(12) 
			sf_ct=sf_ct*b11_ct/b12_ct
		CASE(-12)
			sf_ct=sf_ct*bb11_ct/bb12_ct
		CASE(35) 
			sf_ct=sf_ct*glu11_ct/glu12_ct
		END SELECT

		SELECT CASE(ifl(2))
		CASE(3)  
			sf_ct=sf_ct*u21_ct/u22_ct
		CASE(-3) 
			sf_ct=sf_ct*ub21_ct/ub22_ct
		CASE(4)  
			sf_ct=sf_ct*d21_ct/d22_ct
		CASE(-4) 
			sf_ct=sf_ct*db21_ct/db22_ct
		CASE(7)  
			sf_ct=sf_ct*c21_ct/c22_ct
		CASE(-7) 
			sf_ct=sf_ct*cb21_ct/cb22_ct
		CASE(8)  
			sf_ct=sf_ct*s21_ct/s22_ct
		CASE(-8) 
			sf_ct=sf_ct*sb21_ct/sb22_ct
!       if(ifl(2).eq.11) sf=sf*t2
!       if(ifl(2).eq.-11)sf=sf*tb2
		CASE(12) 
			sf_ct=sf_ct*b21_ct/b22_ct
		CASE(-12)
			sf_ct=sf_ct*bb21_ct/bb22_ct
		CASE(35) 
			sf_ct=sf_ct*glu21_ct/glu22_ct
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

w=w*sf_ct
IF(init.EQ.0)init=1

END SUBROUTINE PDFreweight
END MODULE KT_Clustering
