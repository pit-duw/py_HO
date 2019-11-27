MODULE Helac_pan1
USE Helac_Global
USE Helac_Func_1
USE Helac_Feynman
USE Helac_wavef
USE Helac_SM_FeynRule
USE Helac_pan2
USE Helac_ranmar_mod
!USE Projectors
IMPLICIT NONE
INTEGER,DIMENSION(:),ALLOCATABLE::jxx,kcou     
INTEGER,DIMENSION(:,:),ALLOCATABLE::ih,ih1,icole,list1,ico,ife
INTEGER,DIMENSION(:,:,:),ALLOCATABLE::ipp,list,list_tmp
INTEGER::iv1,iv2
LOGICAL::log1,log2,log3
SAVE ih,ih1,icole,ipp,jxx,kcou,list,list1
SAVE ife,iv1,iv2
INTEGER,DIMENSION(-12:44)::an=(/12,11,10,9,8,7,6,5,4,3,2,1,0,&
                      -1,-2,-3,-4,-5,-6,-7,-8,-9,-10,-11,-12,&
                      13,14,15,16,17,18,19,20,&
                      21,22,23,24,25,26,27,28,29,30,&
                      31,32,34,33,35,36,37,38,39,40,41,42,44,43/)
SAVE an
COMPLEX(KIND=DBL)::zero,zi
SAVE zero,zi
INTEGER::kcount,icc,kmax=0,imax=0,igluon        
SAVE kcount,icc,kmax,imax,igluon
COMPLEX(KIND=DBL),DIMENSION(:,:,:),ALLOCATABLE::zy
SAVE zy
COMPLEX(KIND=DBL),DIMENSION(:,:),ALLOCATABLE::zyq
SAVE zyq
INTEGER::init_pi=0
SAVE init_pi
REAL(KIND=DBL)::pi
SAVE pi
COMPLEX(KIND=DBL),DIMENSION(4)::z1,z2,z3,zp0,zp1,zp2
COMPLEX(KIND=DBL),DIMENSION(0:7,4)::z0,z11,z21,z31
COMPLEX(KIND=DBL),DIMENSION(5)::zq1,zPQQ
INTEGER::m3,m4 
SAVE m3,m4

CONTAINS

SUBROUTINE Helac_pan(nc1)
IMPLICIT NONE
INTEGER,INTENT(IN)::nc1
INTEGER::istat
CALL Helac_mypi(pi)
zero=DCMPLX(dnou(0),dnou(0))
zi=DCMPLX(dnou(0),dnou(1))
ncc=nc1
!WRITE(nunit1,*)'ncc,ngues',ncc,ngues
IF(ALLOCATED(ih))THEN
	DEALLOCATE(ih)
ENDIF
IF(ALLOCATED(ih1))THEN
	DEALLOCATE(ih1)
ENDIF
IF(ALLOCATED(icole))THEN
	DEALLOCATE(icole)
ENDIF
IF(ALLOCATED(ipp))THEN
	DEALLOCATE(ipp)
ENDIF
IF(ALLOCATED(kcou))THEN
	DEALLOCATE(kcou)
ENDIF
IF(ALLOCATED(list))THEN
	DEALLOCATE(list)
ENDIF
IF(ALLOCATED(jxx))THEN
	DEALLOCATE(jxx)
ENDIF
IF(ALLOCATED(list1))THEN
	DEALLOCATE(list1)
ENDIF
IF(ALLOCATED(ife))THEN
	DEALLOCATE(ife)
ENDIF
ALLOCATE( ih(2**n-2,-12:44),ih1(2**n-2,-12:44),&
          icole(2**n-2,2),ipp(0:2**n-2,-12:44,3),&
		  kcou(ncc),list(0:ngues,18,ncc),jxx(ngues),&
          list1(n,3),ife(2**n-2,-12:44),stat=istat)
IF(istat.NE.0)THEN
   WRITE(*,*)'warning: allocation is not working properly in Helac_pan'
   STOP
ENDIF
icc=1
iv1=31
iv2=34
ife(1:2**n-2,-12:44)=0
list(0:ngues,1:18,1:ncc)=0
m3=0
m4=0
END SUBROUTINE Helac_pan
       
! ----------------------------------------------------------------------- 
      
SUBROUTINE Helac_pan_ini()
IMPLICIT NONE
INTEGER::k,ii,ixix       
IF(ncc.GT.1)iv2=35
ih(1:2**n-2,-12:44)=0
ih1(1:2**n-2,-12:44)=0
ife(1:2**n-3,-12:44)=0
icole(1:2**n-2,1:2)=0
DO k=1,n
   ii=2**(k-1)
   icole(ii,1)=icol(k,1)
   icole(ii,2)=icol(k,2)
   ! Incoming
   IF(io(k).EQ.1)ih(ii,ifl(k))=1
   IF(io(k).EQ.1)ife(ii,ifl(k))=1
   ! outgoing, particle<->antiparticle
   IF(io(k).EQ.-1)ih(ii,an(ifl(k)))=1
   IF(io(k).EQ.-1)ife(ii,an(ifl(k)))=1
   ! traceless of T^a
   IF(ifl(k).EQ.35)THEN
      IF(icol(k,1).EQ.icol(k,2))THEN
          icole(ii,1:2)=0
      ENDIF
   ENDIF
ENDDO
! list1(k,1)=2**(k-1)
! list1(k,2)=the flavor number of incoming particles 
! or the anti flavor number of outgoing particles
! list1(k,3)=1,A,g;2,Z,W+,W-;3,incoming fermion;4,outgoing fermion;
! 5,incoming antifermion;6,outgoing antifermion
DO k=1,n
   list1(k,1)=2**(k-1)
   IF(io(k).EQ.1)list1(k,2)=ifl(k)
   IF(io(k).EQ.-1)list1(k,2)=an(ifl(k))
   ! massless vector bosons
   IF(list1(k,2).EQ.31.OR.list1(k,2).EQ.35)list1(k,3)=1
   ! massive vector bosons
   IF(list1(k,2).GE.32.AND.list1(k,2).LE.34)list1(k,3)=2
   IF(list1(k,2).GE.1.AND.list1(k,2).LE.12)THEN
     IF(io(k).EQ.1)list1(k,3)=3
     IF(io(k).EQ.-1)list1(k,3)=6
   ENDIF
   IF(list1(k,2).GE.-12.AND.list1(k,2).LE.-1)THEN
     IF(io(k).EQ.1)list1(k,3)=5
     IF(io(k).EQ.-1)list1(k,3)=4
   ENDIF
   IF(list1(k,2).GE.41.AND.list1(k,2).LE.44)list1(k,3)=7
ENDDO
kcount=0
jxx(1:ngues)=0
END SUBROUTINE Helac_pan_ini

! ----------------------------------------------------------------------- 

SUBROUTINE Helac_v3(i0,le,m0,i1,i2,i3)
IMPLICIT NONE
INTEGER,INTENT(IN)::i1,i2,i3,m0,i0,le
INTEGER::lpol,m1,m2,irou,i23,i_ex
m3=0
lpol=0
DO m1=iv1,iv2
   DO m2=iv1,iv2
      IF(m1.EQ.m2)THEN
         IF(Helac_level(n,i1).GT.Helac_level(n,i2))CYCLE
         IF(Helac_level(n,i1).EQ.Helac_level(n,i2).AND.i1.GT.i2)CYCLE
      ENDIF
	  IF(ih(i1,m1).NE.0)THEN
         IF(ih(i2,m2).NE.0)THEN
           IF(zgv3(m0,m1,m2).NE.zero)THEN
              igluon=0              ! once more
			  ! label 1
	          IF(m1.EQ.35.AND.icc.LE.ncc)THEN
			      IF(icole(i1,1)*icole(i1,2)*icole(i2,1)*icole(i2,2).EQ.0)CYCLE
	              IF(icole(i1,1).EQ.icole(i2,2))THEN
	                   icole(i0,2)=icole(i1,2)
	                   icole(i0,1)=icole(i2,1)
	                   igluon=1
	              ELSEIF(icole(i1,2).EQ.icole(i2,1))THEN
	                   icole(i0,2)=icole(i2,2)
	                   icole(i0,1)=icole(i1,1)
	                   igluon=2
	              ENDIF
                  IF(igluon.EQ.0)CYCLE 
				  ! there is no U(1) gluon between 3 gluon vertex
	              IF(icole(i0,1).EQ.icole(i0,2))THEN
                       icole(i0,1:2)=0
	                   igluon=0
	                   CYCLE
	              ENDIF
	           ENDIF
               ih(i0,m0)=1
               irou=4
 !          include 'list.h' 
               kcount=kcount+1
               list(kcount,1,icc)=irou
               list(kcount,2,icc)=i0
               list(kcount,3,icc)=m0
               list(kcount,7,icc)=i1
               list(kcount,8,icc)=m1
               list(kcount,10,icc)=i2
               list(kcount,11,icc)=m2
               list(kcount,13,icc)=i3
               list(kcount,14,icc)=m3
               list(kcount,16,icc)=lpol
               IF(i3.EQ.0)THEN
                   list(kcount,17,icc)=Helac_ipsgn(i1,i2)
               ELSE
                   i23=i2+i3
                   list(kcount,17,icc)=Helac_ipsgn(i2,i3)*Helac_ipsgn(i1,i23)
               ENDIF
               list(kcount,18,icc)=igluon
               IF(MOD(kcount,100).EQ.0)WRITE(nunit1,*)&
			   'I am calculating the ',kcount,'-th subamplitude'
               i_ex=0
			   ! If we haven't dealt with this situation exactly, we set i_ex=1.
	           IF(kcount.GT.1)THEN
	               IF(list(kcount-1,2,icc).NE.i0)i_ex=1
	               IF(list(kcount-1,3,icc).NE.m0)i_ex=1
	               IF(list(kcount-1,7,icc).NE.i1)i_ex=1
   	               IF(list(kcount-1,8,icc).NE.m1)i_ex=1
                   IF(list(kcount-1,10,icc).NE.i2)i_ex=1
                   IF(list(kcount-1,11,icc).NE.m2)i_ex=1
                   IF(list(kcount-1,13,icc).NE.i3)i_ex=1
                   IF(list(kcount-1,14,icc).NE.m3)i_ex=1
	           ELSEIF(kcount.EQ.1)THEN
                   i_ex=1
	           ENDIF
               IF(igluon.GE.2)i_ex=0
               IF(i_ex.EQ.1)THEN
                   IF(i3.EQ.0)ife(i0,m0)=ife(i0,m0)+ife(i1,m1)*ife(i2,m2)
                   IF(i3.NE.0)ife(i0,m0)=ife(i0,m0)+ife(i1,m1)*ife(i2,m2)*ife(i3,m3)
               ENDIF
! end of list.h 
           ENDIF
         ENDIF
       ENDIF
	   igluon=0              ! once more
    ENDDO
ENDDO
END SUBROUTINE Helac_v3

! ----------------------------------------------------------------------- 
SUBROUTINE Helac_v4(i0,le,m0,i1,i2,i3)
IMPLICIT NONE
INTEGER,INTENT(IN)::i1,i2,i3,m0,i0,le
INTEGER::lpol,m1,m2,irou,i23,i_ex
lpol=0
DO m1=iv1,iv2
   DO m2=iv1,iv2
      DO m3=iv1,iv2 
        IF(m1.EQ.m2)THEN
           IF(Helac_level(n,i1).GT.Helac_level(n,i2))CYCLE
           IF(Helac_level(n,i1).EQ.Helac_level(n,i2).AND.i1.GT.i2)CYCLE
        ENDIF
        IF(m1.EQ.m3)THEN
           IF(Helac_level(n,i1).GT.Helac_level(n,i3))CYCLE
           IF(Helac_level(n,i1).EQ.Helac_level(n,i3).AND.i1.GT.i3)CYCLE
        ENDIF
        IF(m2.EQ.m3)THEN
           IF(Helac_level(n,i2).GT.Helac_level(n,i3))CYCLE
           IF(Helac_level(n,i2).EQ.Helac_level(n,i3).AND.i2.GT.i3)CYCLE
        ENDIF  
        lll1:IF(ih(i1,m1).NE.0)THEN
           lll2:IF(ih(i2,m2).NE.0)THEN
              lll3:IF(ih(i3,m3).NE.0)THEN
!			     PRINT *,zero
                 lll4:IF(zgv4(m0,m1,m2,m3).NE.zero)THEN
                    igluon=0              ! once more
					! Label 2
                    IF(m1.EQ.35.AND.icc.LE.ncc)THEN
                         IF(icole(i1,1)*icole(i1,2)*icole(i2,1)*icole(i2,2)&
                         *icole(i3,1)*icole(i3,2).EQ.0)CYCLE
                         IF(icole(i1,1).EQ.icole(i2,2).AND.icole(i1,2).EQ.icole(i3,1))THEN
                              icole(i0,2)=icole(i3,2)
                              icole(i0,1)=icole(i2,1)
                              igluon=1
                         ELSEIF(icole(i1,1).EQ.icole(i3,2).AND.icole(i1,2).EQ.icole(i2,1))THEN
                              icole(i0,2)=icole(i2,2)
                              icole(i0,1)=icole(i3,1)
                              igluon=2+2
                         ELSEIF(icole(i2,1).EQ.icole(i1,2).AND.icole(i2,2).EQ.icole(i3,1))THEN
                              icole(i0,2)=icole(i3,2)
                              icole(i0,1)=icole(i1,1)
                              igluon=3+2
                         ELSEIF(icole(i2,1).EQ.icole(i3,2).AND.icole(i2,2).EQ.icole(i1,1))THEN
                              icole(i0,2)=icole(i1,2)
                              icole(i0,1)=icole(i3,1)
                              igluon=3+2
                         ELSEIF(icole(i3,1).EQ.icole(i1,2).AND.icole(i3,2).EQ.icole(i2,1))THEN
                              icole(i0,2)=icole(i2,2)
                              icole(i0,1)=icole(i1,1)
                              igluon=4+2
                         ELSEIF(icole(i3,1).EQ.icole(i2,2).AND.icole(i3,2).EQ.icole(i1,1))THEN
                              icole(i0,2)=icole(i1,2)
                              icole(i0,1)=icole(i2,1)
                              igluon=4+2
                         ELSE
                              CYCLE
                         ENDIF
	                     IF(igluon.EQ.0)CYCLE
						 ! there is no U(1) gluon between 4 gluon vertex
	                     IF(icole(i0,1).EQ.icole(i0,2))THEN
                              icole(i0,1:2)=0
	                          igluon=0
	                          CYCLE
	                     ENDIF
                     ENDIF

                     ih(i0,m0)=1
                     irou=5
!           include 'list.h'         
                     kcount=kcount+1
                     list(kcount,1,icc)=irou
                     list(kcount,2,icc)=i0
                     list(kcount,3,icc)=m0
                     list(kcount,7,icc)=i1
                     list(kcount,8,icc)=m1
                     list(kcount,10,icc)=i2
                     list(kcount,11,icc)=m2
                     list(kcount,13,icc)=i3
                     list(kcount,14,icc)=m3
                     list(kcount,16,icc)=lpol
                     IF(i3.EQ.0)THEN
                        list(kcount,17,icc)=Helac_ipsgn(i1,i2)
                     ELSE
                        i23=i2+i3
                        list(kcount,17,icc)=Helac_ipsgn(i2,i3)*Helac_ipsgn(i1,i23)
                     ENDIF
                     list(kcount,18,icc)=igluon
                     IF(MOD(kcount,100).EQ.0)WRITE(nunit1,*)&
			         'I am calculating the ',kcount,'-th subamplitude'
                     i_ex=0
			   ! If we haven't dealt with this situation exactly, we set i_ex=1.
	                 IF(kcount.GT.1)THEN
	                   IF(list(kcount-1,2,icc).NE.i0)i_ex=1
	                   IF(list(kcount-1,3,icc).NE.m0)i_ex=1
	                   IF(list(kcount-1,7,icc).NE.i1)i_ex=1
   	                   IF(list(kcount-1,8,icc).NE.m1)i_ex=1
                       IF(list(kcount-1,10,icc).NE.i2)i_ex=1
                       IF(list(kcount-1,11,icc).NE.m2)i_ex=1
                       IF(list(kcount-1,13,icc).NE.i3)i_ex=1
                       IF(list(kcount-1,14,icc).NE.m3)i_ex=1
	                 ELSEIF(kcount.EQ.1)THEN
                       i_ex=1
	                 ENDIF
                     IF(igluon.GE.2)i_ex=0
                     IF(i_ex.EQ.1)THEN
                        IF(i3.EQ.0)ife(i0,m0)=ife(i0,m0)+ife(i1,m1)*ife(i2,m2)
                        IF(i3.NE.0)ife(i0,m0)=ife(i0,m0)+ife(i1,m1)*ife(i2,m2)*ife(i3,m3)
                     ENDIF
! end of list.h           
	                 igluon=0              ! once more
                 ENDIF  lll4
              ENDIF lll3
           ENDIF lll2
        ENDIF lll1
         
      ENDDO
   ENDDO
ENDDO
END SUBROUTINE Helac_v4

! ----------------------------------------------------------------------- 

SUBROUTINE Helac_vvs(i0,le,m0,i1,i2,i3)
IMPLICIT NONE
INTEGER,INTENT(IN)::i1,i2,i3,m0,i0,le
INTEGER::lpol,m1,m2,irou,i23,i_ex
m3=0
lpol=0
         
DO m1=iv1,iv2
   DO m2=41,44
      IF(zgvvs(m0,m1,m2).NE.zero.AND.ih(i1,m1).NE.0.AND.ih(i2,m2).NE.0)THEN
         ih(i0,m0)=1
         irou=9
         !           include 'list.h'         
         kcount=kcount+1
         list(kcount,1,icc)=irou
         list(kcount,2,icc)=i0
         list(kcount,3,icc)=m0
         list(kcount,7,icc)=i1
         list(kcount,8,icc)=m1
         list(kcount,10,icc)=i2
         list(kcount,11,icc)=m2
         list(kcount,13,icc)=i3
         list(kcount,14,icc)=m3
         list(kcount,16,icc)=lpol
         IF(i3.EQ.0)THEN
             list(kcount,17,icc)=Helac_ipsgn(i1,i2)
         ELSE
             i23=i2+i3
             list(kcount,17,icc)=Helac_ipsgn(i2,i3)*Helac_ipsgn(i1,i23)
         ENDIF
         list(kcount,18,icc)=igluon
         IF(MOD(kcount,100).EQ.0)WRITE(nunit1,*)&
		'I am calculating the ',kcount,'-th subamplitude'
         i_ex=0
			   ! If we haven't dealt with this situation exactly, we set i_ex=1.
	     IF(kcount.GT.1)THEN
	         IF(list(kcount-1,2,icc).NE.i0)i_ex=1
	         IF(list(kcount-1,3,icc).NE.m0)i_ex=1
	         IF(list(kcount-1,7,icc).NE.i1)i_ex=1
   	         IF(list(kcount-1,8,icc).NE.m1)i_ex=1
             IF(list(kcount-1,10,icc).NE.i2)i_ex=1
             IF(list(kcount-1,11,icc).NE.m2)i_ex=1
             IF(list(kcount-1,13,icc).NE.i3)i_ex=1
             IF(list(kcount-1,14,icc).NE.m3)i_ex=1
	     ELSEIF(kcount.EQ.1)THEN
             i_ex=1
	     ENDIF
         IF(igluon.GE.2)i_ex=0
         IF(i_ex.EQ.1)THEN
             IF(i3.EQ.0)ife(i0,m0)=ife(i0,m0)+ife(i1,m1)*ife(i2,m2)
             IF(i3.NE.0)ife(i0,m0)=ife(i0,m0)+ife(i1,m1)*ife(i2,m2)*ife(i3,m3)
         ENDIF
! end of list.h                  
      ENDIF
   ENDDO
ENDDO
END SUBROUTINE Helac_vvs

! ----------------------------------------------------------------------- 

SUBROUTINE Helac_vss(i0,le,m0,i1,i2,i3)
IMPLICIT NONE
INTEGER,INTENT(IN)::i1,i2,i3,m0,i0,le
INTEGER::lpol,m1,m2,irou,i23,i_ex,igo
m3=0
lpol=0
         
DO m1=41,44
   DO m2=41,44
       igo=1
       IF(m1.EQ.m2.AND.i1.GT.i2)igo=-1   !12 
       IF(zgvss(m0,m1,m2).NE.zero.AND.ih(i1,m1).NE.0&
         .AND.ih(i2,m2).NE.0.AND.igo.EQ.1)THEN
         ih(i0,m0)=1
         irou=16
         !           include 'list.h'         
         kcount=kcount+1
         list(kcount,1,icc)=irou
         list(kcount,2,icc)=i0
         list(kcount,3,icc)=m0
         list(kcount,7,icc)=i1
         list(kcount,8,icc)=m1
         list(kcount,10,icc)=i2
         list(kcount,11,icc)=m2
         list(kcount,13,icc)=i3
         list(kcount,14,icc)=m3
         list(kcount,16,icc)=lpol
         IF(i3.EQ.0)THEN
             list(kcount,17,icc)=Helac_ipsgn(i1,i2)
         ELSE
             i23=i2+i3
             list(kcount,17,icc)=Helac_ipsgn(i2,i3)*Helac_ipsgn(i1,i23)
         ENDIF
         list(kcount,18,icc)=igluon
         IF(MOD(kcount,100).EQ.0)WRITE(nunit1,*)&
		'I am calculating the ',kcount,'-th subamplitude'
         i_ex=0
			   ! If we haven't dealt with this situation exactly, we set i_ex=1.
	     IF(kcount.GT.1)THEN
	        IF(list(kcount-1,2,icc).NE.i0)i_ex=1
	        IF(list(kcount-1,3,icc).NE.m0)i_ex=1
	        IF(list(kcount-1,7,icc).NE.i1)i_ex=1
   	        IF(list(kcount-1,8,icc).NE.m1)i_ex=1
            IF(list(kcount-1,10,icc).NE.i2)i_ex=1
            IF(list(kcount-1,11,icc).NE.m2)i_ex=1
            IF(list(kcount-1,13,icc).NE.i3)i_ex=1
            IF(list(kcount-1,14,icc).NE.m3)i_ex=1
	      ELSEIF(kcount.EQ.1)THEN
            i_ex=1
	      ENDIF
          IF(igluon.GE.2)i_ex=0
          IF(i_ex.EQ.1)THEN
             IF(i3.EQ.0)ife(i0,m0)=ife(i0,m0)+ife(i1,m1)*ife(i2,m2)
             IF(i3.NE.0)ife(i0,m0)=ife(i0,m0)+ife(i1,m1)*ife(i2,m2)*ife(i3,m3)
          ENDIF
! end of list.h        
        ENDIF
    ENDDO
ENDDO
END SUBROUTINE Helac_vss

! ----------------------------------------------------------------------- 

SUBROUTINE Helac_vff(i0,le,m0,i1,i2,i3)
IMPLICIT NONE
INTEGER,INTENT(IN)::i1,i2,i3,m0,i0,le
INTEGER::lpol,m1,m2,irou,i23,i_ex
m3=0
lpol=0
! to remove any contribution with a g,Z and W as propagators
!        include 'removeZ.h'
!        include 'removeW.h'
!        include 'removeA.h'
!        include 'removeZW.h'

DO m1=-12,-1
    DO m2=1,12
![ special instructions
!        include 'specialinstvff.h'
!] special instructions
         IF(ih(i1,m1).NE.0.AND.ih(i2,m2).NE.0.AND.&
             (zgvffr(m0,m1,m2).NE.zero.OR.zgvffl(m0,m1,m2).NE.zero))THEN
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   for -12<m<0    icole=(0,x)
!   for   0<m<12   icole=(x,0)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
             igluon=0          ! once more
			 ! label 3
             IF(m0.EQ.35)THEN
			   IF(ColorSinglet2(icole(i2,1),icole(i1,2)))CYCLE ! Quarkonium
               IF(icole(i1,2).EQ.icole(i2,1))THEN
				 igluon=2
				 IF(Helac_level(n,i0).EQ.n-1)igluon=3
				 icole(i0,1:2)=0
               ELSE
                 igluon=1
                 icole(i0,1)=icole(i2,1)
                 icole(i0,2)=icole(i1,2)
               ENDIF
             ELSE
               IF(icole(i1,2).NE.icole(i2,1))CYCLE
               icole(i0,1:2)=0
             ENDIF
 !           include 'largeNclimit.h'
             ih(i0,m0)=1
             irou=1
         !           include 'list.h'         
             kcount=kcount+1
             list(kcount,1,icc)=irou
             list(kcount,2,icc)=i0
             list(kcount,3,icc)=m0
             list(kcount,7,icc)=i1
             list(kcount,8,icc)=m1
             list(kcount,10,icc)=i2
             list(kcount,11,icc)=m2
             list(kcount,13,icc)=i3
             list(kcount,14,icc)=m3
             list(kcount,16,icc)=lpol
             IF(i3.EQ.0)THEN
                list(kcount,17,icc)=Helac_ipsgn(i1,i2)
             ELSE
                i23=i2+i3
                list(kcount,17,icc)=Helac_ipsgn(i2,i3)*Helac_ipsgn(i1,i23)
             ENDIF
             list(kcount,18,icc)=igluon
             IF(MOD(kcount,100).EQ.0)WRITE(nunit1,*)&
		    'I am calculating the ',kcount,'-th subamplitude'
             i_ex=0
			   ! If we haven't dealt with this situation exactly, we set i_ex=1.
	         IF(kcount.GT.1)THEN
	            IF(list(kcount-1,2,icc).NE.i0)i_ex=1
	            IF(list(kcount-1,3,icc).NE.m0)i_ex=1
	            IF(list(kcount-1,7,icc).NE.i1)i_ex=1
   	            IF(list(kcount-1,8,icc).NE.m1)i_ex=1
                IF(list(kcount-1,10,icc).NE.i2)i_ex=1
                IF(list(kcount-1,11,icc).NE.m2)i_ex=1
                IF(list(kcount-1,13,icc).NE.i3)i_ex=1
                IF(list(kcount-1,14,icc).NE.m3)i_ex=1
	         ELSEIF(kcount.EQ.1)THEN
                i_ex=1
	         ENDIF
             IF(igluon.GE.2)i_ex=0
             IF(i_ex.EQ.1)THEN
                IF(i3.EQ.0)ife(i0,m0)=ife(i0,m0)+ife(i1,m1)*ife(i2,m2)
                IF(i3.NE.0)ife(i0,m0)=ife(i0,m0)+ife(i1,m1)*ife(i2,m2)*ife(i3,m3)
             ENDIF
! end of list.h        
         ENDIF
         igluon=0          ! once more
    ENDDO
ENDDO
END SUBROUTINE Helac_vff

! ----------------------------------------------------------------------- 

SUBROUTINE Helac_ffv(i0,le,m0,i1,i2,i3) 
IMPLICIT NONE
INTEGER,INTENT(IN)::i1,i2,i3,m0,i0,le
INTEGER::lpol,m1,m2,irou,i23,i_ex
m3=0
DO m2=1,12
   DO m1=iv1,iv2
      IF(ih(i1,m1).EQ.0.OR.ih(i2,m2).EQ.0)CYCLE
      IF(zgvffr(an(m1),-m0,m2).EQ.zero.AND.zgvffl(an(m1),-m0,m2).EQ.zero)CYCLE
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   for -12<m<0    icole=(0,x)
!   for   0<m<12   icole=(x,0)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      igluon=0            ! once more
      IF(m1.EQ.35)THEN
	  ! label 4
         IF(icole(i1,2).EQ.icole(i2,1))THEN
            igluon=1
            icole(i0,1)=icole(i1,1)
            icole(i0,2)=0
         ELSEIF(icole(i1,1).EQ.0)THEN
            igluon=2
            IF(Helac_level(n,i1).EQ.1)igluon=3
            icole(i0,1:2)=icole(i2,1:2)
         ELSE
            CYCLE
         ENDIF
      ELSE
         icole(i0,1:2)=icole(i2,1:2)
      ENDIF
 !          include 'largeNclimit.h'

      IF(zgvffl(an(m1),-m0,m2).NE.zero)THEN
           ih(i0,m0)=1
           irou=3
           lpol=1
         !           include 'list.h'         
           kcount=kcount+1
           list(kcount,1,icc)=irou
           list(kcount,2,icc)=i0
           list(kcount,3,icc)=m0
           list(kcount,7,icc)=i1
           list(kcount,8,icc)=m1
           list(kcount,10,icc)=i2
           list(kcount,11,icc)=m2
           list(kcount,13,icc)=i3
           list(kcount,14,icc)=m3
           list(kcount,16,icc)=lpol
           IF(i3.EQ.0)THEN
               list(kcount,17,icc)=Helac_ipsgn(i1,i2)
           ELSE
               i23=i2+i3
               list(kcount,17,icc)=Helac_ipsgn(i2,i3)*Helac_ipsgn(i1,i23)
           ENDIF
           list(kcount,18,icc)=igluon
           IF(MOD(kcount,100).EQ.0)WRITE(nunit1,*)&
	       'I am calculating the ',kcount,'-th subamplitude'
           i_ex=0
			   ! If we haven't dealt with this situation exactly, we set i_ex=1.
	      IF(kcount.GT.1)THEN
	          IF(list(kcount-1,2,icc).NE.i0)i_ex=1
	          IF(list(kcount-1,3,icc).NE.m0)i_ex=1
	          IF(list(kcount-1,7,icc).NE.i1)i_ex=1
   	          IF(list(kcount-1,8,icc).NE.m1)i_ex=1
              IF(list(kcount-1,10,icc).NE.i2)i_ex=1
              IF(list(kcount-1,11,icc).NE.m2)i_ex=1
              IF(list(kcount-1,13,icc).NE.i3)i_ex=1
              IF(list(kcount-1,14,icc).NE.m3)i_ex=1
	      ELSEIF(kcount.EQ.1)THEN
              i_ex=1
	      ENDIF
          IF(igluon.GE.2)i_ex=0
          IF(i_ex.EQ.1)THEN
              IF(i3.EQ.0)ife(i0,m0)=ife(i0,m0)+ife(i1,m1)*ife(i2,m2)
              IF(i3.NE.0)ife(i0,m0)=ife(i0,m0)+ife(i1,m1)*ife(i2,m2)*ife(i3,m3)
          ENDIF
! end of list.h 
      ENDIF
      IF(zgvffr(an(m1),-m0,m2).NE.zero)THEN
           ih(i0,m0)=1
           irou=3
           lpol=2
         !           include 'list.h'         
           kcount=kcount+1
           list(kcount,1,icc)=irou
           list(kcount,2,icc)=i0
           list(kcount,3,icc)=m0
           list(kcount,7,icc)=i1
           list(kcount,8,icc)=m1
           list(kcount,10,icc)=i2
           list(kcount,11,icc)=m2
           list(kcount,13,icc)=i3
           list(kcount,14,icc)=m3
           list(kcount,16,icc)=lpol
           IF(i3.EQ.0)THEN
               list(kcount,17,icc)=Helac_ipsgn(i1,i2)
           ELSE
               i23=i2+i3
               list(kcount,17,icc)=Helac_ipsgn(i2,i3)*Helac_ipsgn(i1,i23)
           ENDIF
           list(kcount,18,icc)=igluon
           IF(MOD(kcount,100).EQ.0)WRITE(nunit1,*)&
	       'I am calculating the ',kcount,'-th subamplitude'
           i_ex=0
			   ! If we haven't dealt with this situation exactly, we set i_ex=1.
	      IF(kcount.GT.1)THEN
	          IF(list(kcount-1,2,icc).NE.i0)i_ex=1
	          IF(list(kcount-1,3,icc).NE.m0)i_ex=1
	          IF(list(kcount-1,7,icc).NE.i1)i_ex=1
   	          IF(list(kcount-1,8,icc).NE.m1)i_ex=1
              IF(list(kcount-1,10,icc).NE.i2)i_ex=1
              IF(list(kcount-1,11,icc).NE.m2)i_ex=1
              IF(list(kcount-1,13,icc).NE.i3)i_ex=1
              IF(list(kcount-1,14,icc).NE.m3)i_ex=1
	      ELSEIF(kcount.EQ.1)THEN
              i_ex=1
	      ENDIF
          IF(igluon.GE.2)i_ex=0
          IF(i_ex.EQ.1)THEN
              IF(i3.EQ.0)ife(i0,m0)=ife(i0,m0)+ife(i1,m1)*ife(i2,m2)
              IF(i3.NE.0)ife(i0,m0)=ife(i0,m0)+ife(i1,m1)*ife(i2,m2)*ife(i3,m3)
          ENDIF
! end of list.h 
	  ENDIF  
      igluon=0            ! once more
   ENDDO
ENDDO
END SUBROUTINE Helac_ffv

! ----------------------------------------------------------------------- 

SUBROUTINE Helac_ffs(i0,le,m0,i1,i2,i3) 
IMPLICIT NONE
INTEGER,INTENT(IN)::i1,i2,i3,m0,i0,le
INTEGER::lpol,m1,m2,irou,i23,i_ex
m3=0
DO m2=1,12
   DO m1=41,44
      IF(ih(i1,m1).EQ.0.OR.ih(i2,m2).EQ.0)CYCLE
      IF(zgsffr(an(m1),-m0,m2).EQ.zero.AND.zgsffl(an(m1),-m0,m2).EQ.zero)CYCLE
      icole(i0,1:2)=icole(i2,1:2)
      IF(zgsffl(an(m1),-m0,m2).NE.zero)THEN
          ih(i0,m0)=1
          irou=6
          lpol=1
	           !           include 'list.h'         
          kcount=kcount+1
          list(kcount,1,icc)=irou
          list(kcount,2,icc)=i0
          list(kcount,3,icc)=m0
          list(kcount,7,icc)=i1
          list(kcount,8,icc)=m1
          list(kcount,10,icc)=i2
          list(kcount,11,icc)=m2
          list(kcount,13,icc)=i3
          list(kcount,14,icc)=m3
          list(kcount,16,icc)=lpol
          IF(i3.EQ.0)THEN
              list(kcount,17,icc)=Helac_ipsgn(i1,i2)
          ELSE
              i23=i2+i3
              list(kcount,17,icc)=Helac_ipsgn(i2,i3)*Helac_ipsgn(i1,i23)
          ENDIF
          list(kcount,18,icc)=igluon
          IF(MOD(kcount,100).EQ.0)WRITE(nunit1,*)&
	      'I am calculating the ',kcount,'-th subamplitude'
          i_ex=0
			   ! If we haven't dealt with this situation exactly, we set i_ex=1.
	      IF(kcount.GT.1)THEN
	          IF(list(kcount-1,2,icc).NE.i0)i_ex=1
	          IF(list(kcount-1,3,icc).NE.m0)i_ex=1
	          IF(list(kcount-1,7,icc).NE.i1)i_ex=1
   	          IF(list(kcount-1,8,icc).NE.m1)i_ex=1
              IF(list(kcount-1,10,icc).NE.i2)i_ex=1
              IF(list(kcount-1,11,icc).NE.m2)i_ex=1
              IF(list(kcount-1,13,icc).NE.i3)i_ex=1
              IF(list(kcount-1,14,icc).NE.m3)i_ex=1
	      ELSEIF(kcount.EQ.1)THEN
              i_ex=1
	      ENDIF
          IF(igluon.GE.2)i_ex=0
          IF(i_ex.EQ.1)THEN
              IF(i3.EQ.0)ife(i0,m0)=ife(i0,m0)+ife(i1,m1)*ife(i2,m2)
              IF(i3.NE.0)ife(i0,m0)=ife(i0,m0)+ife(i1,m1)*ife(i2,m2)*ife(i3,m3)
          ENDIF
! end of list.h              
      ENDIF
      IF(zgsffr(an(m1),-m0,m2).NE.zero)THEN
          ih(i0,m0)=1
          irou=6
          lpol=2
	           !           include 'list.h'         
          kcount=kcount+1
          list(kcount,1,icc)=irou
          list(kcount,2,icc)=i0
          list(kcount,3,icc)=m0
          list(kcount,7,icc)=i1
          list(kcount,8,icc)=m1
          list(kcount,10,icc)=i2
          list(kcount,11,icc)=m2
          list(kcount,13,icc)=i3
          list(kcount,14,icc)=m3
          list(kcount,16,icc)=lpol
          IF(i3.EQ.0)THEN
              list(kcount,17,icc)=Helac_ipsgn(i1,i2)
          ELSE
              i23=i2+i3
              list(kcount,17,icc)=Helac_ipsgn(i2,i3)*Helac_ipsgn(i1,i23)
          ENDIF
          list(kcount,18,icc)=igluon
          IF(MOD(kcount,100).EQ.0)WRITE(nunit1,*)&
	      'I am calculating the ',kcount,'-th subamplitude'
          i_ex=0
			   ! If we haven't dealt with this situation exactly, we set i_ex=1.
	      IF(kcount.GT.1)THEN
	          IF(list(kcount-1,2,icc).NE.i0)i_ex=1
	          IF(list(kcount-1,3,icc).NE.m0)i_ex=1
	          IF(list(kcount-1,7,icc).NE.i1)i_ex=1
   	          IF(list(kcount-1,8,icc).NE.m1)i_ex=1
              IF(list(kcount-1,10,icc).NE.i2)i_ex=1
              IF(list(kcount-1,11,icc).NE.m2)i_ex=1
              IF(list(kcount-1,13,icc).NE.i3)i_ex=1
              IF(list(kcount-1,14,icc).NE.m3)i_ex=1
	      ELSEIF(kcount.EQ.1)THEN
              i_ex=1
	      ENDIF
          IF(igluon.GE.2)i_ex=0
          IF(i_ex.EQ.1)THEN
              IF(i3.EQ.0)ife(i0,m0)=ife(i0,m0)+ife(i1,m1)*ife(i2,m2)
              IF(i3.NE.0)ife(i0,m0)=ife(i0,m0)+ife(i1,m1)*ife(i2,m2)*ife(i3,m3)
          ENDIF
! end of list.h      
      ENDIF  
   ENDDO
ENDDO
END SUBROUTINE Helac_ffs
! -----------------------------------------------------------------------      

SUBROUTINE Helac_affv(i0,le,m0,i1,i2,i3) 
IMPLICIT NONE
INTEGER,INTENT(IN)::i1,i2,i3,m0,i0,le
INTEGER::lpol,m1,m2,irou,i23,i_ex
m3=0
DO m2=-12,-1
   DO m1=iv1,iv2
      IF(ih(i1,m1).EQ.0.OR.ih(i2,m2).EQ.0)CYCLE
      IF(zgvffr(an(m1),m2,-m0).EQ.zero.AND.zgvffl(an(m1),m2,-m0).EQ.zero)CYCLE
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   for -12<m<0    icole=(0,x)
!   for   0<m<12   icole=(x,0)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      igluon=0            ! once more
      IF(m1.EQ.35)THEN
	  ! label 5
          IF(icole(i1,1).EQ.icole(i2,2))THEN
             igluon=1
             icole(i0,2)=icole(i1,2)
             icole(i0,1)=0
          ELSEIF(icole(i1,1).EQ.0)THEN
             igluon=2
             IF(Helac_level(n,i1).EQ.1)igluon=3 
             icole(i0,1:2)=icole(i2,1:2)
          ELSE
             CYCLE
          ENDIF
       ELSE
          icole(i0,1:2)=icole(i2,1:2)
       ENDIF
!         include 'largeNclimit.h'
       IF(zgvffl(an(m1),m2,-m0).NE.zero)THEN
	      ih(i0,m0)=1
          irou=2
          lpol=1
	           !           include 'list.h'         
          kcount=kcount+1
          list(kcount,1,icc)=irou
          list(kcount,2,icc)=i0
          list(kcount,3,icc)=m0
          list(kcount,7,icc)=i1
          list(kcount,8,icc)=m1
          list(kcount,10,icc)=i2
          list(kcount,11,icc)=m2
          list(kcount,13,icc)=i3
          list(kcount,14,icc)=m3
          list(kcount,16,icc)=lpol
          IF(i3.EQ.0)THEN
            list(kcount,17,icc)=Helac_ipsgn(i1,i2)
          ELSE
            i23=i2+i3
            list(kcount,17,icc)=Helac_ipsgn(i2,i3)*Helac_ipsgn(i1,i23)
          ENDIF
          list(kcount,18,icc)=igluon
          IF(MOD(kcount,100).EQ.0)WRITE(nunit1,*)&
	      'I am calculating the ',kcount,'-th subamplitude'
          i_ex=0
			   ! If we haven't dealt with this situation exactly, we set i_ex=1.
	      IF(kcount.GT.1)THEN
	         IF(list(kcount-1,2,icc).NE.i0)i_ex=1
	         IF(list(kcount-1,3,icc).NE.m0)i_ex=1
	         IF(list(kcount-1,7,icc).NE.i1)i_ex=1
   	         IF(list(kcount-1,8,icc).NE.m1)i_ex=1
             IF(list(kcount-1,10,icc).NE.i2)i_ex=1
             IF(list(kcount-1,11,icc).NE.m2)i_ex=1
             IF(list(kcount-1,13,icc).NE.i3)i_ex=1
             IF(list(kcount-1,14,icc).NE.m3)i_ex=1
	      ELSEIF(kcount.EQ.1)THEN
             i_ex=1
	      ENDIF
          IF(igluon.GE.2)i_ex=0
          IF(i_ex.EQ.1)THEN
             IF(i3.EQ.0)ife(i0,m0)=ife(i0,m0)+ife(i1,m1)*ife(i2,m2)
             IF(i3.NE.0)ife(i0,m0)=ife(i0,m0)+ife(i1,m1)*ife(i2,m2)*ife(i3,m3)
          ENDIF
! end of list.h   
       ENDIF
       IF(zgvffr(an(m1),m2,-m0).NE.zero)THEN
          ih(i0,m0)=1
          irou=2
          lpol=2
	           !           include 'list.h'         
          kcount=kcount+1
          list(kcount,1,icc)=irou
          list(kcount,2,icc)=i0
          list(kcount,3,icc)=m0
          list(kcount,7,icc)=i1
          list(kcount,8,icc)=m1
          list(kcount,10,icc)=i2
          list(kcount,11,icc)=m2
          list(kcount,13,icc)=i3
          list(kcount,14,icc)=m3
          list(kcount,16,icc)=lpol
          IF(i3.EQ.0)THEN
            list(kcount,17,icc)=Helac_ipsgn(i1,i2)
          ELSE
            i23=i2+i3
            list(kcount,17,icc)=Helac_ipsgn(i2,i3)*Helac_ipsgn(i1,i23)
          ENDIF
          list(kcount,18,icc)=igluon
          IF(MOD(kcount,100).EQ.0)WRITE(nunit1,*)&
	      'I am calculating the ',kcount,'-th subamplitude'
          i_ex=0
			   ! If we haven't dealt with this situation exactly, we set i_ex=1.
	      IF(kcount.GT.1)THEN
	         IF(list(kcount-1,2,icc).NE.i0)i_ex=1
	         IF(list(kcount-1,3,icc).NE.m0)i_ex=1
	         IF(list(kcount-1,7,icc).NE.i1)i_ex=1
   	         IF(list(kcount-1,8,icc).NE.m1)i_ex=1
             IF(list(kcount-1,10,icc).NE.i2)i_ex=1
             IF(list(kcount-1,11,icc).NE.m2)i_ex=1
             IF(list(kcount-1,13,icc).NE.i3)i_ex=1
             IF(list(kcount-1,14,icc).NE.m3)i_ex=1
	      ELSEIF(kcount.EQ.1)THEN
             i_ex=1
	      ENDIF
          IF(igluon.GE.2)i_ex=0
          IF(i_ex.EQ.1)THEN
             IF(i3.EQ.0)ife(i0,m0)=ife(i0,m0)+ife(i1,m1)*ife(i2,m2)
             IF(i3.NE.0)ife(i0,m0)=ife(i0,m0)+ife(i1,m1)*ife(i2,m2)*ife(i3,m3)
          ENDIF
! end of list.h     
       ENDIF
       igluon=0            ! once more
    ENDDO
ENDDO
END SUBROUTINE Helac_affv

! -----------------------------------------------------------------------      

SUBROUTINE Helac_affs(i0,le,m0,i1,i2,i3) 
IMPLICIT NONE
INTEGER,INTENT(IN)::i1,i2,i3,m0,i0,le
INTEGER::lpol,m1,m2,irou,i23,i_ex
m3=0        
DO m2=-12,-1
   DO m1=41,44
      IF(ih(i1,m1).NE.0.AND.ih(i2,m2).NE.0.AND.zgsffl(an(m1),m2,-m0).NE.zero)THEN
          icole(i0,1:2)=icole(i2,1:2)
          ih(i0,m0)=1
          irou=7
          lpol=1
	           !           include 'list.h'         
          kcount=kcount+1
          list(kcount,1,icc)=irou
          list(kcount,2,icc)=i0
          list(kcount,3,icc)=m0
          list(kcount,7,icc)=i1
          list(kcount,8,icc)=m1
          list(kcount,10,icc)=i2
          list(kcount,11,icc)=m2
          list(kcount,13,icc)=i3
          list(kcount,14,icc)=m3
          list(kcount,16,icc)=lpol
          IF(i3.EQ.0)THEN
            list(kcount,17,icc)=Helac_ipsgn(i1,i2)
          ELSE
            i23=i2+i3
            list(kcount,17,icc)=Helac_ipsgn(i2,i3)*Helac_ipsgn(i1,i23)
          ENDIF
          list(kcount,18,icc)=igluon
          IF(MOD(kcount,100).EQ.0)WRITE(nunit1,*)&
	      'I am calculating the ',kcount,'-th subamplitude'
          i_ex=0
			   ! If we haven't dealt with this situation exactly, we set i_ex=1.
	      IF(kcount.GT.1)THEN
	         IF(list(kcount-1,2,icc).NE.i0)i_ex=1
	         IF(list(kcount-1,3,icc).NE.m0)i_ex=1
	         IF(list(kcount-1,7,icc).NE.i1)i_ex=1
   	         IF(list(kcount-1,8,icc).NE.m1)i_ex=1
             IF(list(kcount-1,10,icc).NE.i2)i_ex=1
             IF(list(kcount-1,11,icc).NE.m2)i_ex=1
             IF(list(kcount-1,13,icc).NE.i3)i_ex=1
             IF(list(kcount-1,14,icc).NE.m3)i_ex=1
	      ELSEIF(kcount.EQ.1)THEN
             i_ex=1
	      ENDIF
          IF(igluon.GE.2)i_ex=0
          IF(i_ex.EQ.1)THEN
             IF(i3.EQ.0)ife(i0,m0)=ife(i0,m0)+ife(i1,m1)*ife(i2,m2)
             IF(i3.NE.0)ife(i0,m0)=ife(i0,m0)+ife(i1,m1)*ife(i2,m2)*ife(i3,m3)
          ENDIF
! end of list.h  
      ENDIF
      IF(ih(i1,m1).NE.0.AND.ih(i2,m2).NE.0.AND.zgsffr(an(m1),m2,-m0).NE.zero)THEN
          icole(i0,1:2)=icole(i2,1:2)
          ih(i0,m0)=1
          irou=7
          lpol=2
	           !           include 'list.h'         
          kcount=kcount+1
          list(kcount,1,icc)=irou
          list(kcount,2,icc)=i0
          list(kcount,3,icc)=m0
          list(kcount,7,icc)=i1
          list(kcount,8,icc)=m1
          list(kcount,10,icc)=i2
          list(kcount,11,icc)=m2
          list(kcount,13,icc)=i3
          list(kcount,14,icc)=m3
          list(kcount,16,icc)=lpol
          IF(i3.EQ.0)THEN
            list(kcount,17,icc)=Helac_ipsgn(i1,i2)
          ELSE
            i23=i2+i3
            list(kcount,17,icc)=Helac_ipsgn(i2,i3)*Helac_ipsgn(i1,i23)
          ENDIF
          list(kcount,18,icc)=igluon
          IF(MOD(kcount,100).EQ.0)WRITE(nunit1,*)&
	      'I am calculating the ',kcount,'-th subamplitude'
          i_ex=0
			   ! If we haven't dealt with this situation exactly, we set i_ex=1.
	      IF(kcount.GT.1)THEN
	         IF(list(kcount-1,2,icc).NE.i0)i_ex=1
	         IF(list(kcount-1,3,icc).NE.m0)i_ex=1
	         IF(list(kcount-1,7,icc).NE.i1)i_ex=1
   	         IF(list(kcount-1,8,icc).NE.m1)i_ex=1
             IF(list(kcount-1,10,icc).NE.i2)i_ex=1
             IF(list(kcount-1,11,icc).NE.m2)i_ex=1
             IF(list(kcount-1,13,icc).NE.i3)i_ex=1
             IF(list(kcount-1,14,icc).NE.m3)i_ex=1
	      ELSEIF(kcount.EQ.1)THEN
             i_ex=1
	      ENDIF
          IF(igluon.GE.2)i_ex=0
          IF(i_ex.EQ.1)THEN
             IF(i3.EQ.0)ife(i0,m0)=ife(i0,m0)+ife(i1,m1)*ife(i2,m2)
             IF(i3.NE.0)ife(i0,m0)=ife(i0,m0)+ife(i1,m1)*ife(i2,m2)*ife(i3,m3)
          ENDIF
! end of list.h 
      ENDIF
   ENDDO
ENDDO
END SUBROUTINE Helac_affs

! -----------------------------------------------------------------------      

SUBROUTINE Helac_sff(i0,le,m0,i1,i2,i3)
IMPLICIT NONE
INTEGER,INTENT(IN)::i1,i2,i3,m0,i0,le
INTEGER::lpol,m1,m2,irou,i23,i_ex
m3=0
lpol=0
DO m1=-12,-1
   DO m2=1,12
      IF((zgsffr(m0,m1,m2).NE.zero.OR.zgsffl(m0,m1,m2).NE.zero).AND.ih(i1,m1).NE.0&
       .AND.ih(i2,m2).NE.0)THEN
          IF(icole(i1,1).NE.icole(i2,2).OR.icole(i1,2).NE.icole(i2,1))CYCLE
          ih(i0,m0)=1
          irou=8
	           !           include 'list.h'         
          kcount=kcount+1
          list(kcount,1,icc)=irou
          list(kcount,2,icc)=i0
          list(kcount,3,icc)=m0
          list(kcount,7,icc)=i1
          list(kcount,8,icc)=m1
          list(kcount,10,icc)=i2
          list(kcount,11,icc)=m2
          list(kcount,13,icc)=i3
          list(kcount,14,icc)=m3
          list(kcount,16,icc)=lpol
          IF(i3.EQ.0)THEN
            list(kcount,17,icc)=Helac_ipsgn(i1,i2)
          ELSE
            i23=i2+i3
            list(kcount,17,icc)=Helac_ipsgn(i2,i3)*Helac_ipsgn(i1,i23)
          ENDIF
          list(kcount,18,icc)=igluon
          IF(MOD(kcount,100).EQ.0)WRITE(nunit1,*)&
	      'I am calculating the ',kcount,'-th subamplitude'
          i_ex=0
			   ! If we haven't dealt with this situation exactly, we set i_ex=1.
	      IF(kcount.GT.1)THEN
	         IF(list(kcount-1,2,icc).NE.i0)i_ex=1
	         IF(list(kcount-1,3,icc).NE.m0)i_ex=1
	         IF(list(kcount-1,7,icc).NE.i1)i_ex=1
   	         IF(list(kcount-1,8,icc).NE.m1)i_ex=1
             IF(list(kcount-1,10,icc).NE.i2)i_ex=1
             IF(list(kcount-1,11,icc).NE.m2)i_ex=1
             IF(list(kcount-1,13,icc).NE.i3)i_ex=1
             IF(list(kcount-1,14,icc).NE.m3)i_ex=1
	      ELSEIF(kcount.EQ.1)THEN
             i_ex=1
	      ENDIF
          IF(igluon.GE.2)i_ex=0
          IF(i_ex.EQ.1)THEN
             IF(i3.EQ.0)ife(i0,m0)=ife(i0,m0)+ife(i1,m1)*ife(i2,m2)
             IF(i3.NE.0)ife(i0,m0)=ife(i0,m0)+ife(i1,m1)*ife(i2,m2)*ife(i3,m3)
          ENDIF
! end of list.h 
       ENDIF
    ENDDO
ENDDO
END SUBROUTINE Helac_sff

! -----------------------------------------------------------------------      

SUBROUTINE Helac_s3(i0,le,m0,i1,i2,i3)
IMPLICIT NONE
INTEGER,INTENT(IN)::i1,i2,i3,m0,i0,le
INTEGER::lpol,m1,m2,irou,i23,i_ex,igo
m3=0
lpol=0
DO m1=41,44
   DO m2=41,44
       igo=1
       IF(m1.EQ.m2.AND.i1.GT.i2)igo=-1  !12
       IF(zgs3(m0,m1,m2).NE.zero.AND.ih(i1,m1).NE.0&
         .AND.ih(i2,m2).NE.0.AND.igo.EQ.1)THEN
           ih(i0,m0)=1
           irou=13
	           !           include 'list.h'         
           kcount=kcount+1
           list(kcount,1,icc)=irou
           list(kcount,2,icc)=i0
           list(kcount,3,icc)=m0
           list(kcount,7,icc)=i1
           list(kcount,8,icc)=m1
           list(kcount,10,icc)=i2
           list(kcount,11,icc)=m2
           list(kcount,13,icc)=i3
           list(kcount,14,icc)=m3
           list(kcount,16,icc)=lpol
           IF(i3.EQ.0)THEN
             list(kcount,17,icc)=Helac_ipsgn(i1,i2)
           ELSE
             i23=i2+i3
             list(kcount,17,icc)=Helac_ipsgn(i2,i3)*Helac_ipsgn(i1,i23)
           ENDIF
           list(kcount,18,icc)=igluon
           IF(MOD(kcount,100).EQ.0)WRITE(nunit1,*)&
	      'I am calculating the ',kcount,'-th subamplitude'
           i_ex=0
			   ! If we haven't dealt with this situation exactly, we set i_ex=1.
	       IF(kcount.GT.1)THEN
	         IF(list(kcount-1,2,icc).NE.i0)i_ex=1
	         IF(list(kcount-1,3,icc).NE.m0)i_ex=1
	         IF(list(kcount-1,7,icc).NE.i1)i_ex=1
   	         IF(list(kcount-1,8,icc).NE.m1)i_ex=1
             IF(list(kcount-1,10,icc).NE.i2)i_ex=1
             IF(list(kcount-1,11,icc).NE.m2)i_ex=1
             IF(list(kcount-1,13,icc).NE.i3)i_ex=1
             IF(list(kcount-1,14,icc).NE.m3)i_ex=1
	       ELSEIF(kcount.EQ.1)THEN
             i_ex=1
	       ENDIF
           IF(igluon.GE.2)i_ex=0
           IF(i_ex.EQ.1)THEN
             IF(i3.EQ.0)ife(i0,m0)=ife(i0,m0)+ife(i1,m1)*ife(i2,m2)
             IF(i3.NE.0)ife(i0,m0)=ife(i0,m0)+ife(i1,m1)*ife(i2,m2)*ife(i3,m3)
           ENDIF
! end of list.h 
       ENDIF
   ENDDO
ENDDO
END SUBROUTINE Helac_s3

! -----------------------------------------------------------------------      

SUBROUTINE Helac_s4(i0,le,m0,i1,i2,i3)
IMPLICIT NONE
INTEGER,INTENT(IN)::i1,i2,i3,m0,i0,le
INTEGER::lpol,m1,m2,irou,i23,i_ex,igo
lpol=0
DO m1=41,44
   DO m2=41,44
      DO m3=41,44
         igo=1
         IF(m1.EQ.m2.AND.i1.GT.i2)igo=-1  !12
         IF(m1.EQ.m3.AND.i1.GT.i3)igo=-1  !13
         IF(m2.EQ.m3.AND.i2.GT.i3)igo=-1  !23
         IF(zgs4(m0,m1,m2,m3).NE.zero.AND.ih(i1,m1).NE.0&
          .AND.ih(i2,m2).NE.0.AND.ih(i3,m3).NE.0.AND.igo.EQ.1)THEN
             ih(i0,m0)=1
             irou=14
	           !           include 'list.h'         
             kcount=kcount+1
             list(kcount,1,icc)=irou
             list(kcount,2,icc)=i0
             list(kcount,3,icc)=m0
             list(kcount,7,icc)=i1
             list(kcount,8,icc)=m1
             list(kcount,10,icc)=i2
             list(kcount,11,icc)=m2
             list(kcount,13,icc)=i3
             list(kcount,14,icc)=m3
             list(kcount,16,icc)=lpol
             IF(i3.EQ.0)THEN
               list(kcount,17,icc)=Helac_ipsgn(i1,i2)
             ELSE
               i23=i2+i3
               list(kcount,17,icc)=Helac_ipsgn(i2,i3)*Helac_ipsgn(i1,i23)
             ENDIF
             list(kcount,18,icc)=igluon
             IF(MOD(kcount,100).EQ.0)WRITE(nunit1,*)&
	        'I am calculating the ',kcount,'-th subamplitude'
             i_ex=0
			   ! If we haven't dealt with this situation exactly, we set i_ex=1.
	         IF(kcount.GT.1)THEN
	           IF(list(kcount-1,2,icc).NE.i0)i_ex=1
	           IF(list(kcount-1,3,icc).NE.m0)i_ex=1
	           IF(list(kcount-1,7,icc).NE.i1)i_ex=1
   	           IF(list(kcount-1,8,icc).NE.m1)i_ex=1
               IF(list(kcount-1,10,icc).NE.i2)i_ex=1
               IF(list(kcount-1,11,icc).NE.m2)i_ex=1
               IF(list(kcount-1,13,icc).NE.i3)i_ex=1
               IF(list(kcount-1,14,icc).NE.m3)i_ex=1
	         ELSEIF(kcount.EQ.1)THEN
               i_ex=1
	         ENDIF
             IF(igluon.GE.2)i_ex=0
             IF(i_ex.EQ.1)THEN
               IF(i3.EQ.0)ife(i0,m0)=ife(i0,m0)+ife(i1,m1)*ife(i2,m2)
               IF(i3.NE.0)ife(i0,m0)=ife(i0,m0)+ife(i1,m1)*ife(i2,m2)*ife(i3,m3)
             ENDIF
! end of list.h 
         ENDIF
      ENDDO
   ENDDO
ENDDO
END SUBROUTINE Helac_s4

! -----------------------------------------------------------------------      

SUBROUTINE Helac_svv(i0,le,m0,i1,i2,i3)
IMPLICIT NONE
INTEGER,INTENT(IN)::i1,i2,i3,m0,i0,le
INTEGER::lpol,m1,m2,irou,i23,i_ex,igo
m3=0
lpol=0
DO m1=iv1,iv2
   DO m2=iv1,iv2
      igo=1
      IF(m1.EQ.m2.AND.i1.GT.i2)igo=-1
      IF(zgsvv(m0,m1,m2).NE.zero.AND.ih(i1,m1).NE.0&
       .AND.ih(i2,m2).NE.0.AND.igo.EQ.1)THEN
             ih(i0,m0)=1
             irou=10
	           !           include 'list.h'         
             kcount=kcount+1
             list(kcount,1,icc)=irou
             list(kcount,2,icc)=i0
             list(kcount,3,icc)=m0
             list(kcount,7,icc)=i1
             list(kcount,8,icc)=m1
             list(kcount,10,icc)=i2
             list(kcount,11,icc)=m2
             list(kcount,13,icc)=i3
             list(kcount,14,icc)=m3
             list(kcount,16,icc)=lpol
             IF(i3.EQ.0)THEN
               list(kcount,17,icc)=Helac_ipsgn(i1,i2)
             ELSE
               i23=i2+i3
               list(kcount,17,icc)=Helac_ipsgn(i2,i3)*Helac_ipsgn(i1,i23)
             ENDIF
             list(kcount,18,icc)=igluon
             IF(MOD(kcount,100).EQ.0)WRITE(nunit1,*)&
	        'I am calculating the ',kcount,'-th subamplitude'
             i_ex=0
			   ! If we haven't dealt with this situation exactly, we set i_ex=1.
	         IF(kcount.GT.1)THEN
	           IF(list(kcount-1,2,icc).NE.i0)i_ex=1
	           IF(list(kcount-1,3,icc).NE.m0)i_ex=1
	           IF(list(kcount-1,7,icc).NE.i1)i_ex=1
   	           IF(list(kcount-1,8,icc).NE.m1)i_ex=1
               IF(list(kcount-1,10,icc).NE.i2)i_ex=1
               IF(list(kcount-1,11,icc).NE.m2)i_ex=1
               IF(list(kcount-1,13,icc).NE.i3)i_ex=1
               IF(list(kcount-1,14,icc).NE.m3)i_ex=1
	         ELSEIF(kcount.EQ.1)THEN
               i_ex=1
	         ENDIF
             IF(igluon.GE.2)i_ex=0
             IF(i_ex.EQ.1)THEN
               IF(i3.EQ.0)ife(i0,m0)=ife(i0,m0)+ife(i1,m1)*ife(i2,m2)
               IF(i3.NE.0)ife(i0,m0)=ife(i0,m0)+ife(i1,m1)*ife(i2,m2)*ife(i3,m3)
             ENDIF
! end of list.h 
      ENDIF
   ENDDO
ENDDO
END SUBROUTINE Helac_svv

! -----------------------------------------------------------------------      

SUBROUTINE Helac_ssv(i0,le,m0,i1,i2,i3)
IMPLICIT NONE
INTEGER,INTENT(IN)::i1,i2,i3,m0,i0,le
INTEGER::lpol,m1,m2,irou,i23,i_ex,igo
m3=0
lpol=0
DO m1=41,44
   DO m2=iv1,iv2
      IF(zgssv(m0,m1,m2).NE.zero.AND.ih(i1,m1).NE.0.AND.ih(i2,m2).NE.0)THEN
             ih(i0,m0)=1
             irou=15
	           !           include 'list.h'         
             kcount=kcount+1
             list(kcount,1,icc)=irou
             list(kcount,2,icc)=i0
             list(kcount,3,icc)=m0
             list(kcount,7,icc)=i1
             list(kcount,8,icc)=m1
             list(kcount,10,icc)=i2
             list(kcount,11,icc)=m2
             list(kcount,13,icc)=i3
             list(kcount,14,icc)=m3
             list(kcount,16,icc)=lpol
             IF(i3.EQ.0)THEN
               list(kcount,17,icc)=Helac_ipsgn(i1,i2)
             ELSE
               i23=i2+i3
               list(kcount,17,icc)=Helac_ipsgn(i2,i3)*Helac_ipsgn(i1,i23)
             ENDIF
             list(kcount,18,icc)=igluon
             IF(MOD(kcount,100).EQ.0)WRITE(nunit1,*)&
	        'I am calculating the ',kcount,'-th subamplitude'
             i_ex=0
			   ! If we haven't dealt with this situation exactly, we set i_ex=1.
	         IF(kcount.GT.1)THEN
	           IF(list(kcount-1,2,icc).NE.i0)i_ex=1
	           IF(list(kcount-1,3,icc).NE.m0)i_ex=1
	           IF(list(kcount-1,7,icc).NE.i1)i_ex=1
   	           IF(list(kcount-1,8,icc).NE.m1)i_ex=1
               IF(list(kcount-1,10,icc).NE.i2)i_ex=1
               IF(list(kcount-1,11,icc).NE.m2)i_ex=1
               IF(list(kcount-1,13,icc).NE.i3)i_ex=1
               IF(list(kcount-1,14,icc).NE.m3)i_ex=1
	         ELSEIF(kcount.EQ.1)THEN
               i_ex=1
	         ENDIF
             IF(igluon.GE.2)i_ex=0
             IF(i_ex.EQ.1)THEN
               IF(i3.EQ.0)ife(i0,m0)=ife(i0,m0)+ife(i1,m1)*ife(i2,m2)
               IF(i3.NE.0)ife(i0,m0)=ife(i0,m0)+ife(i1,m1)*ife(i2,m2)*ife(i3,m3)
             ENDIF
! end of list.h 
       ENDIF
   ENDDO
ENDDO
END SUBROUTINE Helac_ssv

! -----------------------------------------------------------------------      

SUBROUTINE Helac_vvss(i0,le,m0,i1,i2,i3)
IMPLICIT NONE
INTEGER,INTENT(IN)::i1,i2,i3,m0,i0,le
INTEGER::lpol,m1,m2,irou,i23,i_ex,igo
lpol=0
DO m1=iv1,iv2
   DO m2=41,44
      DO m3=41,44
           igo=1
           IF(m2.EQ.m3.AND.i2.GT.i3)igo=-1
           IF(zgvvss(m0,m1,m2,m3).NE.zero.AND.ih(i1,m1).NE.0&
           .AND.ih(i2,m2).NE.0.AND.ih(i3,m3).NE.0.AND.igo.EQ.1)THEN
             ih(i0,m0)=1
             irou=11
	           !           include 'list.h'         
             kcount=kcount+1
             list(kcount,1,icc)=irou
             list(kcount,2,icc)=i0
             list(kcount,3,icc)=m0
             list(kcount,7,icc)=i1
             list(kcount,8,icc)=m1
             list(kcount,10,icc)=i2
             list(kcount,11,icc)=m2
             list(kcount,13,icc)=i3
             list(kcount,14,icc)=m3
             list(kcount,16,icc)=lpol
             IF(i3.EQ.0)THEN
               list(kcount,17,icc)=Helac_ipsgn(i1,i2)
             ELSE
               i23=i2+i3
               list(kcount,17,icc)=Helac_ipsgn(i2,i3)*Helac_ipsgn(i1,i23)
             ENDIF
             list(kcount,18,icc)=igluon
             IF(MOD(kcount,100).EQ.0)WRITE(nunit1,*)&
	        'I am calculating the ',kcount,'-th subamplitude'
             i_ex=0
			   ! If we haven't dealt with this situation exactly, we set i_ex=1.
	         IF(kcount.GT.1)THEN
	           IF(list(kcount-1,2,icc).NE.i0)i_ex=1
	           IF(list(kcount-1,3,icc).NE.m0)i_ex=1
	           IF(list(kcount-1,7,icc).NE.i1)i_ex=1
   	           IF(list(kcount-1,8,icc).NE.m1)i_ex=1
               IF(list(kcount-1,10,icc).NE.i2)i_ex=1
               IF(list(kcount-1,11,icc).NE.m2)i_ex=1
               IF(list(kcount-1,13,icc).NE.i3)i_ex=1
               IF(list(kcount-1,14,icc).NE.m3)i_ex=1
	         ELSEIF(kcount.EQ.1)THEN
               i_ex=1
	         ENDIF
             IF(igluon.GE.2)i_ex=0
             IF(i_ex.EQ.1)THEN
               IF(i3.EQ.0)ife(i0,m0)=ife(i0,m0)+ife(i1,m1)*ife(i2,m2)
               IF(i3.NE.0)ife(i0,m0)=ife(i0,m0)+ife(i1,m1)*ife(i2,m2)*ife(i3,m3)
             ENDIF
! end of list.h 
           ENDIF
      ENDDO
   ENDDO
ENDDO
END SUBROUTINE Helac_vvss

! -----------------------------------------------------------------------      

SUBROUTINE Helac_ssvv(i0,le,m0,i1,i2,i3)
IMPLICIT NONE
INTEGER,INTENT(IN)::i1,i2,i3,m0,i0,le
INTEGER::lpol,m1,m2,irou,i23,i_ex,igo
lpol=0
DO m1=41,44
   DO m2=iv1,iv2
      DO m3=iv1,iv2
          igo=1
          IF(m2.EQ.m3.AND.i2.GT.i3)igo=-1
          IF(zgssvv(m0,m1,m2,m3).NE.zero.AND.ih(i1,m1).NE.0&
		  .AND.ih(i2,m2).NE.0.AND.ih(i3,m3).NE.0.AND.igo.EQ.1)THEN
          
             ih(i0,m0)=1
             irou=12
	           !           include 'list.h'         
             kcount=kcount+1
             list(kcount,1,icc)=irou
             list(kcount,2,icc)=i0
             list(kcount,3,icc)=m0
             list(kcount,7,icc)=i1
             list(kcount,8,icc)=m1
             list(kcount,10,icc)=i2
             list(kcount,11,icc)=m2
             list(kcount,13,icc)=i3
             list(kcount,14,icc)=m3
             list(kcount,16,icc)=lpol
             IF(i3.EQ.0)THEN
               list(kcount,17,icc)=Helac_ipsgn(i1,i2)
             ELSE
               i23=i2+i3
               list(kcount,17,icc)=Helac_ipsgn(i2,i3)*Helac_ipsgn(i1,i23)
             ENDIF
             list(kcount,18,icc)=igluon
             IF(MOD(kcount,100).EQ.0)WRITE(nunit1,*)&
	        'I am calculating the ',kcount,'-th subamplitude'
             i_ex=0
			   ! If we haven't dealt with this situation exactly, we set i_ex=1.
	         IF(kcount.GT.1)THEN
	           IF(list(kcount-1,2,icc).NE.i0)i_ex=1
	           IF(list(kcount-1,3,icc).NE.m0)i_ex=1
	           IF(list(kcount-1,7,icc).NE.i1)i_ex=1
   	           IF(list(kcount-1,8,icc).NE.m1)i_ex=1
               IF(list(kcount-1,10,icc).NE.i2)i_ex=1
               IF(list(kcount-1,11,icc).NE.m2)i_ex=1
               IF(list(kcount-1,13,icc).NE.i3)i_ex=1
               IF(list(kcount-1,14,icc).NE.m3)i_ex=1
	         ELSEIF(kcount.EQ.1)THEN
               i_ex=1
	         ENDIF
             IF(igluon.GE.2)i_ex=0
             IF(i_ex.EQ.1)THEN
               IF(i3.EQ.0)ife(i0,m0)=ife(i0,m0)+ife(i1,m1)*ife(i2,m2)
               IF(i3.NE.0)ife(i0,m0)=ife(i0,m0)+ife(i1,m1)*ife(i2,m2)*ife(i3,m3)
             ENDIF
! end of list.h 
          ENDIF
      ENDDO
   ENDDO
ENDDO
END SUBROUTINE Helac_ssvv

! -----------------------------------------------------------------------      

SUBROUTINE Helac_redo(ic1,ino)
IMPLICIT NONE
INTEGER::maxmul=0
INTEGER,INTENT(OUT)::ino
INTEGER,INTENT(INOUT)::ic1
INTEGER::k,j0,n0,j1,n1,j2,n2,j3,n3,k1,kk,jj,ii1,mm1,ip,i,inofeyn,is,is1,istat
PRINT *,'ic2=',icc  
!kmax=0
!imax=0
IF(GLOBALINIT_redo.EQ.0)THEN
	kmax=0
	imax=0
	maxmul=0
ENDIF
ino=0
IF(kcount.EQ.0)THEN
   icc=icc+1 
   ino=1
!ic1=ic1+1
!kcou(ic1)=0
   WRITE(nunit1,*)' '
   WRITE(nunit1,'(A7,I6,A24,I8,A15)')&
      'For the',icc-1,'  colour conf. there are',0,'  subamplitudes'
   RETURN
ENDIF
! begin of note1
DO k=kcount,1,-1
   j0=list(k,2,icc)
   n0=list(k,3,icc)
   j1=list(k,7,icc)
   n1=list(k,8,icc)
   j2=list(k,10,icc)
   n2=list(k,11,icc)
   j3=list(k,13,icc)
   n3=list(k,14,icc)
   ! choose the suitable ones
   IF(j0.EQ.2**n-2.AND.n0.EQ.an(ifl(1)))THEN
      jxx(k)=1
      ih1(j1,n1)=1
      ih1(j2,n2)=1
      IF(j3.NE.0)ih1(j3,n3)=1
   ELSE
      IF(ih1(j0,n0).EQ.1)THEN
         jxx(k)=1
         ih1(j1,n1)=1
         ih1(j2,n2)=1
         IF(j3.NE.0)ih1(j3,n3)=1
      ENDIF
   ENDIF
ENDDO
! end of note1
! begin of note2
k1=0
ic1=ic1+1
DO k=1,kcount
   IF(jxx(k).EQ.1)THEN
       k1=k1+1
       list(k1,1:18,ic1)=list(k,1:18,icc)
   ENDIF
ENDDO
! end of note2
! begin of note3
kcou(ic1)=k1
IF(ic1.EQ.ncc)THEN
   kk=0
   DO jj=1,ncc
      kk=kk+kcou(jj)
   ENDDO
   IF(kk.EQ.0)THEN
      WRITE(*,*)'There are no contributions to this process'
      STOP
   ENDIF
ENDIF
icc=icc+1          
IF(k1.EQ.0)THEN
   ic1=ic1-1
   ino=1
   WRITE(nunit1,'(a7,i6,a24,i8,a15)')&
    'For the',icc-1,'  colour conf. there are',k1,'  subamplitudes'  
   RETURN
ENDIF
! end of note3
! begin of note4
ii1=0
mm1=0
ip=n
ipp(0:2**n-2,-12:44,1:3)=0
DO i=1,n
  ipp(list1(i,1),list1(i,2),1)=i
ENDDO
! end of note4
! begin of note5
inofeyn=0
DO k=1,k1
   is=0
   is1=0
   IF(k.GT.1)THEN
      IF(list(k,2,ic1).EQ.list(k-1,2,ic1))is=is+1
      IF(list(k,3,ic1).EQ.list(k-1,3,ic1))is=is+1     
      IF(list(k,7,ic1).EQ.list(k-1,7,ic1))is=is+1     
      IF(list(k,8,ic1).EQ.list(k-1,8,ic1))is=is+1 
      IF(list(k,10,ic1).EQ.list(k-1,10,ic1))is=is+1   
      IF(list(k,11,ic1).EQ.list(k-1,11,ic1))is=is+1   
      IF(list(k,13,ic1).EQ.list(k-1,13,ic1))is=is+1    
      IF(list(k,14,ic1).EQ.list(k-1,14,ic1))is=is+1       
   ENDIF
        
   IF(is.NE.8.AND.is1.NE.1)& !8)       
   ipp(list(k,2,ic1),list(k,3,ic1),2)=ipp(list(k,2,ic1),list(k,3,ic1),2)+1
   IF(ipp(list(k,2,ic1),list(k,3,ic1),1).EQ.0)THEN
      ip=ip+1
      ipp(list(k,2,ic1),list(k,3,ic1),1)=ip
   ENDIF
   IF(ipp(list(k,2,ic1),list(k,3,ic1),1).NE.0)THEN
      IF(is.NE.8.AND.is1.NE.1)THEN  !8)then
         ipp(list(k,2,ic1),list(k,3,ic1),3)=&
         ipp(list(k,2,ic1),list(k,3,ic1),3)+1
         list(k,5,ic1)=ipp(list(k,2,ic1),list(k,3,ic1),3)
      ENDIF
   ENDIF
   ii1=list(k,2,ic1)
   mm1=list(k,3,ic1)
ENDDO
! end of note5
! begion of note6
ii1=0
mm1=0
DO k=1,k1
   list(k,4,ic1)=ipp(list(k,2,ic1),list(k,3,ic1),1)
   list(k,6,ic1)=ipp(list(k,2,ic1),list(k,3,ic1),2)
   list(k,9,ic1)=ipp(list(k,7,ic1),list(k,8,ic1),1)
   list(k,12,ic1)=ipp(list(k,10,ic1),list(k,11,ic1),1)
   list(k,15,ic1)=ipp(list(k,13,ic1),list(k,14,ic1),1)
   IF(list(k,6,ic1).GT.maxmul)maxmul=list(k,6,ic1)
ENDDO
! end of note6
WRITE(nunit1,'(A7,2I6,A24,I8,A15)')'For the',icc-1,ic1,&
     '  colour conf. there are',k1,'  subamplitudes'        
IF(k1.GT.kmax)kmax=k1                   ! kmax denotes the maximum subamplitudes via colors
IF(list(k1,4,ic1).GT.imax)imax=list(k1,4,ic1)  ! imax denotes the k1's maximum splitting
WRITE(nunit1,*)'-----------------------------------------------------------------------------'
WRITE(nunit1,*)'0 index, 1 rout., 2 i0, 3 m0, 4 n_(i0,m0), 5 s_(i0,m0), 6 s^max_(i0,m0)'
WRITE(nunit1,*)'7 i1, 8 m1, 9 n_(i1,m1), 10 i2, 11 m2, 12 n_(i2,m2)'
WRITE(nunit1,*)'13 i3, 14 m3, 15 n_(i3,m3)'
WRITE(nunit1,*)'16 chiral of fermion (L1R2), 17 Sign(i1,i2,i3), 18 igloun'
WRITE(nunit1,'(19I4)')(k,list(k,1:18,ic1),k=1,k1)
WRITE(nunit1,*)'The number of Feynman graphs = ',ife(2**n-2,an(ifl(1)))
WRITE(nunit1,*)'------------------------------------------------------------------------------'
WRITE(nunit1,*)' '
IF(k1.LE.ngues.AND.inofeyn.EQ.0)THEN
   ALLOCATE(ico(k1,18),STAT=istat)
   ico(1:k1,1:18)=list(1:k1,1:18,ic1)
   CALL Helac_feyn(n,k1,maxmul,ico)
   GLOBALINIT_feyn=1
   DEALLOCATE(ico)
ENDIF
END SUBROUTINE Helac_redo

! -----------------------------------------------------------------------   
           
SUBROUTINE Helac_setlist(newncc)
IMPLICIT NONE
INTEGER,INTENT(IN)::newncc
INTEGER::istat
ncc=newncc
DEALLOCATE(ih,ih1,icole,ipp,jxx,ife)
ALLOCATE(list_tmp(kmax,18,ncc),STAT=istat)
list_tmp(1:kmax,1:18,1:ncc)=list(1:kmax,1:18,1:ncc)
DEALLOCATE(list)
ALLOCATE(list(kmax,18,ncc),stat=istat)
list(1:kmax,1:18,1:ncc)=list_tmp(1:kmax,1:18,1:ncc)
DEALLOCATE(list_tmp)
IF(repeat.EQ.1)THEN
   WRITE(21,*)kmax,ncc,imax
   WRITE(21,*)list1
   WRITE(21,*)list
   WRITE(21,*)kcou(1:ncc)
ENDIF
IF(ALLOCATED(zy))THEN
	DEALLOCATE(zy)
ENDIF
IF(ALLOCATED(zyq))THEN
	DEALLOCATE(zyq)
ENDIF
ALLOCATE(zy(0:2**Pwavenum-1,imax,4),zyq(imax,4),STAT=istat)
END SUBROUTINE Helac_setlist

! -----------------------------------------------------------------------      

SUBROUTINE Helac_getlist()
IMPLICIT NONE
INTEGER::istat
READ(21,*)kmax,ncc,imax
WRITE(*,*)'kmax,ncc,imax',kmax,ncc,imax
IF(ALLOCATED(list))THEN
	DEALLOCATE(list)
ENDIF
IF(ALLOCATED(list1))THEN
	DEALLOCATE(list1)
ENDIF
IF(ALLOCATED(kcou))THEN
	DEALLOCATE(kcou)
ENDIF
ALLOCATE(list(kmax,18,ncc),list1(n,3),kcou(ncc),STAT=istat)
READ(21,*)list1
READ(21,*)list
READ(21,*)kcou
IF(ALLOCATED(zy))THEN
	DEALLOCATE(zy)
ENDIF
IF(ALLOCATED(zyq))THEN
	DEALLOCATE(zyq)
ENDIF
ALLOCATE(zy(0:2**Pwavenum-1,imax,4),zyq(imax,4),stat=istat)
END SUBROUTINE Helac_getlist

! -----------------------------------------------------------------------      

SUBROUTINE Helac_iniqq()
! corrected on 6.7.2006 (Worek) in order to get correct results for repeat=2
IMPLICIT NONE
INTEGER::k,j,ijij,index,index3,kk2,spin
INTEGER,DIMENSION(3)::index2
REAL(KIND=DBL)::r,r1     !,r2,r3,rquark
!COMPLEX(KIND=DBL),DIMENSION(9,4)::zarray1,zarray2
COMPLEX(KIND=DBL),DIMENSION(3,4)::zexex1
COMPLEX(KIND=DBL),DIMENSION(3,3,4)::zexex2
IF(init_pi.EQ.0)THEN
  CALL Helac_mypi(pi)
  zi=DCMPLX(0,1d0)
  init_pi=1
ENDIF
kk2=1
DO k=1,n
   DO j=1,5
        zq1(j)=zq(k,j)       ! common_mom.h
                             !zq(1:n,1)=(p0+pz,pz);zq(1:n,2)=(p0-pz,pt);zq(1:n,3)=(px,py);
                             !zq(1:n,4)=(px,-py);zq(1:n,4)=(m,p)
   ENDDO
   z1(1:4)=0
   z2(1:4)=0
   z3(1:4)=0
   z0(0:2**Pwavenum-1,1:4)=0
!   zarray1(1:9,1:4)=0
!   zarray2(1:9,1:4)=0
!  goto(41,42,43,44,45,46,47)list1(k,3)
   SELECT CASE(list1(k,3))
   CASE(1)    ! A or g      
! 41      continue
! io  1=incoming -1=outgoing
      zq1(2)=DCMPLX(DREAL(zq1(2))*io(k),DIMAG(zq1(2)))! p0-pz+i*pt  incoming -p0+pz+i*pt outgoing
      zq1(1)=zq1(1)*io(k)                  ! po+pz+i*pz incoming -po-pz-i*pz outgoing
      zq1(3)=zq1(3)*io(k)                  ! px+i*py incoming -px-i*py outgoing
      zq1(4)=zq1(4)*io(k)                  ! px-i*py incoming -px+i*py outgoing

      IF(iranhel.EQ.0)THEN                  ! iranhel is a variable in common_new.h
         IF(ipol(kk2).EQ.1)CALL Helac_eminus(z0(0,1:4),zq1)
         IF(ipol(kk2).EQ.2)CALL Helac_eplus(z0(0,1:4),zq1) !ipol(k)=2:rigth-handed 
                                                  !ipol(k)=1:left-handed 
                                                  !ipol(k)=3:long
      ELSE
!        r=Helac_rnmy(0)*2*pi
		IF(Qhelran(kk2).EQ.-1)THEN
			r=Helac_rnmy(0)*2*pi
			Qhelran(kk2)=r
		ELSE
			r=Qhelran(kk2)
		ENDIF
        CALL Helac_eminus(z1,zq1)
        CALL Helac_eplus(z2,zq1)
        z0(0,1:4)=(DCOS(r)+zi*DSIN(r))*z1(1:4)+(DCOS(r)-zi*DSIN(r))*z2(1:4)
      ENDIF
	  kk2=kk2+1
!       goto50
! 42      continue
   CASE(2)   ! W,Z 
      zq1(2)=DCMPLX(DREAL(zq1(2))*io(k),DIMAG(zq1(2)))
      zq1(1)=zq1(1)*io(k)
      zq1(3)=zq1(3)*io(k)
      zq1(4)=zq1(4)*io(k)

      IF(iranhel.EQ.0)THEN
        IF(ipol(kk2).EQ.1)CALL  Helac_eminus(z0(0,1:4),zq1)
        IF(ipol(kk2).EQ.2)CALL  Helac_eplus(z0(0,1:4),zq1)
        IF(ipol(kk2).EQ.3)CALL  Helac_elong(z0(0,1:4),zq1)
      ELSE
!        r=Helac_rnmy(0)*2*pi
		IF(Qhelran(kk2).EQ.-1)THEN
			r=Helac_rnmy(0)*2*pi
			Qhelran(kk2)=r
		ELSE
			r=Qhelran(kk2)
		ENDIF
        CALL Helac_eminus(z1,zq1)
        CALL Helac_eplus(z2,zq1)
        CALL Helac_elong(z3,zq1)
        z0(0,1:4)=(DCOS(r)+zi*DSIN(r))*z1(1:4)+(DCOS(r)-zi*DSIN(r))*z2(1:4)+z3(1:4)
      ENDIF
	  kk2=kk2+1
!        goto50
! 43      continue
   CASE(3)   !incoming fermion      
      IF(iranhel.EQ.0)THEN
        IF(ipol(kk2).EQ.2)CALL Helac_uplus(z0(0,1:4),zq1) 
        IF(ipol(kk2).EQ.1)CALL Helac_uminus(z0(0,1:4),zq1)
      ELSE
!        r=Helac_rnmy(0)*2*pi
		IF(Qhelran(kk2).EQ.-1)THEN
			r=Helac_rnmy(0)*2*pi
			Qhelran(kk2)=r
		ELSE
			r=Qhelran(kk2)
		ENDIF
        CALL Helac_uplus(z1,zq1)
        CALL Helac_uminus(z2,zq1)
        z0(0,1:4)=(DCOS(r)+zi*DSIN(r))*z1(1:4)+(DCOS(r)-zi*DSIN(r))*z2(1:4)
      ENDIF
	  kk2=kk2+1
!        goto50
! 44      continue
    CASE(4)  ! outgoing fermion
	  IF(Quarkonium3(k).EQ.0)THEN 
		IF(iranhel.EQ.0)THEN
			IF(ipol(kk2).EQ.2)CALL Helac_ubplus(z0(0,1:4),zq1) 
			IF(ipol(kk2).EQ.1)CALL Helac_ubminus(z0(0,1:4),zq1) 
		ELSE
!			r=Helac_rnmy(0)*2*pi
			IF(Qhelran(kk2).EQ.-1)THEN
				r=Helac_rnmy(0)*2*pi
				Qhelran(kk2)=r
			ELSE
				r=Qhelran(kk2)
			ENDIF
			CALL Helac_ubplus(z1,zq1)
			CALL Helac_ubminus(z2,zq1)
			z0(0,1:4)=(DCOS(r)+zi*DSIN(r))*z1(1:4)+(DCOS(r)-zi*DSIN(r))*z2(1:4)
		ENDIF
		kk2=kk2+1
	  ! Quarkoniums
	  ELSE
		CALL PlusZq(zPQQ(1:5),zq1(1:5),zq(k+1,1:5))
		! zero derivation
		IF(iranhel.EQ.1.OR.iranhel.EQ.0)THEN
			CALL Qeveryhel(index2,ipol(kk2),kk2)
			IF(index2(1).EQ.2)CALL Helac_Qubplus(z0(0,1:4),zq1,zPQQ,0,index2(2))
			IF(index2(1).EQ.1)CALL Helac_Qubminus(z0(0,1:4),zq1,zPQQ,0,index2(2))
			! one derivation
			IF(QN3PJF(kk2).OR.QN1P1F(kk2))THEN
				IF(index2(1).EQ.2)CALL Helac_Qubplus(z0(2**(PwaveQ(kk2)-1),1:4),zq1,zPQQ,1,index2(2))
				IF(index2(1).EQ.1)CALL Helac_Qubminus(z0(2**(PwaveQ(kk2)-1),1:4),zq1,zPQQ,1,index2(2))
			ENDIF
		ELSEIF(iranhel.GE.2)THEN
			CALL Qeveryhel(index2,ipol(kk2),kk2)
			IF(index2(1).EQ.2)CALL Helac_Qubplus(z0(0,1:4),zq1,zPQQ,0,index2(2))
			IF(index2(1).EQ.1)CALL Helac_Qubminus(z0(0,1:4),zq1,zPQQ,0,index2(2))
			! one derivation
			IF(QN3PJF(kk2).OR.QN1P1F(kk2))THEN
				IF(Qhelran(kk2).EQ.-1)THEN
					r=Helac_rnmy(0)*2*pi
					Qhelran(kk2)=r
				ELSE
					r=Qhelran(kk2)
				ENDIF
				zexex1(1:3,1:4)=0
				IF(index2(1).EQ.2)THEN
					DO index=1,3
						CALL Helac_Qubplus(zexex1(index,1:4),zq1,zPQQ,1,index)
					ENDDO
				ENDIF
				IF(index2(1).EQ.1)THEN
					DO index=1,3
						CALL Helac_Qubminus(zexex1(index,1:4),zq1,zPQQ,1,index)
					ENDDO
				ENDIF
				z0(2**(PwaveQ(kk2)-1),1:4)=(DCOS(r)+zi*DSIN(r))*zexex1(1,1:4)+&
				(DCOS(r)-zi*DSIN(r))*zexex1(2,1:4)+zexex1(3,1:4)
			! zero derivation same with 1P1 and 3PJ
!			ELSEIF(QN3S1F(kk2).OR.QN1S0F(kk2))THEN
!				IF(index2(1).EQ.2)CALL Helac_Qubplus(z0,zq1,zPQQ,0,1)
!				IF(index2(1).EQ.1)CALL Helac_Qubminus(z0,zq1,zPQQ,0,1)
			ENDIF
!		ELSEIF(iranhel.EQ.3)THEN
		ENDIF
!			rquark=r1
			! zero derivation
			! same with zero derivation in 1P1 and 3PJ
!			IF(QN1S0F(kk2).OR.QN3S1F(kk2))THEN
!				CALL Helac_Qubplus(zarray1(1,1:4),zq1,zPQQ,0,1)
!				CALL Helac_Qubminus(zarray2(1,1:4),zq1,zPQQ,0,1)
!				z0(1:4)=(DCOS(rquark)+zi*DSIN(rquark))*zarray1(1,1:4)+&
!				(DCOS(rquark)-zi*DSIN(rquark))*zarray2(1,1:4)
!			ELSEIF(QN1P1F(kk2).OR.QN3PJF(kk2))THEN
			! one derivation
!				DO index=1,3
!					CALL Helac_Qubplus(zarray1(index,1:4),zq1,zPQQ,1,index)
!					CALL Helac_Qubminus(zarray2(index,1:4),zq1,zPQQ,1,index)
!				ENDDO
!				r2=Helac_rnmy(0)*2*pi
!				Qhelran(kk2)=r2
!				zarray1(1,1:4)=(DCOS(r2)+zi*DSIN(r2))*zarray1(1,1:4)+&
!				(DCOS(r2)-zi*DSIN(r2))*zarray1(2,1:4)+zarray1(3,1:4)
!				zarray2(1,1:4)=(DCOS(r2)+zi*DSIN(r2))*zarray2(1,1:4)+&
!				(DCOS(r2)-zi*DSIN(r2))*zarray2(2,1:4)+zarray2(3,1:4)
!				z0(1:4)=(DCOS(rquark)+zi*DSIN(rquark))*zarray1(1,1:4)+&
!				(DCOS(rquark)-zi*DSIN(rquark))*zarray2(1,1:4)
!			ENDIF
!		ENDIF
	  ENDIF
!        goto50
! 45      continue
    CASE(5)  ! incoming antifermion
	  IF(iranhel.EQ.0)THEN
		IF(ipol(kk2).EQ.2)CALL Helac_vbplus(z0(0,1:4),zq1)
		IF(ipol(kk2).EQ.1)CALL Helac_vbminus(z0(0,1:4),zq1)
	  ELSE
!		r=Helac_rnmy(0)*2*pi
		IF(Qhelran(kk2).EQ.-1)THEN
			r=Helac_rnmy(0)*2*pi
			Qhelran(kk2)=r
		ELSE
			r=Qhelran(kk2)
		ENDIF
		CALL Helac_vbplus(z1,zq1)
		CALL Helac_vbminus(z2,zq1)
		z0(0,1:4)=(DCOS(r)+zi*DSIN(r))*z1(1:4)+(DCOS(r)-zi*DSIN(r))*z2(1:4)
	  ENDIF	
	  kk2=kk2+1
!        goto50
! 46      continue
    CASE(6)  ! outgoing antifermion
	  IF(Quarkonium3(k).EQ.0)THEN
		IF(iranhel.EQ.0)THEN
			IF(ipol(kk2).EQ.2)CALL Helac_vplus(z0(0,1:4),zq1) 
			IF(ipol(kk2).EQ.1)CALL Helac_vminus(z0(0,1:4),zq1)
		ELSE
!			r=Helac_rnmy(0)*2*pi
			IF(Qhelran(kk2).EQ.-1)THEN
				r=Helac_rnmy(0)*2*pi
				Qhelran(kk2)=r
			ELSE
				r=Qhelran(kk2)
			ENDIF
			CALL Helac_vplus(z1,zq1)
			CALL Helac_vminus(z2,zq1)
			z0(0,1:4)=(DCOS(r)+zi*DSIN(r))*z1(1:4)+(DCOS(r)-zi*DSIN(r))*z2(1:4)
		ENDIF
	  ELSE
		CALL PlusZq(zPQQ(1:5),zq1(1:5),zq(k-1,1:5))
		IF(QN1S0F(kk2).OR.QN1P1F(kk2))THEN
			spin=0
		ELSEIF(QN3S1F(kk2).OR.QN3PJF(kk2))THEN
			spin=1
		ENDIF
		! zero derivation
		CALL Qeveryhel(index2,ipol(kk2),kk2)
		IF(iranhel.EQ.0.OR.iranhel.EQ.1)THEN
			IF(index2(1).EQ.2)CALL Helac_Qvplus(z0(0,1:4),zq1,zPQQ,spin,0,index2(2),index2(3))
			IF(index2(1).EQ.1)CALL Helac_Qvminus(z0(0,1:4),zq1,zPQQ,spin,0,index2(2),index2(3))
			! one derivation
			IF(QN3PJF(kk2).OR.QN1P1F(kk2))THEN
				IF(index2(1).EQ.2)CALL Helac_Qvplus(z0(2**(PwaveQ(kk2)-1),1:4),zq1,zPQQ,spin,1,index2(2),index2(3))
				IF(index2(1).EQ.1)CALL Helac_Qvminus(z0(2**(PwaveQ(kk2)-1),1:4),zq1,zPQQ,spin,1,index2(2),index2(3))
			ENDIF
		ELSEIF(iranhel.EQ.2.OR.(iranhel.EQ.3.AND.spin.EQ.0))THEN
			IF(index2(1).EQ.2)CALL Helac_Qvplus(z0(0,1:4),zq1,zPQQ,spin,0,index2(2),index2(3))
			IF(index2(1).EQ.1)CALL Helac_Qvminus(z0(0,1:4),zq1,zPQQ,spin,0,index2(2),index2(3))
			IF(QN3PJF(kk2).OR.QN1P1F(kk2))THEN
			! one derivation
				IF(Qhelran(kk2).EQ.-1)THEN
					r=Helac_rnmy(0)*2*pi
					Qhelran(kk2)=r
				ELSE
					r=Qhelran(kk2)
				ENDIF
				IF(index2(1).EQ.2)THEN
					DO index=1,3
						CALL Helac_Qvplus(zexex1(index,1:4),zq1,zPQQ,spin,1,index,index2(3))
					ENDDO
				ENDIF
				IF(index2(1).EQ.1)THEN
					DO index=1,3
						CALL Helac_Qvminus(zexex1(index,1:4),zq1,zPQQ,spin,1,index,index2(3))
					ENDDO
				ENDIF
				z0(2**(PwaveQ(kk2)-1),1:4)=(DCOS(r)+zi*DSIN(r))*zexex1(1,1:4)+&
				(DCOS(r)-zi*DSIN(r))*zexex1(2,1:4)+zexex1(3,1:4)
			! zero derivation same with 1P1 and 3PJ
!			ELSEIF(QN1S0F(kk2).OR.QN3S1F(kk2))THEN
!				IF(index2(1).EQ.2)CALL Helac_Qvplus(z0,zq1,zPQQ,spin,0,1,index2(3))
!				IF(index2(1).EQ.1)CALL Helac_Qvminus(z0,zq1,zPQQ,spin,0,1,index2(3))								
			ENDIF
		ELSEIF(iranhel.EQ.3.AND.spin.EQ.1)THEN
			IF(QSzran(kk2).EQ.-1)THEN
				r1=Helac_rnmy(0)*2*pi
				QSzran(kk2)=r1
			ELSE
				r1=QSzran(kk2)
			ENDIF
			IF(index2(1).EQ.2)THEN
				DO index=1,3
					CALL Helac_Qvplus(zexex1(index,1:4),zq1,zPQQ,spin,0,1,index)
				ENDDO
			ENDIF
			IF(index2(1).EQ.1)THEN
				DO index=1,3
					CALL Helac_Qvminus(zexex1(index,1:4),zq1,zPQQ,spin,0,1,index)
				ENDDO
			ENDIF
			z0(0,1:4)=(DCOS(r1)+zi*DSIN(r1))*zexex1(1,1:4)+&
			(DCOS(r1)-zi*DSIN(r1))*zexex1(2,1:4)+zexex1(3,1:4)
			IF(QN3PJF(kk2).OR.QN1P1F(kk2))THEN
			! one derivation
				IF(Qhelran(kk2).EQ.-1)THEN
					r=Helac_rnmy(0)*2*pi
					Qhelran(kk2)=r
				ELSE
					r=Qhelran(kk2)
				ENDIF
				IF(index2(1).EQ.2)THEN
					DO index=1,3
						DO index3=1,3
							CALL Helac_Qvplus(zexex2(index,index3,1:4),zq1,zPQQ,spin,1,index,index3)
						ENDDO
						zexex1(index,1:4)=(DCOS(r1)+zi*DSIN(r1))*zexex2(index,1,1:4)+&
						(DCOS(r1)-zi*DSIN(r1))*zexex2(index,2,1:4)+zexex2(index,3,1:4)
					ENDDO
				ENDIF
				IF(index2(1).EQ.1)THEN
					DO index=1,3
						DO index3=1,3
							CALL Helac_Qvminus(zexex2(index,index3,1:4),zq1,zPQQ,spin,1,index,index3)
						ENDDO
						zexex1(index,1:4)=(DCOS(r1)+zi*DSIN(r1))*zexex2(index,1,1:4)+&
						(DCOS(r1)-zi*DSIN(r1))*zexex2(index,2,1:4)+zexex2(index,3,1:4)
					ENDDO
				ENDIF
				z0(2**(PwaveQ(kk2)-1),1:4)=(DCOS(r)+zi*DSIN(r))*zexex1(1,1:4)+&
				(DCOS(r)-zi*DSIN(r))*zexex1(2,1:4)+zexex1(3,1:4)
			ENDIF
		ENDIF
			! zero derivation
!			IF(QN1S0F(kk2))THEN
			! same with zero derivation in 1P1
!				CALL Helac_Qvplus(zarray1(1,1:4),zq1,zPQQ,spin,0,1,1)
!				CALL Helac_Qvminus(zarray2(1,1:4),zq1,zPQQ,spin,0,1,1)
!				z0(1:4)=(DCOS(rquark)-zi*DSIN(rquark))*zarray1(1,1:4)+&
!				(DCOS(rquark)+zi*DSIN(rquark))*zarray2(1,1:4)
!			ELSEIF(QN3S1F(kk2))THEN
			! same with zero derivation in 3PJ
!				DO index=1,3
!					CALL Helac_Qvplus(zarray1(index,1:4),zq1,zPQQ,spin,0,1,index)
!					CALL Helac_Qvminus(zarray2(index,1:4),zq1,zPQQ,spin,0,1,index)
!				ENDDO
!				r2=Helac_rnmy(0)*2*pi
!				zarray1(1,1:4)=(DCOS(r2)+zi*DSIN(r2))*zarray1(1,1:4)+&
!				(DCOS(r2)-zi*DSIN(r2))*zarray1(2,1:4)+zarray1(3,1:4)
!				zarray2(1,1:4)=(DCOS(r2)+zi*DSIN(r2))*zarray2(1,1:4)+&
!				(DCOS(r2)-zi*DSIN(r2))*zarray2(2,1:4)+zarray2(3,1:4)
!				z0(1:4)=(DCOS(rquark)-zi*DSIN(rquark))*zarray1(1,1:4)+&
!				(DCOS(rquark)+zi*DSIN(rquark))*zarray2(1,1:4)
!			ELSEIF(QN1P1F(kk2))THEN
			! one derivation
!				DO index=1,3
!					CALL Helac_Qvplus(zarray1(index,1:4),zq1,zPQQ,spin,1,index,1)
!					CALL Helac_Qvminus(zarray2(index,1:4),zq1,zPQQ,spin,1,index,1)
!				ENDDO
!				IF(Qhelran(kk2).EQ.-1)THEN
!					r2=Helac_rnmy(0)*2*pi
!					Qhelran(kk2)=r2
!				ELSE
!					r2=Qhelran(kk2)
!				ENDIF
!				zarray1(1,1:4)=(DCOS(r2)+zi*DSIN(r2))*zarray1(1,1:4)+&
!				(DCOS(r2)-zi*DSIN(r2))*zarray1(2,1:4)+zarray1(3,1:4)
!				zarray2(1,1:4)=(DCOS(r2)+zi*DSIN(r2))*zarray2(1,1:4)+&
!				(DCOS(r2)-zi*DSIN(r2))*zarray2(2,1:4)+zarray2(3,1:4)
!				z0(1:4)=(DCOS(rquark)-zi*DSIN(rquark))*zarray1(1,1:4)+&
!				(DCOS(rquark)+zi*DSIN(rquark))*zarray2(1,1:4)
!			ELSEIF(QN3PJF(kk2))THEN
			! one derivation
!				DO index3=1,3
!				  DO index=1,3
!					CALL Helac_Qvplus(zarray1(3*index3+index-3,1:4),zq1,zPQQ,spin,1,index,index3)
!					CALL Helac_Qvminus(zarray2(3*index3+index-3,1:4),zq1,zPQQ,spin,1,index,index3)
!				  ENDDO
!				ENDDO
!				IF(Qhelran(kk2).EQ.-1)THEN
!					r2=Helac_rnmy(0)*2*pi
!					Qhelran(kk2)=r2
!				ELSE
!					r2=Qhelran(kk2)
!				ENDIF
!				r3=Helac_rnmy(0)*2*pi
				! S_z=1
!				zexex1(1,1:4)=(DCOS(r2)+zi*DSIN(r2))*zarray1(1,1:4)+ &
!				(DCOS(r2)-zi*DSIN(r2))*zarray1(2,1:4)+zarray1(3,1:4)
!				zexex2(1,1:4)=(DCOS(r2)+zi*DSIN(r2))*zarray2(1,1:4)+ &
!				(DCOS(r2)-zi*DSIN(r2))*zarray2(2,1:4)+zarray2(3,1:4)
				! S_z=2
!				zexex1(2,1:4)=(DCOS(r2)+zi*DSIN(r2))*zarray1(4,1:4)+ &
!				(DCOS(r2)-zi*DSIN(r2))*zarray1(5,1:4)+zarray1(6,1:4)
!				zexex2(2,1:4)=(DCOS(r2)+zi*DSIN(r2))*zarray2(4,1:4)+ &
!				(DCOS(r2)-zi*DSIN(r2))*zarray2(5,1:4)+zarray2(6,1:4)
				! S_z=3
!				zexex1(3,1:4)=(DCOS(r2)+zi*DSIN(r2))*zarray1(7,1:4)+ &
!				(DCOS(r2)-zi*DSIN(r2))*zarray1(8,1:4)+zarray1(9,1:4)
!				zexex2(3,1:4)=(DCOS(r2)+zi*DSIN(r2))*zarray2(7,1:4)+ &
!				(DCOS(r2)-zi*DSIN(r2))*zarray2(8,1:4)+zarray2(9,1:4)
				! S_z average
!				zexex1(1,1:4)=(DCOS(r3)+zi*DSIN(r3))*zexex1(1,1:4)+ &
!				(DCOS(r3)-zi*DSIN(r3))*zexex1(2,1:4)+zexex1(3,1:4)
!				zexex2(1,1:4)=(DCOS(r3)+zi*DSIN(r3))*zexex2(1,1:4)+ &
!				(DCOS(r3)-zi*DSIN(r3))*zexex2(2,1:4)+zexex2(3,1:4)
!				z0(1:4)=(DCOS(rquark)-zi*DSIN(rquark))*zexex1(1,1:4)+ &
!				(DCOS(rquark)+zi*DSIN(rquark))*zexex2(1,1:4)
!			ENDIF
!		ENDIF
	  ENDIF
	  kk2=kk2+1
!        goto50
! 47      continue
    CASE(7)  ! scalar
      z0(0,1)=DCMPLX(dnou(1),dnou(0))
      z0(0,2)=DCMPLX(dnou(0),dnou(0))
      z0(0,3)=DCMPLX(dnou(0),dnou(0))
      z0(0,4)=DCMPLX(dnou(0),dnou(0))
	  kk2=kk2+1
	END SELECT
!        goto50
! 50     continue
    zy(0:2**Pwavenum-1,k,1)=z0(0:2**Pwavenum-1,1)
    zy(0:2**Pwavenum-1,k,2)=z0(0:2**Pwavenum-1,2)
    zy(0:2**Pwavenum-1,k,3)=z0(0:2**Pwavenum-1,3)
    zy(0:2**Pwavenum-1,k,4)=z0(0:2**Pwavenum-1,4)

    zyq(k,1)=DCMPLX(DREAL(zq(k,1)),dnou(0))*io(k)
    zyq(k,2)=DCMPLX(DREAL(zq(k,2)),dnou(0))*io(k)
    zyq(k,3)=zq(k,3)*io(k)
    zyq(k,4)=zq(k,4)*io(k)

! GAUGE INVARIANCE CHECK
    IF(k.EQ.n+1)THEN
		DO ijij=0,2**Pwavenum-1
			zy(ijij,k,1:4)=zyq(k,1:4)!/(zy(k,5)+zy(k,6))/2
		ENDDO
    ENDIF
ENDDO

END SUBROUTINE Helac_iniqq

! -----------------------------------------------------------------------      

SUBROUTINE Helac_nextq(kc)
IMPLICIT NONE
INTEGER,INTENT(IN)::kc
INTEGER::i,k,j,j0,j1,j2,j3,ir,ii,ii1,ii2,ii3,mm,m1,m2,m3,l,isgn,deriv  
COMPLEX(KIND=DBL)::zg,zgl,zgr,zmas
COMPLEX(KIND=DBL),DIMENSION(4)::zz01,zz02,zz00
INTEGER,DIMENSION(10)::ip1array,ip2array
DO i=n+1,imax
    DO k=1,4
 	  zy(0:2**Pwavenum-1,i,k)=DCMPLX(dnou(0),dnou(0))
	  zyq(i,k)=DCMPLX(dnou(0),dnou(0))
 	ENDDO
ENDDO

IF(kcou(kc).LT.1)RETURN
DO k=1,kcou(kc)
   ir=list(k,1,kc)           ! the number of routine
   ii=list(k,2,kc)           ! i0
   ii1=list(k,7,kc)          ! i1
   ii2=list(k,10,kc)         ! i2
   ii3=list(k,13,kc)         ! i3

   mm=list(k,3,kc)           ! m0 flavor
   m1=list(k,8,kc)           ! m1
   m2=list(k,11,kc)          ! m2
   m3=list(k,14,kc)          ! m3
   l= list(k,16,kc)          ! chiral of fermion left-1 right-2
   isgn=list(k,17,kc)        ! eps(i1,i2,i3)
   igluon=list(k,18,kc)      ! igluon
!   deriv=list(k,19,kc)       ! the rank of derivation
        
   j0=list(k,4,kc)           ! the label for i0,m0
   j1=list(k,9,kc)           ! the label for i1,m1
   j2=list(k,12,kc)          ! the label for i2,m2
   j3=0
   IF(list(k,13,kc).NE.0)j3=list(k,15,kc)  ! the label for i3,m3
   DO j=1,4
      z0(0:2**Pwavenum-1,j)=0
      z11(0:2**Pwavenum-1,j)=zy(0:2**Pwavenum-1,j1,j)   ! the external wave function of i1
      z21(0:2**Pwavenum-1,j)=zy(0:2**Pwavenum-1,j2,j)   ! the external wave function of i2
      zp0(j)=zyq(j1,j)+zyq(j2,j) ! the sum of momentum of external legs
      IF(list(k,13,kc).NE.0)THEN
          z31(0:2**Pwavenum-1,j)=zy(0:2**Pwavenum-1,j3,j)
          zp0(j)=zp0(j)+zyq(j3,j)
      ENDIF
      zyq(j0,j)=zp0(j)  ! the conservation of momentum
   ENDDO
   ip1array(1:10)=0
   ip2array(1:10)=0
!        goto(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)ir
   SELECT CASE(ir)        
! 4      continue
   CASE(4) !v3
     DO j=1,4
        zp1(j)=zyq(j1,j)
        zp2(j)=zyq(j2,j)
     ENDDO

     zg=zgv3(mm,m1,m2)
   ! the feynman rule in Phys.Rev.D 67,014026 (2003)
     IF(igluon.EQ.2)zg=-zg
     CALL Helac_xv3(z0(0,1:4),zp0,z11(0,1:4),zp1,z21(0,1:4),zp2,mm,zg,isgn)
	! It can be generalized to double or more p-waves
	 IF(Pwavenum.NE.0)THEN
		SELECT CASE(Pwavenum)
		CASE(1)
			CALL MomInclude(ii1,ip1array)
			CALL MomInclude(ii2,ip2array)
			CALL Helac_xv3d(zz00,zp0,z11(0,1:4),zp1,z21(0,1:4),zp2,mm,zg,isgn,1,ip1array,ip2array)
			CALL Helac_xv3(zz01,zp0,z11(1,1:4),zp1,z21(0,1:4),zp2,mm,zg,isgn)
			CALL Helac_xv3(zz02,zp0,z11(0,1:4),zp1,z21(1,1:4),zp2,mm,zg,isgn)
			z0(1,1:4)=zz00(1:4)+zz01(1:4)+zz02(1:4)
		CASE(2)
			CALL MomInclude(ii1,ip1array)
			CALL MomInclude(ii2,ip2array)
			CALL Helac_xv3d(zz00,zp0,z11(0,1:4),zp1,z21(0,1:4),zp2,mm,zg,isgn,1,ip1array,ip2array)
			CALL Helac_xv3(zz01,zp0,z11(1,1:4),zp1,z21(0,1:4),zp2,mm,zg,isgn)
			CALL Helac_xv3(zz02,zp0,z11(0,1:4),zp1,z21(1,1:4),zp2,mm,zg,isgn)
			z0(1,1:4)=zz00(1:4)+zz01(1:4)+zz02(1:4)
			CALL Helac_xv3d(zz00,zp0,z11(0,1:4),zp1,z21(0,1:4),zp2,mm,zg,isgn,2,ip1array,ip2array)
			CALL Helac_xv3(zz01,zp0,z11(2,1:4),zp1,z21(0,1:4),zp2,mm,zg,isgn)
			CALL Helac_xv3(zz02,zp0,z11(0,1:4),zp1,z21(2,1:4),zp2,mm,zg,isgn)
			z0(2,1:4)=zz00(1:4)+zz01(1:4)+zz02(1:4)
			z0(3,1:4)=0
			DO j=0,3
				CALL Helac_xv3(zz00,zp0,z11(3-j,1:4),zp1,z21(j,1:4),zp2,mm,zg,isgn)
				z0(3,1:4)=z0(3,1:4)+zz00(1:4)
			ENDDO
			DO j=1,2
				CALL Helac_xv3d(zz01,zp0,z11(j,1:4),zp1,z21(0,1:4),zp2,mm,zg,isgn,3-j,ip1array,ip2array)
				CALL Helac_xv3d(zz02,zp0,z11(0,1:4),zp1,z21(j,1:4),zp2,mm,zg,isgn,3-j,ip1array,ip2array)
				z0(3,1:4)=z0(3,1:4)+zz01(1:4)+zz02(1:4)
			ENDDO
		CASE DEFAULT
			PRINT *,"Wrong (v3) of Helac_nextq in Helac_pan1.f90! STOP !"
			STOP
		END SELECT
	 ENDIF
!        goto99
! 5      continue
   CASE(5)   ! v4     
     zg=zgv4(mm,m1,m2,m3)
     IF(igluon.LE.1)THEN
		CALL Helac_xv4(z0(0,1:4),zp0,z11(0,1:4),z21(0,1:4),z31(0,1:4),mm,zg,isgn)
		IF(Pwavenum.NE.0)THEN
			SELECT CASE(Pwavenum)
			CASE(1)
				CALL Helac_xv4(zz00(1:4),zp0,z11(0,1:4),z21(0,1:4),z31(1,1:4),mm,zg,isgn)
				CALL Helac_xv4(zz01(1:4),zp0,z11(1,1:4),z21(0,1:4),z31(0,1:4),mm,zg,isgn)
				CALL Helac_xv4(zz02(1:4),zp0,z11(0,1:4),z21(1,1:4),z31(0,1:4),mm,zg,isgn)
				z0(1,1:4)=zz00(1:4)+zz01(1:4)+zz02(1:4)
			CASE(2)
				CALL Helac_xv4(zz00(1:4),zp0,z11(0,1:4),z21(0,1:4),z31(1,1:4),mm,zg,isgn)
				CALL Helac_xv4(zz01(1:4),zp0,z11(1,1:4),z21(0,1:4),z31(0,1:4),mm,zg,isgn)
				CALL Helac_xv4(zz02(1:4),zp0,z11(0,1:4),z21(1,1:4),z31(0,1:4),mm,zg,isgn)
				z0(1,1:4)=zz00(1:4)+zz01(1:4)+zz02(1:4)
				CALL Helac_xv4(zz00(1:4),zp0,z11(0,1:4),z21(0,1:4),z31(2,1:4),mm,zg,isgn)
				CALL Helac_xv4(zz01(1:4),zp0,z11(2,1:4),z21(0,1:4),z31(0,1:4),mm,zg,isgn)
				CALL Helac_xv4(zz02(1:4),zp0,z11(0,1:4),z21(2,1:4),z31(0,1:4),mm,zg,isgn)
				z0(2,1:4)=zz00(1:4)+zz01(1:4)+zz02(1:4)
				z0(3,1:4)=0
				DO j=0,3
					DO i=0,3-j
						CALL Helac_xv4(zz00,zp0,z11(j,1:4),z21(i,1:4),z31(3-i-j,1:4),mm,zg,isgn)
						z0(3,1:4)=z0(3,1:4)+zz00(1:4)
					ENDDO
				ENDDO
			CASE DEFAULT
				PRINT *,"Wrong (v4,1) of Helac_nextq in Helac_pan1.f90! STOP !"
				STOP
			END SELECT
		ENDIF
	 ENDIF
     IF(igluon.EQ.2+2)THEN
		CALL Helac_xv4(z0(0,1:4),zp0,z11(0,1:4),z21(0,1:4),z31(0,1:4),mm,zg,isgn)
		IF(Pwavenum.NE.0)THEN
			SELECT CASE(Pwavenum)
			CASE(1)
				CALL Helac_xv4(zz00(1:4),zp0,z11(0,1:4),z21(0,1:4),z31(1,1:4),mm,zg,isgn)
				CALL Helac_xv4(zz01(1:4),zp0,z11(1,1:4),z21(0,1:4),z31(0,1:4),mm,zg,isgn)
				CALL Helac_xv4(zz02(1:4),zp0,z11(0,1:4),z21(1,1:4),z31(0,1:4),mm,zg,isgn)
				z0(1,1:4)=zz00(1:4)+zz01(1:4)+zz02(1:4)
			CASE(2)
				CALL Helac_xv4(zz00(1:4),zp0,z11(0,1:4),z21(0,1:4),z31(1,1:4),mm,zg,isgn)
				CALL Helac_xv4(zz01(1:4),zp0,z11(1,1:4),z21(0,1:4),z31(0,1:4),mm,zg,isgn)
				CALL Helac_xv4(zz02(1:4),zp0,z11(0,1:4),z21(1,1:4),z31(0,1:4),mm,zg,isgn)
				z0(1,1:4)=zz00(1:4)+zz01(1:4)+zz02(1:4)
				CALL Helac_xv4(zz00(1:4),zp0,z11(0,1:4),z21(0,1:4),z31(2,1:4),mm,zg,isgn)
				CALL Helac_xv4(zz01(1:4),zp0,z11(2,1:4),z21(0,1:4),z31(0,1:4),mm,zg,isgn)
				CALL Helac_xv4(zz02(1:4),zp0,z11(0,1:4),z21(2,1:4),z31(0,1:4),mm,zg,isgn)
				z0(2,1:4)=zz00(1:4)+zz01(1:4)+zz02(1:4)
				z0(3,1:4)=0
				DO j=0,3
					DO i=0,3-j
						CALL Helac_xv4(zz00,zp0,z11(j,1:4),z21(i,1:4),z31(3-i-j,1:4),mm,zg,isgn)
						z0(3,1:4)=z0(3,1:4)+zz00(1:4)
					ENDDO
				ENDDO				
			CASE DEFAULT
				PRINT *,"Wrong (v4,2+2) of Helac_nextq in Helac_pan1.f90! STOP !"
				STOP
			END SELECT
		ENDIF
	 ENDIF
     IF(igluon.EQ.3+2)THEN
		CALL Helac_xv4(z0(0,1:4),zp0,z21(0,1:4),z31(0,1:4),z11(0,1:4),mm,zg,isgn)
		IF(Pwavenum.NE.0)THEN
			SELECT CASE(Pwavenum)
			CASE(1)
				CALL Helac_xv4(zz00,zp0,z21(0,1:4),z31(1,1:4),z11(0,1:4),mm,zg,isgn)
				CALL Helac_xv4(zz01,zp0,z21(0,1:4),z31(0,1:4),z11(1,1:4),mm,zg,isgn)
				CALL Helac_xv4(zz02,zp0,z21(1,1:4),z31(0,1:4),z11(0,1:4),mm,zg,isgn)
				z0(1,1:4)=zz00(1:4)+zz01(1:4)+zz02(1:4)
			CASE(2)
				CALL Helac_xv4(zz00,zp0,z21(0,1:4),z31(1,1:4),z11(0,1:4),mm,zg,isgn)
				CALL Helac_xv4(zz01,zp0,z21(0,1:4),z31(0,1:4),z11(1,1:4),mm,zg,isgn)
				CALL Helac_xv4(zz02,zp0,z21(1,1:4),z31(0,1:4),z11(0,1:4),mm,zg,isgn)
				z0(1,1:4)=zz00(1:4)+zz01(1:4)+zz02(1:4)
				CALL Helac_xv4(zz00,zp0,z21(0,1:4),z31(2,1:4),z11(0,1:4),mm,zg,isgn)
				CALL Helac_xv4(zz01,zp0,z21(0,1:4),z31(0,1:4),z11(2,1:4),mm,zg,isgn)
				CALL Helac_xv4(zz02,zp0,z21(2,1:4),z31(0,1:4),z11(0,1:4),mm,zg,isgn)
				z0(2,1:4)=zz00(1:4)+zz01(1:4)+zz02(1:4)
				z0(3,1:4)=0
				DO j=0,3
					DO i=0,3-j
						CALL Helac_xv4(zz00,zp0,z21(j,1:4),z31(i,1:4),z11(3-i-j,1:4),mm,zg,isgn)
						z0(3,1:4)=z0(3,1:4)+zz00(1:4)
					ENDDO
				ENDDO
			CASE DEFAULT
				PRINT *,"Wrong (v4,3+2) of Helac_nextq in Helac_pan1.f90! STOP !"
				STOP
			END SELECT
		ENDIF
	 ENDIF
     IF(igluon.EQ.4+2)THEN
		CALL Helac_xv4(z0(0,1:4),zp0,z31(0,1:4),z11(0,1:4),z21(0,1:4),mm,zg,isgn)
		IF(Pwavenum.NE.0)THEN
			SELECT CASE(Pwavenum)
			CASE(1)
				CALL Helac_xv4(zz00,zp0,z31(1,1:4),z11(0,1:4),z21(0,1:4),mm,zg,isgn)
				CALL Helac_xv4(zz01,zp0,z31(0,1:4),z11(1,1:4),z21(0,1:4),mm,zg,isgn)
				CALL Helac_xv4(zz02,zp0,z31(0,1:4),z11(0,1:4),z21(1,1:4),mm,zg,isgn)
				z0(1,1:4)=zz00(1:4)+zz01(1:4)+zz02(1:4)
			CASE(2)
				CALL Helac_xv4(zz00,zp0,z31(1,1:4),z11(0,1:4),z21(0,1:4),mm,zg,isgn)
				CALL Helac_xv4(zz01,zp0,z31(0,1:4),z11(1,1:4),z21(0,1:4),mm,zg,isgn)
				CALL Helac_xv4(zz02,zp0,z31(0,1:4),z11(0,1:4),z21(1,1:4),mm,zg,isgn)
				z0(1,1:4)=zz00(1:4)+zz01(1:4)+zz02(1:4)
				CALL Helac_xv4(zz00,zp0,z31(2,1:4),z11(0,1:4),z21(0,1:4),mm,zg,isgn)
				CALL Helac_xv4(zz01,zp0,z31(0,1:4),z11(2,1:4),z21(0,1:4),mm,zg,isgn)
				CALL Helac_xv4(zz02,zp0,z31(0,1:4),z11(0,1:4),z21(2,1:4),mm,zg,isgn)
				z0(2,1:4)=zz00(1:4)+zz01(1:4)+zz02(1:4)
				z0(3,1:4)=0
				DO j=0,3
					DO i=0,3-j
						CALL Helac_xv4(zz00,zp0,z31(j,1:4),z11(i,1:4),z21(3-i-j,1:4),mm,zg,isgn)
						z0(3,1:4)=z0(3,1:4)+zz00(1:4)
					ENDDO
				ENDDO
			CASE DEFAULT
				PRINT *,"Wrong (v4,4+2) of Helac_nextq in Helac_pan1.f90! STOP !"
				STOP
			END SELECT
		ENDIF
	 ENDIF
!        goto99
! 1      continue
    CASE(1)  !vff
      DO j=1,4
        zp1(j)=zyq(j1,j)
        zp2(j)=zyq(j2,j)
      ENDDO
       
      zgl=zgvffl(mm,m1,m2)
      zgr=zgvffr(mm,m1,m2)
      CALL Helac_xvff(z0(0,1:4),zp0,z11(0,1:4),z21(0,1:4),mm,zgl,zgr,isgn)
	  IF(Pwavenum.NE.0)THEN
		SELECT CASE(Pwavenum)
		CASE(1)
			CALL Helac_xvff(zz01,zp0,z11(1,1:4),z21(0,1:4),mm,zgl,zgr,isgn)
			CALL Helac_xvff(zz02,zp0,z11(0,1:4),z21(1,1:4),mm,zgl,zgr,isgn)
			z0(1,1:4)=zz01(1:4)+zz02(1:4)
		CASE(2)
			CALL Helac_xvff(zz01,zp0,z11(1,1:4),z21(0,1:4),mm,zgl,zgr,isgn)
			CALL Helac_xvff(zz02,zp0,z11(0,1:4),z21(1,1:4),mm,zgl,zgr,isgn)
			z0(1,1:4)=zz01(1:4)+zz02(1:4)
			CALL Helac_xvff(zz01,zp0,z11(2,1:4),z21(0,1:4),mm,zgl,zgr,isgn)
			CALL Helac_xvff(zz02,zp0,z11(0,1:4),z21(2,1:4),mm,zgl,zgr,isgn)
			z0(2,1:4)=zz01(1:4)+zz02(1:4)
			z0(3,1:4)=0
			DO j=0,3
				CALL Helac_xvff(zz01,zp0,z11(j,1:4),z21(3-j,1:4),mm,zgl,zgr,isgn)
				z0(3,1:4)=z0(3,1:4)+zz01(1:4)
			ENDDO
		CASE DEFAULT
			PRINT *,"Wrong (vff) of Helac_nextq in Helac_pan1.f90! STOP !"
			STOP
		END SELECT
	  ENDIF
        
      IF(igluon.EQ.2)THEN
        DO i=1,4
           z0(0:2**(Pwavenum)-1,i)=z0(0:2**(Pwavenum)-1,i)/DSQRT(dnou(NCOL))&
		   *DCMPLX(dnou(0),dnou(1))
        ENDDO
      ENDIF 
      IF(igluon.EQ.3)THEN
        DO i=1,4
           z0(0:2**(Pwavenum)-1,i)=-z0(0:2**(Pwavenum)-1,i)/dnou(NCOL)
        ENDDO
      ENDIF
!        goto99
! 3      continue
    CASE(3)  ! ffv
       IF(ii.NE.2**n-2)THEN
          IF(l.EQ.2)THEN   ! right-handed of fermion
             zgr=zgvffr(an(m1),-mm,m2)
             zmas=SQRT(parmas(mm)**2-DCMPLX(0,1)*parmas(mm)*parwid(mm))
             CALL Helac_xffvr(z0(0,1:4),zp0,z11(0,1:4),z21(0,1:4),mm,zgr,isgn,zmas)
			 IF(Pwavenum.NE.0)THEN
				SELECT CASE(Pwavenum)
				CASE(1)
					CALL MomInclude(ii,ip1array)
					CALL Helac_xffvrd(zz00,zp0,z11(0,1:4),z21(0,1:4),mm,zgr,isgn,zmas,1,ip1array)
					CALL Helac_xffvr(zz01,zp0,z11(1,1:4),z21(0,1:4),mm,zgr,isgn,zmas)
					CALL Helac_xffvr(zz02,zp0,z11(0,1:4),z21(1,1:4),mm,zgr,isgn,zmas)
					z0(1,1:4)=zz00(1:4)+zz01(1:4)+zz02(1:4)
				CASE(2)
					CALL MomInclude(ii,ip1array)
					CALL Helac_xffvrd(zz00,zp0,z11(0,1:4),z21(0,1:4),mm,zgr,isgn,zmas,1,ip1array)
					CALL Helac_xffvr(zz01,zp0,z11(1,1:4),z21(0,1:4),mm,zgr,isgn,zmas)
					CALL Helac_xffvr(zz02,zp0,z11(0,1:4),z21(1,1:4),mm,zgr,isgn,zmas)
					z0(1,1:4)=zz00(1:4)+zz01(1:4)+zz02(1:4)
					CALL Helac_xffvrd(zz00,zp0,z11(0,1:4),z21(0,1:4),mm,zgr,isgn,zmas,2,ip1array)
					CALL Helac_xffvr(zz01,zp0,z11(2,1:4),z21(0,1:4),mm,zgr,isgn,zmas)
					CALL Helac_xffvr(zz02,zp0,z11(0,1:4),z21(2,1:4),mm,zgr,isgn,zmas)
					z0(2,1:4)=zz00(1:4)+zz01(1:4)+zz02(1:4)
					z0(3,1:4)=0
					DO j=0,3
						CALL Helac_xffvr(zz01,zp0,z11(j,1:4),z21(3-j,1:4),mm,zgr,isgn,zmas)
						z0(3,1:4)=z0(3,1:4)+zz01(1:4)
					ENDDO
					DO j=1,2
						CALL Helac_xffvrd(zz01,zp0,z11(j,1:4),z21(0,1:4),mm,zgr,isgn,zmas,3-j,ip1array)
						CALL Helac_xffvrd(zz02,zp0,z11(0,1:4),z21(j,1:4),mm,zgr,isgn,zmas,3-j,ip1array)
						z0(3,1:4)=z0(3,1:4)+zz01(1:4)+zz02(1:4)
					ENDDO
				CASE DEFAULT
					PRINT *,"Wrong (ffvr) of Helac_nextq in Helac_pan1.f90 ! STOP !"
					STOP
				END SELECT
			 ENDIF
          ELSEIF(l.EQ.1)THEN ! left-handed of fermion
             zgl=zgvffl(an(m1),-mm,m2)
             zmas=SQRT(parmas(mm)**2-DCMPLX(0,1)*parmas(mm)*parwid(mm))
             CALL Helac_xffvl(z0(0,1:4),zp0,z11(0,1:4),z21(0,1:4),mm,zgl,isgn,zmas)
			 IF(Pwavenum.NE.0)THEN
				SELECT CASE(Pwavenum)
				CASE(1)
					CALL MomInclude(ii,ip1array)
					CALL Helac_xffvld(zz00,zp0,z11(0,1:4),z21(0,1:4),mm,zgl,isgn,zmas,1,ip1array)
					CALL Helac_xffvl(zz01,zp0,z11(1,1:4),z21(0,1:4),mm,zgl,isgn,zmas)
					CALL Helac_xffvl(zz02,zp0,z11(0,1:4),z21(1,1:4),mm,zgl,isgn,zmas)
					z0(1,1:4)=zz00(1:4)+zz01(1:4)+zz02(1:4)
				CASE(2)
					CALL MomInclude(ii,ip1array)
					CALL Helac_xffvld(zz00,zp0,z11(0,1:4),z21(0,1:4),mm,zgl,isgn,zmas,1,ip1array)
					CALL Helac_xffvl(zz01,zp0,z11(1,1:4),z21(0,1:4),mm,zgl,isgn,zmas)
					CALL Helac_xffvl(zz02,zp0,z11(0,1:4),z21(1,1:4),mm,zgl,isgn,zmas)
					z0(1,1:4)=zz00(1:4)+zz01(1:4)+zz02(1:4)
					CALL Helac_xffvld(zz00,zp0,z11(0,1:4),z21(0,1:4),mm,zgl,isgn,zmas,2,ip1array)
					CALL Helac_xffvl(zz01,zp0,z11(2,1:4),z21(0,1:4),mm,zgl,isgn,zmas)
					CALL Helac_xffvl(zz02,zp0,z11(0,1:4),z21(2,1:4),mm,zgl,isgn,zmas)
					z0(2,1:4)=zz00(1:4)+zz01(1:4)+zz02(1:4)
					z0(3,1:4)=0
					DO j=0,3
						CALL Helac_xffvl(zz00,zp0,z11(j,1:4),z21(3-j,1:4),mm,zgl,isgn,zmas)
						z0(3,1:4)=z0(3,1:4)+zz00(1:4)
					ENDDO
					DO j=1,2
						CALL Helac_xffvld(zz01,zp0,z11(j,1:4),z21(0,1:4),mm,zgl,isgn,zmas,3-j,ip1array)
						CALL Helac_xffvld(zz02,zp0,z11(0,1:4),z21(j,1:4),mm,zgl,isgn,zmas,3-j,ip1array)
						z0(3,1:4)=z0(3,1:4)+zz01(1:4)+zz02(1:4)
					ENDDO
				CASE DEFAULT
					PRINT *,"Wrong (ffvl) of Helac_nextq in Helac_pan1.f90 ! STOP !"
					STOP
				END SELECT
			 ENDIF
          ENDIF
       ELSE   ! the last level
          IF(l.EQ.1)THEN
             zgl=zgvffl(an(m1),-mm,m2)
             CALL Helac_xffvl0(z0(0,1:4),z11(0,1:4),z21(0,1:4),zgl,isgn)
			 IF(Pwavenum.NE.0)THEN
				SELECT CASE(Pwavenum)
				CASE(1)
					CALL Helac_xffvl0(zz01,z11(1,1:4),z21(0,1:4),zgl,isgn)
					CALL Helac_xffvl0(zz02,z11(0,1:4),z21(1,1:4),zgl,isgn)
					z0(1,1:4)=zz01(1:4)+zz02(1:4)
				CASE(2)
					CALL Helac_xffvl0(zz01,z11(1,1:4),z21(0,1:4),zgl,isgn)
					CALL Helac_xffvl0(zz02,z11(0,1:4),z21(1,1:4),zgl,isgn)
					z0(1,1:4)=zz01(1:4)+zz02(1:4)
					CALL Helac_xffvl0(zz01,z11(2,1:4),z21(0,1:4),zgl,isgn)
					CALL Helac_xffvl0(zz02,z11(0,1:4),z21(2,1:4),zgl,isgn)
					z0(2,1:4)=zz01(1:4)+zz02(1:4)
					z0(3,1:4)=0
					DO j=0,3
						CALL Helac_xffvl0(zz00,z11(j,1:4),z21(3-j,1:4),zgl,isgn)
						z0(3,1:4)=z0(3,1:4)+zz00(1:4)
					ENDDO
				CASE DEFAULT
					PRINT *,"Wrong (vffl0) of Helac_nextq in Helac_pan1.f90 ! STOP !"
					STOP
				END SELECT
			 ENDIF
          ELSEIF(l.EQ.2)THEN
             zgr=zgvffr(an(m1),-mm,m2)
             CALL Helac_xffvr0(z0(0,1:4),z11(0,1:4),z21(0,1:4),zgr,isgn)
			 IF(Pwavenum.NE.0)THEN
				SELECT CASE(Pwavenum)
				CASE(1)
					CALL Helac_xffvr0(zz01,z11(1,1:4),z21(0,1:4),zgr,isgn)
					CALL Helac_xffvr0(zz02,z11(0,1:4),z21(1,1:4),zgr,isgn)
					z0(1,1:4)=zz01(1:4)+zz02(1:4)
				CASE(2)
					CALL Helac_xffvr0(zz01,z11(1,1:4),z21(0,1:4),zgr,isgn)
					CALL Helac_xffvr0(zz02,z11(0,1:4),z21(1,1:4),zgr,isgn)
					z0(1,1:4)=zz01(1:4)+zz02(1:4)
					CALL Helac_xffvr0(zz01,z11(2,1:4),z21(0,1:4),zgr,isgn)
					CALL Helac_xffvr0(zz02,z11(0,1:4),z21(2,1:4),zgr,isgn)
					z0(2,1:4)=zz01(1:4)+zz02(1:4)
					z0(3,1:4)=0
					DO j=0,3
						CALL Helac_xffvr0(zz00,z11(j,1:4),z21(3-j,1:4),zgr,isgn)
						z0(3,1:4)=z0(3,1:4)+zz00(1:4)
					ENDDO
				CASE DEFAULT
					PRINT *,"Wrong (vffr0) of Helac_nextq in Helac_pan1.f90 ! STOP !"
					STOP
				END SELECT
			 ENDIF
          ENDIF
       ENDIF

       IF(igluon.EQ.2)THEN
          DO i=1,4
             z0(0:2**Pwavenum-1,i)=z0(0:2**Pwavenum-1,i)/DSQRT(dnou(NCOL))&
			 *DCMPLX(dnou(0),dnou(1))
          ENDDO
       ENDIF
       IF(igluon.EQ.3)THEN
          DO i=1,4
             z0(0:2**Pwavenum-1,i)=-z0(0:2**Pwavenum-1,i)/dnou(NCOL)
          ENDDO
       ENDIF
!       goto99
! 2      continue
    CASE(2)  ! affv
       DO j=1,4
           zp1(j)=zyq(j1,j)
           zp2(j)=zyq(j2,j)
       ENDDO
       IF(ii.NE.2**n-2)THEN
           IF(l.EQ.2)THEN
               zgr=zgvffr(an(m1),m2,-mm)
               zmas=SQRT(parmas(mm)**2-DCMPLX(0,1)*parmas(mm)*parwid(mm))
               CALL Helac_xaffvr(z0(0,1:4),zp0,z11(0,1:4),z21(0,1:4),mm,zgr,isgn,zmas)
			   IF(Pwavenum.NE.0)THEN
					SELECT CASE(Pwavenum)
					CASE(1)
						CALL MomInclude(ii,ip1array)
						CALL Helac_xaffvrd(zz00,zp0,z11(0,1:4),z21(0,1:4),mm,zgr,isgn,zmas,1,ip1array)
						CALL Helac_xaffvr(zz01,zp0,z11(1,1:4),z21(0,1:4),mm,zgr,isgn,zmas)
						CALL Helac_xaffvr(zz02,zp0,z11(0,1:4),z21(1,1:4),mm,zgr,isgn,zmas)
						z0(1,1:4)=zz00(1:4)+zz01(1:4)+zz02(1:4)
					CASE(2)
						CALL MomInclude(ii,ip1array)
						CALL Helac_xaffvrd(zz00,zp0,z11(0,1:4),z21(0,1:4),mm,zgr,isgn,zmas,1,ip1array)
						CALL Helac_xaffvr(zz01,zp0,z11(1,1:4),z21(0,1:4),mm,zgr,isgn,zmas)
						CALL Helac_xaffvr(zz02,zp0,z11(0,1:4),z21(1,1:4),mm,zgr,isgn,zmas)
						z0(1,1:4)=zz00(1:4)+zz01(1:4)+zz02(1:4)
						CALL Helac_xaffvrd(zz00,zp0,z11(0,1:4),z21(0,1:4),mm,zgr,isgn,zmas,2,ip1array)
						CALL Helac_xaffvr(zz01,zp0,z11(2,1:4),z21(0,1:4),mm,zgr,isgn,zmas)
						CALL Helac_xaffvr(zz02,zp0,z11(0,1:4),z21(2,1:4),mm,zgr,isgn,zmas)
						z0(2,1:4)=zz00(1:4)+zz01(1:4)+zz02(1:4)
						z0(3,1:4)=0
						DO j=0,3
							CALL Helac_xaffvr(zz00,zp0,z11(j,1:4),z21(3-j,1:4),mm,zgr,isgn,zmas)
							z0(3,1:4)=z0(3,1:4)+zz00(1:4)
						ENDDO
						DO j=1,2
							CALL Helac_xaffvrd(zz01,zp0,z11(j,1:4),z21(0,1:4),mm,zgr,isgn,zmas,3-j,ip1array)
							CALL Helac_xaffvrd(zz02,zp0,z11(0,1:4),z21(j,1:4),mm,zgr,isgn,zmas,3-j,ip1array)
							z0(3,1:4)=z0(3,1:4)+zz01(1:4)+zz02(1:4)
						ENDDO
					CASE DEFAULT
						PRINT *,"Wrong (affvr) of Helac_nextq in Helac_pan1.f90 ! STOP !"
						STOP
					END SELECT
			   ENDIF
           ELSEIF(l.EQ.1)THEN
               zgl=zgvffl(an(m1),m2,-mm)
               zmas=SQRT(parmas(mm)**2-DCMPLX(0,1)*parmas(mm)*parwid(mm))
               CALL Helac_xaffvl(z0(0,1:4),zp0,z11(0,1:4),z21(0,1:4),mm,zgl,isgn,zmas)
			   IF(Pwavenum.NE.0)THEN
					SELECT CASE(Pwavenum)
					CASE(1)
						CALL MomInclude(ii,ip1array)
						CALL Helac_xaffvld(zz00,zp0,z11(0,1:4),z21(0,1:4),mm,zgl,isgn,zmas,1,ip1array)
						CALL Helac_xaffvl(zz01,zp0,z11(1,1:4),z21(0,1:4),mm,zgl,isgn,zmas)
						CALL Helac_xaffvl(zz02,zp0,z11(0,1:4),z21(1,1:4),mm,zgl,isgn,zmas)
						z0(1,1:4)=zz00(1:4)+zz01(1:4)+zz02(1:4)
					CASE(2)
						CALL MomInclude(ii,ip1array)
						CALL Helac_xaffvld(zz00,zp0,z11(0,1:4),z21(0,1:4),mm,zgl,isgn,zmas,1,ip1array)
						CALL Helac_xaffvl(zz01,zp0,z11(1,1:4),z21(0,1:4),mm,zgl,isgn,zmas)
						CALL Helac_xaffvl(zz02,zp0,z11(0,1:4),z21(1,1:4),mm,zgl,isgn,zmas)
						z0(1,1:4)=zz00(1:4)+zz01(1:4)+zz02(1:4)
						CALL Helac_xaffvld(zz00,zp0,z11(0,1:4),z21(0,1:4),mm,zgl,isgn,zmas,2,ip1array)
						CALL Helac_xaffvl(zz01,zp0,z11(2,1:4),z21(0,1:4),mm,zgl,isgn,zmas)
						CALL Helac_xaffvl(zz02,zp0,z11(0,1:4),z21(2,1:4),mm,zgl,isgn,zmas)
						z0(2,1:4)=zz00(1:4)+zz01(1:4)+zz02(1:4)
						z0(3,1:4)=0
						DO j=0,3
							CALL Helac_xaffvl(zz00,zp0,z11(j,1:4),z21(3-j,1:4),mm,zgl,isgn,zmas)
							z0(3,1:4)=z0(3,1:4)+zz00(1:4)
						ENDDO
						DO j=1,2
							CALL Helac_xaffvld(zz01,zp0,z11(j,1:4),z21(0,1:4),mm,zgl,isgn,zmas,3-j,ip1array)
							CALL Helac_xaffvld(zz02,zp0,z11(0,1:4),z21(j,1:4),mm,zgl,isgn,zmas,3-j,ip1array)
							z0(3,1:4)=z0(3,1:4)+zz01(1:4)+zz02(1:4)
						ENDDO
					CASE DEFAULT
						PRINT *,"Wrong (affvl) of Helac_nextq in Helac_pan1.f90 ! STOP !"
						STOP
					END SELECT
			   ENDIF
           ENDIF
       ELSEIF(ii.EQ.2**n-2)THEN
           IF(l.EQ.2)THEN
               zgr=zgvffr(an(m1),m2,-mm)
               CALL Helac_xaffvr0(z0(0,1:4),z11(0,1:4),z21(0,1:4),zgr,isgn)
			   IF(Pwavenum.NE.0)THEN
					SELECT CASE(Pwavenum)
					CASE(1)
						CALL Helac_xaffvr0(zz01,z11(1,1:4),z21(0,1:4),zgr,isgn)
						CALL Helac_xaffvr0(zz02,z11(0,1:4),z21(1,1:4),zgr,isgn)
						z0(1,1:4)=zz01(1:4)+zz02(1:4)
					CASE(2)
						CALL Helac_xaffvr0(zz01,z11(1,1:4),z21(0,1:4),zgr,isgn)
						CALL Helac_xaffvr0(zz02,z11(0,1:4),z21(1,1:4),zgr,isgn)
						z0(1,1:4)=zz01(1:4)+zz02(1:4)
						CALL Helac_xaffvr0(zz01,z11(2,1:4),z21(0,1:4),zgr,isgn)
						CALL Helac_xaffvr0(zz02,z11(0,1:4),z21(2,1:4),zgr,isgn)
						z0(2,1:4)=zz01(1:4)+zz02(1:4)
						z0(3,1:4)=0
						DO j=0,3
							CALL Helac_xaffvr0(zz00,z11(j,1:4),z21(3-j,1:4),zgr,isgn)
							z0(3,1:4)=z0(3,1:4)+zz00(1:4)
						ENDDO
					CASE DEFAULT
						PRINT *,"Wrong (affvr0) of Helac_nextq in Helac_pan1.f90 ! STOP !"
						STOP
					END SELECT
			   ENDIF
           ELSEIF(l.EQ.1)THEN
               zgl=zgvffl(an(m1),m2,-mm)
               CALL Helac_xaffvl0(z0(0,1:4),z11(0,1:4),z21(0,1:4),zgl,isgn)
			   IF(Pwavenum.NE.0)THEN
					SELECT CASE(Pwavenum)
					CASE(1)
						CALL Helac_xaffvl0(zz01,z11(1,1:4),z21(0,1:4),zgl,isgn)
						CALL Helac_xaffvl0(zz02,z11(0,1:4),z21(1,1:4),zgl,isgn)
						z0(1,1:4)=zz01(1:4)+zz02(1:4)
					CASE(2)
						CALL Helac_xaffvl0(zz01,z11(1,1:4),z21(0,1:4),zgl,isgn)
						CALL Helac_xaffvl0(zz02,z11(0,1:4),z21(1,1:4),zgl,isgn)
						z0(1,1:4)=zz01(1:4)+zz02(1:4)
						CALL Helac_xaffvl0(zz01,z11(2,1:4),z21(0,1:4),zgl,isgn)
						CALL Helac_xaffvl0(zz02,z11(0,1:4),z21(2,1:4),zgl,isgn)
						z0(2,1:4)=zz01(1:4)+zz02(1:4)
						z0(3,1:4)=0
						DO j=0,3
							CALL Helac_xaffvl0(zz00,z11(j,1:4),z21(3-j,1:4),zgl,isgn)
							z0(3,1:4)=z0(3,1:4)+zz00(1:4)
						ENDDO
					CASE DEFAULT
						PRINT *,"Wrong (affvl0) of Helac_nextq in Helac_pan1.f90 ! STOP !"
						STOP
					END SELECT
			   ENDIF
           ENDIF
       ENDIF
 
       IF(igluon.EQ.2)THEN
           DO i=1,4
               z0(0:2**Pwavenum-1,i)=z0(0:2**Pwavenum-1,i)/DSQRT(dnou(NCOL))&
			   *DCMPLX(dnou(0),dnou(1))
           ENDDO
       ENDIF
       IF(igluon.EQ.3)THEN
           DO i=1,4
               z0(0:2**Pwavenum-1,i)=-z0(0:2**Pwavenum-1,i)/dnou(NCOL)
           ENDDO
       ENDIF
!  goto99
! 6      continue
    CASE(6)   ! ffs
       IF(ii.NE.2**n-2)THEN
           IF(l.EQ.2)THEN
               zgr=zgsffr(an(m1),-mm,m2)
               CALL Helac_xffsr(z0(0,1:4),zp0,z11(0,1),z21(0,1:4),mm,zgr,isgn,DCMPLX(parmas(mm)))
			   IF(Pwavenum.NE.0)THEN
					SELECT CASE(Pwavenum)
					CASE(1)
						CALL MomInclude(ii,ip1array)
						CALL Helac_xffsrd(zz00,zp0,z11(0,1),z21(0,1:4),mm,zgr,isgn,DCMPLX(parmas(mm)),1,ip1array)
						CALL Helac_xffsr(zz01,zp0,z11(1,1),z21(0,1:4),mm,zgr,isgn,DCMPLX(parmas(mm)))
						CALL Helac_xffsr(zz02,zp0,z11(0,1),z21(1,1:4),mm,zgr,isgn,DCMPLX(parmas(mm)))
						z0(1,1:4)=zz00(1:4)+zz01(1:4)+zz02(1:4)
					CASE(2)
						CALL MomInclude(ii,ip1array)
						CALL Helac_xffsrd(zz00,zp0,z11(0,1),z21(0,1:4),mm,zgr,isgn,DCMPLX(parmas(mm)),1,ip1array)
						CALL Helac_xffsr(zz01,zp0,z11(1,1),z21(0,1:4),mm,zgr,isgn,DCMPLX(parmas(mm)))
						CALL Helac_xffsr(zz02,zp0,z11(0,1),z21(1,1:4),mm,zgr,isgn,DCMPLX(parmas(mm)))
						z0(1,1:4)=zz00(1:4)+zz01(1:4)+zz02(1:4)
						CALL Helac_xffsrd(zz00,zp0,z11(0,1),z21(0,1:4),mm,zgr,isgn,DCMPLX(parmas(mm)),2,ip1array)
						CALL Helac_xffsr(zz01,zp0,z11(2,1),z21(0,1:4),mm,zgr,isgn,DCMPLX(parmas(mm)))
						CALL Helac_xffsr(zz02,zp0,z11(0,1),z21(2,1:4),mm,zgr,isgn,DCMPLX(parmas(mm)))
						z0(2,1:4)=zz00(1:4)+zz01(1:4)+zz02(1:4)
						z0(3,1:4)=0
						DO j=0,3
							CALL Helac_xffsr(zz00,zp0,z11(j,1),z21(3-j,1:4),mm,zgr,isgn,DCMPLX(parmas(mm)))
							z0(3,1:4)=z0(3,1:4)+zz00(1:4)
						ENDDO
						DO j=1,2
							CALL Helac_xffsrd(zz01,zp0,z11(j,1),z21(0,1:4),mm,zgr,isgn,DCMPLX(parmas(mm)),3-j,ip1array)
							CALL Helac_xffsrd(zz02,zp0,z11(0,1),z21(j,1:4),mm,zgr,isgn,DCMPLX(parmas(mm)),3-j,ip1array)
							z0(3,1:4)=z0(3,1:4)+zz01(1:4)+zz02(1:4)
						ENDDO
					CASE DEFAULT
						PRINT *,"Wrong (ffsr) of Helac_nextq in Helac_pan1.f90 ! STOP !"
						STOP
					END SELECT
			   ENDIF
           ELSEIF(l.EQ.1)THEN
               zgl=zgsffl(an(m1),-mm,m2)
               CALL Helac_xffsl(z0(0,1:4),zp0,z11(0,1),z21(0,1:4),mm,zgl,isgn,DCMPLX(parmas(mm)))
			   IF(Pwavenum.NE.0)THEN
					SELECT CASE(Pwavenum)
					CASE(1)
						CALL MomInclude(ii,ip1array)
						CALL Helac_xffsld(zz00,zp0,z11(0,1),z21(0,1:4),mm,zgl,isgn,DCMPLX(parmas(mm)),1,ip1array)
						CALL Helac_xffsl(zz01,zp0,z11(1,1),z21(0,1:4),mm,zgl,isgn,DCMPLX(parmas(mm)))
						CALL Helac_xffsl(zz02,zp0,z11(0,1),z21(1,1:4),mm,zgl,isgn,DCMPLX(parmas(mm)))
						z0(1,1:4)=zz00(1:4)+zz01(1:4)+zz02(1:4)
					CASE(2)
						CALL MomInclude(ii,ip1array)
						CALL Helac_xffsld(zz00,zp0,z11(0,1),z21(0,1:4),mm,zgl,isgn,DCMPLX(parmas(mm)),1,ip1array)
						CALL Helac_xffsl(zz01,zp0,z11(1,1),z21(0,1:4),mm,zgl,isgn,DCMPLX(parmas(mm)))
						CALL Helac_xffsl(zz02,zp0,z11(0,1),z21(1,1:4),mm,zgl,isgn,DCMPLX(parmas(mm)))
						z0(1,1:4)=zz00(1:4)+zz01(1:4)+zz02(1:4)
						CALL Helac_xffsld(zz00,zp0,z11(0,1),z21(0,1:4),mm,zgl,isgn,DCMPLX(parmas(mm)),2,ip1array)
						CALL Helac_xffsl(zz01,zp0,z11(2,1),z21(0,1:4),mm,zgl,isgn,DCMPLX(parmas(mm)))
						CALL Helac_xffsl(zz02,zp0,z11(0,1),z21(2,1:4),mm,zgl,isgn,DCMPLX(parmas(mm)))
						z0(2,1:4)=zz00(1:4)+zz01(1:4)+zz02(1:4)
						z0(3,1:4)=0
						DO j=0,3
							CALL Helac_xffsl(zz00,zp0,z11(j,1),z21(3-j,1:4),mm,zgl,isgn,DCMPLX(parmas(mm)))
							z0(3,1:4)=z0(3,1:4)+zz00(1:4)
						ENDDO
						DO j=1,2
							CALL Helac_xffsld(zz01,zp0,z11(j,1),z21(0,1:4),mm,zgl,isgn,DCMPLX(parmas(mm)),3-j,ip1array)
							CALL Helac_xffsld(zz02,zp0,z11(0,1),z21(j,1:4),mm,zgl,isgn,DCMPLX(parmas(mm)),3-j,ip1array)
							z0(3,1:4)=z0(3,1:4)+zz01(1:4)+zz02(1:4)
						ENDDO
					CASE DEFAULT
						PRINT *,"Wrong (ffsl) of Helac_nextq in Helac_pan1.f90 ! STOP !"
						STOP
					END SELECT
			   ENDIF
           ENDIF
        ELSE
           IF(l.EQ.1)THEN
               zgl=zgsffl(an(m1),-mm,m2)
               CALL Helac_xffsl0(z0(0,1:4),z11(0,1),z21(0,1:4),zgl,isgn)
			   IF(Pwavenum.NE.0)THEN
					SELECT CASE(Pwavenum)
					CASE(1)
						CALL Helac_xffsl0(zz01,z11(1,1),z21(0,1:4),zgl,isgn)
						CALL Helac_xffsl0(zz02,z11(0,1),z21(1,1:4),zgl,isgn)
						z0(1,1:4)=zz01(1:4)+zz02(1:4)
					CASE(2)
						CALL Helac_xffsl0(zz01,z11(1,1),z21(0,1:4),zgl,isgn)
						CALL Helac_xffsl0(zz02,z11(0,1),z21(1,1:4),zgl,isgn)
						z0(1,1:4)=zz01(1:4)+zz02(1:4)
						CALL Helac_xffsl0(zz01,z11(2,1),z21(0,1:4),zgl,isgn)
						CALL Helac_xffsl0(zz02,z11(0,1),z21(2,1:4),zgl,isgn)
						z0(2,1:4)=zz01(1:4)+zz02(1:4)
						z0(3,1:4)=0
						DO j=0,3
							CALL Helac_xffsl0(zz00,z11(j,1),z21(3-j,1:4),zgl,isgn)
							z0(3,1:4)=z0(3,1:4)+zz00(1:4)
						ENDDO
					CASE DEFAULT
						PRINT *,"Wrong (ffsl0) of Helac_nextq in Helac_pan1.f90 ! STOP !"
						STOP
					END SELECT
			   ENDIF
           ELSEIF(l.EQ.2)THEN
               zgr=zgsffr(an(m1),-mm,m2)
               CALL Helac_xffsr0(z0(0,1:4),z11(0,1),z21(0,1:4),zgr,isgn)
			   IF(Pwavenum.NE.0)THEN
					SELECT CASE(Pwavenum)
					CASE(1)
						CALL Helac_xffsr0(zz01,z11(1,1),z21(0,1:4),zgr,isgn)
						CALL Helac_xffsr0(zz02,z11(0,1),z21(1,1:4),zgr,isgn)
						z0(1,1:4)=zz01(1:4)+zz02(1:4)
					CASE(2)
						CALL Helac_xffsr0(zz01,z11(1,1),z21(0,1:4),zgr,isgn)
						CALL Helac_xffsr0(zz02,z11(0,1),z21(1,1:4),zgr,isgn)
						z0(1,1:4)=zz01(1:4)+zz02(1:4)
						CALL Helac_xffsr0(zz01,z11(2,1),z21(0,1:4),zgr,isgn)
						CALL Helac_xffsr0(zz02,z11(0,1),z21(2,1:4),zgr,isgn)
						z0(2,1:4)=zz01(1:4)+zz02(1:4)
						z0(3,1:4)=0
						DO j=0,3
							CALL Helac_xffsr0(zz00,z11(1,1),z21(0,1:4),zgr,isgn)
							z0(3,1:4)=z0(3,1:4)+zz00(1:4)
						ENDDO
					CASE DEFAULT
						PRINT *,"Wrong (sffr0) of Helac_nextq in Helac_pan1.f90 ! STOP !"
						STOP
					END SELECT
			   ENDIF
           ENDIF
        ENDIF
!  goto99
! 7      continue
    CASE(7)  ! affs
        IF(ii.NE.2**n-2)THEN
           IF(l.EQ.1)THEN
               zgl=zgsffl(an(m1),m2,-mm)
               CALL Helac_xaffsl(z0(0,1:4),zp0,z11(0,1),z21(0,1:4),mm,zgl,isgn,DCMPLX(parmas(mm)))
			   IF(Pwavenum.NE.0)THEN
					SELECT CASE(Pwavenum)
					CASE(1)
						CALL MomInclude(ii,ip1array)
						CALL Helac_xaffsld(zz00,zp0,z11(0,1),z21(0,1:4),mm,zgl,isgn,DCMPLX(parmas(mm)),1,ip1array)
						CALL Helac_xaffsl(zz01,zp0,z11(1,1),z21(0,1:4),mm,zgl,isgn,DCMPLX(parmas(mm)))
						CALL Helac_xaffsl(zz02,zp0,z11(0,1),z21(1,1:4),mm,zgl,isgn,DCMPLX(parmas(mm)))
						z0(1,1:4)=zz00(1:4)+zz01(1:4)+zz02(1:4)
					CASE(2)
						CALL MomInclude(ii,ip1array)
						CALL Helac_xaffsld(zz00,zp0,z11(0,1),z21(0,1:4),mm,zgl,isgn,DCMPLX(parmas(mm)),1,ip1array)
						CALL Helac_xaffsl(zz01,zp0,z11(1,1),z21(0,1:4),mm,zgl,isgn,DCMPLX(parmas(mm)))
						CALL Helac_xaffsl(zz02,zp0,z11(0,1),z21(1,1:4),mm,zgl,isgn,DCMPLX(parmas(mm)))
						z0(1,1:4)=zz00(1:4)+zz01(1:4)+zz02(1:4)
						CALL Helac_xaffsld(zz00,zp0,z11(0,1),z21(0,1:4),mm,zgl,isgn,DCMPLX(parmas(mm)),2,ip1array)
						CALL Helac_xaffsl(zz01,zp0,z11(2,1),z21(0,1:4),mm,zgl,isgn,DCMPLX(parmas(mm)))
						CALL Helac_xaffsl(zz02,zp0,z11(0,1),z21(2,1:4),mm,zgl,isgn,DCMPLX(parmas(mm)))
						z0(2,1:4)=zz00(1:4)+zz01(1:4)+zz02(1:4)
						z0(3,1:4)=0
						DO j=0,3
							CALL Helac_xaffsl(zz00,zp0,z11(j,1),z21(3-j,1:4),mm,zgl,isgn,DCMPLX(parmas(mm)))
							z0(3,1:4)=z0(3,1:4)+zz00(1:4)
						ENDDO
						DO j=1,2
							CALL Helac_xaffsld(zz01,zp0,z11(j,1),z21(0,1:4),mm,zgl,isgn,DCMPLX(parmas(mm)),3-j,ip1array)
							CALL Helac_xaffsld(zz02,zp0,z11(0,1),z21(j,1:4),mm,zgl,isgn,DCMPLX(parmas(mm)),3-j,ip1array)
							z0(3,1:4)=z0(3,1:4)+zz01(1:4)+zz02(1:4)
						ENDDO
					CASE DEFAULT
						PRINT *,"Wrong (affsl) of Helac_nextq in Helac_pan1.f90 ! STOP !"
						STOP
					END SELECT
			   ENDIF
           ELSEIF(l.EQ.2)THEN
               zgr=zgsffr(an(m1),m2,-mm)
               CALL Helac_xaffsr(z0(0,1:4),zp0,z11(0,1),z21(0,1:4),mm,zgr,isgn,DCMPLX(parmas(mm)))
			   IF(Pwavenum.NE.0)THEN
					SELECT CASE(Pwavenum)
					CASE(1)
						CALL MomInclude(ii,ip1array)
						CALL Helac_xaffsrd(zz00,zp0,z11(0,1),z21(0,1:4),mm,zgr,isgn,DCMPLX(parmas(mm)),1,ip1array)
						CALL Helac_xaffsr(zz01,zp0,z11(1,1),z21(0,1:4),mm,zgr,isgn,DCMPLX(parmas(mm)))
						CALL Helac_xaffsr(zz02,zp0,z11(0,1),z21(1,1:4),mm,zgr,isgn,DCMPLX(parmas(mm)))
						z0(1,1:4)=zz00(1:4)+zz01(1:4)+zz02(1:4)
					CASE(2)
						CALL MomInclude(ii,ip1array)
						CALL Helac_xaffsrd(zz00,zp0,z11(0,1),z21(0,1:4),mm,zgr,isgn,DCMPLX(parmas(mm)),1,ip1array)
						CALL Helac_xaffsr(zz01,zp0,z11(1,1),z21(0,1:4),mm,zgr,isgn,DCMPLX(parmas(mm)))
						CALL Helac_xaffsr(zz02,zp0,z11(0,1),z21(1,1:4),mm,zgr,isgn,DCMPLX(parmas(mm)))
						z0(1,1:4)=zz00(1:4)+zz01(1:4)+zz02(1:4)
						CALL Helac_xaffsrd(zz00,zp0,z11(0,1),z21(0,1:4),mm,zgr,isgn,DCMPLX(parmas(mm)),2,ip1array)
						CALL Helac_xaffsr(zz01,zp0,z11(2,1),z21(0,1:4),mm,zgr,isgn,DCMPLX(parmas(mm)))
						CALL Helac_xaffsr(zz02,zp0,z11(0,1),z21(2,1:4),mm,zgr,isgn,DCMPLX(parmas(mm)))
						z0(2,1:4)=zz00(1:4)+zz01(1:4)+zz02(1:4)
						z0(3,1:4)=0
						DO j=0,3
							CALL Helac_xaffsr(zz00,zp0,z11(j,1),z21(3-j,1:4),mm,zgr,isgn,DCMPLX(parmas(mm)))
							z0(3,1:4)=z0(3,1:4)+zz00(1:4)
						ENDDO
						DO j=1,2
							CALL Helac_xaffsrd(zz01,zp0,z11(j,1),z21(0,1:4),mm,zgr,isgn,DCMPLX(parmas(mm)),3-j,ip1array)
							CALL Helac_xaffsrd(zz02,zp0,z11(0,1),z21(j,1:4),mm,zgr,isgn,DCMPLX(parmas(mm)),3-j,ip1array)
							z0(3,1:4)=z0(3,1:4)+zz01(1:4)+zz02(1:4)
						ENDDO
					CASE DEFAULT
						PRINT *,"Wrong (affsr) of Helac_nextq in Helac_pan1.f90 ! STOP"
						STOP
					END SELECT
			   ENDIF
           ENDIF  
        ELSE
           IF(l.EQ.2)THEN
               zgr=zgsffr(an(m1),m2,-mm)
               CALL Helac_xaffsr0(z0(0,1:4),z11(0,1),z21(0,1:4),zgr,isgn)
			   IF(Pwavenum.NE.0)THEN
					SELECT CASE(Pwavenum)
					CASE(1)
						CALL Helac_xaffsr0(zz01,z11(1,1),z21(0,1:4),zgr,isgn)
						CALL Helac_xaffsr0(zz02,z11(0,1),z21(1,1:4),zgr,isgn)
						z0(1,1:4)=zz01(1:4)+zz02(1:4)
					CASE(2)
						CALL Helac_xaffsr0(zz01,z11(1,1),z21(0,1:4),zgr,isgn)
						CALL Helac_xaffsr0(zz02,z11(0,1),z21(1,1:4),zgr,isgn)
						z0(1,1:4)=zz01(1:4)+zz02(1:4)
						CALL Helac_xaffsr0(zz01,z11(2,1),z21(0,1:4),zgr,isgn)
						CALL Helac_xaffsr0(zz02,z11(0,1),z21(2,1:4),zgr,isgn)
						z0(2,1:4)=zz01(1:4)+zz02(1:4)
						z0(3,1:4)=0
						DO j=0,3
							CALL Helac_xaffsr0(zz00,z11(j,1),z21(3-j,1:4),zgr,isgn)
							z0(3,1:4)=z0(3,1:4)+zz00(1:4)
						ENDDO
					CASE DEFAULT
						PRINT *,"Wrong (affsr0) of Helac_nextq in Helac_pan1.f90 ! STOP !"
						STOP
					END SELECT
			   ENDIF
           ELSEIF(l.EQ.1)THEN
               zgl=zgsffl(an(m1),m2,-mm)
               CALL Helac_xaffsl0(z0(0,1:4),z11(0,1),z21(0,1:4),zgl,isgn)
			   IF(Pwavenum.NE.0)THEN
					SELECT CASE(Pwavenum)
					CASE(1)
						CALL Helac_xaffsl0(zz01,z11(1,1),z21(0,1:4),zgl,isgn)
						CALL Helac_xaffsl0(zz02,z11(0,1),z21(1,1:4),zgl,isgn)
						z0(1,1:4)=zz01(1:4)+zz02(1:4)
					CASE(2)
						CALL Helac_xaffsl0(zz01,z11(1,1),z21(0,1:4),zgl,isgn)
						CALL Helac_xaffsl0(zz02,z11(0,1),z21(1,1:4),zgl,isgn)
						z0(1,1:4)=zz01(1:4)+zz02(1:4)
						CALL Helac_xaffsl0(zz01,z11(2,1),z21(0,1:4),zgl,isgn)
						CALL Helac_xaffsl0(zz02,z11(0,1),z21(2,1:4),zgl,isgn)
						z0(2,1:4)=zz01(1:4)+zz02(1:4)
						z0(3,1:4)=0
						DO j=0,3
							CALL Helac_xaffsl0(zz00,z11(j,1),z21(3-j,1:4),zgl,isgn)
							z0(3,1:4)=z0(3,1:4)+zz00(1:4)
						ENDDO
					CASE DEFAULT
						PRINT *,"Wrong (affsl0) of Helac_nextq in Helac_pan1.f90 ! STOP !"
						STOP
					END SELECT
			   ENDIF
           ENDIF
        ENDIF
!        goto99
! 8      continue
    CASE(8)  ! sff      
        zgl=zgsffl(mm,m1,m2)
        zgr=zgsffr(mm,m1,m2) 
        CALL Helac_xsff(z0(0,1),zp0,z11(0,1:4),z21(0,1:4),mm,zgl,zgr,isgn)
		IF(Pwavenum.NE.0)THEN
			SELECT CASE(Pwavenum)
			CASE(1)
				CALL Helac_xsff(zz01(1),zp0,z11(1,1:4),z21(0,1:4),mm,zgl,zgr,isgn)
				CALL Helac_xsff(zz02(1),zp0,z11(0,1:4),z21(1,1:4),mm,zgl,zgr,isgn)
				z0(1,1)=zz01(1)+zz02(1)
			CASE(2)
				CALL Helac_xsff(zz01(1),zp0,z11(1,1:4),z21(0,1:4),mm,zgl,zgr,isgn)
				CALL Helac_xsff(zz02(1),zp0,z11(0,1:4),z21(1,1:4),mm,zgl,zgr,isgn)
				z0(1,1)=zz01(1)+zz02(1)
				CALL Helac_xsff(zz01(1),zp0,z11(2,1:4),z21(0,1:4),mm,zgl,zgr,isgn)
				CALL Helac_xsff(zz02(1),zp0,z11(0,1:4),z21(2,1:4),mm,zgl,zgr,isgn)
				z0(2,1)=zz01(1)+zz02(1)
				z0(3,1:4)=0
				DO j=0,3
					CALL Helac_xsff(zz00(1),zp0,z11(j,1:4),z21(3-j,1:4),mm,zgl,zgr,isgn)
					z0(3,1)=z0(3,1)+zz00(1)
				ENDDO
			CASE DEFAULT
				PRINT *,"Wrong (sff) of Helac_nextq in Helac_pan1.f90 ! STOP !"
				STOP
			END SELECT
		ENDIF
!        goto99
! 9      continue
    CASE(9)  !vvs       
        zg=zgvvs(mm,m1,m2)
        CALL Helac_xvvs(z0(0,1:4),zp0,z11(0,1:4),z21(0,1),mm,zg,isgn)
		IF(Pwavenum.NE.0)THEN
			SELECT CASE(Pwavenum)
			CASE(1)
				CALL Helac_xvvs(zz01,zp0,z11(1,1:4),z21(0,1),mm,zg,isgn)
				CALL Helac_xvvs(zz02,zp0,z11(0,1:4),z21(1,1),mm,zg,isgn)
				z0(1,1:4)=zz01(1:4)+zz02(1:4)
			CASE(2)
				CALL Helac_xvvs(zz01,zp0,z11(1,1:4),z21(0,1),mm,zg,isgn)
				CALL Helac_xvvs(zz02,zp0,z11(0,1:4),z21(1,1),mm,zg,isgn)
				z0(1,1:4)=zz01(1:4)+zz02(1:4)
				CALL Helac_xvvs(zz01,zp0,z11(2,1:4),z21(0,1),mm,zg,isgn)
				CALL Helac_xvvs(zz02,zp0,z11(0,1:4),z21(2,1),mm,zg,isgn)
				z0(2,1:4)=zz01(1:4)+zz02(1:4)
				z0(3,1:4)=0
				DO j=0,3
					CALL Helac_xvvs(zz00,zp0,z11(j,1:4),z21(3-j,1),mm,zg,isgn)
					z0(3,1:4)=z0(3,1:4)+zz00(1:4)
				ENDDO
			CASE DEFAULT
				PRINT *,"Wrong (vvs) of Helac_nextq in Helac_pan1.f90 ! STOP !"
				STOP
			END SELECT
		ENDIF
!       goto99
! 10     continue
    CASE(10)  ! svv
        zg=zgsvv(mm,m1,m2)
        CALL Helac_xsvv(z0(0,1),zp0,z11(0,1:4),z21(0,1:4),mm,zg,isgn)
		IF(Pwavenum.NE.0)THEN
			SELECT CASE(Pwavenum)
			CASE(1)
				CALL Helac_xsvv(zz01(1),zp0,z11(1,1:4),z21(0,1:4),mm,zg,isgn)
				CALL Helac_xsvv(zz02(1),zp0,z11(0,1:4),z21(1,1:4),mm,zg,isgn)
				z0(1,1)=zz01(1)+zz02(1)
			CASE(2)
				CALL Helac_xsvv(zz01(1),zp0,z11(1,1:4),z21(0,1:4),mm,zg,isgn)
				CALL Helac_xsvv(zz02(1),zp0,z11(0,1:4),z21(1,1:4),mm,zg,isgn)
				z0(1,1)=zz01(1)+zz02(1)
				CALL Helac_xsvv(zz01(1),zp0,z11(2,1:4),z21(0,1:4),mm,zg,isgn)
				CALL Helac_xsvv(zz02(1),zp0,z11(0,1:4),z21(2,1:4),mm,zg,isgn)
				z0(2,1)=zz01(1)+zz02(1)
				z0(3,1:4)=0
				DO j=0,3
					CALL Helac_xsvv(zz00(1),zp0,z11(j,1:4),z21(3-j,1:4),mm,zg,isgn)
					z0(3,1)=z0(3,1)+zz00(1)
				ENDDO
			CASE DEFAULT
				PRINT *,"Wrong (svv) of Helac_nextq in Helac_pan1.f90 ! STOP !"
				STOP
			END SELECT
		ENDIF
!        goto99
! 11     continue
    CASE(11) ! vvss      
        zg=zgvvss(mm,m1,m2,m3)
        CALL Helac_xvvss(z0(0,1:4),zp0,z11(0,1:4),z21(0,1),z31(0,1),mm,zg,isgn)
		IF(Pwavenum.NE.0)THEN
			SELECT CASE(Pwavenum)
			CASE(1)
				CALL Helac_xvvss(zz00,zp0,z11(1,1:4),z21(0,1),z31(0,1),mm,zg,isgn)
				CALL Helac_xvvss(zz01,zp0,z11(0,1:4),z21(1,1),z31(0,1),mm,zg,isgn)
				CALL Helac_xvvss(zz02,zp0,z11(0,1:4),z21(0,1),z31(1,1),mm,zg,isgn)
				z0(1,1:4)=zz00(1:4)+zz01(1:4)+zz02(1:4)
			CASE(2)
				CALL Helac_xvvss(zz00,zp0,z11(1,1:4),z21(0,1),z31(0,1),mm,zg,isgn)
				CALL Helac_xvvss(zz01,zp0,z11(0,1:4),z21(1,1),z31(0,1),mm,zg,isgn)
				CALL Helac_xvvss(zz02,zp0,z11(0,1:4),z21(0,1),z31(1,1),mm,zg,isgn)
				z0(1,1:4)=zz00(1:4)+zz01(1:4)+zz02(1:4)
				CALL Helac_xvvss(zz00,zp0,z11(2,1:4),z21(0,1),z31(0,1),mm,zg,isgn)
				CALL Helac_xvvss(zz01,zp0,z11(0,1:4),z21(2,1),z31(0,1),mm,zg,isgn)
				CALL Helac_xvvss(zz02,zp0,z11(0,1:4),z21(0,1),z31(2,1),mm,zg,isgn)
				z0(2,1:4)=zz00(1:4)+zz01(1:4)+zz02(1:4)
				z0(3,1:4)=0
				DO j=0,3
					DO i=0,3-j
						CALL Helac_xvvss(zz00,zp0,z11(j,1:4),z21(i,1),z31(3-i-j,1),mm,zg,isgn)
						z0(3,1:4)=z0(3,1:4)+zz00(1:4)
					ENDDO
				ENDDO
			CASE DEFAULT
				PRINT *,"Wrong (vvss) of Helac_nextq in Helac_pan1.f90 ! STOP !"
				STOP
			END SELECT
		ENDIF
!        goto99
! 12     continue
    CASE(12) ! ssvv       
        zg=zgssvv(mm,m1,m2,m3)
        CALL Helac_xssvv(z0(0,1),zp0,z11(0,1),z21(0,1:4),z31(0,1:4),mm,zg,isgn)
		IF(Pwavenum.NE.0)THEN
			SELECT CASE(Pwavenum)
			CASE(1)
				CALL Helac_xssvv(zz00(1),zp0,z11(1,1),z21(0,1:4),z31(0,1:4),mm,zg,isgn)
				CALL Helac_xssvv(zz01(1),zp0,z11(0,1),z21(1,1:4),z31(0,1:4),mm,zg,isgn)
				CALL Helac_xssvv(zz02(1),zp0,z11(0,1),z21(0,1:4),z31(1,1:4),mm,zg,isgn)
				z0(1,1)=zz00(1)+zz01(1)+zz02(1)
			CASE(2)
				CALL Helac_xssvv(zz00(1),zp0,z11(1,1),z21(0,1:4),z31(0,1:4),mm,zg,isgn)
				CALL Helac_xssvv(zz01(1),zp0,z11(0,1),z21(1,1:4),z31(0,1:4),mm,zg,isgn)
				CALL Helac_xssvv(zz02(1),zp0,z11(0,1),z21(0,1:4),z31(1,1:4),mm,zg,isgn)
				z0(1,1)=zz00(1)+zz01(1)+zz02(1)
				CALL Helac_xssvv(zz00(1),zp0,z11(2,1),z21(0,1:4),z31(0,1:4),mm,zg,isgn)
				CALL Helac_xssvv(zz01(1),zp0,z11(0,1),z21(2,1:4),z31(0,1:4),mm,zg,isgn)
				CALL Helac_xssvv(zz02(1),zp0,z11(0,1),z21(0,1:4),z31(2,1:4),mm,zg,isgn)
				z0(2,1)=zz00(1)+zz01(1)+zz02(1)
				z0(3,1:4)=0
				DO j=0,3
					DO i=0,3-j
						CALL Helac_xssvv(zz00(1),zp0,z11(j,1),z21(i,1:4),z31(3-i-j,1:4),mm,zg,isgn)
						z0(3,1)=z0(3,1)+zz00(1)
					ENDDO
				ENDDO
			CASE DEFAULT
				PRINT *,"Wrong (ssvv) of Helac_nextq in Helac_pan1.f90 ! STOP !"
				STOP
			END SELECT
		ENDIF
!        goto99 
! 13     continue
    CASE(13) ! s3        
        zg=zgs3(mm,m1,m2)
        CALL Helac_xs3(z0(0,1),zp0,z11(0,1),z21(0,1),mm,zg,isgn)
		IF(Pwavenum.NE.0)THEN
			SELECT CASE(Pwavenum)
			CASE(1)
				CALL Helac_xs3(zz01(1),zp0,z11(1,1),z21(0,1),mm,zg,isgn)
				CALL Helac_xs3(zz02(1),zp0,z11(0,1),z21(1,1),mm,zg,isgn)
				z0(1,1)=zz01(1)+zz02(1)
			CASE(2)
				CALL Helac_xs3(zz01(1),zp0,z11(1,1),z21(0,1),mm,zg,isgn)
				CALL Helac_xs3(zz02(1),zp0,z11(0,1),z21(1,1),mm,zg,isgn)
				z0(1,1)=zz01(1)+zz02(1)
				CALL Helac_xs3(zz01(1),zp0,z11(2,1),z21(0,1),mm,zg,isgn)
				CALL Helac_xs3(zz02(1),zp0,z11(0,1),z21(2,1),mm,zg,isgn)
				z0(2,1)=zz01(1)+zz02(1)
				z0(3,1:4)=0
				DO j=0,3
					CALL Helac_xs3(zz00(1),zp0,z11(j,1),z21(3-j,1),mm,zg,isgn)
					z0(3,1)=z0(3,1)+zz00(1)
				ENDDO
			CASE DEFAULT
				PRINT *,"Wrong (s3) of Helac_nextq in Helac_pan1.f90 ! STOP !"
				STOP
			END SELECT
		ENDIF
!        goto99
! 14     continue
    CASE(14) ! s4      
        zg=zgs4(mm,m1,m2,m3)
        CALL Helac_xs4(z0(0,1),zp0,z11(0,1),z21(0,1),z31(0,1),mm,zg,isgn)
		IF(Pwavenum.NE.0)THEN
			SELECT CASE(Pwavenum)
			CASE(1)
				CALL Helac_xs4(zz00(1),zp0,z11(1,1),z21(0,1),z31(0,1),mm,zg,isgn)
				CALL Helac_xs4(zz01(1),zp0,z11(0,1),z21(1,1),z31(0,1),mm,zg,isgn)
				CALL Helac_xs4(zz02(1),zp0,z11(0,1),z21(0,1),z31(1,1),mm,zg,isgn)
				z0(1,1)=zz00(1)+zz01(1)+zz02(1)
			CASE(2)
				CALL Helac_xs4(zz00(1),zp0,z11(1,1),z21(0,1),z31(0,1),mm,zg,isgn)
				CALL Helac_xs4(zz01(1),zp0,z11(0,1),z21(1,1),z31(0,1),mm,zg,isgn)
				CALL Helac_xs4(zz02(1),zp0,z11(0,1),z21(0,1),z31(1,1),mm,zg,isgn)
				z0(1,1)=zz00(1)+zz01(1)+zz02(1)
				CALL Helac_xs4(zz00(1),zp0,z11(2,1),z21(0,1),z31(0,1),mm,zg,isgn)
				CALL Helac_xs4(zz01(1),zp0,z11(0,1),z21(2,1),z31(0,1),mm,zg,isgn)
				CALL Helac_xs4(zz02(1),zp0,z11(0,1),z21(0,1),z31(2,1),mm,zg,isgn)
				z0(2,1)=zz00(1)+zz01(1)+zz02(1)
				z0(3,1:4)=0
				DO j=0,3
					DO i=0,3-j
						CALL Helac_xs4(zz00(1),zp0,z11(j,1),z21(i,1),z31(3-i-j,1),mm,zg,isgn)
						z0(3,1)=z0(3,1)+zz00(1)
					ENDDO
				ENDDO
			CASE DEFAULT
				PRINT *,"Wrong (s4) of Helac_nextq in Helac_pan1.f90 ! STOP !"
				STOP
			END SELECT
		ENDIF
!        goto99 
! 15     continue
    CASE(15) ! ssv       
        DO j=1,4
           zp1(j)=zyq(j1,j)
        ENDDO
        zg=zgssv(mm,m1,m2)
        CALL Helac_xssv(z0(0,1),zp0,z11(0,1),zp1,z21(0,1:4),mm,zg,isgn)
		IF(Pwavenum.NE.0)THEN
			SELECT CASE(Pwavenum)
			CASE(1)
				CALL MomInclude(ii,ip1array)
				CALL MomInclude(ii1,ip2array)
				CALL Helac_xssvd(zz00(1),zp0,z11(0,1),zp1,z21(0,1:4),mm,zg,isgn,1,ip1array,ip2array)
				CALL Helac_xssv(zz01(1),zp0,z11(1,1),zp1,z21(0,1:4),mm,zg,isgn)
				CALL Helac_xssv(zz02(1),zp0,z11(0,1),zp1,z21(1,1:4),mm,zg,isgn)
				z0(1,1)=zz00(1)+zz01(1)+zz02(1)
			CASE(2)
				CALL MomInclude(ii,ip1array)
				CALL MomInclude(ii1,ip2array)
				CALL Helac_xssvd(zz00(1),zp0,z11(0,1),zp1,z21(0,1:4),mm,zg,isgn,1,ip1array,ip2array)
				CALL Helac_xssv(zz01(1),zp0,z11(1,1),zp1,z21(0,1:4),mm,zg,isgn)
				CALL Helac_xssv(zz02(1),zp0,z11(0,1),zp1,z21(1,1:4),mm,zg,isgn)
				z0(1,1)=zz00(1)+zz01(1)+zz02(1)
				CALL Helac_xssvd(zz00(1),zp0,z11(0,1),zp1,z21(0,1:4),mm,zg,isgn,2,ip1array,ip2array)
				CALL Helac_xssv(zz01(1),zp0,z11(2,1),zp1,z21(0,1:4),mm,zg,isgn)
				CALL Helac_xssv(zz02(1),zp0,z11(0,1),zp1,z21(2,1:4),mm,zg,isgn)
				z0(2,1)=zz00(1)+zz01(1)+zz02(1)
				z0(3,1:4)=0
				DO j=0,3
					CALL Helac_xssv(zz00(1),zp0,z11(j,1),zp1,z21(3-j,1:4),mm,zg,isgn)
					z0(3,1)=z0(3,1)+zz00(1)
				ENDDO
				DO j=1,2
					CALL Helac_xssvd(zz01(1),zp0,z11(j,1),zp1,z21(0,1:4),mm,zg,isgn,3-j,ip1array,ip2array)
					CALL Helac_xssvd(zz02(1),zp0,z11(0,1),zp1,z21(j,1:4),mm,zg,isgn,3-j,ip1array,ip2array)
					z0(3,1)=z0(3,1)+zz01(1)+zz02(1)
				ENDDO
			CASE DEFAULT
				PRINT *,"Wrong (ssv) of Helac_nextq in Helac_pan1.f90 ! STOP !"
				STOP
			END SELECT
		ENDIF
!        goto99
! 16     continue
    CASE(16) !vss     
        DO j=1,4
           zp1(j)=zyq(j1,j)
           zp2(j)=zyq(j2,j)
        ENDDO
        zg=zgvss(mm,m1,m2)
        CALL Helac_xvss(z0(0,1:4),zp0,z11(0,1),zp1,z21(0,1),zp2,mm,zg,isgn)
		IF(Pwavenum.NE.0)THEN
			SELECT CASE(Pwavenum)
			CASE(1)
				CALL MomInclude(ii1,ip1array)
				CALL MomInclude(ii2,ip2array)
				CALL Helac_xvssd(zz00,zp0,z11(0,1),zp1,z21(0,1),zp2,mm,zg,isgn,1,ip1array,ip2array)
				CALL Helac_xvss(zz01,zp0,z11(1,1),zp1,z21(0,1),zp2,mm,zg,isgn)
				CALL Helac_xvss(zz02,zp0,z11(0,1),zp1,z21(1,1),zp2,mm,zg,isgn)
				z0(1,1:4)=zz00(1:4)+zz01(1:4)+zz02(1:4)
			CASE(2)
				CALL MomInclude(ii1,ip1array)
				CALL MomInclude(ii2,ip2array)
				CALL Helac_xvssd(zz00,zp0,z11(0,1),zp1,z21(0,1),zp2,mm,zg,isgn,1,ip1array,ip2array)
				CALL Helac_xvss(zz01,zp0,z11(1,1),zp1,z21(0,1),zp2,mm,zg,isgn)
				CALL Helac_xvss(zz02,zp0,z11(0,1),zp1,z21(1,1),zp2,mm,zg,isgn)
				z0(1,1:4)=zz00(1:4)+zz01(1:4)+zz02(1:4)
				CALL Helac_xvssd(zz00,zp0,z11(0,1),zp1,z21(0,1),zp2,mm,zg,isgn,2,ip1array,ip2array)
				CALL Helac_xvss(zz01,zp0,z11(2,1),zp1,z21(0,1),zp2,mm,zg,isgn)
				CALL Helac_xvss(zz02,zp0,z11(0,1),zp1,z21(2,1),zp2,mm,zg,isgn)
				z0(2,1:4)=zz00(1:4)+zz01(1:4)+zz02(1:4)
				z0(3,1:4)=0
				DO j=0,3
					CALL Helac_xvss(zz00,zp0,z11(j,1),zp1,z21(3-j,1),zp2,mm,zg,isgn)
					z0(3,1:4)=z0(3,1:4)+zz00(1:4)
				ENDDO
				DO j=1,2
					CALL Helac_xvssd(zz01,zp0,z11(j,1),zp1,z21(0,1),zp2,mm,zg,isgn,3-j,ip1array,ip2array)
					CALL Helac_xvssd(zz02,zp0,z11(0,1),zp1,z21(j,1),zp2,mm,zg,isgn,3-j,ip1array,ip2array)
					z0(3,1:4)=z0(3,1:4)+zz01(1:4)+zz02(1:4)
				ENDDO
			CASE DEFAULT
				PRINT *,"Wrong (vss) of Helac_nextq in Helac_pan1.f90 ! STOP !"
				STOP
			END SELECT
		END IF
    END SELECT
!        goto99
! 99     continue
    DO j=1,4
        zy(0:2**Pwavenum-1,j0,j)=zy(0:2**Pwavenum-1,j0,j)+z0(0:2**Pwavenum-1,j)
    ENDDO

    log1=(list(k,5,kc).EQ.list(k,6,kc).OR.list(k,5,kc).EQ.0)
    IF(k.LT.kcou(kc))THEN
        log2=list(k,4,kc).NE.list(k+1,4,kc)
    ELSE
        log2=.TRUE.
    ENDIF
    log3=k.LE.kcou(kc)
    IF(log1.AND.log2.AND.log3)THEN
        DO j=1,4
            z0(0:2**Pwavenum-1,j)=zy(0:2**Pwavenum-1,j0,j)
        ENDDO
        IF(ii.NE.2**n-2)THEN
			zz00(1:4)=z0(0,1:4)
            CALL Helac_propag(z0(0,1:4),mm,zp0)
			IF(Pwavenum.NE.0)THEN
				SELECT CASE(Pwavenum)
				CASE(1)
					zz01(1:4)=z0(1,1:4)
					CALL Helac_propag(zz01,mm,zp0)
					CALL MomInclude(ii,ip1array)
					CALL Helac_propagd(zz00,mm,zp0,1,ip1array)
					z0(1,1:4)=zz01(1:4)+zz00(1:4)
				CASE(2)
					CALL MomInclude(ii,ip1array)
					CALL Helac_propag(z0(3,1:4),mm,zp0)
					DO j=1,2
						zz01(1:4)=z0(3-j,1:4)
						CALL Helac_propagd(zz01,mm,zp0,j,ip1array)
						z0(3,1:4)=z0(3,1:4)+zz01(1:4)
					ENDDO
					zz01(1:4)=zz00(1:4)
					CALL Helac_propagd(zz01,mm,zp0,3,ip1array)
					z0(3,1:4)=z0(3,1:4)+zz01(1:4)
					zz01(1:4)=z0(1,1:4)
					CALL Helac_propag(zz01,mm,zp0)
					zz02(1:4)=zz00(1:4)
					CALL Helac_propagd(zz02,mm,zp0,1,ip1array)
					z0(1,1:4)=zz01(1:4)+zz02(1:4)
					zz01(1:4)=z0(2,1:4)
					CALL Helac_propag(zz01,mm,zp0)
					zz02(1:4)=zz00(1:4)
					CALL Helac_propagd(zz02,mm,zp0,2,ip1array)
					z0(2,1:4)=zz01(1:4)+zz02(1:4)
				CASE DEFAULT
					PRINT *,"Wrong (propag) of Helac_nextq in Helac_pan1.f90! STOP !"
					STOP
				END SELECT
			ENDIF
        ENDIF

        DO j=1,4
            zy(0:2**Pwavenum-1,j0,j)=z0(0:2**Pwavenum-1,j)
        ENDDO 
    ENDIF  
ENDDO
END SUBROUTINE Helac_nextq
! -----------------------------------------------------------------------      
       
SUBROUTINE Helac_ampq(k2,zamp)
IMPLICIT NONE 
INTEGER,INTENT(IN)::k2      
COMPLEX(KIND=DBL),INTENT(OUT)::zamp
INTEGER::i,ijij        
zamp=DCMPLX(dnou(0),dnou(0))
IF(kcou(k2).LT.1)RETURN
IF(ifl(1).LT.30)THEN
! for fermions
   DO i=1,4
     zamp=zamp+zy(2**(Pwavenum)-1,list(kcou(k2),4,k2),i)*zy(0,1,i)
   ENDDO
ELSEIF(ifl(1).GE.31.AND.ifl(1).LE.40)THEN
! for v-bosons
   zamp=Helac_zprod(zy(2**(Pwavenum)-1,list(kcou(k2),4,k2),1:4),zy(0,1,1:4))
ELSEIF(ifl(1).GE.41)THEN
! for scalars
   zamp=zy(2**(Pwavenum)-1,list(kcou(k2),4,k2),1)*zy(0,1,1)
ENDIF  
END SUBROUTINE Helac_ampq
END MODULE Helac_pan1
