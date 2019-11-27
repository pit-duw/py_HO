MODULE Helac_Feynman
USE Helac_Global
USE Helac_Func_1
IMPLICIT NONE
CONTAINS
  SUBROUTINE Helac_feyn(n1,k1,maxmul,ico)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!                                                                     
!             Feynman graph finding algorithm                         
!                                                                     
!                                                                     
!    The subroutine redo results to ico(k1,1:18) elements             
!    that have the following description :                            
!                                                                     
!    Take an example  from w- w+ -> w- w+ a a                         
!    which has 35 subamplitudes coming out form 156 contributions     
!    resulting to 192 Feynman graphs (unitary-gauge)                  
!                                                                     
!  156  5  62  33  35  1  44  56  34  28  2 33  2  4 33  3 0 1 0        
!                                                                     
!                                                                     
!    1/18 : the routine that has been used to produce it  (5)         
!    2/18 : the momentum i0  (62)                                    
!    3/18 : the flavour  m0  (33)                                     
!    4/18 : the rank of the subamplitude (35 ;the 35-th subamplitude) 
!    5/18 : the rank of the multiplicity (1 ;the 1st contribution)    
!    6/18 : the max of the multiplicity (44; 44 contributions)        
!    7/18 : the momentum i1 (56)                                      
!    8/18 : the flavour  m1 (34)                                      
!    9/18 : the rank of its subamplitude (28;it comes form the 28-th) 
!   10/18 : the momentum i2 (2)                                       
!   11/18 : the flavour  m2 (33)                                      
!   12/18 : the rank of its subamplitude (2;it comes form the 2-nd)   
!   13/18 : the momentum i3 (4)                                       
!   14/18 : the flavour  m3 (33)                                      
!   15/18 : the rank of its subamplitude (4;it comes form the 4-th)   
!   16/18 : the chirality (0); relevant for fermions only  (1/2)      
!   17/18 : the isign (1); from fermions permutations                 
!   18/18 : igluon                                                    
!                                                                     
!                                                                     
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  IMPLICIT NONE
! include 'common_feyn.h'
! include 'common_print.h'
  INTEGER,INTENT(IN)::k1
  INTEGER,DIMENSION(k1,18),INTENT(INOUT)::ico
  INTEGER,DIMENSION(:,:),ALLOCATABLE::icoo
  INTEGER,DIMENSION(:),ALLOCATABLE::ik
  INTEGER,DIMENSION(:),ALLOCATABLE::ipa
  INTEGER,DIMENSION(:),ALLOCATABLE::m,f,isani
  INTEGER,INTENT(IN)::maxmul,n1
  ! there is some problem with is (allocate?)
  INTEGER::nn,istat,k,i,i1,j,j1,ix,l,ixch,iv,flag1,flag2,flag0,istop,il2,initinit3=0 !,tgtg
  SAVE ik,icoo,ipa,m,f,isani,initinit3
  WRITE(*,*)' '
  !WRITE(*,*)'k1=',k1,ico(k1,4),ico(k1,6),maxmul,n1
  WRITE(*,*)'The number of still alived Green functions(subamplitudes) k1=',k1
  WRITE(*,*)'The number of still alived (i,m) with level 1 ico(k1,4)=',ico(k1,4)
  WRITE(*,*)'The number of the last (i0.m0) spllitings ico(k1,6)=',ico(k1,6) 
  WRITE(*,*)'The upper limit number of splitting for any (i,m) maxmul =',maxmul
  WRITE(*,*)'The number of external legs n=',n1
  WRITE(*,*)' '
  nn=n1
  IF(GLOBALINIT_feyn.EQ.0)THEN
	initinit3=0
  ENDIF
  IF(initinit3.EQ.0)THEN
	ng=0
	is(1:ISMAX,1:11,1:8)=0
	initinit3=1
  ENDIF
  ALLOCATE(icoo(0:ico(k1,4),0:maxmul),ik(0:ico(k1,4)),ipa(n1),&
          m(nn),f(nn),isani(2:2**(n1-1)),STAT=istat)
  IF(istat.NE.0)THEN
       WRITE(*,*)'warning: allocation is not working properly in Helac_feyn'
       STOP
  ENDIF
  ipa(1:n1)=0
  m(1:nn)=0
  f(1:nn)=0
  icoo(0:ico(k1,4),0:maxmul)=0
  ik(0:ico(k1,4))=1

  IF(maxmul.GT.0)THEN
     k=1
     DO i1=1,k1                   ! delete the same configuration of splitting
       IF(ico(k,5).EQ.0)THEN
          ico(k:k1-1,1:18)=ico(k+1:k1,1:18)
          k=k-1
       ENDIF
       k=k+1
     ENDDO
      
     k=k-1                        ! the number of different splitting configuarations
     IF(k.NE.0)THEN               ! k=0 means there is no contribution to this color configuration
       DO j=1,k
          ico(j,5)=ico(j,6)-ico(j,5)+1
          icoo(ico(j,4),ico(j,5))=j
       ENDDO
       ik(0:ico(k,4))=1
       j1=k
       ix=1
       l=1
       ixch=n1-1
       
	   
       levelm1:DO
	      flag0=0                          !11
          ipa(ix)=icoo(ico(j1,4),ik(ico(j1,4)))
          flag0=1
          
          level0:DO
! 10    continue
             flag1=0
             level1:DO
               IF(ixch.EQ.0.AND.flag0.NE.1)THEN   !goto13
		          flag1=1
		          EXIT
		       ELSE
			      IF(flag0.NE.1)THEN
                    j1=ipa(ix)
! 12    continue
                  ENDIF
			      flag0=0
                  IF(ix.EQ.ixch)THEN
                     ik(ico(j1,4))=ik(ico(j1,4))+1
                     IF(ik(ico(j1,4)).LE.ico(j1,6))THEN
                        ipa(ix)=icoo(ico(j1,4),ik(ico(j1,4)))
                     ENDIF
                  ENDIF
                  IF(ik(ico(j1,4)).GT.ico(j1,6))THEN
                     ik(ico(j1,4))=1
                     ix=1
                     ixch=ixch-1
                     l=1
                     f(1:nn)=0
                     m(1:nn)=0
		          ELSE
		             EXIT
!        goto10
                  ENDIF
		        ENDIF
		     ENDDO level1
		
             if1:IF(flag1.EQ.1)THEN
		        EXIT
		     ELSE if1
                j1=icoo(ico(j1,4),ik(ico(j1,4)))
                iv=1
                m(l)  =ico(j1,7)
                f(l)  =ico(j1,9)
                m(l+1)=ico(j1,10)
                f(l+1)=ico(j1,12)
                m(l+2)=ico(j1,13)
                f(l+2)=ico(j1,15)
                IF(Helac_level(n1,m(l+2)).LE.1.OR.icoo(f(l+2),ik(f(l+2))).EQ.0)THEN
                   m(l+2:nn-1)=m(l+3:nn)
                   m(nn)=0
                   f(l+2:nn-1)=f(l+3:nn)
                   f(nn)=0
                ENDIF
                IF(Helac_level(n1,m(l+1)).LE.1.OR.icoo(f(l+1),ik(f(l+1))).EQ.0)THEN
                   m(l+1:nn-1)=m(l+2:nn)
                   m(nn)=0
                   f(l+1:nn-1)=f(l+2:nn)
                   f(nn)=0
                ENDIF
                IF(Helac_level(n1,m(l)).LE.1.OR.icoo(f(l),ik(f(l))).EQ.0)THEN
                   m(l:nn-1)=m(l+1:nn)
                   m(nn)=0
                   f(l:nn-1)=f(l+1:nn)
                   f(nn)=0
                ENDIF
                j=0
                istop=0
                DO WHILE(istop.EQ.0)
                   j=j+1
                   IF(m(j).EQ.0)istop=1
                END DO
                l=j-1                ! l=0 means we get the level 2 or can't go backward.

                IF (l.GT.0)THEN
                   ix=ix+iv
                   IF(iv.EQ.2)ipa(ix-1)=0
                ENDIF
                IF(ix.LE.0.OR.ix.GT.n1-1) WRITE(*,*)'warning in Helac_feyn',ix

                if2:IF(l.EQ.0)THEN
                   l=1
                   ng=ng+1
				   IF(ng.GT.ISMAX)THEN
					  PRINT *,'The number of Feynman diagrams is lager than ISMAX=',ISMAX
					  PRINT *,'Please enlarge the ISMAX in Feynman_Helac.f90'
					  STOP 
				   ENDIF
                   j=0
                   i=0
	               isani(2:2**(n1-1))=0
		           DO
                     j=j+1               !101
	                 i=i+1
                     IF(i.GT.ix)THEN    !goto102
			            EXIT
			         ELSE
                        IF(ipa(i).NE.0)THEN
                          IF(ico(ipa(i),18).LE.1)THEN  ! start
                             is(ng,j,1:2)=ico(ipa(i),2:3)
                             is(ng,j,3:4)=ico(ipa(i),7:8)
                             is(ng,j,5:6)=ico(ipa(i),10:11)
                             is(ng,j,7:8)=ico(ipa(i),13:14)
	                         IF(Helac_level(n1,is(ng,j,1)).EQ.1)isani(is(ng,j,1))=isani(is(ng,j,1))+1
	                         IF(Helac_level(n1,is(ng,j,3)).EQ.1)isani(is(ng,j,3))=isani(is(ng,j,3))+1
	                         IF(Helac_level(n1,is(ng,j,5)).EQ.1)isani(is(ng,j,5))=isani(is(ng,j,5))+1
	                         IF(Helac_level(n1,is(ng,j,7)).EQ.1)isani(is(ng,j,7))=isani(is(ng,j,7))+1 
				  !        goto101
				          ELSE
				             EXIT                   
                          ENDIF                    ! end
	                    ELSEIF(ipa(i).EQ.0)THEN
	                      j=j-1
                                          !goto101
			            ELSE
			              EXIT
	                    ENDIF
			          ENDIF
		           ENDDO
! 102    continue
                   flag2=0
                   DO il2=2,n1
	                  IF(isani(2**(il2-1)).NE.1)THEN
                         is(ng,1:8,1:8)=0
                         ng=ng-1
                         ixch=ix
                         ix=1
                         f(1:nn)=0
                         m(1:nn)=0
					     flag2=1
                         EXIT                     ! goto10
                      ENDIF
	               ENDDO
				   IF(flag2.EQ.0)THEN

! 1001   continue
                      ixch=ix
                      ix=1
                      f(1:nn)=0
                      m(1:nn)=0
					  !        goto10
			       ENDIF
!         endif
                 ELSE if2
			       EXIT
                 ENDIF if2
             ENDIF if1
	  ENDDO level0
      IF(flag1.NE.1)THEN
        j1=icoo(f(l),ik(f(l)))                 ! goto 11
      ELSE
	     EXIT
	  ENDIF
	 ENDDO levelm1
	ENDIF
  ENDIF
  WRITE(nunit1,*)' '
  WRITE(nunit1,*) '(Feynman_Helac.f90) the number of Feynman graphs is:',ng    !13
!  READ(*,*)tgtg
  WRITE(nunit1,*)' '
      
  DEALLOCATE(icoo,ik,ipa,m,f,isani)      
  END SUBROUTINE Helac_feyn
END MODULE Helac_Feynman