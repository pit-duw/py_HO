MODULE MC_PARNI_Weight
!USE Helac_Global  ! use DBL definition and tmp_dir in this module
IMPLICIT NONE
!********************************************************************
!                                                                   *
!                           This is PARNI                           *
!             (Practical Adaptive Random Number Idealizer)          *
!                                                                   *
! author: Andre van Hameren                                         *
!   date: Sat Apr 19 2003                                           *
!********************************************************************
!
!********************************************************************
! usage:
!      call initparni(ndim,readchan)
!
!      do iev=1,nev
!         
!        generate_ndim_random_numbers(rho,ndim)
!        generate_one_random_number(rch)
!         
!        call gnrtparni(xx,rho,rch)
!
!        evaluate_your_integrand(ff,xx)
!
!        call wghtparni(wght,xx)
!
!        call adapparni(ff*wght,2d0)
!
!        apply_your_statistics(ff*wght)
!
!      enddo
!
! initparni, input:
!     ndim = the number of dimensions of the hypercube you want to 
!            generate random vectors in
! readchan = a logical, default is .false.. 
!            If you put readchan=.true. then parni will load channels 
!            from a file "channels.parni" and use them as starting 
!            point. The subroutine "writparni" is available to print 
!            the channels to this file.
! 
! gnrtparni, input:
!      rho = array of ndim (random) numbers between 0 and 1
!      rch = (random) numbers between 0 and 1
!           output:
!       xx = array of ndim (random) numbers between 0 and 1
!
! wghtparni, input: 
!       xx = array of ndim (random) numbers between 0 and 1
!           output:
!     wght = wght coming with xx ( wght=1/density(xx) )
!           
! adapparni, input:
!  ff*wght = total weight, so integrand times "wght" from parni
!    expon = 2d0. Exponent in adaptive multi-channeling, should be
!            put to 1d0 for the creation of a unitary probability
!            decomposition.
!
! For now, I suggest the user to play around with the parameters in 
! the subroutine "initparni" below in order to get an optimal result.
!********************************************************************
REAL(KIND(1d0)),DIMENSION(5555,37),PUBLIC::bxm,bxp ! the lower and upper corrdi for every channel
REAL(KIND(1d0)),DIMENSION(5555),PUBLIC::vol   ! the volume for every channel
REAL(KIND(1d0)),DIMENSION(5555),PUBLIC::ave
INTEGER,PUBLIC:: ndat,ndatmin,nchmax,niter,nmerge,istp,nstp
REAL(KIND(1d0)),DIMENSION(5555),PUBLIC::wch ! the weight for every channel
INTEGER,PUBLIC::nch                     ! the number of actually channels
INTEGER,DIMENSION(5555),PUBLIC::lch     ! label for every channel
INTEGER,PUBLIC::nnd                     ! the number of dimensions (for public)
REAL(KIND(1d0)),PUBLIC::sumdns
INTEGER,PUBLIC::nchloc
REAL(KIND(1d0)),DIMENSION(5555),PUBLIC::loc
LOGICAL,PUBLIC::debug
INTEGER,PUBLIC::NPRN=1                  ! printed or not from MC_Vegas
                                        ! If NPRN<0 don't print
INTEGER,PUBLIC::NDEV=6                  ! device number for output
SAVE NPRN,NDEV
CONTAINS
SUBROUTINE INITPARNI(nndI,readchan)
!  **************************************************************
!  **************************************************************
IMPLICIT NONE
INTEGER,INTENT(IN)::nndI                ! the number of dimensions of the hypercube you want to 
                                        ! generate random vectors in
REAL(KIND(1d0)),DIMENSION(5555)::arr
INTEGER::iid,ich,iter
LOGICAL,INTENT(IN)::readchan       ! a logical, default is .false.. 
                                   ! If you put readchan=.true. then parni will load channels 
                                   ! from a file "channels.parni" and use them as starting 
                                   ! point. The subroutine "writparni" is available to print 
                                   ! the channels to this file.
! put input parameters
nnd = nndI ! dimension of hypercube
! set further parameters
nchmax  = 1000 ! maximum number of channels (=boxes)
ndatmin = 1000 ! batch-size of data for weight optimization
niter   = 40   ! # boxes that is splitted in one optimization step
nmerge  = 20   ! # channels that is discard in one optimization step
nstp    = 1    ! # weight optimization steps before channel optimization
IF(NPRN.GE.0)WRITE(NDEV,101) nchmax,ndatmin,niter,nmerge,nstp 
101 FORMAT('##### MC_PARNI #####'&
     ,' nch=',I4,', ndat=',I5,', niter=',I3&
     ,', nmerge=',I3,', nstp=',I2,' #') 
!initialize lch, the labels of the channels
DO ich=1,nchmax
   lch(ich) = ich
ENDDO
! initialize statistics
istp = 0
ndat = 0
DO ich=1,nchmax
   ave(ich) = 0d0
ENDDO
! initialize channels
IF (readchan) THEN
   CALL READPARNI    ! we rarely use it
ELSE
   nch = 1           ! the actual number of boxes
   wch(nch) = 1d0    ! the weight w_k for each boxes
   DO iid=1,nnd
      bxp(nch,iid) = 1d0  ! the upper corrdinates for ~nch box in iid dimension
      bxm(nch,iid) = 0d0  ! the lower corrdinates for ~nch box in iid dimension
   ENDDO
   vol(nch) = 1d0       ! the volume of every box
! For starting with more than 2*nnd channels use iter=1,1+2*nnd  etc
! after this we get the (2*nnd)**3 boxes with equal wch=1d0/((2*nnd)**3)
   DO iter=1,1+2*nnd +(2*nnd)**2
! order channels: wch(lch(1)) is going to be largest
      DO ich=1,nch
          arr(ich) = -wch(lch(ich))
      ENDDO
      CALL SORTPARNI(arr,lch,nch)
! create new channels by dissecting the largest
      CALL dissectfirst    ! split the box
   ENDDO
ENDIF
END SUBROUTINE INITPARNI

SUBROUTINE ADAPPARNI(wght,expon)
!  **************************************************************
!  **************************************************************
IMPLICIT NONE
REAL(KIND(1d0)),INTENT(IN)::wght,expon
REAL(KIND(1d0))::factor,sum
REAL(KIND(1d0)),DIMENSION(5555)::arr
INTEGER::ich,iid,nless,nmore ,iter
LOGICAL::adaptwg,adaptch 
INTEGER::index
! update ave,  loc  labels channels with non-zero density-value
factor = wght**expon/sumdns
DO ich=1,nchloc
   index=loc(ich)
   ave(index) = ave(index) + factor/vol(index)
ENDDO
ndat = ndat + 1
! check if enough data have been collected
adaptwg = (ndat.EQ.ndatmin)
! adapt channel-weights
IF(adaptwg)THEN
   istp = istp+1
   sum = 0d0
   DO ich=1,nch
      ave(lch(ich)) = ave(lch(ich))/DBLE(ndat)
      wch(lch(ich)) = wch(lch(ich))*ave(lch(ich))**(1d0/expon)
      ave(lch(ich)) = 0d0
      sum = sum + wch(lch(ich))
   ENDDO
! normalize wch
   DO ich=1,nch
      wch(lch(ich)) = wch(lch(ich))/sum 
   ENDDO  
! re-initialize ndat
   ndat = 0
ENDIF
! adapt channels
adaptch = (istp.EQ.nstp)
IF(adaptch) THEN
   istp = 0
! order channels: wch(lch(1)) is going to be largest
   DO ich=1,nch
      arr(ich) = -wch(lch(ich))
   ENDDO
   CALL SORTPARNI(arr,lch,nch)
! merge smallest  nmerge  channels 
   IF(nch.GT.nmerge.AND.nmerge.GT.1) CALL mergelastfew(nmerge-1)
! create new channels by dissecting the largest
   CALL dissectfirst
   DO iter=1,niter-1
! order channels: wch(lch(1)) is going to be largest
      DO ich=1,nch
         arr(ich) = -wch(lch(ich))
      ENDDO
      CALL SORTPARNI(arr,lch,nch)
! create new channels by dissecting the largest
      CALL dissectfirst
   ENDDO
ENDIF
END SUBROUTINE ADAPPARNI

SUBROUTINE mergelastfew(nless)
!  **************************************************************
!  * merge  nless+1  smallest channels to 1 channel
!  **************************************************************
IMPLICIT NONE
REAL(KIND(1d0))::wch1
REAL(KIND(1d0)),DIMENSION(37)::bxp1,bxm1
INTEGER,INTENT(IN)::nless
INTEGER::ich,iid
! initialize result of merging
wch1 = 0d0
DO iid=1,nnd
   bxp1(iid) = 0d0
   bxm1(iid) = 1d0
ENDDO
! merge
DO ich=nch-nless,nch 
   wch1 = wch1 + wch(lch(ich))
   DO iid=1,nnd
      IF(bxp(lch(ich),iid).GT.bxp1(iid)) &
           bxp1(iid) = bxp(lch(ich),iid)
      IF(bxm(lch(ich),iid).LT.bxm1(iid)) &
           bxm1(iid) = bxm(lch(ich),iid)
   ENDDO
ENDDO
! put result of merging
ich = lch(nch-nless)
wch(ich) = wch1
vol(ich) = 1d0
DO iid=1,nnd
   bxp(ich,iid) = bxp1(iid)
   bxm(ich,iid) = bxm1(iid)
   vol(ich) = vol(ich)*( bxp(ich,iid)-bxm(ich,iid) )
ENDDO
! update nch
nch = nch-nless
END SUBROUTINE mergelastfew

SUBROUTINE dissectfirst
!  **************************************************************
!  * Replace channel lch(1) by its 2*nnd pieces after dissection 
!  **************************************************************
IMPLICIT NONE
REAL(KIND(1d0))::wch1,vol1
REAL(KIND(1d0)),DIMENSION(37)::bxp1,bxm1
INTEGER::nmore,ich,jch,iid,nmany
INTEGER,DIMENSION(74)::kch
! number of extra channels
nmore = 2*nnd-1
! merge smallest  nmany+1  channels if necessary
nmany = nch+nmore-nchmax
IF(nmany.GT.0) CALL mergelastfew(nmany)
! store the first channel
wch1 = wch(lch(1))/DBLE(2*nnd)
DO iid=1,nnd
   bxp1(iid) = bxp(lch(1),iid)
   bxm1(iid) = bxm(lch(1),iid)
ENDDO
vol1 = vol(lch(1))/2d0
! store last  nmore  entries of  lch
DO ich=1,nmore
  kch(ich) = lch(nch+ich)
ENDDO
! shift lch with  nmore  entries
DO ich=nch,1,-1
  lch(ich+nmore) = lch(ich)
ENDDO
nch = nch+nmore
! re-fill first  nmore  entries of  lch
DO ich=1,nmore
  lch(ich) = kch(ich)
ENDDO
! copy the previously first channel  2*nnd  times
DO ich=1,2*nnd
  wch(lch(ich)) = wch1
  DO iid=1,nnd
      bxp(lch(ich),iid) = bxp1(iid)
      bxm(lch(ich),iid) = bxm1(iid)
  ENDDO
  vol(lch(ich)) = vol1
ENDDO
! cut the copies in two, each one along another direction
DO iid=1,nnd
  ich = 2*iid
  bxp(lch(ich),iid) = ( bxp(lch(ich),iid)& 
   +bxm(lch(ich),iid) )/2d0
  ich = ich-1
  bxm(lch(ich),iid) = ( bxp(lch(ich),iid)&
   +bxm(lch(ich),iid) )/2d0
ENDDO
END SUBROUTINE dissectfirst

SUBROUTINE GNTRPARNI(xx,rho,rch)
!  **************************************************************
!  * "rch" is used to choose a channel (box), 
!  * "rho" is used to construct the point "xx" in the chosen box 
!  **************************************************************
IMPLICIT NONE
REAL(KIND(1d0)),DIMENSION(37),INTENT(OUT)::xx
REAL(KIND(1d0)),DIMENSION(37),INTENT(IN)::rho
REAL(KIND(1d0)),INTENT(IN)::rch
REAL(KIND(1d0))::wchsum,xm,xp
INTEGER::ich,iid,flag=0
CHARACTER(len=24),PARAMETER::tmp_dir="./tmp/"
! choose channel
wchsum = 0d0
DO ich=1,nch
   wchsum = wchsum + wch(lch(ich))
   IF(rch.LT.wchsum)THEN
      flag=1
	  EXIT           ! goto 100
   ENDIF
ENDDO

IF(flag.EQ.0)THEN
  WRITE(NDEV,*) 'Error in GNRTPARNI: no channel chosen.'
  WRITE(NDEV,*) 'Channels dumped in file "dumped.parni"'
  OPEN(99,FILE=TRIM(tmp_dir)//'dumped.parni')
  51 FORMAT('channel',18x,'weight',13x,'sum weights')
  52 FORMAT(I7,2d24.16)
  WRITE(99,51)
  xm = 0d0
  DO ich=1,nch
      xm = xm + wch(lch(ich))
      WRITE(99,52) lch(ich),wch(lch(ich)),xm
  ENDDO
  CLOSE(99)
  STOP
ENDIF
      
!  100 continue
! generate  xx  in box  ich 
DO iid=1,nnd
    xm = bxm(lch(ich),iid)
    xp = bxp(lch(ich),iid)
    xx(iid) = xm + (xp-xm)*rho(iid)
ENDDO
END SUBROUTINE GNTRPARNI

SUBROUTINE WGHTPARNI(wght,xx)
!  **************************************************************   
!  **************************************************************
IMPLICIT NONE
REAL(KIND(1d0)),INTENT(OUT)::wght
REAL(KIND(1d0)),DIMENSION(37),INTENT(IN)::xx
INTEGER::ich,iid,jch,flag
sumdns = 0d0
jch = 0
DO ich=1,nch
   flag=0
   DO iid=1,nnd
      IF(xx(iid).LE.bxm(lch(ich),iid).OR.bxp(lch(ich),iid).LT.xx(iid) )THEN !goto 100
	      flag=1
		  EXIT
	  ENDIF
   ENDDO
   IF(flag.EQ.0)THEN 
      jch = jch+1
      loc(jch) = lch(ich)
      sumdns = sumdns + wch(lch(ich))/vol(lch(ich))
   ENDIF
!  100   continue
ENDDO
nchloc = jch
wght = 1d0/sumdns
END SUBROUTINE WGHTPARNI

 
SUBROUTINE SORTPARNI(arr,lbl,nn)
!  **************************************************************
!  * heapsort:
!  * order the array "arr" of double precisions from small to 
!  * large, and re-order the array "lbl" of integers with it.
!  * So after 
!  *      do ii=1,nn
!  *        arr(ii) = brr(lbl(ii))
!  *      enddo
!  *      call sortparni(arr,lbl,nn)
!  * brr(lbl(ii)) will be ordered from small to large, and after
!  *      do ii=1,nn
!  *        arr(ii) = -brr(lbl(ii))
!  *      enddo
!  *      call sortparni(arr,lbl,nn)
!  * brr(lbl(ii)) will be ordered from large to small.
!  **************************************************************
IMPLICIT NONE 
REAL(KIND(1d0)),DIMENSION(5555),INTENT(INOUT)::arr
REAL(KIND(1d0))::arrb 
INTEGER,DIMENSION(5555),INTENT(INOUT)::lbl
INTEGER,INTENT(IN)::nn
INTEGER::lblb,ll,ii,jj,ir,flag=0
IF(nn.GT.1)THEN   !goto 30     
   ll = nn/2 + 1
   ir = nn
   DO
!  10  continue
      IF (ll.GT.1)THEN 
          ll = ll - 1
          arrb = arr(ll)
          lblb = lbl(ll)
      ELSE
          arrb = arr(ir)
          lblb = lbl(ir)
          arr(ir) = arr(1)
          lbl(ir) = lbl(1)
          ir = ir - 1
          IF(ir.EQ.1)THEN
            arr(1) = arrb
            lbl(1) = lblb
			EXIT
!            goto 30
          ENDIF
       ENDIF
       ii = ll
       jj = ll + ll
!  20    continue
       DO
         IF(jj.LE.ir)THEN
            IF(jj.LT.ir)THEN
               IF(arr(jj).LT.arr(jj+1))jj = jj + 1
            ENDIF
            IF(arrb.LT.arr(jj))THEN
                arr(ii) = arr(jj)
                lbl(ii) = lbl(jj)
                ii = jj
                jj = jj + jj
            ELSE
                jj = ir + 1
            ENDIF
		  ELSE
		     EXIT
!        goto 20
          ENDIF
	    ENDDO
        arr(ii) = arrb
        lbl(ii) = lblb
    ENDDO
ENDIF
!  30  continue
END SUBROUTINE SORTPARNI
      
      
SUBROUTINE WRITPARNI
!  **************************************************************   
!  **************************************************************
IMPLICIT NONE
INTEGER::ich,iid
CHARACTER(len=24),PARAMETER::tmp_dir="./tmp/"
OPEN(99,FILE=TRIM(tmp_dir)//'channels.parni')
WRITE(99,101) nch,nnd
DO ich=1,nch
   WRITE(99,102) lch(ich),wch(lch(ich))&
   ,(bxm(lch(ich),iid),iid=1,nnd),(bxp(lch(ich),iid),iid=1,nnd)
ENDDO
CLOSE(99)
101 FORMAT(I4,I3)
102 FORMAT(I4,75d24.16)
END SUBROUTINE WRITPARNI

SUBROUTINE READPARNI
!  **************************************************************   
!  **************************************************************
IMPLICIT NONE
INTEGER::ich,iid
CHARACTER(len=24),PARAMETER::tmp_dir="./tmp/"
OPEN(99,FILE=TRIM(tmp_dir)//'channels.parni')
READ(99,101) nch,nnd
DO ich=1,nch
    READ(99,102) lch(ich),wch(lch(ich))&
    ,(bxm(lch(ich),iid),iid=1,nnd),(bxp(lch(ich),iid),iid=1,nnd)
ENDDO
CLOSE(99)
101 FORMAT(I4,I3)
102 FORMAT(I4,75d24.16)
END SUBROUTINE READPARNI
END MODULE MC_PARNI_Weight
