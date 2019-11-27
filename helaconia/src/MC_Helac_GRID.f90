MODULE MC_Helac_GRID
USE Helac_ranmar_mod
IMPLICIT NONE
!**********************************************************************
!*                                                                    *
!*                     Adaptive grids for Helac                       *
!*                                                                    *
!* author: Andreas van Hameren <hamerenREMOVETHIS@ifj.edu.pl>         *
!*   date: 11-05-2010                                                 *
!*                                                                    *
!* The first module "avh_helac_grid" defines the grid-type and the    *
!* routines acting on it.                                             *
!* The second module "helac_grids" contains an allocatable array of   *
!* "grids" to be used in Helac.                                       *
!* At the end of this file are the routines to be actually called in  *
!* Helac:                                                             *
!*                                                                    *
!* subroutine helac_grid_alloc_2( nsize )                             *
!* subroutine helac_grid_init_2( ID ,nbatch,nchmax,frac ,igraph,iinv )*
!* subroutine helac_grid_gnrt_2( ID ,xx )                             *
!* subroutine helac_grid_wght_2( ID ,ww ,xx )                         *
!* subroutine helac_grid_collect_2( ww )                              *
!* subroutine helac_grid_plot_2( iunit )                              *
!*                                                                    *
!* Since this file contians fortran90 modules, you need to make sure  *
!* the makefile contains lines like                                   *
!*                                                                    *
!#%.o : %.f90
!# 	$(FC) $(FFLAGS) -c  $*.f90 -o $*.o
!*                                                                    *
!**********************************************************************
PRIVATE
PUBLIC::Helac_grid_type,Helac_grid_init &
        ,Helac_grid_close,Helac_grid_gnrt &
        ,Helac_grid_wght,Helac_grid_collect &
        ,Helac_grid_plot,Helac_grid_random,&
		grid,active,nsize,NCHM,Helac_grid_alloc_2&
		,Helac_grid_init_2,Helac_grid_gnrt_2&
		,Helac_grid_wght_2,Helac_grid_collect_2&
		,Helac_grid_plot_2,Helac_grid_close_2,&
		Helac_grid_remove_2
! Character lenght label
INTEGER,PARAMETER:: nlbl = 8
! Message unit
INTEGER,PARAMETER:: nunit = 6
INTEGER,PARAMETER::NCHM=10000
TYPE :: Helac_grid_type
    PRIVATE
    REAL(KIND(1d0)),DIMENSION(0:NCHM):: xri,sch
	REAL(KIND(1d0)),DIMENSION(1:NCHM)::wch
!	REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE:: wch
!	REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE:: sch
!    ALLOCATABLE::xri,wch,sch
    INTEGER:: idat,ndat,nch,nchmax
    REAL(KIND(1d0)):: frac
    INTEGER:: iw,ic
    CHARACTER(nlbl):: lbl
END TYPE Helac_grid_type
TYPE(Helac_grid_type),DIMENSION(:),ALLOCATABLE,SAVE:: grid
LOGICAL,DIMENSION(:),ALLOCATABLE,SAVE:: active
INTEGER,SAVE:: nsize=0
CONTAINS
SUBROUTINE Helac_grid_random( xx )
!**********************************************************************
! The random number generator
!**********************************************************************
IMPLICIT NONE
REAL(KIND(1d0)) ,INTENT(OUT)::xx
xx =Helac_rnmy(0)
END SUBROUTINE Helac_grid_random
!
!**********************************************************************
!**********************************************************************
!*
!* Now following defining the grid-type and the routines
!* acting on it.
!*
!**********************************************************************
!**********************************************************************
SUBROUTINE Helac_grid_init( grid1 ,nbatch,nchmax,frac ,lbl )
!********************************************************************
!* nbatch = number of events to be collected before adaptation step
!* nchmax = the maximal number of bins/channels
!*   frac = fraction of the number of channels/bins to be merged in
!*          an adaptation step
!*    lbl = a character label for the grid
!********************************************************************
IMPLICIT NONE
TYPE(Helac_grid_type),INTENT(INOUT)::grid1
INTEGER,INTENT(IN):: nbatch,nchmax
REAL(KIND(1d0)),INTENT(IN):: frac
CHARACTER(nlbl),INTENT(IN):: lbl
INTEGER:: ii
grid1%lbl  = lbl
grid1%ndat = nbatch
grid1%idat = 0
grid1%nchmax = nchmax
grid1%frac   = frac
grid1%nch = 1
grid1%iw = 0
grid1%ic = 0
IF(grid1%nchmax.LE.1)THEN
  WRITE(*,*) 'ERROR in Helac_grid_init: nchmax should be at least 2'
  STOP
ENDIF
!ALLOCATE( grid1%xri(0:nchmax) )
!ALLOCATE( grid1%wch(1:nchmax) )
!ALLOCATE( grid1%sch(0:nchmax) )
grid1%xri(0) = 0d0
grid1%sch(0) = 0d0
DO ii=1,grid1%nch
   grid1%xri(ii) = DBLE(ii)/DBLE(grid1%nch) ! the coordinate of the division
   grid1%wch(ii) = 0d0                    ! weightting for every channel
   grid1%sch(ii) = DBLE(ii)/DBLE(grid1%nch)   ! the importance weighting for every channel
                                              ! i.e.when grid1%sch(i)-grid1%sch(i-1) is larger
											  ! ,then the probability for choosing the grid1%xri(i)
											  ! to grid1%xri(i-1) larger.
ENDDO
END SUBROUTINE Helac_grid_init

SUBROUTINE Helac_grid_close( grid1 )
!********************************************************************
!********************************************************************
IMPLICIT NONE
TYPE(Helac_grid_type),INTENT(INOUT)::grid1
!IF(ALLOCATED(grid1%xri))DEALLOCATE(grid1%xri)
!IF(ALLOCATED(grid1%wch))DEALLOCATE(grid1%wch)
!IF(ALLOCATED(grid1%sch))DEALLOCATE(grid1%sch )
RETURN
END SUBROUTINE Helac_grid_close

SUBROUTINE Helac_grid_gnrt( grid1 ,xx )
!********************************************************************
! xx is output
!********************************************************************
IMPLICIT NONE
TYPE(Helac_grid_type),INTENT(INOUT) :: grid1
REAL(KIND(1d0)),INTENT(OUT):: xx
INTEGER:: i0,i1,ii
CALL Helac_grid_random( xx )
i0 = 0
i1 = grid1%nch
DO WHILE(i1-i0.GT.1)
   ii = (i0+i1)/2
   IF(xx.LE.grid1%sch(ii))THEN
      i1 = ii
   ELSE
      i0 = ii
   ENDIF
ENDDO
ii = i1
grid1%iw = ii
grid1%ic = ii
CALL Helac_grid_random( xx )
xx = ( grid1%xri(ii)-grid1%xri(ii-1) )*xx + grid1%xri(ii-1)
END SUBROUTINE Helac_grid_gnrt

SUBROUTINE Helac_grid_wght( grid1 ,weight ,xx )
!********************************************************************
! Put xx=-1d0 to use the latest generated, and stored, channel
!********************************************************************
IMPLICIT NONE
TYPE(Helac_grid_type),INTENT(INOUT)::grid1
REAL(KIND(1d0)),INTENT(OUT):: weight
REAL(KIND(1d0)),INTENT(IN):: xx
INTEGER:: ii
IF(xx.LT.0d0)THEN
    ii = grid1%iw 
ELSE
    CALL Helac_findch( grid1 ,ii ,xx )
    IF(ii.EQ.0) THEN
         IF(nunit.GT.0)WRITE(nunit,*) &
        'ERROR in Helac_grid_wght: no channel found, putting weight=0'
         weight = 0d0
         RETURN
    ENDIF
    grid1%ic = ii
ENDIF
weight = ( grid1%xri(ii) - grid1%xri(ii-1) ) &
         / ( grid1%sch(ii) - grid1%sch(ii-1) )
END SUBROUTINE Helac_grid_wght

SUBROUTINE Helac_grid_collect( grid1 ,weight ,xx )
!********************************************************************
! Put xx=-1d0 to update the latest active, and stored, channel
!********************************************************************
IMPLICIT NONE
TYPE(Helac_grid_type),INTENT(INOUT) :: grid1
REAL(KIND(1d0)),INTENT(IN):: weight,xx
INTEGER:: ii,jj
LOGICAL:: update

! Find the bin into which xx belongs.
! If xx<0, use the latest active (by gnrt or wght) bin.
IF(xx.LT.0d0)THEN
    ii = grid1%ic
ELSE
    CALL Helac_findch( grid1 ,ii ,xx )
ENDIF

! Make sure the weight coming with a channel is used only once
IF(ii.EQ.0)RETURN
grid1%ic = 0

! Update the weights of the bins
IF(weight.NE.0d0)THEN
    grid1%wch(ii) = grid1%wch(ii) + weight    ! entropy
!   grid%wch(ii) = grid%wch(ii) + weight**2 ! variance
    grid1%idat = grid1%idat+1
ENDIF

! Adapt the bins
update = .FALSE.
IF(grid1%idat.EQ.grid1%ndat)THEN
    CALL Helac_splitch(grid1)
    grid1%idat = 0
    update = .TRUE.
ELSEIF(grid1%idat.EQ.grid1%ndat/2)THEN
    CALL Helac_mergech(grid1)
    update = .TRUE.
ENDIF

! Activate the updated bins for gnrt/wght
IF(update)THEN
    DO jj=1,grid1%nch
      grid1%sch(jj) = grid1%sch(jj-1) + grid1%wch(jj)        ! entropy
!     grid%sch(jj) = grid%sch(jj-1) + dsqrt(grid%wch(jj)) ! variance
    ENDDO
    DO jj=1,grid1%nch
      grid1%sch(jj) = grid1%sch(jj)/grid1%sch(grid1%nch)
    ENDDO
    IF(nunit.GT.0)WRITE(nunit,*) &
      'MESSAGE from MC_Helac_GRID ',grid1%lbl,': nch =',grid1%nch
ENDIF
END SUBROUTINE Helac_grid_collect


SUBROUTINE Helac_splitch( grid1 )
!********************************************************************
! split channel
! Eq(10) in arXiv:0710.2448v2[hep-ph]
!********************************************************************
IMPLICIT NONE
TYPE(Helac_grid_type),INTENT(INOUT) :: grid1
REAL(KIND(1d0)):: ineff,ineff_old,wchtmp,wchmax,xx
INTEGER:: nch,ii,nsame
INTEGER,DIMENSION(grid1%nchmax)::lch
REAL(KIND(1d0)),PARAMETER:: ex=0.0d0

ineff_old = 1d300
DO WHILE(grid1%nch.LT.grid1%nchmax)
    nch = grid1%nch
    wchmax = -1d0
    DO ii=1,nch
      wchtmp = grid1%wch(ii) !**(1d0-ex) * (grid%xri(ii)-grid%xri(ii-1))**ex
      IF(wchtmp.GT.wchmax)THEN
        wchmax = wchtmp
        lch(1) = ii
        nsame = 1
      ELSEIF(wchtmp.EQ.wchmax)THEN
        nsame = nsame+1
        lch(nsame) = ii
      ENDIF
    ENDDO

    ineff = wchmax*nch
    IF(ineff.GE.ineff_old)THEN
      EXIT
    ELSE
      ineff_old = ineff
    ENDIF

    IF(nsame.GT.1)THEN
      CALL Helac_grid_random( xx )
      ii = lch(1+INT(nsame*xx))
    ELSE
      ii = lch(1)
    ENDIF

    IF(grid1%xri(ii)-grid1%xri(ii-1).LE.EPSILON(1d0))EXIT 
    grid1%xri(ii+1:nch+1) = grid1%xri(ii:nch)
    grid1%wch(ii+1:nch+1) = grid1%wch(ii:nch)
    grid1%nch = nch+1
    grid1%xri(ii)   = ( grid1%xri(ii-1)+grid1%xri(ii) )/2
    grid1%wch(ii)   = grid1%wch(ii)/2
    grid1%wch(ii+1) = grid1%wch(ii)
ENDDO

END SUBROUTINE Helac_splitch


SUBROUTINE Helac_mergech( grid1 )
!********************************************************************
! merge channels
!********************************************************************
IMPLICIT NONE
TYPE(Helac_grid_type),INTENT(INOUT) :: grid1
REAL(KIND(1d0)) :: wchtmp,wchmin,xx
INTEGER:: nch,ii,nsame,newnch
INTEGER,DIMENSION(grid1%nchmax)::lch
REAL(KIND(1d0)),PARAMETER:: ex=0.0d0

newnch = 1+INT( (1d0-grid1%frac)*DBLE(grid1%nch) )
DO WHILE(grid1%nch.GT.newnch)
    nch = grid1%nch
    wchmin = 1d300
    DO ii=1,nch-1
!      wchtmp = grid1%wch(ii)+grid1%wch(ii+1)
      wchtmp = grid1%xri(ii+1)-grid1%xri(ii-1)
      IF(wchtmp.LT.wchmin)THEN
        wchmin = wchtmp
        lch(1) = ii
        nsame = 1
      ELSEIF(wchtmp.EQ.wchmin)THEN
        nsame = nsame+1
        lch(nsame) = ii
      ENDIF
    ENDDO

    IF(nsame.GT.1)THEN
      CALL Helac_grid_random( xx )
      ii = lch(1+INT(nsame*xx))
    ELSE
      ii = lch(1)
    ENDIF

    grid1%xri(ii) = grid1%xri(ii+1)
    grid1%wch(ii) = grid1%wch(ii)+grid1%wch(ii+1)
    grid1%xri(ii+1:nch-1) = grid1%xri(ii+2:nch)
    grid1%wch(ii+1:nch-1) = grid1%wch(ii+2:nch)
    grid1%nch = nch-1
ENDDO

END SUBROUTINE Helac_mergech


SUBROUTINE Helac_findch( grid1 ,ii ,xx )
!********************************************************************
!********************************************************************
IMPLICIT NONE
TYPE(Helac_grid_type),INTENT(IN):: grid1
INTEGER,INTENT(OUT):: ii
REAL(KIND(1d0)),INTENT(IN):: xx
INTEGER::i0,i1

IF(xx.LT.0d0.OR.xx.GT.1d0)THEN
    ii = 0
    RETURN
ENDIF

i0 = 0
i1 = grid1%nch
DO WHILE(i1-i0.GT.1)
    ii = (i0+i1)/2
    IF(xx.LE.grid1%xri(ii))THEN
      i1 = ii
    ELSE
      i0 = ii
    ENDIF
ENDDO
ii = i1
END SUBROUTINE Helac_findch


SUBROUTINE Helac_grid_plot( grid1  ,iunit )
!********************************************************************
!********************************************************************
IMPLICIT NONe
TYPE(Helac_grid_type),INTENT(IN) :: grid1
INTEGER,INTENT(IN) :: iunit
CHARACTER(4+nlbl) :: filename
REAL(KIND(1d0)),PARAMETER::o0 = 0d0
REAL(KIND(1d0)):: h1,h2
INTEGER :: ii
IF(iunit.LE.0)RETURN
IF(grid1%nchmax.LE.1)RETURN
filename = 'grid'//grid1%lbl
OPEN(UNIT=iunit,FILE=TRIM(tmp_dir)//filename,STATUS='REPLACE' )
DO ii=1,grid1%nch
   h1 = ( grid1%sch(ii) - grid1%sch(ii-1) ) &
       / ( grid1%xri(ii) - grid1%xri(ii-1) )
   h2 = grid1%sch(ii)
   WRITE(iunit,'(99e16.8)') grid1%xri(ii-1),o0,o0
   WRITE(iunit,'(99e16.8)') grid1%xri(ii-1),h1,h2
   WRITE(iunit,'(99e16.8)') grid1%xri(ii  ),h1,h2
   WRITE(iunit,'(99e16.8)') grid1%xri(ii  ),o0,o0
ENDDO
CLOSE(UNIT=iunit )
END SUBROUTINE Helac_grid_plot
 
!**********************************************************************
!**********************************************************************
!*
!* Now follow the routines to be called in Helac
!*
!**********************************************************************
!**********************************************************************


SUBROUTINE Helac_grid_alloc_2( nsize_in )
!********************************************************************
!********************************************************************
IMPLICIT NONE
INTEGER,INTENT(IN) :: nsize_in
nsize = nsize_in
IF(ALLOCATED(grid))THEN
	DEALLOCATE(grid)
ENDIF
IF(ALLOCATED(active))THEN
	DEALLOCATE(active)
ENDIF
ALLOCATE( grid(1:nsize) )
ALLOCATE( active(1:nsize) )
active(1:nsize) = .FALSE.
END SUBROUTINE Helac_grid_alloc_2


SUBROUTINE Helac_grid_init_2( ID ,nbatch,nchmax,frac ,igraph,iinv )
!********************************************************************
!* "igraph" and "iinv" are just used to give a character label to
!* the grid. If "helac_grid_plot_2" is called, the grid will be put
!* in the file "grid<igraph>i<iinv>"
!********************************************************************
IMPLICIT NONE
INTEGER,INTENT(IN) :: ID,nbatch,nchmax,igraph,iinv
REAL(KIND(1d0)),INTENT(IN)::frac
CHARACTER(8)::aa
IF(ID.GT.nsize)THEN
    WRITE(*,*) 'ERROR in Helac_grid_init_2: ID=',ID,', while nsize=',nsize &
              ,'. You have to allocate at a higher value of nsize'
    STOP
ENDIF
active(ID) = .TRUE.
aa = '0000i000'
WRITE(aa(1:4),'(i4.4)') igraph
WRITE(aa(6:8),'(i3.3)') iinv
CALL Helac_grid_init( grid(ID) ,nbatch,nchmax,frac ,aa )
END SUBROUTINE Helac_grid_init_2


SUBROUTINE Helac_grid_gnrt_2( ID ,xx )
!********************************************************************
!********************************************************************
IMPLICIT NONE
INTEGER,INTENT(IN)  :: ID
REAL(KIND(1d0)),INTENT(OUT) :: xx
IF(active(ID))THEN
    CALL Helac_grid_gnrt( grid(ID) ,xx )
ELSE
    WRITE(*,*) 'ERROR in Helac_grid_gnrt_2: grid(',ID,') was not initialized'
    STOP
ENDIF
END SUBROUTINE

  
SUBROUTINE Helac_grid_wght_2( ID ,ww ,xx )
!********************************************************************
!********************************************************************
IMPLICIT NONE
INTEGER,INTENT(IN)  :: ID
REAL(KIND(1d0)),INTENT(OUT) :: ww
REAL(KIND(1d0)),INTENT(IN) :: xx
IF(active(ID))THEN
    CALL Helac_grid_wght( grid(ID) ,ww ,xx )
ELSE
    WRITE(*,*) 'ERROR in Helac_grid_wght_2: grid(',ID,') was not initialized'
    STOP
ENDIF
END SUBROUTINE Helac_grid_wght_2


SUBROUTINE Helac_grid_collect_2( ww )
!********************************************************************
!********************************************************************
IMPLICIT NONE
REAL(KIND(1d0)),INTENT(IN) :: ww
INTEGER:: ii
DO ii=1,nsize
   IF(active(ii)) CALL Helac_grid_collect( grid(ii) ,ww ,-1d0 )
ENDDO
END SUBROUTINE Helac_grid_collect_2

  
SUBROUTINE Helac_grid_plot_2( iunit )
!********************************************************************
!********************************************************************
IMPLICIT NONE
INTEGER,INTENT(IN) :: iunit
INTEGER:: ii
IF(iunit.LE.0)RETURN
DO ii=1,nsize
   IF(active(ii))CALL Helac_grid_plot( grid(ii) ,iunit )
ENDDO
END SUBROUTINE Helac_grid_plot_2


SUBROUTINE Helac_grid_close_2
!********************************************************************
!********************************************************************
IMPLICIT NONE
INTEGER:: ii
DO ii=1,nsize
   IF(active(ii))CALL Helac_grid_close( grid(ii) )
ENDDO
IF(ALLOCATED(grid))DEALLOCATE( grid )
IF(ALLOCATED(active))DEALLOCATE( active )
nsize = 0
END SUBROUTINE Helac_grid_close_2


SUBROUTINE Helac_grid_remove_2( ii )
!********************************************************************
!********************************************************************
IMPLICIT NONE
INTEGER,INTENT(IN) :: ii
IF(active(ii))THEN
    active(ii)=.FALSE.
    CALL Helac_grid_close( grid(ii) )
ENDIF
END SUBROUTINE Helac_grid_remove_2
END MODULE MC_Helac_GRID
