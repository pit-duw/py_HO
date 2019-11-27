MODULE FO_plot
  ! Wrapper routines for the fixed order analyses
  USE HELAC_Global
  USE plot_user
  USE Helac_Func_1
  IMPLICIT NONE
  INTEGER::max_weight=1 ! max_weight=maxscales*maxscales+maxpdfs+1
  SAVE
  !LOGICAL::do_rwgt_scale = .FALSE. ! now don't do reweight for scales
  !LOGICAL::do_rwgt_pdf = .FALSE. ! now don't do reweight for pdf
  !INTEGER::numscales=3
  !REAL(KIND(1d0))::rw_Fscale_up=1d0,rw_Fscale_down=1d0
  !REAL(KIND(1d0))::rw_Rscale_up=1d0,rw_Rscale_down=1d0
  !INTEGER::pdf_set_max,pdf_set_min,npdfs=1
  !INTERFACE outfun
  !   MODULE PROCEDURE outfun_c
  !   MODULE PROCEDURE outfun_u
  !END INTERFACE outfun
CONTAINS
  SUBROUTINE initplot
    IMPLICIT NONE
!    include 'run.inc'
!    include 'reweight0.inc'
    INTEGER::nwgt
    !INTEGER,PARAMETER::max_weight=3 ! max_weight=maxscales*maxscales+maxpdfs+1
    CHARACTER(len=15),DIMENSION(:),ALLOCATABLE::weights_info
    INTEGER::i
    max_weight=1+ho_nscale+ho_npdf
    IF(ALLOCATED(weights_info))THEN
       DEALLOCATE(weights_info)
    ENDIF
    ALLOCATE(weights_info(max_weight))
    nwgt=1
    weights_info(nwgt)="central value  "
    IF(reweight_scale)THEN
       nwgt=nwgt+ho_nscale
       IF(ho_nscale.NE.8)THEN
          WRITE(*,*) 'ERROR #1 in initplot:',ho_nscale
          STOP
       ENDIF
       WRITE(weights_info(nwgt-7),'(a4,f3.1,x,a4,f3.1)')&
            "muR=",1.0,"muF=",rw_Fscale_up
       WRITE(weights_info(nwgt-6),'(a4,f3.1,x,a4,f3.1)')&
            "muR=",1.0,"muF=",rw_Fscale_down
       WRITE(weights_info(nwgt-5),'(a4,f3.1,x,a4,f3.1)')&
            "muR=",rw_Rscale_up,"muF=",1d0
       WRITE(weights_info(nwgt-4),'(a4,f3.1,x,a4,f3.1)')&
            "muR=",rw_Rscale_up,"muF=",rw_Fscale_up
       WRITE(weights_info(nwgt-3),'(a4,f3.1,x,a4,f3.1)')&
            "muR=",rw_Rscale_up,"muF=",rw_Fscale_down
       WRITE(weights_info(nwgt-2),'(a4,f3.1,x,a4,f3.1)')&
            "muR=",rw_Rscale_down,"muF=",1d0
       WRITE(weights_info(nwgt-1),'(a4,f3.1,x,a4,f3.1)')&
            "muR=",rw_Rscale_down,"muF=",rw_Fscale_up
       WRITE(weights_info(nwgt  ),'(a4,f3.1,x,a4,f3.1)')&
            "muR=",rw_Rscale_down,"muF=",rw_Fscale_down
    END IF
    IF(reweight_pdf)THEN
       DO i=1,ho_npdf
          WRITE(weights_info(nwgt+i),'(a4,i8,a3)')&
               'PDF=',ipdf_set_min-1+i,'   '
       ENDDO
       nwgt=nwgt+ho_npdf
    ENDIF
    ! output plot files
    CALL plot_begin(nwgt,max_weight,weights_info)
    RETURN
  END SUBROUTINE initplot


  SUBROUTINE plotout
    USE MC_VEGAS
  ! output of plot files
    IMPLICIT NONE
    LOGICAL::usexinteg=.FALSE.,mint=.FALSE. ! I have not implemented MINT
!    common/cusexinteg/usexinteg,mint
    !INTEGER::itmax=1,ncall
!    common/citmax/itmax,ncall
    LOGICAL::useitmax=.TRUE.
!    common/cuseitmax/useitmax
    REAL(KIND(1d0))::xnorm

    IF(usexinteg.AND..NOT.mint.AND.gener.EQ.3)THEN
       xnorm=1.d0/float(itmx)
    ELSEIF(mint)THEN
       xnorm=1.d0/float(ntotps)
    ELSEIF(gener.NE.3)THEN ! not vegas
       xnorm=1.d0/float(ntotps)
    ELSE
       xnorm=1d0
    ENDIF
    IF(useitmax.AND.gener.EQ.3)xnorm=xnorm/float(itmx) ! use vegas with itmax
    ! output plot files
    CALL plot_end(xnorm)
    RETURN                
  END SUBROUTINE plotout


  SUBROUTINE outfun(www)
    IMPLICIT NONE
    !include 'run.inc'
    !include 'genps.inc'
    !include 'reweight.inc'
    !include 'reweightNLO.inc'
    REAL(KIND(1d0)),INTENT(IN)::www
    INTEGER::i,j
    !REAL(KIND(1d0)),DIMENSION(3)::xd=(/0d0,0d0,1d0/)
    !INTEGER,PARAMETER::maxflow=999
    INTEGER::nwgt
    REAL(KIND(1d0))::alphasc,alphasu,alphasd
    REAL(KIND(1d0)),DIMENSION(:),ALLOCATABLE::wgts !,wgtden
    INTEGER::init=0
    SAVE init,wgts
    ! for the future use of estimating scale/pdf uncertainties
    IF(init.EQ.0)THEN
       max_weight=1+ho_nscale+ho_npdf
       IF(ALLOCATED(wgts))THEN
          DEALLOCATE(wgts)
       ENDIF
       ALLOCATE(wgts(max_weight))
       init=1
    ENDIF
    nwgt=1
    wgts(1)=www
    IF(reweight_scale)THEN
       DO i=1,3
          DO j=1,3
             IF(i.EQ.1.AND.j.EQ.1)CYCLE
             nwgt=nwgt+1
             wgts(nwgt)=www*wgtxsecmu(i,j)
          ENDDO
       ENDDO
    ENDIF
    IF(reweight_pdf)THEN
       DO i=1,ho_npdf
          nwgt=nwgt+1
          wgts(nwgt)=www*wgtxsecpdf(i)
       ENDDO
    ENDIF
    ! output plot file
    CALL plot_fill(wgts)
    RETURN
  END SUBROUTINE outfun
  
END MODULE FO_plot
