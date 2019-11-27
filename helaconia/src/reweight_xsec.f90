MODULE reweight_xsec
  USE Helac_Global
  USE Helac_Func_1
  USE Structf_PDFs
  USE setscale
  IMPLICIT NONE
CONTAINS
  SUBROUTINE reweight_xsec_scale(wsf0)
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::wsf0
    REAL(KIND(1d0))::alphasc=1d0,alphasu=1d0,alphasd=1d0
    REAL(KIND(1d0)),DIMENSION(3)::R_alphas,R_pdf
    REAL(KIND(1d0))::wsfw
    INTEGER::i,j
    IF(.NOT.reweight_scale)THEN
       WRITE(*,*)"ERROR:Please set reweight_scale to be T first !"
       STOP
    ENDIF
    ! reweighting for alpha_s
    IF(alphas_power.EQ.0)THEN
       alphasc=1d0
       alphasu=1d0
       alphasd=1d0
    ELSE
       CALL qcdscale(scale)
       alphasc=ALPHAS(scale)
       alphasu=ALPHAS(scale*rw_Rscale_up)
       alphasd=ALPHAS(scale*rw_Rscale_down)
    ENDIF
    R_alphas(1)=1d0
    R_alphas(2)=(alphasu/alphasc)**alphas_power
    R_alphas(3)=(alphasd/alphasc)**alphas_power
    R_pdf(1)=1d0
    IF(istruc.EQ.1)THEN
       reweight_Fscale_phase=1
       CALL strf_pdf(wsfw)
       R_pdf(2)=wsfw/wsf0
       reweight_Fscale_phase=2
       CALL strf_pdf(wsfw)
       R_pdf(3)=wsfw/wsf0
       ! recover it
       reweight_Fscale_phase=0
    ELSE
       R_pdf(2:3)=1d0
    ENDIF
    wgtxsecmu(1:3,1:3)=1d0
    DO i=1,3
       wgtxsecmu(i,1:3)=R_alphas(i)
       DO j=1,3
          wgtxsecmu(i,j)=wgtxsecmu(i,j)*R_pdf(j)
       ENDDO
    END DO
    RETURN
  END SUBROUTINE reweight_xsec_scale

  SUBROUTINE reweight_xsec_pdf(wsf0)
    IMPLICIT NONE
    INTEGER::init=0,i
    REAL(KIND(1d0)),INTENT(IN)::wsf0
    REAL(KIND(1d0))::wsfw
    SAVE init
    IF(.NOT.reweight_pdf)THEN
       WRITE(*,*)"ERROR:Please set reweight_pdf to be T first !"
       STOP
    ENDIF
    IF(init.EQ.0)THEN
       IF(MOD(ho_npdf,2).NE.0)THEN
          WRITE(*,*)'The number of error sets must be even',ho_npdf
          STOP
       ENDIF
       IF(ALLOCATED(wgtxsecpdf))THEN
          DEALLOCATE(wgtxsecpdf)
       ENDIF
       ALLOCATE(wgtxsecpdf(ho_npdf))
       IF(ALLOCATED(idpdf))THEN
          DEALLOCATE(idpdf)
       ENDIF
       ALLOCATE(idpdf(0:ho_npdf))
       idpdf(0)=iPDFSUP1
       idpdf(1)=ipdf_set_min
       IF(ipdf_set_max-ipdf_set_min+1.NE.ho_npdf)THEN
          WRITE(*,*)'Incorrect the number of error pdf',ho_npdf,ipdf_set_min,ipdf_set_max
          STOP
       ENDIF
       DO i=2,ho_npdf
          idpdf(i)=idpdf(1)+i-1
       ENDDO
       WRITE(*,*) 'Doing PDF reweight:'
       WRITE(*,*) 'Central set id: ',idpdf(0)
       WRITE(*,*) 'Min error set id: ',idpdf(1)
       WRITE(*,*) 'Max error set id: ',idpdf(ho_npdf)
       init=1
    ENDIF
    DO i=1,ho_npdf
       INCLUDE "../input/lhapdf/call_initpdf"
       CALL strf_pdf(wsfw)
       wgtxsecpdf(i)=wsfw/wsf0
    ENDDO
    ! recover it
    INCLUDE "../input/lhapdf/call_initpdf0"
    RETURN
  END SUBROUTINE reweight_xsec_pdf
END MODULE reweight_xsec
