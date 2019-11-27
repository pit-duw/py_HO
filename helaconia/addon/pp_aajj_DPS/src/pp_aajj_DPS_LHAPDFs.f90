SUBROUTINE strf_lhapdf_pp_aajj_dps(ipip,id1,id2,wsf)
  USE Helac_Global
  USE pp_aajj_dps_global
  USE pp_aajj_dps_func
  INTEGER,INTENT(IN)::ipip
  INTEGER::ipp=1  ! ipp=1 ppbar ; ipp=2 pp
  INTEGER::ih=1   ! ih=1 no photon PDF,ih=2, photon from proton/anti-proton,ih=3 photon from electron/positron 
  REAL(KIND=DBL),INTENT(OUT)::wsf          
  CHARACTER(len=20),DIMENSION(20)::parm
  REAL(KIND=DBL),DIMENSION(20)::val
  REAL(KIND=DBL),DIMENSION(-7:7)::pdflist
  INTEGER::init=0
  REAL(KIND(1d0))::xpp1,xpp2
  SAVE init,ipp,ih
  REAL(KIND=DBL)::glu1_ct,glu2_ct,u1_ct,u2_ct,d1_ct,d2_ct,s1_ct,s2_ct,c1_ct,c2_ct,b1_ct,b2_ct,&
       ub1_ct,ub2_ct,db1_ct,db2_ct,sb1_ct,sb2_ct,cb1_ct,cb2_ct,bb1_ct,bb2_ct,sf_ct
  
  IF(init.EQ.0)THEN
     IF(.NOT.lhapdfwrapinit)THEN
        CALL lhapdfwrap
        lhapdfwrapinit=.TRUE.
     ENDIF
     IF(Coll_Type.EQ.1)THEN
        ipp=2
     ELSE
        ipp=1
     ENDIF
     init=1
  ENDIF

  CALL setscale_dps(ipip,scale)
  xpp1=xpi(ipip,1)
  xpp2=xpi(ipip,2)
  IF(strf_pdf_saveq(ipip))THEN
     wsf=strf_pdf_save(ipip,1,id1)*strf_pdf_save(ipip,2,id2)/xpp1/xpp2
     IF(init.EQ.0)init=1
     RETURN
  ENDIF

  CALL pftopdglha(ih,xpp1,scale,pdflist(-7:7))
  glu1_ct = xpp1*pdflist(0)
  u1_ct   = xpp1*pdflist(2)
  d1_ct   = xpp1*pdflist(1)
  s1_ct   = xpp1*pdflist(3)
  c1_ct   = xpp1*pdflist(4)
  b1_ct   = xpp1*pdflist(5)
  ub1_ct  = xpp1*pdflist(-2)
  db1_ct  = xpp1*pdflist(-1)
  sb1_ct  = xpp1*pdflist(-3)
  cb1_ct  = xpp1*pdflist(-4)
  bb1_ct  = xpp1*pdflist(-5)
  CALL pftopdglha(ih,xpp2,scale,pdflist(-7:7))
  glu2_ct = xpp2*pdflist(0)
  u2_ct   = xpp2*pdflist(2)
  d2_ct   = xpp2*pdflist(1)
  s2_ct   = xpp2*pdflist(3)
  c2_ct   = xpp2*pdflist(4)
  b2_ct   = xpp2*pdflist(5)
  ub2_ct  = xpp2*pdflist(-2)
  db2_ct  = xpp2*pdflist(-1)
  sb2_ct  = xpp2*pdflist(-3)
  cb2_ct  = xpp2*pdflist(-4)
  bb2_ct  = xpp2*pdflist(-5)

  strf_pdf_saveq(ipip)=.TRUE.
  strf_pdf_save(ipip,1,0)=glu1_ct
  strf_pdf_save(ipip,1,1)=d1_ct
  strf_pdf_save(ipip,1,2)=u1_ct
  strf_pdf_save(ipip,1,3)=s1_ct
  strf_pdf_save(ipip,1,4)=c1_ct
  strf_pdf_save(ipip,1,5)=b1_ct
  strf_pdf_save(ipip,1,-1)=db1_ct
  strf_pdf_save(ipip,1,-2)=ub1_ct
  strf_pdf_save(ipip,1,-3)=sb1_ct
  strf_pdf_save(ipip,1,-4)=cb1_ct
  strf_pdf_save(ipip,1,-5)=bb1_ct
  IF(ipp.EQ.2)THEN
     ! p p
     strf_pdf_save(ipip,2,0)=glu2_ct
     strf_pdf_save(ipip,2,1)=d2_ct
     strf_pdf_save(ipip,2,2)=u2_ct
     strf_pdf_save(ipip,2,3)=s2_ct
     strf_pdf_save(ipip,2,4)=c2_ct
     strf_pdf_save(ipip,2,5)=b2_ct
     strf_pdf_save(ipip,2,-1)=db2_ct
     strf_pdf_save(ipip,2,-2)=ub2_ct
     strf_pdf_save(ipip,2,-3)=sb2_ct
     strf_pdf_save(ipip,2,-4)=cb2_ct
     strf_pdf_save(ipip,2,-5)=bb2_ct
  ELSE
     ! p p~
     strf_pdf_save(ipip,2,0)=glu2_ct
     strf_pdf_save(ipip,2,-1)=d2_ct
     strf_pdf_save(ipip,2,-2)=u2_ct
     strf_pdf_save(ipip,2,-3)=s2_ct
     strf_pdf_save(ipip,2,-4)=c2_ct
     strf_pdf_save(ipip,2,-5)=b2_ct
     strf_pdf_save(ipip,2,1)=db2_ct
     strf_pdf_save(ipip,2,2)=ub2_ct
     strf_pdf_save(ipip,2,3)=sb2_ct
     strf_pdf_save(ipip,2,4)=cb2_ct
     strf_pdf_save(ipip,2,5)=bb2_ct
  ENDIF
  wsf=strf_pdf_save(ipip,1,id1)*strf_pdf_save(ipip,2,id2)/xpp1/xpp2
  IF(init.EQ.0)init=1
  RETURN     
END SUBROUTINE strf_lhapdf_pp_aajj_dps
