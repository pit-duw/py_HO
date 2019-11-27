MODULE plot_fit_cb
  IMPLICIT NONE
  INTEGER,PARAMETER::MXBIN=2000,MXSET=20,MXFILE=20
  REAL(KIND(1d0)),DIMENSION(MXSET,MXBIN)::FIT_TH_HIST,FIT_TH_ERR,FIT_EX_HIST,FIT_EX_ERR
  REAL(KIND(1d0)),DIMENSION(MXSET,MXBIN)::FIT_XHIS
  REAL(KIND(1d0)),DIMENSION(MXSET)::FIT_HMAX,FIT_HMIN
  INTEGER,DIMENSION(MXSET)::FIT_NBIN
  INTEGER::FIT_NMAX
  CHARACTER*100,DIMENSION(MXSET)::FIT_LABEL
CONTAINS
  SUBROUTINE FGNUPLOT(N,itype)
    IMPLICIT NONE
    REAL(KIND(1d0))::xpos,ypos,ymin,tmp
    INTEGER,INTENT(IN)::N,itype
    CHARACTER*100::psfile
    SAVE psfile
    INTEGER::J
    CHARACTER*15::title
    CHARACTER*18::datafile
    IF(N.EQ.1)THEN
       psfile="fit_results.ps"
       OPEN(UNIT=9797,FILE="./output/fit_results.gnu")
       WRITE(9797,102)TRIM(psfile)
102    FORMAT(/1X,&
            ' reset',/1X,&
            ' set lmargin 10',/1X,&
            ' set rmargin 0',/1X,&
            ' set terminal postscript col enhanced',/1X,&
            ' set output ','"',A,'"',/1X,&
            ' set style data histeps',/1X,&
            ' set key off',/1X,&
            ' set style line 1 lt 1 lc rgb "#006D4F" lw 1.8',/1X,&
            ' set style line 2 lt 1 lc rgb "#B90091" lw 1.8')
       !WRITE(9797,*)' set multiplot'
       !WRITE(9797,*)' set style data histeps'
    ELSE
       WRITE(9797,*)' unset label'
    ENDIF
    xpos=FIT_HMAX(N)+(FIT_HMAX(N)-FIT_HMIN(N))*0.05d0/2.2d0
    ypos=0d0
    ymin=1d99
    DO J=1,FIT_NBIN(N)
       tmp=MAX(FIT_TH_HIST(N,J),FIT_EX_HIST(N,J))
       IF(ypos.LT.tmp)ypos=tmp
    ENDDO
    WRITE(title, "(A13,I2)") "SET NUMBER = ",N
    SELECT CASE(itype)
    CASE(2)
       WRITE(9797,101)TRIM(FIT_LABEL(N)),title,"P_T [GeV]","d^2{/Symbol s}/dP_Tdx_F",&
            FIT_HMIN(N),FIT_HMAX(N)
    CASE DEFAULT
       WRITE(9797,101)TRIM(FIT_LABEL(N)),title,"P_T [GeV]","d^2{/Symbol s}/dP_Tdy",&
            FIT_HMIN(N),FIT_HMAX(N)
    END SELECT
101 FORMAT( /1x,&
         ' set label ','"',A,'" font ",14" at graph 0.1, graph 0.94',/1X,&
         ' set title ','"',A,' distribution" font "Helvetica, 20"',/1X,&
         ' set xlabel ','"',A,'" font "Helvetica, 20"',/1X,&
         ' set ylabel ','"',A,&
         ' [nb/GeV]" font "Helvetica, 20"',/1X,&
         ' set label front ','"','HELAC-ONIA','" font "Helvetica, 15" rotate by 90 at graph 1.02, graph 0.4',/1X,&
         !' rotate by 90 at 'G13.6,',',G13.6,/1X,&
         ' set xrange [ ',F10.5,':',F10.5,']',/1X,&
         ' set format y "10^{%T}"',/1X,&
         ' set key horizontal noreverse maxcols 1 width -4 at graph 0.92, graph 0.9')
    WRITE(9797,*)' set logscale y'
    WRITE(9797,*)' plot \'
    IF(N.LE.9)THEN
       WRITE(datafile,"(A11,I1,A4)")"comparison_",N,".dat"
    ELSEIF(N.LE.99)THEN
       WRITE(datafile,"(A11,I2,A4)")"comparison_",N,".dat"
    ELSEIF(N.LE.999)THEN
       WRITE(datafile,"(A11,I3,A4)")"comparison_",N,".dat"
    ELSE
       WRITE(*,*)"ERROR:Wrong 1 !"
       STOP
    ENDIF
    WRITE(9797,*)' "./'//TRIM(datafile)//'" using 1:($4):($4+$5):($4-$5) with yerror ls 1 title "EX",\'
    WRITE(9797,*)' "./'//TRIM(datafile)//'" using 1:($2) with histeps ls 2 title "TH"'
107 FORMAT( /1x,&
         ' set title ','"',A,' distribution" font "Helvetica, 20"',/1X,&
         ' set xlabel ','"',A,'" font "Helvetica, 20"',/1X,&
         ' set ylabel ','"',A,&
         ' [nb/GeV]" font "Helvetica, 20"',/1X,&
         ' set label ','"','HELAC-ONIA','" font "Helvetica, 15"',&
         ' rotate by 90 at 'G13.6,',',G13.6,/1X,&
         ' set xrange [ ',F10.5,':',F10.5,']',/1X,&
         ' set format y "10^{%T}"')
    OPEN(UNIT=29943,FILE="./output/"//TRIM(datafile))
    DO J=1,FIT_NBIN(N)
       WRITE(29943,'(3X,G13.6,4(2X,G13.6))')FIT_XHIS(N,J),FIT_TH_HIST(N,J),FIT_TH_ERR(N,J),FIT_EX_HIST(N,J),FIT_EX_ERR(N,J)
    ENDDO
    CLOSE(UNIT=29943)
    !WRITE(9797,200)
200 FORMAT(' e')
    IF(N.EQ.FIT_NMAX)THEN
       CLOSE(UNIT=9797)
    ENDIF
  END SUBROUTINE FGNUPLOT
END MODULE plot_fit_cb
