MODULE unweight_lhe
IMPLICIT NONE
CONTAINS
SUBROUTINE Generate_lhe(n1,nevent,process,icase)
!USE Helac_Global
!USE CONSTANTS
IMPLICIT NONE
INTEGER,INTENT(IN)::n1,nevent,icase
CHARACTER(len=400),INTENT(IN)::process
CHARACTER(len=24),PARAMETER::tmp_dir="./tmp/",output_dir="./output/"
INTEGER::i,nunit4,nunit5,nunit3
INTEGER::istop,k
REAL(KIND(1d0))::p0,px,py,pz,SPINUP,EBMUP1,EBMUP2,XSECUP1,XSECUP2,XERRUP1,XMAXUP1,&
XWGTUP,VTIMUP,SCALUP,PM0,AQEDUP,AQCDUP,XWGTUP2
INTEGER::IDBMUP1,IDBMUP2,IDWTUP,NPRUP,IDPRUP,NUP,IDUP,ISTUP,IMOTHUP1,IMOTHUP2,&
ICOLUP1,ICOLUP2,iPDFGUP1,iPDFGUP2,iPDFSUP1,iPDFSUP2,LPRUP1
nunit3=30
CLOSE(nunit3)
OPEN(nunit3,FILE=TRIM(tmp_dir)//'even_'//TRIM(process(2:10*n1+1))//'.out',FORM='unformatted')
nunit4=31
CLOSE(nunit4)
OPEN(nunit4,FILE=TRIM(tmp_dir)//'sample'//TRIM(process(2:10*n1+1))//'.init')
nunit5=32
CLOSE(nunit5)
OPEN(nunit5,FILE=TRIM(output_dir)//'sample'//TRIM(process(2:10*n1+1))//'.lhe')
WRITE(nunit5,'(A)') '<LesHouchesEvents version="1.0">'
WRITE(nunit5,'(A)') '<!--'
WRITE(nunit5,'(A)') 'File generated with HELAC-ONIA '
WRITE(nunit5,'(A)') '-->'
WRITE(nunit5,'(A)') '<init>' 
READ(nunit4,*)IDBMUP1,IDBMUP2,EBMUP1,EBMUP2,iPDFGUP1,iPDFGUP2,iPDFSUP1,iPDFSUP2,IDWTUP,NPRUP
READ(nunit4,*)XSECUP1,XERRUP1,XMAXUP1,LPRUP1
WRITE(nunit5,5000)IDBMUP1,IDBMUP2,EBMUP1,EBMUP2,iPDFGUP1,iPDFGUP2,iPDFSUP1,iPDFSUP2,IDWTUP,NPRUP
WRITE(nunit5,5100)XSECUP1,XERRUP1,XMAXUP1,LPRUP1       
WRITE(nunit5,'(A)') '</init>'
istop=1
k=0
XWGTUP2=XSECUP1/nevent
DO WHILE(istop.EQ.1)
	k=k+1
	READ(nunit3,END=100) NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP
        IF(icase.EQ.1)XWGTUP=XWGTUP2
	WRITE(nunit5,'(A)') '<event>' 
	WRITE(nunit5,5200) NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP
	DO i=1,NUP
		READ(nunit3)IDUP,ISTUP,iMOTHUP1,iMOTHUP2,ICOLUP1,ICOLUP2,px,py,pz,p0,pm0,VTIMUP,SPINUP
		WRITE(nunit5,5300)IDUP,ISTUP,iMOTHUP1,iMOTHUP2,ICOLUP1,ICOLUP2,px,py,pz,p0,pm0,VTIMUP,SPINUP
	ENDDO
	WRITE(nunit5,'(A)') '</event>' 
ENDDO
100  CONTINUE
IF(icase.EQ.1.AND.k.NE.nevent+1)THEN
   WRITE(*,*)"WARNING:mismatching of the unweighted lhe events number ",k-1,nevent
ENDIF
WRITE(nunit5,'(A)') '</LesHouchesEvents>' 
CLOSE(nunit3,STATUS='delete')
CLOSE(nunit4,STATUS='delete')
5200 FORMAT(1P,2I6,4E14.6)
5300 FORMAT(1P,I8,5I5,5E18.10,E14.6,E12.4)
5000 FORMAT(1P,2I8,2E14.6,6I8)
5100 FORMAT(1P,3E20.10,I6)
CLOSE(nunit5,STATUS='keep')
END SUBROUTINE Generate_lhe

SUBROUTINE FILL_MC_MASS(MonteCarlo)
  IMPLICIT NONE
  INTEGER::i
  CHARACTER*10,INTENT(IN)::MonteCarlo
  REAL(KIND(1d0)),DIMENSION(-16:21)::MCMASS
  DO i=-16,21
     MCMASS(i)=-1d10
  ENDDO
  IF(MonteCarlo.EQ.'HERWIG6')THEN
     INCLUDE '../shower/HERWIG6/MCmasses_HERWIG6.inc'
  ELSEIF(MonteCarlo.EQ.'HERWIGPP')THEN
     INCLUDE '../shower/HERWIGPP/MCmasses_HERWIGPP.inc'
  ELSEIF(MonteCarlo.EQ.'PYTHIA6Q')THEN
     INCLUDE '../shower/PYTHIA6/MCmasses_PYTHIA6Q.inc'
  ELSEIF(MonteCarlo.EQ.'PYTHIA6PT')THEN
     INCLUDE '../shower/PYTHIA6/MCmasses_PYTHIA6PT.inc'
  ELSEIF(MonteCarlo.EQ.'PYTHIA8')THEN
     INCLUDE '../shower/PYTHIA8/MCmasses_PYTHIA8.inc'
  ELSE
     WRITE(*,*)'ERROR:Unknown MC = '//TRIM(MonteCarlo)//' in fill_MC_MASS'
     STOP
  ENDIF
  DO i=-5,-1
     MCMASS(i)=MCMASS(-i)
  ENDDO
  DO i=-16,-11
     MCMASS(i)=MCMASS(-i)
  ENDDO
  RETURN
END SUBROUTINE FILL_MC_MASS
END MODULE unweight_lhe
