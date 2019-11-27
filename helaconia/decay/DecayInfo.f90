MODULE DecayInfo
  USE Helac_Global
  IMPLICIT NONE
CONTAINS
  SUBROUTINE ReadDecayInfo
    IMPLICIT NONE
    LOGICAL::lexist,needed
    INTEGER::iounit,iostate
    CHARACTER(100)::ctmp
    INTEGER::k,i,jjj,lll,kkk
    REAL(KIND(1d0))::BR
    INTEGER,DIMENSION(0:5)::decay
    INTEGER::udecay=301127
    INQUIRE(FILE=TRIM(input_dir)//"decay_user.inp",EXIST=lexist)
    IF(.NOT.lexist)THEN
       PRINT *,"ERROR: the file decay_user.inp does not exist ! STOP !"
       STOP
    ENDIF
    INQUIRE(FILE=TRIM(input_dir)//"decay_user.inp",OPENED=lexist)
    IF(lexist)THEN
       INQUIRE(FILE=TRIM(input_dir)//"decay_user.inp",NUMBER=iounit)
       IF(iounit.NE.udecay)THEN
          PRINT *,"WARNING: the decay_user.inp has been linked with another unit ! Close and reopen !"
          CLOSE(UNIT=iounit)
          OPEN(UNIT=udecay,FILE=TRIM(input_dir)//"decay_user.inp")
       ENDIF
    ELSE
       OPEN(UNIT=udecay,FILE=TRIM(input_dir)//"decay_user.inp")
    ENDIF
    iostate=0
    jjj=0
    kkk=0
    DecayChains(1:MAX_DecayChain,-1:5)=0
    DecayBR(1:MAX_DecayChain)=0d0
    iflh2DecayChains(1:20,0:MAX_DecayChain)=0
    NDecayChains=0
    DO WHILE(iostate.NE.-1.AND.iostate.LE.0)
       READ(udecay,"(A)",IOSTAT=iostate)ctmp
       IF(iostate.LT.0)EXIT
       IF(ctmp(1:1)=="#")CYCLE
       IF(ctmp(1:11)=="Decay Chain")THEN
          READ(udecay,*,IOSTAT=iostate)k,BR
          IF(.NOT.(iostate.NE.-1.AND.iostate.LE.0))THEN
             PRINT *,"ERROR:Wrong format (number,BR) in decay_user.inp"
             STOP
          ENDIF
          IF(k.GT.5)THEN
             PRINT *,"ERROR:HELAC-Onia only supports up to 5-body decays"
             STOP
          ENDIF
          READ(udecay,*,IOSTAT=iostate)(decay(i),i=0,k)
          IF(.NOT.(iostate.NE.-1.AND.iostate.LE.0))THEN
             PRINT *,"ERROR:Wrong format (decay chain) in decay_user.inp"
             STOP
          ENDIF
          IF(BR.LE.0d0)CYCLE
          ! only select the needed ones
          needed=.FALSE.
          DO i=1,nhad
             IF(decay(0).EQ.iflh(i))THEN
                needed=.TRUE.
                lll=iflh2DecayChains(i,0)
                iflh2DecayChains(i,0)=lll+1
                iflh2DecayChains(i,lll+1)=jjj+1
             ENDIF
          ENDDO
          IF(needed)kkk=kkk+1
          jjj=jjj+1
          IF(jjj.GT.MAX_DecayChain)THEN
             PRINT *,"ERROR:Find too many decay chains !"
             STOP
          ENDIF
          DecayBR(jjj)=BR
          DecayChains(jjj,0:k)=decay(0:k)
          DecayChains(jjj,-1)=k
       ENDIF
    END DO
    NDecayChains=kkk
    AllNDecayChains=jjj
    NDecayIflh=0
    DO i=1,nhad
       jjj=iflh2DecayChains(i,0)
       IF(jjj.GT.1)NDecayIflh=NDecayIflh+1
    ENDDO
    CLOSE(UNIT=udecay)
    !    CALL SelectDecayInfo
    RETURN
  END SUBROUTINE ReadDecayInfo
  
END MODULE DecayInfo
