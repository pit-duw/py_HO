MODULE Phegas_Choice       
IMPLICIT NONE
INTEGER,PRIVATE::i1_c,m1_c,i2_c,m2_c,i3_c,m3_c
SAVE i1_c,m1_c,i2_c,m2_c,i3_c,m3_c
CONTAINS
SUBROUTINE init_choice(i1i,m1i,i2i,m2i,i3i,m3i)
INTEGER,INTENT(IN)::i1i,m1i,i2i,m2i,i3i,m3i
i1_c=i1i
m1_c=m1i
i2_c=i2i
m2_c=m2i
i3_c=i3i
m3_c=m3i
END SUBROUTINE init_choice

SUBROUTINE fini_choice(i1i,m1i,i2i,m2i,i3i,m3i)
INTEGER,INTENT(OUT)::i1i,m1i,i2i,m2i,i3i,m3i
i1i=i1_c
m1i=m1_c
i2i=i2_c
m2i=m2_c
i3i=i3_c
m3i=m3_c
END SUBROUTINE fini_choice

SUBROUTINE symm_choice(icase,j1,n1,j2,n2,j3,n3)
INTEGER,INTENT(IN)::icase
INTEGER,INTENT(OUT)::j1,n1,j2,n2,j3,n3
IF(i3_c.EQ.0)THEN
   IF(icase.EQ.1)THEN
        j1=i1_c
        n1=m1_c
        j2=i2_c
        n2=m2_c
   ELSEIF(icase.EQ.2)THEN
        j1=i2_c
        n1=m2_c
        j2=i1_c
        n2=m1_c
   ELSE
        WRITE(*,*)'something is wrong in Phegas_Choice',icase
        STOP
   ENDIF
ELSE
   IF(icase.EQ.1)THEN
        j1=i1_c
        n1=m1_c
        j2=i2_c
        n2=m2_c
        j3=i3_c
        n3=m3_c
   ELSEIF(icase.EQ.2)THEN
        j1=i1_c
        n1=m1_c
        j2=i3_c
        n2=m3_c
        j3=i2_c
        n3=m2_c
   ELSEIF(icase.eq.3)THEN
        j1=i2_c
        n1=m2_c
        j2=i1_c
        n2=m1_c
        j3=i3_c
        n3=m3_c
   ELSEIF(icase.EQ.4)THEN
        j1=i2_c
        n1=m2_c
        j2=i3_c
        n2=m3_c
        j3=i1_c
        n3=m1_c
   ELSEIF(icase.EQ.5)THEN
        j1=i3_c
        n1=m3_c
        j2=i1_c
        n2=m1_c
        j3=i2_c
        n3=m2_c
   ELSEIF(icase.EQ.6)THEN
        j1=i3_c
        n1=m3_c
        j2=i2_c
        n2=m2_c
        j3=i1_c
        n3=m1_c
   ELSE
        WRITE(*,*)'something is wrong in Phegas_Choice',icase
        STOP
   ENDIF
ENDIF
END SUBROUTINE symm_choice
END MODULE Phegas_Choice