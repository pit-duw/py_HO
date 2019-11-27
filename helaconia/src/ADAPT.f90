MODULE Adapt
USE Helac_ranmar_mod
USE Helac_Global
IMPLICIT NONE
INTEGER,PARAMETER::maxn=100,maxv=100
REAL(KIND(1d0)),DIMENSION(maxn,maxv)::alfa_adapt,alfao_adapt,low_adapt,upp_adapt,wei_adapt,we0_adapt
REAL(KIND(1d0)),DIMENSION(0:maxn,maxv)::beta_adapt
REAL(KIND(1d0)),DIMENSION(maxv)::lo_adapt,up_adapt
REAL(KIND(1d0))::ex_adapt,ey_adapt
INTEGER::ii,kount_g,kount_p,maxadap_adapt,iadap_adapt
INTEGER,DIMENSION(maxv)::ich_adapt,init_adapt=(/(0,ii=1,maxv)/),nch_adapt,kount0_adapt,kount1_adapt
SAVE
CONTAINS

SUBROUTINE adapt_gen(nchi,lo1,up1,x,iv)
IMPLICIT NONE
INTEGER,INTENT(IN)::nchi
INTEGER,INTENT(IN)::iv
REAL(KIND(1d0)),INTENT(IN)::lo1,up1
INTEGER::i
REAL(KIND(1d0))::r
REAL(KIND(1d0)),INTENT(OUT)::x
nch_adapt(iv)=nchi
lo_adapt(iv)=lo1
up_adapt(iv)=up1
IF(GLOBALINIT_adapt_gen.EQ.0)init_adapt(1:maxv)=0

IF(init_adapt(iv).EQ.0)THEN
![1
	kount_g=0
	kount_p=0

	wei_adapt(1:maxn,1:maxv)=0
	we0_adapt(1:maxn,1:maxv)=0
	maxadap_adapt=MaxAdapt
	iadap_adapt=0
	kount0_adapt(iv)=250 !50*nch(iv)
	kount1_adapt(iv)=250 !50*nch(iv)
	ex_adapt=0.1d0
	ey_adapt=1d0

	beta_adapt(0,iv)=0
	DO i=1,nch_adapt(iv)
		alfa_adapt(i,iv)=1d0/nch_adapt(iv)     
		beta_adapt(i,iv)=beta_adapt(i-1,iv)+alfa_adapt(i,iv)
		low_adapt(i,iv)=lo_adapt(iv)+(up_adapt(iv)-lo_adapt(iv))/nch_adapt(iv)*(i-1)
		upp_adapt(i,iv)=lo_adapt(iv)+(up_adapt(iv)-lo_adapt(iv))/nch_adapt(iv)*i
	ENDDO
	init_adapt(iv)=1
!]1  
ENDIF

r=Helac_rnmy(0)
DO i=1,nch_adapt(iv)
	IF(r.LT.beta_adapt(i,iv))EXIT  ! goto 1  ! the probability changed by alfa
ENDDO
! 1
ich_adapt(iv)=i

r=Helac_rnmy(0)
x=(1-r)*low_adapt(ich_adapt(iv),iv)+r*upp_adapt(ich_adapt(iv),iv) 
IF(iv.EQ.1)kount_g=kount_g+1

END SUBROUTINE adapt_gen

SUBROUTINE adapt_weo(wo,iv)
IMPLICIT NONE
INTEGER,INTENT(IN)::iv
REAL(KIND(1d0)),INTENT(OUT)::wo      
wo=(upp_adapt(ich_adapt(iv),iv)-low_adapt(ich_adapt(iv),iv))/alfa_adapt(ich_adapt(iv),iv)
END SUBROUTINE adapt_weo

SUBROUTINE adapt_wei(w,iv,niv,lgo)
IMPLICIT NONE
REAL(KIND(1d0)),INTENT(IN)::w
INTEGER,INTENT(IN)::iv,niv
LOGICAL,INTENT(IN)::lgo
INTEGER::kcrit,i
REAL(KIND(1d0))::ww,alfamin,tota,ave,avem,da
wei_adapt(ich_adapt(iv),iv)=wei_adapt(ich_adapt(iv),iv)+w/alfa_adapt(ich_adapt(iv),iv)
! ww=w**2->w in SingProd.f90 ww/alfa=w**2/alfa is the same with Eq.(6,7) in hep-ph/9405257  
we0_adapt(1,iv)=we0_adapt(1,iv)+1d0
IF(iv.EQ.1.AND.w.GT.0)THEN
	kount_p=kount_p+1
ENDIF

kcrit=MAX(kount_p,1)
! iadap_adapt means the number of adapt steps
IF(iadap_adapt.LT.maxadap_adapt.AND.lgo)THEN
![3
	alfao_adapt(1:nch_adapt(iv),iv)=alfa_adapt(1:nch_adapt(iv),iv)
	IF(iv.EQ.niv)iadap_adapt=iadap_adapt+1
	kount1_adapt(iv)=kount1_adapt(iv) ! *1.2
	kount0_adapt(iv)=kount0_adapt(iv)+kount1_adapt(iv)

	DO i=1,nch_adapt(iv)
		ww=wei_adapt(i,iv)/we0_adapt(1,iv)
		alfa_adapt(i,iv)= alfa_adapt(i,iv)**ey_adapt * ww**ex_adapt  ! optimize with multichannel
	ENDDO

	alfamin=1d-3
	tota=0
	DO i=1,nch_adapt(iv)
		IF(alfa_adapt(i,iv).LT.alfamin)alfa_adapt(i,iv)=alfamin
		tota=tota+alfa_adapt(i,iv)
	ENDDO

	DO i=1,nch_adapt(iv)
		alfa_adapt(i,iv)=alfa_adapt(i,iv)/tota
		beta_adapt(i,iv)=beta_adapt(i-1,iv)+alfa_adapt(i,iv)
	ENDDO

	ave=0
	avem=0
	DO i=1,nch_adapt(iv)
		da=0
		IF(alfao_adapt(i,iv).GT.0)da=DABS(alfa_adapt(i,iv)-alfao_adapt(i,iv))/alfao_adapt(i,iv)
        ave=ave+da/nch_adapt(iv)
        avem=MAX(avem,da)
	ENDDO
	PRINT *,'ADAPT.F90 params: ',iv,ave,avem,iadap_adapt,kount0_adapt(iv)

	wei_adapt(1:nch_adapt(iv),iv)=0
	we0_adapt(1:nch_adapt(iv),iv)=0
!]3
ENDIF
END SUBROUTINE adapt_wei
END MODULE Adapt
