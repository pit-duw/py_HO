!=================================================================================
! Package MINT
!=================================================================================
! Integrator Package for POWHEG
! Written by Paolo Nason (arXiv:0709.2085)
! Modified by Hua-Sheng Shao (23 May 2013)
!      subroutine mint(fun,ndim,ncalls0,nitmax,imode,
! ndim=number of dimensions
! ncalls0=# of calls per iteration
! nitmax =# of iterations
! fun(xx,www,ifirst): returns the function to be integrated multiplied by www;
!                     xx(1:ndim) are the variables of integration
!                     ifirst=0: normal behaviour
! imode: integer flag
!
! imode=0:
! When called with imode=0 the routine integrates the absolute value of the function
! and sets up a grid xgrid(0:50,ndim) such that in each ndim-1 dimensional slice
! (i.e. xgrid(m-1,n)<xx(n)<xgrid(m,n)) the contribution of the integral is the same
! the array xgrid is setup at this stage; ans and err are the integral and its error 
!
! imode=1 (in fact #0)
! When called with imode=1, the routine performs the integral of the function fun
! using the grid xgrid. If some number in the array ifold, (say, ifold(n))
! is different from 1, it must be a divisor of 50, and the 50 intervals xgrid(0:50,n)
! are grouped into ifold(n) groups, each group containing 50/ifold(n) nearby
! intervals. For example, if ifold(1)=5, the 50 intervals for the first dimension
! are divided in 5 groups of 10. The integral is then performed by folding on top
! of each other these 5 groups. Suppose, for example, that we choose a random point
! in xx(1) = xgrid(2,1)+x*(xgrid(3,1)-xgrid(2,1)), in the group of the first 5 interval.
! we sum the contribution of this point to the contributions of points
! xgrid(2+m*10,1)+x*(xgrid(3+m*10,1)-xgrid(2+m*10,1)), with m=1,...,4.
! In the sequence of calls to the
! function fun, the call for the first point is performed with ifirst=0, and that for
! all subsequent points with ifirst=1, so that the function can avoid to compute
! quantities that only depend upon dimensions that have ifold=1, and do not change
! in each group of folded call. The values returned by fun in a sequence of folded
! calls with ifirst=0 and ifirst=1 are not used. The function itself must accumulate
! the values, and must return them when called with ifirst=2.
! 
MODULE mint_integrator
INTEGER,PARAMETER::ndimmax=6
INTEGER,PARAMETER::nintervals=50
! specify how many times each dimension is folded
INTEGER,DIMENSION(ndimmax)::ifold
SAVE ifold
CONTAINS
  SUBROUTINE mint(fun,ndim,ncalls0,nitmax,imode,&
    xgrid,xint,ymax,ans,err)
! imode=0: integrate and adapt the grid
! imode=1: frozen grid, compute the integral and the upper bounds
! others: same as 1 (for now)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::ncalls0,ndim,nitmax,imode
    REAL(KIND(1d0)),DIMENSION(0:nintervals,ndim),INTENT(INOUT)::xgrid
    REAL(KIND(1d0)),DIMENSION(nintervals,ndim),INTENT(INOUT)::ymax ! the upper bounding envelope(all same index cell in every group in kdim-th) of the folded function 
    REAL(KIND(1d0)),DIMENSION(ndimmax)::x,rand,dx
    REAL(KIND(1d0)),INTENT(INOUT)::xint
    REAL(KIND(1d0)),EXTERNAL::fun ! the user function fun=f(h(y))*dhk(yk)/dyk,where dhk(yk)/dyk=w
    REAL(KIND(1d0)),INTENT(OUT)::ans,err
    REAL(KIND(1d0))::vol,f,vtot,etot,prod
    REAL(KIND(1d0)),DIMENSION(0:nintervals,ndimmax)::xacc ! accumulated value of f in each cell
    ! ncell(kdim) is the integer index(1-nintervals/ifold(kdim)) of cell/interval that chosen in each group (1-ifold(kdim))
    ! kfold(kdim) is the index of the chosen group (1-ifold(kdim))
    ! icell(kdim) is the integer index (1-nintervals) of cell/interval that chosen from ncell(kdim) and kfold(kdim)
    INTEGER,DIMENSION(ndimmax)::icell,ncell,kfold
    INTEGER,DIMENSION(1:nintervals,ndimmax)::nhits ! it records the number of PS point happens in each cell
    INTEGER::kdim,kint,kpoint,nit,ncalls,ibin,iret,nintcurr,ifirst
    INTEGER::flag10,flag1,flag101
    IF(imode.EQ.0)THEN
       ! performs the integration of the absolute value of the function, 
       ! finds the optimal grid xgrid, stores the answer in xint and ans
       ! the error in err
       DO kdim=1,ndim
          ifold(kdim)=1 ! first set each dimension is 1 folded.
          ! divide the grid in each dimension as equal intervals
          DO kint=0,nintervals
             xgrid(kint,kdim)=DBLE(kint)/nintervals
          ENDDO
       ENDDO
    ELSEIF(imode.EQ.1)THEN
       ! perfoms the folded integration (folded integration is defined in Section 3 in arXiv:0709.2085) with grid fixed
       ! fbar(z1,z2,...,zn)=(p1*p2*...*pn)^(-1)*Sum^{p1-1}_{l1=0}*Sum^{p2-1}_{l2=0}*...*Sum^{pn-1}_{ln=0}f(h(y))Prod^n_{k=1}dhk/dyk|_{y=y(z,l)}
       ! where yk(zk,lk)=(lk+zk)/pk, lk/pk<yk<(lk+1)/pk with zk spans the [0,1] interval.
       ! int f(x)d^nx=int f(h(y))Prod^n_{k=1}dhk(yk)/dyk*dyk=int fbar(z1,z2,...,zn)dz1dz2...dzn
       ! here pk is just ifold(k)
       ! use the frozen grid,i.e. iflod has been assigned.
       DO kdim=1,ndim
          nintcurr=nintervals/ifold(kdim) ! how many intervals in each group of each dimension kdim
          IF(nintcurr*ifold(kdim).NE.nintervals)THEN
             WRITE(*,*)'mint: the values in the ifold array shoud be divisors of',nintervals
             STOP
          ENDIF
          ! assign the initial value of ymax
          ! xins is just the ans at least in imode=0,i.e. calculated before.
          DO kint=1,nintcurr
             ymax(kint,kdim)=xint**(1d0/ndim)
          ENDDO
       ENDDO
    ENDIF
    ncalls=ncalls0 ! maximum number of calls per iteration
    nit=0 ! records the number of iteration
    ans=0
    err=0
    flag10=0
    DO WHILE(flag10.EQ.0)
! 10   continue
       nit=nit+1
       IF(nit.GT.nitmax)THEN
          IF(imode.EQ.0)xint=ans
          RETURN
       ENDIF
       IF(imode.EQ.0)THEN
          DO kdim=1,ndim
             DO kint=0,nintervals
                xacc(kint,kdim)=0
                IF(kint.GT.0)THEN
                   nhits(kint,kdim)=0
                ENDIF
             ENDDO
          ENDDO
       ENDIF
       vtot=0
       etot=0
       DO kpoint=1,ncalls
! find random x, and its random cell
          DO kdim=1,ndim
             kfold(kdim)=1
             ! randomly choose the cell
             ! a cell means the intervals in each group (1-ifold(kdim)) in kdim
             ncell(kdim)=nintervals/ifold(kdim)*random()+1
             ! random x in this cell
             rand(kdim)=random()
          ENDDO
          f=0
          ifirst=0
! 1       continue
          flag1=0
          DO WHILE(flag1.EQ.0)
             vol=1 ! the volum of the choosen (cell(s)*cell number in each group) in ndim dimnesions
             ! i.e. the volum of the choosen group (if we assume the width of the group is the width of the
             ! choosen cell* cell number in this group) in ndim dimensions
             DO kdim=1,ndim
                ! the number of cells/intervals in each group
                ! kfold records the chosen group
                nintcurr=nintervals/ifold(kdim)
                icell(kdim)=ncell(kdim)+(kfold(kdim)-1)*nintcurr
                ibin=icell(kdim) ! the chosen cell index
                dx(kdim)=xgrid(icell(kdim),kdim)-xgrid(icell(kdim)-1,kdim) ! the width of the chosen cell
                vol=vol*dx(kdim)*nintcurr
                ! the x in this chosen cell.
                x(kdim)=xgrid(icell(kdim)-1,kdim)+rand(kdim)*dx(kdim)
                IF(imode.EQ.0) nhits(ibin,kdim)=nhits(ibin,kdim)+1
             ENDDO
! contribution to integral
             flag1=1
             IF(imode.EQ.0)THEN
                f=ABS(fun(x,vol,ifirst))+f
             ELSE
! this accumulated value will not be used
! accumulated value from kfold group to ifold group (last group)
                f=fun(x,vol,ifirst)+f
                ifirst=1
                ! next group
                CALL nextlexi(ndim,ifold,kfold,iret)
                IF(iret.EQ.0)flag1=0 ! goto 1
! closing call: accumulated value with correct sign
                IF(flag1.EQ.1)f=fun(x,vol,2)
             ENDIF
          ENDDO ! END OF DO WHILE(flag1.EQ.0)
!
          IF(imode.EQ.0)THEN
! accumulate the function in xacc(icell(kdim),kdim) to adjust the grid later
             DO kdim=1,ndim
                xacc(icell(kdim),kdim)=xacc(icell(kdim),kdim)+f
             ENDDO
          ELSE
! update the upper bounding envelope
             prod=1
             DO kdim=1,ndim
                prod=prod*ymax(ncell(kdim),kdim)
             ENDDO
             prod=(f/prod)
             IF(prod.GT.1)THEN
! This guarantees a 10% increase of the upper bound in this cell
                prod=1+0.1d0/ndim
                DO kdim=1,ndim
                   ymax(ncell(kdim),kdim)=ymax(ncell(kdim),kdim)*prod
                ENDDO
             ENDIF
          ENDIF
          vtot=vtot+f/ncalls ! Eq.(2.4) in arXiv:0709.2085
          etot=etot+f**2/ncalls
       ENDDO ! END OF DO kpoint=1,ncalls
       IF(imode.EQ.0)THEN
! iteration is finished; now rearrange the grid
          DO kdim=1,ndim
             CALL regrid(xacc(0,kdim),xgrid(0,kdim),nhits(1,kdim),nintervals,nit)
          ENDDO
       ENDIF
! the abs is to avoid tiny negative values
       etot=SQRT(ABS(etot-vtot**2)/ncalls) ! standard deviation
       WRITE(*,*) vtot,etot
       IF(nit.EQ.1)THEN
          ans=vtot
          err=etot
       ELSE
! prevent annoying division by zero for nearly zero
! integrands
          flag101=0
          IF(etot.EQ.0.AND.err.EQ.0)THEN
             IF(ans.EQ.vtot)THEN
                flag101=1
!                goto 10
             ELSE
                err=ABS(vtot-ans)
                etot=ABS(vtot-ans)
             ENDIF
          ELSEIF(etot.EQ.0)THEN
             etot=err
          ELSEIF(err.EQ.0)THEN
             err=etot
          ENDIF
          IF(flag101.EQ.0)THEN
             ! the formula for the combining different interations
             ! which is the same as in VEGAS
             ans=(ans/err+vtot/etot)/(1/err+1/etot)
             err=1/SQRT(1/err**2+1/etot**2)
          ENDIF
       ENDIF
!    goto 10
    ENDDO ! END OF DO WHILE(flag10.EQ.0)
  END SUBROUTINE mint

  FUNCTION random()
! a warper for the function Helac_rnmy in HELAC-Onia
    USE Helac_ranmar_mod
    REAL(KIND(1d0))::random
    random=Helac_rnmy(0)
  END FUNCTION random

  SUBROUTINE regrid(xacc,xgrid,nhits,nint,nit)
! rearrange the grid for each dimension
    IMPLICIT NONE
    ! nint=nintervals
    INTEGER,INTENT(IN)::nint,nit ! nit is the index for the current interation
    INTEGER,DIMENSION(nint),INTENT(IN)::nhits
    REAL(KIND(1d0)),DIMENSION(0:nint),INTENT(INOUT)::xacc,xgrid
    REAL(KIND(1d0)),DIMENSION(100)::xn
    REAL(KIND(1d0))::r
    INTEGER::kint,jint,flag11
    ! xacc(0)=0
    DO kint=1,nint
! xacc (xerr) already containe a factor equal to the interval size
! Thus the integral of rho is performed by summing up
       IF(nhits(kint).NE.0)THEN
          xacc(kint)= xacc(kint-1)+ABS(xacc(kint))/nhits(kint)
       ELSE
          xacc(kint)=xacc(kint-1)
       ENDIF
    ENDDO
    ! determine the fraction
    DO kint=1,nint
       xacc(kint)=xacc(kint)/xacc(nint)
    ENDDO
    WRITE(11,*) 'set limits x 0 1 y 0 1'
    WRITE(11,*) 0, 0
    DO kint=1,nint
       WRITE(11,*) xgrid(kint),xacc(kint)
    ENDDO
    WRITE(11,*) 'join 0'

    DO kint=1,nint
       r=DBLE(kint)/nint

       WRITE(11,*) 0, r
       WRITE(11,*) 1, r
       WRITE(11,*) ' join'
       flag11=0
       DO jint=1,nint
          IF(r.LT.xacc(jint))THEN
             ! main regriding function
             ! the larger value the cell width is smaller
             xn(kint)=xgrid(jint-1)+(r-xacc(jint-1))&
                  /(xacc(jint)-xacc(jint-1))*(xgrid(jint)-xgrid(jint-1))
             flag11=1
             EXIT
!             goto 11
          ENDIF
       ENDDO
       IF(flag11.EQ.0)THEN
          IF(jint.NE.nint+1.AND.kint.NE.nint)THEN
             WRITE(*,*) ' error',jint,nint
             STOP
          ENDIF
          xn(nint)=1
       ENDIF
! 11      continue
    ENDDO
    DO kint=1,nint
       xgrid(kint)=xn(kint)
!         xgrid(kint)=(xn(kint)+2*xgrid(kint))/3
!         xgrid(kint)=(xn(kint)+xgrid(kint)*log(dble(nit)))
!     #        /(log(dble(nit))+1)
       WRITE(11,*) xgrid(kint), 0
       WRITE(11,*) xgrid(kint), 1
       WRITE(11,*) ' join'
    ENDDO
    WRITE(11,*) ' newplot'
  END SUBROUTINE regrid

  SUBROUTINE nextlexi(ndim,iii,kkk,iret)
! kkk: array of integers 1 <= kkk(j) <= iii(j), j=1,ndim
! at each call iii is increased lexicographycally.
! for example, starting from ndim=3, kkk=(1,1,1), iii=(2,3,2)
! subsequent calls to nextlexi return
!         kkk(1)      kkk(2)      kkk(3)    iret
! 0 calls   1           1           1       0
! 1         1           1           2       0    
! 2         1           2           1       0
! 3         1           2           2       0
! 4         1           3           1       0
! 5         1           3           2       0
! 6         2           1           1       0
! 7         2           1           2       0
! 8         2           2           1       0
! 9         2           2           2       0
! 10        2           3           1       0
! 11        2           3           2       0
! 12        2           3           2       1
    IMPLICIT NONE
    INTEGER,INTENT(IN)::ndim
    INTEGER,INTENT(OUT)::iret
    INTEGER,DIMENSION(ndim),INTENT(INOUT)::kkk
    INTEGER,DIMENSION(ndim),INTENT(IN)::iii
    INTEGER::k,flag1
    k=ndim
    flag1=0
! 1    continue
    DO WHILE(flag1.EQ.0)
       IF(kkk(k).LT.iii(k))THEN
          kkk(k)=kkk(k)+1
          iret=0
          RETURN
       ELSE
          kkk(k)=1
          k=k-1
          IF(k.EQ.0)THEN
             iret=1
             RETURN
          ENDIF
!         goto 1
       ENDIF
    ENDDO
  END SUBROUTINE nextlexi


  SUBROUTINE gen(fun,ndim,xgrid,ymax,imode,x)
! generate events with the distribution of the folded function
! Usage:
!    imode=0
!    call gen(fun,ndim,xgrid,ymax,imode,x)
!    imode=1
!    do j=1,10000
!    call gen(fun,ndim,xgrid,ymax,imode,x)
!    ...
!    enddo
!    imode=3
!    call gen(fun,ndim,xgrid,ymax,imode,x)
! End of Usage
! imode=0 to initialize the generation
! imode=1 to generate the ndim-tuples,see Section 4 in arXiv:0709.2085
! imode=3 store generation efficiency in x(1), prints out generation statistics
    IMPLICIT NONE
    INTEGER,INTENT(IN)::ndim,imode
    REAL(KIND(1d0)),EXTERNAL::fun
    REAL(KIND(1d0)),DIMENSION(0:nintervals,ndim),INTENT(OUT)::xgrid
    REAL(KIND(1d0)),DIMENSION(nintervals,ndim),INTENT(OUT)::ymax
    REAL(KIND(1d0)),DIMENSION(ndim),INTENT(OUT)::x
    REAL(KIND(1d0)),DIMENSION(ndimmax)::dx
    INTEGER,DIMENSION(ndimmax)::icell,ncell,kfold
    REAL(KIND(1d0))::r,f,ubound,vol
    REAL(KIND(1d0)),DIMENSION(ndimmax)::rand
    REAL(KIND(1d0)),DIMENSION(nintervals,ndimmax)::xmmm
    INTEGER::icalls,mcalls,kdim,kint,nintcurr,iret,ifirst
    SAVE icalls,mcalls,xmmm
    INTEGER::flag10,flag5
    IF(imode.EQ.0)THEN
       DO kdim=1,ndim
          nintcurr=nintervals/ifold(kdim)
          xmmm(1,kdim)=ymax(1,kdim)
          DO kint=2,nintcurr
             xmmm(kint,kdim)=xmmm(kint-1,kdim)+ymax(kint,kdim)
          ENDDO
          DO kint=1,nintcurr
             xmmm(kint,kdim)=xmmm(kint,kdim)/xmmm(nintcurr,kdim)
          ENDDO
       ENDDO
       icalls=0
       mcalls=0
       RETURN
    ELSEIF(imode.EQ.3)THEN
       IF(icalls.GT.0)THEN
          x(1)=DBLE(mcalls)/icalls
       ELSE
          x(1)=-1
       ENDIF
       RETURN
    ENDIF
    mcalls=mcalls+1
    flag10=0
! generate one unweighted event
    DO WHILE(flag10.EQ.0)
! 10   continue
       DO kdim=1,ndim
          nintcurr=nintervals/ifold(kdim)
          r=random()
          ! select the cell index in each group ncell following the probability of xmmm,i.e.ymax
          DO kint=1,nintcurr
             IF(r.LT.xmmm(kint,kdim))THEN
                ncell(kdim)=kint
                EXIT
!                goto 1
             ENDIF
          ENDDO
! 1       continue
          rand(kdim)=random()
       ENDDO
       ! the bound value for this cell i.e. ncell in every group
       ubound=1
       DO kdim=1,ndim
          ubound=ubound*ymax(ncell(kdim),kdim)
       ENDDO
       DO kdim=1,ndim
          kfold(kdim)=1 ! take the first group in each dimension
       ENDDO
       f=0
       ifirst=0
       flag5=0
       ! exhaust all groups
       DO WHILE(flag5.EQ.0)
! 5    continue
          vol=1
          DO kdim=1,ndim
             nintcurr=nintervals/ifold(kdim)
             icell(kdim)=ncell(kdim)+(kfold(kdim)-1)*nintcurr
             dx(kdim)=xgrid(icell(kdim),kdim)-xgrid(icell(kdim)-1,kdim)
             vol=vol*dx(kdim)*nintervals/ifold(kdim)
             x(kdim)=xgrid(icell(kdim)-1,kdim)+rand(kdim)*dx(kdim)
          ENDDO
          f=f+fun(x,vol,ifirst)
          ifirst=1
          ! next group
          CALL nextlexi(ndim,ifold,kfold,iret)
          flag5=1
          IF(iret.EQ.0)flag5=0   !  goto 5
       ENDDO ! END OF DO WHILE(flag5.EQ.0)
! get final value (x and vol not used in this call)
       f=fun(x,vol,2)
       IF(f.LT.0)THEN
          WRITE(*,*) 'gen: non positive function'
          STOP
       ENDIF
       IF(f.GT.ubound)THEN
          CALL increasecnt('upper bound failure in inclusive cross section')
       ENDIF
       ubound=ubound*random() ! this is the unweighted procedure
       icalls=icalls+1
       flag10=1
       IF(ubound.GT.f)THEN
          ! at this time, the event should *not* be recorded
          CALL increasecnt('vetoed calls in inclusive cross section')
          flag10=0
!       goto 10
       ENDIF
       ! otherwise, the event should be recorded
    ENDDO ! END OF DO WHILE(flag10.EQ.0)
  END SUBROUTINE gen

  SUBROUTINE increasecnt(str)
    IMPLICIT NONE
    CHARACTER(len=*),INTENT(IN)::str
    WRITE(*,*)str
  END SUBROUTINE increasecnt

END MODULE mint_integrator
