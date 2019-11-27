c
c Example analysis for "p p > t t~ h, (t > w+ b, w+ > lj lj),(t > w- b~, w- > lj lj), h > a a" process.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine analysis_begin(nwgt,weights_info)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer idim,n_plot
      parameter (idim=7,n_plot=20)
C      common/plot_analy/idim,n_plot
      integer nwgt
      character*(*) weights_info(*)
      integer i,kk,l,nwgt_analysis
      common/c_analysis/nwgt_analysis
      character*5 cc(idim)
      ! hadhad,lephad,leplep,hadtau,leptau,tautau 
      data cc/'pp LO','pp hh','pp lh','pp ll','pp ht','pp lt','pp tt'/
      double precision pi
      PARAMETER (PI=3.14159265358979312D0)
      include 'dbook.inc'
      call inihist
      nwgt_analysis=nwgt
      if (nwgt_analysis*idim*n_plot.gt.nplots/4) then
         write (*,*) 'error in analysis_begin: '/
     &        /'too many histograms, increase NPLOTS to',
     &        nwgt_analysis*idim*n_plot*4
         stop 1
      endif
      do kk=1,nwgt_analysis
      do i=1,idim
        l=(kk-1)*idim*n_plot+(i-1)*n_plot
        call bookup(l+ 1,'tt pt            '
     &       //weights_info(kk)//cc(i),10.d0,0.d0,5000.d0)
        call bookup(l+ 2,'del R min  '
     &       //weights_info(kk)//cc(i),pi/20.d0,0.d0,3*pi)
        call bookup(l+ 3,'ba min del R         '
     &       //weights_info(kk)//cc(i),pi/20.d0,0.d0,3*pi)
        call bookup(l+ 4,'aa del R    '
     &       //weights_info(kk)//cc(i),pi/20.d0,0.d0,3*pi)
        call bookup(l+ 5,'wa min del R         '
     &       //weights_info(kk)//cc(i),pi/20.d0,0.d0,3*pi)
        call bookup(l+ 6,'tb pt            '
     &       //weights_info(kk)//cc(i),10.d0,0.d0,5000.d0)
        call bookup(l+ 7,'minv aa       '
     &       //weights_info(kk)//cc(i),0.01d0,86.d0,96.d0)
        call bookup(l+ 8,'t pt             '
     &       //weights_info(kk)//cc(i),10.d0,0.d0,5000.d0)
        call bookup(l+ 9,'h pt        '
     &       //weights_info(kk)//cc(i),10.d0,0.d0,5000.d0)
        call bookup(l+10,'tt delta eta     '
     &       //weights_info(kk)//cc(i),0.2d0,-4.d0,4.d0)
        call bookup(l+11,'y_tt             '
     &       //weights_info(kk)//cc(i),0.1d0,-4.d0,4.d0)
        call bookup(l+12,'delta y          '
     &       //weights_info(kk)//cc(i),0.2d0,-4.d0,4.d0)
        call bookup(l+13,'tt azimt         '
     &       //weights_info(kk)//cc(i),pi/60.d0,2*pi/3,pi)
        call bookup(l+14,'tt del R         '
     &       //weights_info(kk)//cc(i),pi/60.d0,2*pi/3,4*pi/3)
        call bookup(l+15,'delta(y_tt,h)    '
     &       //weights_info(kk)//cc(i),0.1d0,-4.d0,4.d0)
        call bookup(l+16,'y_t              '
     &       //weights_info(kk)//cc(i),0.1d0,-4.d0,4.d0)
        call bookup(l+17,'tt log[pi-azimt] '
     &       //weights_info(kk)//cc(i),0.05d0,-4.d0,0.1d0)
        call bookup(l+18,'h pt            '
     &       //weights_info(kk)//cc(i),10.d0,0.d0,5000.d0)
        call bookup(l+19,'minv tot            '
     &       //weights_info(kk)//cc(i),10.d0,0.d0,5000.d0)
        call bookup(l+20,'minv tx h             '
     &       //weights_info(kk)//cc(i),10.d0,0.d0,5000.d0)
      enddo
      enddo
      return
      end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine analysis_end(xnorm)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer idim,n_plot
      parameter (idim=7,n_plot=20)
C      common/plot_analy/idim,n_plot
      character*14 ytit
      double precision xnorm
      integer i
      integer kk,l,nwgt_analysis
      common/c_analysis/nwgt_analysis
      include 'dbook.inc'
      call open_topdrawer_file
      call mclear
      do i=1,NPLOTS
         call mopera(i,'+',i,i,xnorm,0.d0)
         call mfinal(i)
      enddo
      ytit='sigma per bin '
      do kk=1,nwgt_analysis
      do i=1,idim
        l=(kk-1)*idim*n_plot+(i-1)*n_plot
        call multitop(l+ 1,2,3,'tt pt',' ','LOG')
        call multitop(l+ 2,2,3,'del R min',' ','LOG')
        call multitop(l+ 3,2,3,'ba min del R',' ','LOG')
        call multitop(l+ 4,2,3,'aa del R',' ','LOG')
        call multitop(l+ 5,2,3,'wa min del R',' ','LOG')
        call multitop(l+ 6,2,3,'tb pt',' ','LOG')
        call multitop(l+ 7,2,3,'minv aa',' ','LOG')
        call multitop(l+ 8,2,3,'t pt',' ','LOG')
        call multitop(l+ 9,2,3,'pt h',' ','LOG')
        call multitop(l+10,2,3,'tt Delta eta',' ','LOG')
        call multitop(l+11,2,3,'y_tt',' ','LOG')
        call multitop(l+12,2,3,'tt Delta y',' ','LOG')
        call multitop(l+13,2,3,'tt azimt',' ','LOG')
        call multitop(l+14,2,3,'tt del R',' ','LOG')
        call multitop(l+15,2,3,'delta(y_tt,h)',' ','LOG')
        call multitop(l+16,2,3,'t y',' ','LOG')
        call multitop(l+17,2,3,'tt log[pi-azimt]',' ','LOG')
        call multitop(l+18,2,3,'tt pth',' ','LOG')
        call multitop(l+19,2,3,'minv tot',' ','LOG')
        call multitop(l+20,2,3,'minv tx h',' ','LOG')
      enddo
      enddo
      call close_topdrawer_file
      return                
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine analysis_fill(p,nexternal,ipdg,wgts)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer idim,n_plot
      parameter (idim=7,n_plot=20)
C      common/plot_analy/idim,n_plot
      integer nexternal
      integer iPDG(nexternal)
      double precision p(0:4,nexternal)
      double precision wgts(*)
      integer ibody,order
      double precision wgt,var
      integer i,ii,kk,l,nwgt_analysis
      common/c_analysis/nwgt_analysis
      double precision www,xptq(0:3),xptb(0:3),yptqtb(0:3),ycut,ptcut
     $     ,ptq1,ptq2,ptg,yq1,yq2,etaq1,etaq2,azi,azinorm,qqm,dr,yqq,yh,
     $     xpta1(0:3),xpta2(0:3),yptaa(0:3),ya1,ya2,drtamin,azita,drta,
     $     pta1,pta2,ptamin,etaamax,drmin,draa,xpwp(0:3),xpwm(0:3),
     $     xpbq(0:3),xpbb(0:3),drba,drbamin,drwa,drwamin,ywp,ywm,ybq,ybb
      double precision xpth(0:3),yptqtbth(0:3),yptqth(0:3),ptamax,
     $ yptbth(0:3),Ht,minvtot,minvth,minvtxh,pth,minvaa,minvwpb,minvwmbx
     $ ,yptwpb(0:3),yptwmbx(0:3)    
      logical siq1flag,siq2flag,siq3flag,cutflag,ggflag,uuflag,ddflag,
     . ugflag,dgflag,cutflag2,leplep(6)
      integer nlep,ntau
      double precision getrapidity,getpseudorap,getinvm,getdelphi
      external getrapidity,getpseudorap,getinvm,getdelphi
      double precision pi
      PARAMETER (PI=3.14159265358979312D0)
c      print*, "i pdg quando nn si blocca sono",ipdg(1),ipdg(2)
c      if (nexternal.ne.10) then
c         print *, nexternal
c         write (*,*) 'error #1 in analysis_fill: '/
c     &        /'only for process "p p > t t~ > w+ w- b b~ a a"'
c         stop 1
c      endif
      if (.not. (abs(ipdg(1)).le.5 .or. ipdg(1).eq.21)) then
         write (*,*) 'error #2 in analysis_fill: '/
     &        /'only for process "p p > t t~ > w+ w- b b~ a a"'
         stop 1
      endif
      if (.not. (abs(ipdg(2)).le.5 .or. ipdg(2).eq.21)) then
         write (*,*) 'error #3 in analysis_fill: '/
     &        /'only for process "p p > t t~ > w+ w- b b~ a a"'
         stop 1
      endif
c      if (abs(ipdg(3)).ne.6) then
c         write (*,*) 'error #4 in analysis_fill: '/
c     &        /'only for process "p p > t t~ > w+ w- b b~ a a"'
c         stop 1
c      endif
c      if (abs(ipdg(4)).ne.6) then
c         write (*,*) 'error #5 in analysis_fill: '/
c     &        /'only for process "p p > t t~ > w+ w- b b~ a a"'
c         stop 1
c      endif
c      if (ipdg(5).ne.24) then
c         write (*,*) 'error #6 in analysis_fill: '/
c     &        /'only for process "p p > t t~ > w+ w- b b~ a a"'
c         stop 1
c      endif
c      if (ipdg(6).ne.-24) then
c         write (*,*) 'error #7 in analysis_fill: '/
c     &        /'only for process "p p > t t~ > w+ w- b b~ a a"'
c         stop 1
c      endif
c      if (ipdg(7).ne.5) then
c         write (*,*) 'error #8 in analysis_fill: '/
c     &        /'only for process "p p > t t~ > w+ w- b b~ a a"'
c         stop 1
c      endif
c      if (ipdg(8).ne.-5) then
c         write (*,*) 'error #9 in analysis_fill: '/
c     &        /'only for process "p p > t t~ > w+ w- b b~ a a"'
c         stop 1
c      endif
      if (abs(ipdg(nexternal-1)).ne.11.and.abs(ipdg(nexternal-1)).ne.13) then
          write (*,*) 'error #10 in analysis_fill: '/
     &        /'only for process "p p > t t~ > w+ w- b b~ a a"'
         stop 1
      endif
      if (abs(ipdg(nexternal)).ne.11.and.abs(ipdg(nexternal)).ne.13) then
          write (*,*) 'error #11 in analysis_fill: '/
     &        /'only for process "p p > t t~ > w+ w- b b~ a a"'
         stop 1
      endif

      nlep=0 ! mu+,e+,mu-,e-
      ntau=0 ! ta+,ta-
      do ii=1,6
         leplep(ii)=.FALSE.
      enddo
      do ii=3,nexternal-2
         if(abs(ipdg(ii)).EQ.11.OR.abs(ipdg(ii)).EQ.13)nlep=nlep+1
         if(abs(ipdg(ii)).EQ.15)ntau=ntau+1
      enddo
      ! hadhad,lephad,leplep,hadtau,leptau,tautau
      SELECT CASE(nlep)
      CASE(0)
         SELECT CASE(ntau)
         CASE(0)
            leplep(1)=.TRUE.
         CASE(1)
            leplep(4)=.TRUE.
         CASE(2)
            leplep(6)=.TRUE.
         END SELECT
      CASE(1)
         SELECT CASE(ntau)
         CASE(0)
            leplep(2)=.TRUE.
         CASE(1)
            leplep(5)=.TRUE.
         END SELECT
      CASE(2)
         leplep(3)=.TRUE.
      END SELECT

      do i=0,3
c        xpth(i)=p(i,3)
        xptq(i)=0d0
        xptb(i)=0d0
        xpta1(i)=p(i,nexternal-1)
        xpta2(i)=p(i,nexternal)
        xpwp(i)=0d0
        xpwm(i)=0d0
        xpbq(i)=0d0
        xpbb(i)=0d0
        do ii=3,nexternal
           if(ipdg(ii).EQ.6)xptq(i)=p(i,ii)
           if(ipdg(ii).EQ.-6)xptb(i)=p(i,ii)
           if(ipdg(ii).EQ.24)xpwp(i)=p(i,ii)
           if(ipdg(ii).EQ.-24)xpwm(i)=p(i,ii)
           if(ipdg(ii).EQ.5)xpbq(i)=p(i,ii)
           if(ipdg(ii).EQ.-5)xpbb(i)=p(i,ii)
        enddo
        xpth(i)=xpta1(i)+xpta2(i)
        yptqtb(i)=xptq(i)+xptb(i)
        yptqtbth(i)=xptq(i)+xptb(i)+xpth(i)
        yptbth(i)=xptb(i)+xpth(i)
        yptqth(i)=xptq(i)+xpth(i)
        yptaa(i)=xpta1(i)+xpta2(i)
        yptwpb(i)=xpwp(i)+xpbq(i)
        yptwmbx(i)=xpwm(i)+xpbb(i)
      enddo
      PTCUT=200.D0
C
      ptq1 = dsqrt(xptq(1)**2+xptq(2)**2)
      ptq2 = dsqrt(xptb(1)**2+xptb(2)**2)
      pth = dsqrt(xpth(1)**2+xpth(2)**2)
      ptg = dsqrt(yptqtb(1)**2+yptqtb(2)**2)
      pta1 = dsqrt(xpta1(1)**2+xpta1(2)**2)
      pta2 = dsqrt(xpta2(1)**2+xpta2(2)**2)
      ptamin = min(pta1,pta2)
      ptamax = min(pta1,pta2)
      yq1=getrapidity(xptq(0),xptq(3))
      yq2=getrapidity(xptb(0),xptb(3))
      ywp=getrapidity(xpwp(0),xpwp(3))
      ywm=getrapidity(xpwm(0),xpwm(3))
      ybq=getrapidity(xpbq(0),xpbq(3))
      ybb=getrapidity(xpbb(0),xpbb(3))
      yh=getrapidity(xpth(0),xpth(3))
      ya1=getrapidity(xpta1(0),xpta1(3))
      ya2=getrapidity(xpta2(0),xpta2(3))
      etaamax=max(abs(ya1),abs(ya2))
C     a1 a2
      azita=getdelphi(xpta1(1),xpta1(2),xpta2(1),xpta2(2))
      draa=dsqrt(azita**2+(ya1-ya2)**2)
C     w+ a1
      azita=getdelphi(xpwp(1),xpwp(2),xpta1(1),xpta1(2))
      drwa = dsqrt(azita**2+(ya1-ywp)**2)
      drwamin = drwa
C     w- a1
      azita=getdelphi(xpwm(1),xpwm(2),xpta1(1),xpta1(2))
      drwa = dsqrt(azita**2+(ya1-ywm)**2)
      if(drwa.lt.drwamin)drwamin=drwa
C     w+ a2
      azita=getdelphi(xpwp(1),xpwp(2),xpta2(1),xpta2(2))
      drwa = dsqrt(azita**2+(ya2-ywp)**2)
      if(drwa.lt.drwamin)drwamin=drwa
C     w- a2
      azita=getdelphi(xpwm(1),xpwm(2),xpta2(1),xpta2(2))
      drwa = dsqrt(azita**2+(ya2-ywm)**2)
      if(drwa.lt.drwamin)drwamin=drwa
C     b a1
      azita=getdelphi(xpbq(1),xpbq(2),xpta1(1),xpta1(2))
      drba = dsqrt(azita**2+(ya1-ybq)**2)
      drbamin = drba
C     b~ a1
      azita=getdelphi(xpbb(1),xpbb(2),xpta1(1),xpta1(2))
      drba = dsqrt(azita**2+(ya1-ybb)**2)
      if(drba.lt.drbamin)drbamin=drba
C     b a2
      azita=getdelphi(xpbq(1),xpbq(2),xpta2(1),xpta2(2))
      drba = dsqrt(azita**2+(ya2-ybq)**2)
      if(drba.lt.drbamin)drbamin=drba
C     b~ a2
      azita=getdelphi(xpbb(1),xpbb(2),xpta2(1),xpta2(2))
      drba = dsqrt(azita**2+(ya2-ybb)**2)
      if(drba.lt.drbamin)drbamin=drba

C     t t~
      azita=getdelphi(xptq(1),xptq(2),xptb(1),xptb(2))
      drta = dsqrt(azita**2+(yq1-yq2)**2)
      drmin=drta
C     t h
      azita=getdelphi(xptq(1),xptq(2),xpth(1),xpth(2))
      drta = dsqrt(azita**2+(yq1-yh)**2)
      if(drta.lt.drmin)drmin=drta
C     t~ h
      azita=getdelphi(xptb(1),xptb(2),xpth(1),xpth(2))
      drta = dsqrt(azita**2+(yq2-yh)**2)
      if(drta.lt.drmin)drmin=drta

      etaq1=getpseudorap(xptq(0),xptq(1),xptq(2),xptq(3))
      etaq2=getpseudorap(xptb(0),xptb(1),xptb(2),xptb(3))
      azi=getdelphi(xptq(1),xptq(2),xptb(1),xptb(2))
      azinorm = (pi-azi)/pi
      qqm=getinvm(yptqtb(0),yptqtb(1),yptqtb(2),yptqtb(3))
      dr  = dsqrt(azi**2+(etaq1-etaq2)**2)
      yqq=getrapidity(yptqtb(0),yptqtb(3))
      Ht=ptq1+ptq2+pth
      minvtot=getinvm(yptqtbth(0),yptqtbth(1),yptqtbth(2),yptqtbth(3))
      minvth=getinvm(yptqth(0),yptqth(1),yptqth(2),yptqth(3))
      minvtxh=getinvm(yptbth(0),yptbth(1),yptbth(2),yptbth(3))
      minvaa=getinvm(yptaa(0),yptaa(1),yptaa(2),yptaa(3))
      minvwpb=getinvm(yptwpb(0),yptwpb(1),yptwpb(2),yptwpb(3))
      minvwmbx=getinvm(yptwmbx(0),yptwmbx(1),yptwmbx(2),yptwmbx(3))
c-------------------------------------------------------------
C      siq1flag=ptq1.gt.ptcut!.and.abs(yq1).lt.ycut
C      siq2flag=ptq2.gt.ptcut!.and.abs(yq2).lt.ycut
C      siq3flag=pth.gt.ptcut
C      cutflag=siq1flag.and.siq2flag.and.siq3flag
      cutflag=ptamin.GT.20d0.AND.etaamax.LT.2.5d0
     &     .AND.minvaa.GE.86d0.AND.minvaa.LE.96d0
     &     .AND.(ptamax.GT.minvaa/2d0)
      cutflag2=ptamin.GT.20d0.AND.etaamax.LT.2.5d0
     &     .AND.minvaa.GE.86d0.AND.minvaa.LE.96d0
c      ggflag=(ipdg(1).eq.21) .and. (ipdg(2).eq.21)
c      uuflag=(abs(ipdg(1)).eq.2) .and.(abs(ipdg(2)).eq.2)
c      ddflag=(abs(ipdg(1)).eq.1) .and.(abs(ipdg(2)).eq.1)
c      ugflag=((abs(ipdg(1)).eq.2).and.(ipdg(2).eq.21)).or.
c     .((ipdg(1).eq.21) .and.(abs(ipdg(2)).eq.2))
c      dgflag=((abs(ipdg(1)).eq.1).and.(ipdg(2).eq.21)).or.
c     .((ipdg(1).eq.21) .and.(abs(ipdg(2)).eq.1))
c      if((.not.ggflag).and.(.not.uuflag).and.(.not.ddflag).and.(.not.ugflag).and.(.not.dgflag)) then
c      print*, "not uu dd gg ug dg",ipdg(1),ipdg(2)
c      stop 
c      endif
c-------------------------------------------------------------
c      double precision xpth(0:3),yptqtbth(0:3),yptqth(0:3),
c     $ yptbth(0:3),Ht,minvtot,minvth,minvtxh
c      double precision xpth(0:3),yptqtbth(0:3),yptqth(0:3),
c     $ yptbth(0:3),Ht,minvtot,minvth,minvtxh,pth

      do kk=1,nwgt_analysis
         www=wgts(kk)
         l=(kk-1)*idim*n_plot

         if(cutflag2)then
c         call mfill(l+1,minvaa,WWW)
         call mfill(l+1,ptg,WWW)
         call mfill(l+18,pth,WWW)
         call mfill(l+2,drmin,WWW)
         call mfill(l+3,drbamin,WWW)
         call mfill(l+4,draa,WWW)
         call mfill(l+13,azi,WWW)
         if(azinorm.gt.0)call mfill(l+17,log10(azinorm),WWW)
         call mfill(l+5,drwamin,WWW)
         call mfill(l+14,dr,WWW)
         call mfill(l+10,etaq1-etaq2,WWW)
         call mfill(l+11,yqq,WWW)
         call mfill(l+12,yq1-yq2,WWW)
         call mfill(l+6,ptq2,WWW)
         call mfill(l+19,minvtot,WWW)
         call mfill(l+7,minvaa,WWW)
         call mfill(l+15,yqq-yh,WWW)
         call mfill(l+8,ptq1,WWW)
         call mfill(l+20,minvtxh,WWW)
         call mfill(l+9,ptamin,WWW)
         call mfill(l+16,yq1,WWW)
         else
            print *,ptamin,etaamax
            stop
         endif
c
c***************************************************** with cuts
c
         do ii=1,idim-1
         l=l+n_plot
c
         if(leplep(ii))then
c         call mfill(l+1,minvaa,WWW)
         call mfill(l+1,ptg,WWW)
         call mfill(l+18,pth,WWW)
         call mfill(l+2,drmin,WWW)
         call mfill(l+3,drbamin,WWW)
         call mfill(l+4,draa,WWW)
         call mfill(l+13,azi,WWW)
         if(azinorm.gt.0)call mfill(l+17,log10(azinorm),WWW)
         call mfill(l+5,drwamin,WWW)
         call mfill(l+14,dr,WWW)
         call mfill(l+10,etaq1-etaq2,WWW)
         call mfill(l+11,yqq,WWW)
         call mfill(l+12,yq1-yq2,WWW)
         call mfill(l+6,ptq2,WWW)
         call mfill(l+19,minvtot,WWW)
         call mfill(l+7,minvaa,WWW)
         call mfill(l+15,yqq-yh,WWW)
         call mfill(l+8,ptq1,WWW)
         call mfill(l+20,minvtxh,WWW)
         call mfill(l+9,ptamin,WWW)
         call mfill(l+16,yq1,WWW)
         endif
      enddo
      enddo
 999  return      
      end


      SUBROUTINE open_topdrawer_file
      IMPLICIT NONE
      OPEN(unit=99,FILE='plot_lhe.top',status='unknown')
      END SUBROUTINE open_topdrawer_file

      SUBROUTINE close_topdrawer_file
      IMPLICIT NONE
      CLOSE(99)
      RETURN
      END SUBROUTINE close_topdrawer_file

      function getrapidity(en,pl)
      implicit none
      real*8 getrapidity,en,pl,tiny,xplus,xminus,y
      parameter (tiny=1.d-8)
      xplus=en+pl
      xminus=en-pl
      if(xplus.gt.tiny.and.xminus.gt.tiny)then
         if( (xplus/xminus).gt.tiny.and.(xminus/xplus).gt.tiny)then
            y=0.5d0*log( xplus/xminus  )
         else
            y=sign(1.d0,pl)*1.d8
         endif
      else 
         y=sign(1.d0,pl)*1.d8
      endif
      getrapidity=y
      return
      end


      function getpseudorap(en,ptx,pty,pl)
      implicit none
      real*8 getpseudorap,en,ptx,pty,pl,tiny,pt,eta,th
      parameter (tiny=1.d-5)
c
      pt=sqrt(ptx**2+pty**2)
      if(pt.lt.tiny.and.abs(pl).lt.tiny)then
        eta=sign(1.d0,pl)*1.d8
      else
        th=atan2(pt,pl)
        eta=-log(tan(th/2.d0))
      endif
      getpseudorap=eta
      return
      end


      function getinvm(en,ptx,pty,pl)
      implicit none
      real*8 getinvm,en,ptx,pty,pl,tiny,tmp
      parameter (tiny=1.d-5)
c
      tmp=en**2-ptx**2-pty**2-pl**2
      if(tmp.gt.0.d0)then
        tmp=sqrt(tmp)
      elseif(tmp.gt.-tiny)then
        tmp=0.d0
      else
         tmp=-1d0
c        write(*,*)'Attempt to compute a negative mass'
c        stop
      endif
      getinvm=tmp
      return
      end


      function getdelphi(ptx1,pty1,ptx2,pty2)
      implicit none
      real*8 getdelphi,ptx1,pty1,ptx2,pty2,tiny,pt1,pt2,tmp
      parameter (tiny=1.d-5)
c
      pt1=sqrt(ptx1**2+pty1**2)
      pt2=sqrt(ptx2**2+pty2**2)
      if(pt1.ne.0.d0.and.pt2.ne.0.d0)then
        tmp=ptx1*ptx2+pty1*pty2
        tmp=tmp/(pt1*pt2)
        if(abs(tmp).gt.1.d0+tiny)then
          write(*,*)'Cosine larger than 1'
          stop
        elseif(abs(tmp).ge.1.d0)then
          tmp=sign(1.d0,tmp)
        endif
        tmp=acos(tmp)
      else
        tmp=1.d8
      endif
      getdelphi=tmp
      return
      end
