c
c Example analysis for "p p > t t~ h, (t > w+ b, w+ > lj lj),(t > w- b~, w- > lj lj), h > b b~" process.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine analysis_begin(nwgt,weights_info)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer idim,n_plot
      parameter (idim=7,n_plot=21)
C      common/plot_analy/idim,n_plot
      integer nwgt
      character*(*) weights_info(*)
      integer i,kk,l,nwgt_analysis
      common/c_analysis/nwgt_analysis
      character*5 cc(idim)
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
        call bookup(l+ 1,'t pt            '
     &       //weights_info(kk)//cc(i),10.d0,0.d0,5000.d0)
        call bookup(l+ 2,'h pt  '
     &       //weights_info(kk)//cc(i),10.d0,0.d0,5000.d0)
        call bookup(l+ 3,'j pt  '
     &       //weights_info(kk)//cc(i),1.d0,0.d0,500.d0)
        call bookup(l+ 4,'b pt   '
     &       //weights_info(kk)//cc(i),1.d0,0.d0,500.d0)
        call bookup(l+ 5,'a pt         '
     &       //weights_info(kk)//cc(i),1.d0,0.d0,500.d0)
        call bookup(l+ 6,'l pt         '
     &       //weights_info(kk)//cc(i),1.d0,0.d0,500.d0)
        call bookup(l+ 7,'minv aa       '
     &       //weights_info(kk)//cc(i),0.02d0,115.d0,135.d0)
        call bookup(l+ 8,'pt miss       '
     &       //weights_info(kk)//cc(i),1.d0,0.d0,500.d0)
        call bookup(l+ 9,'j eta        '
     &       //weights_info(kk)//cc(i),0.01d0,0.d0,5.d0)
        call bookup(l+10,'b eta     '
     &       //weights_info(kk)//cc(i),0.01d0,0.d0,5.d0)
        call bookup(l+11,'a eta      '
     &       //weights_info(kk)//cc(i),0.01d0,0.d0,5.d0)
        call bookup(l+12,'l eta          '
     &       //weights_info(kk)//cc(i),0.01d0,0.d0,5.d0)
        call bookup(l+13,'aa del R       '
     &       //weights_info(kk)//cc(i),pi/20.d0,0.d0,3*pi)
        call bookup(l+14,'jj del R         '
     &       //weights_info(kk)//cc(i),pi/20.d0,0.d0,3*pi)
        call bookup(l+15,'bb del R    '
     &       //weights_info(kk)//cc(i),pi/20.d0,0.d0,3*pi)
        call bookup(l+16,'bj del R    '
     &       //weights_info(kk)//cc(i),pi/20.d0,0.d0,3*pi)
        call bookup(l+17,'al del R '
     &       //weights_info(kk)//cc(i),pi/20.d0,0.d0,3*pi)
        call bookup(l+18,'htt pt min            '
     &       //weights_info(kk)//cc(i),10.d0,0.d0,5000.d0)
        call bookup(l+19,'aj del R       '
     &       //weights_info(kk)//cc(i),pi/20.d0,0.d0,3*pi)
        call bookup(l+20,'ab del R             '
     &       //weights_info(kk)//cc(i),pi/20.d0,0.d0,3*pi)
        call bookup(l+21,'htt del R min '
     &       //weights_info(kk)//cc(i),pi/20.d0,0.d0,3*pi)
      enddo
      enddo
      return
      end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine analysis_end(xnorm)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer idim,n_plot
      parameter (idim=7,n_plot=21)
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
        call multitop(l+ 1,2,3,'t pt',' ','LOG')
        call multitop(l+ 2,2,3,'h pt',' ','LOG')
        call multitop(l+ 3,2,3,'j pt',' ','LOG')
        call multitop(l+ 4,2,3,'b pt',' ','LOG')
        call multitop(l+ 5,2,3,'a pt',' ','LOG')
        call multitop(l+ 6,2,3,'l pt',' ','LOG')
        call multitop(l+ 7,2,3,'minv aa',' ','LOG')
        call multitop(l+ 8,2,3,'pt miss',' ','LOG')
        call multitop(l+ 9,2,3,'j eta',' ','LOG')
        call multitop(l+10,2,3,'b eta',' ','LOG')
        call multitop(l+11,2,3,'a eta',' ','LOG')
        call multitop(l+12,2,3,'l eta',' ','LOG')
        call multitop(l+13,2,3,'aa del R',' ','LOG')
        call multitop(l+14,2,3,'jj del R',' ','LOG')
        call multitop(l+15,2,3,'bb del R',' ','LOG')
        call multitop(l+16,2,3,'bj del R',' ','LOG')
        call multitop(l+17,2,3,'al del R',' ','LOG')
        call multitop(l+18,2,3,'htt pt min',' ','LOG')
        call multitop(l+19,2,3,'aj del R',' ','LOG')
        call multitop(l+20,2,3,'ab del R',' ','LOG')
        call multitop(l+21,2,3,'htt del R min',' ','LOG')
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
      parameter (idim=7,n_plot=21)
C      common/plot_analy/idim,n_plot
      integer nexternal
      integer iPDG(nexternal)
      double precision p(0:4,nexternal)
      double precision wgts(*)
      integer ibody,order
      double precision wgt,var
      integer i,ii,kk,l,nwgt_analysis
      common/c_analysis/nwgt_analysis
      double precision www,yptqtb(0:3),ycut,ptcut
     $     ,ptq1,ptq2,ptg,yq1,yq2,etaq1,etaq2,azi,azinorm,qqm,dr,yqq,yh,
     $     xpta1(0:3),xpta2(0:3),ya1,ya2,drtamin,azita,drta,
     $     pta1,pta2,ptamin,etaamax,drmin,draa,
     $     drba,drbamin,drwa,drwamin,ywp,ywm,ybq,ybb
      double precision xpth(0:3),yptqtbth(0:3),yptqth(0:3),ptamax,
     $ yptbth(0:3),Ht,minvtot,minvth,minvtxh,pth,minvwpb,minvwmbx
      double precision toppt,hpt,jpt,bpt,apt,lpt,maa,ptmiss,jeta,beta,
     $     aeta,leta,aadelR,jjdelR,bbdelR,bjdelR,aldelR,pttthmin,
     $     ajdelR,abdelR,httdelRmin
      double precision xplep(nexternal-2,0:3),xpjet(nexternal-2,0:3),
     $     xpbjet(nexternal-2,0:3),xptau(nexternal-2,0:3),
     $     xpneu(nexternal-2,0:3),xpgam(nexternal-2,0:3),xptop(2,0:3)
      logical cutflag(idim-1),cutflag2,leplep(6)
      integer nlep,ntau,njet,nbjet,nneu,ngam,ntop
      double precision getrapidity,getpseudorap,getinvm,getdelphi,
     $     getdelR
      external getrapidity,getpseudorap,getinvm,getdelphi,getdelR
      logical isjet,isbjet,islep,istau,isleptau,isneu,isgam
      external isjet,isbjet,islep,istau,isleptau,isneu,isgam
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
      if (ipdg(nexternal-1).ne.22) then
          write (*,*) 'error #10 in analysis_fill: '/
     &        /'only for process "p p > t t~ > w+ w- b b~ a a"'
         stop 1
      endif
      if (ipdg(nexternal).ne.22) then
          write (*,*) 'error #11 in analysis_fill: '/
     &        /'only for process "p p > t t~ > w+ w- b b~ a a"'
         stop 1
      endif

      nlep=0 ! mu+,e+,mu-,e-                                                            
      ntau=0 ! ta+,ta-
      njet=0 ! g,u,d,s,c
      nbjet=0 ! b
      nneu=0 ! ve,vm,vt
      ngam=0 ! a
      ntop=0 ! top
      do ii=1,6
         leplep(ii)=.FALSE.
      enddo
      do ii=3,nexternal-2
         if(islep(ipdg(ii)))then
            nlep=nlep+1
            do i=0,3
               xplep(nlep,i)=p(i,ii)
            enddo
         endif
         if(istau(ipdg(ii)))then
            ntau=ntau+1
            do i=0,3
               xptau(ntau,i)=p(i,ii)
            enddo
         endif
         if(isjet(ipdg(ii)))then
            njet=njet+1
            do i=0,3
               xpjet(njet,i)=p(i,ii)
            enddo
         endif
         if(isbjet(ipdg(ii)))then
            nbjet=nbjet+1
            do i=0,3
               xpbjet(nbjet,i)=p(i,ii)
            enddo
         endif
         if(isneu(ipdg(ii)))then
            nneu=nneu+1
            do i=0,3
               xpneu(nneu,i)=p(i,ii)
            enddo
         endif
         if(isgam(ipdg(ii)))then
            ngam=ngam+1
            do i=0,3
               xpgam(ngam,i)=p(i,ii)
            enddo
         endif
         if(abs(ipdg(ii)).eq.6)then
            ntop=ntop+1
            do i=0,3
               xptop(ntop,i)=p(i,ii)
            enddo
 1       endif
      enddo
      nbjet=nbjet+2
      do i=0,3
         xpbjet(nbjet-1,i)=p(i,nexternal-1)
         xpbjet(nbjet,i)=p(i,nexternal)
         xpta1(i)=p(i,nexternal-1)
         xpta2(i)=p(i,nexternal)
         xpth(i)=xpta1(i)+xpta2(i)
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

      toppt=1d9
      do i=1,ntop
         toppt=min(toppt,dsqrt(xptop(i,1)**2+xptop(i,2)**2))
      enddo
      hpt=dsqrt(xpth(1)**2+xpth(2)**2)
      pttthmin=min(toppt,hpt)
      jpt=1d9
      do i=1,njet
         jpt=min(jpt,dsqrt(xpjet(i,1)**2+xpjet(i,2)**2))
      enddo
      bpt=1d9
      do i=1,nbjet
         bpt=min(bpt,dsqrt(xpbjet(i,1)**2+xpbjet(i,2)**2))
      enddo
      apt=1d9
      do i=1,ngam
         apt=min(apt,dsqrt(xpgam(i,1)**2+xpgam(i,2)**2))
      enddo
      lpt=1d9
      do i=1,nlep
         lpt=min(lpt,dsqrt(xplep(i,1)**2+xplep(i,2)**2))
      enddo
      do i=1,ntau
         lpt=min(lpt,dsqrt(xptau(i,1)**2+xptau(i,2)**2))
      enddo
      maa=getinvm(xpth(0),xpth(1),xpth(2),xpth(3))
      if(nneu.EQ.0)then
         ptmiss=0d0
      else
         yptqtb(0:3)=0d0
         do i=1,nneu
            do ii=0,3
               yptqtb(ii)=yptqtb(ii)+xpneu(i,ii)
            enddo
         enddo
         ptmiss=dsqrt(yptqtb(1)**2+yptqtb(2)**2)
      endif
      jeta=-1d0
      do i=1,njet
         jeta=max(jeta,
     $        abs(getpseudorap(xpjet(i,0),xpjet(i,1),
     $        xpjet(i,2),xpjet(i,3))))
      enddo
      beta=-1d0
      do i=1,nbjet
         beta=max(beta,
     $        abs(getpseudorap(xpbjet(i,0),xpbjet(i,1),
     $        xpbjet(i,2),xpbjet(i,3))))
      enddo
      aeta=-1d0
      do i=1,ngam
         aeta=max(aeta,
     $        abs(getpseudorap(xpgam(i,0),xpgam(i,1),
     $        xpgam(i,2),xpgam(i,3))))
      enddo
      leta=-1d0
      do i=1,nlep
         leta=max(leta,
     $        abs(getpseudorap(xplep(i,0),xplep(i,1),
     $        xplep(i,2),xplep(i,3))))
      enddo
      do i=1,ntau
         leta=max(leta,
     $        abs(getpseudorap(xptau(i,0),xptau(i,1),
     $        xptau(i,2),xptau(i,3))))
      enddo
      aadelR=1d9
      do i=1,ngam-1
         do ii=i+1,ngam
            aadelR=min(aadelR,getdelR(xpgam(i,0:3),xpgam(ii,0:3)))
         enddo
      enddo
      jjdelR=1d9
      do i=1,njet-1
         do ii=i+1,njet
            jjdelR=min(jjdelR,getdelR(xpjet(i,0:3),xpjet(ii,0:3)))
         enddo
      enddo
      bbdelR=1d9
      do i=1,nbjet-1
         do ii=i+1,nbjet
            bbdelR=min(bbdelR,getdelR(xpbjet(i,0:3),xpbjet(ii,0:3)))
         enddo
      enddo
      bjdelR=1d9
      do i=1,nbjet
         do ii=1,njet
            bjdelR=min(bjdelR,getdelR(xpbjet(i,0:3),xpjet(ii,0:3)))
         enddo
      enddo
      aldelR=1d9
      do i=1,ngam
         do ii=1,nlep
            aldelR=min(aldelR,getdelR(xpgam(i,0:3),xplep(ii,0:3)))
         enddo
         do ii=1,ntau
            aldelR=min(aldelR,getdelR(xpgam(i,0:3),xptau(ii,0:3)))
         enddo
      enddo
      ajdelR=1d9
      do i=1,ngam
         do ii=1,njet
            ajdelR=min(ajdelR,getdelR(xpgam(i,0:3),xpjet(ii,0:3)))
         enddo
      enddo
      abdelR=1d9
      do i=1,ngam
         do ii=1,nbjet
            abdelR=min(abdelR,getdelR(xpgam(i,0:3),xpbjet(ii,0:3)))
         enddo
      enddo
      httdelRmin=1d9
      do i=1,ntop
         httdelRmin=min(httdelRmin,getdelR(xpth(0:3),xptop(i,0:3)))
      enddo
      ! basic cuts
      ptamin=dsqrt(xpta1(1)**2+xpta1(2)**2)
      ptamax=max(ptamin,dsqrt(xpta2(1)**2+xpta2(2)**2))
      ptamin=min(ptamin,dsqrt(xpta2(1)**2+xpta2(2)**2))
      etaamax=dabs(getpseudorap(xpta1(0),xpta1(1),xpta1(2),xpta1(3)))
      etaamax=max(etaamax,abs(getpseudorap(xpta2(0),xpta2(1),
     $     xpta2(2),xpta2(3))))

      cutflag2=ptamin.GT.10d0.AND.etaamax.LT.5d0
     &     .AND.maa.GE.120d0.AND.maa.LE.130d0
      do i=1,idim-1
         cutflag(i)=cutflag2.AND.leplep(i)
      enddo

      do kk=1,nwgt_analysis
         www=wgts(kk)
         l=(kk-1)*idim*n_plot

         if(cutflag2)then

         call mfill(l+1,toppt,WWW)
         call mfill(l+2,hpt,WWW)
         call mfill(l+3,jpt,WWW)
         call mfill(l+4,bpt,WWW)
         call mfill(l+5,apt,WWW)
         call mfill(l+6,lpt,WWW)
         call mfill(l+7,maa,WWW)
         call mfill(l+8,ptmiss,WWW)
         call mfill(l+9,jeta,WWW)
         call mfill(l+10,beta,WWW)
         call mfill(l+11,aeta,WWW)
         call mfill(l+12,leta,WWW)
         call mfill(l+13,aadelR,WWW)
         call mfill(l+14,jjdelR,WWW)
         call mfill(l+15,bbdelR,WWW)
         call mfill(l+16,bjdelR,WWW)
         call mfill(l+17,aldelR,WWW)
         call mfill(l+18,pttthmin,WWW)
         call mfill(l+19,ajdelR,WWW)
         call mfill(l+20,abdelR,WWW)
         call mfill(l+21,httdelRmin,WWW)
         endif
c
c***************************************************** with cuts
c
         do ii=1,idim-1
         l=l+n_plot
c
         if(cutflag(ii))then
         call mfill(l+1,toppt,WWW)
         call mfill(l+2,hpt,WWW)
         call mfill(l+3,jpt,WWW)
         call mfill(l+4,bpt,WWW)
         call mfill(l+5,apt,WWW)
         call mfill(l+6,lpt,WWW)
         call mfill(l+7,maa,WWW)
         call mfill(l+8,ptmiss,WWW)
         call mfill(l+9,jeta,WWW)
         call mfill(l+10,beta,WWW)
         call mfill(l+11,aeta,WWW)
         call mfill(l+12,leta,WWW)
         call mfill(l+13,aadelR,WWW)
         call mfill(l+14,jjdelR,WWW)
         call mfill(l+15,bbdelR,WWW)
         call mfill(l+16,bjdelR,WWW)
         call mfill(l+17,aldelR,WWW)
         call mfill(l+18,pttthmin,WWW)
         call mfill(l+19,ajdelR,WWW)
         call mfill(l+20,abdelR,WWW)
         call mfill(l+21,httdelRmin,WWW)
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

      function isjet(ipdg)
      implicit none
      logical isjet
      integer ipdg
      if(abs(ipdg).le.4.or.ipdg.eq.21)then
         isjet=.TRUE.
      else
         isjet=.FALSE.
      endif
      end

      function isbjet(ipdg)
      implicit none
      logical isbjet
      integer ipdg
      if(abs(ipdg).eq.5)then
         isbjet=.TRUE.
      else
         isbjet=.FALSE.
      endif
      end

      function isgam(ipdg)
      implicit none
      logical isgam
      integer ipdg
      if(ipdg.eq.22)then
         isgam=.TRUE.
      else
         isgam=.FALSE.
      endif
      end
      
      function islep(ipdg)
      implicit none
      logical islep
      integer ipdg
      if(abs(ipdg).eq.11.or.abs(ipdg).eq.13)then
         islep=.TRUE.
      else
         islep=.FALSE.
      endif
      end

      function istau(ipdg)
      implicit none
      logical istau
      integer ipdg
      if(abs(ipdg).eq.15)then
         istau=.TRUE.
      else
         istau=.FALSE.
      endif
      end

      function isleptau(ipdg)
      implicit none
      logical isleptau
      integer ipdg
      logical istau,islep
      external istau,islep
      isleptau=istau(ipdg).OR.islep(ipdg)
      end

      function isneu(ipdg)
      implicit none
      logical isneu
      integer ipdg
      if(abs(ipdg).eq.12.or.abs(ipdg).eq.14.or.abs(ipdg).eq.16)then
         isneu=.TRUE.
      else
         isneu=.FALSE.
      endif
      end

      function getdelR(pmom1,pmom2)
      implicit none
      double precision pmom1(0:3),pmom2(0:3),getdelR
      double precision y1,y2,dphi
      double precision getrapidity,getdelphi
      external getrapidity,getdelphi
      dphi=getdelphi(pmom1(1),pmom1(2),pmom2(1),pmom2(2))
      y1=getrapidity(pmom1(0),pmom1(3))
      y2=getrapidity(pmom2(0),pmom2(3))
      getdelR=dsqrt(dphi**2+(y1-y2)**2)
      end

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
