c
c Example analysis for "p p > h t t~, h > a a" process.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine analysis_begin(nwgt,weights_info)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer idim,n_plot
      parameter (idim=2,n_plot=1)
C      common/plot_analy/idim,n_plot
      integer nwgt
      character*(*) weights_info(*)
      integer i,kk,l,nwgt_analysis
      common/c_analysis/nwgt_analysis
      character*5 cc(idim)
      data cc/'pp LO','cutsL'/
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
        call bookup(l+ 1,'maa            '
     &       //weights_info(kk)//cc(i),0.1d0,70.d0,170.d0)
      enddo
      enddo
      return
      end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine analysis_end(xnorm)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer idim,n_plot
      parameter (idim=2,n_plot=1)
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
        call multitop(l+ 1,2,3,'maa',' ','LOG')
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
      parameter (idim=2,n_plot=1)
C      common/plot_analy/idim,n_plot
      integer nexternal
      integer iPDG(nexternal)
      double precision p(0:4,nexternal)
      double precision wgts(*)
      integer ibody,order
      double precision wgt,var
      integer i,kk,l,nwgt_analysis
      common/c_analysis/nwgt_analysis
      double precision www,xptq(0:3),xptb(0:3),yptqtb(0:3),ycut,ptcut
     $     ,ptq1,ptq2,ptg,yq1,yq2,etaq1,etaq2,azi,azinorm,qqm,dr,yqq,yh,
     $     xpta1(0:3),xpta2(0:3),yptaa(0:3)
      double precision xpth(0:3),yptqtbth(0:3),yptqth(0:3),
     $ yptbth(0:3),Ht,minvtot,minvth,minvtxh,pth,minvaa
      logical siq1flag,siq2flag,siq3flag,cutflag,ggflag,uuflag,ddflag,
     . ugflag,dgflag
      double precision getrapidity,getpseudorap,getinvm,getdelphi
      external getrapidity,getpseudorap,getinvm,getdelphi
      double precision pi
      PARAMETER (PI=3.14159265358979312D0)
c      print*, "i pdg quando nn si blocca sono",ipdg(1),ipdg(2)
      if (nexternal.ne.7) then
         write (*,*) 'error #1 in analysis_fill: '/
     &        /'only for process "p p > h t t~, h > a a"'
         stop 1
      endif
      if (.not. (abs(ipdg(1)).le.5 .or. ipdg(1).eq.21)) then
         write (*,*) 'error #2 in analysis_fill: '/
     &        /'only for process "p p > h t t~, h > a a"'
         stop 1
      endif
      if (.not. (abs(ipdg(2)).le.5 .or. ipdg(2).eq.21)) then
         write (*,*) 'error #3 in analysis_fill: '/
     &        /'only for process "p p > h t t~, h > a a"'
         stop 1
      endif
      if (ipdg(3).ne.25) then
         write (*,*) 'error #4 in analysis_fill: '/
     &        /'only for process "p p > h t t~, h > a a"'
         stop 1
      endif
      if (ipdg(4).ne.6) then
         write (*,*) 'error #5 in analysis_fill: '/
     &        /'only for process "p p > h t t~, h > a a"'
         stop 1
      endif
      if (ipdg(5).ne.-6) then
         write (*,*) 'error #6 in analysis_fill: '/
     &        /'only for process "p p > h t t~, h > a a"'
         stop 1
      endif
      if (ipdg(6).ne.22) then
         write (*,*) 'error #7 in analysis_fill: '/
     &        /'only for process "p p > h t t~, h > a a"'
         stop 1
      endif
      if (ipdg(7).ne.22) then
         write (*,*) 'error #8 in analysis_fill: '/
     &        /'only for process "p p > h t t~, h > a a"'
         stop 1
      endif

      do i=0,3
        xpth(i)=p(i,3)
        xptq(i)=p(i,4)
        xptb(i)=p(i,5)
        xpta1(i)=p(i,6)
        xpta2(i)=p(i,7)
C        yptqtb(i)=xptq(i)+xptb(i)
C        yptqtbth(i)=xptq(i)+xptb(i)+xpth(i)
C        yptbth(i)=xptb(i)+xpth(i)
C        yptqth(i)=xptq(i)+xpth(i)
        yptaa(i)=xpta1(i)+xpta2(i)
      enddo
      PTCUT=200.D0
C
C      ptq1 = dsqrt(xptq(1)**2+xptq(2)**2)
C      ptq2 = dsqrt(xptb(1)**2+xptb(2)**2)
C      pth = dsqrt(xpth(1)**2+xpth(2)**2)
C      ptg = dsqrt(yptqtb(1)**2+yptqtb(2)**2)
C      yq1=getrapidity(xptq(0),xptq(3))
C      yq2=getrapidity(xptb(0),xptb(3))
C      yh=getrapidity(xpth(0),xpth(3))
C      etaq1=getpseudorap(xptq(0),xptq(1),xptq(2),xptq(3))
C      etaq2=getpseudorap(xptb(0),xptb(1),xptb(2),xptb(3))
C      azi=getdelphi(xptq(1),xptq(2),xptb(1),xptb(2))
C      azinorm = (pi-azi)/pi
C      qqm=getinvm(yptqtb(0),yptqtb(1),yptqtb(2),yptqtb(3))
C      dr  = dsqrt(azi**2+(etaq1-etaq2)**2)
C      yqq=getrapidity(yptqtb(0),yptqtb(3))
C      Ht=ptq1+ptq2+pth
C      minvtot=getinvm(yptqtbth(0),yptqtbth(1),yptqtbth(2),yptqtbth(3))
C      minvth=getinvm(yptqth(0),yptqth(1),yptqth(2),yptqth(3))
C      minvtxh=getinvm(yptbth(0),yptbth(1),yptbth(2),yptbth(3))
      minvaa=getinvm(yptaa(0),yptaa(1),yptaa(2),yptaa(3))
c-------------------------------------------------------------
C      siq1flag=ptq1.gt.ptcut!.and.abs(yq1).lt.ycut
C      siq2flag=ptq2.gt.ptcut!.and.abs(yq2).lt.ycut
C      siq3flag=pth.gt.ptcut
C      cutflag=siq1flag.and.siq2flag.and.siq3flag
      cutflag=minvaa.LT.100d0
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

         call mfill(l+1,minvaa,WWW)
c
c***************************************************** with cuts
c
         l=l+n_plot
c
         if(cutflag)then
         call mfill(l+1,minvaa,WWW)
         endif
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
        write(*,*)'Attempt to compute a negative mass'
        stop
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
