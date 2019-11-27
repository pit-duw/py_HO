c
c Example analysis for "p p > psi + X" process.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine analysis_begin(nwgt,weights_info)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer nwgt
      character*(*) weights_info(*)
      integer i,kk,l,nwgt_analysis
      common/c_analysis/nwgt_analysis
      character*5 cc(2)
      data cc/'pp LO','cutsL'/
      double precision pi
      PARAMETER (PI=3.14159265358979312D0)
      include 'dbook.inc'
      call inihist
      nwgt_analysis=nwgt
      if (nwgt_analysis*40.gt.nplots/4) then
         write (*,*) 'error in analysis_begin: '/
     &        /'too many histograms, increase NPLOTS to',
     &        nwgt_analysis*40*4
         stop 1
      endif
      do kk=1,nwgt_analysis
      do i=1,2
        l=(kk-1)*4+(i-1)*10
        call bookup(l+ 1,'psi pt            '
     &       //weights_info(kk)//cc(i),0.2d0,0.d0,20.d0)
        call bookup(l+ 2,'psi y             '
     &       //weights_info(kk)//cc(i),0.1d0,0.d0,10d0)
        call bookup(l+ 3,'psi dpt dy [0,1]  '
     &       //weights_info(kk)//cc(i),0.2d0,0.d0,20.d0)
        call bookup(l+ 4,'psi dpt dy [1,2]  '
     &       //weights_info(kk)//cc(i),0.2d0,0.d0,20.d0)
        call bookup(l+ 5,'psi dpt dy [2,2.5]'
     &       //weights_info(kk)//cc(i),0.2d0,0.d0,20.d0)
        call bookup(l+ 6,'psi dpt dy [2.5,3]'
     &       //weights_info(kk)//cc(i),0.2d0,0.d0,20.d0)
        call bookup(l+ 7,'psi dpt dy [3,3.5]'
     &       //weights_info(kk)//cc(i),0.2d0,0.d0,20.d0)
        call bookup(l+ 8,'psi dpt dy [3.5,4]'
     &       //weights_info(kk)//cc(i),0.2d0,0.d0,20.d0)
        call bookup(l+ 9,'psi dpt dy [4,4.5]'
     &       //weights_info(kk)//cc(i),0.2d0,0.d0,20.d0)
        call bookup(l+ 10,'psi dpt dy [4.5,5]'
     &       //weights_info(kk)//cc(i),0.2d0,0.d0,20.d0)
      enddo
      enddo
      return
      end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine analysis_end(xnorm)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
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
      do i=1,2
        l=(kk-1)*4+(i-1)*10
        call multitop(l+ 1,2,3,'psi pt',' ','LOG')
        call multitop(l+ 2,2,3,'psi y',' ','LOG')
        call multitop(l+ 3,2,3,'psi dpt dy [0,1]',' ','LOG')
        call multitop(l+ 4,2,3,'psi dpt dy [1,2]',' ','LOG')
        call multitop(l+ 5,2,3,'psi dpt dy [2,2.5]',' ','LOG')
        call multitop(l+ 6,2,3,'psi dpt dy [2.5,3]',' ','LOG')
        call multitop(l+ 7,2,3,'psi dpt dy [3,3.5]',' ','LOG')
        call multitop(l+ 8,2,3,'psi dpt dy [3.5,4]',' ','LOG')
        call multitop(l+ 9,2,3,'psi dpt dy [4,4.5]',' ','LOG')
        call multitop(l+ 10,2,3,'psi dpt dy [4.5,5]',' ','LOG')
      enddo
      enddo
      call close_topdrawer_file
      return                
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine analysis_fill(p,nexternal,ipdg,wgts)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer nexternal
      integer iPDG(nexternal)
      double precision p(0:4,nexternal)
      double precision wgts(*)
      integer ibody,order
      double precision wgt,var
      integer i,kk,l,nwgt_analysis
      common/c_analysis/nwgt_analysis
      double precision www,xptq(0:3),xptb(0:3),yptqtb(0:3),ycut,ptcut
     $     ,ptq1,ptq2,ptg,yq1,yq2,etaq1,etaq2,azi,azinorm,qqm,dr,yqq,yh
      double precision xpth(0:3),yptqtbth(0:3),yptqth(0:3),
     $ yptbth(0:3),Ht,minvtot,minvth,minvtxh,pth
      logical siq1flag,siq2flag,siq3flag,cutflag,ggflag,uuflag,ddflag,
     . ugflag,dgflag
      double precision getrapidity,getpseudorap,getinvm,getdelphi
      external getrapidity,getpseudorap,getinvm,getdelphi
      double precision pi
      PARAMETER (PI=3.14159265358979312D0)
      if (nexternal.ne.4) then
         write (*,*) 'error #1 in analysis_fill: '/
     &        /'only for process "p p > psi + X"'
         stop 1
      endif
      if (.not. (abs(ipdg(1)).le.5 .or. ipdg(1).eq.21)) then
         write (*,*) 'error #2 in analysis_fill: '/
     &        /'only for process "p p > psi + X"'
         stop 1
      endif
      if (.not. (abs(ipdg(2)).le.5 .or. ipdg(2).eq.21)) then
         write (*,*) 'error #3 in analysis_fill: '/
     &        /'only for process "p p > psi + X"'
         stop 1
      endif
      if (ipdg(3).ne.443) then
         write (*,*) 'error #4 in analysis_fill: '/
     &        /'only for process "p p > psi + X"'
         stop 1
      endif
      if (.not.(abs(ipdg(4)).le.5.or.ipdg(4).eq.21)) then
         write (*,*) 'error #5 in analysis_fill: '/
     &        /'only for process "p p > psi + X"'
         stop 1
      endif
c
      do i=0,3
        xpth(i)=p(i,3)
        xptq(i)=p(i,4)
      enddo
      PTCUT=5D0
C
      pth = dsqrt(xpth(1)**2+xpth(2)**2)
      yh=getrapidity(xpth(0),xpth(3))
      yh=abs(yh)
c-------------------------------------------------------------
      siq3flag=pth.gt.ptcut
      cutflag=siq3flag
c-------------------------------------------------------------
c      double precision xpth(0:3),yptqtbth(0:3),yptqth(0:3),
c     $ yptbth(0:3),Ht,minvtot,minvth,minvtxh
c      double precision xpth(0:3),yptqtbth(0:3),yptqth(0:3),
c     $ yptbth(0:3),Ht,minvtot,minvth,minvtxh,pth

      do kk=1,nwgt_analysis
         www=wgts(kk)
         l=(kk-1)*4
      
         call mfill(l+1,pth,WWW)
         call mfill(l+2,yh,WWW)
         if(yh.le.1d0)then
            call mfill(l+3,pth,WWW/0.2d0/2d0)
         elseif(yh.gt.1d0.and.yh.le.2d0)then
            call mfill(l+4,pth,WWW/0.2d0/2d0)
         elseif(yh.gt.2d0.and.yh.le.2.5d0)then
            call mfill(l+5,pth,WWW/0.2d0/1d0)
         elseif(yh.gt.2.5d0.and.yh.le.3d0)then
            call mfill(l+6,pth,WWW/0.2d0/1d0)
         elseif(yh.gt.3d0.and.yh.le.3.5d0)then
            call mfill(l+7,pth,WWW/0.2d0/1d0)
         elseif(yh.gt.3.5d0.and.yh.le.4d0)then
            call mfill(l+8,pth,WWW/0.2d0/1d0)
         elseif(yh.gt.4d0.and.yh.le.4.5d0)then
            call mfill(l+9,pth,WWW/0.2d0/1d0)
         elseif(yh.gt.4.5d0.and.yh.le.5d0)then
            call mfill(l+10,pth,WWW/0.2d0/1d0)
         endif
c
c***************************************************** with cuts
c
         l=l+2
c
         if(cutflag)then
         call mfill(l+1,pth,WWW)
         call mfill(l+2,yh,WWW)
         if(yh.le.1d0)then
            call mfill(l+3,pth,WWW/0.2d0/2d0)
         elseif(yh.gt.1d0.and.yh.le.2d0)then
            call mfill(l+4,pth,WWW/0.2d0/2d0)
         elseif(yh.gt.2d0.and.yh.le.2.5d0)then
            call mfill(l+5,pth,WWW/0.2d0/1d0)
         elseif(yh.gt.2.5d0.and.yh.le.3d0)then
            call mfill(l+6,pth,WWW/0.2d0/1d0)
         elseif(yh.gt.3d0.and.yh.le.3.5d0)then
            call mfill(l+7,pth,WWW/0.2d0/1d0)
         elseif(yh.gt.3.5d0.and.yh.le.4d0)then
            call mfill(l+8,pth,WWW/0.2d0/1d0)
         elseif(yh.gt.4d0.and.yh.le.4.5d0)then
            call mfill(l+9,pth,WWW/0.2d0/1d0)
         elseif(yh.gt.4.5d0.and.yh.le.5d0)then
            call mfill(l+10,pth,WWW/0.2d0/1d0)
         endif
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
