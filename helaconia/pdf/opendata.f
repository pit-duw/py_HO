      Integer Function NextUnopen()
c********************************************************************
C     Returns an unallocated FORTRAN i/o unit.
c********************************************************************

      Logical EX
C
      Do 10 N = 10, 300
         INQUIRE (UNIT=N, OPENED=EX)
         If (.NOT. EX) then
            NextUnopen = N
            Return
         Endif
 10   Continue
      Stop ' There is no available I/O unit. '
C               *************************
      End



      subroutine OpenData(Tablefile)
c********************************************************************
c generic subroutine to open the table files in the right directories
c********************************************************************
      implicit none
c
      Character Tablefile*40,up*3,lib*4,cteq_dir*9,mrs_dir*8,
     $     tempname*100
      data up,lib,cteq_dir,mrs_dir/'../','lib/','pdf/cteq/',
     $     'pdf/mrs/'/
      Integer IU,NextUnopen,i
      External NextUnopen
      common/IU/IU
c
c--   start
c
      IU=NextUnopen()

c     first try in the current directory (for cluster use)
      tempname=Tablefile
      open(IU,file=tempname,status='old',ERR=10)
      return

 10   tempname=up//Tablefile
      open(IU,file=tempname,status='old',ERR=20)
      return

c     then try data directories
 20   tempname=cteq_dir//Tablefile
      open(IU,file=tempname,status='old',ERR=30)
      return


 30   tempname=mrs_dir//Tablefile
      open(IU,file=tempname,status='old',ERR=40)
      return

 40   tempname=lib//tempname
      open(IU,file=tempname,status='old',ERR=50)

 50   continue
      do i=0,6
         open(IU,file=tempname,status='old',ERR=60)
         return
 60      tempname=up//tempname
         if (i.eq.6)then
            write(*,*) 'Error: PDF file ',Tablefile,' not found'
            stop
         endif
      enddo

      print*,'table for the pdf NOT found!!!'
      
      return
      end

