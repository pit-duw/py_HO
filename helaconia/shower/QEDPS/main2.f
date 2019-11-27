         Program main


* Test program for QEDPS
*        Differential cross section of ee --> ZH

         IMPLICIT REAL*8(A-H,O-Z)
         COMMON/QPLIST/PLPTN(10,1000),NLPTN(10,1000),NPTCL

*         NPTCL     :  Number of particles
*         PLPTN(1, ):  The lightcone fraction
*         PLPTN(2, ):  virtual mass squared
*         PLPTN(3, ):  Pt**2
*         PLPTN(4, ):  x-component of the four momentum
*         PLPTN(5, ):  y-component of the four momentum
*         PLPTN(6, ):  z-component of the four momentum
*         PLPTN(7, ):  E-component of the four momentum
*         PLPTN(8, ):  + component of the lightcone momentum
*         PLPTN(9, ):  - component of the lightcone momentum
*
*         NLPTN(1, ): particle ID
*                      electron:11, positron:-11
*                      photon :22
*                     (PDG standard)
*         NLPTN(2, ): relative address of the parent
*         NLPTN(3, ): number of the children
*         NLPTN(4, ): relative address of the first child
*         NLPTN(5, ): status of the process
*         NLPTN(6, ): spacelike(-1) or timelike(+1)
c--------------------------------------------------------------------
* Differential cross section of e+e- -> Z Higgs(for test)
*  xxx:CM-energy square
*  ccc:cos(theta),theta:polar angle of Higgs boson w.r.t. incident e-
          sig_zh(xxx,ccc)=
     .      gzzh2*alpha/4.d0/64.d0/pi/amz**2*beta(amz**2/xxx,amh**2/xxx)
     .    /((xxx-amz**2)**2+amz**2*agz**2)
     .    *(1-4*xw+8*xw**2)/xw/(1-xw)*(2*amz**2+ph**2*(1-ccc**2))
c---------------------------------------------------------------------
        lu=10
        open(lu,file='QEDPS.out',status='unknown')

*  physical constants
        alpha =  1.0D0/137.0359895D0
        pi    = acos(-1.d0)
        gevpb = 0.38937966d9

        amel  = 0.511d-3
        amz   = 91.1888D0
	agz   = 2.49D0
        amw   = 80.23D0
        amh   = amw

        xw    = 1-amw**2/amz**2
        gzzh2 = amz**2*4*pi*alpha/xw/(1-xw)

* Initialization
* Q-square maximum
        q2max=500.d0**2
* Mass of incident particles
        ainmas=amel
* Random number seed (if negative, use default value)
        iseed=-999

        CALL QPINIT(q2max,ainmas,iseed)

* Q-square minimum
        q2min=(amh+amz)**2

* angular cut on Higgs boson in the LAB-frame
	coscut=0.9d0

* Number of total events
        NEVEN=50000

* Initialization of the corrected cross sections.
        sigqp=0

* Event DO-loop
        DO 1 I=1,NEVEN
*  QEDPS -----------------------------------------------------------
* Q2IN : (CM enegy)**2 before radiative photon emission.
* Q2OUT: (CM enegy)**2 after  radiative photon emission.
*   Q2IN is NOT constant for future linear colliders due to beamstrahlung.
          q2in=q2max
          CALL QPGEN(q2in, q2out)
* q-square minimum check   
	  if(q2out.lt.q2min) goto 1

* set 4-momentum in CM-system after phton emission ***********************
* Higgs (PDG ID = 23)
          ID    = 23
	  cosh =QPRAND(1)*2-1
	  sinh =sqrt(1.d0-cosh**2)
	                      ajacob=2
          phih =QPRAND(1)*2*pi
	                      ajacob=ajacob*2*pi
          eh=(q2out+amh**2-amz**2)/2.d0/sqrt(q2out)
	  ph  = sqrt( (eh-amh)*(eh+amh) )
	  px_h= ph*sinh*cos(phih)
	  py_h= ph*sinh*sin(phih)
	  pz_h= ph*cosh

* IFRAG = 0 : first call of QPSET
	  IFRAG = 0

* set necessary information in 'COMMON/QPLIST/'
	  CALL QPSET(ID,eh,px_h,py_h,pz_h,IFLAG)

* Z  (PDG ID =25)
          ID    = 25
          ez  = (q2out-amh**2+amz**2)/2.d0/sqrt(q2out)
	  px_z=-px_h
	  py_z=-py_h
	  pz_z=-pz_h

* IFRAG = 1 : last call of QPSET
	  IFRAG = 1

* set necessary information in 'COMMON/QPLIST/' and boost them to LAB frame
	  CALL QPSET(ID,ez,px_z,py_z,pz_z,IFRAG)

* Cut on Higgs angle in LAB-frame (CM-system of incident beams)  
          do 2 j=1,nptcl
	   if(nlptn(1,j).eq.23) then
	    ph_lab=sqrt(plptn(4,j)**2+plptn(5,j)**2+plptn(6,j)**2)
	    coshlb=plptn(6,j)/ph_lab
	    if(abs(coshlb).gt.coscut) goto 1 
           end if
2         continue

* summing up cross section
          sigqp=sigqp+sig_zh(q2out,cosh)*ajacob
*--------------------------------------------------------------------
*  to look shower development ....
          if(i.le.10) then
              write(6,*)' CM energy  =  ',sqrt(q2out),'  GeV'
              CALL QPDUMP
          end if
  1     CONTINUE

        write(6,*)' total cross section (ee->zh) with cos cut= ',coscut
        write(6,*)'  w/  QEDPS :',sigqp/neven*gevpb
        write(lu,*)' total cross section (ee->zh) with cos cut= ',coscut
        write(lu,*)'  w/  QEDPS :',sigqp/neven*gevpb

        STOP
        END
*-----------------------------------------------------------------------
      function beta(z1,z2)
      implicit real* 8(a-h,o-z)
      if(1-2*(z1+z2)+(z1-z2)**2.lt.0)then
       beta = 0
      else
       beta=sqrt(1-2*(z1+z2)+(z1-z2)**2)
      end if
      return
      end

