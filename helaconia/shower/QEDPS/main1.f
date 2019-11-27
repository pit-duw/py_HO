         Program main
   
* Test program for QEDPS
*        total cross section of ee --> ZH
*        comparison between QEDPS and the structure function
 
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
 
*--------------------------------------------------------------------
* total cross section of e+e- -> Z Higgs(for test)
*   xxx:CM-energy square
          sig_zh(xxx)=
     .      gzzh2*alpha/96.d0/amz**2*beta(amz**2/xxx,amh**2/xxx)
     .    *(ph**2+3*amz**2)/((xxx-amz**2)**2+amz**2*agz**2)
     .    *(1-4*xw+8*xw**2)/xw/(1-xw)
*---------------------------------------------------------------------
        lu=10
        open(lu,file='QEDPS.out',status='unknown')
  
* physical constants
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
* Mass of incident particles which annihilate
        ainmas=amel
* Random number seed (if negative, use default value)
        iseed=-999

        CALL QPINIT(q2max,ainmas,iseed)
 
* Q-square minimum
        q2min=(amh+amz)**2

* Number of total events
        NEVEN=50000

* Initialization of the corrected cross sections.
	sigqp=0
	sigsf=0

* Event DO-loop
        DO 1 I = 1, NEVEN
  
* QEDPS -----------------------------------------------------------
   
* Q2IN : (CM enegy)**2 before radiative photon emission.
* Q2OUT: (CM enegy)**2 after  radiative photon emission.
* Q2IN is NOT constant for future linear colliders due to beamstrahlung.
   
          q2in=q2max
          CALL QPGEN(q2in, q2out)
	  if(q2out.lt.q2min) goto 2

*  to look shower development ....
	  if(i.le.10) then
              write(6,*)' CM energy  =  ',sqrt(q2out),'  GeV'
              CALL QPDUMP
          end if
  
          eh=(q2out+amh**2-amz**2)/2.d0/sqrt(q2out)
	  ph=sqrt( (eh-amh)*(eh+amh) )
* answer 
          sigqp=sigqp+sig_zh(q2out)
   
2         continue
  
* structure function ------------------------------------------------
          al     = log(q2in/ainmas**2)
          abeta   = 2.d0*alpha/pi*(al-1.d0)
          xx1max = 1.d0-q2min/q2in
          xx1min = 0.d0
          yy1max = xx1max**abeta
          yy1    = yy1max*qprand(1)
          xx1    = yy1**(1.d0/abeta)
          radi   = hradi(xx1,abeta)
          radi   = radi*yy1max
          q2out=q2in*(1.d0-xx1)
          eh=(q2out+amh**2-amz**2)/2.d0/sqrt(q2out)
	  ph=sqrt( (eh-amh)*(eh+amh) )
* answer 
          sigsf=sigsf+sig_zh(q2out)*radi
  
*--------------------------------------------------------------------
  1     CONTINUE

          eh=(q2max+amh**2-amz**2)/2.d0/sqrt(q2max)
	  ph=sqrt( (eh-amh)*(eh+amh) )

* print out the results
        write(6,*)' total cross section (ee->zh) '
        write(6,*)'  w/o ISR   :',sig_zh(q2max)*gevpb  
        write(6,*)'  w/  QEDPS :',sigqp/neven*gevpb
        write(6,*)'  w/  SF    :',sigsf/neven*gevpb
        write(lu,*)' total cross section (ee->zh) '
        write(lu,*)'  w/o ISR   :',sig_zh(q2max)*gevpb  
        write(lu,*)'  w/  QEDPS :',sigqp/neven*gevpb
        write(lu,*)'  w/  SF    :',sigsf/neven*gevpb

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

*-----------------------------------------------------------------------
      function hradi(x1,abeta)
      implicit real* 8(a-h,o-z)
*--------------------------------------------
      pi2=acos(-1.d0)**2
      hradi =1.d0+0.75d0*abeta+abeta**2/4.d0*(9.d0/8.d0-pi2/3.d0)
     .+(-abeta*(1.d0-x1/2.d0)
     .+abeta*abeta/8.d0*(-4.d0*(2.d0-x1)*log(x1)-(1.d0+3.d0*(1.d0-x1)**2
     .)/x1*log(1.d0-x1)-6+x1)
     .)*x1**(1.d0-abeta)/abeta
      return
      end

