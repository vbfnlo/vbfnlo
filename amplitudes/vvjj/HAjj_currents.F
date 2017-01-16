      subroutine Get_HAjj_couplings()

      IMPLICIT NONE

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/scales.inc"
#include "VBFNLO/utilities/coupl.inc"

      integer mu
      double complex GBF(2,3,3:4)
      
      COMMON /HAJJCOUPL/ GBF
      
c     get vector boson fermion couplings GBF(iVecBos,iflav)
      do mu=1,2
         GBF(mu,1,3)=GAU(mu)
         GBF(mu,1,4)=GAD(mu)
         GBF(mu,2,3)=GZU(mu)
         GBF(mu,2,4)=GZD(mu)
         GBF(mu,3,3)=GWF(mu)
         GBF(mu,3,4)=GWF(mu)
      enddo

      return
      end

      subroutine Get_HAjj_currents(pbar,vbar,nlo)

      IMPLICIT NONE

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/scales.inc"
#include "VBFNLO/utilities/coupl.inc"
      
      integer nlo, iflav,isigA,isig,isigF,mu,i
      double precision pbar(0:3,max_p), vbar(0:3,max_v), pA(0:3), ph(0:4),zero,
     $     c2o4pi,alsUpperc2o4pi,alsLowerc2o4pi,alsULcvirtc2o4pi,pIntBosu(0:4),
     $     pIntBosl(0:4),fac,scale(2),stpi

  
      DOUBLE PRECISION CLR,XM2,XMG,B,V,A
      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),V(4,5),A(4,5)

      double complex GBF(2,3,3:4),SpinorI(6,4,-1:1),SpinorO(6,4,-1:1),
     $     Photon(6,-1:1),Higgs(6),Curr12(6,3:4,-1:1,2:3),
     $     Curr21(6,3:4,-1:1,2:3),Curr34(6,3:4,-1:1,2:3),
     $     Curr43(6,3:4,-1:1,2:3),SpinorIA(6,4,3:4,-1:1,-1:1),
     $     SpinorOA(6,4,3:4,-1:1,-1:1),Curr1A2(6,3:4,-1:1,2,-1:1,2:3),
     $     Curr2A1(6,3:4,-1:1,2,-1:1,2:3),Curr3A4(6,3:4,-1:1,3:4,-1:1,2:3),
     $     Curr4A3(6,3:4,-1:1,3:4,-1:1,2:3),Curr12tbv(6,3:4,-1:1,-1:1),Curr21tbv(6,3:4,-1:1,-1:1),
     $     Curr34tbv(6,3:4,-1:1,-1:1),Curr43tbv(6,3:4,-1:1,-1:1),
     $     Curr12H(6,3:4,-1:1,2:3),Curr21H(6,3:4,-1:1,2:3),Curr34H(6,3:4,-1:1,2:3),
     $     Curr43H(6,3:4,-1:1,2:3)

      COMMON /HAJJCOUPL/ GBF

      COMMON /HAJJCURRENTS/ Curr12,Curr21,Curr34,Curr43,Curr1A2,Curr2A1,Curr3A4,
     $     Curr4A3,Curr12tbv,Curr21tbv,Curr34tbv,Curr43tbv,Higgs,fac


      COMMON /HAJJCURRENTSVIRT/ SpinorI,SpinorO,Photon,Curr12H,Curr21H,Curr34H,Curr43H,
     $     pIntBosu,pIntBosl,scale,alsUpperc2o4pi,alsLowerc2o4pi,alsULcvirtc2o4pi

      double precision dotrr
      EXTERNAL dotrr
      parameter(zero=0d0, c2o4pi=1d0/3d0/pi, stpi=16d0*pi)

      do mu=0,3
         ph(mu) = 0d0
         do i=1,n_v-1
            ph(mu) = ph(mu)+vbar(mu,i) ! Higgs momentum
         enddo
         pA(mu) = vbar(mu,n_v)
      enddo
      
      ph(4)=dotrr(ph,ph)

c     incoming and outgoing particles with dependency of their helicities:
c     SpinorI/O(1,particle#,sigma_particle)
      do isig=-1,1,2
         CALL IXXXXX(pbar(0,1),zero,isig, 1,SpinorI(1,1,isig))
         CALL OXXXXX(pbar(0,3),zero,isig, 1,SpinorO(1,2,isig))
         CALL IXXXXX(pbar(0,2),zero,isig, 1,SpinorI(1,3,isig))
         CALL OXXXXX(pbar(0,4),zero,isig, 1,SpinorO(1,4,isig))

         CALL OXXXXX(pbar(0,1),zero,isig,-1,SpinorO(1,1,isig))
         CALL IXXXXX(pbar(0,3),zero,isig,-1,SpinorI(1,2,isig))
         CALL OXXXXX(pbar(0,2),zero,isig,-1,SpinorO(1,3,isig))
         CALL IXXXXX(pbar(0,4),zero,isig,-1,SpinorI(1,4,isig))
         
         CALL VXXXXX(pA,zero,isig,1,Photon(1,isig))

      enddo

      CALL SXXXXX(ph,1,Higgs)

c     attatch photon to each of the external fermion legs and get
c     SpinorI/OA(1,iins,iflav,isigF,isigA) 
c     with isigF/A=helicity index of fermion/photon and iflav=fermion flavour
      do iflav=3,4
         do isigF=-1,1,2
            do isigA=-1,1,2
               CALL FVIXXX(SpinorI(1,1,isigF),Photon(1,isigA),GBF(1,1,iflav),zero,zero,SpinorIA(1,1,iflav,isigF,isigA))
               CALL FVIXXX(SpinorI(1,2,isigF),Photon(1,isigA),GBF(1,1,iflav),zero,zero,SpinorIA(1,2,iflav,isigF,isigA))
               CALL FVIXXX(SpinorI(1,3,isigF),Photon(1,isigA),GBF(1,1,iflav),zero,zero,SpinorIA(1,3,iflav,isigF,isigA))
               CALL FVIXXX(SpinorI(1,4,isigF),Photon(1,isigA),GBF(1,1,iflav),zero,zero,SpinorIA(1,4,iflav,isigF,isigA))
               CALL FVOXXX(SpinorO(1,1,isigF),Photon(1,isigA),GBF(1,1,iflav),zero,zero,SpinorOA(1,1,iflav,isigF,isigA))
               CALL FVOXXX(SpinorO(1,2,isigF),Photon(1,isigA),GBF(1,1,iflav),zero,zero,SpinorOA(1,2,iflav,isigF,isigA))
               CALL FVOXXX(SpinorO(1,3,isigF),Photon(1,isigA),GBF(1,1,iflav),zero,zero,SpinorOA(1,3,iflav,isigF,isigA))
               CALL FVOXXX(SpinorO(1,4,isigF),Photon(1,isigA),GBF(1,1,iflav),zero,zero,SpinorOA(1,4,iflav,isigF,isigA))
            enddo
         enddo
      enddo

c     unmodified currents for the no-emission lines: 
c     Curr..(1,in-flavor,isig,iVecBos)
c     CurrXY means fermion number flow from X to Y
c     first: Z Boson
      do iflav=3,4
         do isigF=-1,1,2   
            CALL JIOXXX(SpinorI(1,1,isigF),SpinorO(1,2,isigF),GBF(1,2,iflav),ZMASS,ZWIDTH,Curr12(1,iflav,isigF,2))
            CALL JIOXXX(SpinorI(1,2,isigF),SpinorO(1,1,isigF),GBF(1,2,iflav),ZMASS,ZWIDTH,Curr21(1,iflav,isigF,2))
            CALL JIOXXX(SpinorI(1,3,isigF),SpinorO(1,4,isigF),GBF(1,2,iflav),ZMASS,ZWIDTH,Curr34(1,iflav,isigF,2))
            CALL JIOXXX(SpinorI(1,4,isigF),SpinorO(1,3,isigF),GBF(1,2,iflav),ZMASS,ZWIDTH,Curr43(1,iflav,isigF,2))
         enddo
c     next: W Boson
         CALL JIOXXX(SpinorI(1,1,-1),SpinorO(1,2,-1),GBF(1,3,iflav),WMASS,WWIDTH,Curr12(1,iflav,-1,3))
         CALL JIOXXX(SpinorI(1,2, 1),SpinorO(1,1, 1),GBF(1,3,iflav),WMASS,WWIDTH,Curr21(1,iflav, 1,3))
         CALL JIOXXX(SpinorI(1,3,-1),SpinorO(1,4,-1),GBF(1,3,iflav),WMASS,WWIDTH,Curr34(1,iflav,-1,3))
         CALL JIOXXX(SpinorI(1,4, 1),SpinorO(1,3, 1),GBF(1,3,iflav),WMASS,WWIDTH,Curr43(1,iflav, 1,3))
         
         
c     currents for the photon emission lines: 
c     Curr..A(1,in-flavor,isigF,iinsA,isigA,iWeakBos)
c     CurrXYA means fermion number flow from X to Y           !needs SpinorI/OA(1,iins,iflav,isigF,isigA) 
c     first: Z Boson
         do isigA=-1,1,2
            do isigF=-1,1,2
               CALL JIOXXX(SpinorIA(1,1,iflav,isigF,isigA),SpinorO(1,2,isigF),GBF(1,2,iflav),
     $              ZMASS,ZWIDTH,Curr1A2(1,iflav,isigF,1,isigA,2))
               CALL JIOXXX(SpinorI(1,1,isigF),SpinorOA(1,2,iflav,isigF,isigA),GBF(1,2,iflav),
     $              ZMASS,ZWIDTH,Curr1A2(1,iflav,isigF,2,isigA,2))
               CALL JIOXXX(SpinorI(1,2,isigF),SpinorOA(1,1,iflav,isigF,isigA),GBF(1,2,iflav),
     $              ZMASS,ZWIDTH,Curr2A1(1,iflav,isigF,1,isigA,2))
               CALL JIOXXX(SpinorIA(1,2,iflav,isigF,isigA),SpinorO(1,1,isigF),GBF(1,2,iflav),
     $              ZMASS,ZWIDTH,Curr2A1(1,iflav,isigF,2,isigA,2))

               CALL JIOXXX(SpinorIA(1,3,iflav,isigF,isigA),SpinorO(1,4,isigF),GBF(1,2,iflav),
     $              ZMASS,ZWIDTH,Curr3A4(1,iflav,isigF,3,isigA,2))
               CALL JIOXXX(SpinorI(1,3,isigF),SpinorOA(1,4,iflav,isigF,isigA),GBF(1,2,iflav),
     $              ZMASS,ZWIDTH,Curr3A4(1,iflav,isigF,4,isigA,2))
               CALL JIOXXX(SpinorI(1,4,isigF),SpinorOA(1,3,iflav,isigF,isigA),GBF(1,2,iflav),
     $              ZMASS,ZWIDTH,Curr4A3(1,iflav,isigF,3,isigA,2))
               CALL JIOXXX(SpinorIA(1,4,iflav,isigF,isigA),SpinorO(1,3,isigF),GBF(1,2,iflav),
     $              ZMASS,ZWIDTH,Curr4A3(1,iflav,isigF,4,isigA,2))
               
            enddo
c     next: W Boson
            CALL JIOXXX(SpinorIA(1,1,  iflav,-1,isigA),SpinorO(1,2,-1),GBF(1,3,iflav),WMASS,WWIDTH,Curr1A2(1,iflav,-1,1,isigA,3))
            CALL JIOXXX(SpinorIA(1,2,7-iflav, 1,isigA),SpinorO(1,1, 1),GBF(1,3,iflav),WMASS,WWIDTH,Curr2A1(1,iflav, 1,2,isigA,3))
            CALL JIOXXX(SpinorIA(1,3,  iflav,-1,isigA),SpinorO(1,4,-1),GBF(1,3,iflav),WMASS,WWIDTH,Curr3A4(1,iflav,-1,3,isigA,3))
            CALL JIOXXX(SpinorIA(1,4,7-iflav, 1,isigA),SpinorO(1,3, 1),GBF(1,3,iflav),WMASS,WWIDTH,Curr4A3(1,iflav, 1,4,isigA,3))
            
            CALL JIOXXX(SpinorI(1,1,-1),SpinorOA(1,2,7-iflav,-1,isigA),GBF(1,3,iflav),WMASS,WWIDTH,Curr1A2(1,iflav,-1,2,isigA,3))
            CALL JIOXXX(SpinorI(1,2, 1),SpinorOA(1,1,  iflav, 1,isigA),GBF(1,3,iflav),WMASS,WWIDTH,Curr2A1(1,iflav, 1,1,isigA,3))
            CALL JIOXXX(SpinorI(1,3,-1),SpinorOA(1,4,7-iflav,-1,isigA),GBF(1,3,iflav),WMASS,WWIDTH,Curr3A4(1,iflav,-1,4,isigA,3))
            CALL JIOXXX(SpinorI(1,4, 1),SpinorOA(1,3,  iflav, 1,isigA),GBF(1,3,iflav),WMASS,WWIDTH,Curr4A3(1,iflav, 1,3,isigA,3))
         
c     currents with three boson vertex attatched            !needs: Curr..(1,in-flavor,isig,iVecBos=3)
c     Curr..tbv(1,in-flavor,isigF,isigA)
         
            CALL JVVXXX(Photon(1,isigA),Curr12(1,iflav,-1,3),( 7-2*iflav)*GWWA,WMASS,WWIDTH,Curr12tbv(1,iflav,-1,isigA))
            CALL JVVXXX(Photon(1,isigA),Curr21(1,iflav, 1,3),(-7+2*iflav)*GWWA,WMASS,WWIDTH,Curr21tbv(1,iflav, 1,isigA))
            CALL JVVXXX(Curr34(1,iflav,-1,3),Photon(1,isigA),(-7+2*iflav)*GWWA,WMASS,WWIDTH,Curr34tbv(1,iflav,-1,isigA))
            CALL JVVXXX(Curr43(1,iflav, 1,3),Photon(1,isigA),( 7-2*iflav)*GWWA,WMASS,WWIDTH,Curr43tbv(1,iflav, 1,isigA))
         enddo
      enddo

      if(nlo.eq.1)then
c     for the virtual contributions, we also need the currents with the Higgs Boson attatched
c     Curr..H(1,iflav,isigF,iWeakBos)
         do iflav=3,4
            do isigA=-1,1,2
               do isigF=-1,1,2
                  CALL JVSXXX(Curr12(1,iflav,isigF,2),Higgs,GZZH,ZMASS,ZWIDTH,Curr12H(1,iflav,isigF,2))
                  CALL JVSXXX(Curr21(1,iflav,isigF,2),Higgs,GZZH,ZMASS,ZWIDTH,Curr21H(1,iflav,isigF,2))
                  CALL JVSXXX(Curr34(1,iflav,isigF,2),Higgs,GZZH,ZMASS,ZWIDTH,Curr34H(1,iflav,isigF,2))
                  CALL JVSXXX(Curr43(1,iflav,isigF,2),Higgs,GZZH,ZMASS,ZWIDTH,Curr43H(1,iflav,isigF,2))
                  
                  CALL JVSXXX(Curr12(1,iflav,isigF,3),Higgs,GWWH,WMASS,WWIDTH,Curr12H(1,iflav,isigF,3))
                  CALL JVSXXX(Curr21(1,iflav,isigF,3),Higgs,GWWH,WMASS,WWIDTH,Curr21H(1,iflav,isigF,3))
                  CALL JVSXXX(Curr34(1,iflav,isigF,3),Higgs,GWWH,WMASS,WWIDTH,Curr34H(1,iflav,isigF,3))
                  CALL JVSXXX(Curr43(1,iflav,isigF,3),Higgs,GWWH,WMASS,WWIDTH,Curr43H(1,iflav,isigF,3))
               enddo
            enddo
         enddo

         alsUpperc2o4pi=als(1,1)*c2o4pi
         alsLowerc2o4pi=als(2,1)*c2o4pi
         alsULcvirtc2o4pi=(alsUpperc2o4pi+alsLowerc2o4pi)*cvirt

         scale(1)=2d0*dotrr(pbar(0,1),pbar(0,3))
         scale(2)=2d0*dotrr(pbar(0,2),pbar(0,4))

         pIntBosu(0)=-dreal(Curr34H(5,3,1,2))
         pIntBosu(1)=-dreal(Curr34H(6,3,1,2))
         pIntBosu(2)=-dimag(Curr34H(6,3,1,2))
         pIntBosu(3)=-dimag(Curr34H(5,3,1,2))
         pIntBosu(4)=dotrr(pIntBosu,pIntBosu)
         pIntBosl(0)=-dreal(Curr12H(5,3,1,2))
         pIntBosl(1)=-dreal(Curr12H(6,3,1,2))
         pIntBosl(2)=-dimag(Curr12H(6,3,1,2))
         pIntBosl(3)=-dimag(Curr12H(5,3,1,2))
         pIntBosl(4)=dotrr(pIntBosl,pIntBosl)

      endif

c now get the coupling*propagator factors for subprocess k, helicity
c combination isig and sum the contributions from the 4 helicity 
c combinations for all subprocesses. Multiply by the factor corresponding to
c    SUM_pol |D_H(p_6-p_5) * psi-bar(6) psi(5)|^2 * color factor for H->bb
c
c multiply propagator by 16 pi*M Gamma to compensate for the phase 
c space integral in LIPSN
c
c   Int dq^2/2pi * dOmega/(32 pi^2) * 1/(q^2-M*2)^2+(MGamma)^2) 
c
c   = 1/(16 pi MGamma)
c
c 
      if (n_v.eq.1+1) then ! only H+A
        fac = 1d0
      else
        fac = stpi*xmg(6) / ( (ph(4)-xm2(6))**2 + xmg(6)**2 )
      endif

      fac=fac*9d0

      return
      end


      subroutine Get_HAjj_subtr_currents(pbar,vbar,kin)

c     This subroutine is in principle the same as Get_HAjj_currents, but has an additional index in
c     its outputs to identify the dipole kinematic.

      IMPLICIT NONE

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/scales.inc"
#include "VBFNLO/utilities/coupl.inc"
      
      integer iflav,isigA,isig,isigF,mu,kin,i
      double precision pbar(0:3,max_p,max_kin), vbar(0:3,max_v), pA(0:3), ph(0:4),
     $     zero,fac,stpi

 
      DOUBLE PRECISION CLR,XM2,XMG,B,V,A
      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),V(4,5),A(4,5)

      double complex GBF(2,3,3:4),SpinorI(6,4,-1:1),SpinorO(6,4,-1:1),
     $     Photon(6,-1:1),Higgs(6),Curr12(6,2:3,3:4,-1:1,2:3),
     $     Curr21(6,2:3,3:4,-1:1,2:3),Curr34(6,2:3,3:4,-1:1,2:3),
     $     Curr43(6,2:3,3:4,-1:1,2:3),SpinorIA(6,4,3:4,-1:1,-1:1),
     $     SpinorOA(6,4,3:4,-1:1,-1:1),Curr1A2(6,2:3,3:4,-1:1,2,-1:1,2:3),
     $     Curr2A1(6,2:3,3:4,-1:1,2,-1:1,2:3),Curr3A4(6,2:3,3:4,-1:1,3:4,-1:1,2:3),
     $     Curr4A3(6,2:3,3:4,-1:1,3:4,-1:1,2:3),Curr12tbv(6,2:3,3:4,-1:1,-1:1),Curr21tbv(6,2:3,3:4,-1:1,-1:1),
     $     Curr34tbv(6,2:3,3:4,-1:1,-1:1),Curr43tbv(6,2:3,3:4,-1:1,-1:1)

      COMMON /HAJJCOUPL/ GBF

      COMMON /HAJJSUBTRCURRENTS/ Curr12,Curr21,Curr34,Curr43,Curr1A2,Curr2A1,Curr3A4,
     $     Curr4A3,Curr12tbv,Curr21tbv,Curr34tbv,Curr43tbv,Higgs,fac

      parameter(zero=0d0, stpi=16d0*pi)
      
      double precision dotrr
      EXTERNAL dotrr

      do mu=0,3
         ph(mu) = 0d0
         do i=1,n_v-1
            ph(mu) = ph(mu)+vbar(mu,i) ! Higgs momentum
         enddo
         pA(mu) = vbar(mu,n_v)
      enddo

      ph(4)=dotrr(ph,ph)
      
      CALL VXXXXX(pA,zero,-1,1,Photon(1,-1))
      CALL VXXXXX(pA,zero,1,1,Photon(1,1))

      CALL SXXXXX(ph,1,Higgs)
      

c     incoming and outgoing particles with dependency of their helicities:
c     SpinorI/O(1,particle#,sigma_particle)
      do isig=-1,1,2
         CALL IXXXXX(pbar(0,1,kin),zero,isig, 1,SpinorI(1,1,isig))
         CALL OXXXXX(pbar(0,3,kin),zero,isig, 1,SpinorO(1,2,isig)) !exchanged 2<->3 in particle IDs
         CALL IXXXXX(pbar(0,2,kin),zero,isig, 1,SpinorI(1,3,isig))
         CALL OXXXXX(pbar(0,4,kin),zero,isig, 1,SpinorO(1,4,isig))

         CALL OXXXXX(pbar(0,1,kin),zero,isig,-1,SpinorO(1,1,isig))
         CALL IXXXXX(pbar(0,3,kin),zero,isig,-1,SpinorI(1,2,isig))
         CALL OXXXXX(pbar(0,2,kin),zero,isig,-1,SpinorO(1,3,isig))
         CALL IXXXXX(pbar(0,4,kin),zero,isig,-1,SpinorI(1,4,isig))
      enddo

c     attatch photon to each of the external fermion legs and get
c     SpinorI/OA(1,iins,iflav,isigF,isigA) 
c     with isigF/A=helicity index of fermion/photon and iflav=fermion flavour
      do iflav=3,4
         do isigF=-1,1,2
            do isigA=-1,1,2
               CALL FVIXXX(SpinorI(1,1,isigF),Photon(1,isigA),GBF(1,1,iflav),zero,zero,SpinorIA(1,1,iflav,isigF,isigA))
               CALL FVIXXX(SpinorI(1,2,isigF),Photon(1,isigA),GBF(1,1,iflav),zero,zero,SpinorIA(1,2,iflav,isigF,isigA))
               CALL FVIXXX(SpinorI(1,3,isigF),Photon(1,isigA),GBF(1,1,iflav),zero,zero,SpinorIA(1,3,iflav,isigF,isigA))
               CALL FVIXXX(SpinorI(1,4,isigF),Photon(1,isigA),GBF(1,1,iflav),zero,zero,SpinorIA(1,4,iflav,isigF,isigA))
               CALL FVOXXX(SpinorO(1,1,isigF),Photon(1,isigA),GBF(1,1,iflav),zero,zero,SpinorOA(1,1,iflav,isigF,isigA))
               CALL FVOXXX(SpinorO(1,2,isigF),Photon(1,isigA),GBF(1,1,iflav),zero,zero,SpinorOA(1,2,iflav,isigF,isigA))
               CALL FVOXXX(SpinorO(1,3,isigF),Photon(1,isigA),GBF(1,1,iflav),zero,zero,SpinorOA(1,3,iflav,isigF,isigA))
               CALL FVOXXX(SpinorO(1,4,isigF),Photon(1,isigA),GBF(1,1,iflav),zero,zero,SpinorOA(1,4,iflav,isigF,isigA))
            enddo
         enddo
      enddo

c     unmodified currents for the no-emission lines: 
c     Curr..(1,in-flavor,isig,iVecBos)
c     CurrXY means fermion number flow from X to Y
c     first: Z Boson
      do iflav=3,4
         do isigF=-1,1,2   
            CALL JIOXXX(SpinorI(1,1,isigF),SpinorO(1,2,isigF),GBF(1,2,iflav),ZMASS,ZWIDTH,Curr12(1,kin,iflav,isigF,2))
            CALL JIOXXX(SpinorI(1,2,isigF),SpinorO(1,1,isigF),GBF(1,2,iflav),ZMASS,ZWIDTH,Curr21(1,kin,iflav,isigF,2))
            CALL JIOXXX(SpinorI(1,3,isigF),SpinorO(1,4,isigF),GBF(1,2,iflav),ZMASS,ZWIDTH,Curr34(1,kin,iflav,isigF,2))
            CALL JIOXXX(SpinorI(1,4,isigF),SpinorO(1,3,isigF),GBF(1,2,iflav),ZMASS,ZWIDTH,Curr43(1,kin,iflav,isigF,2))
         enddo
c     next: W Boson
         CALL JIOXXX(SpinorI(1,1,-1),SpinorO(1,2,-1),GBF(1,3,iflav),WMASS,WWIDTH,Curr12(1,kin,iflav,-1,3))
         CALL JIOXXX(SpinorI(1,2, 1),SpinorO(1,1, 1),GBF(1,3,iflav),WMASS,WWIDTH,Curr21(1,kin,iflav, 1,3))
         CALL JIOXXX(SpinorI(1,3,-1),SpinorO(1,4,-1),GBF(1,3,iflav),WMASS,WWIDTH,Curr34(1,kin,iflav,-1,3))
         CALL JIOXXX(SpinorI(1,4, 1),SpinorO(1,3, 1),GBF(1,3,iflav),WMASS,WWIDTH,Curr43(1,kin,iflav, 1,3))            
         
c     currents for the photon emission lines: 
c     Curr..A(1,in-flavor,isigF,iinsA,isigA,iWeakBos)
c     CurrXYA means fermion number flow from X to Y           !needs SpinorI/OA(1,iins,iflav,isigF,isigA) 
c     first: Z Boson
         do isigA=-1,1,2
            do isigF=-1,1,2
               CALL JIOXXX(SpinorIA(1,1,iflav,isigF,isigA),SpinorO(1,2,isigF),GBF(1,2,iflav),
     $              ZMASS,ZWIDTH,Curr1A2(1,kin,iflav,isigF,1,isigA,2))
               CALL JIOXXX(SpinorI(1,1,isigF),SpinorOA(1,2,iflav,isigF,isigA),GBF(1,2,iflav),
     $              ZMASS,ZWIDTH,Curr1A2(1,kin,iflav,isigF,2,isigA,2))
               CALL JIOXXX(SpinorI(1,2,isigF),SpinorOA(1,1,iflav,isigF,isigA),GBF(1,2,iflav),
     $              ZMASS,ZWIDTH,Curr2A1(1,kin,iflav,isigF,1,isigA,2))
               CALL JIOXXX(SpinorIA(1,2,iflav,isigF,isigA),SpinorO(1,1,isigF),GBF(1,2,iflav),
     $              ZMASS,ZWIDTH,Curr2A1(1,kin,iflav,isigF,2,isigA,2))

               CALL JIOXXX(SpinorIA(1,3,iflav,isigF,isigA),SpinorO(1,4,isigF),GBF(1,2,iflav),
     $              ZMASS,ZWIDTH,Curr3A4(1,kin,iflav,isigF,3,isigA,2))
               CALL JIOXXX(SpinorI(1,3,isigF),SpinorOA(1,4,iflav,isigF,isigA),GBF(1,2,iflav)
     $              ,ZMASS,ZWIDTH,Curr3A4(1,kin,iflav,isigF,4,isigA,2))
               CALL JIOXXX(SpinorI(1,4,isigF),SpinorOA(1,3,iflav,isigF,isigA),GBF(1,2,iflav),
     $              ZMASS,ZWIDTH,Curr4A3(1,kin,iflav,isigF,3,isigA,2))
               CALL JIOXXX(SpinorIA(1,4,iflav,isigF,isigA),SpinorO(1,3,isigF),GBF(1,2,iflav),
     $              ZMASS,ZWIDTH,Curr4A3(1,kin,iflav,isigF,4,isigA,2))
               
            enddo
c     next: W Boson
            CALL JIOXXX(SpinorIA(1,1,  iflav,-1,isigA),SpinorO(1,2,-1),GBF(1,3,iflav),
     $              WMASS,WWIDTH,Curr1A2(1,kin,iflav,-1,1,isigA,3))
            CALL JIOXXX(SpinorIA(1,2,7-iflav, 1,isigA),SpinorO(1,1, 1),GBF(1,3,iflav)
     $              ,WMASS,WWIDTH,Curr2A1(1,kin,iflav, 1,2,isigA,3))
            CALL JIOXXX(SpinorIA(1,3,  iflav,-1,isigA),SpinorO(1,4,-1),GBF(1,3,iflav),
     $              WMASS,WWIDTH,Curr3A4(1,kin,iflav,-1,3,isigA,3))
            CALL JIOXXX(SpinorIA(1,4,7-iflav, 1,isigA),SpinorO(1,3, 1),GBF(1,3,iflav),
     $              WMASS,WWIDTH,Curr4A3(1,kin,iflav, 1,4,isigA,3))
            
            CALL JIOXXX(SpinorI(1,1,-1),SpinorOA(1,2,7-iflav,-1,isigA),GBF(1,3,iflav),
     $              WMASS,WWIDTH,Curr1A2(1,kin,iflav,-1,2,isigA,3))
            CALL JIOXXX(SpinorI(1,2, 1),SpinorOA(1,1,  iflav, 1,isigA),GBF(1,3,iflav),
     $              WMASS,WWIDTH,Curr2A1(1,kin,iflav, 1,1,isigA,3))
            CALL JIOXXX(SpinorI(1,3,-1),SpinorOA(1,4,7-iflav,-1,isigA),GBF(1,3,iflav)
     $              ,WMASS,WWIDTH,Curr3A4(1,kin,iflav,-1,4,isigA,3))
            CALL JIOXXX(SpinorI(1,4, 1),SpinorOA(1,3,  iflav, 1,isigA),GBF(1,3,iflav),
     $              WMASS,WWIDTH,Curr4A3(1,kin,iflav, 1,3,isigA,3))
            
c     currents with three boson vertex attatched            !needs: Curr..(1,in-flavor,isig,iVecBos=3)
c     Curr..tbv(1,in-flavor,isigF,isigA)
            
            CALL JVVXXX(Photon(1,isigA),Curr12(1,kin,iflav,-1,3),( 7-2*iflav)*GWWA,WMASS,WWIDTH,Curr12tbv(1,kin,iflav,-1,isigA))
            CALL JVVXXX(Photon(1,isigA),Curr21(1,kin,iflav, 1,3),(-7+2*iflav)*GWWA,WMASS,WWIDTH,Curr21tbv(1,kin,iflav, 1,isigA))
            CALL JVVXXX(Curr34(1,kin,iflav,-1,3),Photon(1,isigA),(-7+2*iflav)*GWWA,WMASS,WWIDTH,Curr34tbv(1,kin,iflav,-1,isigA))
            CALL JVVXXX(Curr43(1,kin,iflav, 1,3),Photon(1,isigA),( 7-2*iflav)*GWWA,WMASS,WWIDTH,Curr43tbv(1,kin,iflav, 1,isigA))

         enddo
      enddo

c now get the coupling*propagator factors for subprocess k, helicity
c combination isig and sum the contributions from the 4 helicity 
c combinations for all subprocesses. Multiply by the factor corresponding to
c    SUM_pol |D_H(p_6-p_5) * psi-bar(6) psi(5)|^2 * color factor for H->bb
c


c
c multiply propagator by 16 pi*M Gamma to compensate for the phase 
c space integral in LIPSN
c
c   Int dq^2/2pi * dOmega/(32 pi^2) * 1/(q^2-M*2)^2+(MGamma)^2) 
c
c   = 1/(16 pi MGamma)
c
c 
      if (n_v.eq.1+1) then ! only H+A
        fac = 1d0
      else
        fac = stpi*xmg(6) / ( (ph(4)-xm2(6))**2 + xmg(6)**2 )
      endif
      
      fac=fac*9d0
      
      return
      end


      subroutine Get_HAjjj_currents(pbar,vbar)

      IMPLICIT NONE

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/scales.inc"
#include "VBFNLO/utilities/coupl.inc"
      
      integer iflav,isigA,isigB,isigG,isig,isigF,mu,i
      double precision pbar(0:3,max_p), vbar(0:3,max_v), pA(0:3), ph(0:4),zero,fpi,stpi,fac

      DOUBLE PRECISION CLR,XM2,XMG,B,V,A
      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),V(4,5),A(4,5)


      double complex GBF(2,3,3:4),GGupper(2),GGlower(2),
     $     SpinorI(6,5,-1:1),SpinorO(6,5,-1:1),
     $     Photon(6,-1:1),Higgs(6),GluonI(6,3,-1:1),GluonO(6,-1:1),
     $     SpinorIA(6,4,3:4,-1:1,-1:1),SpinorOA(6,5,3:4,-1:1,-1:1),
     $     SpinorOGO(6,4,-1:1,-1:1),SpinorIGO(6,4,-1:1,-1:1),
     $     SpinorOGI(6,1:3,-1:1,-1:1),SpinorIGI(6,2:4,-1:1,-1:1),
     $     SpinorOAGO(6,4,3:4,-1:1,-1:1,-1:1),
     $     SpinorIAGO(6,4,3:4,-1:1,-1:1,-1:1),
     $     SpinorOAGI(6,4,3:4,-1:1,-1:1,-1:1),
     $     SpinorIAGI(6,4,3:4,-1:1,-1:1,-1:1),
     $     SpinorOGOA(6,4,3:4,-1:1,-1:1,-1:1),
     $     SpinorIGOA(6,4,3:4,-1:1,-1:1,-1:1),
     $     SpinorOGIA(6,4,3:4,-1:1,-1:1,-1:1),
     $     SpinorIGIA(6,4,3:4,-1:1,-1:1,-1:1),
     $     Curr12(6,3:4,-1:1,2:3),Curr21(6,3:4,-1:1,2:3),
     $     Curr34(6,3:4,-1:1,2:3),Curr43(6,3:4,-1:1,2:3),
     $     Curr12A(6,3:4,-1:1,2,-1:1,2:3),Curr21A(6,3:4,-1:1,2,-1:1,2:3),
     $     Curr34A(6,3:4,-1:1,3:4,-1:1,2:3),Curr43A(6,3:4,-1:1,3:4,-1:1,2:3),
     $     Curr12GO(6,3:4,-1:1,4,-1:1,2:3),Curr21GO(6,3:4,-1:1,4,-1:1,2:3),
     $     Curr34GO(6,3:4,-1:1,4,-1:1,2:3),Curr43GO(6,3:4,-1:1,4,-1:1,2:3),
     $     Curr21GI(6,3:4,-1:1,4,-1:1,2:3),Curr43GI(6,3:4,-1:1,4,-1:1,2:3),
     $     Curr12GOt1Z(6,3:4,-1:1,-1:1,-1:1,2,2),
     $     Curr21GOt1Z(6,3:4,-1:1,-1:1,-1:1,2,2),
     $     Curr34GOt1Z(6,3:4,-1:1,-1:1,-1:1,3:4,3:4),
     $     Curr43GOt1Z(6,3:4,-1:1,-1:1,-1:1,3:4,3:4),
     $     Curr12GOt1W(6,3:4,-1:1,-1:1,-1:1,2,2),
     $     Curr21GOt1W(6,3:4,-1:1,-1:1,-1:1,2,2),
     $     Curr34GOt1W(6,3:4,-1:1,-1:1,-1:1,3:4,3:4),
     $     Curr43GOt1W(6,3:4,-1:1,-1:1,-1:1,3:4,3:4),
     $     Curr12GOt2Z(6,3:4,-1:1,-1:1,-1:1,2,2),
     $     Curr21GOt2Z(6,3:4,-1:1,-1:1,-1:1,2,2),
     $     Curr34GOt2Z(6,3:4,-1:1,-1:1,-1:1,3:4,2),
     $     Curr43GOt2Z(6,3:4,-1:1,-1:1,-1:1,3:4,2),
     $     Curr12GOt2W(6,3:4,-1:1,-1:1,-1:1,2,2),
     $     Curr21GOt2W(6,3:4,-1:1,-1:1,-1:1,2,2),
     $     Curr34GOt2W(6,3:4,-1:1,-1:1,-1:1,3:4,2),
     $     Curr43GOt2W(6,3:4,-1:1,-1:1,-1:1,3:4,2),
     $     Curr21GIt1Z(6,3:4,-1:1,-1:1,-1:1,2,2),
     $     Curr43GIt1Z(6,3:4,-1:1,-1:1,-1:1,3:4,3:4),
     $     Curr21GIt1W(6,3:4,-1:1,-1:1,-1:1,2,2),
     $     Curr43GIt1W(6,3:4,-1:1,-1:1,-1:1,3:4,3:4),
     $     Curr21GIt2Z(6,3:4,-1:1,-1:1,-1:1,2,2),
     $     Curr43GIt2Z(6,3:4,-1:1,-1:1,-1:1,3:4,2),
     $     Curr21GIt2W(6,3:4,-1:1,-1:1,-1:1,2,2),
     $     Curr43GIt2W(6,3:4,-1:1,-1:1,-1:1,3:4,2),
     $     Curr12tbv(6,3:4,-1:1,-1:1),Curr21tbv(6,3:4,-1:1,-1:1),
     $     Curr34tbv(6,3:4,-1:1,-1:1),Curr43tbv(6,3:4,-1:1,-1:1),
     $     Curr12GOtbv(6,3:4,-1:1,-1:1,2,-1:1),
     $     Curr21GOtbv(6,3:4,-1:1,-1:1,2,-1:1),
     $     Curr34GOtbv(6,3:4,-1:1,-1:1,3:4,-1:1),
     $     Curr43GOtbv(6,3:4,-1:1,-1:1,3:4,-1:1),
     $     Curr21GItbv(6,3:4,-1:1,-1:1,2,-1:1),
     $     Curr43GItbv(6,3:4,-1:1,-1:1,3:4,-1:1)
      
      COMMON /HAJJCOUPL/ GBF

      COMMON /HAJJJCURRENTS/ fac,Higgs,Curr12,Curr21,Curr34,Curr43,
     $     Curr12A,Curr21A,Curr34A,Curr43A,
     $     Curr12GO,Curr21GO,Curr34GO,Curr43GO,Curr21GI,Curr43GI,
     $     Curr12GOt1Z,Curr21GOt1Z,Curr34GOt1Z,Curr43GOt1Z,Curr21GIt1Z,Curr43GIt1Z,
     $     Curr12GOt2Z,Curr21GOt2Z,Curr34GOt2Z,Curr43GOt2Z,Curr21GIt2Z,Curr43GIt2Z,
     $     Curr12GOt1W,Curr21GOt1W,Curr34GOt1W,Curr43GOt1W,Curr21GIt1W,Curr43GIt1W,
     $     Curr12GOt2W,Curr21GOt2W,Curr34GOt2W,Curr43GOt2W,Curr21GIt2W,Curr43GIt2W,
     $     Curr12tbv,Curr21tbv,Curr34tbv,Curr43tbv,
     $     Curr12GOtbv,Curr21GOtbv,Curr34GOtbv,Curr43GOtbv,Curr21GItbv,Curr43GItbv
      
      double precision dotrr
      EXTERNAL dotrr

      parameter(zero=0d0)
      parameter(fpi=4d0*pi, stpi=16d0*pi)

      do mu=1,2
         GGupper(mu)=SQRT(fpi*als(1,1))
         GGlower(mu)=SQRT(fpi*als(2,1))
      enddo

      do mu=0,3
         ph(mu) = 0d0
         do i=1,n_v-1
            ph(mu) = ph(mu)+vbar(mu,i) ! Higgs momentum
         enddo
         pA(mu) = vbar(mu,n_v)
      enddo

      ph(4)=dotrr(ph,ph)

c     incoming and outgoing particles with dependency of their helicities:
c     SpinorI/O(1,particle#,sigma_particle)
      do isig=-1,1,2
         CALL IXXXXX(pbar(0,1),zero,isig, 1,SpinorI(1,1,isig))
         CALL OXXXXX(pbar(0,3),zero,isig, 1,SpinorO(1,2,isig))
         CALL IXXXXX(pbar(0,2),zero,isig, 1,SpinorI(1,3,isig))
         CALL OXXXXX(pbar(0,4),zero,isig, 1,SpinorO(1,4,isig))
         CALL OXXXXX(pbar(0,5),zero,-isig,1,SpinorO(1,5,isig))

         CALL OXXXXX(pbar(0,1),zero,isig,-1,SpinorO(1,1,isig))
         CALL IXXXXX(pbar(0,3),zero,isig,-1,SpinorI(1,2,isig))
         CALL OXXXXX(pbar(0,2),zero,isig,-1,SpinorO(1,3,isig))
         CALL IXXXXX(pbar(0,4),zero,isig,-1,SpinorI(1,4,isig))
         CALL IXXXXX(pbar(0,5),zero,-isig,1,SpinorI(1,5,isig))
         
         CALL VXXXXX(pA,zero,isig,1,Photon(1,isig))
         CALL VXXXXX(pbar(0,1),zero,isig,-1,GluonI(1,1,isig))
         CALL VXXXXX(pbar(0,2),zero,isig,-1,GluonI(1,3,isig))
         CALL VXXXXX(pbar(0,5),zero,isig, 1,GluonO(1,isig))
      enddo
      CALL SXXXXX(ph,1,Higgs)

c     attatch photon to each of the external fermion legs and get
c     SpinorI/OA(1,iins,iflav,isigF,isigA) 
c     with isigF/B=helicity index of fermion/photon and iflav=fermion flavour
      do isigF=-1,1,2
         do isigB=-1,1,2
            do iflav=3,4
               CALL FVIXXX(SpinorI(1,1,isigF),Photon(1,isigB),GBF(1,1,iflav),zero,zero,SpinorIA(1,1,iflav,isigF,isigB))
               CALL FVIXXX(SpinorI(1,2,isigF),Photon(1,isigB),GBF(1,1,iflav),zero,zero,SpinorIA(1,2,iflav,isigF,isigB))
               CALL FVIXXX(SpinorI(1,3,isigF),Photon(1,isigB),GBF(1,1,iflav),zero,zero,SpinorIA(1,3,iflav,isigF,isigB))
               CALL FVIXXX(SpinorI(1,4,isigF),Photon(1,isigB),GBF(1,1,iflav),zero,zero,SpinorIA(1,4,iflav,isigF,isigB))
               CALL FVOXXX(SpinorO(1,1,isigF),Photon(1,isigB),GBF(1,1,iflav),zero,zero,SpinorOA(1,1,iflav,isigF,isigB))
               CALL FVOXXX(SpinorO(1,2,isigF),Photon(1,isigB),GBF(1,1,iflav),zero,zero,SpinorOA(1,2,iflav,isigF,isigB))
               CALL FVOXXX(SpinorO(1,3,isigF),Photon(1,isigB),GBF(1,1,iflav),zero,zero,SpinorOA(1,3,iflav,isigF,isigB))
               CALL FVOXXX(SpinorO(1,4,isigF),Photon(1,isigB),GBF(1,1,iflav),zero,zero,SpinorOA(1,4,iflav,isigF,isigB))

               CALL FVOXXX(SpinorO(1,5,isigF),Photon(1,isigB),GBF(1,1,iflav),zero,zero,SpinorOA(1,5,iflav,isigF,isigB))
            enddo

c     attatch gluon to each of the external fermion legs and get
c     SpinorI/OGI/O(1,iins,isigF,isigA) 
c     with isigF/B=helicity index of fermion/gluon
            CALL FVIXXX(SpinorI(1,1,isigF),GluonO(1,isigB),GGupper,zero,zero,SpinorIGO(1,1,isigF,isigB))
            CALL FVIXXX(SpinorI(1,2,isigF),GluonO(1,isigB),GGupper,zero,zero,SpinorIGO(1,2,isigF,isigB))
            CALL FVIXXX(SpinorI(1,3,isigF),GluonO(1,isigB),GGlower,zero,zero,SpinorIGO(1,3,isigF,isigB))
            CALL FVIXXX(SpinorI(1,4,isigF),GluonO(1,isigB),GGlower,zero,zero,SpinorIGO(1,4,isigF,isigB))
            CALL FVOXXX(SpinorO(1,1,isigF),GluonO(1,isigB),GGupper,zero,zero,SpinorOGO(1,1,isigF,isigB))
            CALL FVOXXX(SpinorO(1,2,isigF),GluonO(1,isigB),GGupper,zero,zero,SpinorOGO(1,2,isigF,isigB))
            CALL FVOXXX(SpinorO(1,3,isigF),GluonO(1,isigB),GGlower,zero,zero,SpinorOGO(1,3,isigF,isigB))
            CALL FVOXXX(SpinorO(1,4,isigF),GluonO(1,isigB),GGlower,zero,zero,SpinorOGO(1,4,isigF,isigB))

            CALL FVOXXX(SpinorO(1,5,isigF),GluonI(1,1,isigB),GGupper,zero,zero,SpinorOGI(1,1,isigF,isigB))
            CALL FVOXXX(SpinorO(1,5,isigF),GluonI(1,3,isigB),GGlower,zero,zero,SpinorOGI(1,3,isigF,isigB))
            CALL FVIXXX(SpinorI(1,2,isigF),GluonI(1,1,isigB),GGupper,zero,zero,SpinorIGI(1,2,isigF,isigB))
            CALL FVIXXX(SpinorI(1,4,isigF),GluonI(1,3,isigB),GGlower,zero,zero,SpinorIGI(1,4,isigF,isigB))

c     attatch gluon to each external leg with photon emission, resulting in
c     SpinorAGI/GO(1,iins,iflav,isigF,isigA,isigG)

            do iflav=3,4
               do isigG=-1,1,2
                  CALL FVIXXX(SpinorIA(1,1,iflav,isigF,isigB),GluonO(1,isigG),GGupper,
     $              zero,zero,SpinorIAGO(1,1,iflav,isigF,isigB,isigG))
                  CALL FVIXXX(SpinorIA(1,2,iflav,isigF,isigB),GluonO(1,isigG),GGupper,
     $              zero,zero,SpinorIAGO(1,2,iflav,isigF,isigB,isigG))
                  CALL FVIXXX(SpinorIA(1,3,iflav,isigF,isigB),GluonO(1,isigG),GGlower,
     $              zero,zero,SpinorIAGO(1,3,iflav,isigF,isigB,isigG))
                  CALL FVIXXX(SpinorIA(1,4,iflav,isigF,isigB),GluonO(1,isigG),GGlower,
     $              zero,zero,SpinorIAGO(1,4,iflav,isigF,isigB,isigG))
                  CALL FVOXXX(SpinorOA(1,1,iflav,isigF,isigB),GluonO(1,isigG),GGupper,
     $              zero,zero,SpinorOAGO(1,1,iflav,isigF,isigB,isigG))
                  CALL FVOXXX(SpinorOA(1,2,iflav,isigF,isigB),GluonO(1,isigG),GGupper,
     $              zero,zero,SpinorOAGO(1,2,iflav,isigF,isigB,isigG))
                  CALL FVOXXX(SpinorOA(1,3,iflav,isigF,isigB),GluonO(1,isigG),GGlower,
     $              zero,zero,SpinorOAGO(1,3,iflav,isigF,isigB,isigG))
                  CALL FVOXXX(SpinorOA(1,4,iflav,isigF,isigB),GluonO(1,isigG),GGlower,
     $              zero,zero,SpinorOAGO(1,4,iflav,isigF,isigB,isigG))
                  
                  CALL FVOXXX(SpinorOA(1,5,iflav,isigF,isigB),GluonI(1,1,isigG),GGupper,
     $              zero,zero,SpinorOAGI(1,1,iflav,isigF,isigB,isigG))
                  CALL FVOXXX(SpinorOA(1,5,iflav,isigF,isigB),GluonI(1,3,isigG),GGlower,
     $              zero,zero,SpinorOAGI(1,3,iflav,isigF,isigB,isigG))
                  CALL FVIXXX(SpinorIA(1,2,iflav,isigF,isigB),GluonI(1,1,isigG),GGupper,
     $              zero,zero,SpinorIAGI(1,2,iflav,isigF,isigB,isigG))
                  CALL FVIXXX(SpinorIA(1,4,iflav,isigF,isigB),GluonI(1,3,isigG),GGlower,
     $              zero,zero,SpinorIAGI(1,4,iflav,isigF,isigB,isigG))
               enddo
c the same with photon attatched to a spinor with gluon emission
            
               do isigA=-1,1,2
                  CALL FVIXXX(SpinorIGO(1,1,isigF,isigB),Photon(1,isigA),GBF(1,1,iflav)
     $              ,zero,zero,SpinorIGOA(1,1,iflav,isigF,isigA,isigB))
                  CALL FVIXXX(SpinorIGO(1,2,isigF,isigB),Photon(1,isigA),GBF(1,1,iflav),
     $              zero,zero,SpinorIGOA(1,2,iflav,isigF,isigA,isigB))
                  CALL FVIXXX(SpinorIGO(1,3,isigF,isigB),Photon(1,isigA),GBF(1,1,iflav),
     $              zero,zero,SpinorIGOA(1,3,iflav,isigF,isigA,isigB))
                  CALL FVIXXX(SpinorIGO(1,4,isigF,isigB),Photon(1,isigA),GBF(1,1,iflav),
     $              zero,zero,SpinorIGOA(1,4,iflav,isigF,isigA,isigB))
                  CALL FVOXXX(SpinorOGO(1,1,isigF,isigB),Photon(1,isigA),GBF(1,1,iflav)
     $              ,zero,zero,SpinorOGOA(1,1,iflav,isigF,isigA,isigB))
                  CALL FVOXXX(SpinorOGO(1,2,isigF,isigB),Photon(1,isigA),GBF(1,1,iflav),
     $              zero,zero,SpinorOGOA(1,2,iflav,isigF,isigA,isigB))
                  CALL FVOXXX(SpinorOGO(1,3,isigF,isigB),Photon(1,isigA),GBF(1,1,iflav),
     $              zero,zero,SpinorOGOA(1,3,iflav,isigF,isigA,isigB))
                  CALL FVOXXX(SpinorOGO(1,4,isigF,isigB),Photon(1,isigA),GBF(1,1,iflav),
     $              zero,zero,SpinorOGOA(1,4,iflav,isigF,isigA,isigB))
                  
                  CALL FVOXXX(SpinorOGI(1,1,isigF,isigB),Photon(1,isigA),GBF(1,1,iflav),
     $              zero,zero,SpinorOGIA(1,1,iflav,isigF,isigA,isigB))
                  CALL FVOXXX(SpinorOGI(1,3,isigF,isigB),Photon(1,isigA),GBF(1,1,iflav)
     $              ,zero,zero,SpinorOGIA(1,3,iflav,isigF,isigA,isigB))
                  CALL FVIXXX(SpinorIGI(1,2,isigF,isigB),Photon(1,isigA),GBF(1,1,iflav),
     $              zero,zero,SpinorIGIA(1,2,iflav,isigF,isigA,isigB))
                  CALL FVIXXX(SpinorIGI(1,4,isigF,isigB),Photon(1,isigA),GBF(1,1,iflav),
     $              zero,zero,SpinorIGIA(1,4,iflav,isigF,isigA,isigB))
               enddo
            enddo
         enddo
      enddo

c     unmodified currents for the no-emission lines: 
c     Curr..(1,in-flavor,isig,iVecBos)
c     CurrXY means fermion number flow from X to Y
c     first: Z Boson
      do iflav=3,4
         do isigF=-1,1,2   
            CALL JIOXXX(SpinorI(1,1,isigF),SpinorO(1,2,isigF),GBF(1,2,iflav),ZMASS,ZWIDTH,Curr12(1,iflav,isigF,2))
            CALL JIOXXX(SpinorI(1,2,isigF),SpinorO(1,1,isigF),GBF(1,2,iflav),ZMASS,ZWIDTH,Curr21(1,iflav,isigF,2))
            CALL JIOXXX(SpinorI(1,3,isigF),SpinorO(1,4,isigF),GBF(1,2,iflav),ZMASS,ZWIDTH,Curr34(1,iflav,isigF,2))
            CALL JIOXXX(SpinorI(1,4,isigF),SpinorO(1,3,isigF),GBF(1,2,iflav),ZMASS,ZWIDTH,Curr43(1,iflav,isigF,2))
         enddo
c     next: W Boson
         CALL JIOXXX(SpinorI(1,1,-1),SpinorO(1,2,-1),GBF(1,3,iflav),WMASS,WWIDTH,Curr12(1,iflav,-1,3))
         CALL JIOXXX(SpinorI(1,2, 1),SpinorO(1,1, 1),GBF(1,3,iflav),WMASS,WWIDTH,Curr21(1,iflav, 1,3))
         CALL JIOXXX(SpinorI(1,3,-1),SpinorO(1,4,-1),GBF(1,3,iflav),WMASS,WWIDTH,Curr34(1,iflav,-1,3))
         CALL JIOXXX(SpinorI(1,4, 1),SpinorO(1,3, 1),GBF(1,3,iflav),WMASS,WWIDTH,Curr43(1,iflav, 1,3))
         
         
c     currents for the photon emission lines: 
c     Curr..A(1,in-flavor,isigF,iinsA,isigA,iWeakBos)
c     CurrXYA means fermion number flow from X to Y           !needs SpinorI/OA(1,iins,iflav,isigF,isigA) 
c     first: Z Boson
         do isigB=-1,1,2
            do isigF=-1,1,2
               CALL JIOXXX(SpinorIA(1,1,iflav,isigF,isigB),SpinorO(1,2,isigF),GBF(1,2,iflav),
     $              ZMASS,ZWIDTH,Curr12A(1,iflav,isigF,1,isigB,2))
               CALL JIOXXX(SpinorI(1,1,isigF),SpinorOA(1,2,iflav,isigF,isigB),GBF(1,2,iflav),
     $              ZMASS,ZWIDTH,Curr12A(1,iflav,isigF,2,isigB,2))
               CALL JIOXXX(SpinorI(1,2,isigF),SpinorOA(1,1,iflav,isigF,isigB),GBF(1,2,iflav),
     $              ZMASS,ZWIDTH,Curr21A(1,iflav,isigF,1,isigB,2))
               CALL JIOXXX(SpinorIA(1,2,iflav,isigF,isigB),SpinorO(1,1,isigF),GBF(1,2,iflav),
     $              ZMASS,ZWIDTH,Curr21A(1,iflav,isigF,2,isigB,2))

               CALL JIOXXX(SpinorIA(1,3,iflav,isigF,isigB),SpinorO(1,4,isigF),GBF(1,2,iflav)
     $              ,ZMASS,ZWIDTH,Curr34A(1,iflav,isigF,3,isigB,2))
               CALL JIOXXX(SpinorI(1,3,isigF),SpinorOA(1,4,iflav,isigF,isigB),GBF(1,2,iflav),
     $              ZMASS,ZWIDTH,Curr34A(1,iflav,isigF,4,isigB,2))
               CALL JIOXXX(SpinorI(1,4,isigF),SpinorOA(1,3,iflav,isigF,isigB),GBF(1,2,iflav)
     $              ,ZMASS,ZWIDTH,Curr43A(1,iflav,isigF,3,isigB,2))
               CALL JIOXXX(SpinorIA(1,4,iflav,isigF,isigB),SpinorO(1,3,isigF),GBF(1,2,iflav)
     $              ,ZMASS,ZWIDTH,Curr43A(1,iflav,isigF,4,isigB,2))
               
            enddo
c     next: W Boson
            CALL JIOXXX(SpinorIA(1,1,  iflav,-1,isigB),SpinorO(1,2,-1),GBF(1,3,iflav),WMASS,WWIDTH,Curr12A(1,iflav,-1,1,isigB,3))
            CALL JIOXXX(SpinorIA(1,2,7-iflav, 1,isigB),SpinorO(1,1, 1),GBF(1,3,iflav),WMASS,WWIDTH,Curr21A(1,iflav, 1,2,isigB,3))
            CALL JIOXXX(SpinorIA(1,3,  iflav,-1,isigB),SpinorO(1,4,-1),GBF(1,3,iflav),WMASS,WWIDTH,Curr34A(1,iflav,-1,3,isigB,3))
            CALL JIOXXX(SpinorIA(1,4,7-iflav, 1,isigB),SpinorO(1,3, 1),GBF(1,3,iflav),WMASS,WWIDTH,Curr43A(1,iflav, 1,4,isigB,3))
            
            CALL JIOXXX(SpinorI(1,1,-1),SpinorOA(1,2,7-iflav,-1,isigB),GBF(1,3,iflav),WMASS,WWIDTH,Curr12A(1,iflav,-1,2,isigB,3))
            CALL JIOXXX(SpinorI(1,2, 1),SpinorOA(1,1,  iflav, 1,isigB),GBF(1,3,iflav),WMASS,WWIDTH,Curr21A(1,iflav, 1,1,isigB,3))
            CALL JIOXXX(SpinorI(1,3,-1),SpinorOA(1,4,7-iflav,-1,isigB),GBF(1,3,iflav),WMASS,WWIDTH,Curr34A(1,iflav,-1,4,isigB,3))
            CALL JIOXXX(SpinorI(1,4, 1),SpinorOA(1,3,  iflav, 1,isigB),GBF(1,3,iflav),WMASS,WWIDTH,Curr43A(1,iflav, 1,3,isigB,3))

c     currents for the gluon emission lines: 
c     Curr...(1,in-flavor,isigF,iinsG,isigG,iWeakBos)
c     CurrXYZ means fermion number flow from X to Y to Z where one of the three indices denotes a gluon
c     first: Z Boson
            do isigF=-1,1,2
               CALL JIOXXX(SpinorIGO(1,1,isigF,isigB),SpinorO(1,2,isigF),GBF(1,2,iflav),
     $              ZMASS,ZWIDTH,Curr12GO(1,iflav,isigF,1,isigB,2))
               CALL JIOXXX(SpinorI(1,1,isigF),SpinorOGO(1,2,isigF,isigB),GBF(1,2,iflav),
     $              ZMASS,ZWIDTH,Curr12GO(1,iflav,isigF,2,isigB,2))
               CALL JIOXXX(SpinorI(1,2,isigF),SpinorOGO(1,1,isigF,isigB),GBF(1,2,iflav),
     $              ZMASS,ZWIDTH,Curr21GO(1,iflav,isigF,1,isigB,2))
               CALL JIOXXX(SpinorIGO(1,2,isigF,isigB),SpinorO(1,1,isigF),GBF(1,2,iflav),
     $              ZMASS,ZWIDTH,Curr21GO(1,iflav,isigF,2,isigB,2))

               CALL JIOXXX(SpinorIGO(1,3,isigF,isigB),SpinorO(1,4,isigF),GBF(1,2,iflav),
     $              ZMASS,ZWIDTH,Curr34GO(1,iflav,isigF,3,isigB,2))
               CALL JIOXXX(SpinorI(1,3,isigF),SpinorOGO(1,4,isigF,isigB),GBF(1,2,iflav),
     $              ZMASS,ZWIDTH,Curr34GO(1,iflav,isigF,4,isigB,2))
               CALL JIOXXX(SpinorI(1,4,isigF),SpinorOGO(1,3,isigF,isigB),GBF(1,2,iflav),
     $              ZMASS,ZWIDTH,Curr43GO(1,iflav,isigF,3,isigB,2))
               CALL JIOXXX(SpinorIGO(1,4,isigF,isigB),SpinorO(1,3,isigF),GBF(1,2,iflav),
     $              ZMASS,ZWIDTH,Curr43GO(1,iflav,isigF,4,isigB,2))
               
               CALL JIOXXX(SpinorI(1,2,isigF),SpinorOGI(1,1,isigF,isigB),GBF(1,2,iflav)
     $              ,ZMASS,ZWIDTH,Curr21GI(1,iflav,isigF,1,isigB,2))
               CALL JIOXXX(SpinorIGI(1,2,isigF,isigB),SpinorO(1,5,isigF),GBF(1,2,iflav),
     $              ZMASS,ZWIDTH,Curr21GI(1,iflav,isigF,2,isigB,2))
               CALL JIOXXX(SpinorI(1,4,isigF),SpinorOGI(1,3,isigF,isigB),GBF(1,2,iflav),
     $              ZMASS,ZWIDTH,Curr43GI(1,iflav,isigF,3,isigB,2))
               CALL JIOXXX(SpinorIGI(1,4,isigF,isigB),SpinorO(1,5,isigF),GBF(1,2,iflav),
     $              ZMASS,ZWIDTH,Curr43GI(1,iflav,isigF,4,isigB,2))


            enddo
c     next: W Boson
            CALL JIOXXX(SpinorIGO(1,1,-1,isigB),SpinorO(1,2,-1),GBF(1,3,iflav),WMASS,WWIDTH,Curr12GO(1,iflav,-1,1,isigB,3))
            CALL JIOXXX(SpinorIGO(1,2, 1,isigB),SpinorO(1,1, 1),GBF(1,3,iflav),WMASS,WWIDTH,Curr21GO(1,iflav, 1,2,isigB,3))
            CALL JIOXXX(SpinorIGO(1,3,-1,isigB),SpinorO(1,4,-1),GBF(1,3,iflav),WMASS,WWIDTH,Curr34GO(1,iflav,-1,3,isigB,3))
            CALL JIOXXX(SpinorIGO(1,4, 1,isigB),SpinorO(1,3, 1),GBF(1,3,iflav),WMASS,WWIDTH,Curr43GO(1,iflav, 1,4,isigB,3))
            
            CALL JIOXXX(SpinorI(1,1,-1),SpinorOGO(1,2,-1,isigB),GBF(1,3,iflav),WMASS,WWIDTH,Curr12GO(1,iflav,-1,2,isigB,3))
            CALL JIOXXX(SpinorI(1,2, 1),SpinorOGO(1,1, 1,isigB),GBF(1,3,iflav),WMASS,WWIDTH,Curr21GO(1,iflav, 1,1,isigB,3))
            CALL JIOXXX(SpinorI(1,3,-1),SpinorOGO(1,4,-1,isigB),GBF(1,3,iflav),WMASS,WWIDTH,Curr34GO(1,iflav,-1,4,isigB,3))
            CALL JIOXXX(SpinorI(1,4, 1),SpinorOGO(1,3, 1,isigB),GBF(1,3,iflav),WMASS,WWIDTH,Curr43GO(1,iflav, 1,3,isigB,3))

            CALL JIOXXX(SpinorI(1,2, 1),SpinorOGI(1,1, 1,isigB),GBF(1,3,iflav),WMASS,WWIDTH,Curr21GI(1,iflav, 1,1,isigB,3))
            CALL JIOXXX(SpinorIGI(1,2, 1,isigB),SpinorO(1,5, 1),GBF(1,3,iflav),WMASS,WWIDTH,Curr21GI(1,iflav, 1,2,isigB,3))
            CALL JIOXXX(SpinorI(1,4, 1),SpinorOGI(1,3, 1,isigB),GBF(1,3,iflav),WMASS,WWIDTH,Curr43GI(1,iflav, 1,3,isigB,3))
            CALL JIOXXX(SpinorIGI(1,4, 1,isigB),SpinorO(1,5, 1),GBF(1,3,iflav),WMASS,WWIDTH,Curr43GI(1,iflav, 1,4,isigB,3))
         enddo

c     Get currents with photon and gluon emissions.
c     There are two types of those:
c     -one emission on each side of the FFB vertex which emits the weak boson (Curr21t1)
c     -two emissions on one side of the FFB vertex which emits the weak boson (Curr21t2)
c     The indices are labeled in the following way:
c     Curr..t1(GI/GO)(WeakBos)(1,iflav-in,isigF,isigA,isigG,iinsA,iinsG)
c     Curr..t2(GI/GO)(WeakBos)(1,iflav-in,isigF,isigA,isigG,iins,iorder)
c     where 
c     iorder=1 means a call JIO(Spinor AG somewhere) and
c     iorder=2 means a call JIO(Spinor GA somewhere)

         do isigA=-1,1,2
            do isigG=-1,1,2
c     first: Z Boson
               do isigF=-1,1,2
c     upper current
                  CALL JIOXXX(SpinorIA(1,1,iflav,isigF,isigA),SpinorOGO(1,2,isigF,isigG),GBF(1,2,iflav)
     $              ,ZMASS,ZWIDTH,Curr12GOt1Z(1,iflav,isigF,isigA,isigG,1,2))
                  CALL JIOXXX(SpinorIAGO(1,1,iflav,isigF,isigA,isigG),SpinorO(1,2,isigF),GBF(1,2,iflav)
     $              ,ZMASS,ZWIDTH,Curr12GOt2Z(1,iflav,isigF,isigA,isigG,1,1))
                  CALL JIOXXX(SpinorI(1,1,isigF),SpinorOAGO(1,2,iflav,isigF,isigA,isigG),GBF(1,2,iflav)
     $              ,ZMASS,ZWIDTH,Curr12GOt2Z(1,iflav,isigF,isigA,isigG,2,1))

                  CALL JIOXXX(SpinorIGO(1,1,isigF,isigG),SpinorOA(1,2,iflav,isigF,isigA),GBF(1,2,iflav),
     $              ZMASS,ZWIDTH,Curr12GOt1Z(1,iflav,isigF,isigA,isigG,2,1))
                  CALL JIOXXX(SpinorIGOA(1,1,iflav,isigF,isigA,isigG),SpinorO(1,2,isigF),GBF(1,2,iflav),
     $              ZMASS,ZWIDTH,Curr12GOt2Z(1,iflav,isigF,isigA,isigG,1,2))
                  CALL JIOXXX(SpinorI(1,1,isigF),SpinorOGOA(1,2,iflav,isigF,isigA,isigG),GBF(1,2,iflav),
     $              ZMASS,ZWIDTH,Curr12GOt2Z(1,iflav,isigF,isigA,isigG,2,2))

                  CALL JIOXXX(SpinorIA(1,2,iflav,isigF,isigA),SpinorOGO(1,1,isigF,isigG),GBF(1,2,iflav),
     $              ZMASS,ZWIDTH,Curr21GOt1Z(1,iflav,isigF,isigA,isigG,2,1))
                  CALL JIOXXX(SpinorIAGO(1,2,iflav,isigF,isigA,isigG),SpinorO(1,1,isigF),GBF(1,2,iflav)
     $              ,ZMASS,ZWIDTH,Curr21GOt2Z(1,iflav,isigF,isigA,isigG,2,1))
                  CALL JIOXXX(SpinorI(1,2,isigF),SpinorOAGO(1,1,iflav,isigF,isigA,isigG),GBF(1,2,iflav),
     $              ZMASS,ZWIDTH,Curr21GOt2Z(1,iflav,isigF,isigA,isigG,1,1))

                  CALL JIOXXX(SpinorIGO(1,2,isigF,isigG),SpinorOA(1,1,iflav,isigF,isigA),GBF(1,2,iflav)
     $              ,ZMASS,ZWIDTH,Curr21GOt1Z(1,iflav,isigF,isigA,isigG,1,2))
                  CALL JIOXXX(SpinorIGOA(1,2,iflav,isigF,isigA,isigG),SpinorO(1,1,isigF),GBF(1,2,iflav),
     $              ZMASS,ZWIDTH,Curr21GOt2Z(1,iflav,isigF,isigA,isigG,2,2))
                  CALL JIOXXX(SpinorI(1,2,isigF),SpinorOGOA(1,1,iflav,isigF,isigA,isigG),GBF(1,2,iflav),
     $              ZMASS,ZWIDTH,Curr21GOt2Z(1,iflav,isigF,isigA,isigG,1,2))

                  CALL JIOXXX(SpinorIA(1,2,iflav,isigF,isigA),SpinorOGI(1,1,isigF,isigG),GBF(1,2,iflav),
     $              ZMASS,ZWIDTH,Curr21GIt1Z(1,iflav,isigF,isigA,isigG,2,1))
                  CALL JIOXXX(SpinorIAGI(1,2,iflav,isigF,isigA,isigG),SpinorO(1,5,isigF),GBF(1,2,iflav),
     $              ZMASS,ZWIDTH,Curr21GIt2Z(1,iflav,isigF,isigA,isigG,2,1))
                  CALL JIOXXX(SpinorI(1,2,isigF),SpinorOAGI(1,1,iflav,isigF,isigA,isigG),GBF(1,2,iflav),
     $              ZMASS,ZWIDTH,Curr21GIt2Z(1,iflav,isigF,isigA,isigG,1,1))

                  CALL JIOXXX(SpinorIGI(1,2,isigF,isigG),SpinorOA(1,5,iflav,isigF,isigA),GBF(1,2,iflav),
     $              ZMASS,ZWIDTH,Curr21GIt1Z(1,iflav,isigF,isigA,isigG,1,2)) !
                  CALL JIOXXX(SpinorIGIA(1,2,iflav,isigF,isigA,isigG),SpinorO(1,5,isigF),GBF(1,2,iflav),
     $              ZMASS,ZWIDTH,Curr21GIt2Z(1,iflav,isigF,isigA,isigG,2,2))
                  CALL JIOXXX(SpinorI(1,2,isigF),SpinorOGIA(1,1,iflav,isigF,isigA,isigG),GBF(1,2,iflav),
     $              ZMASS,ZWIDTH,Curr21GIt2Z(1,iflav,isigF,isigA,isigG,1,2))

c     lower current
                  CALL JIOXXX(SpinorIA(1,3,iflav,isigF,isigA),SpinorOGO(1,4,isigF,isigG),GBF(1,2,iflav),
     $              ZMASS,ZWIDTH,Curr34GOt1Z(1,iflav,isigF,isigA,isigG,3,4))
                  CALL JIOXXX(SpinorIAGO(1,3,iflav,isigF,isigA,isigG),SpinorO(1,4,isigF),GBF(1,2,iflav),
     $              ZMASS,ZWIDTH,Curr34GOt2Z(1,iflav,isigF,isigA,isigG,3,1))
                  CALL JIOXXX(SpinorI(1,3,isigF),SpinorOAGO(1,4,iflav,isigF,isigA,isigG),GBF(1,2,iflav),
     $              ZMASS,ZWIDTH,Curr34GOt2Z(1,iflav,isigF,isigA,isigG,4,1))

                  CALL JIOXXX(SpinorIGO(1,3,isigF,isigG),SpinorOA(1,4,iflav,isigF,isigA),GBF(1,2,iflav),
     $              ZMASS,ZWIDTH,Curr34GOt1Z(1,iflav,isigF,isigA,isigG,4,3))
                  CALL JIOXXX(SpinorIGOA(1,3,iflav,isigF,isigA,isigG),SpinorO(1,4,isigF),GBF(1,2,iflav),
     $              ZMASS,ZWIDTH,Curr34GOt2Z(1,iflav,isigF,isigA,isigG,3,2))
                  CALL JIOXXX(SpinorI(1,3,isigF),SpinorOGOA(1,4,iflav,isigF,isigA,isigG),GBF(1,2,iflav),
     $              ZMASS,ZWIDTH,Curr34GOt2Z(1,iflav,isigF,isigA,isigG,4,2))


                  CALL JIOXXX(SpinorIA(1,4,iflav,isigF,isigA),SpinorOGO(1,3,isigF,isigG),GBF(1,2,iflav),
     $              ZMASS,ZWIDTH,Curr43GOt1Z(1,iflav,isigF,isigA,isigG,4,3))
                  CALL JIOXXX(SpinorIAGO(1,4,iflav,isigF,isigA,isigG),SpinorO(1,3,isigF),GBF(1,2,iflav),
     $              ZMASS,ZWIDTH,Curr43GOt2Z(1,iflav,isigF,isigA,isigG,4,1))
                  CALL JIOXXX(SpinorI(1,4,isigF),SpinorOAGO(1,3,iflav,isigF,isigA,isigG),GBF(1,2,iflav),
     $              ZMASS,ZWIDTH,Curr43GOt2Z(1,iflav,isigF,isigA,isigG,3,1))

                  CALL JIOXXX(SpinorIGO(1,4,isigF,isigG),SpinorOA(1,3,iflav,isigF,isigA),GBF(1,2,iflav),
     $              ZMASS,ZWIDTH,Curr43GOt1Z(1,iflav,isigF,isigA,isigG,3,4))
                  CALL JIOXXX(SpinorIGOA(1,4,iflav,isigF,isigA,isigG),SpinorO(1,3,isigF),GBF(1,2,iflav),
     $              ZMASS,ZWIDTH,Curr43GOt2Z(1,iflav,isigF,isigA,isigG,4,2))
                  CALL JIOXXX(SpinorI(1,4,isigF),SpinorOGOA(1,3,iflav,isigF,isigA,isigG),GBF(1,2,iflav),
     $              ZMASS,ZWIDTH,Curr43GOt2Z(1,iflav,isigF,isigA,isigG,3,2))


                  CALL JIOXXX(SpinorIA(1,4,iflav,isigF,isigA),SpinorOGI(1,3,isigF,isigG),GBF(1,2,iflav),
     $              ZMASS,ZWIDTH,Curr43GIt1Z(1,iflav,isigF,isigA,isigG,4,3))
                  CALL JIOXXX(SpinorIAGI(1,4,iflav,isigF,isigA,isigG),SpinorO(1,5,isigF),GBF(1,2,iflav),
     $              ZMASS,ZWIDTH,Curr43GIt2Z(1,iflav,isigF,isigA,isigG,4,1))
                  CALL JIOXXX(SpinorI(1,4,isigF),SpinorOAGI(1,3,iflav,isigF,isigA,isigG),GBF(1,2,iflav),
     $              ZMASS,ZWIDTH,Curr43GIt2Z(1,iflav,isigF,isigA,isigG,3,1))

                  CALL JIOXXX(SpinorIGI(1,4,isigF,isigG),SpinorOA(1,5,iflav,isigF,isigA),GBF(1,2,iflav),
     $              ZMASS,ZWIDTH,Curr43GIt1Z(1,iflav,isigF,isigA,isigG,3,4))
                  CALL JIOXXX(SpinorIGIA(1,4,iflav,isigF,isigA,isigG),SpinorO(1,5,isigF),GBF(1,2,iflav),
     $              ZMASS,ZWIDTH,Curr43GIt2Z(1,iflav,isigF,isigA,isigG,4,2))
                  CALL JIOXXX(SpinorI(1,4,isigF),SpinorOGIA(1,3,iflav,isigF,isigA,isigG),GBF(1,2,iflav),
     $              ZMASS,ZWIDTH,Curr43GIt2Z(1,iflav,isigF,isigA,isigG,3,2))
               enddo


c     next: W Boson
c     upper current
               CALL JIOXXX(SpinorIA(1,1,  iflav,-1,isigA),SpinorOGO(1,2,-1,isigG),GBF(1,3,iflav),
     $              WMASS,WWIDTH,Curr12GOt1W(1,iflav,-1,isigA,isigG,1,2))
               CALL JIOXXX(SpinorIAGO(1,1,  iflav,-1,isigA,isigG),SpinorO(1,2,-1),GBF(1,3,iflav),
     $              WMASS,WWIDTH,Curr12GOt2W(1,iflav,-1,isigA,isigG,1,1))
               CALL JIOXXX(SpinorI(1,1,-1),SpinorOAGO(1,2,7-iflav,-1,isigA,isigG),GBF(1,3,iflav),
     $              WMASS,WWIDTH,Curr12GOt2W(1,iflav,-1,isigA,isigG,2,1))

               CALL JIOXXX(SpinorIGO(1,1,-1,isigG),SpinorOA(1,2,7-iflav,-1,isigA),GBF(1,3,iflav),
     $              WMASS,WWIDTH,Curr12GOt1W(1,iflav,-1,isigA,isigG,2,1))
               CALL JIOXXX(SpinorIGOA(1,1,  iflav,-1,isigA,isigG),SpinorO(1,2,-1),GBF(1,3,iflav),
     $              WMASS,WWIDTH,Curr12GOt2W(1,iflav,-1,isigA,isigG,1,2))
               CALL JIOXXX(SpinorI(1,1,-1),SpinorOGOA(1,2,7-iflav,-1,isigA,isigG),GBF(1,3,iflav),
     $              WMASS,WWIDTH,Curr12GOt2W(1,iflav,-1,isigA,isigG,2,2))

               CALL JIOXXX(SpinorIA(1,2,7-iflav, 1,isigA),SpinorOGO(1,1, 1,isigG),GBF(1,3,iflav),
     $              WMASS,WWIDTH,Curr21GOt1W(1,iflav, 1,isigA,isigG,2,1))
               CALL JIOXXX(SpinorIAGO(1,2,7-iflav, 1,isigA,isigG),SpinorO(1,1, 1),GBF(1,3,iflav),
     $              WMASS,WWIDTH,Curr21GOt2W(1,iflav, 1,isigA,isigG,2,1))
               CALL JIOXXX(SpinorI(1,2, 1),SpinorOAGO(1,1,  iflav, 1,isigA,isigG),GBF(1,3,iflav),
     $              WMASS,WWIDTH,Curr21GOt2W(1,iflav, 1,isigA,isigG,1,1))

               CALL JIOXXX(SpinorIGO(1,2, 1,isigG),SpinorOA(1,1,  iflav, 1,isigA),GBF(1,3,iflav),
     $              WMASS,WWIDTH,Curr21GOt1W(1,iflav, 1,isigA,isigG,1,2))
               CALL JIOXXX(SpinorIGOA(1,2,7-iflav, 1,isigA,isigG),SpinorO(1,1, 1),GBF(1,3,iflav),
     $              WMASS,WWIDTH,Curr21GOt2W(1,iflav, 1,isigA,isigG,2,2))
               CALL JIOXXX(SpinorI(1,2, 1),SpinorOGOA(1,1,  iflav, 1,isigA,isigG),GBF(1,3,iflav),
     $              WMASS,WWIDTH,Curr21GOt2W(1,iflav, 1,isigA,isigG,1,2))


               CALL JIOXXX(SpinorIA(1,2,7-iflav, 1,isigA),SpinorOGI(1,1, 1,isigG),GBF(1,3,iflav),
     $              WMASS,WWIDTH,Curr21GIt1W(1,iflav, 1,isigA,isigG,2,1))
               CALL JIOXXX(SpinorIAGI(1,2,7-iflav, 1,isigA,isigG),SpinorO(1,5, 1),GBF(1,3,iflav),
     $              WMASS,WWIDTH,Curr21GIt2W(1,iflav, 1,isigA,isigG,2,1))
               CALL JIOXXX(SpinorI(1,2, 1),SpinorOAGI(1,1,  iflav, 1,isigA,isigG),GBF(1,3,iflav),
     $              WMASS,WWIDTH,Curr21GIt2W(1,iflav, 1,isigA,isigG,1,1))

               CALL JIOXXX(SpinorIGI(1,2, 1,isigG),SpinorOA(1,5,  iflav, 1,isigA),GBF(1,3,iflav),
     $              WMASS,WWIDTH,Curr21GIt1W(1,iflav, 1,isigA,isigG,1,2))
               CALL JIOXXX(SpinorIGIA(1,2,7-iflav, 1,isigA,isigG),SpinorO(1,5, 1),GBF(1,3,iflav),
     $              WMASS,WWIDTH,Curr21GIt2W(1,iflav, 1,isigA,isigG,2,2))
               CALL JIOXXX(SpinorI(1,2, 1),SpinorOGIA(1,1,  iflav, 1,isigA,isigG),GBF(1,3,iflav),
     $              WMASS,WWIDTH,Curr21GIt2W(1,iflav, 1,isigA,isigG,1,2))

c     lower current
               CALL JIOXXX(SpinorIA(1,3,  iflav,-1,isigA),SpinorOGO(1,4,-1,isigG),GBF(1,3,iflav),
     $              WMASS,WWIDTH,Curr34GOt1W(1,iflav,-1,isigA,isigG,3,4))
               CALL JIOXXX(SpinorIAGO(1,3,  iflav,-1,isigA,isigG),SpinorO(1,4,-1),GBF(1,3,iflav),
     $              WMASS,WWIDTH,Curr34GOt2W(1,iflav,-1,isigA,isigG,3,1))
               CALL JIOXXX(SpinorI(1,3,-1),SpinorOAGO(1,4,7-iflav,-1,isigA,isigG),GBF(1,3,iflav),
     $              WMASS,WWIDTH,Curr34GOt2W(1,iflav,-1,isigA,isigG,4,1))

               CALL JIOXXX(SpinorIGO(1,3,-1,isigG),SpinorOA(1,4,7-iflav,-1,isigA),GBF(1,3,iflav),
     $              WMASS,WWIDTH,Curr34GOt1W(1,iflav,-1,isigA,isigG,4,3))
               CALL JIOXXX(SpinorIGOA(1,3,  iflav,-1,isigA,isigG),SpinorO(1,4,-1),GBF(1,3,iflav),
     $              WMASS,WWIDTH,Curr34GOt2W(1,iflav,-1,isigA,isigG,3,2))
               CALL JIOXXX(SpinorI(1,3,-1),SpinorOGOA(1,4,7-iflav,-1,isigA,isigG),GBF(1,3,iflav),
     $              WMASS,WWIDTH,Curr34GOt2W(1,iflav,-1,isigA,isigG,4,2))

               CALL JIOXXX(SpinorIA(1,4,7-iflav, 1,isigA),SpinorOGO(1,3, 1,isigG),GBF(1,3,iflav),
     $              WMASS,WWIDTH,Curr43GOt1W(1,iflav, 1,isigA,isigG,4,3))
               CALL JIOXXX(SpinorIAGO(1,4,7-iflav, 1,isigA,isigG),SpinorO(1,3, 1),GBF(1,3,iflav),
     $              WMASS,WWIDTH,Curr43GOt2W(1,iflav, 1,isigA,isigG,4,1))
               CALL JIOXXX(SpinorI(1,4, 1),SpinorOAGO(1,3,  iflav, 1,isigA,isigG),GBF(1,3,iflav),
     $              WMASS,WWIDTH,Curr43GOt2W(1,iflav, 1,isigA,isigG,3,1))

               CALL JIOXXX(SpinorIGO(1,4, 1,isigG),SpinorOA(1,3,  iflav, 1,isigA),GBF(1,3,iflav),
     $              WMASS,WWIDTH,Curr43GOt1W(1,iflav, 1,isigA,isigG,3,4))
               CALL JIOXXX(SpinorIGOA(1,4,7-iflav, 1,isigA,isigG),SpinorO(1,3, 1),GBF(1,3,iflav),
     $              WMASS,WWIDTH,Curr43GOt2W(1,iflav, 1,isigA,isigG,4,2))
               CALL JIOXXX(SpinorI(1,4, 1),SpinorOGOA(1,3,  iflav, 1,isigA,isigG),GBF(1,3,iflav),
     $              WMASS,WWIDTH,Curr43GOt2W(1,iflav, 1,isigA,isigG,3,2))


               CALL JIOXXX(SpinorIA(1,4,7-iflav, 1,isigA),SpinorOGI(1,3, 1,isigG),GBF(1,3,iflav),
     $              WMASS,WWIDTH,Curr43GIt1W(1,iflav, 1,isigA,isigG,4,3))
               CALL JIOXXX(SpinorIAGI(1,4,7-iflav, 1,isigA,isigG),SpinorO(1,5, 1),GBF(1,3,iflav),
     $              WMASS,WWIDTH,Curr43GIt2W(1,iflav, 1,isigA,isigG,4,1))
               CALL JIOXXX(SpinorI(1,4, 1),SpinorOAGI(1,3,  iflav, 1,isigA,isigG),GBF(1,3,iflav),
     $              WMASS,WWIDTH,Curr43GIt2W(1,iflav, 1,isigA,isigG,3,1))

               CALL JIOXXX(SpinorIGI(1,4, 1,isigG),SpinorOA(1,5,  iflav, 1,isigA),GBF(1,3,iflav),
     $              WMASS,WWIDTH,Curr43GIt1W(1,iflav, 1,isigA,isigG,3,4))
               CALL JIOXXX(SpinorIGIA(1,4,7-iflav, 1,isigA,isigG),SpinorO(1,5, 1),GBF(1,3,iflav),
     $              WMASS,WWIDTH,Curr43GIt2W(1,iflav, 1,isigA,isigG,4,2))
               CALL JIOXXX(SpinorI(1,4, 1),SpinorOGIA(1,3,  iflav, 1,isigA,isigG),GBF(1,3,iflav),
     $              WMASS,WWIDTH,Curr43GIt2W(1,iflav, 1,isigA,isigG,3,2))  
               
c     currents with three boson vertex attatched            !needs: Curr..(1,in-flavor,isig,iVecBos=3)
c     Curr..tbv(1,in-flavor,isigF,isigB)
               
               CALL JVVXXX(Photon(1,isigA),Curr12(1,iflav,-1,3),( 7-2*iflav)*GWWA,WMASS,WWIDTH,Curr12tbv(1,iflav,-1,isigA))
               CALL JVVXXX(Photon(1,isigA),Curr21(1,iflav, 1,3),(-7+2*iflav)*GWWA,WMASS,WWIDTH,Curr21tbv(1,iflav, 1,isigA))
               CALL JVVXXX(Curr34(1,iflav,-1,3),Photon(1,isigA),(-7+2*iflav)*GWWA,WMASS,WWIDTH,Curr34tbv(1,iflav,-1,isigA))
               CALL JVVXXX(Curr43(1,iflav, 1,3),Photon(1,isigA),( 7-2*iflav)*GWWA,WMASS,WWIDTH,Curr43tbv(1,iflav, 1,isigA))
               
               CALL JVVXXX(Photon(1,isigA),Curr12GO(1,iflav,-1,1,isigG,3),( 7-2*iflav)*GWWA,
     $              WMASS,WWIDTH,Curr12GOtbv(1,iflav,-1,isigA,1,isigG))
               CALL JVVXXX(Photon(1,isigA),Curr12GO(1,iflav,-1,2,isigG,3),( 7-2*iflav)*GWWA,
     $              WMASS,WWIDTH,Curr12GOtbv(1,iflav,-1,isigA,2,isigG))
               CALL JVVXXX(Photon(1,isigA),Curr21GO(1,iflav, 1,1,isigG,3),(-7+2*iflav)*GWWA,
     $              WMASS,WWIDTH,Curr21GOtbv(1,iflav, 1,isigA,1,isigG))
               CALL JVVXXX(Photon(1,isigA),Curr21GO(1,iflav, 1,2,isigG,3),(-7+2*iflav)*GWWA,
     $              WMASS,WWIDTH,Curr21GOtbv(1,iflav, 1,isigA,2,isigG))
               CALL JVVXXX(Curr34GO(1,iflav,-1,3,isigG,3),Photon(1,isigA),(-7+2*iflav)*GWWA,
     $              WMASS,WWIDTH,Curr34GOtbv(1,iflav,-1,isigA,3,isigG))
               CALL JVVXXX(Curr34GO(1,iflav,-1,4,isigG,3),Photon(1,isigA),(-7+2*iflav)*GWWA,
     $              WMASS,WWIDTH,Curr34GOtbv(1,iflav,-1,isigA,4,isigG))
               CALL JVVXXX(Curr43GO(1,iflav, 1,3,isigG,3),Photon(1,isigA),( 7-2*iflav)*GWWA,
     $              WMASS,WWIDTH,Curr43GOtbv(1,iflav, 1,isigA,3,isigG))
               CALL JVVXXX(Curr43GO(1,iflav, 1,4,isigG,3),Photon(1,isigA),( 7-2*iflav)*GWWA,
     $              WMASS,WWIDTH,Curr43GOtbv(1,iflav, 1,isigA,4,isigG))         
               
               CALL JVVXXX(Photon(1,isigA),Curr21GI(1,iflav, 1,1,isigG,3),(-7+2*iflav)*GWWA,
     $              WMASS,WWIDTH,Curr21GItbv(1,iflav, 1,isigA,1,isigG))
               CALL JVVXXX(Photon(1,isigA),Curr21GI(1,iflav, 1,2,isigG,3),(-7+2*iflav)*GWWA,
     $              WMASS,WWIDTH,Curr21GItbv(1,iflav, 1,isigA,2,isigG))
               CALL JVVXXX(Curr43GI(1,iflav, 1,3,isigG,3),Photon(1,isigA),( 7-2*iflav)*GWWA,
     $              WMASS,WWIDTH,Curr43GItbv(1,iflav, 1,isigA,3,isigG))
               CALL JVVXXX(Curr43GI(1,iflav, 1,4,isigG,3),Photon(1,isigA),( 7-2*iflav)*GWWA,
     $              WMASS,WWIDTH,Curr43GItbv(1,iflav, 1,isigA,4,isigG))
            enddo
         enddo
      enddo

c
c multiply propagator by 16 pi*M Gamma to compensate for the phase 
c space integral in LIPSN
c
c   Int dq^2/2pi * dOmega/(32 pi^2) * 1/(q^2-M*2)^2+(MGamma)^2) 
c
c   = 1/(16 pi MGamma)
c
c 
      if (n_v.eq.1+1) then ! only H+A
        fac = 1d0
      else
        fac = stpi*xmg(6) / ( (ph(4)-xm2(6))**2 + xmg(6)**2 )
      endif
      
      fac=12d0*fac

      return
      end
