c*************************************************************************  
c     This file contains subroutines which call the different
c     functions and subroutines to calculate matrix elements squared.
c
c*************************************************************************  
c   LIST OF ALL FUNCTIONS AND SUBROUTINES IN THIS FILE:
c
c     DOUBLE PRECISION FUNCTION Amplitude(rand, p, x, v, nlo, ps_number)
c     SUBROUTINE RE_Amplitude(rand, p, x, v, lokt, xuz, m2s, ps_number)
c     SUBROUTINE QED_Amplitude(rand, p, x, v, lokt, xuz, m2s, ps_number)
c
c*************************************************************************  

c*************************************************************************
      DOUBLE PRECISION FUNCTION Amplitude(rand, p, x, v, nlo, ps_number)
c*************************************************************************
c     This is the manager of all matrixelement routines and the (only) 
c     interface of those to the main program.  Depending on the process ID
c     the right matrix element is called with the correct calling 
c     convention, which might well be individually different.
c     Once many processes are included in the program, this subroutine
c     will be quite complex, but it should actually not be much more than a
c     big select-case statement.
c
c*************************************************************************
c     INPUT
c     rand      : Array with random numbers driven by monaco. 
c     p         : 4-momenta of the partons involved in the process.
c                 The first 2 entries are the incoming partons.
c     x         : The feynman x of the incoming partons 
c     v         : The 4-momenta of the other particles
c     nlo       : switch for NLO processes
c     ps_number : Number specifying the desired matrix element squared
c     
c     OUTPUT
c     Amplitude : The amplitude square as returned by the matrix element code.
c*************************************************************************
          use globalvars, only: ldoblha

      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/process.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"

      double precision x(nx) 
      double precision rand(max_PS_dim)
      double precision p(0:3,max_p,max_kin)
      double precision v(0:3,max_v,max_kin)

      double precision rn, rn2(2), rn4(4), xuz(2,2:max_kin), m2s(0:max_kin)
      double precision m2s_dec
      logical lokt(max_kin)
      integer NLO,i, ps_number

c     to convert "NLO" to box- and pentagon call conventions
c     last entry is for non-interfering closed fermion loops 
c     (GF contributions)
      integer VBFNLOswitch(0:11)
      data VBFNLOswitch / 0, 4,-4,-5,-6,-8,-999,-999,-999,-999,-999,-8 /

c     masses and widths
      double precision XM2s(6),XMGs(6)
      COMMON /BKOPOUshort/ XM2s,XMGs


c     external matrixelement functions
      double precision m2s_ubtd
      double precision m2s_wbfvv
      double precision m2s_wbfw
      double precision m2s_ggfh
      double precision m2s_QCDzz2j
      double precision m2s_QCDwz2j
      double precision m2s_QCDww2j
      double precision RandomNumber
      double precision m2s_qqVV
      double precision m2s_qqZZ
      double precision m2s_qqZA
      double precision m2s_qqAA
      double precision m2s_qqVVj2
      double precision m2s_qqZZj2
      double precision m2s_qqWWZ
      double precision m2s_qqZZZ
      double precision m2s_qqZZW
      double precision m2s_qqWWW
      double precision m2s_qqWWA
      double precision m2s_qqZZA
      double precision m2s_qqWZA
      double precision m2s_qqWAA
      double precision m2s_qqZAA
      double precision m2s_qqAAA
      double precision m2s_qqWAAj2
      double precision m2s_qqwgaglu_hel
      double precision m2s_qqwzglu_hel 
      double precision m2s_qqwpgaglu_hel
      double precision m2s_qqwpzglu_hel
      double precision m2s_qqWZ
      double precision m2s_qqWA
      double precision m2s_wbfwpp
      double precision m2s_qqWH, m2s_qqWHj2 


      integer mu
      external m2s_ubtd
      external m2s_wbfvv
      external m2s_wbfw
      external m2s_ggfh
      external m2s_QCDwz2j
      external m2s_QCDww2j
      external m2s_QCDzz2j
      external RandomNumber
      external m2s_qqVV
      external m2s_qqZZ
      external m2s_qqZA
      external m2s_qqAA
      external m2s_qqVVj2
      external m2s_qqZZj2
      external m2s_qqWWZ
      external m2s_qqZZZ
      external m2s_qqZZW
      external m2s_qqWWW
      external m2s_qqWWA
      external m2s_qqZZA
      external m2s_qqWZA
      external m2s_qqWAA
      external m2s_qqZAA
      external m2s_qqAAA
      external m2s_qqWAAj2
      external m2s_qqwgaglu_hel
      external m2s_qqwzglu_hel 
      external m2s_qqwpgaglu_hel
      external m2s_qqwpzglu_hel
      external m2s_qqWZ
      external m2s_qqWA
      external m2s_wbfwpp
      external m2s_qqWH, m2s_qqWHj2

c     auxiliary 
      double precision amp, Q(0:3)

      double precision mass2
      external mass2

c initialize
      do mu=1,2
        do i=2,max_kin
          xuz(mu,i) = 0d0
        enddo
      enddo
      if (ldoblha) VBFNLOswitch(1) = 1

      SELECT CASE(process)

#ifdef WITH_VBF
**********************
*
*     VBF : Hjj
*
**********************
c     pp -> Hjj
      CASE(Hjj)
         rn = RandomNumber()
         Amplitude = m2s_wbfvv(6,nlo,x,p(0,1,1),v(0,1,1),rn, ps_number)
c     pp -> Hjj -> gamma gamma jj
      CASE(Hjj_AA)
         rn = RandomNumber()
         Amplitude = m2s_wbfvv(6,nlo,x,p(0,1,1),v(0,1,1),rn, ps_number)
         AMplitude = Amplitude*BHGAM
c     pp -> spin-2 (or Higgs) jj -> aajj
      CASE(AAjj)
         rn = RandomNumber()
         Amplitude = m2s_wbfvv(11,VBFNLOswitch(nlo),x,p(0,1,1),v(0,1,1),rn, ps_number)
c     pp -> spin-2 jj -> WWjj
      CASE(Sp2jj_WW)
         rn = RandomNumber()
         Amplitude = m2s_wbfvv(734,VBFNLOswitch(nlo),x,p(0,1,1),v(0,1,1),rn, ps_number)
c     pp -> spin-2 jj -> ZZjj
      CASE(Sp2jj_ZZ_ll)
         rn = RandomNumber()
         Amplitude = m2s_wbfvv(722,VBFNLOswitch(nlo),x,p(0,1,1),v(0,1,1),rn, ps_number)
      CASE(Sp2jj_ZZ_lnu)
         rn = RandomNumber()
         Amplitude = m2s_wbfvv(721,VBFNLOswitch(nlo),x,p(0,1,1),v(0,1,1),rn, ps_number)
c     pp -> Hjj -> mu+mu- jj
      CASE(Hjj_mu)
         rn = RandomNumber()
         Amplitude = m2s_wbfvv(6,nlo,x,p(0,1,1),v(0,1,1),rn, ps_number)
         Amplitude = Amplitude*BHMU
c     pp -> Hjj -> tau+tau- jj
      CASE(Hjj_tau)
         rn = RandomNumber()
         Amplitude = m2s_wbfvv(6,nlo,x,p(0,1,1),v(0,1,1),rn, ps_number)
         Amplitude = Amplitude*BHTAU
c     pp -> Hjj -> b bbar jj
      CASE(Hjj_bbar)
         rn = RandomNumber()
         Amplitude = m2s_wbfvv(6,nlo,x,p(0,1,1),v(0,1,1),rn, ps_number)
         Amplitude = AMplitude*BHBB
c     pp -> Hjj -> W W jj
      CASE(Hjj_WW)
         rn = RandomNumber()
         Amplitude = m2s_wbfvv(6,nlo,x,p(0,1,1),v(0,1,1),rn, ps_number)
         call m2s_VVsum(v(0,1,1),3,4,1,2,m2s_dec)
         m2s_dec = m2s_dec/(16.0d0 * pi * xmgs(6))
         Amplitude = Amplitude*m2s_dec
      CASE(Hjj_WPhadWM, Hjj_WPWMhad)
         rn = RandomNumber()
         Amplitude = m2s_wbfvv(6,nlo,x,p(0,1,1),v(0,1,1),rn, ps_number)
         if (process.eq.Hjj_WPhadWM) then
            call m2s_VVsum_had(v(0,1,1),3,4,3,2,m2s_dec)
         else
            call m2s_VVsum_had(v(0,1,1),3,4,1,4,m2s_dec)
         endif
         m2s_dec = m2s_dec/(16.0d0 * pi * xmgs(6))
         Amplitude = Amplitude*m2s_dec
c     pp -> Hjj -> Z Z jj 
      CASE(Hjj_ZZ_ll)
         rn = RandomNumber()
         Amplitude = m2s_wbfvv(6,nlo,x,p(0,1,1),v(0,1,1),rn, ps_number)
         call m2s_VVsum(v(0,1,1),2,2,2,2,m2s_dec)
         m2s_dec = m2s_dec/(16.0d0 * pi * xmgs(6))
         Amplitude = Amplitude*m2s_dec
      CASE(Hjj_ZZhad)
         rn = RandomNumber()
         Amplitude = m2s_wbfvv(6,nlo,x,p(0,1,1),v(0,1,1),rn, ps_number)
         call m2s_VVsum_had(v(0,1,1),2,2,5,2,m2s_dec)
         m2s_dec = m2s_dec/(16.0d0 * pi * xmgs(6))
         Amplitude = Amplitude*m2s_dec
      CASE(Hjj_ZZ_lnu)
         rn = RandomNumber()
         Amplitude = m2s_wbfvv(6,nlo,x,p(0,1,1),v(0,1,1),rn, ps_number)
         call m2s_VVsum(v(0,1,1),2,2,2,1,m2s_dec)
         m2s_dec = m2s_dec/(16.0d0 * pi * xmgs(6))
         Amplitude = Amplitude*m2s_dec  

**********************
*
*     VBF : HHjj
*
**********************
c     pp -> HHjj
      CASE(HHjj)
         rn = RandomNumber()
         Amplitude = m2s_wbfvv(66,nlo,x,p(0,1,1),v(0,1,1),rn, ps_number)
c     pp -> HHjj -> b bbar tau+tau- jj
      CASE(HHjj_bbtau)
         rn = RandomNumber()
         Amplitude = m2s_wbfvv(66,nlo,x,p(0,1,1),v(0,1,1),rn)
         Amplitude = Amplitude*BHBB*BHTAU
c     pp -> HHjj -> b bbar gamma gamma jj
      CASE(HHjj_bbAA)
         rn = RandomNumber()
         Amplitude = m2s_wbfvv(66,nlo,x,p(0,1,1),v(0,1,1),rn)
         Amplitude = Amplitude*BHBB*BHGAM
c     pp -> HHjjj
      CASE(HHjjjLO)
         rn = RandomNumber()
         call VBFptildemake(3,p(0,1,1),xuz,v(0,1,1))
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         call m2s_wbfzhg( 66, 0, lokt, 
     1        x, p(0,1,1), v(0,1,1), rn, xuz, m2s)
         Amplitude=m2s(0)
c     pp -> HHjjj -> b bbar tau+tau- jj
      CASE(HHjjjLO_bbtau)
         rn = RandomNumber()
         call VBFptildemake(3,p(0,1,1),xuz,v(0,1,1))
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         call m2s_wbfzhg( 66, 0, lokt, 
     1        x, p(0,1,1), v(0,1,1), rn, xuz, m2s)
         Amplitude=m2s(0)
         Amplitude=Amplitude*BHBB*BHTAU
c     pp -> HHjjj -> b bbar gamma gamma jj
      CASE(HHjjjLO_bbAA)
         rn = RandomNumber()
         call VBFptildemake(3,p(0,1,1),xuz,v(0,1,1))
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         call m2s_wbfzhg( 66, 0, lokt, 
     1        x, p(0,1,1), v(0,1,1), rn, xuz, m2s)
         Amplitude=m2s(0)
         Amplitude=Amplitude*BHBB*BHGAM

**********************
*
*     VBF : Hjjj
*
**********************
c     pp -> Hjjj
      CASE(HjjjLO)
         rn = RandomNumber()
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         call m2s_wbfzhg( 6, 0, lokt,
     1        x, p(0,1,1), v(0,1,1), rn, xuz, m2s)
         Amplitude=m2s(0)
c     pp -> Hjjj -> gamma gamma jjj
      CASE(HjjjLO_AA)
         rn = RandomNumber()
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         call m2s_wbfzhg( 6, 0, lokt,
     1        x, p(0,1,1), v(0,1,1), rn, xuz, m2s)
         Amplitude=m2s(0)
         Amplitude = Amplitude*BHGAM
c     pp -> spin-2 jjj -> gamma gamma jjj
      CASE(AAjjjLO)
         rn = RandomNumber()
         call VBFptildemake(3,p(0,1,1),xuz,v(0,1,1))
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         call m2s_wbfzhg( 11, VBFNLOswitch(nlo), lokt, 
     1                    x, p(0,1,1), v(0,1,1), rn, xuz, m2s)
         Amplitude=m2s(0)
c     pp -> spin-2 jjj -> WW jjj
      CASE(Sp2jjjLO_WW)
         rn = RandomNumber()
         call VBFptildemake(3,p(0,1,1),xuz,v(0,1,1))
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         call m2s_wbfzhg( 734, VBFNLOswitch(nlo), lokt, 
     1                    x, p(0,1,1), v(0,1,1), rn, xuz, m2s)
         Amplitude=m2s(0)
c     pp -> spin-2 jjj -> ZZ jjj
      CASE(Sp2jjjLO_ZZ_ll)
         rn = RandomNumber()
         call VBFptildemake(3,p(0,1,1),xuz,v(0,1,1))
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         call m2s_wbfzhg( 722, VBFNLOswitch(nlo), lokt, 
     1                    x, p(0,1,1), v(0,1,1), rn, xuz, m2s)
         Amplitude=m2s(0)

      CASE(Sp2jjjLO_ZZ_lnu)
         rn = RandomNumber()
         call VBFptildemake(3,p(0,1,1),xuz,v(0,1,1))
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         call m2s_wbfzhg( 721, VBFNLOswitch(nlo), lokt, 
     1                    x, p(0,1,1), v(0,1,1), rn, xuz, m2s)
         Amplitude=m2s(0)
c     pp -> Hjjj -> mu+mu- jjj
      CASE(HjjjLO_mu)
         rn = RandomNumber()
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         call m2s_wbfzhg( 6, 0, lokt,
     1        x, p(0,1,1), v(0,1,1), rn, xuz, m2s)
         Amplitude=m2s(0)
         Amplitude = Amplitude*BHMU
c     pp -> Hjjj -> tau+tau- jjj
      CASE(HjjjLO_tau)
         rn = RandomNumber()
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         call m2s_wbfzhg( 6, 0, lokt,
     1        x, p(0,1,1), v(0,1,1), rn, xuz, m2s)
         Amplitude=m2s(0)
         Amplitude = Amplitude*BHTAU
c     pp -> Hjjj -> b bbar jjj
      CASE(HjjjLO_bbar)
         rn = RandomNumber()
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         call m2s_wbfzhg( 6, 0, lokt,
     1        x, p(0,1,1), v(0,1,1), rn, xuz, m2s)
         Amplitude=m2s(0)
         Amplitude = AMplitude*BHBB
c     pp -> Hjjj -> W W jjj
      CASE(HjjjLO_WW)
         rn = RandomNumber()
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         call m2s_wbfzhg( 6, 0, lokt,
     1        x, p(0,1,1), v(0,1,1), rn, xuz, m2s)
         Amplitude=m2s(0)
         call m2s_VVsum(v(0,1,1),3,4,1,2,m2s_dec)
         m2s_dec = m2s_dec/(16.0d0 * pi * xmgs(6))
         Amplitude = Amplitude*m2s_dec
      CASE(HjjjLO_WPhadWM, HjjjLO_WPWMhad)
         rn = RandomNumber()
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         call m2s_wbfzhg( 6, 0, lokt,
     1        x, p(0,1,1), v(0,1,1), rn, xuz, m2s)
         Amplitude=m2s(0)
         if (process.eq.HjjjLO_WPhadWM) then
            call m2s_VVsum_had(v(0,1,1),3,4,3,2,m2s_dec)
         else
            call m2s_VVsum_had(v(0,1,1),3,4,1,4,m2s_dec)
         endif
         m2s_dec = m2s_dec/(16.0d0 * pi * xmgs(6))
         Amplitude = Amplitude*m2s_dec
c     pp -> Hjjj -> Z Z jjj
      CASE(HjjjLO_ZZ_ll)
         rn = RandomNumber()
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         call m2s_wbfzhg( 6,0, lokt,
     1        x, p(0,1,1), v(0,1,1), rn, xuz, m2s)
         Amplitude=m2s(0)
         call m2s_VVsum(v(0,1,1),2,2,2,2,m2s_dec)
         m2s_dec = m2s_dec/(16.0d0 * pi * xmgs(6))
         Amplitude = Amplitude*m2s_dec
      CASE(HjjjLO_ZZhad)
         rn = RandomNumber()
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         call m2s_wbfzhg( 6,0, lokt,
     1        x, p(0,1,1), v(0,1,1), rn, xuz, m2s)
         Amplitude=m2s(0)
         call m2s_VVsum_had(v(0,1,1),2,2,5,2,m2s_dec)
         m2s_dec = m2s_dec/(16.0d0 * pi * xmgs(6))
         Amplitude = Amplitude*m2s_dec
      CASE(HjjjLO_ZZ_lnu)
         rn = RandomNumber()
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         call m2s_wbfzhg( 6, 0, lokt,
     1        x, p(0,1,1), v(0,1,1), rn, xuz, m2s)
         Amplitude=m2s(0)
         call m2s_VVsum(v(0,1,1),2,2,2,1,m2s_dec)
         m2s_dec = m2s_dec/(16.0d0 * pi * xmgs(6))
         Amplitude = Amplitude*m2s_dec
#endif

#ifdef WITH_HJJJ
********************
*
*     VBF : Hjjj
*
********************
c     pp -> Hjjj
      CASE(Hjjj)
         rn = RandomNumber()
         call m2s_wbfhj(6, nlo, x, p(0,1,1), v(0,1,1), rn, m2s(0))
         Amplitude = m2s(0)
c         pp -> Hjjj -> gamma gamma jjj
      CASE(Hjjj_AA)
         rn = RandomNumber()
         call m2s_wbfhj(6, nlo, x, p(0,1,1), v(0,1,1), rn, m2s(0))
         Amplitude = m2s(0)
         Amplitude = Amplitude*BHGAM
c     pp -> Hjjj -> mu+mu- jjj
      CASE(Hjjj_mu)
         rn = RandomNumber()
         call m2s_wbfhj(6, nlo, x, p(0,1,1), v(0,1,1), rn, m2s(0))
         Amplitude = m2s(0)
         Amplitude = Amplitude*BHMU
c     pp -> Hjjj -> tau+tau- jjj
      CASE(Hjjj_tau)
         rn = RandomNumber()
         call m2s_wbfhj(6, nlo, x, p(0,1,1), v(0,1,1), rn, m2s(0))
         Amplitude = m2s(0)
         Amplitude = Amplitude*BHTAU
c     pp -> Hjjj -> b bbar jjj
      CASE(Hjjj_bbar)
         rn = RandomNumber()
         call m2s_wbfhj(6, nlo, x, p(0,1,1), v(0,1,1), rn, m2s(0))
         Amplitude = m2s(0)
         Amplitude = AMplitude*BHBB
c     pp -> Hjjj -> W W jjj
      CASE(Hjjj_WW)
         rn = RandomNumber()
         call m2s_wbfhj(6, nlo, x, p(0,1,1), v(0,1,1), rn, m2s(0))
         Amplitude = m2s(0)
         call m2s_VVsum(v(0,1,1),3,4,2,1,m2s_dec)
         m2s_dec = m2s_dec/(16.0d0 * pi * xmgs(6))
         Amplitude = Amplitude*m2s_dec
c     pp -> Hjjj -> Z Z jjj 
      CASE(Hjjj_ZZ_ll)
         rn = RandomNumber()
         call m2s_wbfhj(6, nlo, x, p(0,1,1), v(0,1,1), rn, m2s(0))
         Amplitude = m2s(0)
         call m2s_VVsum(v(0,1,1),2,2,2,2,m2s_dec)
         m2s_dec = m2s_dec/(16.0d0 * pi * xmgs(6))
         Amplitude = Amplitude*m2s_dec
      CASE(Hjjj_ZZ_lnu)
         rn = RandomNumber()
         call m2s_wbfhj(6, nlo, x, p(0,1,1), v(0,1,1), rn, m2s(0))
         Amplitude = m2s(0)
         call m2s_VVsum(v(0,1,1),2,2,2,1,m2s_dec)
         m2s_dec = m2s_dec/(16.0d0 * pi * xmgs(6))
         Amplitude = Amplitude*m2s_dec

*********************
*
*     VBF : Hjjjj
*
*********************
c     pp -> Hjjjj
      CASE(HjjjjLO)
         do i=1,n_kin
            lokt(i)=.false.
         enddo
         lokt(1)=.true.
         call m2s_wbfh4j(6, 0, lokt, x, p(0,1,1), v(0,1,1), xuz, m2s(0))
         Amplitude=m2s(0)
c     pp -> Hjjjj-> gamma gamma jjjj
      CASE(HjjjjLO_AA)
         do i=1,n_kin
            lokt(i)=.false.
         enddo
         lokt(1)=.true.
         call m2s_wbfh4j(6, 0, lokt, x, p(0,1,1), v(0,1,1), xuz, m2s(0))
         Amplitude=m2s(0) 
         Amplitude = Amplitude*BHGAM
c     pp -> Hjjjj -> mu+mu- jjjj
      CASE(HjjjjLO_mu)
         do i=1,n_kin
            lokt(i)=.false.
         enddo
         lokt(1)=.true.
         call m2s_wbfh4j(6, 0, lokt, x, p(0,1,1), v(0,1,1), xuz, m2s(0))
         Amplitude=m2s(0)
         Amplitude = Amplitude*BHMU
c     pp -> Hjjjj -> tau+tau- jjjj
      CASE(HjjjjLO_tau)
         do i=1,n_kin
            lokt(i)=.false.
         enddo
         lokt(1)=.true.
         call m2s_wbfh4j(6, 0, lokt, x, p(0,1,1), v(0,1,1), xuz, m2s(0))
         Amplitude=m2s(0)
         Amplitude = Amplitude*BHTAU
c     pp -> Hjjjj -> b bbar jjjj
      CASE(HjjjjLO_bbar)
         do i=1,n_kin
            lokt(i)=.false.
         enddo
         lokt(1)=.true.
         call m2s_wbfh4j(6, 0, lokt, x, p(0,1,1), v(0,1,1), xuz, m2s(0))
         Amplitude=m2s(0)
         Amplitude = Amplitude*BHBB
c     pp -> Hjjjj -> W W jjjj
      CASE(HjjjjLO_WW)
         do i=1,n_kin
            lokt(i)=.false.
         enddo
         lokt(1)=.true.
         call m2s_wbfh4j(6, 0, lokt, x, p(0,1,1), v(0,1,1), xuz, m2s(0))
         Amplitude = m2s(0)
         call m2s_VVsum(v(0,1,1),3,4,2,1,m2s_dec)
         m2s_dec = m2s_dec/(16.0d0 * pi * xmgs(6))
         Amplitude = Amplitude*m2s_dec
c     pp -> Hjjjj -> Z Z jjjj 
      CASE(HjjjjLO_ZZ_ll)
         do i=1,n_kin
            lokt(i)=.false.
         enddo
         lokt(1)=.true.
         call m2s_wbfh4j(6, 0, lokt, x, p(0,1,1), v(0,1,1), xuz, m2s(0))
         Amplitude = m2s(0)
         call m2s_VVsum(v(0,1,1),2,2,2,2,m2s_dec)
         m2s_dec = m2s_dec/(16.0d0 * pi * xmgs(6))
         Amplitude = Amplitude*m2s_dec
      CASE(HjjjjLO_ZZ_lnu)
         do i=1,n_kin
            lokt(i)=.false.
         enddo
         lokt(1)=.true.
         call m2s_wbfh4j(6, 0, lokt, x, p(0,1,1), v(0,1,1), xuz, m2s(0))
         Amplitude = m2s(0)
         call m2s_VVsum(v(0,1,1),2,2,2,1,m2s_dec)
         m2s_dec = m2s_dec/(16.0d0 * pi * xmgs(6))
         Amplitude = Amplitude*m2s_dec 
         
#endif

#ifdef WITH_VBF
***********************
*
*     VBF : Wjj, Zjj, Ajj
*
***********************
c     pp -> Zjj
      CASE(Zjj_l, Zjj_nu) 
         rn = RandomNumber()
         Amplitude = m2s_wbfvv(2,nlo,x,p(0,1,1),v(0,1,1),rn, ps_number)
c     pp -> W+jj
      CASE(WPjj) 
         rn = RandomNumber()
         Amplitude = m2s_wbfw(3,nlo,x,p(0,1,1),v(0,1,1),rn)
c     pp -> W-jj
      CASE(WMjj) 
         rn = RandomNumber()
         Amplitude = m2s_wbfw(4,nlo,x,p(0,1,1),v(0,1,1),rn)
c     pp -> Ajj
      CASE(Ajj)
         rn = RandomNumber()
         Amplitude = m2s_wbfvv(1,nlo,x,p(0,1,1),v(0,1,1),rn, ps_number)

**************************
*
*     VBF : Wjjj, Zjjj, Ajjj
*
**************************
c     pp -> Zjjj
      CASE(ZjjjLO_l, ZjjjLO_nu) 
         rn = RandomNumber()
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         call m2s_wbfzhg( 2, 0, lokt, 
     1        x, p(0,1,1), v(0,1,1), rn, xuz, m2s)
         Amplitude=m2s(0)
c     pp -> W+jjj
      CASE(WPjjjLO)
         rn = RandomNumber()
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         call m2s_wbfw3j( 3, 0, lokt, 
     1        x, p(0,1,1), v(0,1,1), rn, xuz, m2s)
         Amplitude=m2s(0)
c     pp -> W-jjj
      CASE(WMjjjLO)
         rn = RandomNumber()
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         call m2s_wbfw3j( 4, 0, lokt, 
     1        x, p(0,1,1), v(0,1,1), rn, xuz, m2s)
         Amplitude=m2s(0)
c     pp ->Ajjj
      CASE(AjjjLO)
         rn = RandomNumber()
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         call m2s_wbfzhg( 1, 0, lokt,
     1        x, p(0,1,1), v(0,1,1), rn, xuz, m2s)
         Amplitude=m2s(0)

************************************
*
*     VBF : WWjj, ZZjj, WZjj, WAjj, ZAjj
*
************************************
c     pp -> WWjj 
      CASE(WPWMjj, WPhadWMjj, WPWMhadjj)
         rn = RandomNumber()
         Amplitude = m2s_wbfvv(34,VBFNLOswitch(nlo),x,p(0,1,1),v(0,1,1),rn, ps_number)
c     pp -> W+W+jj 
      CASE(WPWPjj,WMWMjj, WPhadWPjj, WMhadWMjj)  ! note that W-W-jj is obtained from W+W+jj
         rn = RandomNumber()
         Amplitude = m2s_wbfwpp(33,VBFNLOswitch(nlo),x,p(0,1,1),v(0,1,1),rn)
c     pp -> ZZjj 
      CASE(ZZjj_ll, ZZhadjj)
         rn = RandomNumber()
         Amplitude = m2s_wbfvv(22,VBFNLOswitch(nlo),x,p(0,1,1),v(0,1,1),rn, ps_number)
      CASE(ZZjj_lnu)
         rn = RandomNumber()
         Amplitude = m2s_wbfvv(21,VBFNLOswitch(nlo),x,p(0,1,1),v(0,1,1),rn, ps_number)
c     pp -> W+Zjj
      CASE(WPZjj, WPhadZjj, WPZhadjj)
         rn = RandomNumber()
         Amplitude = m2s_wbfw(32,VBFNLOswitch(nlo),x,p(0,1,1),v(0,1,1),rn,ps_number)
c     pp -> W-Zjj
      CASE(WMZjj, WMhadZjj, WMZhadjj)
         rn = RandomNumber()
         Amplitude = m2s_wbfw(42,VBFNLOswitch(nlo),x,p(0,1,1),v(0,1,1),rn,ps_number)
c     pp -> W+Ajj
      CASE(WPAjj)
         rn = RandomNumber()
         Amplitude = m2s_wbfw(31,VBFNLOswitch(nlo),x,p(0,1,1),v(0,1,1),rn,ps_number)
c     pp -> W-Ajj
      CASE(WMAjj)
         rn = RandomNumber()
         Amplitude = m2s_wbfw(41,VBFNLOswitch(nlo),x,p(0,1,1),v(0,1,1),rn,ps_number)
c     pp -> ZAjj
      CASE(ZAjj)
         rn = RandomNumber()
         Amplitude = m2s_wbfvv(211,VBFNLOswitch(nlo),x,p(0,1,1),v(0,1,1),rn,ps_number)
      CASE(ZAjj_n)
         rn = RandomNumber()
         Amplitude = m2s_wbfvv(212,VBFNLOswitch(nlo),x,p(0,1,1),v(0,1,1),rn,ps_number)

***************************************
*
*     VBF : WWjjj, ZZjjj, WZjjj, WAjjj
*
***************************************
c     pp -> WW jjj
      CASE(WPWMjjjLO, WPhadWMjjjLO, WPWMhadjjjLO)
         rn = RandomNumber()
         call VBFptildemake(3,p(0,1,1),xuz,v(0,1,1))
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         call m2s_wbfzhg( 34, VBFNLOswitch(nlo), lokt, 
     1                    x, p(0,1,1), v(0,1,1), rn, xuz, m2s)
         Amplitude=m2s(0)
c     pp -> ZZ jjj
      CASE(ZZjjjLO_ll, ZZhadjjjLO)
         rn = RandomNumber()
         call VBFptildemake(3,p(0,1,1),xuz,v(0,1,1))
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         call m2s_wbfzhg( 22, VBFNLOswitch(nlo), lokt, 
     1                    x, p(0,1,1), v(0,1,1), rn, xuz, m2s)
         Amplitude=m2s(0)
      CASE(ZZjjjLO_lnu)
         rn = RandomNumber()
         call VBFptildemake(3,p(0,1,1),xuz,v(0,1,1))
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         call m2s_wbfzhg( 21, VBFNLOswitch(nlo), lokt, 
     1                    x, p(0,1,1), v(0,1,1), rn, xuz, m2s)
         Amplitude=m2s(0)
c     pp -> W+Z jjj
      CASE(WPZjjjLO, WPhadZjjjLO, WPZhadjjjLO)
         rn = RandomNumber()
         call VBFptildemake(3,p(0,1,1),xuz,v(0,1,1))
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         call m2s_wbfw3j( 32,0,  lokt, 
     1        x, p(0,1,1), v(0,1,1), rn, xuz, m2s)
         Amplitude=m2s(0)
c     pp -> W-Z jjj
      CASE(WMZjjjLO, WMhadZjjjLO, WMZhadjjjLO)
         rn = RandomNumber()
         call VBFptildemake(3,p(0,1,1),xuz,v(0,1,1))
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         call m2s_wbfw3j( 42,0, lokt, 
     1        x, p(0,1,1), v(0,1,1), rn, xuz, m2s)
         Amplitude=m2s(0)
c     pp -> W+A jjj
      CASE(WPAjjjLO)
         rn = RandomNumber()
         call VBFptildemake(3,p(0,1,1),xuz,v(0,1,1))
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         call m2s_wbfw3j( 31,0,  lokt, 
     1        x, p(0,1,1), v(0,1,1), rn, xuz, m2s)
         Amplitude=m2s(0)
c     pp -> W-A jjj
      CASE(WMAjjjLO)
         rn = RandomNumber()
         call VBFptildemake(3,p(0,1,1),xuz,v(0,1,1))
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         call m2s_wbfw3j( 41,0, lokt, 
     1        x, p(0,1,1), v(0,1,1), rn, xuz, m2s)
         Amplitude=m2s(0)
c     pp -> W+W+ jjj
      CASE(WPWPjjjLO, WMWMjjjLO, WPhadWPjjjLO, WMhadWMjjjLO)
         rn = RandomNumber()
         call VBFptildemake(3,p(0,1,1),xuz,v(0,1,1))
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         call m2s_wbfwppg( 33,0, lokt, 
!         call m2s_wbfwppg( 33,VBFNLOswitch(nlo), lokt, 
     1                    x, p(0,1,1), v(0,1,1), rn, xuz, m2s) 
         Amplitude=m2s(0)
c     pp -> ZA jjj
      CASE(ZAjjjLO)
         rn = RandomNumber()
         call VBFptildemake(3,p(0,1,1),xuz,v(0,1,1))
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         call m2s_wbfzhg( 211, 0, lokt, 
     1                    x, p(0,1,1), v(0,1,1), rn, xuz, m2s)
         Amplitude=m2s(0)
c     pp -> ZA jjj
      CASE(ZAjjjLO_n)
         rn = RandomNumber()
         call VBFptildemake(3,p(0,1,1),xuz,v(0,1,1))
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         call m2s_wbfzhg( 212, 0, lokt, 
     1                    x, p(0,1,1), v(0,1,1), rn, xuz, m2s)
         Amplitude=m2s(0)

***********************
*
*     VBF : HAjj
*
***********************
c     pp -> HAjj
      CASE(HAjj)
         rn = RandomNumber()
         Amplitude = m2s_wbfvv(61,nlo,x,p(0,1,1),v(0,1,1),rn, ps_number)
c     pp -> HAjj -> A gamma gamma jj
      CASE(HAjj_AA)
         rn = RandomNumber()
         Amplitude = m2s_wbfvv(61,nlo,x,p(0,1,1),v(0,1,1),rn, ps_number)
         AMplitude = Amplitude*BHGAM
c     pp -> HAjj -> A mu+mu- jj
      CASE(HAjj_mu)
         rn = RandomNumber()
         Amplitude = m2s_wbfvv(61,nlo,x,p(0,1,1),v(0,1,1),rn, ps_number)
         Amplitude = Amplitude*BHMU
c     pp -> HAjj -> A tau+tau- jj
      CASE(HAjj_tau)
         rn = RandomNumber()
         Amplitude = m2s_wbfvv(61,nlo,x,p(0,1,1),v(0,1,1),rn, ps_number)
         Amplitude = Amplitude*BHTAU
c     pp -> Hjj -> A b bbar jj
      CASE(HAjj_bbar)
         rn = RandomNumber()
         Amplitude = m2s_wbfvv(61,nlo,x,p(0,1,1),v(0,1,1),rn, ps_number)
         Amplitude = Amplitude*BHBB
c     pp -> HAjj -> A W W jj
      CASE(HAjj_WW)
         rn = RandomNumber()
         Amplitude = m2s_wbfvv(61,nlo,x,p(0,1,1),v(0,1,1),rn, ps_number)
         call m2s_VVsum(v(0,1,1),3,4,2,1,m2s_dec)
         m2s_dec = m2s_dec/(16.0d0 * pi * xmgs(6))
         Amplitude = Amplitude*m2s_dec
c     pp -> HAjj -> A Z Z jj 
      CASE(HAjj_ZZ_ll)
         rn = RandomNumber()
         Amplitude = m2s_wbfvv(61,nlo,x,p(0,1,1),v(0,1,1),rn, ps_number)
         call m2s_VVsum(v(0,1,1),2,2,2,2,m2s_dec)
         m2s_dec = m2s_dec/(16.0d0 * pi * xmgs(6))
         Amplitude = Amplitude*m2s_dec
      CASE(HAjj_ZZ_lnu)
         rn = RandomNumber()
         Amplitude = m2s_wbfvv(61,nlo,x,p(0,1,1),v(0,1,1),rn, ps_number)
         call m2s_VVsum(v(0,1,1),2,2,2,1,m2s_dec)
         m2s_dec = m2s_dec/(16.0d0 * pi * xmgs(6))
         Amplitude = Amplitude*m2s_dec  


***********************
*
*     VBF : HAjjj
*
***********************
c     pp -> HAjjj
      CASE(HAjjjLO)
         rn = RandomNumber()
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         call m2s_wbfzhg( 16, 0, lokt,
     1        x, p(0,1,1), v(0,1,1), rn, xuz, m2s)
         Amplitude=m2s(0)
c     pp -> HAjjjLO -> A gamma gamma jj
      CASE(HAjjjLO_AA)
         rn = RandomNumber()
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         call m2s_wbfzhg( 16, 0, lokt,
     1        x, p(0,1,1), v(0,1,1), rn, xuz, m2s)
         Amplitude=m2s(0)
         AMplitude = Amplitude*BHGAM
c     pp -> HAjjjLO -> A mu+mu- jj
      CASE(HAjjjLO_mu)
         rn = RandomNumber()
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         call m2s_wbfzhg( 16, 0, lokt,
     1        x, p(0,1,1), v(0,1,1), rn, xuz, m2s)
         Amplitude=m2s(0)
         Amplitude = Amplitude*BHMU
c     pp -> HAjjjLO -> A tau+tau- jj
      CASE(HAjjjLO_tau)
         rn = RandomNumber()
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         call m2s_wbfzhg( 16, 0, lokt,
     1        x, p(0,1,1), v(0,1,1), rn, xuz, m2s)
         Amplitude=m2s(0)
         Amplitude = Amplitude*BHTAU
c     pp -> Hjj -> A b bbar jj
      CASE(HAjjjLO_bbar)
         rn = RandomNumber()
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         call m2s_wbfzhg( 16, 0, lokt,
     1        x, p(0,1,1), v(0,1,1), rn, xuz, m2s)
         Amplitude=m2s(0)
         Amplitude = Amplitude*BHBB
c     pp -> HAjjjLO -> A W W jj
      CASE(HAjjjLO_WW)
         rn = RandomNumber()
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         call m2s_wbfzhg( 16, 0, lokt,
     1        x, p(0,1,1), v(0,1,1), rn, xuz, m2s)
         Amplitude=m2s(0)
         call m2s_VVsum(v(0,1,1),3,4,2,1,m2s_dec)
         m2s_dec = m2s_dec/(16.0d0 * pi * xmgs(6))
         Amplitude = Amplitude*m2s_dec
c     pp -> HAjjjLO -> A Z Z jj 
      CASE(HAjjjLO_ZZ_ll)
         rn = RandomNumber()
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         call m2s_wbfzhg( 16, 0, lokt,
     1        x, p(0,1,1), v(0,1,1), rn, xuz, m2s)
         Amplitude=m2s(0)
         call m2s_VVsum(v(0,1,1),2,2,2,2,m2s_dec)
         m2s_dec = m2s_dec/(16.0d0 * pi * xmgs(6))
         Amplitude = Amplitude*m2s_dec
      CASE(HAjjjLO_ZZ_lnu)
         rn = RandomNumber()
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         call m2s_wbfzhg( 16, 0, lokt,
     1        x, p(0,1,1), v(0,1,1), rn, xuz, m2s)
         Amplitude=m2s(0)
         call m2s_VVsum(v(0,1,1),2,2,2,1,m2s_dec)
         Amplitude=m2s(0)
         m2s_dec = m2s_dec/(16.0d0 * pi * xmgs(6))
         Amplitude = Amplitude*m2s_dec  

#endif

#ifdef WITH_DIBOSON
*********************
*
*     DIBOSON : WW, WZ, ZZ, WA, ZA, AA, WH
*
*********************
c pp -> W+ W- -> ve e+ e- ve~ 
      CASE(WW,WPhadWMlep,WPlepWMhad)
         rn = RandomNumber()
         Amplitude = m2s_qqvv( 34, VBFNLOswitch(nlo), x, p, v, rn)
      CASE(GFWW,GFWPhadWMlep,GFWPlepWMhad)
         rn = RandomNumber()
         Amplitude = m2s_qqvv( 34, -8, x, p, v, rn)
c pp -> W+ Z -> ve e+ mu- mu+

      CASE(WPZ,WPhadZ,WPZhad)
         rn = RandomNumber()
         Amplitude = m2s_qqWZ( 32, nlo, x, p, v, rn)
c pp -> W- Z -> ve e- mu- mu+

      CASE(WMZ,WMhadZ,WMZhad)
         rn = RandomNumber()
         Amplitude = m2s_qqWZ( 42, nlo, x, p, v, rn)

c pp -> Z Z -> e+ e- mu+ mu-
      CASE(ZZ,ZZhad)
         rn = RandomNumber()
         Amplitude = m2s_qqzz( 22, VBFNLOswitch(nlo), x, p, v, rn)
      CASE(GFZZ,GFZZhad)
         rn = RandomNumber()
         Amplitude = m2s_qqzz( 22, -8, x, p, v, rn)

c pp -> W+ A -> ve e+ mu- mu+

      CASE(WPA)
         rn = RandomNumber()
         Amplitude = m2s_qqWA( 31, nlo, x, p, v, rn)

c pp -> W- A -> ve e+ mu- mu+

      CASE(WMA)
         rn = RandomNumber()
         Amplitude = m2s_qqWA( 41, nlo, x, p, v, rn)

c pp -> Z A -> e+ e- A
      CASE(ZA)
         rn = RandomNumber()
         Amplitude = m2s_qqza( 21, VBFNLOswitch(nlo), x, p, v, rn)
      CASE(GFZA)
         rn = RandomNumber()
         Amplitude = m2s_qqza( 21, -8, x, p, v, rn)

c pp -> A A 
      CASE(AA)
         rn = RandomNumber()
         Amplitude = m2s_qqAA( 11, VBFNLOswitch(nlo), x, p, v, rn)
      CASE(GFAA)
         rn = RandomNumber()
         Amplitude = m2s_qqAA( 11, -8, x, p, v, rn)

c pp -> W -> ve e
      CASE(WP_only)
         rn = RandomNumber()
         Amplitude = m2s_qqWH( 3, nlo, x, p, v, rn)

      CASE(WM_only)
         rn = RandomNumber()
         Amplitude = m2s_qqWH( 4, nlo, x, p, v, rn)

c     pp -> W H
      CASE(WPH)
         rn = RandomNumber()
         Amplitude = m2s_qqWH( 36, nlo, x, p, v, rn)
c     pp -> W H -> W gamma gamma
      CASE(WPH_AA)
         rn = RandomNumber()
         Amplitude = m2s_qqWH( 36, nlo, x, p, v, rn)
         Amplitude = Amplitude*BHGAM
c     pp -> W H -> W mu+mu-
      CASE(WPH_mu)
         rn = RandomNumber()
         Amplitude = m2s_qqWH( 36, nlo, x, p, v, rn)
         Amplitude = Amplitude*BHMU
c     pp -> W H -> W tau+tau-
      CASE(WPH_tau)
         rn = RandomNumber()
         Amplitude = m2s_qqWH( 36, nlo, x, p, v, rn)
         Amplitude = Amplitude*BHTAU
c     pp -> W H -> W b bbar
      CASE(WPH_bbar)
         rn = RandomNumber()
         Amplitude = m2s_qqWH( 36, nlo, x, p, v, rn)
         Amplitude = Amplitude*BHBB
c     pp -> W H -> W W W
      CASE(WPH_WW)
         rn = RandomNumber()
         Amplitude = m2s_qqWH( 36, nlo, x, p, v, rn)
         call m2s_VVsum(v(0,3,1),3,4,2,1,m2s_dec) !use v3..v6 for H->VV
         m2s_dec = m2s_dec/(16.0d0 * pi * xmgs(6))
         Amplitude = Amplitude*m2s_dec
!         print*, 'decfac', m2s_dec
c     pp -> W H -> W Z Z 
      CASE(WPH_ZZ_ll)
         rn = RandomNumber()
         Amplitude = m2s_qqWH( 36, nlo, x, p, v, rn)
         call m2s_VVsum(v(0,3,1),2,2,2,2,m2s_dec)
         m2s_dec = m2s_dec/(16.0d0 * pi * xmgs(6))
         Amplitude = Amplitude*m2s_dec
      CASE(WPH_ZZ_lnu)
         rn = RandomNumber()
         Amplitude = m2s_qqWH( 36, nlo, x, p, v, rn)
         call m2s_VVsum(v(0,3,1),2,2,2,1,m2s_dec)
         m2s_dec = m2s_dec/(16.0d0 * pi * xmgs(6))
         Amplitude = Amplitude*m2s_dec  


c     pp -> W H
      CASE(WMH)
         rn = RandomNumber()
         Amplitude = m2s_qqWH( 46, nlo, x, p, v, rn)
c     pp -> W H -> W gamma gamma
      CASE(WMH_AA)
         rn = RandomNumber()
         Amplitude = m2s_qqWH( 46, nlo, x, p, v, rn)
         Amplitude = Amplitude*BHGAM
c     pp -> W H -> W mu+mu-
      CASE(WMH_mu)
         rn = RandomNumber()
         Amplitude = m2s_qqWH( 46, nlo, x, p, v, rn)
         Amplitude = Amplitude*BHMU
c     pp -> W H -> W tau+tau-
      CASE(WMH_tau)
         rn = RandomNumber()
         Amplitude = m2s_qqWH( 46, nlo, x, p, v, rn)
         Amplitude = Amplitude*BHTAU
c     pp -> W H -> W b bbar
      CASE(WMH_bbar)
         rn = RandomNumber()
         Amplitude = m2s_qqWH( 46, nlo, x, p, v, rn)
         Amplitude = Amplitude*BHBB
c     pp -> W H -> W W W
      CASE(WMH_WW)
         rn = RandomNumber()
         Amplitude = m2s_qqWH( 46, nlo, x, p, v, rn)
         call m2s_VVsum(v(0,3,1),3,4,2,1,m2s_dec) !use v3..v6 for H->VV
         m2s_dec = m2s_dec/(16.0d0 * pi * xmgs(6))
         Amplitude = Amplitude*m2s_dec
!         print*, 'decfac', m2s_dec
c     pp -> W H -> W Z Z 
      CASE(WMH_ZZ_ll)
         rn = RandomNumber()
         Amplitude = m2s_qqWH( 46, nlo, x, p, v, rn)
         call m2s_VVsum(v(0,3,1),2,2,2,2,m2s_dec)
         m2s_dec = m2s_dec/(16.0d0 * pi * xmgs(6))
         Amplitude = Amplitude*m2s_dec
      CASE(WMH_ZZ_lnu)
         rn = RandomNumber()
         Amplitude = m2s_qqWH( 46, nlo, x, p, v, rn)
         call m2s_VVsum(v(0,3,1),2,2,2,1,m2s_dec)
         m2s_dec = m2s_dec/(16.0d0 * pi * xmgs(6))
         Amplitude = Amplitude*m2s_dec  


**********************
*
*     DIBOSON : WWj, WZj, ZZj, WAj, ZAj, AAj
*
**********************
c pp -> W+ W- j -> ve e+ e- ve~ j
      CASE(WWjLO,WPhadWMlepjLO,WPlepWMhadjLO)
         rn = RandomNumber()
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         call m2s_qqVVj( 34, 0, lokt,x, p, v, rn, xuz(1,2), m2s, ps_number)
         Amplitude = m2s(0)
       CASE(WPZjLOb,WPhadZjLO,WPZhadjLO)
         rn = RandomNumber()
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         call m2s_qqWZj( 32, 0, lokt,x, p, v, rn, xuz(1,2), m2s)
         Amplitude = m2s(0)
      CASE(WMZjLOb,WMhadZjLO,WMZhadjLO)
         rn = RandomNumber()
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         call m2s_qqWZj( 42, 0, lokt,x, p, v, rn, xuz(1,2), m2s)
         Amplitude = m2s(0)

c pp -> Z Z j -> e- e+ mu- mu+ j
      CASE(ZZjLO,ZZhadjLO)
         rn = RandomNumber()
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         call m2s_qqZZj( 22, 0, lokt,x, p, v, rn, xuz(1,2), m2s, ps_number)
         Amplitude = m2s(0)

      CASE(WPAjLOb)
         rn = RandomNumber()
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         call m2s_qqWAj( 31, 0, lokt,x, p, v, rn, xuz(1,2), m2s)
         Amplitude = m2s(0)
      CASE(WMAjLOb)
         rn = RandomNumber()
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         call m2s_qqWAj( 41, 0, lokt,x, p, v, rn, xuz(1,2), m2s)
         Amplitude = m2s(0)      

c pp -> W j -> e ve j

      CASE(WPjLOb)
         rn = RandomNumber()
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         call m2s_qqWHj( 3, 0, lokt,x, p, v, rn, xuz(1,2), m2s)
         Amplitude = m2s(0)
      CASE(WMjLOb)
         rn = RandomNumber()
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         call m2s_qqWHj( 4, 0, lokt,x, p, v, rn, xuz(1,2), m2s)
         Amplitude = m2s(0)      

c pp -> W H j -> e ve H j

c     pp -> WPHj
      CASE(WPHjLOb)
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         rn = RandomNumber()
         call m2s_qqWHj( 36, 0, lokt, x, p, v, rn, xuz(1,2), m2s)
         Amplitude = m2s(0)
c     pp -> WPHj -> gamma gamma j
      CASE(WPHjLOb_AA)
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         rn = RandomNumber()
         call m2s_qqWHj( 36, 0, lokt, x, p, v, rn, xuz(1,2), m2s)
         Amplitude = m2s(0)*BHGAM
c     pp -> WPHj -> mu+mu- j
      CASE(WPHjLOb_mu)
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         rn = RandomNumber()
         call m2s_qqWHj( 36, 0, lokt, x, p, v, rn, xuz(1,2), m2s)
         Amplitude = m2s(0)*BHMU
c     pp -> WPHj -> tau+tau- j
      CASE(WPHjLOb_tau)
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         rn = RandomNumber()
         call m2s_qqWHj( 36, 0, lokt, x, p, v, rn, xuz(1,2), m2s)
         Amplitude = m2s(0)*BHTAU
c     pp -> WPHj -> b bbar j
      CASE(WPHjLOb_bbar)
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         rn = RandomNumber()
         call m2s_qqWHj( 36, 0, lokt, x, p, v, rn, xuz(1,2), m2s)
         Amplitude = m2s(0)*BHBB
c     pp -> WPHj -> W W j
      CASE(WPHjLOb_WW)
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         rn = RandomNumber()
         call m2s_qqWHj( 36, 0, lokt, x, p, v, rn, xuz(1,2), m2s)
         call m2s_VVsum(v(0,3,1),3,4,2,1,m2s_dec)
         m2s_dec = m2s_dec/(16.0d0 * pi * xmgs(6))
         Amplitude = m2s(0)*m2s_dec
c     pp -> WPHj -> Z Z j
      CASE(WPHjLOb_ZZ_ll)
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         rn = RandomNumber()
         call m2s_qqWHj( 36, 0, lokt, x, p, v, rn, xuz(1,2), m2s)
         call m2s_VVsum(v(0,3,1),2,2,2,2,m2s_dec)
         m2s_dec = m2s_dec/(16.0d0 * pi * xmgs(6))
         Amplitude = m2s(0)*m2s_dec
      CASE(WPHjLOb_ZZ_lnu)
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         rn = RandomNumber()
         call m2s_qqWHj( 36, 0, lokt, x, p, v, rn, xuz(1,2), m2s)
         call m2s_VVsum(v(0,3,1),2,2,2,1,m2s_dec)
         m2s_dec = m2s_dec/(16.0d0 * pi * xmgs(6))
         Amplitude = m2s(0)*m2s_dec


c     pp -> WMHj
      CASE(WMHjLOb)
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         rn = RandomNumber()
         call m2s_qqWHj( 46, 0, lokt, x, p, v, rn, xuz(1,2), m2s)
         Amplitude = m2s(0)
c     pp -> WMHj -> gamma gamma j
      CASE(WMHjLOb_AA)
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         rn = RandomNumber()
         call m2s_qqWHj( 46, 0, lokt, x, p, v, rn, xuz(1,2), m2s)
         Amplitude = m2s(0)*BHGAM
c     pp -> WMHj -> mu+mu- j
      CASE(WMHjLOb_mu)
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         rn = RandomNumber()
         call m2s_qqWHj( 46, 0, lokt, x, p, v, rn, xuz(1,2), m2s)
         Amplitude = m2s(0)*BHMU
c     pp -> WMHj -> tau+tau- j
      CASE(WMHjLOb_tau)
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         rn = RandomNumber()
         call m2s_qqWHj( 46, 0, lokt, x, p, v, rn, xuz(1,2), m2s)
         Amplitude = m2s(0)*BHTAU
c     pp -> WMHj -> b bbar j
      CASE(WMHjLOb_bbar)
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         rn = RandomNumber()
         call m2s_qqWHj( 46, 0, lokt, x, p, v, rn, xuz(1,2), m2s)
         Amplitude = m2s(0)*BHBB
c     pp -> WMHj -> W W j
      CASE(WMHjLOb_WW)
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         rn = RandomNumber()
         call m2s_qqWHj( 46, 0, lokt, x, p, v, rn, xuz(1,2), m2s)
         call m2s_VVsum(v(0,3,1),3,4,2,1,m2s_dec)
         m2s_dec = m2s_dec/(16.0d0 * pi * xmgs(6))
         Amplitude = m2s(0)*m2s_dec
c     pp -> WMHj -> Z Z j
      CASE(WMHjLOb_ZZ_ll)
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         rn = RandomNumber()
         call m2s_qqWHj( 46, 0, lokt, x, p, v, rn, xuz(1,2), m2s)
         call m2s_VVsum(v(0,3,1),2,2,2,2,m2s_dec)
         m2s_dec = m2s_dec/(16.0d0 * pi * xmgs(6))
         Amplitude = m2s(0)*m2s_dec
      CASE(WMHjLOb_ZZ_lnu)
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         rn = RandomNumber()
         call m2s_qqWHj( 46, 0, lokt, x, p, v, rn, xuz(1,2), m2s)
         call m2s_VVsum(v(0,3,1),2,2,2,1,m2s_dec)
         m2s_dec = m2s_dec/(16.0d0 * pi * xmgs(6))
         Amplitude = m2s(0)*m2s_dec

!      CASE(WPHjLOb)
!         rn = RandomNumber()
!         lokt(1) = .true.
!         lokt(2) = .false.
!         lokt(3) = .false.
!         call m2s_qqWHj( 36, 0, lokt,x, p, v, rn, xuz(1,2), m2s)
!         Amplitude = m2s(0)
!      CASE(WMHjLOb)
!         rn = RandomNumber()
!         lokt(1) = .true.
!         lokt(2) = .false.
!         lokt(3) = .false.
!         call m2s_qqWHj( 46, 0, lokt,x, p, v, rn, xuz(1,2), m2s)
!         Amplitude = m2s(0)      

c pp -> Z A j -> e- e+ A j
      CASE(ZAjLO)
         rn = RandomNumber()
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         call m2s_qqZAj( 21, 0, lokt,x, p, v, rn, xuz(1,2), m2s, ps_number)
         Amplitude = m2s(0)

c pp -> A A j 
      CASE(AAjLO)
         rn = RandomNumber()
         lokt(1) = .true.
         do i=2,n_kin
           lokt(i) = .false.
         enddo
         call m2s_qqAAj( 11, 0, lokt,x, p, v, rn, xuz(1,2), m2s, ps_number)
         Amplitude = m2s(0)

#endif          

#ifdef WITH_DIBOSONJET
**********************
*
*     DIBOSON PLUS JET : WWj, WZj, ZZj, WAj
*
**********************
      CASE(WWJ,WPHADWMJ,WPWMHADJ)
         rn = RandomNumber()
         Amplitude = m2s_qqvvj2( 34, VBFNLOswitch(nlo), x, p, v, rn, ps_number)
      CASE(GFWWJ,GFWPHADWMJ,GFWPWMHADJ)
         rn = RandomNumber()
         Amplitude = m2s_qqvvj2( 34, -8, x, p, v, rn, ps_number)
      CASE(WMAJ)
         if (nlo.eq.0) then
            Amplitude = m2s_qqwgaglu_hel(x,p(0,1,1),v(0,1,1),nlo,.false.)
         else 
#ifdef WITH_NLO
            Call m2s_ppwgajvirt(x,rn,.false.,p(0,1,1),v(0,1,1), m2s(0),nlo)
            Amplitude=m2s(0)
#endif
         endif
      CASE(WPAJ)
         if (nlo.eq.0) then
            Amplitude = m2s_qqwpgaglu_hel(x,p(0,1,1),v(0,1,1),nlo,.true.)
         else 
#ifdef WITH_NLO
            Call m2s_ppwpgajvirt(x,rn,.true.,p(0,1,1),v(0,1,1), m2s(0),nlo)
            Amplitude=m2s(0)
#endif
         endif
      CASE(WPZJ, WPhadZJ, WPZhadJ)
         if (nlo.eq.0) then
            Amplitude = m2s_qqwpzglu_hel(x,p(0,1,1),v(0,1,1),nlo,.true.)
         else 
#ifdef WITH_NLO
            Call m2s_ppwpzjvirt(x,rn,.true.,p(0,1,1),v(0,1,1), m2s(0),nlo)
            Amplitude=m2s(0)
#endif
         endif
      CASE(WMZJ, WMhadZJ, WMZhadJ)
         if (nlo.eq.0) then
            Amplitude = m2s_qqwzglu_hel(x,p(0,1,1),v(0,1,1),nlo,.false.)
         else 
#ifdef WITH_NLO
            Call m2s_ppwzjvirt(x,rn,.false.,p(0,1,1),v(0,1,1), m2s(0),nlo)
            Amplitude=m2s(0)
#endif
         endif
      CASE(ZZJ)
         rn = RandomNumber()
         Amplitude = m2s_qqzzj2( 22, VBFNLOswitch(nlo), x, p, v, rn, ps_number)
      CASE(GFZZJ)
         rn = RandomNumber()
         Amplitude = m2s_qqzzj2( 22, -8, x, p, v, rn, ps_number)

      CASE(WWjjLO,WPHADWMJJLO,WPWMHADJJLO)
         lokt(1) = .true.
         do i=2,max_kin
         lokt(i) = .false.
         enddo
         call m2s_qqVVjj( 34, 0, lokt,x, p, v, rn, xuz, m2s)
         Amplitude = m2s(0)
      CASE(WMAJLO)
         lokt(1) = .true.
         do i=2,max_kin
         lokt(i) = .false.
         enddo
         CALL m2s_qqwgagluj(xuz,0,.false.,lokt,x,p,v,m2s)
         Amplitude = m2s(0)
      CASE(WPAJLO)
         lokt(1) = .true.
         do i=2,max_kin
         lokt(i) = .false.
         enddo
         CALL m2s_qqwpgagluj(xuz,0,.true.,lokt,x,p,v,m2s)
         Amplitude = m2s(0)
      CASE(WMZJLO, WMhadZJjLO, WMZhadJjLO)
         lokt(1) = .true.
         do i=2,max_kin
         lokt(i) = .false.
         enddo
         CALL m2s_qqwzgluj(xuz,0,.false.,lokt,x,p,v,m2s)
         Amplitude = m2s(0)
      CASE(WPZJLO, WPhadZJjLO, WPZhadJjLO)
         lokt(1) = .true.
         do i=2,max_kin
         lokt(i) = .false.
         enddo
         CALL m2s_qqwpzgluj(xuz,0,.true.,lokt,x,p,v,m2s)
         Amplitude = m2s(0)
      CASE(ZZjjLO)
         lokt(1) = .true.
         do i=2,max_kin
         lokt(i) = .false.
         enddo
         call m2s_qqZZjj( 22, 0, lokt,x, p, v, rn, xuz, m2s)
         Amplitude = m2s(0)
#endif

#ifdef WITH_TRIBOSON
**********************************************************************
*
*     TRIBOSON : WWZ, ZZZ, ZZW, WWW, WWA, ZZA, WZA, WAA, ZAA, AAA
*
**********************************************************************         
c pp -> W+ W- Z -> ve e+ e- ve~ e- e+
      CASE(WWZ, WPhadWMZ, WPWMhadZ, WWZhad)
         rn = RandomNumber()
         Amplitude = m2s_qqWWZ(342,VBFNLOswitch(nlo),x,p,v,rn,narrowwidth)

c pp -> Z Z Z -> e- e+ e- e+ e- e+
      CASE(ZZZ, ZZZhad)
         rn = RandomNumber()
         Amplitude = m2s_qqZZZ(222,VBFNLOswitch(nlo),x,p,v,rn,
     1                                       narrowwidth,ps_number)

c pp -> Z Z W+ -> e- e+ e- e+ ve e+
      CASE(ZZWP, ZZhadWP, ZZWPhad)
         rn = RandomNumber()
         Amplitude = m2s_qqZZW(223,VBFNLOswitch(nlo),x,p,v,rn,narrowwidth)

c pp -> Z Z W- -> e- e+ e- e+ e- ve~
      CASE(ZZWM, ZZhadWM, ZZWMhad)
         rn = RandomNumber()
         Amplitude = m2s_qqZZW(224,VBFNLOswitch(nlo),x,p,v,rn,narrowwidth)

c pp -> W+W-W+ -> ve e+ e- ve~ ve e+
      CASE(WWWP, WPhadWMWP, WPWMhadWP)
         rn = RandomNumber()
         Amplitude = m2s_qqWWW( 343,VBFNLOswitch(nlo),x,p,v,rn,narrowwidth)

c pp -> W+W-W- -> e- ve~ ve e+ e- ve~
      CASE(WWWM, WMhadWPWM, WMWPhadWM)
         rn = RandomNumber()
         Amplitude = m2s_qqWWW( 434,VBFNLOswitch(nlo),x,p,v,rn,narrowwidth)

c pp -> W+ W- A -> ve e+ e- ve~ A
      CASE(WWA, WPhadWMA, WPWMhadA)
         rn = RandomNumber()
         Amplitude = m2s_qqWWA( 341, VBFNLOswitch(nlo), x, p, v, rn)

c pp -> Z Z A -> e- e+ e- e+ A
      CASE(ZZA, ZZhadA)
         rn = RandomNumber()
         Amplitude = m2s_qqZZA( 221, VBFNLOswitch(nlo), x, p, v, rn)

c pp -> Z Z A -> e- e+ ve ve~ A
      CASE(ZZnA)
         rn = RandomNumber()
         Amplitude = m2s_qqZZA( 251, VBFNLOswitch(nlo), x, p, v, rn)

c pp -> W+ Z A -> e+ ve e- e+ A
      CASE(WPZA, WPhadZA, WPZhadA)
         rn = RandomNumber()
         Amplitude = m2s_qqWZA( 321, VBFNLOswitch(nlo), x, p, v, rn, ps_number)

c pp -> W- Z A -> e- ve~ e- e+ A
      CASE(WMZA, WMhadZA, WMZhadA)
         rn = RandomNumber()
         Amplitude = m2s_qqWZA( 421, VBFNLOswitch(nlo), x, p, v, rn, ps_number)

c pp -> W+ A A -> e+ ve A A
      CASE(WPAA)
         rn = RandomNumber()
         Amplitude = m2s_qqWAA( 311, VBFNLOswitch(nlo), x, p, v, rn, 
     &        ps_number, narrowwidth)

c pp -> W- A A -> e- ve~ A A
      CASE(WMAA)
         rn = RandomNumber()
         Amplitude = m2s_qqWAA( 411, VBFNLOswitch(nlo), x, p, v, rn, 
     &        ps_number, narrowwidth)

c pp -> Z A A -> e+ e- A A
      CASE(ZAA)
         rn = RandomNumber()
         Amplitude = m2s_qqZAA( 211, VBFNLOswitch(nlo), x, p, v, rn, 
     &        ps_number, narrowwidth)

c pp -> Z A A -> ve~ ve A A
      CASE(ZnAA)
         rn = RandomNumber()
         Amplitude = m2s_qqZAA( 511, VBFNLOswitch(nlo), x, p, v, rn, 
     &        ps_number, narrowwidth)

c pp -> A A A 
      CASE(AAA)
         rn = RandomNumber()
         Amplitude = m2s_qqAAA( 111, VBFNLOswitch(nlo), x, p, v, rn, ps_number)

*****************************************************************************
*
*     TRIBOSON : WWZj, ZZZj, ZZWj, WWWj, WWAj, ZZAj, WZAj, WAAj, ZAAj, AAAj
*
*****************************************************************************
c pp -> W+ W- Z j -> ve e+ e- ve~ e- e+ j
      CASE(WWZjLO, WPhadWMZjLO, WPWMhadZjLO, WWZhadjLO)
         rn = RandomNumber()
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         call m2s_qqWWZj(342,0,lokt,x,p,v,rn,xuz(1,2),m2s,narrowwidth)
         Amplitude = m2s(0)

c pp -> Z Z Z j -> e- e+ e- e+ e- e+ j
      CASE(ZZZjLO, ZZZhadjLO)
         rn = RandomNumber()
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         call m2s_qqZZZj(222, 0, lokt, x, p, v, rn, xuz(1,2), m2s,
     1                       narrowwidth, ps_number)
         Amplitude = m2s(0)

c pp -> Z Z W+ j -> e- e+ e- e+ ve e+ j
      CASE(ZZWPjLO, ZZhadWPjLO, ZZWPhadjLO)
         rn = RandomNumber()
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         call m2s_qqZZWj(223,0,lokt,x,p,v,rn,xuz(1,2),m2s,narrowwidth)
         Amplitude = m2s(0)

c pp -> Z Z W- j -> e- e+ e- e+ e- ve~ j
      CASE(ZZWMjLO, ZZhadWMjLO, ZZWMhadjLO)
         rn = RandomNumber()
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         call m2s_qqZZWj(224,0,lokt,x,p,v,rn,xuz(1,2),m2s,narrowwidth)
         Amplitude = m2s(0)

c pp -> W+W-W+j -> ve e+ e- ve~ vmu mu+ j
      CASE(WWWPjLO, WPhadWMWPjLO, WPWMhadWPjLO)
         rn = RandomNumber()
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         call m2s_qqWWWj(343,0,lokt,x,p,v,rn,xuz(1,2),m2s,narrowwidth)
         Amplitude = m2s(0)

c pp -> W+W-W-j -> e- ve~ ve e+ e- ve~ j
      CASE(WWWMjLO, WMhadWPWMjLO, WMWPhadWMjLO)
         rn = RandomNumber()
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         call m2s_qqWWWj(434,0,lokt,x,p,v,rn,xuz(1,2),m2s,narrowwidth)
         Amplitude = m2s(0)

c pp -> W+ W- A j -> ve e+ e- ve~ A j
      CASE(WWAjLO, WPhadWMAjLO, WPWMhadAjLO)
         rn = RandomNumber()
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         call m2s_qqVVV1j( 341, 0, lokt,x, p, v, rn, xuz(1,2), m2s)
         Amplitude = m2s(0)

c pp -> Z Z A j -> e- e+ e- e+ A j
      CASE(ZZAjLO, ZZhadAjLO)
         rn = RandomNumber()
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         call m2s_qqZZAj( 221, 0, lokt,x, p, v, rn, xuz(1,2), m2s)
         Amplitude = m2s(0)

c pp -> Z Z A j -> e- e+ ve ve~ A j
      CASE(ZZnAjLO)
         rn = RandomNumber()
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         call m2s_qqZZAj( 251, 0, lokt,x, p, v, rn, xuz(1,2), m2s)
         Amplitude = m2s(0)

c pp -> W+ Z A j -> e+ ve e- e+ A j
      CASE(WPZAjLO, WPhadZAjLO, WPZhadAjLO)
         rn = RandomNumber()
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         call m2s_qqWZAj( 321, 0, lokt,x, p, v, rn, xuz(1,2), m2s, ps_number)
         Amplitude = m2s(0)

c pp -> W- Z A j -> e- ve~ e- e+ A j
      CASE(WMZAjLO, WMhadZAjLO, WMZhadAjLO)
         rn = RandomNumber()
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         call m2s_qqWZAj( 421, 0, lokt,x, p, v, rn, xuz(1,2), m2s, ps_number)
         Amplitude = m2s(0)

c pp -> W+ A A j -> e+ ve A A j
      CASE(WPAAjLO)
         rn = RandomNumber()
         lokt(1) = .true.
         do i = 2, 11
            lokt(i) = .false.
         end do
         call m2s_qqWAAj( 311, 0, lokt,x, p, v, rn, xuz(1,2), m2s, 
     -        ps_number, narrowwidth)
         Amplitude = m2s(0)

c pp -> W- A A j -> e- ve~ A A j
      CASE(WMAAjLO)
         rn = RandomNumber()
         lokt(1) = .true.
         do i = 2, 11
            lokt(i) = .false.
         end do
         call m2s_qqWAAj( 411, 0, lokt,x, p, v, rn, xuz(1,2), m2s, 
     -        ps_number, narrowwidth)
         Amplitude = m2s(0)

c pp -> Z A A j -> e+ e- A A j
      CASE(ZAAjLO)
         rn = RandomNumber()
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         call m2s_qqZAAj( 211, 0, lokt,x, p, v, rn, xuz(1,2), m2s, 
     -        ps_number, narrowwidth)
         Amplitude = m2s(0)

c pp -> Z A A j -> ve~ ve A A j
      CASE(ZnAAjLO)
         rn = RandomNumber()
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         call m2s_qqZAAj( 511, 0, lokt,x, p, v, rn, xuz(1,2), m2s, 
     &        ps_number, narrowwidth)
         Amplitude = m2s(0)

c pp -> A A A j 
      CASE(AAAjLO)
         rn = RandomNumber()
         lokt(1) = .true.
         lokt(2) = .false.
         lokt(3) = .false.
         call m2s_qqAAAj( 111, 0, lokt,x, p, v, rn, xuz(1,2), m2s, ps_number)
         Amplitude = m2s(0)

#endif

**************************************
*
*     VVVj 
*
**************************************
#ifdef WITH_TRIBOSONJET
c pp -> W+ A A j -> e+ ve A A j
      CASE(WPAAj)
         rn = RandomNumber()
         Amplitude = m2s_qqWAAj2( 311, VBFNLOswitch(nlo), x, p, v, rn, 
     &        ps_number, narrowwidth)

c pp -> W- A A j -> e- ve~ A A j
      CASE(WMAAj)
         rn = RandomNumber()
         Amplitude = m2s_qqWAAj2( 411, VBFNLOswitch(nlo), x, p, v, rn,
     &        ps_number, narrowwidth)

c pp -> W+ A A j j -> e+ ve A A j j
      CASE(WPAAjjLO)
         rn = RandomNumber()
         lokt( 1) = .true.
         do i = 2, 11
            lokt(i) = .false.
         end do
         call m2s_qqWAAjj( 311, 0, lokt,x, p, v, rn, xuz(1,2), 
     &        ps_number, narrowwidth, m2s)
         Amplitude = m2s(0)

c pp -> W- A A j j -> e- ve~ A A j j
      CASE(WMAAjjLO)
         rn = RandomNumber()
         lokt( 1) = .true.
         do i = 2, 11
            lokt(i) = .false.
         end do
         call m2s_qqWAAjj( 411, 0, lokt,x, p, v, rn, xuz(1,2), 
     &        ps_number, narrowwidth, m2s)
         Amplitude = m2s(0)
#endif


#ifdef WITH_DIBOSONJET
c pp -> W+ j -> e+ ve j
      CASE(WPj)
         rn = RandomNumber()
         Amplitude = m2s_qqWHj2( 3, VBFNLOswitch(nlo), x, p, v, rn, 
     &        ps_number, narrowwidth)

c pp -> W- j -> e- ve~ j
      CASE(WMj)
         rn = RandomNumber()
         Amplitude = m2s_qqWHj2( 4, VBFNLOswitch(nlo), x, p, v, rn,
     &        ps_number, narrowwidth)

c pp -> W+ j j -> e+ ve j j
      CASE(WPjjLO)
         rn = RandomNumber()
         lokt( 1) = .true.
         do i = 2, 11
            lokt(i) = .false.
         end do
         call m2s_qqWHjj( 3, 0, lokt,x, p, v, rn, xuz(1,2), 
     &        narrowwidth, 0, m2s)
         Amplitude = m2s(0)

c pp -> W- j j -> e- ve~ j j
      CASE(WMjjLO)
         rn = RandomNumber()
         lokt( 1) = .true.
         do i = 2, 11
            lokt(i) = .false.
         end do
         call m2s_qqWHjj( 4, 0, lokt,x, p, v, rn, xuz(1,2), 
     &        narrowwidth, 0, m2s)
         Amplitude = m2s(0)


c     pp -> W HJ
      CASE(WPHJ)
         rn = RandomNumber()
         Amplitude = m2s_qqWHj2( 36, VBFNLOswitch(nlo), x, p, v, rn, narrowwidth)
c     pp -> W HJ -> W gamma gamma
      CASE(WPHJ_AA)
         rn = RandomNumber()
         Amplitude = m2s_qqWHj2( 36, VBFNLOswitch(nlo), x, p, v, rn, narrowwidth)
         Amplitude = Amplitude*BHGAM
c     pp -> W HJ -> W mu+mu-
      CASE(WPHJ_mu)
         rn = RandomNumber()
         Amplitude = m2s_qqWHj2( 36, VBFNLOswitch(nlo), x, p, v, rn, narrowwidth)
         Amplitude = Amplitude*BHMU
c     pp -> W HJ -> W tau+tau-
      CASE(WPHJ_tau)
         rn = RandomNumber()
         Amplitude = m2s_qqWHj2( 36, VBFNLOswitch(nlo), x, p, v, rn, narrowwidth)
         Amplitude = Amplitude*BHTAU
c     pp -> W HJ -> W b bbar
      CASE(WPHJ_bbar)
         rn = RandomNumber()
         Amplitude = m2s_qqWHj2( 36, VBFNLOswitch(nlo), x, p, v, rn, narrowwidth)
         Amplitude = Amplitude*BHBB
c     pp -> W HJ -> W W W
      CASE(WPHJ_WW)
         rn = RandomNumber()
         Amplitude = m2s_qqWHj2( 36, VBFNLOswitch(nlo), x, p, v, rn, narrowwidth)
         call m2s_VVsum(v(0,3,1),3,4,2,1,m2s_dec) !use v3..v6 for H->VV
         m2s_dec = m2s_dec/(16.0d0 * pi * xmgs(6))
         Amplitude = Amplitude*m2s_dec
!         print*, 'decfac', m2s_dec
c     pp -> W HJ -> W Z Z 
      CASE(WPHJ_ZZ_ll)
         rn = RandomNumber()
         Amplitude = m2s_qqWHj2( 36, VBFNLOswitch(nlo), x, p, v, rn, narrowwidth)
         call m2s_VVsum(v(0,3,1),2,2,2,2,m2s_dec)
         m2s_dec = m2s_dec/(16.0d0 * pi * xmgs(6))
         Amplitude = Amplitude*m2s_dec
      CASE(WPHJ_ZZ_lnu)
         rn = RandomNumber()
         Amplitude = m2s_qqWHj2( 36, VBFNLOswitch(nlo), x, p, v, rn, narrowwidth)
         call m2s_VVsum(v(0,3,1),2,2,2,1,m2s_dec)
         m2s_dec = m2s_dec/(16.0d0 * pi * xmgs(6))
         Amplitude = Amplitude*m2s_dec  


c     pp -> W HJ
      CASE(WMHJ)
         rn = RandomNumber()
         Amplitude = m2s_qqWHj2( 46, VBFNLOswitch(nlo), x, p, v, rn, narrowwidth)
c     pp -> W HJ -> W gamma gamma
      CASE(WMHJ_AA)
         rn = RandomNumber()
         Amplitude = m2s_qqWHj2( 46, VBFNLOswitch(nlo), x, p, v, rn, narrowwidth)
         Amplitude = Amplitude*BHGAM
c     pp -> W HJ -> W mu+mu-
      CASE(WMHJ_mu)
         rn = RandomNumber()
         Amplitude = m2s_qqWHj2( 46, VBFNLOswitch(nlo), x, p, v, rn, narrowwidth)
         Amplitude = Amplitude*BHMU
c     pp -> W HJ -> W tau+tau-
      CASE(WMHJ_tau)
         rn = RandomNumber()
         Amplitude = m2s_qqWHj2( 46, VBFNLOswitch(nlo), x, p, v, rn, narrowwidth)
         Amplitude = Amplitude*BHTAU
c     pp -> W HJ -> W b bbar
      CASE(WMHJ_bbar)
         rn = RandomNumber()
         Amplitude = m2s_qqWHj2( 46, VBFNLOswitch(nlo), x, p, v, rn, narrowwidth)
         Amplitude = Amplitude*BHBB
c     pp -> W HJ -> W W W
      CASE(WMHJ_WW)
         rn = RandomNumber()
         Amplitude = m2s_qqWHj2( 46, VBFNLOswitch(nlo), x, p, v, rn, narrowwidth)
         call m2s_VVsum(v(0,3,1),3,4,2,1,m2s_dec) !use v3..v6 for H->VV
         m2s_dec = m2s_dec/(16.0d0 * pi * xmgs(6))
         Amplitude = Amplitude*m2s_dec
!         print*, 'decfac', m2s_dec
c     pp -> W HJ -> W Z Z 
      CASE(WMHJ_ZZ_ll)
         rn = RandomNumber()
         Amplitude = m2s_qqWHj2( 46, VBFNLOswitch(nlo), x, p, v, rn, narrowwidth)
         call m2s_VVsum(v(0,3,1),2,2,2,2,m2s_dec)
         m2s_dec = m2s_dec/(16.0d0 * pi * xmgs(6))
         Amplitude = Amplitude*m2s_dec
      CASE(WMHJ_ZZ_lnu)
         rn = RandomNumber()
         Amplitude = m2s_qqWHj2( 46, VBFNLOswitch(nlo), x, p, v, rn, narrowwidth)
         call m2s_VVsum(v(0,3,1),2,2,2,1,m2s_dec)
         m2s_dec = m2s_dec/(16.0d0 * pi * xmgs(6))
         Amplitude = Amplitude*m2s_dec


c     pp -> W HJ
      CASE(WPHJJLO)
         rn = RandomNumber()
         lokt( 1) = .true.
         do i = 2, 11
            lokt(i) = .false. 
         end do
         call m2s_qqWHjj( 36, 0, lokt,x, p, v, rn, xuz(1,2), 
     &                           narrowwidth, 0, m2s)
         Amplitude=m2s(0)
c     pp -> W HJ -> W gamma gamma
      CASE(WPHJJLO_AA)
         rn = RandomNumber()
         lokt( 1) = .true.
         do i = 2, 11
            lokt(i) = .false. 
         end do
         call m2s_qqWHjj( 36, 0, lokt,x, p, v, rn, xuz(1,2), 
     &                           narrowwidth, 0, m2s)
         Amplitude=m2s(0)
         Amplitude = Amplitude*BHGAM
c     pp -> W HJ -> W mu+mu-
      CASE(WPHJJLO_mu)
         rn = RandomNumber()
         lokt( 1) = .true.
         do i = 2, 11
            lokt(i) = .false. 
         end do
         call m2s_qqWHjj( 36, 0, lokt,x, p, v, rn, xuz(1,2), 
     &                           narrowwidth, 0, m2s)
         Amplitude=m2s(0)
         Amplitude = Amplitude*BHMU
c     pp -> W HJ -> W tau+tau-
      CASE(WPHJJLO_tau)
         rn = RandomNumber()
         lokt( 1) = .true.
         do i = 2, 11
            lokt(i) = .false. 
         end do
         call m2s_qqWHjj( 36, 0, lokt,x, p, v, rn, xuz(1,2), 
     &                           narrowwidth, 0, m2s)
         Amplitude=m2s(0)
         Amplitude = Amplitude*BHTAU
c     pp -> W HJ -> W b bbar
      CASE(WPHJJLO_bbar)
         rn = RandomNumber()
         lokt( 1) = .true.
         do i = 2, 11
            lokt(i) = .false. 
         end do
         call m2s_qqWHjj( 36, 0, lokt,x, p, v, rn, xuz(1,2), 
     &                           narrowwidth, 0, m2s)
         Amplitude=m2s(0)
         Amplitude = Amplitude*BHBB
c     pp -> W HJ -> W W W
      CASE(WPHJJLO_WW)
         rn = RandomNumber()
         lokt( 1) = .true.
         do i = 2, 11
            lokt(i) = .false. 
         end do
         call m2s_qqWHjj( 36, 0, lokt,x, p, v, rn, xuz(1,2), 
     &                           narrowwidth, 0, m2s)
         Amplitude=m2s(0)
         call m2s_VVsum(v(0,3,1),3,4,2,1,m2s_dec) !use v3..v6 for H->VV
         m2s_dec = m2s_dec/(16.0d0 * pi * xmgs(6))
         Amplitude = Amplitude*m2s_dec
!         print*, 'decfac', m2s_dec
c     pp -> W HJ -> W Z Z 
      CASE(WPHJJLO_ZZ_ll)
         rn = RandomNumber()
         lokt( 1) = .true.
         do i = 2, 11
            lokt(i) = .false. 
         end do
         call m2s_qqWHjj( 36, 0, lokt,x, p, v, rn, xuz(1,2), 
     &                           narrowwidth, 0, m2s)
         Amplitude=m2s(0)
         call m2s_VVsum(v(0,3,1),2,2,2,2,m2s_dec)
         m2s_dec = m2s_dec/(16.0d0 * pi * xmgs(6))
         Amplitude = Amplitude*m2s_dec
      CASE(WPHJJLO_ZZ_lnu)
         rn = RandomNumber()
         lokt( 1) = .true.
         do i = 2, 11
            lokt(i) = .false. 
         end do
         call m2s_qqWHjj( 36, 0, lokt,x, p, v, rn, xuz(1,2), 
     &                           narrowwidth, 0, m2s)
         Amplitude=m2s(0)
         call m2s_VVsum(v(0,3,1),2,2,2,1,m2s_dec)
         m2s_dec = m2s_dec/(16.0d0 * pi * xmgs(6))
         Amplitude = Amplitude*m2s_dec  


c     pp -> W HJ
      CASE(WMHJJLO)
         rn = RandomNumber()
         lokt( 1) = .true.
         do i = 2, 11
            lokt(i) = .false. 
         end do
         call m2s_qqWHjj( 46, 0, lokt,x, p, v, rn, xuz(1,2), 
     &                           narrowwidth, 0, m2s)
         Amplitude=m2s(0)
c     pp -> W HJ -> W gamma gamma
      CASE(WMHJJLO_AA)
         rn = RandomNumber()
         lokt( 1) = .true.
         do i = 2, 11
            lokt(i) = .false. 
         end do
         call m2s_qqWHjj( 46, 0, lokt,x, p, v, rn, xuz(1,2), 
     &                           narrowwidth, 0, m2s)
         Amplitude=m2s(0)
         Amplitude = Amplitude*BHGAM
c     pp -> W HJ -> W mu+mu-
      CASE(WMHJJLO_mu)
         rn = RandomNumber()
         lokt( 1) = .true.
         do i = 2, 11
            lokt(i) = .false. 
         end do
         call m2s_qqWHjj( 46, 0, lokt,x, p, v, rn, xuz(1,2), 
     &                           narrowwidth, 0, m2s)
         Amplitude=m2s(0)
         Amplitude = Amplitude*BHMU
c     pp -> W HJ -> W tau+tau-
      CASE(WMHJJLO_tau)
         rn = RandomNumber()
         lokt( 1) = .true.
         do i = 2, 11
            lokt(i) = .false. 
         end do
         call m2s_qqWHjj( 46, 0, lokt,x, p, v, rn, xuz(1,2), 
     &                           narrowwidth, 0, m2s)
         Amplitude=m2s(0)
         Amplitude = Amplitude*BHTAU
c     pp -> W HJ -> W b bbar
      CASE(WMHJJLO_bbar)
         rn = RandomNumber()
         lokt( 1) = .true.
         do i = 2, 11
            lokt(i) = .false. 
         end do
         call m2s_qqWHjj( 46, 0, lokt,x, p, v, rn, xuz(1,2), 
     &                           narrowwidth, 0, m2s)
         Amplitude=m2s(0)
         Amplitude = Amplitude*BHBB
c     pp -> W HJ -> W W W
      CASE(WMHJJLO_WW)
         rn = RandomNumber()
         lokt( 1) = .true.
         do i = 2, 11
            lokt(i) = .false. 
         end do
         call m2s_qqWHjj( 46, 0, lokt,x, p, v, rn, xuz(1,2), 
     &                           narrowwidth, 0, m2s)
         Amplitude=m2s(0)
         call m2s_VVsum(v(0,3,1),3,4,2,1,m2s_dec) !use v3..v6 for H->VV
         m2s_dec = m2s_dec/(16.0d0 * pi * xmgs(6))
         Amplitude = Amplitude*m2s_dec
!         print*, 'decfac', m2s_dec
c     pp -> W HJ -> W Z Z 
      CASE(WMHJJLO_ZZ_ll)
         rn = RandomNumber()
         lokt( 1) = .true.
         do i = 2, 11
            lokt(i) = .false. 
         end do
         call m2s_qqWHjj( 46, 0, lokt,x, p, v, rn, xuz(1,2), 
     &                           narrowwidth, 0, m2s)
         Amplitude=m2s(0)
         call m2s_VVsum(v(0,3,1),2,2,2,2,m2s_dec)
         m2s_dec = m2s_dec/(16.0d0 * pi * xmgs(6))
         Amplitude = Amplitude*m2s_dec
      CASE(WMHJJLO_ZZ_lnu)
         rn = RandomNumber()
         lokt( 1) = .true.
         do i = 2, 11
            lokt(i) = .false. 
         end do
         call m2s_qqWHjj( 46, 0, lokt,x, p, v, rn, xuz(1,2), 
     &                           narrowwidth, 0, m2s)
         Amplitude=m2s(0)
         call m2s_VVsum(v(0,3,1),2,2,2,1,m2s_dec)
         m2s_dec = m2s_dec/(16.0d0 * pi * xmgs(6))
         Amplitude = Amplitude*m2s_dec  

!c pp -> W+ H j -> e+ ve H j
!      CASE(WPHj)
!         rn = RandomNumber()
!         Amplitude = m2s_qqWHj2( 36, VBFNLOswitch(nlo), x, p, v, rn, 
!     &        ps_number, narrowwidth)
!
!c pp -> W- H j -> e- ve~ H j
!      CASE(WMHj)
!         rn = RandomNumber()
!         Amplitude = m2s_qqWHj2( 46, VBFNLOswitch(nlo), x, p, v, rn,
!     &        ps_number, narrowwidth)
!


!c pp -> W+ H j j -> e+ ve H j j
!      CASE(WPHjjLO)
!         rn = RandomNumber()
!         lokt( 1) = .true.
!         do i = 2, 11
!            lokt(i) = .false.
!         end do
!         call m2s_qqWHjj( 36, 0, lokt,x, p, v, rn, xuz(1,2), 
!     &        narrowwidth, m2s)
!         Amplitude = m2s(0)
!
!c pp -> W- H j j -> e- ve~ H j j
!      CASE(WMHjjLO)
!         rn = RandomNumber()
!         lokt( 1) = .true.
!         do i = 2, 11
!            lokt(i) = .false.
!         end do
!         call m2s_qqWHjj( 46, 0, lokt,x, p, v, rn, xuz(1,2), 
!     &        narrowwidth, m2s)
!         Amplitude = m2s(0)
#endif      

#ifdef WITH_GGF
**************************************
*
*     GF : Hjj 
*
**************************************
c     GF pp -> Hjj
      CASE(GFHjj)
         do mu = 0,3 !     sum the 2 lepton momenta up to Higgs momentum
            q(mu) = v(mu,1,1)+v(mu,2,1)
         enddo
         Amplitude = m2s_ggfh( x, p(0,1,1),q,rand )
*SEGorig         Amplitude = m2s_ggfh( x, p(0,1,1),v(0,1,1),rand )
c     GF pp -> Hjj -> gamma gamma jj
      CASE(GFHjj_AA)
         do mu = 0,3 !     sum the 2 lepton momenta up to Higgs momentum
            q(mu) = v(mu,1,1)+v(mu,2,1)
         enddo
         m2s_dec = q(0)*q(0) - q(1)*q(1) - q(2)*q(2) - q(3)*q(3)
         amplitude = m2s_ggfh( x, p(0,1,1),q,rand ) ! pp-> Hjj
         AMplitude = Amplitude*BHGAM
         Amplitude = Amplitude*16D0*pi*XMGS(6)/
     -        ((m2s_dec - XM2S(6))**2 + XMGS(6)**2)
c     GF pp -> Hjj -> mu+ mu- jj
      CASE(GFHjj_mu)  
         do mu = 0,3 !     sum the 2 lepton momenta up to Higgs momentum
            q(mu) = v(mu,1,1)+v(mu,2,1)
         enddo
         m2s_dec = q(0)*q(0) - q(1)*q(1) - q(2)*q(2) - q(3)*q(3)
         amplitude = m2s_ggfh( x, p(0,1,1),q,rand ) ! pp-> Hjj
         Amplitude = Amplitude*BHMU
         Amplitude = Amplitude*16D0*pi*XMGS(6)/
     -        ((m2s_dec - XM2S(6))**2 + XMGS(6)**2)
c     GF pp -> Hjj -> tau+ tau- jj
      CASE(GFHjj_tau)
         do mu = 0,3 !     sum the 2 lepton momenta up to Higgs momentum
            q(mu) = v(mu,1,1)+v(mu,2,1)
         enddo
         m2s_dec = q(0)*q(0) - q(1)*q(1) - q(2)*q(2) - q(3)*q(3)
         amplitude = m2s_ggfh( x, p(0,1,1),q,rand ) ! pp-> Hjj
         Amplitude = Amplitude*BHTAU
         Amplitude = Amplitude*16D0*pi*XMGS(6)/
     -        ((m2s_dec - XM2S(6))**2 + XMGS(6)**2)
c     GF pp -> Hjj -> b bbar jj
      CASE(GFHjj_bbar)
         do mu = 0,3 !     sum the 2 lepton momenta up to Higgs momentum
            q(mu) = v(mu,1,1)+v(mu,2,1)
         enddo
         m2s_dec = q(0)*q(0) - q(1)*q(1) - q(2)*q(2) - q(3)*q(3)
         amplitude = m2s_ggfh( x, p(0,1,1),q,rand ) ! pp-> Hjj
         Amplitude = amplitude*BHBB
         Amplitude = Amplitude*16D0*pi*XMGS(6)/
     -        ((m2s_dec - XM2S(6))**2 + XMGS(6)**2)
c     GF pp -> Hjj -> WWjj -> llnunu jj
      CASE(GFHjj_WW)
         do mu = 0,3 !     sum the 4 lepton momenta up to Higgs momentum
            q(mu) = v(mu,1,1)+v(mu,2,1)+v(mu,3,1)+v(mu,4,1)
         enddo
         call m2s_VVsum(v(0,1,1),3,4,2,1,amp)
         amp = amp / ( (xm2s(6) - mass2(q) )**2  +  xmgs(6)**2 ) 
         Amplitude = amp * m2s_ggfh( x, p(0,1,1),q,rand ) ! pp-> Hjj
c     GF pp -> Hjj -> Z Z jj 
      CASE(GFHjj_ZZ_ll)
         do mu = 0,3 !     sum the 4 lepton momenta up to Higgs momentum
            q(mu) = v(mu,1,1)+v(mu,2,1)+v(mu,3,1)+v(mu,4,1)
         enddo
         call m2s_VVsum(v(0,1,1),2,2,2,2,amp)
         amp = amp / ( (xm2s(6) - mass2(q) )**2  +  xmgs(6)**2 ) 
         Amplitude = amp * m2s_ggfh( x, p(0,1,1),q,rand ) ! pp-> Hjj
      CASE(GFHjj_ZZ_lnu)
         do mu = 0,3 !     sum the 4 lepton momenta up to Higgs momentum
            q(mu) = v(mu,1,1)+v(mu,2,1)+v(mu,3,1)+v(mu,4,1)
         enddo
         call m2s_VVsum(v(0,1,1),2,2,2,1,amp)
         amp = amp / ( (xm2s(6) - mass2(q) )**2  +  xmgs(6)**2 ) 
         Amplitude = amp * m2s_ggfh( x, p(0,1,1),q,rand ) ! pp-> Hjj
#endif

#ifdef WITH_QCDVV
**************************************
*
*     QCD : VVjj
*
**************************************
c     QCD pp -> W+ Z jj -> e+ ve mu+ mu- jj
      CASE(QCDWPZjj) 
         rn4(1) = RandomNumber()
         rn4(2) = rand(16) !RandomNumber()
         rn4(3) = RandomNumber()
         Amplitude = m2s_QCDwz2j( x, p(0,1,1), v(0,1,1) ,rn4,32,nlo)
c     QCD pp -> W+ Z jjj -> e+ ve mu+ mu- jjj
      CASE(QCDWPZjjjLO)
         do i=2,max_kin
            lokt(i) = .false.
         enddo
         lokt(1) = .true.
         rn2(1) = RandomNumber()
         rn2(2) = RandomNumber()
         call m2s_QCDwz3j( x, p(0,1,1), v(0,1,1) ,32,.false.,lokt,xuz,rn2,m2s)
         Amplitude = m2s(0)
c     QCD pp -> W- Z jj -> e- ve~ mu- mu+ jj
      CASE(QCDWMZjj) 
         rn4(1) = RandomNumber()
         rn4(2) = RandomNumber()
         rn4(3) = RandomNumber()
         Amplitude = m2s_QCDwz2j( x, p(0,1,1), v(0,1,1) ,rn4,42,nlo)
         rn=RandomNumber()         
c     QCD pp -> W- Z jjj -> e- ve~ mu- mu+ jjj
      CASE(QCDWMZjjjLO)
         do i=2,max_kin
            lokt(i) = .false.
         enddo
         lokt(1) = .true.
         rn = RandomNumber()
         call m2s_QCDwz3j( x, p(0,1,1), v(0,1,1) ,42,.false.,lokt,xuz,rn,m2s)
         Amplitude = m2s(0)
c     QCD pp -> W+ A jj -> e+ ve mu+ mu- jj
      CASE(QCDWPAjj) 
         rn4(1) = RandomNumber()
         rn4(2) = rand(13) !RandomNumber()
         rn4(3) = RandomNumber()
         Amplitude = m2s_QCDwz2j( x, p(0,1,1), v(0,1,1) ,rn4,31,nlo)
c     QCD pp -> W+ A jjj -> e+ ve mu+ mu- jjj
      CASE(QCDWPAjjjLO)
         do i=2,max_kin
            lokt(i) = .false.
         enddo
         lokt(1) = .true.
         rn2(1) = RandomNumber()
         rn2(2) = RandomNumber()
         call m2s_QCDwz3j( x, p(0,1,1), v(0,1,1) ,31,.false.,lokt,xuz,rn2,m2s)
         Amplitude = m2s(0)
c     QCD pp -> W- A jj -> e- ve~ mu- mu+ jj
      CASE(QCDWMAjj) 
         rn4(1) = RandomNumber()
         rn4(2) = RandomNumber()
         rn4(3) = RandomNumber()
         Amplitude = m2s_QCDwz2j( x, p(0,1,1), v(0,1,1) ,rn4,41,nlo)
         rn=RandomNumber()         
c     QCD pp -> W- A jjj -> e- ve~ mu- mu+ jjj
      CASE(QCDWMAjjjLO)
         do i=2,max_kin
            lokt(i) = .false.
         enddo
         lokt(1) = .true.
         rn = RandomNumber()
         call m2s_QCDwz3j( x, p(0,1,1), v(0,1,1) ,41,.false.,lokt,xuz,rn,m2s)
         Amplitude = m2s(0)
c     QCD pp -> Z A jj
      CASE(QCDZAjj_l) 
         rn4(1) = RandomNumber()
         rn4(2) = RandomNumber()
         rn4(3) = RandomNumber()
         rn4(4) = RandomNumber()
         Amplitude = m2s_QCDzz2j( x, p(0,1,1), v(0,1,1) ,rn4,20,nlo)        
      CASE(QCDZAjj_n) 
         rn4(1) = RandomNumber()
         rn4(2) = RandomNumber()
         rn4(3) = RandomNumber()
         rn4(4) = RandomNumber()
         Amplitude = m2s_QCDzz2j( x, p(0,1,1), v(0,1,1) ,rn4,23,nlo)        
      CASE(QCDZAJJJLO_l)
         do i=2,max_kin
            lokt(i) = .false.
         enddo
         lokt(1) = .true.
         rn4(1) = RandomNumber()
         rn4(2) = RandomNumber()
         rn4(3) = RandomNumber()
         call m2s_QCDzz3j(x,p,v,20,.false.,lokt,xuz,rn4,m2s)
         Amplitude=m2s(0)
      CASE(QCDZAJJJLO_n)
         do i=2,max_kin
            lokt(i) = .false.
         enddo
         lokt(1) = .true.
         rn4(1) = RandomNumber()
         rn4(2) = RandomNumber()
         rn4(3) = RandomNumber()
         call m2s_QCDzz3j(x,p,v,23,.false.,lokt,xuz,rn4,m2s)
         Amplitude=m2s(0)
      CASE(QCDAAjj) 
         rn4(1) = RandomNumber()
         rn4(2) = RandomNumber()
         rn4(3) = RandomNumber()
         rn4(4) = RandomNumber()
         Amplitude = m2s_QCDzz2j( x, p(0,1,1), v(0,1,1) ,rn4,24,nlo)     
      CASE(QCDZZjj_ll) 
         rn4(1) = RandomNumber()
         rn4(2) = RandomNumber()
         rn4(3) = RandomNumber()
         rn4(4) = RandomNumber()
         Amplitude = m2s_QCDzz2j( x, p(0,1,1), v(0,1,1) ,rn4,22,nlo)        
      CASE(QCDZZjj_lnu) 
         rn4(1) = RandomNumber()
         rn4(2) = RandomNumber()
         rn4(3) = RandomNumber()
         rn4(4) = RandomNumber()
         Amplitude = m2s_QCDzz2j( x, p(0,1,1), v(0,1,1) ,rn4,21,nlo)         
      CASE(QCDZZJJJLO_ll)
         do i=2,max_kin
            lokt(i) = .false.
         enddo
         lokt(1) = .true.
         rn4(1) = RandomNumber()
         rn4(2) = RandomNumber()
         rn4(3) = RandomNumber()
         call m2s_QCDzz3j(x,p,v,22,.false.,lokt,xuz,rn4,m2s)
         Amplitude=m2s(0)
      CASE(QCDAAjjjLO)
         do i=2,max_kin
            lokt(i) = .false.
         enddo
         lokt(1) = .true.
         rn4(1) = RandomNumber()
         rn4(2) = RandomNumber()
         rn4(3) = RandomNumber()
         call m2s_QCDzz3j(x,p,v,24,.false.,lokt,xuz,rn4,m2s)
         Amplitude=m2s(0)
      CASE(QCDZZJJJLO_lnu)
         do i=2,max_kin
            lokt(i) = .false.
         enddo
         lokt(1) = .true.
         rn4(1) = RandomNumber()
         rn4(2) = RandomNumber()
         rn4(3) = RandomNumber()
         call m2s_QCDzz3j(x,p,v,21,.false.,lokt,xuz,rn4,m2s)
         Amplitude=m2s(0)
c     QCD pp -> W+ W+ jj
      CASE(QCDWPWPjj) 
         rn2(1) = RandomNumber()
         rn2(2) = RandomNumber()
         Amplitude = m2s_QCDww2j( x, p(0,1,1), v(0,1,1) ,rn2,33,nlo)
c     QCD pp -> W+ W+ jjj
      CASE(QCDWPWPjjjLO)
         do i=2,max_kin
            lokt(i) = .false.
         enddo
         lokt(1) = .true.
         rn = RandomNumber()
         call m2s_QCDww3j( x, p(0,1,1), v(0,1,1) ,33,.false.,lokt,xuz,rn,m2s)
         Amplitude = m2s(0)
c     QCD pp -> W- W- jj
      CASE(QCDWMWMjj) 
         rn2(1) = RandomNumber()
         rn2(2) = RandomNumber()
         Amplitude = m2s_QCDww2j( x, p(0,1,1), v(0,1,1) ,rn2,44,nlo)
c     QCD pp -> W- W- jjj
      CASE(QCDWMWMjjjLO)
         do i=2,max_kin
            lokt(i) = .false.
         enddo
         lokt(1) = .true.
         rn = RandomNumber()
         call m2s_QCDww3j( x, p(0,1,1), v(0,1,1) ,44,.false.,lokt,xuz,rn,m2s)
         Amplitude = m2s(0)
#endif

#ifdef WITH_QCDV
**************************************
*
*     QCD : Zjj
*
**************************************
c     QCD pp -> Z jj -> ll jj
      CASE(QCDZjj_l)
         rn4(1) = RandomNumber()
         rn4(2) = RandomNumber()
         rn4(3) = RandomNumber()
         rn4(4) = RandomNumber()
         Amplitude = m2s_QCDzz2j( x, p(0,1,1), v(0,1,1) ,rn4,2,nlo)
c     QCD pp -> Z jj -> nunu jj
      CASE(QCDZjj_nu)
         rn4(1) = RandomNumber()
         rn4(2) = RandomNumber()
         rn4(3) = RandomNumber()
         rn4(4) = RandomNumber()
         Amplitude = m2s_QCDzz2j( x, p(0,1,1), v(0,1,1) ,rn4,-2,nlo)
      CASE(QCDZJJJLO_l)
         do i=2,max_kin
            lokt(i) = .false.
         enddo
         lokt(1) = .true.
         rn4(1) = RandomNumber()
         rn4(2) = RandomNumber()
         rn4(3) = RandomNumber()
         call m2s_QCDzz3j(x,p,v,2,.false.,lokt,xuz,rn4,m2s)
         Amplitude=m2s(0)
      CASE(QCDZJJJLO_nu)
         do i=2,max_kin
            lokt(i) = .false.
         enddo
         lokt(1) = .true.
         rn4(1) = RandomNumber()
         rn4(2) = RandomNumber()
         rn4(3) = RandomNumber()
         call m2s_QCDzz3j(x,p,v,-2,.false.,lokt,xuz,rn4,m2s)
         Amplitude=m2s(0)
c     QCD pp -> W+ jj 
      CASE(QCDWPjj)
         rn4(1) = RandomNumber()
         rn4(2) = RandomNumber()
         rn4(3) = RandomNumber()
         Amplitude = m2s_QCDwz2j( x, p(0,1,1), v(0,1,1) ,rn4,3,nlo)
c     QCD pp -> W- jj 
      CASE(QCDWMjj)
         rn4(1) = RandomNumber()
         rn4(2) = RandomNumber()
         rn4(3) = RandomNumber()
         Amplitude = m2s_QCDwz2j( x, p(0,1,1), v(0,1,1) ,rn4,4,nlo)
c     QCD pp -> W+ jjj 
      CASE(QCDWPjjjLO)
         do i=2,max_kin
            lokt(i) = .false.
         enddo
         lokt(1) = .true.
         rn2(1) = RandomNumber()
         rn2(2) = RandomNumber()
         call m2s_QCDwz3j( x, p(0,1,1), v(0,1,1) ,3,.false.,lokt,xuz,rn2,m2s)
         Amplitude = m2s(0)
c     QCD pp -> W- jjj 
      CASE(QCDWMjjjLO)
         do i=2,max_kin
            lokt(i) = .false.
         enddo
         lokt(1) = .true.
         rn2(1) = RandomNumber()
         rn2(2) = RandomNumber()
         call m2s_QCDwz3j( x, p(0,1,1), v(0,1,1) ,4,.false.,lokt,xuz,rn2,m2s)
         Amplitude = m2s(0)
#endif

      CASE DEFAULT
         print *,"Error: Amplitude for process not implemented, ID = ",procID
         print *,"Check process-ID and compiler flags"
         STOP
      END SELECT
      

      RETURN
      END


#ifdef WITH_NLO
c*************************************************************************
      SUBROUTINE RE_Amplitude(rand, p, x, v, lokt, xuz, m2s, ps_number)
c*************************************************************************
c     Amplitude for real emissions.
c*************************************************************************
      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/process.inc"

      double precision x(nx) 
      double precision rand(max_PS_dim)
      double precision p(0:3,max_p,max_kin)
      double precision v(0:3,max_v,max_kin)
      double precision rn, rn3(3), xuz(2,2:max_kin), m2s(0:max_kin), m2s_dec
      logical lokt(max_kin)
      integer L, ps_number

c     masses and widths
      double precision XM2s(6),XMGs(6)
      COMMON /BKOPOUshort/ XM2s,XMGs
      double precision BWNE,BWUD,BWTB,BZNN,BZEE,BZUU,BZDD,BZTT,
     &                BHWW,BHZZ,BHGG,BHTT,BHBB,BHCC,BHTAU,BHMU,
     &                BHGAM, BHGAMZ, XGW, XGZ, XGH, GAMT
      COMMON /BRANCH/ BWNE,BWUD,BWTB,BZNN,BZEE,BZUU,BZDD,BZTT,
     &                BHWW,BHZZ,BHGG,BHTT,BHBB,BHCC,BHTAU,BHMU,
     &                BHGAM, BHGAMZ, XGW, XGZ, XGH, GAMT

      double precision RandomNumber
      external RandomNumber
   

      SELECT CASE(process)
#ifdef WITH_VBF
**************************************
*
*     VBF : Hjjj, WWjjj, ZZjjj, WZjjj
*
**************************************
c     pp -> Hjjj
      CASE(Hjj)
         rn = RandomNumber()
         call m2s_wbfzhg( 6, 1, lokt,
     1        x, p(0,1,1), v(0,1,1), rn, xuz, m2s)
c     pp -> Hjjj -> gamma gamma jjj
      CASE(Hjj_AA)
         rn = RandomNumber()
         call m2s_wbfzhg( 6, 1, lokt,
     1        x, p(0,1,1), v(0,1,1), rn, xuz, m2s)
         do L = 0,3
            m2s(L) = m2s(L)*BHGAM
         enddo
c     pp -> spin-2 jjj -> gamma gamma jjj
      CASE(AAjj)
         rn = RandomNumber()
         call m2s_wbfzhg( 11, 1, lokt, 
     1                    x, p(0,1,1), v(0,1,1), rn, xuz, m2s) 
c     pp -> spin-2 jjj -> WW jjj
      CASE(Sp2jj_WW)
         rn = RandomNumber()
         call m2s_wbfzhg( 734, 1, lokt, 
     1                    x, p(0,1,1), v(0,1,1), rn, xuz, m2s)
c     pp -> spin-2 jjj -> ZZ jjj
      CASE(Sp2jj_ZZ_ll)
         rn = RandomNumber()
         call m2s_wbfzhg( 722, 1, lokt, 
     1                    x, p(0,1,1), v(0,1,1), rn, xuz, m2s)
      CASE(Sp2jj_ZZ_lnu)
         rn = RandomNumber()
         call m2s_wbfzhg( 721, 1, lokt, 
     1                    x, p(0,1,1), v(0,1,1), rn, xuz, m2s)
c     pp -> Hjjj -> mu+mu- jjj
      CASE(Hjj_mu)
         rn = RandomNumber()
         call m2s_wbfzhg( 6, 1, lokt,
     1        x, p(0,1,1), v(0,1,1), rn, xuz, m2s)
         do L = 0,3
            m2s(L) = m2s(L)*BHMU
         enddo
c     pp -> Hjjj -> tau+tau- jjj
      CASE(Hjj_tau)
         rn = RandomNumber()
         call m2s_wbfzhg( 6, 1, lokt,
     1        x, p(0,1,1), v(0,1,1), rn, xuz, m2s)
         do L = 0,3
            m2s(L) = m2s(L)*BHTAU
         enddo
c     pp -> Hjjj -> b bbar jjj
      CASE(Hjj_bbar)
         rn = RandomNumber()
         call m2s_wbfzhg( 6, 1, lokt,
     1        x, p(0,1,1), v(0,1,1), rn, xuz, m2s)
         do L = 0,3
            m2s(L) = m2s(L)*BHBB
         enddo
c     pp -> Hjjj -> W W jjj
      CASE(Hjj_WW)
         rn = RandomNumber()
         call m2s_wbfzhg(6, 1, lokt,
     1        x, p(0,1,1), v(0,1,1), rn, xuz, m2s)
         call m2s_VVsum(v(0,1,1),3,4,1,2,m2s_dec)
         m2s_dec = m2s_dec/(16.0d0 * pi * xmgs(6))
         do L = 0,3
            m2s(L) = m2s(L)*m2s_dec
         enddo
      CASE(Hjj_WPhadWM, Hjj_WPWMhad)
         rn = RandomNumber()
         call m2s_wbfzhg(6, 1, lokt,
     1        x, p(0,1,1), v(0,1,1), rn, xuz, m2s)
         if (process.eq.Hjj_WPhadWM) then
            call m2s_VVsum_had(v(0,1,1),3,4,3,2,m2s_dec)
         else
            call m2s_VVsum_had(v(0,1,1),3,4,1,4,m2s_dec)
         endif
         m2s_dec = m2s_dec/(16.0d0 * pi * xmgs(6))
         do L = 0,3
            m2s(L) = m2s(L)*m2s_dec
         enddo
c     pp -> Hjjj -> Z Z jjj
      CASE(Hjj_ZZ_ll)
         rn = RandomNumber()
         call m2s_wbfzhg(6, 1, lokt,
     1        x, p(0,1,1), v(0,1,1), rn, xuz, m2s)
         call m2s_VVsum(v(0,1,1),2,2,2,2,m2s_dec)
         m2s_dec = m2s_dec/(16.0d0 * pi * xmgs(6))
         do L = 0,3
            m2s(L) = m2s(L)*m2s_dec
         enddo
      CASE(Hjj_ZZhad)
         rn = RandomNumber()
         call m2s_wbfzhg(6, 1, lokt,
     1        x, p(0,1,1), v(0,1,1), rn, xuz, m2s)
         call m2s_VVsum_had(v(0,1,1),2,2,5,2,m2s_dec)
         m2s_dec = m2s_dec/(16.0d0 * pi * xmgs(6))
         do L = 0,3
            m2s(L) = m2s(L)*m2s_dec
         enddo
      CASE(Hjj_ZZ_lnu)
         rn = RandomNumber()
         call m2s_wbfzhg(6, 1, lokt,
     1        x, p(0,1,1), v(0,1,1), rn, xuz, m2s)
         call m2s_VVsum(v(0,1,1),2,2,2,1,m2s_dec)
         m2s_dec = m2s_dec/(16.0d0 * pi * xmgs(6))
         do L = 0,3
            m2s(L) = m2s(L)*m2s_dec
         enddo
**************************************
*
*     VBF : HHjjj
*
**************************************
c     pp -> HHjjj
      CASE(HHjj)
         rn = RandomNumber()
         call m2s_wbfzhg(66, 1, lokt,x, p(0,1,1), v(0,1,1), rn, xuz, m2s)
      CASE(HHjj_bbtau)
         rn = RandomNumber()
         call m2s_wbfzhg(66, 1, lokt,x, p(0,1,1), v(0,1,1), rn, xuz, m2s)
         do L = 0,3
            m2s(L) = m2s(L)*BHBB*BHTAU
         enddo
      CASE(HHjj_bbAA)
         rn = RandomNumber()
         call m2s_wbfzhg(66, 1, lokt,x, p(0,1,1), v(0,1,1), rn, xuz, m2s)
         do L = 0,3
            m2s(L) = m2s(L)*BHBB*BHGAM
         enddo
         
#endif
         
#ifdef WITH_HJJJ
**************************************
*
*     VBF : Hjjjj
*
**************************************
c     pp -> Hjjjj
      CASE(Hjjj)
         call m2s_wbfh4j(6, 1, lokt, x, p(0,1,1), v(0,1,1), xuz, m2s)
c     pp -> Hjjjj -> gamma gamma jjjj
      CASE(Hjjj_AA)
         call m2s_wbfh4j(6, 1, lokt, x, p(0,1,1), v(0,1,1), xuz, m2s)
         do L = 0,n_kin
            m2s(L) = m2s(L)*BHGAM
         enddo
c     pp -> Hjjjj -> mu+mu- jjjj
      CASE(Hjjj_mu)
         call m2s_wbfh4j(6, 1, lokt, x, p(0,1,1), v(0,1,1), xuz, m2s)
         do L = 0,n_kin
            m2s(L) = m2s(L)*BHMU
         enddo
c     pp -> Hjjjj -> tau+tau- jjjj
      CASE(Hjjj_tau)
         call m2s_wbfh4j(6, 1, lokt, x, p(0,1,1), v(0,1,1), xuz, m2s)
         do L = 0,n_kin
            m2s(L) = m2s(L)*BHTAU
         enddo
c     pp -> Hjjjj -> b bbar jjjj
      CASE(Hjjj_bbar)
         call m2s_wbfh4j(6, 1, lokt, x, p(0,1,1), v(0,1,1), xuz, m2s)
         do L = 0,n_kin
            m2s(L) = m2s(L)*BHBB
         enddo
c     pp -> Hjjjj -> W W jjjj
      CASE(Hjjj_WW)
         call m2s_wbfh4j(6, 1, lokt, x, p(0,1,1), v(0,1,1), xuz, m2s)
         call m2s_VVsum(v(0,1,1),3,4,2,1,m2s_dec)
         m2s_dec = m2s_dec/(16.0d0 * pi * xmgs(6))
         do L = 0,n_kin
            m2s(L) = m2s(L)*m2s_dec
         enddo
c     pp -> Hjjjj -> Z Z jjjj
      CASE(Hjjj_ZZ_ll)
         call m2s_wbfh4j(6, 1, lokt, x, p(0,1,1), v(0,1,1), xuz, m2s)
         call m2s_VVsum(v(0,1,1),2,2,2,2,m2s_dec)
         m2s_dec = m2s_dec/(16.0d0 * pi * xmgs(6))
         do L = 0,n_kin
            m2s(L) = m2s(L)*m2s_dec
         enddo

      CASE(Hjjj_ZZ_lnu)
         call m2s_wbfh4j(6, 1, lokt, x, p(0,1,1), v(0,1,1), xuz, m2s)
         call m2s_VVsum(v(0,1,1),2,2,2,1,m2s_dec)
         m2s_dec = m2s_dec/(16.0d0 * pi * xmgs(6))
         do L = 0,n_kin
            m2s(L) = m2s(L)*m2s_dec
         enddo
#endif 

#ifdef WITH_VBF
**************************************
*
*     VBF : Wjjj, Zjjj
*
**************************************
c     pp -> Ajjj 
      CASE(Ajj)
         rn = RandomNumber()
         call m2s_wbfzhg( 1, 1, lokt,
     1        x, p(0,1,1), v(0,1,1), rn, xuz, m2s)
c     pp -> Zjjj 
      CASE(Zjj_l, Zjj_nu) 
         rn = RandomNumber()
         call m2s_wbfzhg( 2, 1, lokt, 
     1        x, p(0,1,1), v(0,1,1), rn, xuz, m2s)
c     pp -> W+ jjj
      CASE(WPjj) 
         rn = RandomNumber()
         call m2s_wbfw3j( 3, 1, lokt, 
     1        x, p(0,1,1), v(0,1,1), rn, xuz, m2s)
c     pp -> W- jjj
      CASE(WMjj)
         rn = RandomNumber()
         call m2s_wbfw3j( 4, 1, lokt,
     1        x, p(0,1,1), v(0,1,1), rn, xuz, m2s)


**************************************
*
*     VBF : WWjjj, ZZjjj, WZjjj, WAjjj
*
**************************************
c     pp -> WW jjj
      CASE(WPWMjj, WPhadWMjj, WPWMhadjj)
         rn = RandomNumber()
         call m2s_wbfzhg( 34, 1, lokt, 
     1                    x, p(0,1,1), v(0,1,1), rn, xuz, m2s) 
c     pp -> ZZ jjj
      CASE(ZZjj_ll, ZZhadjj)
         rn = RandomNumber()
         call m2s_wbfzhg( 22, 1, lokt, 
     1                    x, p(0,1,1), v(0,1,1), rn, xuz, m2s) 
      CASE(ZZjj_lnu)
         rn = RandomNumber()
         call m2s_wbfzhg( 21, 1, lokt, 
     1                    x, p(0,1,1), v(0,1,1), rn, xuz, m2s) 
c     pp -> W+Z jjj
      CASE(WPZjj, WPhadZjj, WPZhadjj) 
         rn = RandomNumber()
         call m2s_wbfw3j( 32, 1, lokt, 
     1        x, p(0,1,1), v(0,1,1), rn, xuz, m2s)
c     pp -> W-Z jjj
      CASE(WMZjj, WMhadZjj, WMZhadjj) 
         rn = RandomNumber()
         call m2s_wbfw3j( 42, 1, lokt, 
     1        x, p(0,1,1), v(0,1,1), rn, xuz, m2s)
c     pp -> W+A jjj
      CASE(WPAjj) 
         rn = RandomNumber()
         call m2s_wbfw3j( 31, 1, lokt, 
     1        x, p(0,1,1), v(0,1,1), rn, xuz, m2s)
c     pp -> W-A jjj
      CASE(WMAjj) 
         rn = RandomNumber()
         call m2s_wbfw3j( 41, 1, lokt, 
     1        x, p(0,1,1), v(0,1,1), rn, xuz, m2s)
      CASE(WPWPjj, WMWMjj, WPhadWPjj, WMhadWMjj)
         rn = RandomNumber()
         call m2s_wbfwppg( 33,1, lokt, 
     1                    x, p(0,1,1), v(0,1,1), rn, xuz, m2s) 
c     pp -> ZA jjj
      CASE(ZAjj)
         rn = RandomNumber()
         call m2s_wbfzhg( 211, 1, lokt, 
     1          x, p(0,1,1), v(0,1,1), rn, xuz, m2s) 
c     pp -> ZnA jjj
      CASE(ZAjj_n)
         rn = RandomNumber()
         call m2s_wbfzhg( 212, 1, lokt, 
     1         x, p(0,1,1), v(0,1,1), rn, xuz, m2s) 

**************************************
*
*     VBF : HAjjj
*
**************************************
c     pp -> HAjjj
      CASE(HAjj)
         rn = RandomNumber()
         call m2s_wbfzhg( 61, 1, lokt,
     1        x, p(0,1,1), v(0,1,1), rn, xuz, m2s)
c     pp -> HAjjj -> A gamma gamma jjj
      CASE(HAjj_AA)
         rn = RandomNumber()
         call m2s_wbfzhg( 61, 1, lokt,
     1        x, p(0,1,1), v(0,1,1), rn, xuz, m2s)
         do L = 0,3
            m2s(L) = m2s(L)*BHGAM
         enddo
c     pp -> HAjjj -> A mu+mu- jjj
      CASE(HAjj_mu)
         rn = RandomNumber()
         call m2s_wbfzhg( 61, 1, lokt,
     1        x, p(0,1,1), v(0,1,1), rn, xuz, m2s)
         do L = 0,3
            m2s(L) = m2s(L)*BHMU
         enddo
c     pp -> HAjjj -> A tau+tau- jjj
      CASE(HAjj_tau)
         rn = RandomNumber()
         call m2s_wbfzhg( 61, 1, lokt,
     1        x, p(0,1,1), v(0,1,1), rn, xuz, m2s)
         do L = 0,3
            m2s(L) = m2s(L)*BHTAU
         enddo
c     pp -> HAjjj -> A b bbar jjj
      CASE(HAjj_bbar)
         rn = RandomNumber()
         call m2s_wbfzhg( 61, 1, lokt,
     1        x, p(0,1,1), v(0,1,1), rn, xuz, m2s)
         do L = 0,3
            m2s(L) = m2s(L)*BHBB
         enddo
c     pp -> HAjjj -> A W W jjj
      CASE(HAjj_WW)
         rn = RandomNumber()
         call m2s_wbfzhg(61, 1, lokt,
     1        x, p(0,1,1), v(0,1,1), rn, xuz, m2s)
         call m2s_VVsum(v(0,1,1),3,4,2,1,m2s_dec)
         m2s_dec = m2s_dec/(16.0d0 * pi * xmgs(6))
         do L = 0,3
            m2s(L) = m2s(L)*m2s_dec
         enddo
c     pp -> HAjjj -> A Z Z jjj
      CASE(HAjj_ZZ_ll)
         rn = RandomNumber()
         call m2s_wbfzhg(61, 1, lokt,
     1        x, p(0,1,1), v(0,1,1), rn, xuz, m2s)
         call m2s_VVsum(v(0,1,1),2,2,2,2,m2s_dec)
         m2s_dec = m2s_dec/(16.0d0 * pi * xmgs(6))
         do L = 0,3
            m2s(L) = m2s(L)*m2s_dec
         enddo
      CASE(HAjj_ZZ_lnu)
         rn = RandomNumber()
         call m2s_wbfzhg(61, 1, lokt,
     1        x, p(0,1,1), v(0,1,1), rn, xuz, m2s)
         call m2s_VVsum(v(0,1,1),2,2,2,1,m2s_dec)
         m2s_dec = m2s_dec/(16.0d0 * pi * xmgs(6))
         do L = 0,3
            m2s(L) = m2s(L)*m2s_dec
         enddo
#endif


#ifdef WITH_DIBOSON
**************************************
*
*     DIBOSON : WWj, WZj, ZZj, WAj, ZAj, AAj
*
**************************************         
c pp -> W+ W- (j) -> ve e+ e- ve~ (j)
      CASE(WW,WPhadWMlep,WPlepWMhad)
         rn = RandomNumber()
         call m2s_qqVVj( 34, 1, lokt, x, p, v, rn, xuz(1,2), m2s, ps_number)
      CASE(WPZ,WPhadZ,WPZhad)
         rn = RandomNumber()
         call m2s_qqWZj( 32, 1, lokt, x, p, v, rn, xuz(1,2), m2s)
      CASE(WMZ,WMhadZ,WMZhad)
         rn = RandomNumber()
         call m2s_qqWZj( 42, 1, lokt, x, p, v, rn, xuz(1,2), m2s)
c pp -> Z Z (j) -> e- e+ mu- mu+ (j)
      CASE(ZZ,ZZhad)
         rn = RandomNumber()
         call m2s_qqZZj( 22, 1, lokt, x, p, v, rn, xuz(1,2), m2s, ps_number)
      CASE(WPA)
         rn = RandomNumber()
         call m2s_qqWAj( 31, 1, lokt, x, p, v, rn, xuz(1,2), m2s)
      CASE(WMA)
         rn = RandomNumber()
         call m2s_qqWAj( 41, 1, lokt, x, p, v, rn, xuz(1,2), m2s)
c pp -> Z A (j) -> e- e+ A (j)
      CASE(ZA)
         rn = RandomNumber()
         call m2s_qqZAj( 21, 1, lokt, x, p, v, rn, xuz(1,2), m2s, ps_number)
c pp -> A A (j) 
      CASE(AA)
         rn = RandomNumber()
         call m2s_qqAAj( 11, 1, lokt, x, p, v, rn, xuz(1,2), m2s, ps_number)


c     pp -> WPHj
      CASE(WPH)
         rn = RandomNumber()
         call m2s_qqWHj( 36, 1, lokt, x, p, v, rn, xuz(1,2), m2s)
c     pp -> WPHj -> gamma gamma j
      CASE(WPH_AA)
         rn = RandomNumber()
         call m2s_qqWHj( 36, 1, lokt, x, p, v, rn, xuz(1,2), m2s)
         do L = 0,3
            m2s(L) = m2s(L)*BHGAM
         enddo
c     pp -> WPHj -> mu+mu- j
      CASE(WPH_mu)
         rn = RandomNumber()
         call m2s_qqWHj( 36, 1, lokt, x, p, v, rn, xuz(1,2), m2s)
         do L = 0,3
            m2s(L) = m2s(L)*BHMU
         enddo
c     pp -> WPHj -> tau+tau- j
      CASE(WPH_tau)
         rn = RandomNumber()
         call m2s_qqWHj( 36, 1, lokt, x, p, v, rn, xuz(1,2), m2s)
         do L = 0,3
            m2s(L) = m2s(L)*BHTAU
         enddo
c     pp -> WPHj -> b bbar j
      CASE(WPH_bbar)
         rn = RandomNumber()
         call m2s_qqWHj( 36, 1, lokt, x, p, v, rn, xuz(1,2), m2s)
         do L = 0,3
            m2s(L) = m2s(L)*BHBB
         enddo
c     pp -> WPHj -> W W j
      CASE(WPH_WW)
         rn = RandomNumber()
         call m2s_qqWHj( 36, 1, lokt, x, p, v, rn, xuz(1,2), m2s)
         call m2s_VVsum(v(0,3,1),3,4,2,1,m2s_dec)
         m2s_dec = m2s_dec/(16.0d0 * pi * xmgs(6))
         do L = 0,3
            m2s(L) = m2s(L)*m2s_dec
         enddo
c     pp -> WPHj -> Z Z j
      CASE(WPH_ZZ_ll)
         rn = RandomNumber()
         call m2s_qqWHj( 36, 1, lokt, x, p, v, rn, xuz(1,2), m2s)
         call m2s_VVsum(v(0,3,1),2,2,2,2,m2s_dec)
         m2s_dec = m2s_dec/(16.0d0 * pi * xmgs(6))
         do L = 0,3
            m2s(L) = m2s(L)*m2s_dec
         enddo
      CASE(WPH_ZZ_lnu)
         rn = RandomNumber()
         call m2s_qqWHj( 36, 1, lokt, x, p, v, rn, xuz(1,2), m2s)
         call m2s_VVsum(v(0,3,1),2,2,2,1,m2s_dec)
         m2s_dec = m2s_dec/(16.0d0 * pi * xmgs(6))
         do L = 0,3
            m2s(L) = m2s(L)*m2s_dec
         enddo


c     pp -> WMHj
      CASE(WMH)
         rn = RandomNumber()
         call m2s_qqWHj( 46, 1, lokt, x, p, v, rn, xuz(1,2), m2s)
c     pp -> WMHj -> gamma gamma j
      CASE(WMH_AA)
         rn = RandomNumber()
         call m2s_qqWHj( 46, 1, lokt, x, p, v, rn, xuz(1,2), m2s)
         do L = 0,3
            m2s(L) = m2s(L)*BHGAM
         enddo
c     pp -> WMHj -> mu+mu- j
      CASE(WMH_mu)
         rn = RandomNumber()
         call m2s_qqWHj( 46, 1, lokt, x, p, v, rn, xuz(1,2), m2s)
         do L = 0,3
            m2s(L) = m2s(L)*BHMU
         enddo
c     pp -> WMHj -> tau+tau- j
      CASE(WMH_tau)
         rn = RandomNumber()
         call m2s_qqWHj( 46, 1, lokt, x, p, v, rn, xuz(1,2), m2s)
         do L = 0,3
            m2s(L) = m2s(L)*BHTAU
         enddo
c     pp -> WMHj -> b bbar j
      CASE(WMH_bbar)
         rn = RandomNumber()
         call m2s_qqWHj( 46, 1, lokt, x, p, v, rn, xuz(1,2), m2s)
         do L = 0,3
            m2s(L) = m2s(L)*BHBB
         enddo
c     pp -> WMHj -> W W j
      CASE(WMH_WW)
         rn = RandomNumber()
         call m2s_qqWHj( 46, 1, lokt, x, p, v, rn, xuz(1,2), m2s)
         call m2s_VVsum(v(0,3,1),3,4,2,1,m2s_dec)
         m2s_dec = m2s_dec/(16.0d0 * pi * xmgs(6))
         do L = 0,3
            m2s(L) = m2s(L)*m2s_dec
         enddo
c     pp -> WMHj -> Z Z j
      CASE(WMH_ZZ_ll)
         rn = RandomNumber()
         call m2s_qqWHj( 46, 1, lokt, x, p, v, rn, xuz(1,2), m2s)
         call m2s_VVsum(v(0,3,1),2,2,2,2,m2s_dec)
         m2s_dec = m2s_dec/(16.0d0 * pi * xmgs(6))
         do L = 0,3
            m2s(L) = m2s(L)*m2s_dec
         enddo
      CASE(WMH_ZZ_lnu)
         rn = RandomNumber()
         call m2s_qqWHj( 46, 1, lokt, x, p, v, rn, xuz(1,2), m2s)
         call m2s_VVsum(v(0,3,1),2,2,2,1,m2s_dec)
         m2s_dec = m2s_dec/(16.0d0 * pi * xmgs(6))
         do L = 0,3
            m2s(L) = m2s(L)*m2s_dec
         enddo

c pp -> ve e      
      CASE(WP_only)
         rn = RandomNumber()
         call m2s_qqWHj( 3, 1, lokt, x, p, v, rn, xuz(1,2), m2s)
      CASE(WM_only)
         rn = RandomNumber()
         call m2s_qqWHj( 4, 1, lokt, x, p, v, rn, xuz(1,2), m2s)

#endif


#ifdef WITH_DIBOSONJET
**************************************
*
*     DIBOSON PLUS JET : WWjj,WZjj,ZZjj,WAjj
*
**************************************         
      CASE(WWJ,WPHADWMJ,WPWMHADJ)
         call m2s_qqVVjj( 34, 1, lokt, x, p, v, rn, xuz, m2s)
c pp -> W+-A (jj) -> ve e+- a (jj)
      CASE(WMAJ)
         call m2s_qqwgagluj(xuz,1,.false.,lokt,x,p,v,m2s)
      CASE(WPAJ)
         call m2s_qqwpgagluj(xuz,1,.true.,lokt,x,p,v,m2s)
c pp -> W+-Z (jj) -> ve e+- mu+ mu- (jj)
      CASE(WMZJ, WMhadZJ, WMZhadJ)
         call m2s_qqwzgluj(xuz,1,.false.,lokt,x,p,v,m2s)
      CASE(WPZJ, WPhadZJ, WPZhadJ)
         call m2s_qqwpzgluj(xuz,1,.true.,lokt,x,p,v,m2s)       
c pp -> ZZ (jj) -> e+ e- mu+ mu- (jj)
      CASE(ZZJ)
         call m2s_qqZZjj( 22, 1, lokt, x, p, v, rn, xuz, m2s)
c     pp -> W HJ
      CASE(WPHj)
         rn = RandomNumber()
         call m2s_qqWHjj( 36, 1, lokt, x, p, v, rn, xuz(1,2), 
     &                           narrowwidth, 0, m2s)
c     pp -> W HJ -> W gamma gamma
      CASE(WPHj_AA)
         rn = RandomNumber()
         call m2s_qqWHjj( 36, 1, lokt, x, p, v, rn, xuz(1,2), 
     &                           narrowwidth, 0, m2s)
         m2s = m2s*BHGAM
c     pp -> W HJ -> W mu+mu-
      CASE(WPHj_mu)
         rn = RandomNumber()
         call m2s_qqWHjj( 36, 1, lokt, x, p, v, rn, xuz(1,2), 
     &                           narrowwidth, 0, m2s)
         m2s = m2s*BHMU
c     pp -> W HJ -> W tau+tau-
      CASE(WPHj_tau)
         rn = RandomNumber()
         call m2s_qqWHjj( 36, 1, lokt, x, p, v, rn, xuz(1,2), 
     &                           narrowwidth, 0, m2s)
         m2s = m2s*BHTAU
c     pp -> W HJ -> W b bbar
      CASE(WPHj_bbar)
         rn = RandomNumber()
         call m2s_qqWHjj( 36, 1, lokt, x, p, v, rn, xuz(1,2), 
     &                           narrowwidth, 0, m2s)
         m2s = m2s*BHBB
c     pp -> W HJ -> W W W
      CASE(WPHj_WW)
         rn = RandomNumber()
         call m2s_qqWHjj( 36, 1, lokt, x, p, v, rn, xuz(1,2), 
     &                           narrowwidth, 0, m2s)
         call m2s_VVsum(v(0,3,1),3,4,2,1,m2s_dec) !use v3..v6 for H->VV
         m2s_dec = m2s_dec/(16.0d0 * pi * xmgs(6))
         m2s = m2s*m2s_dec ! implicit loop over array m2s 
c     pp -> W HJ -> W Z Z 
      CASE(WPHj_ZZ_ll)
         rn = RandomNumber()
         call m2s_qqWHjj( 36, 1, lokt, x, p, v, rn, xuz(1,2), 
     &                           narrowwidth, 0, m2s)
         call m2s_VVsum(v(0,3,1),2,2,2,2,m2s_dec)
         m2s_dec = m2s_dec/(16.0d0 * pi * xmgs(6))
         m2s = m2s*m2s_dec
      CASE(WPHj_ZZ_lnu)
         rn = RandomNumber()
         call m2s_qqWHjj( 36, 1, lokt, x, p, v, rn, xuz(1,2), 
     &                           narrowwidth, 0, m2s)
         call m2s_VVsum(v(0,3,1),2,2,2,1,m2s_dec)
         m2s_dec = m2s_dec/(16.0d0 * pi * xmgs(6))
         m2s = m2s*m2s_dec  


c     pp -> W HJ
      CASE(WMHj)
         rn = RandomNumber()
         call m2s_qqWHjj( 46, 1, lokt, x, p, v, rn, xuz(1,2), 
     &                           narrowwidth, 0, m2s)
c     pp -> W HJ -> W gamma gamma
      CASE(WMHj_AA)
         rn = RandomNumber()
         call m2s_qqWHjj( 46, 1, lokt, x, p, v, rn, xuz(1,2), 
     &                           narrowwidth, 0, m2s)
         m2s = m2s*BHGAM
c     pp -> W HJ -> W mu+mu-
      CASE(WMHj_mu)
         rn = RandomNumber()
         call m2s_qqWHjj( 46, 1, lokt, x, p, v, rn, xuz(1,2), 
     &                           narrowwidth, 0, m2s)
         m2s = m2s*BHMU
c     pp -> W HJ -> W tau+tau-
      CASE(WMHj_tau)
         rn = RandomNumber()
         call m2s_qqWHjj( 46, 1, lokt, x, p, v, rn, xuz(1,2), 
     &                           narrowwidth, 0, m2s)
         m2s = m2s*BHTAU
c     pp -> W HJ -> W b bbar
      CASE(WMHj_bbar)
         rn = RandomNumber()
         call m2s_qqWHjj( 46, 1, lokt, x, p, v, rn, xuz(1,2), 
     &                           narrowwidth, 0, m2s)
         m2s = m2s*BHBB
c     pp -> W HJ -> W W W
      CASE(WMHj_WW)
         rn = RandomNumber()
         call m2s_qqWHjj( 46, 1, lokt, x, p, v, rn, xuz(1,2), 
     &                           narrowwidth, 0, m2s)
         call m2s_VVsum(v(0,3,1),3,4,2,1,m2s_dec) !use v3..v6 for H->VV
         m2s_dec = m2s_dec/(16.0d0 * pi * xmgs(6))
         m2s = m2s*m2s_dec
!         print*, 'decfac', m2s_dec
c     pp -> W HJ -> W Z Z 
      CASE(WMHj_ZZ_ll)
         rn = RandomNumber()
         call m2s_qqWHjj( 46, 1, lokt, x, p, v, rn, xuz(1,2), 
     &                           narrowwidth, 0, m2s)
         call m2s_VVsum(v(0,3,1),2,2,2,2,m2s_dec)
         m2s_dec = m2s_dec/(16.0d0 * pi * xmgs(6))
         m2s = m2s*m2s_dec
      CASE(WMHj_ZZ_lnu)
         rn = RandomNumber()
         call m2s_qqWHjj( 46, 1, lokt, x, p, v, rn, xuz(1,2), 
     &                           narrowwidth, 0, m2s)
         call m2s_VVsum(v(0,3,1),2,2,2,1,m2s_dec)
         m2s_dec = m2s_dec/(16.0d0 * pi * xmgs(6))
         m2s = m2s*m2s_dec  


c pp -> W+ j -> e+ ve j
      CASE(WPj)
         rn = RandomNumber()
         call m2s_qqWHjj( 3, 1, lokt, x, p, v, rn, xuz(1,2),
     &      narrowwidth, 0, m2s)

c pp -> W- j -> e- ve~ j
      CASE(WMj)
         rn = RandomNumber()
         call m2s_qqWHjj( 4, 1, lokt, x, p, v, rn, xuz(1,2),
     &      narrowwidth, 0, m2s)
#endif

#ifdef WITH_QCDVV
**************************************
*
*     QCD VVjj
*
**************************************        
c pp -> Z A jj(j) -> e+ e- A jj(j)
      CASE(QCDZAJJ_l)
         rn3(1) = RandomNumber()
         rn3(2) = RandomNumber()
         rn3(3) = RandomNumber()
         call m2s_QCDzz3j(x,p,v,20,.true.,lokt,xuz,rn3,m2s)
c pp -> Z A jj(j) -> vm~ vm A jj(j)
      CASE(QCDZAJJ_n)
         rn3(1) = RandomNumber()
         rn3(2) = RandomNumber()
         rn3(3) = RandomNumber()
         call m2s_QCDzz3j(x,p,v,23,.true.,lokt,xuz,rn3,m2s)
c pp -> A A jj(j) 
      CASE(QCDAAjj)
         rn3(1) = RandomNumber()
         rn3(2) = RandomNumber()
         rn3(3) = RandomNumber()
         call m2s_QCDzz3j(x,p,v,24,.true.,lokt,xuz,rn3,m2s)
c pp -> Z Z jj(j) -> e+ e- mu+ mu- jj(j)
      CASE(QCDZZJJ_ll)
         rn3(1) = RandomNumber()
         rn3(2) = RandomNumber()
         rn3(3) = RandomNumber()
         call m2s_QCDzz3j(x,p,v,22,.true.,lokt,xuz,rn3,m2s)
c pp -> Z Z jj(j) -> e+ e- vm~ vm jj(j)
      CASE(QCDZZJJ_lnu)
         rn3(1) = RandomNumber()
         rn3(2) = RandomNumber()
         rn3(3) = RandomNumber()
         call m2s_QCDzz3j(x,p,v,21,.true.,lokt,xuz,rn3,m2s)
c pp -> W+ Z jj(j) -> e+ ve mu+ mu- jj(j)
      CASE(QCDWPZjj)
         rn3(1) = RandomNumber()
         rn3(2) = RandomNumber()
         call m2s_QCDwz3j(x,p,v,32,.true.,lokt,xuz,rn3(1),m2s)
c pp -> W+ Z jj(j) -> e- ve~ mu- mu+ jj(j)
      CASE(QCDWMZjj)
         rn3(1) = RandomNumber()
         rn3(2) = RandomNumber()
         call m2s_QCDwz3j(x,p,v,42,.true.,lokt,xuz,rn3(1),m2s)
c pp -> W+ Z jj(j) -> e+ ve mu+ mu- jj(j)
      CASE(QCDWPAjj)
         rn3(1) = RandomNumber()
         rn3(2) = RandomNumber()
         call m2s_QCDwz3j(x,p,v,31,.true.,lokt,xuz,rn3(1),m2s)
c pp -> W+ Z jj(j) -> e- ve~ mu- mu+ jj(j)
      CASE(QCDWMAjj)
         rn3(1) = RandomNumber()
         rn3(2) = RandomNumber()
         call m2s_QCDwz3j(x,p,v,41,.true.,lokt,xuz,rn3(1),m2s)
c pp -> W+ W+ jj(j) -> e+ ve mu+ vm jj(j)
      CASE(QCDWPWPjj)
         rn = RandomNumber()
         call m2s_QCDww3j(x,p,v,33,.true.,lokt,xuz,rn,m2s)
c pp -> W- W- jj(j) -> e- ve mu- vm jj(j)
      CASE(QCDWMWMjj)
         rn = RandomNumber()
         call m2s_QCDww3j(x,p,v,44,.true.,lokt,xuz,rn,m2s)
#endif


#ifdef WITH_QCDV
**************************************
*
*     QCD Vjj
*
**************************************   
c pp -> W+ jj(j) -> e+ ve mu+ mu- jj(j)
      CASE(QCDWPjj)
         rn3(1) = RandomNumber()
         rn3(2) = RandomNumber()
         call m2s_QCDwz3j(x,p,v,3,.true.,lokt,xuz,rn3(1),m2s)
c pp -> W- jj(j) -> e- ve mu- mu+ jj(j)
      CASE(QCDWMjj)
         rn3(1) = RandomNumber()
         rn3(2) = RandomNumber()
         call m2s_QCDwz3j(x,p,v,4,.true.,lokt,xuz,rn3(1),m2s)
c pp -> Z jj(j) -> e+ e- jj(j)
      CASE(QCDZjj_l)
         rn3(1) = RandomNumber()
         rn3(2) = RandomNumber()
         rn3(3) = RandomNumber()
         call m2s_QCDzz3j(x,p,v,2,.true.,lokt,xuz,rn3,m2s)
c pp -> Z jj(j) -> ve ve jj(j)
      CASE(QCDZjj_nu)
         rn3(1) = RandomNumber()
         rn3(2) = RandomNumber()
         rn3(3) = RandomNumber()
         call m2s_QCDzz3j(x,p,v,-2,.true.,lokt,xuz,rn3,m2s)
#endif


#ifdef WITH_TRIBOSON
*****************************************************************************
*
*     TRIBOSON : WWZj, ZZZj, ZZWj, WWWj, WWAj, ZZAj, WZAj, WAAj, ZAAj, AAAj
*
*****************************************************************************
c pp -> W+ W- Z -> ve e+ e- ve~ e- e+
      CASE(WWZ, WPhadWMZ, WPWMhadZ, WWZhad)
         rn = RandomNumber()
         call m2s_qqWWZj(342,1,lokt,x,p,v,rn,xuz(1,2),m2s,narrowwidth)

c pp -> Z Z Z -> e- e+ e- e+ e- e+
      CASE(ZZZ, ZZZhad)
         rn = RandomNumber()
         call m2s_qqZZZj(222, 1, lokt, x, p, v, rn, xuz(1,2), m2s,
     1                                       narrowwidth, ps_number)

c pp -> Z Z W+ -> e- e+ e- e+ ve e+
      CASE(ZZWP, ZZhadWP, ZZWPhad)
         rn = RandomNumber()
         call m2s_qqZZWj(223,1,lokt,x,p,v,rn,xuz(1,2),m2s,narrowwidth)

c pp -> Z Z W- -> e- e+ e- e+ e- ve~
      CASE(ZZWM, ZZhadWM, ZZWMhad)
         rn = RandomNumber()
         call m2s_qqZZWj(224,1,lokt,x,p,v,rn,xuz(1,2),m2s,narrowwidth)

c pp -> W+ W- W+ -> ve e+ e- ve~ ve e+
      CASE(WWWP, WPhadWMWP, WPWMhadWP)
         rn = RandomNumber()
         call m2s_qqWWWj(343,1,lokt,x,p,v,rn,xuz(1,2),m2s,narrowwidth)

c pp -> W- W+ W- -> e- ve~ ve e+ e- ve~
      CASE(WWWM, WMhadWPWM, WMWPhadWM)
         rn = RandomNumber()
         call m2s_qqWWWj(434,1,lokt,x,p,v,rn,xuz(1,2),m2s,narrowwidth)

c pp -> W+ W- A -> ve e+ e- ve~ A
      CASE(WWA, WPhadWMA, WPWMhadA)
         rn = RandomNumber()
         call m2s_qqVVV1j( 341, 1, lokt, x, p, v, rn, xuz(1,2), m2s)

c pp -> Z Z A -> e- e+ e- e+ A
      CASE(ZZA, ZZhadA)
         rn = RandomNumber()
         call m2s_qqZZAj( 221, 1, lokt, x, p, v, rn, xuz(1,2), m2s)

c pp -> Z Z A -> e- e+ ve ve~ A
      CASE(ZZnA)
         rn = RandomNumber()
         call m2s_qqZZAj( 251, 1, lokt, x, p, v, rn, xuz(1,2), m2s)

c pp -> W+ Z A -> e+ ve e- e+ A
      CASE(WPZA,WPhadZA, WPZhadA)
         rn = RandomNumber()
         call m2s_qqWZAj( 321, 1, lokt, x, p, v, rn, xuz(1,2), m2s, ps_number)

c pp -> W- Z A -> e- ve~ e- e+ A
      CASE(WMZA, WMhadZA, WMZhadA)
         rn = RandomNumber()
         call m2s_qqWZAj( 421, 1, lokt, x, p, v, rn, xuz(1,2), m2s, ps_number)

c pp -> W+ A A -> e+ ve A A
      CASE(WPAA)
         rn = RandomNumber()
         call m2s_qqWAAj( 311, 1, lokt, x, p, v, rn, xuz(1,2), m2s, 
     -        ps_number, narrowwidth)

c pp -> W- A A -> e- ve~ A A
      CASE(WMAA)
         rn = RandomNumber()
         call m2s_qqWAAj( 411, 1, lokt, x, p, v, rn, xuz(1,2), m2s, 
     -        ps_number, narrowwidth)

c pp -> Z A A -> e+ e- A A
      CASE(ZAA)
         rn = RandomNumber()
         call m2s_qqZAAj( 211, 1, lokt, x, p, v, rn, xuz(1,2), m2s, 
     &        ps_number, narrowwidth)

c pp -> Z A A -> ve~ ve A A
      CASE(ZnAA)
         rn = RandomNumber()
         call m2s_qqZAAj( 511, 1, lokt, x, p, v, rn, xuz(1,2), m2s, 
     &        ps_number, narrowwidth)

c pp -> A A A 
      CASE(AAA)
         rn = RandomNumber()
         call m2s_qqAAAj( 111, 1, lokt, x, p, v, rn, xuz(1,2), m2s, 
     &        ps_number)

#endif

#ifdef WITH_TRIBOSONJET
*****************************************************************************
*
*     TRIBOSON PLUS JET: WAAjj
*
*****************************************************************************
c pp -> W+ A A j -> e+ ve A A j
      CASE(WPAAj)
         rn = RandomNumber()
         call m2s_qqWAAjj( 311, 1, lokt, x, p, v, rn, xuz(1,2), ps_number
     &        ,narrowwidth, m2s)

c pp -> W- A A j -> e- ve~ A A j
      CASE(WMAAj)
         rn = RandomNumber()
         call m2s_qqWAAjj( 411, 1, lokt, x, p, v, rn, xuz(1,2), ps_number
     &        ,narrowwidth , m2s)
#endif

      CASE DEFAULT
         print *,"Error: Real emissions Amplitude ",
     &           "for process not implemented, ID = ",procID
         STOP
      END SELECT

      RETURN
      END


c*************************************************************************
      SUBROUTINE QED_Amplitude(rand, p, x, v, lokt, xuz, m2s, ps_number)
c*************************************************************************
c     Amplitude for real qed emissions.
c*************************************************************************
      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/process.inc"

      double precision x(nx) 
      double precision rand(max_PS_dim)
      double precision p(0:3,max_p,max_kin)
      double precision v(0:3,max_v,max_kin)
      double precision rn, xuz(2,2:max_kin), m2s(0:max_kin), m2s_dec
      logical lokt(max_kin)
      integer L, ps_number

c     masses and widths
      double precision XM2s(6),XMGs(6)
      COMMON /BKOPOUshort/ XM2s,XMGs
      double precision BWNE,BWUD,BWTB,BZNN,BZEE,BZUU,BZDD,BZTT,
     &                BHWW,BHZZ,BHGG,BHTT,BHBB,BHCC,BHTAU,BHMU,
     &                BHGAM, BHGAMZ, XGW, XGZ, XGH, GAMT
      COMMON /BRANCH/ BWNE,BWUD,BWTB,BZNN,BZEE,BZUU,BZDD,BZTT,
     &                BHWW,BHZZ,BHGG,BHTT,BHBB,BHCC,BHTAU,BHMU,
     &                BHGAM, BHGAMZ, XGW, XGZ, XGH, GAMT

      double precision RandomNumber
      external RandomNumber

      SELECT CASE(process)

#ifdef WITH_VBF
**************************************
*
*     VBF : HjjA
*
**************************************
c     pp -> HjjA
      CASE(Hjj)
         rn = RandomNumber()
         call m2s_wbfHA( 6,1,lokt, x, p, v, xuz, m2s)
c     pp -> HjjA -> gamma gamma jjA
      CASE(Hjj_AA)
         rn = RandomNumber()
         call m2s_wbfHA( 6,1,lokt, x, p, v, xuz, m2s)
         do L = 0,n_qed
            m2s(L) = m2s(L)*BHGAM
         enddo
c     pp -> HjjA -> mu+mu- jjA
      CASE(Hjj_mu)
         rn = RandomNumber()
         call m2s_wbfHA( 6,1,lokt, x, p, v, xuz, m2s)
         do L = 0,n_qed
            m2s(L) = m2s(L)*BHMU
         enddo
c     pp -> HjjA -> tau+tau- jjA
      CASE(Hjj_tau)
         rn = RandomNumber()
         call m2s_wbfHA( 6,1,lokt, x, p, v, xuz, m2s)
         do L = 0,n_qed
            m2s(L) = m2s(L)*BHTAU
         enddo
c     pp -> HjjA -> b bbar jjA
      CASE(Hjj_bbar)
         rn = RandomNumber()
         call m2s_wbfHA( 6,1,lokt, x, p, v, xuz, m2s)
         do L = 0,n_qed
            m2s(L) = m2s(L)*BHBB
         enddo
c     pp -> HjjA -> W W jjA
      CASE(Hjj_WW)
         rn = RandomNumber()
         call m2s_wbfHA( 6,1,lokt, x, p, v, xuz, m2s)
         call m2s_VVsum(v(0,1,1),3,4,2,1,m2s_dec)
         m2s_dec = m2s_dec/(16.0d0 * pi * xmgs(6))
         do L = 0,n_qed
            m2s(L) = m2s(L)*m2s_dec
         enddo
      CASE(Hjj_WPhadWM, Hjj_WPWMhad)
         rn = RandomNumber()
         call m2s_wbfHA( 6,1,lokt, x, p, v, xuz, m2s)
         if (process.eq.Hjj_WPhadWM) then
            call m2s_VVsum_had(v(0,1,1),3,4,3,2,m2s_dec)
         else
            call m2s_VVsum_had(v(0,1,1),3,4,1,4,m2s_dec)
         endif
         m2s_dec = m2s_dec/(16.0d0 * pi * xmgs(6))
         do L = 0,n_qed
            m2s(L) = m2s(L)*m2s_dec
         enddo
c     pp -> HjjA -> Z Z jjA
      CASE(Hjj_ZZ_ll)
         rn = RandomNumber()
         call m2s_wbfHA( 6,1,lokt, x, p, v, xuz, m2s)
         call m2s_VVsum(v(0,1,1),2,2,2,2,m2s_dec)
         m2s_dec = m2s_dec/(16.0d0 * pi * xmgs(6))
         do L = 0,n_qed
            m2s(L) = m2s(L)*m2s_dec
         enddo
      CASE(Hjj_ZZhad)
         rn = RandomNumber()
         call m2s_wbfHA( 6,1,lokt, x, p, v, xuz, m2s)
         call m2s_VVsum_had(v(0,1,1),2,2,5,2,m2s_dec)
         m2s_dec = m2s_dec/(16.0d0 * pi * xmgs(6))
         do L = 0,n_qed
            m2s(L) = m2s(L)*m2s_dec
         enddo
      CASE(Hjj_ZZ_lnu)
         rn = RandomNumber()
         call m2s_wbfHA( 6,1,lokt, x, p, v, xuz, m2s)
         call m2s_VVsum(v(0,1,1),2,2,2,1,m2s_dec)
         m2s_dec = m2s_dec/(16.0d0 * pi * xmgs(6))
         do L = 0,n_qed
            m2s(L) = m2s(L)*m2s_dec
         enddo
#endif
      CASE DEFAULT
         print *,"Error: QED emissions Amplitude ",
     &           "for process not implemented, ID = ",procID
         STOP
      END SELECT

      RETURN
      END

c with_nlo
#endif

