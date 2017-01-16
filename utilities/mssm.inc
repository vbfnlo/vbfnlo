!****************************************************
!**                                               ***
!**  Declaration of model parameters.             ***
!**  sophy@particle.uni-karlsruhe.de              ***
!**                                               ***
!****************************************************

!      double precision Pi
!      parameter (Pi = 3.1415926535897932384626433832795029D0)

!* complex i
	double complex cI
	parameter (cI = (0D0,1D0))


!* em coupling constant, electrix charge, fermi constant
      double precision AlfaQED, EL, gfermi, AlfaQED2, delAlfa
!      parameter (gfermi = 0.0000116637D0)


!* Electroweak mixing angle and fermion and guage boson masses and their squares
!* as well as fermion mass array
      double precision SW, SW2, CW, CW2
      double precision MW, MZ, MW2, MZ2

      double precision MsU, MD, MS, MC, MB, MT
      double precision MU2, MD2, MS2, MC2, MB2, MT2

      double precision ME, MM, ML
      double precision ME2, MM2, ML2

      double precision Mf(4,3), Mf2(4,3)


!* Higgs mixing angles, masses, propagator matrix and width
      double precision TB, TB2, CB2, CB, SB2, SB, S2B, C2B
      double complex SA, CA, SA2, CA2, S2A, C2A
      double precision CAB, SAB, CBA, SBA
      double precision MHiggs(4), MHiggs2(4), MHiggs_loop(4), MHiggs_tree(4)
      double precision MA0, MA02, MHp, MHp2
      double precision Mh0, Mh02, MHH, MHH2
      double complex Hgs(3,3), ZHiggs(3,3)
      double precision MH, MH2

!* Tree level higgs mixing angles
      double precision SBAtree, CBAtree, SABtree, CABtree
      double precision TAtree, SAtree, CAtree, SA2tree, CA2tree
      double precision S2Atree, C2Atree

!* Effective higgs mixing angles
      double precision SBAeff, CBAeff, SABeff, CABeff
      double precision TAeff, SAeff, CAeff, SA2eff, CA2eff
      double precision S2Aeff, C2Aeff


!* Squark masses and mixing parameters
      double precision MSf(2,4,3), MSf2(2,4,3)
      double complex USf(2,2,4,3), USfC(2,2,4,3)

!* SUSY breaking parameters
      double complex Af(4,3,3), MUE, AfC(4,3,3), MUEC

!* Chargino and neutralino parameters
      double precision MCha(2), MNeu(4)
      double complex UCha(2,2), VCha(2,2), ZNeu(4,4)
      double precision MCha2(2), MNeu2(4)
      double complex UChaC(2,2), VChaC(2,2), ZNeuC(4,4)

!* The soft-SUSY breaking parameters in the sfermion sector
      double precision M3SL, M2SL, M1SL
      double precision M3SE, M2SE, M1SE
      double precision M3SQ, M2SQ, M1SQ
      double precision M3SU, M2SU, M1SU
      double precision M3SD, M2SD, M1SD

!* The gaugino mass parameters
      double complex M_1, M_2, M_3


!* Coupling factors: the values that the charges and left/right axial vector 
!* couplings can take
      double precision Qu, Qd, gL3, gL4, gR3, gR4, dgL3, dgL4, dgR3, dgR4
      double precision dgL3fm, dgL4fm, dgR3fm, dgR4fm

!* neutralino couplings: values
!* labelling: (neu1, gen)
      double complex SNeu1u(4,2), SNeu1d(4,2)
      double precision SNeu2u, SNeu2d, SNeu3u, SNeu3d

!* chargino couplings: values
!* labelling (cha1,gen); (cha1)
      double complex SCha1u(2,2), SCha1d(2,2), SCha2u(2), SCha2d(2)

!* couplings chargino-neutralino-W: values
!* labelling (cnw, cha1, neu1); cnw = 1-4
      double complex cnw(4,2,4)

!* couplings chargino-chargino-Z: values
!* labelling (cha1,cha2,zcc); zcc = 1-2
      double complex zcc(2,2,2)

!* couplings higgs-squark-squark
!* labelling (gen)
      double complex Shqku(2), Shqkd(2), Shqk1u, Shqk1d
      double complex SHHu(2), SHHd(2), SHH1u, SHH1d
      double complex SAqku, SAqkd


!* Renormalisation parameters: mass^2 of photon, mass of quark, UV divergence 
      double precision massG, massQ, Divergence


!* Renormalization constants - for definition, see notebook 3, page 77
      double complex dZfLu, dZfLd, dZfLc, dZfLs
      double complex dZfRu, dZfRd, dZfRc, dZfRs
      double complex dMWsq1, dMZsq1, dSW1, dZet1, dZh0
      double complex dZH1, dZHH, dZh0HH, dZW1, dZZZ1, dTB1, dMAZsq1
      double complex dZZA1, dZAA1, dZAZ1


!* renormalization constants - fermion contributions (used if we're calculating
!* ew-loop^2 contributions)
      double complex dMWsq1fm, dMZsq1fm, dSW1fm, dZet1fm, dZh0fm
      double complex dZH1fm, dZHHfm, dZh0HHfm, dZW1fm, dZZZ1fm, dTB1fm
      double complex dZZA1fm, dZAA1fm, dZAZ1fm, dMAZsq1fm


!* Renormalization constants - full MSSM, used when Higgs-VV vertex is 
!* calculated in the full MSSM, but otherwise only SM-type and sfermion
!* contributions are included
      double complex dMWsq1full, dMZsq1full, dSW1full, dZet1full
      double complex dZZA1full, dZW1full, dZZZ1full
      double complex dZHHfull, dZh0full, dZh0HHfull, dTB1full


!* massCur stores all of the external particle masses for each current - the 
!* first entry is particle number i, the second is current label cur.  See 
!* notebook 7, page 1 for notation. chargeCur stores charges for each current 
!* when qk = 1.  Label = (cur, particle i)
      double precision massCur(5,6), chargeCur(6,5)
      double precision qqVmass(2,2,2), qqVcharge(2,2)


!* Flag for benchmark
      integer bench
!* Flag for Higgs propagator corns
      integer propLevel
!* FeynHiggs flag for higgs-mixing
      integer higgsmix


      common /susyparamsD/ AlfaQED, EL, AlfaQED2, gfermi, delAlfa

      common /susyparamsD/ SW, SW2, CW, CW2
      common /susyparamsD/ MW, MZ, MW2, MZ2

      common /susyparamsD/ MsU, MD, MS, MC, MB, MT
      common /susyparamsD/ MU2, MD2, MS2, MC2, MB2, MT2
      common /susyparamsD/ ME, MM, ML, ME2, MM2, ML2
      common /susyparamsD/ Mf, Mf2


      common /susyparamsD/ TB, TB2, CB2, CB, SB2, SB, S2B, C2B
      common /susyparamsC/ SA, CA, SA2, CA2, S2A, C2A
      common /susyparamsD/ CAB, SAB, CBA, SBA
      common /susyparamsD/ MHiggs, MHiggs2, MHiggs_loop, MHiggs_tree
      common /susyparamsD/ MA0, MA02, MHp, MHp2
      common /susyparamsD/ Mh0, Mh02, MHH, MHH2
      common /susyparamsC/ Hgs, ZHiggs
      common /susyparamsD/ MH, MH2

      common /susyparamsD/ SBAtree, CBAtree, SABtree, CABtree
      common /susyparamsD/ TAtree, SAtree, CAtree, SA2tree, CA2tree
      common /susyparamsD/ S2Atree, C2Atree

      common /susyparamsD/ SBAeff, CBAeff, SABeff, CABeff
      common /susyparamsD/ TAeff, SAeff, CAeff, SA2eff, CA2eff
      common /susyparamsD/ S2Aeff, C2Aeff

      common /susyparamsD/ MSf, MSf2
      common /susyparamsC/ USf, USfC, Af, MUE, AfC, MUEC

      common /susyparamsD/ MCha, MNeu
      common /susyparamsC/ UCha, VCha, ZNeu
      common /susyparamsD/ MCha2, MNeu2
      common /susyparamsC/ UChaC, VChaC, ZNeuC

      common /susyparamsD/ M3SL, M2SL, M1SL
      common /susyparamsD/ M3SE, M2SE, M1SE
      common /susyparamsD/ M3SQ, M2SQ, M1SQ
      common /susyparamsD/ M3SU, M2SU, M1SU
      common /susyparamsD/ M3SD, M2SD, M1SD

      common /susyparamsC/ M_1, M_2, M_3


      common /susyparamsD/ Qu, Qd, gL3, gL4, gR3, gR4, dgL3, dgL4, dgR3, dgR4
      common /susyparamsD/ massCur, chargeCur, qqVmass, qqVcharge
      common /susyparamsC/ SNeu1u, SNeu1d
      common /susyparamsD/ SNeu2u, SNeu2d, SNeu3u, SNeu3d
      common /susyparamsC/ SCha1u, SCha1d, SCha2u, SCha2d
      common /susyparamsC/ cnw, zcc
      common /susyparamsC/ Shqku, Shqkd, Shqk1u, Shqk1d
      common /susyparamsC/ SHHu, SHHd, SHH1u, SHH1d
      common /susyparamsC/ SAqku, SAqkd

      common /susyparamsD/ massG, massQ, Divergence

      common /susyparamsC/ dMWsq1, dMZsq1, dSW1, dZet1, dZh0
      common /susyparamsC/ dZH1, dZHH, dZh0HH, dZW1, dZZZ1, dTB1, dMAZsq1
      common /susyparamsC/ dZZA1, dZAA1, dZAZ1
      common /susyparamsC/ dZfLu, dZfLd, dZfLc, dZfLs
      common /susyparamsC/ dZfRu, dZfRd, dZfRc, dZfRs

      common /susyparamsC/ dMWsq1full, dMZsq1full, dSW1full, dZet1full
      common /susyparamsC/ dZZA1full, dZW1full, dZZZ1full
      common /susyparamsC/ dZHHfull, dZh0full, dZh0HHfull, dTB1full

      common /susyparamsI/ bench, propLevel, higgsmix


      common /susyparamsFM/ dMWsq1fm, dMZsq1fm, dSW1fm, dZet1fm, dZh0fm
      common /susyparamsFM/ dZH1fm, dZHHfm, dZh0HHfm, dZW1fm, dZZZ1fm, dTB1fm
      common /susyparamsFM/ dZZA1fm, dZAA1fm, dZAZ1fm, dMAZsq1fm

      common /susyparamsFMd/ dgL3fm, dgL4fm, dgR3fm, dgR4fm



