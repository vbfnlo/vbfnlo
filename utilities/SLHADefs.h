#if 0
	SLHADefs.h
		declarations for SLHALib data
		generated 15 Sep 2009 11:04
#endif

#ifndef SLHADEFS_H
#define SLHADEFS_H

#define invalid (-999)

#define OffsetModSel 0
#define LengthModSel 12
#define BlockModSel(i) SlhaData(i)
#define ModSel_Model Slhadata(1)
#define ModSel_Content Slhadata(2)
#define ModSel_RPV Slhadata(3)
#define ModSel_CPV Slhadata(4)
#define ModSel_FV Slhadata(5)
#define ModSel_GridPts Slhadata(6)
#define ModSel_Qmax Slhadata(7)
#define ModSel_PDG(n) Slhadata(7+n)

#define OffsetSMInputs 12
#define LengthSMInputs 16
#define BlockSMInputs(i) SlhaData(12+i)
#define SMInputs_invAlfaMZ Slhadata(13)
#define SMInputs_GF Slhadata(14)
#define SMInputs_AlfasMZ Slhadata(15)
#define SMInputs_MZ Slhadata(16)
#define SMInputs_Mf(t,g) Slhadata(12+t+4*(g))
#define SMInputs_MfFlat(i) Slhadata(16+i)
#define   SMInputs_Mnu1 SMInputs_Mf(1,1)
#define   SMInputs_Me SMInputs_Mf(2,1)
#define   SMInputs_Mu SMInputs_Mf(3,1)
#define   SMInputs_Md SMInputs_Mf(4,1)
#define   SMInputs_Mnu2 SMInputs_Mf(1,2)
#define   SMInputs_Mmu SMInputs_Mf(2,2)
#define   SMInputs_Mc SMInputs_Mf(3,2)
#define   SMInputs_Ms SMInputs_Mf(4,2)
#define   SMInputs_Mnu3 SMInputs_Mf(1,3)
#define   SMInputs_Mtau SMInputs_Mf(2,3)
#define   SMInputs_Mt SMInputs_Mf(3,3)
#define   SMInputs_Mb SMInputs_Mf(4,3)

#define OffsetMinPar 28
#define LengthMinPar 6
#define BlockMinPar(i) SlhaData(28+i)
#define MinPar_M0 Slhadata(29)
#define   MinPar_Lambda MinPar_M0
#define MinPar_M12 Slhadata(30)
#define   MinPar_Mmess MinPar_M12
#define   MinPar_M32 MinPar_M12
#define MinPar_TB Slhadata(31)
#define MinPar_signMUE Slhadata(32)
#define MinPar_A Slhadata(33)
#define   MinPar_N5 MinPar_A
#define MinPar_cgrav Slhadata(34)

#define OffsetExtPar 34
#define LengthExtPar 42
#define BlockExtPar(i) SlhaData(34+i)
#define ExtPar_Q SlhaData(35)
#define ExtPar_M1 Slhadata(36)
#define ExtPar_M2 Slhadata(37)
#define ExtPar_M3 Slhadata(38)
#define ExtPar_Af(t) Slhadata(37+t)
#define   ExtPar_Atau ExtPar_Af(2)
#define   ExtPar_At ExtPar_Af(3)
#define   ExtPar_Ab ExtPar_Af(4)
#define ExtPar_MHu2 Slhadata(42)
#define ExtPar_MHd2 Slhadata(43)
#define ExtPar_MUE Slhadata(44)
#define ExtPar_MA02 Slhadata(45)
#define ExtPar_TB Slhadata(46)
#define ExtPar_MA0 Slhadata(47)
#define ExtPar_MHp Slhadata(48)
#define ExtPar_MSS(g,q) Slhadata(45+g+3*(q))
#define   ExtPar_MSL(g) ExtPar_MSS(g,1)
#define   ExtPar_MSE(g) ExtPar_MSS(g,2)
#define   ExtPar_MSQ(g) ExtPar_MSS(g,3)
#define   ExtPar_MSU(g) ExtPar_MSS(g,4)
#define   ExtPar_MSD(g) ExtPar_MSS(g,5)
#define ExtPar_N5(g) Slhadata(63+g)
#define ExtPar_lambda Slhadata(67)
#define ExtPar_kappa Slhadata(68)
#define ExtPar_Alambda Slhadata(69)
#define ExtPar_Akappa Slhadata(70)
#define ExtPar_lambdaS Slhadata(71)
#define ExtPar_xiF Slhadata(72)
#define ExtPar_xiS Slhadata(73)
#define ExtPar_MUEprime Slhadata(74)
#define ExtPar_mS2prime Slhadata(75)
#define ExtPar_mS2 Slhadata(76)

#define OffsetQExtPar 76
#define LengthQExtPar 16
#define BlockQExtPar(i) SlhaData(76+i)
#define QExtPar_QM1 Slhadata(77)
#define QExtPar_QM2 Slhadata(78)
#define QExtPar_QM3 Slhadata(79)
#define QExtPar_QAf(t) Slhadata(78+t)
#define   QExtPar_QAtau QExtPar_QAf(2)
#define   QExtPar_QAt QExtPar_QAf(3)
#define   QExtPar_QAb QExtPar_QAf(4)
#define QExtPar_QMHu2 Slhadata(83)
#define QExtPar_QMHd2 Slhadata(84)
#define QExtPar_QMUE Slhadata(85)
#define QExtPar_QMA02 Slhadata(86)
#define QExtPar_QTB Slhadata(87)
#define QExtPar_QMSS(q) Slhadata(87+q)
#define   QExtPar_QMSL QExtPar_QMSS(1)
#define   QExtPar_QMSE QExtPar_QMSS(2)
#define   QExtPar_QMSQ QExtPar_QMSS(3)
#define   QExtPar_QMSU QExtPar_QMSS(4)
#define   QExtPar_QMSD QExtPar_QMSS(5)

#define OffsetNMSSMRun 92
#define LengthNMSSMRun 11
#define BlockNMSSMRun(i) SlhaData(92+i)
#define NMSSMRun_Q SlhaData(93)
#define NMSSMRun_lambda Slhadata(94)
#define NMSSMRun_kappa Slhadata(95)
#define NMSSMRun_Alambda Slhadata(96)
#define NMSSMRun_Akappa Slhadata(97)
#define NMSSMRun_lambdaS Slhadata(98)
#define NMSSMRun_xiF Slhadata(99)
#define NMSSMRun_xiS Slhadata(100)
#define NMSSMRun_MUEprime Slhadata(101)
#define NMSSMRun_mS2prime Slhadata(102)
#define NMSSMRun_mS2 Slhadata(103)

#define OffsetMass 103
#define LengthMass 53
#define BlockMass(i) SlhaData(103+i)
#define Mass_Mf(t,g) Slhadata(99+t+4*(g))
#define Mass_MfFlat(i) Slhadata(103+i)
#define Mass_MSf(s,t,g) Slhadata(105+s+8*(g)+2*(t))
#define Mass_MSfFlat(i) Slhadata(115+i)
#define Mass_MZ Slhadata(140)
#define Mass_MW Slhadata(141)
#define Mass_Mh0 Slhadata(142)
#define Mass_MHH Slhadata(143)
#define Mass_MA0 Slhadata(144)
#define Mass_MHp Slhadata(145)
#define   Mass_MH1 Mass_Mh0
#define   Mass_MH2 Mass_MHH
#define Mass_MH3 Slhadata(146)
#define   Mass_MA1 Mass_MA0
#define Mass_MA2 Slhadata(147)
#define Mass_MNeu(n) Slhadata(147+n)
#define Mass_MCha(c) Slhadata(152+c)
#define Mass_MGl Slhadata(155)
#define Mass_MGrav Slhadata(156)

#define OffsetDMass 156
#define LengthDMass 4
#define BlockDMass(i) SlhaData(156+i)
#define DMass_DeltaMh0 Slhadata(157)
#define DMass_DeltaMHH Slhadata(158)
#define DMass_DeltaMA0 Slhadata(159)
#define DMass_DeltaMHp Slhadata(160)

#define OffsetNMix 160
#define LengthNMix 16
#define BlockNMix(i) SlhaData(160+i)
#define NMix_ZNeu(n1,n2) Slhadata(156+n1+4*(n2))
#define NMix_ZNeuFlat(i) Slhadata(160+i)

#define OffsetUMix 176
#define LengthUMix 4
#define BlockUMix(i) SlhaData(176+i)
#define UMix_UCha(c1,c2) Slhadata(174+c1+2*(c2))
#define UMix_UChaFlat(i) Slhadata(176+i)

#define OffsetVMix 180
#define LengthVMix 4
#define BlockVMix(i) SlhaData(180+i)
#define VMix_VCha(c1,c2) Slhadata(178+c1+2*(c2))
#define VMix_VChaFlat(i) Slhadata(180+i)

#define OffsetSfMix 184
#define LengthSfMix 12
#define BlockSfMix(i) SlhaData(184+i)
#define SfMix_USf(s1,s2,t) Slhadata(174+s1+2*(s2)+4*(t))
#define SfMix_USfFlat(i,t) Slhadata(176+i+4*(t))

#define OffsetStauMix 184
#define LengthStauMix 4
#define BlockStauMix(i) SlhaData(184+i)
#define   StauMix_USf(s1,s2) SfMix_USf(s1,s2,2)
#define   StauMix_USfFlat(i) SfMix_USfFlat(i,2)

#define OffsetStopMix 188
#define LengthStopMix 4
#define BlockStopMix(i) SlhaData(188+i)
#define   StopMix_USf(s1,s2) SfMix_USf(s1,s2,3)
#define   StopMix_USfFlat(i) SfMix_USfFlat(i,3)

#define OffsetSbotMix 192
#define LengthSbotMix 4
#define BlockSbotMix(i) SlhaData(192+i)
#define   SbotMix_USf(s1,s2) SfMix_USf(s1,s2,4)
#define   SbotMix_USfFlat(i) SfMix_USfFlat(i,4)

#define OffsetAlpha 196
#define LengthAlpha 1
#define BlockAlpha(i) SlhaData(196+i)
#define Alpha_Alpha Slhadata(197)

#define OffsetDAlpha 197
#define LengthDAlpha 1
#define BlockDAlpha(i) SlhaData(197+i)
#define DAlpha_DeltaAlpha Slhadata(198)

#define OffsetHMix 198
#define LengthHMix 5
#define BlockHMix(i) SlhaData(198+i)
#define HMix_Q SlhaData(199)
#define HMix_MUE Slhadata(200)
#define HMix_TB Slhadata(201)
#define HMix_VEV Slhadata(202)
#define HMix_MA02 Slhadata(203)

#define OffsetGauge 203
#define LengthGauge 4
#define BlockGauge(i) SlhaData(203+i)
#define Gauge_Q SlhaData(204)
#define Gauge_g1 Slhadata(205)
#define Gauge_g2 Slhadata(206)
#define Gauge_g3 Slhadata(207)

#define OffsetMSoft 207
#define LengthMSoft 21
#define BlockMSoft(i) SlhaData(207+i)
#define MSoft_Q SlhaData(208)
#define MSoft_M1 Slhadata(209)
#define MSoft_M2 Slhadata(210)
#define MSoft_M3 Slhadata(211)
#define MSoft_MHu2 Slhadata(212)
#define MSoft_MHd2 Slhadata(213)
#define MSoft_MSS(g,q) Slhadata(210+g+3*(q))
#define   MSoft_MSL(g) MSoft_MSS(g,1)
#define   MSoft_MSE(g) MSoft_MSS(g,2)
#define   MSoft_MSQ(g) MSoft_MSS(g,3)
#define   MSoft_MSU(g) MSoft_MSS(g,4)
#define   MSoft_MSD(g) MSoft_MSS(g,5)

#define OffsetAf 228
#define LengthAf 30
#define BlockAf(i) SlhaData(228+i)
#define Af_Q(t) Slhadata(209+10*(t))
#define Af_Af(g1,g2,t) Slhadata(206+g1+3*(g2)+10*(t))
#define Af_AfFlat(i,t) Slhadata(209+i+10*(t))

#define OffsetAe 228
#define LengthAe 11
#define BlockAe(i) SlhaData(228+i)
#define   Ae_Q Af_Q(2)
#define   Ae_Af(g1,g2) Af_Af(g1,g2,2)
#define   Ae_AfFlat(i) Af_AfFlat(i,2)
#define   Ae_Atau Ae_Af(3,3)

#define OffsetAu 239
#define LengthAu 11
#define BlockAu(i) SlhaData(239+i)
#define   Au_Q Af_Q(3)
#define   Au_Af(g1,g2) Af_Af(g1,g2,3)
#define   Au_AfFlat(i) Af_AfFlat(i,3)
#define   Au_At Au_Af(3,3)

#define OffsetAd 250
#define LengthAd 11
#define BlockAd(i) SlhaData(250+i)
#define   Ad_Q Af_Q(4)
#define   Ad_Af(g1,g2) Af_Af(g1,g2,4)
#define   Ad_AfFlat(i) Af_AfFlat(i,4)
#define   Ad_Ab Ad_Af(3,3)

#define OffsetYf 261
#define LengthYf 30
#define BlockYf(i) SlhaData(261+i)
#define Yf_Q(t) Slhadata(242+10*(t))
#define Yf_Yf(g1,g2,t) Slhadata(239+g1+3*(g2)+10*(t))
#define Yf_YfFlat(i,t) Slhadata(242+i+10*(t))

#define OffsetYe 261
#define LengthYe 11
#define BlockYe(i) SlhaData(261+i)
#define   Ye_Q Yf_Q(2)
#define   Ye_Yf(g1,g2) Yf_Yf(g1,g2,2)
#define   Ye_YfFlat(i) Yf_YfFlat(i,2)
#define   Ye_Ytau Ye_Yf(3,3)

#define OffsetYu 272
#define LengthYu 11
#define BlockYu(i) SlhaData(272+i)
#define   Yu_Q Yf_Q(3)
#define   Yu_Yf(g1,g2) Yf_Yf(g1,g2,3)
#define   Yu_YfFlat(i) Yf_YfFlat(i,3)
#define   Yu_Yt Yu_Yf(3,3)

#define OffsetYd 283
#define LengthYd 11
#define BlockYd(i) SlhaData(283+i)
#define   Yd_Q Yf_Q(4)
#define   Yd_Yf(g1,g2) Yf_Yf(g1,g2,4)
#define   Yd_YfFlat(i) Yf_YfFlat(i,4)
#define   Yd_Yb Yd_Yf(3,3)

#define OffsetRVLamLLEIn 294
#define LengthRVLamLLEIn 27
#define BlockRVLamLLEIn(i) SlhaData(294+i)
#define RVLamLLEIn_lamLLE(i,j,k) Slhadata(282+i+3*(j)+9*(k))
#define RVLamLLEIn_lamLLEFlat(i) Slhadata(294+i)

#define OffsetRVLamLQDIn 321
#define LengthRVLamLQDIn 27
#define BlockRVLamLQDIn(i) SlhaData(321+i)
#define RVLamLQDIn_lamLQD(i,j,k) Slhadata(309+i+3*(j)+9*(k))
#define RVLamLQDIn_lamLQDFlat(i) Slhadata(321+i)

#define OffsetRVLamUDDIn 348
#define LengthRVLamUDDIn 27
#define BlockRVLamUDDIn(i) SlhaData(348+i)
#define RVLamUDDIn_lamUDD(i,j,k) Slhadata(336+i+3*(j)+9*(k))
#define RVLamUDDIn_lamUDDFlat(i) Slhadata(348+i)

#define OffsetRVLamLLE 375
#define LengthRVLamLLE 28
#define BlockRVLamLLE(i) SlhaData(375+i)
#define RVLamLLE_Q SlhaData(376)
#define RVLamLLE_lamLLE(i,j,k) Slhadata(364+i+3*(j)+9*(k))
#define RVLamLLE_lamLLEFlat(i) Slhadata(376+i)

#define OffsetRVLamLQD 403
#define LengthRVLamLQD 28
#define BlockRVLamLQD(i) SlhaData(403+i)
#define RVLamLQD_Q SlhaData(404)
#define RVLamLQD_lamLQD(i,j,k) Slhadata(392+i+3*(j)+9*(k))
#define RVLamLQD_lamLQDFlat(i) Slhadata(404+i)

#define OffsetRVLamUDD 431
#define LengthRVLamUDD 28
#define BlockRVLamUDD(i) SlhaData(431+i)
#define RVLamUDD_Q SlhaData(432)
#define RVLamUDD_lamUDD(i,j,k) Slhadata(420+i+3*(j)+9*(k))
#define RVLamUDD_lamUDDFlat(i) Slhadata(432+i)

#define OffsetRVTLLEIn 459
#define LengthRVTLLEIn 27
#define BlockRVTLLEIn(i) SlhaData(459+i)
#define RVTLLEIn_TLLE(i,j,k) Slhadata(447+i+3*(j)+9*(k))
#define RVTLLEIn_TLLEFlat(i) Slhadata(459+i)

#define OffsetRVTLQDIn 486
#define LengthRVTLQDIn 27
#define BlockRVTLQDIn(i) SlhaData(486+i)
#define RVTLQDIn_TLQD(i,j,k) Slhadata(474+i+3*(j)+9*(k))
#define RVTLQDIn_TLQDFlat(i) Slhadata(486+i)

#define OffsetRVTUDDIn 513
#define LengthRVTUDDIn 27
#define BlockRVTUDDIn(i) SlhaData(513+i)
#define RVTUDDIn_TUDD(i,j,k) Slhadata(501+i+3*(j)+9*(k))
#define RVTUDDIn_TUDDFlat(i) Slhadata(513+i)

#define OffsetRVTLLE 540
#define LengthRVTLLE 28
#define BlockRVTLLE(i) SlhaData(540+i)
#define RVTLLE_Q SlhaData(541)
#define RVTLLE_TLLE(i,j,k) Slhadata(529+i+3*(j)+9*(k))
#define RVTLLE_TLLEFlat(i) Slhadata(541+i)

#define OffsetRVTLQD 568
#define LengthRVTLQD 28
#define BlockRVTLQD(i) SlhaData(568+i)
#define RVTLQD_Q SlhaData(569)
#define RVTLQD_TLQD(i,j,k) Slhadata(557+i+3*(j)+9*(k))
#define RVTLQD_TLQDFlat(i) Slhadata(569+i)

#define OffsetRVTUDD 596
#define LengthRVTUDD 28
#define BlockRVTUDD(i) SlhaData(596+i)
#define RVTUDD_Q SlhaData(597)
#define RVTUDD_TUDD(i,j,k) Slhadata(585+i+3*(j)+9*(k))
#define RVTUDD_TUDDFlat(i) Slhadata(597+i)

#define OffsetRVKappaIn 624
#define LengthRVKappaIn 3
#define BlockRVKappaIn(i) SlhaData(624+i)
#define RVKappaIn_kappa(i) Slhadata(624+i)

#define OffsetRVKappa 627
#define LengthRVKappa 4
#define BlockRVKappa(i) SlhaData(627+i)
#define RVKappa_Q SlhaData(628)
#define RVKappa_kappa(i) Slhadata(628+i)

#define OffsetRVDIn 631
#define LengthRVDIn 3
#define BlockRVDIn(i) SlhaData(631+i)
#define RVDIn_D(i) Slhadata(631+i)

#define OffsetRVD 634
#define LengthRVD 4
#define BlockRVD(i) SlhaData(634+i)
#define RVD_Q SlhaData(635)
#define RVD_D(i) Slhadata(635+i)

#define OffsetRVSnVEVIn 638
#define LengthRVSnVEVIn 3
#define BlockRVSnVEVIn(i) SlhaData(638+i)
#define RVSnVEVIn_VEV(i) Slhadata(638+i)

#define OffsetRVSnVEV 641
#define LengthRVSnVEV 4
#define BlockRVSnVEV(i) SlhaData(641+i)
#define RVSnVEV_Q SlhaData(642)
#define RVSnVEV_VEV(i) Slhadata(642+i)

#define OffsetRVM2LH1In 645
#define LengthRVM2LH1In 3
#define BlockRVM2LH1In(i) SlhaData(645+i)
#define RVM2LH1In_M2LH1(i) Slhadata(645+i)

#define OffsetRVM2LH1 648
#define LengthRVM2LH1 4
#define BlockRVM2LH1(i) SlhaData(648+i)
#define RVM2LH1_Q SlhaData(649)
#define RVM2LH1_M2LH1(i) Slhadata(649+i)

#define OffsetRVNMix 652
#define LengthRVNMix 49
#define BlockRVNMix(i) SlhaData(652+i)
#define RVNMix_ZNeu(n1,n2) Slhadata(645+n1+7*(n2))
#define RVNMix_ZNeuFlat(i) Slhadata(652+i)

#define OffsetRVUMix 701
#define LengthRVUMix 25
#define BlockRVUMix(i) SlhaData(701+i)
#define RVUMix_UCha(c1,c2) Slhadata(696+c1+5*(c2))
#define RVUMix_UChaFlat(i) Slhadata(701+i)

#define OffsetRVVMix 726
#define LengthRVVMix 25
#define BlockRVVMix(i) SlhaData(726+i)
#define RVVMix_VCha(c1,c2) Slhadata(721+c1+5*(c2))
#define RVVMix_VChaFlat(i) Slhadata(726+i)

#define OffsetRVHMix 751
#define LengthRVHMix 25
#define BlockRVHMix(i) SlhaData(751+i)
#define RVHMix_UH(h1,h2) Slhadata(746+h1+5*(h2))
#define RVHMix_UHFlat(i) Slhadata(751+i)

#define OffsetRVAMix 776
#define LengthRVAMix 25
#define BlockRVAMix(i) SlhaData(776+i)
#define RVAMix_UA(h1,h2) Slhadata(771+h1+5*(h2))
#define RVAMix_UAFlat(i) Slhadata(776+i)

#define OffsetRVLMix 801
#define LengthRVLMix 64
#define BlockRVLMix(i) SlhaData(801+i)
#define RVLMix_CLep(l1,l2) Slhadata(793+l1+8*(l2))
#define RVLMix_CLepFlat(i) Slhadata(801+i)

#define OffsetVCKMIn 865
#define LengthVCKMIn 4
#define BlockVCKMIn(i) SlhaData(865+i)
#define VCKMIn_lambda Slhadata(866)
#define VCKMIn_A Slhadata(867)
#define VCKMIn_rhobar Slhadata(868)
#define VCKMIn_etabar Slhadata(869)

#define OffsetVCKM 869
#define LengthVCKM 10
#define BlockVCKM(i) SlhaData(869+i)
#define VCKM_Q SlhaData(870)
#define VCKM_VCKM(g1,g2) Slhadata(867+g1+3*(g2))
#define VCKM_VCKMFlat(i) Slhadata(870+i)

#define OffsetUPMNSIn 879
#define LengthUPMNSIn 6
#define BlockUPMNSIn(i) SlhaData(879+i)
#define UPMNSIn_theta12 Slhadata(880)
#define UPMNSIn_theta23 Slhadata(881)
#define UPMNSIn_theta13 Slhadata(882)
#define UPMNSIn_delta13 Slhadata(883)
#define UPMNSIn_alpha1 Slhadata(884)
#define UPMNSIn_alpha2 Slhadata(885)

#define OffsetUPMNS 885
#define LengthUPMNS 10
#define BlockUPMNS(i) SlhaData(885+i)
#define UPMNS_Q SlhaData(886)
#define UPMNS_UPMNS(g1,g2) Slhadata(883+g1+3*(g2))
#define UPMNS_UPMNSFlat(i) Slhadata(886+i)

#define OffsetMSS2In 895
#define LengthMSS2In 45
#define BlockMSS2In(i) SlhaData(895+i)
#define MSS2In_MSS2(g1,g2,q) Slhadata(883+g1+3*(g2)+9*(q))
#define MSS2In_MSS2Flat(i,q) Slhadata(886+i+9*(q))

#define OffsetMSL2In 895
#define LengthMSL2In 9
#define BlockMSL2In(i) SlhaData(895+i)
#define   MSL2In_MSL2(g1,g2) MSS2In_MSS2(g1,g2,1)
#define   MSL2In_MSL2Flat(i) MSS2In_MSS2Flat(i,1)

#define OffsetMSE2In 904
#define LengthMSE2In 9
#define BlockMSE2In(i) SlhaData(904+i)
#define   MSE2In_MSE2(g1,g2) MSS2In_MSS2(g1,g2,2)
#define   MSE2In_MSE2Flat(i) MSS2In_MSS2Flat(i,2)

#define OffsetMSQ2In 913
#define LengthMSQ2In 9
#define BlockMSQ2In(i) SlhaData(913+i)
#define   MSQ2In_MSQ2(g1,g2) MSS2In_MSS2(g1,g2,3)
#define   MSQ2In_MSQ2Flat(i) MSS2In_MSS2Flat(i,3)

#define OffsetMSU2In 922
#define LengthMSU2In 9
#define BlockMSU2In(i) SlhaData(922+i)
#define   MSU2In_MSU2(g1,g2) MSS2In_MSS2(g1,g2,4)
#define   MSU2In_MSU2Flat(i) MSS2In_MSS2Flat(i,4)

#define OffsetMSD2In 931
#define LengthMSD2In 9
#define BlockMSD2In(i) SlhaData(931+i)
#define   MSD2In_MSD2(g1,g2) MSS2In_MSS2(g1,g2,5)
#define   MSD2In_MSD2Flat(i) MSS2In_MSS2Flat(i,5)

#define OffsetMSS2 940
#define LengthMSS2 50
#define BlockMSS2(i) SlhaData(940+i)
#define MSS2_Q(q) Slhadata(931+10*(q))
#define MSS2_MSS2(g1,g2,q) Slhadata(928+g1+3*(g2)+10*(q))
#define MSS2_MSS2Flat(i,q) Slhadata(931+i+10*(q))

#define OffsetMSL2 940
#define LengthMSL2 10
#define BlockMSL2(i) SlhaData(940+i)
#define   MSL2_Q MSS2_Q(1)
#define   MSL2_MSL2(g1,g2) MSS2_MSS2(g1,g2,1)
#define   MSL2_MSL2Flat(i) MSS2_MSS2Flat(i,1)

#define OffsetMSE2 950
#define LengthMSE2 10
#define BlockMSE2(i) SlhaData(950+i)
#define   MSE2_Q MSS2_Q(2)
#define   MSE2_MSE2(g1,g2) MSS2_MSS2(g1,g2,2)
#define   MSE2_MSE2Flat(i) MSS2_MSS2Flat(i,2)

#define OffsetMSQ2 960
#define LengthMSQ2 10
#define BlockMSQ2(i) SlhaData(960+i)
#define   MSQ2_Q MSS2_Q(3)
#define   MSQ2_MSQ2(g1,g2) MSS2_MSS2(g1,g2,3)
#define   MSQ2_MSQ2Flat(i) MSS2_MSS2Flat(i,3)

#define OffsetMSU2 970
#define LengthMSU2 10
#define BlockMSU2(i) SlhaData(970+i)
#define   MSU2_Q MSS2_Q(4)
#define   MSU2_MSU2(g1,g2) MSS2_MSS2(g1,g2,4)
#define   MSU2_MSU2Flat(i) MSS2_MSS2Flat(i,4)

#define OffsetMSD2 980
#define LengthMSD2 10
#define BlockMSD2(i) SlhaData(980+i)
#define   MSD2_Q MSS2_Q(5)
#define   MSD2_MSD2(g1,g2) MSS2_MSS2(g1,g2,5)
#define   MSD2_MSD2Flat(i) MSS2_MSS2Flat(i,5)

#define OffsetTfIn 990
#define LengthTfIn 27
#define BlockTfIn(i) SlhaData(990+i)
#define TfIn_Tf(g1,g2,t) Slhadata(969+g1+3*(g2)+9*(t))
#define TfIn_TfFlat(i,t) Slhadata(972+i+9*(t))

#define OffsetTeIn 990
#define LengthTeIn 9
#define BlockTeIn(i) SlhaData(990+i)
#define   TeIn_Tf(g1,g2) TfIn_Tf(g1,g2,2)
#define   TeIn_TfFlat(i) TfIn_TfFlat(i,2)

#define OffsetTuIn 999
#define LengthTuIn 9
#define BlockTuIn(i) SlhaData(999+i)
#define   TuIn_Tf(g1,g2) TfIn_Tf(g1,g2,3)
#define   TuIn_TfFlat(i) TfIn_TfFlat(i,3)

#define OffsetTdIn 1008
#define LengthTdIn 9
#define BlockTdIn(i) SlhaData(1008+i)
#define   TdIn_Tf(g1,g2) TfIn_Tf(g1,g2,4)
#define   TdIn_TfFlat(i) TfIn_TfFlat(i,4)

#define OffsetTf 1017
#define LengthTf 30
#define BlockTf(i) SlhaData(1017+i)
#define Tf_Q(t) Slhadata(998+10*(t))
#define Tf_Tf(g1,g2,t) Slhadata(995+g1+3*(g2)+10*(t))
#define Tf_TfFlat(i,t) Slhadata(998+i+10*(t))

#define OffsetTe 1017
#define LengthTe 10
#define BlockTe(i) SlhaData(1017+i)
#define   Te_Q Tf_Q(2)
#define   Te_Tf(g1,g2) Tf_Tf(g1,g2,2)
#define   Te_TfFlat(i) Tf_TfFlat(i,2)

#define OffsetTu 1027
#define LengthTu 10
#define BlockTu(i) SlhaData(1027+i)
#define   Tu_Q Tf_Q(3)
#define   Tu_Tf(g1,g2) Tf_Tf(g1,g2,3)
#define   Tu_TfFlat(i) Tf_TfFlat(i,3)

#define OffsetTd 1037
#define LengthTd 10
#define BlockTd(i) SlhaData(1037+i)
#define   Td_Q Tf_Q(4)
#define   Td_Tf(g1,g2) Tf_Tf(g1,g2,4)
#define   Td_TfFlat(i) Tf_TfFlat(i,4)

#define OffsetASfMix 1047
#define LengthASfMix 144
#define BlockASfMix(i) SlhaData(1047+i)
#define ASfMix_UASf(s1,s2,t) Slhadata(1005+s1+6*(s2)+36*(t))
#define ASfMix_UASfFlat(i,t) Slhadata(1011+i+36*(t))

#define OffsetSnuMix 1047
#define LengthSnuMix 36
#define BlockSnuMix(i) SlhaData(1047+i)
#define   SnuMix_UASf(s1,s2) ASfMix_UASf(s1,s2,1)
#define   SnuMix_UASfFlat(i) ASfMix_UASfFlat(i,1)

#define OffsetSelMix 1083
#define LengthSelMix 36
#define BlockSelMix(i) SlhaData(1083+i)
#define   SelMix_UASf(s1,s2) ASfMix_UASf(s1,s2,2)
#define   SelMix_UASfFlat(i) ASfMix_UASfFlat(i,2)

#define OffsetUSqMix 1119
#define LengthUSqMix 36
#define BlockUSqMix(i) SlhaData(1119+i)
#define   USqMix_UASf(s1,s2) ASfMix_UASf(s1,s2,3)
#define   USqMix_UASfFlat(i) ASfMix_UASfFlat(i,3)

#define OffsetDSqMix 1155
#define LengthDSqMix 36
#define BlockDSqMix(i) SlhaData(1155+i)
#define   DSqMix_UASf(s1,s2) ASfMix_UASf(s1,s2,4)
#define   DSqMix_UASfFlat(i) ASfMix_UASfFlat(i,4)

#define OffsetSnSMix 1191
#define LengthSnSMix 9
#define BlockSnSMix(i) SlhaData(1191+i)
#define SnSMix_US(g1,g2) Slhadata(1188+g1+3*(g2))
#define SnSMix_USFlat(i) Slhadata(1191+i)

#define OffsetSnAMix 1200
#define LengthSnAMix 9
#define BlockSnAMix(i) SlhaData(1200+i)
#define SnAMix_UA(g1,g2) Slhadata(1197+g1+3*(g2))
#define SnAMix_UAFlat(i) Slhadata(1200+i)

#define OffsetCVHMix 1209
#define LengthCVHMix 16
#define BlockCVHMix(i) SlhaData(1209+i)
#define CVHMix_UH(h1,h2) Slhadata(1205+h1+4*(h2))
#define CVHMix_UHFlat(i) Slhadata(1209+i)

#define OffsetNMNMix 1225
#define LengthNMNMix 25
#define BlockNMNMix(i) SlhaData(1225+i)
#define NMNMix_ZNeu(n1,n2) Slhadata(1220+n1+5*(n2))
#define NMNMix_ZNeuFlat(i) Slhadata(1225+i)

#define OffsetNMHMix 1250
#define LengthNMHMix 9
#define BlockNMHMix(i) SlhaData(1250+i)
#define NMHMix_UH(h1,h2) Slhadata(1247+h1+3*(h2))
#define NMHMix_UHFlat(i) Slhadata(1250+i)

#define OffsetNMAMix 1259
#define LengthNMAMix 9
#define BlockNMAMix(i) SlhaData(1259+i)
#define NMAMix_UA(h1,h2) Slhadata(1256+h1+3*(h2))
#define NMAMix_UAFlat(i) Slhadata(1259+i)

#define OffsetPrecObs 1268
#define LengthPrecObs 13
#define BlockPrecObs(i) SlhaData(1268+i)
#define PrecObs_DeltaRho Slhadata(1269)
#define PrecObs_MWMSSM Slhadata(1270)
#define PrecObs_MWSM Slhadata(1271)
#define PrecObs_SW2effMSSM Slhadata(1272)
#define PrecObs_SW2effSM Slhadata(1273)
#define PrecObs_gminus2mu Slhadata(1274)
#define PrecObs_EDMeTh Slhadata(1275)
#define PrecObs_EDMn Slhadata(1276)
#define PrecObs_EDMHg Slhadata(1277)
#define PrecObs_bsgammaMSSM Slhadata(1278)
#define PrecObs_bsgammaSM Slhadata(1279)
#define PrecObs_DeltaMsMSSM Slhadata(1280)
#define PrecObs_DeltaMsSM Slhadata(1281)

#define OffsetSPInfo 1281
#define LengthSPInfo 92
#define BlockSPInfo(i) SlhaData(1281+i)
#define SPInfo_NLines SlhaData(1282)
#define SPInfo_Severity SlhaData(1283)
#define SPInfo_Code(n) SlhaData(1283+n)
#define SPInfo_Text(i,n) SlhaData(1293+i+5*(n))
#define SPInfo_TextFlat(i) SlhaData(1298+i)
#define   SPInfo_Len 80

#define OffsetDCInfo 1373
#define LengthDCInfo 92
#define BlockDCInfo(i) SlhaData(1373+i)
#define DCInfo_NLines SlhaData(1374)
#define DCInfo_Severity SlhaData(1375)
#define DCInfo_Code(n) SlhaData(1375+n)
#define DCInfo_Text(i,n) SlhaData(1385+i+5*(n))
#define DCInfo_TextFlat(i) SlhaData(1390+i)
#define   DCInfo_Len 80

#define OffsetDecays 1465
#define LengthDecays 4096
#define BlockDecays(i) SlhaData(1465+i)
#define Decays_Data(n) Slhadata(1465+n)

#define nslhadata 5561

#endif
