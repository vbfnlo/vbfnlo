      function virt_tri_box(q1s,q2s,pq1,pq2,ppp,ppq1,ppq2,q1q2)
      implicit none
      complex * 16 virt_tri_box
      real * 8 q1s,q2s,pq1,pq2,ppp,ppq1,ppq2,q1q2
      complex * 16 E0_12345, Eij_12345(4,46)
      complex * 16 D0_2345,D0_1345,D0_1245,D0_1235,D0_1234
      complex * 16 Dij_2345(3,13),Dij_1345(3,13),Dij_1245(3,13)
      complex * 16 Dij_1235(3,13),Dij_1234(3,13)
      complex * 16 C0_345,C0_245,C0_235,C0_234,C0_145
      complex * 16 C0_135,C0_134,C0_125,C0_124,C0_123
      complex * 16 Cij_345(2,4),Cij_245(2,4),Cij_235(2,4),Cij_234(2,4)
      complex * 16 Cij_145(2,4),Cij_135(2,4),Cij_134(2,4),Cij_125(2,4)
      complex * 16 Cij_124(2,4),Cij_123(2,4)
      complex * 16 B0_13,B0_14,B0_24,B0_34,B0_23,B0_45,B0_35
      common/E_functions/E0_12345, Eij_12345
      common/D_functions/D0_2345,D0_1345,D0_1245,D0_1235,D0_1234,
     #     Dij_2345,Dij_1345,Dij_1245,Dij_1235,Dij_1234
      common/C_functions/
     #     C0_345,C0_245,C0_235,C0_234,C0_145,
     #     C0_135,C0_134,C0_125,C0_124,C0_123,
     #     Cij_345,Cij_245,Cij_235,Cij_234,Cij_145,
     #     Cij_135,Cij_134,Cij_125,Cij_124,Cij_123
      common/B_functions/B0_13,B0_14,B0_24,B0_34,B0_23,B0_45,B0_35
      complex * 16 dpmu1,dpmu2,dpmu3,dqmu1,dqmu2,dqmu3
      complex * 16 dlmu1,dlmu2,dlmu3,dppmu1,dppmu2,dppmu3
      complex * 16 dmu1mu2,dmu1mu3,dmu2mu3
      common/dotp_pol_vec/dpmu1,dpmu2,dpmu3,dqmu1,dqmu2,dqmu3,
     #     dlmu1,dlmu2,dlmu3,dppmu1,dppmu2,dppmu3,
     #     dmu1mu2,dmu1mu3,dmu2mu3
      complex * 16 gp_mu3_l_mu2_q_mu1,gp_l,gp_mu2,gp_mu3_l_q
      complex * 16 gp_mu2_q_mu1,gp_mu3_mu2_q,gp_q,gp_mu3_l_mu2
      complex * 16 gp_l_mu2_q,gp_mu3_mu2_mu1,gp_mu3_l_mu1,gp_mu3
      complex * 16 gp_mu1,gp_l_q_mu1,gp_mu3_q_mu1,gp_l_mu2_mu1
      common/gammaprods/gp_mu3_l_mu2_q_mu1,gp_l,gp_mu2,gp_mu3_l_q,
     #     gp_mu2_q_mu1,gp_mu3_mu2_q,gp_q,gp_mu3_l_mu2,
     #     gp_l_mu2_q,gp_mu3_mu2_mu1,gp_mu3_l_mu1,gp_mu3,
     #     gp_mu1,gp_l_q_mu1,gp_mu3_q_mu1,gp_l_mu2_mu1
      complex * 16 t1
      complex * 16 t2
      complex * 16 t3
      complex * 16 t4
      complex * 16 t6
      complex * 16 t7
      complex * 16 t8
      complex * 16 t9
      complex * 16 t12
      complex * 16 t13
      complex * 16 t14
      complex * 16 t15
      complex * 16 t16
      complex * 16 t18
      complex * 16 t19
      complex * 16 t20
      complex * 16 t21
      complex * 16 t24
      complex * 16 t25
      complex * 16 t28
      complex * 16 t33
      complex * 16 t34
      complex * 16 t46
      complex * 16 t48
      complex * 16 t50
      complex * 16 t51
      complex * 16 t52
      complex * 16 t53
      complex * 16 t54
      complex * 16 t55
      complex * 16 t60
      complex * 16 t61
      complex * 16 t62
      complex * 16 t63
      complex * 16 t65
      complex * 16 t66
      complex * 16 t68
      complex * 16 t69
      complex * 16 t75
      complex * 16 t82
      complex * 16 t86
      complex * 16 t90
      complex * 16 t91
      complex * 16 t92
      complex * 16 t95
      complex * 16 t96
      complex * 16 t99
      complex * 16 t100
      complex * 16 t105
      complex * 16 t108
      complex * 16 t109
      complex * 16 t110
      complex * 16 t111
      complex * 16 t112
      complex * 16 t113
      complex * 16 t114
      complex * 16 t115
      complex * 16 t116
      complex * 16 t118
      complex * 16 t119
      complex * 16 t121
      complex * 16 t122
      complex * 16 t124
      complex * 16 t125
      complex * 16 t126
      complex * 16 t128
      complex * 16 t129
      complex * 16 t130
      complex * 16 t131
      complex * 16 t135
      complex * 16 t138
      complex * 16 t139
      complex * 16 t143
      complex * 16 t144
      complex * 16 t146
      complex * 16 t148
      complex * 16 t149
      complex * 16 t151
      complex * 16 t152
      complex * 16 t157
      complex * 16 t160
      complex * 16 t161
      complex * 16 t162
      complex * 16 t163
      complex * 16 t164
      complex * 16 t165
      complex * 16 t166
      complex * 16 t167
      complex * 16 t168
      complex * 16 t169
      complex * 16 t170
      complex * 16 t171
      complex * 16 t172
      complex * 16 t173
      complex * 16 t174
      complex * 16 t175
      complex * 16 t176
      complex * 16 t179
      complex * 16 t180
      complex * 16 t181
      complex * 16 t182
      complex * 16 t183
      complex * 16 t185
      complex * 16 t186
      complex * 16 t187
      complex * 16 t188
      complex * 16 t189
      complex * 16 t190
      complex * 16 t191
      complex * 16 t192
      complex * 16 t194
      complex * 16 t195
      complex * 16 t196
      complex * 16 t197
      complex * 16 t198
      complex * 16 t199
      complex * 16 t200
      complex * 16 t201
      complex * 16 t202
      complex * 16 t203
      complex * 16 t205
      complex * 16 t209
      complex * 16 t210
      complex * 16 t212
      complex * 16 t214
      complex * 16 t215
      complex * 16 t216
      complex * 16 t217
      complex * 16 t219
      complex * 16 t220
      complex * 16 t221
      complex * 16 t223
      complex * 16 t225
      complex * 16 t226
      complex * 16 t228
      complex * 16 t229
      complex * 16 t230
      complex * 16 t231
      complex * 16 t239
      complex * 16 t245
      complex * 16 t253
      complex * 16 t256
      complex * 16 t258
      complex * 16 t260
      complex * 16 t261
      complex * 16 t262
      complex * 16 t264
      complex * 16 t266
      complex * 16 t267
      complex * 16 t268
      complex * 16 t270
      complex * 16 t273
      complex * 16 t277
      complex * 16 t278
      complex * 16 t279
      complex * 16 t280
      complex * 16 t281
      complex * 16 t285
      complex * 16 t288
      complex * 16 t299
      complex * 16 t300
      complex * 16 t301
      complex * 16 t303
      complex * 16 t307
      complex * 16 t309
      complex * 16 t313
      complex * 16 t314
      complex * 16 t320
      complex * 16 t323
      complex * 16 t330
      complex * 16 t331
      complex * 16 t332
      complex * 16 t334
      complex * 16 t336
      complex * 16 t337
      complex * 16 t340
      complex * 16 t342
      complex * 16 t358
      complex * 16 t361
      complex * 16 t362
      complex * 16 t363
      complex * 16 t364
      complex * 16 t365
      complex * 16 t366
      complex * 16 t367
      complex * 16 t368
      complex * 16 t369
      complex * 16 t370
      complex * 16 t371
      complex * 16 t374
      complex * 16 t375
      complex * 16 t377
      complex * 16 t379
      complex * 16 t380
      complex * 16 t381
      complex * 16 t383
      complex * 16 t385
      complex * 16 t386
      complex * 16 t388
      complex * 16 t390
      complex * 16 t395
      complex * 16 t396
      complex * 16 t397
      complex * 16 t398
      complex * 16 t405
      complex * 16 t406
      complex * 16 t407
      complex * 16 t408
      complex * 16 t409
      complex * 16 t415
      complex * 16 t419
      complex * 16 t428
      complex * 16 t429
      complex * 16 t430
      complex * 16 t434
      complex * 16 t447
      complex * 16 t450
      complex * 16 t451
      complex * 16 t456
      complex * 16 t459
      complex * 16 t460
      complex * 16 t461
      complex * 16 t463
      complex * 16 t474
      complex * 16 t481
      complex * 16 t483
      complex * 16 t496
      complex * 16 t501
      complex * 16 t502
      complex * 16 t503
      complex * 16 t504
      complex * 16 t505
      complex * 16 t506
      complex * 16 t507
      complex * 16 t508
      complex * 16 t509
      complex * 16 t510
      complex * 16 t511
      complex * 16 t512
      complex * 16 t513
      complex * 16 t514
      complex * 16 t515
      complex * 16 t516
      complex * 16 t518
      complex * 16 t519
      complex * 16 t520
      complex * 16 t521
      complex * 16 t522
      complex * 16 t523
      complex * 16 t524
      complex * 16 t525
      complex * 16 t526
      complex * 16 t527
      complex * 16 t528
      complex * 16 t529
      complex * 16 t530
      complex * 16 t531
      complex * 16 t532
      complex * 16 t533
      complex * 16 t534
      complex * 16 t535
      complex * 16 t536
      complex * 16 t537
      complex * 16 t538
      complex * 16 t557
      complex * 16 t566
      complex * 16 t568
      complex * 16 t572
      complex * 16 t574
      complex * 16 t575
      complex * 16 t576
      complex * 16 t592
      complex * 16 t593
      complex * 16 t594
      complex * 16 t597
      complex * 16 t598
      complex * 16 t599
      complex * 16 t600
      complex * 16 t601
      complex * 16 t602
      complex * 16 t605
      complex * 16 t608
      complex * 16 t609
      complex * 16 t612
      complex * 16 t614
      complex * 16 t615
      complex * 16 t616
      complex * 16 t618
      complex * 16 t620
      complex * 16 t622
      complex * 16 t623
      complex * 16 t624
      complex * 16 t634
      complex * 16 t635
      complex * 16 t648
      complex * 16 t660
      complex * 16 t663
      complex * 16 t666
      complex * 16 t667
      complex * 16 t668
      complex * 16 t670
      complex * 16 t671
      complex * 16 t672
      complex * 16 t682
      complex * 16 t683
      complex * 16 t684
      complex * 16 t689
      complex * 16 t703
      complex * 16 t706
      complex * 16 t707
      complex * 16 t708
      complex * 16 t709
      complex * 16 t710
      complex * 16 t711
      complex * 16 t712
      complex * 16 t714
      complex * 16 t715
      complex * 16 t716
      complex * 16 t720
      complex * 16 t728
      complex * 16 t733
      complex * 16 t734
      complex * 16 t735
      complex * 16 t741
      complex * 16 t748
      complex * 16 t751
      complex * 16 t753
      complex * 16 t755
      complex * 16 t756
      complex * 16 t758
      complex * 16 t759
      complex * 16 t760
      complex * 16 t761
      complex * 16 t763
      complex * 16 t766
      complex * 16 t770
      complex * 16 t771
      complex * 16 t772
      complex * 16 t773
      complex * 16 t786
      complex * 16 t787
      complex * 16 t801
      complex * 16 t804
      complex * 16 t809
      complex * 16 t812
      complex * 16 t814
      complex * 16 t820
      complex * 16 t821
      complex * 16 t822
      complex * 16 t824
      complex * 16 t830
      complex * 16 t831
      complex * 16 t832
      complex * 16 t833
      complex * 16 t834
      complex * 16 t835
      complex * 16 t836
      complex * 16 t837
      complex * 16 t838
      complex * 16 t839
      complex * 16 t847
      complex * 16 t851
      complex * 16 t865
      complex * 16 t866
      complex * 16 t875
      complex * 16 t886
      complex * 16 t893
      complex * 16 t902
      complex * 16 t904
      complex * 16 t906
      complex * 16 t907
      complex * 16 t908
      complex * 16 t909
      complex * 16 t911
      complex * 16 t913
      complex * 16 t915
      complex * 16 t916
      complex * 16 t918
      complex * 16 t919
      complex * 16 t923
      complex * 16 t926
      complex * 16 t927
      complex * 16 t928
      complex * 16 t929
      complex * 16 t930
      complex * 16 t931
      complex * 16 t932
      complex * 16 t933
      complex * 16 t934
      complex * 16 t935
      complex * 16 t936
      complex * 16 t937
      complex * 16 t938
      complex * 16 t939
      complex * 16 t946
      complex * 16 t955
      complex * 16 t967
      complex * 16 t972
      complex * 16 t977
      complex * 16 t981
      complex * 16 t982
      complex * 16 t987
      complex * 16 t988
      complex * 16 t990
      complex * 16 t991
      complex * 16 t998
      complex * 16 t999
      complex * 16 t1006
      complex * 16 t1011
      complex * 16 t1014
      complex * 16 t1020
      complex * 16 t1023
      complex * 16 t1034
      complex * 16 t1036
      complex * 16 t1039
      complex * 16 t1040
      complex * 16 t1041
      complex * 16 t1042
      complex * 16 t1043
      complex * 16 t1044
      complex * 16 t1045
      complex * 16 t1046
      complex * 16 t1047
      complex * 16 t1048
      complex * 16 t1049
      complex * 16 t1050
      complex * 16 t1051
      complex * 16 t1052
      complex * 16 t1053
      complex * 16 t1054
      complex * 16 t1055
      complex * 16 t1056
      complex * 16 t1058
      complex * 16 t1060
      complex * 16 t1061
      complex * 16 t1062
      complex * 16 t1063
      complex * 16 t1064
      complex * 16 t1065
      complex * 16 t1066
      complex * 16 t1067
      complex * 16 t1077
      complex * 16 t1078
      complex * 16 t1081
      complex * 16 t1093
      complex * 16 t1094
      complex * 16 t1095
      complex * 16 t1096
      complex * 16 t1097
      complex * 16 t1105
      complex * 16 t1110
      complex * 16 t1111
      complex * 16 t1112
      complex * 16 t1113
      complex * 16 t1114
      complex * 16 t1115
      complex * 16 t1116
      complex * 16 t1117
      complex * 16 t1118
      complex * 16 t1119
      complex * 16 t1120
      complex * 16 t1121
      complex * 16 t1122
      complex * 16 t1123
      complex * 16 t1124
      complex * 16 t1125
      complex * 16 t1126
      complex * 16 t1127
      complex * 16 t1128
      complex * 16 t1130
      complex * 16 t1131
      complex * 16 t1133
      complex * 16 t1135
      complex * 16 t1137
      complex * 16 t1138
      complex * 16 t1139
      complex * 16 t1140
      complex * 16 t1149
      complex * 16 t1150
      complex * 16 t1151
      complex * 16 t1152
      complex * 16 t1154
      complex * 16 t1161
      complex * 16 t1164
      complex * 16 t1166
      complex * 16 t1175
      complex * 16 t1177
      complex * 16 t1178
      complex * 16 t1183
      complex * 16 t1184
      complex * 16 t1194
      complex * 16 t1196
      complex * 16 t1200
      complex * 16 t1205
      complex * 16 t1209
      complex * 16 t1226
      complex * 16 t1229
      complex * 16 t1236
      complex * 16 t1252
      complex * 16 t1254
      complex * 16 t1255
      complex * 16 t1264
      complex * 16 t1269
      complex * 16 t1275
      complex * 16 t1276
      complex * 16 t1279
      complex * 16 t1286
      complex * 16 t1305
      complex * 16 t1316
      complex * 16 t1324
      complex * 16 t1325
      complex * 16 t1326
      complex * 16 t1327
      complex * 16 t1328
      complex * 16 t1329
      complex * 16 t1334
      complex * 16 t1336
      complex * 16 t1338
      complex * 16 t1344
      complex * 16 t1382
      complex * 16 t1383
      complex * 16 t1392
      complex * 16 t1394
      complex * 16 t1396
      complex * 16 t1398
      complex * 16 t1405
      complex * 16 t1408
      complex * 16 t1424
      complex * 16 t1426
      complex * 16 t1429
      complex * 16 t1430
      complex * 16 t1431
      complex * 16 t1433
      complex * 16 t1434
      complex * 16 t1435
      complex * 16 t1437
      complex * 16 t1438
      complex * 16 t1439
      complex * 16 t1441
      complex * 16 t1442
      complex * 16 t1443
      complex * 16 t1445
      complex * 16 t1446
      complex * 16 t1449
      complex * 16 t1450
      complex * 16 t1451
      complex * 16 t1456
      complex * 16 t1458
      complex * 16 t1460
      complex * 16 t1461
      complex * 16 t1467
      complex * 16 t1482
      complex * 16 t1486
      complex * 16 t1488
      complex * 16 t1491
      complex * 16 t1494
      complex * 16 t1495
      complex * 16 t1506
      complex * 16 t1518
      complex * 16 t1536
      complex * 16 t1537
      complex * 16 t1541
      complex * 16 t1542
      complex * 16 t1546
      complex * 16 t1557
      complex * 16 t1559
      complex * 16 t1560
      complex * 16 t1567
      complex * 16 t1569
      complex * 16 t1573
      complex * 16 t1574
      complex * 16 t1575
      complex * 16 t1593
      complex * 16 t1597
      complex * 16 t1602
      complex * 16 t1619
      complex * 16 t1622
      complex * 16 t1648
      complex * 16 t1651
      complex * 16 t1658
      complex * 16 t1660
      complex * 16 t1675
      complex * 16 t1695
      complex * 16 t1696
      complex * 16 t1699
      complex * 16 t1700
      complex * 16 t1703
      complex * 16 t1705
      complex * 16 t1706
      complex * 16 t1707
      complex * 16 t1711
      complex * 16 t1724
      complex * 16 t1741
      complex * 16 t1751
      complex * 16 t1761
      complex * 16 t1765
      complex * 16 t1770
      complex * 16 t1790
      complex * 16 t1815
      complex * 16 t1864
      complex * 16 t1866
      complex * 16 t1879
      complex * 16 t1883
      complex * 16 t1890
      complex * 16 t1898
      complex * 16 t1931
      complex * 16 t1935
      complex * 16 t1942
      complex * 16 t1955
      complex * 16 t1978
      complex * 16 t1981
      complex * 16 t1989
      complex * 16 t2008
      complex * 16 t2022
      complex * 16 t2053
      complex * 16 s1
      complex * 16 s2
      complex * 16 s3
      complex * 16 s4
      t1 = 2*dqmu3
      t2 = 2*dpmu3
      t3 = 2*dlmu3
      t4 = t1-dppmu3+t2+t3
      t6 = 2*pq1
      t7 = q1s+t6
      t8 = 1/t7
      t9 = t8*gp_mu2
      t12 = 2*q1q2
      t13 = 2*pq2
      t14 = t6+q1s-ppq2-ppp-ppq1+t12+q2s+t13
      t15 = t14*dpmu1
      t16 = dqmu2+dpmu2
      t18 = q1s+t12+t6+q2s+t13
      t19 = 1/t18
      t20 = t19*t8
      t21 = t20*gp_mu3
      t24 = t14*t19
      t25 = t8*gp_mu3_l_mu2_q_mu1
      t28 = t20*gp_mu3_l_mu2
      t33 = t14*t16
      t34 = t20*gp_mu3_q_mu1
      t46 = q1s*dmu1mu2
      t48 = pq1*dmu1mu2
      t50 = dqmu1*dpmu2
      t51 = 4*t50
      t52 = dqmu1*dqmu2
      t53 = 4*t52
      t54 = dpmu1*dqmu2
      t55 = 4*t54
      t60 = t19*gp_mu3_l_mu2_q_mu1
      t61 = 12*t60
      t62 = t19*dqmu1
      t63 = t62*gp_mu3_l_mu2
      t65 = t7*t19
      t66 = t65*gp_mu3_mu2_mu1
      t68 = 6*dpmu2
      t69 = 4*dqmu2
      t75 = 6*dpmu1
      t82 = 16*t19*dmu1mu2*gp_mu3_l_q
      t86 = dqmu3+dlmu3+dpmu3
      t90 = dpmu1*t16
      t91 = t8*gp_mu3
      t92 = t90*t91
      t95 = t8*dpmu1
      t96 = t95*gp_mu3_l_mu2
      t99 = t16*t8
      t100 = t99*gp_mu3_q_mu1
      t105 = 8*t86*dpmu1*t9-8*t92-2*t25-4*t96+2*gp_mu3_mu2_mu1-4*t100+4*
     #t86*t8*gp_mu2_q_mu1
      t108 = q1s*dpmu2
      t109 = q1q2*dpmu2
      t110 = 2*t109
      t111 = q2s*dpmu2
      t112 = pq1*dqmu2
      t113 = 2*t112
      t114 = pq2*dqmu2
      t115 = 2*t114
      t116 = t108+t110+t111-t113-t115
      t118 = t19*gp_mu3
      t119 = dpmu1*t116*t118
      t121 = pq2+pq1
      t122 = t121*t19
      t124 = 4*t122*gp_mu3_l_mu2_q_mu1
      t125 = gp_mu3_l_mu2*dpmu1
      t126 = t122*t125
      t128 = pq2*q1s
      t129 = pq1*q1q2
      t130 = 2*t129
      t131 = pq1*q2s
      t135 = 4*(-t128+t130+t131)*t19*gp_mu3_mu2_mu1
      t138 = 4*t116*t19*gp_mu3_q_mu1
      t139 = dpmu1*gp_mu3_mu2_q
      t143 = q2s*dpmu3
      t144 = dqmu3*q2s
      t146 = 2*pq2*dlmu3
      t148 = 2*q1q2*dlmu3
      t149 = dppmu3*q2s
      t151 = 2*ppq2*dppmu3
      t152 = -t143-t144+t146+t148+t149-t151
      t157 = 4*dlmu2+4*dpmu2+4*dqmu2-4*dppmu2
      t160 = 4*dlmu3*gp_l_mu2_mu1
      t161 = q1q2*dmu2mu3
      t162 = 8*t161
      t163 = dlmu2*dlmu3
      t164 = 16*t163
      t165 = dlmu3*dpmu2
      t166 = 8*t165
      t167 = dlmu3*dqmu2
      t168 = 8*t167
      t169 = pq2*dmu2mu3
      t170 = 8*t169
      t171 = dmu2mu3*q2s
      t172 = 8*t171
      t173 = dlmu2*dqmu3
      t174 = 8*t173
      t175 = dlmu2*dpmu3
      t176 = 8*t175
      t179 = dppmu2*dqmu3
      t180 = dppmu2*dpmu3
      t181 = dppmu2*dppmu3
      t182 = dpmu2*dpmu3
      t183 = 2*t182
      t185 = dpmu2*dqmu3
      t186 = 2*t185
      t187 = 2*t173
      t188 = 2*t163
      t189 = dlmu2*dppmu3
      t190 = 2*t175
      t191 = dqmu2*dpmu3
      t192 = 2*t191
      t194 = dqmu3*dqmu2
      t195 = 2*t194
      t196 = ppq1*dmu2mu3
      t197 = dmu2mu3*q1s
      t198 = dmu2mu3*ppp
      t199 = pq1*dmu2mu3
      t200 = 2*t199
      t201 = 2*t161
      t202 = 2*t169
      t203 = t179+t180+t181-t183-3*t165-t186-t187-t188-t189-t190-t192-3*
     #t167-t195-t196+t171+t197-t198+t200+t201+t202
      t205 = t8*gp_l
      t209 = 2*q1q2*dppmu2
      t210 = q2s*dppmu2
      t212 = 2*pq2*dppmu2
      t214 = 4*t109
      t215 = dpmu2*pq2
      t216 = 4*t215
      t217 = q1s*dlmu2
      t219 = 2*ppp*dlmu2
      t220 = dlmu2*pq1
      t221 = 2*t220
      t223 = 2*dlmu2*ppq1
      t225 = 2*dlmu2*ppq2
      t226 = q2s*dqmu2
      t228 = dqmu2*q1q2
      t229 = 4*t228
      t230 = 4*t114
      t231 = -t209-t210-t212+3*t111+t214+t216-t217+t219-t221+t223+t225+3
     #*t226+t229+t230
      t239 = (-ppq1+q2s+q1s-ppp+t6+t12+t13)*t8
      t245 = 4*pq2+4*q2s+4*q1q2
      t253 = -8*t152*dpmu1*t9+t157*gp_mu3_l_mu1-t160+(-t162+t164+t166+t1
     #68-t170-t172+t174+t176)*gp_mu1-16*dpmu1*t203*t205-8*dpmu1*t231*t91
     #-8*t203*t8*gp_l_q_mu1-4*t239*gp_mu3_l_mu2_q_mu1-8*t239*t125+t245*g
     #p_mu3_mu2_mu1-4*t231*t8*gp_mu3_q_mu1-4*t152*t8*gp_mu2_q_mu1
      t256 = t90*t21
      t258 = t20*gp_mu3_l_mu2_q_mu1
      t260 = dpmu1*t19
      t261 = t8*gp_mu3_l_mu2
      t262 = t260*t261
      t264 = t19*gp_mu3_mu2_mu1
      t266 = t16*t19
      t267 = t8*gp_mu3_q_mu1
      t268 = t266*t267
      t270 = -8*t256-2*t258-4*t262+2*t264-4*t268
      t273 = t19*dlmu2
      t277 = 2*q1q2*dlmu2
      t278 = dlmu2*q2s
      t279 = pq2*dlmu2
      t280 = 2*t279
      t281 = -t111+t277+t278+t280-t226
      t285 = q2s*t19
      t288 = q2s*dpmu1
      t299 = 8*t163
      t300 = 4*t171
      t301 = -t299+t300
      t303 = -t165-t175-t173-t167+t169+t161
      t307 = t99*gp_mu3
      t309 = 8*t288*t307
      t313 = pq2+q1q2
      t314 = t313*t8
      t320 = 2*q2s*gp_mu3_mu2_mu1
      t323 = 4*q2s*t16*t267
      t330 = 24*t96
      t331 = dpmu1*dlmu2
      t332 = t331*t91
      t334 = dmu2mu3*dpmu1
      t336 = 32*t334*t205
      t337 = 12*t25
      t340 = 16*dmu2mu3*t8*gp_l_q_mu1
      t342 = dlmu2*t8*gp_mu3_q_mu1
      t358 = dpmu3+dqmu3-dppmu3
      t361 = 8*t182
      t362 = dpmu2*dppmu3
      t363 = 16*t362
      t364 = dqmu2*dppmu3
      t365 = 16*t364
      t366 = 4*t197
      t367 = 8*t199
      t368 = 8*t194
      t369 = 8*t185
      t370 = 8*t191
      t371 = -t168-t361+t162-t363-t365-t166+t366+t367+t170-t368-t369-t17
     #4-t370-t176
      t374 = 8*t334*gp_l
      t375 = q1s*dppmu2
      t377 = pq1*dppmu2
      t379 = dpmu2*pq1
      t380 = 2*t379
      t381 = ppp*dpmu2
      t383 = ppq1*dpmu2
      t385 = dqmu2*q1s
      t386 = ppp*dqmu2
      t388 = ppq1*dqmu2
      t390 = -2*t375-4*t377+t108+t380+4*t381+4*t383+t217+t221+t385+t113+
     #4*t386+4*t388
      t395 = 4*dmu2mu3*gp_l_q_mu1
      t396 = 2*ppq1
      t397 = 2*ppp
      t398 = q1s+t6-t396-t397
      t405 = 4*pq1
      t406 = 4*pq2
      t407 = 4*q1q2
      t408 = 2*q1s
      t409 = 4*ppq2
      t415 = 4*dlmu3
      t419 = 8*(dlmu3+2*dppmu3)*dpmu1*gp_mu2+(-4*dppmu2+4*dqmu2+4*dpmu2)
     #*gp_mu3_l_mu1+4*t358*gp_l_mu2_mu1+t371*gp_mu1-t374+8*dpmu1*t390*t9
     #1-t395-2*t398*t8*gp_mu3_l_mu2_q_mu1-4*t398*dpmu1*t261+(-t405-t406-
     #t407-t408+t409)*gp_mu3_mu2_mu1+4*t390*t8*gp_mu3_q_mu1+(t415+8*dppm
     #u3)*gp_mu2_q_mu1
      t428 = 4*t266*gp_mu3_l_mu1
      t429 = t12+q2s+t13
      t430 = t429*t19
      t434 = t430*gp_mu3_mu2_mu1
      t447 = 16*t256+4*t258+8*t262-4*t264+8*t268
      t450 = 2*dqmu2
      t451 = 2*dpmu2
      t456 = t226+t217+t221+t111
      t459 = 8*dpmu1*t456*t21
      t460 = 6*pq1
      t461 = 3*q1s
      t463 = (t460+t461+t13+t12)*t19
      t474 = 4*t456*t19*t267
      t481 = 2*t60
      t483 = 2*t66
      t496 = 2*t215
      t501 = pq2*dmu1mu2
      t502 = t501*q1q2
      t503 = 2*t502
      t504 = t501*q2s
      t505 = pq2**2
      t506 = dmu1mu2*t505
      t507 = 2*t506
      t508 = t48*q1q2
      t509 = 2*t508
      t510 = t48*q2s
      t511 = t48*pq2
      t512 = 2*t511
      t513 = q1q2*dqmu1
      t514 = t513*dpmu2
      t515 = 2*t514
      t516 = t50*q2s
      t518 = 2*t50*pq2
      t519 = t288*dpmu2
      t520 = q1q2*dpmu1
      t521 = t520*dpmu2
      t522 = 4*t521
      t523 = dlmu1*dpmu2
      t524 = t523*q1q2
      t525 = 2*t524
      t526 = dlmu1*q2s
      t527 = t526*dpmu2
      t528 = t523*pq2
      t529 = 2*t528
      t530 = q1s*dpmu1
      t531 = t530*dlmu2
      t532 = t331*q2s
      t533 = t331*q1q2
      t534 = 2*t533
      t535 = dpmu1*pq2
      t536 = t535*dqmu2
      t537 = 4*t536
      t538 = t503+t504+t507+t509+t510+t512-t515-t516-t518-t519-t522-t525
     #-t527-t529-t531-t532-t534+t537
      t557 = dpmu1*dpmu2
      s1 = (8*t4*dpmu1*t9-16*t15*t16*t21-4*t24*t25-8*t15*t28+4*t24*gp_mu
     #3_mu2_mu1-8*t33*t34+4*t4*t8*gp_mu2_q_mu1)*Cij_145(2,3)+(-8*t19*dqm
     #u2*gp_mu3_l_mu1+4*(4*t46+8*t48-t51-t53-t55)*t19*gp_mu3+t61-8*t63-1
     #2*t66+4*(t68+t69)*t19*gp_mu3_q_mu1+4*(2*dqmu1+t75)*t19*gp_mu3_mu2_
     #q+t82)*Dij_1234(3,12)+t105*Cij_145(1,2)+(-8*t119+t124+8*t126+t135-
     #t138-4*t139)*D0_1234+t253*Dij_1345(2,6)+t270*B0_14/2+(4*t273*gp_mu
     #3_l_mu1-8*dpmu1*t281*t21+2*t285*t25+4*t288*t28-2*t285*gp_mu3_mu2_m
     #u1-4*t281*t19*t267)*Cij_134(2,2)+(t301*gp_mu1+16*dpmu1*t303*t205+t
     #309+8*t303*t8*gp_l_q_mu1+4*t314*gp_mu3_l_mu2_q_mu1+8*t314*t125-t32
     #0+t323)*Dij_1345(3,6)
      s2 = s1+(-16*dpmu1*dlmu3*t9+t330-16*t332+t336+t337+t340-8*t342-8*d
     #lmu3*t8*gp_mu2_q_mu1)*Dij_1345(3,12)+t419*Dij_1345(1,1)+(-8*t119+8
     #*t126-8*t122*t139)*Dij_1234(2,1)+(t428+4*t430*t96-8*t332-2*t434+2*
     #t430*t25-4*t342)*Cij_134(1,2)
      t566 = s2+t447*Cij_145(2,4)+(-4*(dlmu2+t450+t451)*t19*gp_mu3_l_mu1
     #-t459+2*t463*t25+4*t463*t96+2*(-t405-t408+q2s)*t19*gp_mu3_mu2_mu1-
     #t474)*Cij_134(1,1)+(-8*dqmu1*t16*t118+t481-4*t63-t483+4*t266*gp_mu
     #3_q_mu1+4*(dqmu1+dpmu1)*t19*gp_mu3_mu2_q)*Cij_123(1,2)-t105*Cij_14
     #5(2,1)+(4*(t496-t110-t108+t280+t221+t230+t113)*t19*gp_mu3_l_mu1+8*
     #t538*t19*gp_mu3-t124+4*(t12+t13+q2s-t6)*t19*t125-4*t121*t429*t264+
     #4*(t111+t496+t280+t221+t115)*t19*gp_mu3_q_mu1-4*(q2s+t406+t12)*t19
     #*t139-8*(t501+t48-t523-t50-3*t557-t54)*t19*gp_mu3_l_q)*Dij_1234(2,
     #5)
      t568 = -2*t358
      t572 = dmu2mu3*gp_mu1
      t574 = 2*dppmu2
      t575 = 4*dpmu2
      t576 = -t574+t575+t69
      t592 = q1q2+pq2+q2s-ppq2
      t593 = t592*dppmu3
      t594 = t95*gp_mu2
      t597 = dlmu3*dppmu2
      t598 = 8*t597
      t599 = 8*t189
      t600 = ppq2*dmu2mu3
      t601 = 8*t600
      t602 = t598-t166+t599-t164+t162+t172-t601+t170-t174-t176-t168
      t605 = 2*t181
      t608 = 2*t362
      t609 = 4*t165
      t612 = 6*t163
      t614 = 2*t180-t605+2*t179+4*t597-t183+t608-t186-t609+4*t189-4*t175
     #-t612-4*t173
      t615 = 2*t364
      t616 = 4*t167
      t618 = 3*t171
      t620 = 2*t196
      t622 = 2*t198
      t623 = -t192+t615-t195-t616+t197+t200+4*t169+t618-4*t600-t620+4*t1
     #61-t622
      t624 = t614+t623
      t634 = 3*q2s
      t635 = q1s+t6+t406+t634-t409-t396+t407-t397
      t648 = t8*gp_mu2_q_mu1
      t660 = q1q2*dmu1mu2
      t663 = q2s*dmu1mu2
      t666 = 4*t557
      t667 = 4*t523
      t668 = 4*t331
      t670 = dlmu1*dqmu2
      t671 = 4*t670
      t672 = 4*t501+4*t660+4*t48+2*t663+2*t46-t51+t666-t667-t53+t668+8*t
     #54-t671
      t682 = 2*q2s
      t683 = 8*pq1
      t684 = 4*q1s
      t689 = 2*dlmu2
      t703 = dpmu1*dppmu3*gp_mu2
      t706 = 16*t173
      t707 = 16*t165
      t708 = 16*t169
      t709 = 16*t161
      t710 = 16*t175
      t711 = 16*t167
      t712 = -t598-t599-t172+t601+t706+t707-t708+t164-t709+t710+t711
      t714 = 2*t165
      t715 = 2*t167
      t716 = t180+t179-t186-t714+t362-t183-t190-t187-t195-t715+t364-t192
     #+t197-t196+t202+t201+t200-t198
      t720 = q2s+t13+t12-ppq2
      t728 = (q1s-ppq1+t13+t12+t6-ppp)*t8
      t733 = 8*q1q2
      t734 = 8*pq2
      t735 = 4*q2s
      t741 = dppmu3*gp_mu2_q_mu1
      t748 = (-ppp+pq2+q1q2+q1s-ppq1+t6)*t8
      t751 = 24*t182
      t753 = 24*t191
      t755 = 24*t199
      t756 = -t751+t709-t707-16*t196-t711-t753+t363-t706+t365+t599+16*t1
     #80+t755
      t758 = 8*t181
      t759 = 24*t185
      t760 = 12*t197
      t761 = 24*t194
      t763 = -16*t198-t299-t601-t710+t598-t758-t759+t760+t708+t300-t761+
     #16*t179
      t766 = t180+t179-t183-t186-t165+t362-t175-t173-t192-t195-t167+t364
     #-t198+t169+t161+t197-t196+t200
      t770 = 4*ppp
      t771 = 2*ppq2
      t772 = 4*ppq1
      t773 = t407-t770+t460+q2s+t406-t771-t772+t461
      t786 = 12*pq1
      t787 = 6*q1s
      t801 = pq1*t19
      t804 = q1s*t19
      t809 = pq1*t16
      t812 = t260*gp_mu3_mu2_q
      t814 = -8*t530*t16*t21+4*t801*t25-4*t804*t96-4*t801*gp_mu3_mu2_mu1
     #+8*t809*t34+4*t812
      t820 = t260*dpmu2*gp_mu3_l_q
      t821 = dpmu1*pq1
      t822 = t19*gp_mu3_l_mu2
      t824 = t19*gp_mu3_mu2_q
      t830 = t48*q1s
      t831 = pq1**2
      t832 = dmu1mu2*t831
      t833 = 2*t832
      t834 = t50*q1s
      t835 = t50*pq1
      t836 = 2*t835
      t837 = t557*q1s
      t838 = t112*dpmu1
      t839 = 2*t838
      t847 = 4*t530*t822
      t851 = t19*gp_mu3_q_mu1
      t865 = 2*t52
      t866 = t46-t865
      t875 = q1s*t16
      t886 = t19*dpmu2
      t893 = t260*gp_mu3_l_mu2
      t902 = 2*ppq2*dpmu3
      t904 = 2*ppq2*dqmu3
      t906 = 2*ppp*dlmu3
      t907 = pq1*dlmu3
      t908 = 2*t907
      t909 = q1s*dlmu3
      t911 = 2*dlmu3*ppq1
      t913 = 2*q1q2*dppmu3
      t915 = 2*pq2*dppmu3
      t916 = dppmu3*q1s
      t918 = 2*pq1*dppmu3
      t919 = t143-t902+t144-t904+t906-t148-t908-t909+t911-t146+t913+t915
     #+t151+t916+t918
      t923 = dlmu3+dpmu3-dppmu3+dqmu3
      t926 = 8*t196
      t927 = 8*t197
      t928 = 16*t199
      t929 = 8*t198
      t930 = t709-t926+t927+t928-t929+t708-t601+t172+t599+t598-t164
      t931 = 8*t179
      t932 = 8*t180
      t933 = 8*t362
      t934 = 8*t364
      t935 = 16*t182
      t936 = 16*t185
      t937 = 16*t194
      t938 = 16*t191
      t939 = t931+t932+t933+t934-t710-t711-t707-t706-t935-t936-t937-t938
      t946 = t605-4*t185-t609+t608-4*t182-t190-t187-t188-4*t194-t616+t61
     #5-4*t191+t197+t201+t200+t171+t202
      t955 = t18*t8
      t967 = -8*t919*dpmu1*t9+4*t923*gp_l_mu2_mu1+(t930+t939)*gp_mu1+8*d
     #pmu1*t946*t205+16*t15*t307+4*t946*t8*gp_l_q_mu1+2*t955*gp_mu3_l_mu
     #2_q_mu1+4*t955*t125+(t770-t683+t409-t734-t733-t684+t772-t735)*gp_m
     #u3_mu2_mu1+8*t33*t267-4*t919*t8*gp_mu2_q_mu1
      t972 = -t143-t144+t146+t148+t913+t149+t915
      t977 = -4*t16
      t981 = 2*t189
      t982 = -t608-t183-t186-t714-t981-t615-t192-t195-t715+t197+t200
      t987 = 2*ppq2*dpmu2
      t988 = 2*t228
      t990 = 2*ppq2*dqmu2
      t991 = -t209-t210-t212+t110+t111+t496+t987-t217+t219-t221+t223+t22
     #5+t988+t226+t115+t990
      t998 = 2*gp_mu3_l_mu2_q_mu1
      t999 = 4*t125
      t1006 = 8*t972*dpmu1*t9+t977*gp_mu3_l_mu1+t160+(-t299+t300-t599)*g
     #p_mu1+8*dpmu1*t982*t205+8*dpmu1*t991*t91+4*t982*t8*gp_l_q_mu1+t998
     #+t999-t320+4*t991*t8*gp_mu3_q_mu1+4*t972*t8*gp_mu2_q_mu1
      t1011 = -t183-t186-t192-t195+t197+t200
      t1014 = 8*dpmu1*t1011*t205
      t1020 = 4*t1011*t8*gp_l_q_mu1
      t1023 = -4*t313*gp_mu3_mu2_mu1
      t1034 = dqmu2+dpmu2+dlmu2
      t1036 = t19*gp_mu3_l_mu1
      t1039 = t46*q2s
      t1040 = t663*q1q2
      t1041 = 2*t1040
      t1042 = q1q2**2
      t1043 = dmu1mu2*t1042
      t1044 = 4*t1043
      t1045 = 4*t502
      t1046 = 2*t510
      t1047 = dlmu2*dqmu1
      t1048 = t1047*q1q2
      t1049 = 4*t1048
      t1050 = t278*dqmu1
      t1051 = 2*t1050
      t1052 = t279*dqmu1
      t1053 = 4*t1052
      t1054 = 4*t533
      t1055 = dlmu1*dlmu2
      t1056 = t1055*q1s
      t1058 = t1055*pq1
      t1060 = t670*q1q2
      t1061 = 4*t1060
      t1062 = t526*dqmu2
      t1063 = 2*t1062
      t1064 = dlmu1*pq2
      t1065 = t1064*dqmu2
      t1066 = 4*t1065
      t1067 = t1039+t1041+t1044+t1045+t1046-t1049-t1051-t1053+t1054-2*t1
     #056-4*t1058-t1061-t1063-t1066
      t1077 = q1s*q2s
      t1078 = q2s*q1q2
      t1081 = pq2*q1q2
      t1093 = 4*t285*t139
      t1094 = 2*t660
      t1095 = 2*t1047
      t1096 = 2*t1055
      t1097 = 2*t670
      t1105 = 2*t108
      t1110 = t52*pq2
      t1111 = 4*t1110
      t1112 = t670*q1s
      t1113 = 4*t511
      t1114 = t670*pq1
      t1115 = 2*t1114
      t1116 = 2*t1043
      t1117 = t523*pq1
      t1118 = 2*t1117
      t1119 = 2*t1065
      t1120 = 2*t1060
      t1121 = 2*t1052
      t1122 = 2*t1048
      t1123 = t52*q2s
      t1124 = t523*q1s
      t1125 = t52*q1q2
      t1126 = 4*t1125
      t1127 = -t1111-t1112+t1113-t1115+t1116-t1118-t1119-t1062-t1120-t11
     #21-t1050-t1122+t1040-t1123-t1124-t1126
      t1128 = t1047*pq1
      t1130 = t1047*q1s
      t1131 = t501*q1s
      t1133 = t660*q1s
      t1135 = 4*t508
      t1137 = t288*dqmu2
      t1138 = 2*t536
      t1139 = 2*t521
      t1140 = -2*t1128-t1130+3*t1131+3*t1133+t1135+t1039-t515-2*t1058-t1
     #056+t503-t1137-t1138+t1139+t534-t518+t510
      t1149 = q2s*dqmu1
      t1150 = 2*t821
      t1151 = 2*t530
      t1152 = q1s*dlmu1
      t1154 = 2*pq1*dlmu1
      t1161 = q1q2*q1s
      t1164 = 4*pq1*pq2
      t1166 = 4*t129
      t1175 = 2*t535
      t1177 = 2*q1q2*dlmu1
      t1178 = 2*t1064
      t1183 = 2*t523
      t1184 = 2*t54
      s3 = (8*t568*dpmu1*t9+8*t572+8*dpmu1*t576*t91+4*t25+8*t96-8*gp_mu3
     #_mu2_mu1+4*t576*t8*gp_mu3_q_mu1+4*t568*t8*gp_mu2_q_mu1)*Dij_1345(2
     #,7)+(-16*t593*t594+t602*gp_mu1+8*dpmu1*t624*t205+16*t592*dpmu1*t30
     #7+4*t624*t8*gp_l_q_mu1+2*t635*t8*gp_mu3_l_mu2_q_mu1+4*t635*dpmu1*t
     #261-4*t592*gp_mu3_mu2_mu1+8*t592*t16*t267-8*t593*t648)*Dij_1345(3,
     #9)
      s2 = s3+(4*(2*dlmu2+2*dqmu2)*t19*gp_mu3_l_mu1+4*t672*t19*gp_mu3+4*
     #t60-4*(2*dqmu1-2*dpmu1+2*dlmu1)*t19*gp_mu3_l_mu2-2*(t407+t406+t682
     #+t683+t684)*t19*gp_mu3_mu2_mu1+4*(t451+t689+t69)*t19*gp_mu3_q_mu1-
     #4*(-2*dqmu1-2*dlmu1)*t19*gp_mu3_mu2_q)*Dij_1234(2,7)+(16*t314*t703
     #+t712*gp_mu1-16*dpmu1*t716*t205-16*t720*dpmu1*t307-8*t716*t8*gp_l_
     #q_mu1-4*t728*gp_mu3_l_mu2_q_mu1-8*t728*t125+(t733-t409+t734+t735)*
     #gp_mu3_mu2_mu1-8*t720*t16*t267+8*t314*t741)*Dij_1345(3,10)
      s1 = s2+(-16*t748*t703+(t756+t763)*gp_mu1+16*dpmu1*t766*t205+8*t77
     #3*dpmu1*t307+8*t766*t8*gp_l_q_mu1+4*t748*gp_mu3_l_mu2_q_mu1+8*t748
     #*t125+(t409-t682-t734-t733+8*ppq1+8*ppp-t786-t787)*gp_mu3_mu2_mu1+
     #4*t773*t16*t267-8*t748*t741)*Dij_1345(3,7)+t814*C0_123+(8*dpmu1*(-
     #t108+t113)*t118-8*t820+8*t821*t822-8*t821*t824)*Dij_1234(3,4)+(8*(
     #t830+t833-t834-t836+t837-t839)*t19*gp_mu3+4*t801*gp_mu3_l_mu2_q_mu
     #1+t847-4*pq1*t7*t264+8*t809*t851+4*(-q1s+t6)*t19*t139+8*(-t50-t54+
     #t48)*t19*gp_mu3_l_q)*Dij_1234(3,6)+(4*t7*t866*t118+2*t804*gp_mu3_l
     #_mu2_q_mu1-2*q1s*t7*t264+4*t875*t851+4*t530*t824+4*t866*t19*gp_mu3
     #_l_q)*Dij_1234(3,2)
      s2 = s1+(-8*t886*gp_mu3_l_mu1+8*dpmu1*(t451+t69)*t118+16*t893-8*t8
     #86*gp_mu3_q_mu1-16*t812)*Dij_1234(3,11)+t967*Dij_1345(2,3)+t814*Ci
     #j_123(2,3)+t1006*Dij_1345(1,2)
      s3 = s2+(8*t303*gp_mu1+t1014+16*t313*dpmu1*t307+t1020+t998+t999+t1
     #023+8*t313*t16*t267)*Dij_1345(3,4)+(-t428+4*t893-t483+t481)*Cij_13
     #4(2,1)
      s4 = s3+(8*q1q2*t1034*t1036+4*t1067*t19*gp_mu3+2*(-t12+q2s)*t19*gp
     #_mu3_l_mu2_q_mu1-8*t520*t822-2*(t1077+2*t1078+4*t1042+4*t1081+2*t1
     #31)*t19*gp_mu3_mu2_mu1+4*(t111+t277+t226)*t19*gp_mu3_q_mu1+t1093+4
     #*(-t1094+t663+t1095-t1096+t1097)*t19*gp_mu3_l_q)*Dij_1234(3,9)
      t1194 = s4+(4*(t110+t380+t1105+t277+t217+t221+t988+t385)*t19*gp_mu
     #3_l_mu1+8*(t1127+t1140)*t19*gp_mu3-4*(-pq2+pq1+q1s)*t19*gp_mu3_l_m
     #u2_q_mu1+4*(t1149-t1150-t1151+t1152+t1154)*t19*gp_mu3_l_mu2-4*(2*t
     #1081+2*t1042+t1078+3*t1161+t1077+t1164+3*t128+t1166+t131)*t19*gp_m
     #u3_mu2_mu1+4*(t496+t110+t277+t217+t221+t230+t229+t226)*t19*gp_mu3_
     #q_mu1+4*(-t1149+t1175+t1177+t526+t1178)*t19*gp_mu3_mu2_q-8*(-t501+
     #t48+t46+t1183-t50-t865+t1055+t331-t1184+t670)*t19*gp_mu3_l_q)*Dij_
     #1234(2,6)+t447*Cij_123(2,4)
      t1196 = t143+t144-t146-t148-t915-t913
      t1200 = -dpmu2+dlmu2-dqmu2
      t1205 = -t165-t183-t362-t186-t167-t192-t364-t195+t197+t200
      t1209 = -t212-t209+t111+t214+t216+t987-t217-t221+t223+t219+t226+t2
     #29+t230+t990
      t1226 = -8*t1196*dpmu1*t9+4*t1200*gp_mu3_l_mu1+t160+(-t599-t176+t1
     #72-t164-t166-t174-t168+t162+t170)*gp_mu1+16*dpmu1*t1205*t205+8*dpm
     #u1*t1209*t91+8*t1205*t8*gp_l_q_mu1+4*gp_mu3_l_mu2_q_mu1+8*t125-t24
     #5*gp_mu3_mu2_mu1+4*t1209*t8*gp_mu3_q_mu1-4*t1196*t8*gp_mu2_q_mu1
      t1229 = dpmu1*gp_mu3
      t1236 = pq1*gp_mu3_mu2_mu1
      t1252 = 2*t1110
      t1254 = -t1252+t512+t1116-t1119-t1062-t1120-t1121-t1050-t1122+t104
     #0-2*t1125+t1131+t1133+t509
      t1255 = t507-t515-t525+t504-t527+t1045+t1137+t1138-t1139-t534-t532
     #-t531-t518-t529
      t1264 = 4*t279
      t1269 = 2*t520
      t1275 = 2*t48
      t1276 = 2*t501
      t1279 = t46+t1094+t1275+t663+t1276-2*t50-t666-t1095-t865-t55-2*t33
     #1
      t1286 = t171-t188
      t1305 = (q1q2+pq2+q1s+t6)*t19
      t1316 = 8*t703
      t1324 = 2*t381
      t1325 = 2*t383
      t1326 = 2*t377
      t1327 = 2*t386
      t1328 = 2*t388
      t1329 = -t1324-t1325+t375+t1326-t1327-t1328
      t1334 = (ppp+ppq1)*t8
      t1336 = 4*t1334*gp_mu3_l_mu2_q_mu1
      t1338 = 8*t1334*t125
      t1344 = 4*t741
      t1382 = q2s**2
      t1383 = dmu1mu2*t1382
      t1392 = 4*t1055*q1q2
      t1394 = 2*t1055*q2s
      t1396 = 4*t1055*pq2
      t1398 = -t1049-t1051-t1053-t1054-t1392-t1394-t1396+4*t1137-t1061-t
     #1063-t1066
      t1405 = t313*t19
      t1408 = t429**2
      t1424 = 4*t112
      t1426 = (t380+t217+t221+t1424+t385)*t19
      t1429 = q1s**2
      t1430 = dmu1mu2*t1429
      t1431 = 2*t1133
      t1433 = 2*t1131
      t1434 = 4*t832
      t1435 = 2*t834
      t1437 = 4*t835
      t1438 = t1430+t1431+4*t830+t1039+t1433+t1135+t1434+t1046+t1113-t14
     #35-2*t516-t1437
      t1439 = 2*t1124
      t1441 = 2*t1130
      t1442 = 4*t1128
      t1443 = t52*q1s
      t1445 = 2*t1123
      t1446 = t52*pq1
      t1449 = 2*t1112
      t1450 = 4*t1114
      t1451 = t522-t1439-4*t1117-t1441-t1442-2*t1443-t1445-4*t1446-2*t11
     #37-t537-t1449-t1450
      t1456 = dqmu1*q1s
      t1458 = 2*pq1*dqmu1
      t1460 = 2*pq2*dqmu1
      t1461 = 2*t513
      t1467 = (-t408-t405)*gp_mu3_mu2_mu1
      t1482 = dppmu3+dlmu3
      t1486 = -t977
      t1488 = dpmu3+dqmu3
      t1491 = -t176-t935-t938-t168+t928-t166-t937-t933+t927-t936-t934-t1
     #74+t170+t162
      t1494 = 2*t385
      t1495 = -t375-t1326+4*t379+t1324+t1325+t1105+t217+t221+t1424+t1327
     #+t1328+t1494
      t1506 = 8*t1482*dpmu1*gp_mu2+t1486*gp_mu3_l_mu1+4*t1488*gp_l_mu2_m
     #u1+t1491*gp_mu1-t374+8*dpmu1*t1495*t91-t395-t998-t999+(-t407-t406-
     #t683-t684)*gp_mu3_mu2_mu1+4*t1495*t8*gp_mu3_q_mu1+4*t1482*gp_mu2_q
     #_mu1
      t1518 = t1456+t1458-t530
      t1536 = 4*dppmu3
      t1537 = t2+t1+t3+t1536
      t1541 = 16*t572
      t1542 = t574+t575-t689+t69
      t1546 = 12*gp_mu3_mu2_mu1
      t1557 = 2*ppq1*dppmu3
      t1559 = 2*dppmu3*ppp
      t1560 = -t143-t144+t909+t908+t146+t148+t913+t149+t915+t916+t918-t1
     #557-t1559
      t1567 = -t162+t934+t370+t361+t369+t368-t300-t366-t367+t168-t170-t7
     #58+t174+t176+t933+t166+t299+t599
      t1569 = t182+t165+t185+t362+t189+t191+t167+t194+t364
      t1573 = -t375-t209-t1326-t210-t212+t380+t111+t1324+t496+t987+t1325
     #+t110
      t1574 = t108+t225+t219+t223+t113+t226+t1327+t115+t990+t1328+t988+t
     #385
      t1575 = t1573+t1574
      t1593 = dlmu3+dpmu3+dqmu3+dppmu3
      t1597 = t6+q2s+t397+t13+t771+t396+t12+q1s
      t1602 = t1597*t19
      t1619 = (-t108+t380+t221+t1424)*t19
      t1622 = t830+t509+t833+t510+t512-t516-t515-t836-t834+t1139-t837-t1
     #118+t1252-t531+t839-t1138-t1137-t1115
      s2 = t1226*Dij_1345(2,4)+3*(-4*t1229*dqmu2-4*t1229*dpmu2-gp_mu3_l_
     #mu2_q_mu1-2*t125+q1s*gp_mu3_mu2_mu1+2*t1236-2*gp_mu3_q_mu1*dqmu2-2
     #*gp_mu3_q_mu1*dpmu2)*t19*t8+t270*B0_13/2+(4*(t496+t280-t217+t226+t
     #230+t988)*t19*gp_mu3_l_mu1+8*(t1254+t1255)*t19*gp_mu3-t998+4*(t114
     #9+t288-t530-t1150+t1152+t1154)*t19*gp_mu3_l_mu2+t1023+4*(t496+t278
     #+t221+t277+t1264+t226+t230+t988)*t19*gp_mu3_q_mu1-4*(t1149+t1269+t
     #1175+2*t288-t1177-t526-t1178)*t19*gp_mu3_mu2_q-4*t1279*t19*gp_mu3_
     #l_q)*Dij_1234(1,3)
      s1 = s2+(8*dpmu1*t1286*t205+4*t1286*t8*gp_l_q_mu1+4*t288*t261+2*q2
     #s*t8*gp_mu3_l_mu2_q_mu1)*Dij_1345(3,2)+(-4*t1034*t19*gp_mu3_l_mu1-
     #t459+4*t1305*t25+8*t1305*t96+2*(q2s-q1s-t6)*t19*gp_mu3_mu2_mu1-t47
     #4)*C0_134+(t1316-4*dppmu2*gp_mu3_l_mu1-4*dppmu3*gp_l_mu2_mu1-8*dpp
     #mu3*t16*gp_mu1-8*dpmu1*t1329*t91+t1336+t1338+4*ppq2*gp_mu3_mu2_mu1
     #-4*t1329*t8*gp_mu3_q_mu1+t1344)*D0_1345+(4*(t68+4*dlmu2+6*dqmu2)*t
     #19*gp_mu3_l_mu1+4*(8*t660+4*t663+8*t501-t667+8*t331-t671)*t19*gp_m
     #u3-t61-4*(t75+2*dlmu1)*t19*gp_mu3_l_mu2-t82+8*dlmu1*t19*gp_mu3_mu2
     #_q-12*t434+16*t273*gp_mu3_q_mu1)*Dij_1234(3,13)
      s3 = s1+(4*(t110+t111+t496+t278+t280-t221-t217+t988+t115+2*t226)*t
     #19*gp_mu3_l_mu1+4*(t1044+4*t1040+8*t502+t1383+4*t504+4*t506+2*t519
     #-4*t524-2*t527-4*t528+t1398)*t19*gp_mu3-2*t430*gp_mu3_l_mu2_q_mu1-
     #8*t1405*t125-2*t1408*t19*gp_mu3_mu2_mu1+4*(t277+t1264+2*t278+t226)
     #*t19*gp_mu3_q_mu1-t1093-4*(t1094+t663+t1276-t1183-t1095-t668-t1096
     #-t1097)*t19*gp_mu3_l_q)*Dij_1234(2,3)
      s2 = s3+(4*t1426*gp_mu3_l_mu1+4*(t1438+t1451)*t19*gp_mu3+4*(t1456+
     #t1458+t1460+t1461+t1269+t1150+t1151-t1152-t1154)*t19*gp_mu3_l_mu2+
     #t1467+4*t1426*gp_mu3_q_mu1+4*(t1149+t288+t1175-t530+t1152+t1154)*t
     #19*gp_mu3_mu2_q-8*(t50+2*t557+t1047+t52+t331+t1184)*t19*gp_mu3_l_q
     #)*Dij_1234(1,2)+t1506*Dij_1345(2,1)+(-16*t119+t124+16*t126+t135-t1
     #38-4*(t405+t12+q1s+t406+q2s)*t19*t139)*Dij_1234(1,1)
      s3 = s2+(-8*t16*t1518*t21+2*t804*t25-4*t1518*t19*t261-2*t804*gp_mu
     #3_mu2_mu1+4*t875*t34+4*t62*gp_mu3_mu2_q)*Cij_123(2,2)+(8*t1537*dpm
     #u1*t9-t1541-t336-8*dpmu1*t1542*t91-t340-t337-t330+t1546-4*t1542*t8
     #*gp_mu3_q_mu1+4*t1537*t8*gp_mu2_q_mu1)*Dij_1345(3,13)
      t1648 = s3+(-8*t1560*dpmu1*t9-4*t86*gp_l_mu2_mu1+t1567*gp_mu1+16*d
     #pmu1*t1569*t205-8*dpmu1*t1575*t91+8*t1569*t8*gp_l_q_mu1+(t408+t407
     #+t405+t682+t406)*gp_mu3_mu2_mu1-4*t1575*t8*gp_mu3_q_mu1-4*t1560*t8
     #*gp_mu2_q_mu1)*Dij_1345(1,3)+(-8*t1593*dpmu1*t9+8*t1597*dpmu1*t16*
     #t21+2*t1602*t25+4*t1602*t96-2*t1602*gp_mu3_mu2_mu1+4*t1597*t16*t34
     #-4*t1593*t8*gp_mu2_q_mu1)*Cij_145(1,1)+(4*t1619*gp_mu3_l_mu1+8*t16
     #22*t19*gp_mu3+4*(t1460+t1269+4*t821+3*t530-t1154)*t19*gp_mu3_l_mu2
     #-4*t1236+4*t1619*gp_mu3_q_mu1+4*(-t1460+t1175+t288-t1150-t1151+t11
     #54)*t19*gp_mu3_mu2_q-8*dpmu1*(3*dpmu2+t450+dlmu2)*t19*gp_mu3_l_q)*
     #Dij_1234(2,4)
      t1651 = 4*q1s*t1034*t1036
      t1658 = 2*t531
      t1660 = 6*t830+2*t1430+t1431+t1039+t1434+t1433+4*t514-t1437-t1435+
     #2*t837-t1439-t1441-t1442-8*t1446-4*t1443-t1445-t1111+t1658-4*t838-
     #t1449
      t1675 = 2*t128
      t1695 = q1s+q2s-t771+t13+t12-t396-t397+t6
      t1696 = t1695*dppmu3
      t1699 = -t934+t166+t168+t361+t758+t369+t299-t366+t174-t933+t368-t5
     #99
      t1700 = t176-t162-t170-t932-t931+t929-t300+t601-t367+t370-t598+t92
     #6
      t1703 = t180+t597+t179-t181-t182-t165-t185+t362-t175-t163-t173+t18
     #9
      t1705 = 2*t600
      t1706 = -t192-t715-t195+t615+t197+t171-t1705+t202+t201-t620-t622+t
     #200
      t1707 = 2*t1703+t1706
      t1711 = t1695*dpmu1
      t1724 = t1695*t16
      t1741 = pq2*t19
      t1751 = dlmu2*gp_mu3_q_mu1
      t1761 = t8*dppmu3
      t1765 = ppq2+ppp+ppq1
      t1770 = t1765*t19
      t1790 = t1695*t19
      t1815 = t1039+4*t1133+t1433+t1135-t1441-t1442-t1126-t1445-t1111+t1
     #658-t1449-t1450
      t1864 = t174+t166-t170+t926-t755-t162-t932+t929+t751+t761+t759+t17
     #6-t931+t753-t760+t168-t934-t933
      t1866 = -t397+t13+t460-t396+t12+t461
      t1879 = 8*t149*t594
      t1883 = 2*t597-t714-t190+t981-t187-t612-t715+t202+t201+t618-t1705
      t1890 = t13+t12+t634-t771
      t1898 = 4*t149*t648
      t1931 = -t714-t190-t188-t187-t981-t715+t201+t171+t202
      t1935 = -t210+t111+t225+t226
      t1942 = t429*t8
      t1955 = t143-t902+t144-t904+t906-t146-t148-2*t909+t911-4*t907+t155
     #9+t1557
      t1978 = t164-t931-t932+24*t175+24*t167+24*t165+24*t173+32*t182+32*
     #t185+32*t194+32*t191
      t1981 = -t180-t179-t183-t165-t186-t192-t167-t195+t198+t196
      t1989 = -t375-t1326+t214+t111+t216+4*t108+8*t379+t217+t221+t229+t2
     #26+t230+4*t385+8*t112
      t2008 = 8*t1955*dpmu1*t9-t157*gp_mu3_l_mu1+(-8*dpmu3-8*dqmu3-t415+
     #t1536)*gp_l_mu2_mu1+(-24*t161+t926-16*t197-32*t199+t929-24*t169+t6
     #01-t172-t758-t598+t1978)*gp_mu1-16*dpmu1*t1981*t205-8*dpmu1*t1989*
     #t91-8*t1981*t8*gp_l_q_mu1-t1336-t1338+(12*q1q2-t770+16*pq1+t735+12
     #*pq2+8*q1s-t409-t772)*gp_mu3_mu2_mu1-4*t1989*t8*gp_mu3_q_mu1+4*t19
     #55*t8*gp_mu2_q_mu1
      t2022 = t1131+t510+t1113+t509-t515-t516-t518+t1139-t1124-t1118-t53
     #1-t1137-t1138
      s3 = (t1651+4*t1660*t19*gp_mu3+2*t65*gp_mu3_l_mu2_q_mu1+4*(t1461+t
     #1456+t1458-t1152)*t19*gp_mu3_l_mu2-2*(6*pq1*q1s+2*t1429+2*t1161+t1
     #077+4*t831+t1675)*t19*gp_mu3_mu2_mu1+4*(t380+t1105+t217+t113+t1494
     #)*t19*gp_mu3_q_mu1+4*(t1149+t1460+t530+t1150+t1152)*t19*gp_mu3_mu2
     #_q+4*(t46+t1275-t51-t1095-t53-t1184)*t19*gp_mu3_l_q)*Dij_1234(2,2)
     #+(8*t1696*t594+(t1699+t1700)*gp_mu1-8*dpmu1*t1707*t205-8*t1711*t30
     #7-4*t1707*t8*gp_l_q_mu1-2*t1695*t8*gp_mu3_l_mu2_q_mu1-4*t1711*t261
     #+(t405+t682-t772+t407-t770-t409+t408+t406)*gp_mu3_mu2_mu1-4*t1724*
     #t267+4*t1696*t648)*Dij_1345(3,3)
      s2 = s3+t814*Cij_123(1,1)+(8*pq2*t1034*t1036+8*(t503+t504+t507+t51
     #9-t525-t527-t529-t532-t534+t1137)*t19*gp_mu3-4*t1741*gp_mu3_l_mu2_
     #q_mu1+4*(-t13+q2s)*t19*t125-4*pq2*t429*t264+8*t1741*t1751-t1093-8*
     #(t501-t523-t331)*t19*gp_mu3_l_q)*Dij_1234(3,7)
      s1 = s2+(-8*t1761*dpmu1*gp_mu2+16*t1765*dpmu1*t16*t21+4*t1770*t25+
     #8*t1770*t96-4*t1770*gp_mu3_mu2_mu1+8*t1765*t16*t34-4*t1761*gp_mu2_
     #q_mu1)*C0_145+(-8*t923*dpmu1*t9+8*t1711*t16*t21+2*t1790*t25+4*t171
     #1*t28-2*t1790*gp_mu3_mu2_mu1+4*t1724*t34-4*t923*t8*gp_mu2_q_mu1)*C
     #ij_145(2,2)+(-8*dpmu1*(t111+t110-t115)*t118+8*t820+8*t535*t822-8*t
     #535*t824)*Dij_1234(3,5)+(t1651+4*t1815*t19*gp_mu3-2*(-t12+q1s)*t19
     #*gp_mu3_l_mu2_q_mu1-t847-2*(t1077+4*t1161+t1675+t1166)*t19*gp_mu3_
     #mu2_mu1+4*(t110+t217+t988)*t19*gp_mu3_q_mu1+8*t520*t824-4*(-t1094+
     #t46-t865+t1095+t1097)*t19*gp_mu3_l_q)*Dij_1234(3,8)+(-16*t1488*dpm
     #u1*t9+t1541+32*t92-t1546+16*t100-8*t1488*t8*gp_mu2_q_mu1)*Dij_1345
     #(3,11)
      s2 = s1+(-4*t1200*t19*gp_mu3_l_mu1-t459+4*t1405*t25+8*t1405*t96-4*
     #t1405*gp_mu3_mu2_mu1-t474)*Cij_134(2,3)+(t1316+t1864*gp_mu1-t1014-
     #8*t1866*dpmu1*t307-t1020-t998-t999+(t786+t787+t407+t406-t772-t770)
     #*gp_mu3_mu2_mu1-4*t1866*t16*t267+t1344)*Dij_1345(3,5)+(t1879-t301*
     #gp_mu1-8*dpmu1*t1883*t205-t309-4*t1883*t8*gp_l_q_mu1-2*t1890*t8*gp
     #_mu3_l_mu2_q_mu1-4*t1890*dpmu1*t261+t320-t323+t1898)*Dij_1345(3,8)
     #+(4*q2s*t1034*t1036+4*(t1041+t1383+2*t504+2*t532-t1392-t1394-t1396
     #)*t19*gp_mu3-2*t285*gp_mu3_l_mu2_q_mu1-4*t285*t125-2*q2s*t429*t264
     #+4*t285*t1751-4*(t663-t1096)*t19*gp_mu3_l_q)*Dij_1234(3,3)
      t2053 = s2+t447*Cij_134(2,4)+(t1879-4*dlmu2*gp_mu3_l_mu1+8*dpmu1*t
     #1931*t205+8*dpmu1*t1935*t91+4*t1931*t8*gp_l_q_mu1+2*t1942*gp_mu3_l
     #_mu2_q_mu1+4*t1942*t125+4*t1935*t8*gp_mu3_q_mu1+t1898)*Dij_1345(2,
     #2)+t2008*Dij_1345(2,5)+((-t370-t368-t361+t367+t366-t369)*gp_mu1+8*
     #t90*gp_mu3+t1467+t1486*gp_mu3_q_mu1)*Dij_1345(3,1)+(8*pq1*t1034*t1
     #036+8*t2022*t19*gp_mu3-4*(pq1-pq2)*t19*gp_mu3_l_mu2_q_mu1-8*(pq1-q
     #1q2)*t19*t125-4*(t128+t131+t1164+t130)*t19*gp_mu3_mu2_mu1+8*(t215+
     #t220+t114)*t19*gp_mu3_q_mu1+8*(-q1q2+pq2)*t19*t139-8*(t48-t501-t50
     #+t523+t331-t54)*t19*gp_mu3_l_q)*Dij_1234(3,10)
      virt_tri_box = t566+t1194+t1648+t2053
      end
