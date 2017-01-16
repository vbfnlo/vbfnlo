      function virt_pent(q1s,q2s,pq1,pq2,ppp,ppq1,ppq2,q1q2)
      implicit none
      complex * 16 virt_pent
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
      complex * 16 t5
      complex * 16 t7
      complex * 16 t8
      complex * 16 t12
      complex * 16 t18
      complex * 16 t20
      complex * 16 t21
      complex * 16 t23
      complex * 16 t26
      complex * 16 t28
      complex * 16 t34
      complex * 16 t36
      complex * 16 t37
      complex * 16 t38
      complex * 16 t39
      complex * 16 t44
      complex * 16 t45
      complex * 16 t46
      complex * 16 t47
      complex * 16 t54
      complex * 16 t56
      complex * 16 t57
      complex * 16 t58
      complex * 16 t59
      complex * 16 t60
      complex * 16 t61
      complex * 16 t62
      complex * 16 t63
      complex * 16 t64
      complex * 16 t69
      complex * 16 t70
      complex * 16 t71
      complex * 16 t72
      complex * 16 t73
      complex * 16 t79
      complex * 16 t83
      complex * 16 t84
      complex * 16 t85
      complex * 16 t86
      complex * 16 t87
      complex * 16 t88
      complex * 16 t89
      complex * 16 t92
      complex * 16 t93
      complex * 16 t94
      complex * 16 t96
      complex * 16 t99
      complex * 16 t100
      complex * 16 t101
      complex * 16 t102
      complex * 16 t103
      complex * 16 t106
      complex * 16 t111
      complex * 16 t115
      complex * 16 t120
      complex * 16 t122
      complex * 16 t124
      complex * 16 t131
      complex * 16 t132
      complex * 16 t133
      complex * 16 t134
      complex * 16 t136
      complex * 16 t137
      complex * 16 t138
      complex * 16 t142
      complex * 16 t143
      complex * 16 t146
      complex * 16 t147
      complex * 16 t148
      complex * 16 t149
      complex * 16 t150
      complex * 16 t156
      complex * 16 t157
      complex * 16 t158
      complex * 16 t159
      complex * 16 t160
      complex * 16 t161
      complex * 16 t162
      complex * 16 t165
      complex * 16 t174
      complex * 16 t177
      complex * 16 t187
      complex * 16 t191
      complex * 16 t192
      complex * 16 t193
      complex * 16 t194
      complex * 16 t195
      complex * 16 t203
      complex * 16 t204
      complex * 16 t206
      complex * 16 t209
      complex * 16 t215
      complex * 16 t217
      complex * 16 t218
      complex * 16 t219
      complex * 16 t220
      complex * 16 t222
      complex * 16 t233
      complex * 16 t243
      complex * 16 t244
      complex * 16 t245
      complex * 16 t246
      complex * 16 t247
      complex * 16 t248
      complex * 16 t249
      complex * 16 t250
      complex * 16 t251
      complex * 16 t252
      complex * 16 t253
      complex * 16 t254
      complex * 16 t255
      complex * 16 t256
      complex * 16 t257
      complex * 16 t260
      complex * 16 t261
      complex * 16 t262
      complex * 16 t263
      complex * 16 t264
      complex * 16 t265
      complex * 16 t268
      complex * 16 t269
      complex * 16 t271
      complex * 16 t272
      complex * 16 t273
      complex * 16 t274
      complex * 16 t275
      complex * 16 t278
      complex * 16 t279
      complex * 16 t282
      complex * 16 t283
      complex * 16 t284
      complex * 16 t285
      complex * 16 t287
      complex * 16 t288
      complex * 16 t289
      complex * 16 t290
      complex * 16 t292
      complex * 16 t301
      complex * 16 t303
      complex * 16 t306
      complex * 16 t307
      complex * 16 t308
      complex * 16 t309
      complex * 16 t310
      complex * 16 t311
      complex * 16 t312
      complex * 16 t315
      complex * 16 t316
      complex * 16 t317
      complex * 16 t318
      complex * 16 t321
      complex * 16 t322
      complex * 16 t325
      complex * 16 t332
      complex * 16 t333
      complex * 16 t338
      complex * 16 t339
      complex * 16 t342
      complex * 16 t348
      complex * 16 t349
      complex * 16 t357
      complex * 16 t358
      complex * 16 t359
      complex * 16 t360
      complex * 16 t361
      complex * 16 t363
      complex * 16 t369
      complex * 16 t378
      complex * 16 t380
      complex * 16 t386
      complex * 16 t387
      complex * 16 t395
      complex * 16 t398
      complex * 16 t399
      complex * 16 t403
      complex * 16 t408
      complex * 16 t410
      complex * 16 t411
      complex * 16 t412
      complex * 16 t413
      complex * 16 t414
      complex * 16 t415
      complex * 16 t416
      complex * 16 t417
      complex * 16 t418
      complex * 16 t419
      complex * 16 t424
      complex * 16 t425
      complex * 16 t426
      complex * 16 t427
      complex * 16 t428
      complex * 16 t429
      complex * 16 t432
      complex * 16 t433
      complex * 16 t434
      complex * 16 t436
      complex * 16 t437
      complex * 16 t438
      complex * 16 t441
      complex * 16 t442
      complex * 16 t443
      complex * 16 t444
      complex * 16 t445
      complex * 16 t446
      complex * 16 t447
      complex * 16 t448
      complex * 16 t451
      complex * 16 t452
      complex * 16 t457
      complex * 16 t460
      complex * 16 t461
      complex * 16 t462
      complex * 16 t463
      complex * 16 t469
      complex * 16 t470
      complex * 16 t471
      complex * 16 t473
      complex * 16 t479
      complex * 16 t480
      complex * 16 t483
      complex * 16 t484
      complex * 16 t495
      complex * 16 t496
      complex * 16 t498
      complex * 16 t499
      complex * 16 t504
      complex * 16 t508
      complex * 16 t510
      complex * 16 t512
      complex * 16 t513
      complex * 16 t518
      complex * 16 t520
      complex * 16 t521
      complex * 16 t522
      complex * 16 t525
      complex * 16 t526
      complex * 16 t527
      complex * 16 t528
      complex * 16 t531
      complex * 16 t532
      complex * 16 t533
      complex * 16 t534
      complex * 16 t535
      complex * 16 t538
      complex * 16 t539
      complex * 16 t540
      complex * 16 t543
      complex * 16 t544
      complex * 16 t545
      complex * 16 t546
      complex * 16 t547
      complex * 16 t555
      complex * 16 t556
      complex * 16 t557
      complex * 16 t560
      complex * 16 t563
      complex * 16 t564
      complex * 16 t565
      complex * 16 t567
      complex * 16 t568
      complex * 16 t569
      complex * 16 t571
      complex * 16 t574
      complex * 16 t577
      complex * 16 t593
      complex * 16 t595
      complex * 16 t596
      complex * 16 t597
      complex * 16 t599
      complex * 16 t605
      complex * 16 t606
      complex * 16 t608
      complex * 16 t610
      complex * 16 t611
      complex * 16 t613
      complex * 16 t614
      complex * 16 t615
      complex * 16 t616
      complex * 16 t617
      complex * 16 t618
      complex * 16 t619
      complex * 16 t621
      complex * 16 t622
      complex * 16 t623
      complex * 16 t624
      complex * 16 t625
      complex * 16 t626
      complex * 16 t628
      complex * 16 t632
      complex * 16 t634
      complex * 16 t635
      complex * 16 t636
      complex * 16 t640
      complex * 16 t641
      complex * 16 t642
      complex * 16 t645
      complex * 16 t648
      complex * 16 t655
      complex * 16 t658
      complex * 16 t659
      complex * 16 t660
      complex * 16 t661
      complex * 16 t662
      complex * 16 t666
      complex * 16 t667
      complex * 16 t671
      complex * 16 t672
      complex * 16 t682
      complex * 16 t687
      complex * 16 t689
      complex * 16 t696
      complex * 16 t701
      complex * 16 t707
      complex * 16 t708
      complex * 16 t709
      complex * 16 t710
      complex * 16 t713
      complex * 16 t714
      complex * 16 t715
      complex * 16 t718
      complex * 16 t719
      complex * 16 t720
      complex * 16 t723
      complex * 16 t726
      complex * 16 t727
      complex * 16 t729
      complex * 16 t733
      complex * 16 t734
      complex * 16 t742
      complex * 16 t743
      complex * 16 t745
      complex * 16 t747
      complex * 16 t761
      complex * 16 t765
      complex * 16 t766
      complex * 16 t767
      complex * 16 t768
      complex * 16 t775
      complex * 16 t777
      complex * 16 t778
      complex * 16 t784
      complex * 16 t785
      complex * 16 t789
      complex * 16 t793
      complex * 16 t794
      complex * 16 t795
      complex * 16 t798
      complex * 16 t799
      complex * 16 t800
      complex * 16 t801
      complex * 16 t803
      complex * 16 t808
      complex * 16 t809
      complex * 16 t815
      complex * 16 t817
      complex * 16 t819
      complex * 16 t820
      complex * 16 t823
      complex * 16 t825
      complex * 16 t826
      complex * 16 t827
      complex * 16 t831
      complex * 16 t838
      complex * 16 t840
      complex * 16 t842
      complex * 16 t851
      complex * 16 t858
      complex * 16 t871
      complex * 16 t874
      complex * 16 t876
      complex * 16 t878
      complex * 16 t885
      complex * 16 t886
      complex * 16 t888
      complex * 16 t889
      complex * 16 t894
      complex * 16 t896
      complex * 16 t897
      complex * 16 t901
      complex * 16 t903
      complex * 16 t913
      complex * 16 t919
      complex * 16 t921
      complex * 16 t923
      complex * 16 t929
      complex * 16 t931
      complex * 16 t932
      complex * 16 t936
      complex * 16 t938
      complex * 16 t939
      complex * 16 t940
      complex * 16 t941
      complex * 16 t943
      complex * 16 t948
      complex * 16 t949
      complex * 16 t950
      complex * 16 t951
      complex * 16 t953
      complex * 16 t957
      complex * 16 t962
      complex * 16 t968
      complex * 16 t972
      complex * 16 t973
      complex * 16 t977
      complex * 16 t981
      complex * 16 t985
      complex * 16 t988
      complex * 16 t989
      complex * 16 t993
      complex * 16 t994
      complex * 16 t995
      complex * 16 t996
      complex * 16 t998
      complex * 16 t1000
      complex * 16 t1001
      complex * 16 t1007
      complex * 16 t1008
      complex * 16 t1022
      complex * 16 t1023
      complex * 16 t1024
      complex * 16 t1026
      complex * 16 t1031
      complex * 16 t1032
      complex * 16 t1033
      complex * 16 t1034
      complex * 16 t1035
      complex * 16 t1040
      complex * 16 t1041
      complex * 16 t1042
      complex * 16 t1043
      complex * 16 t1044
      complex * 16 t1046
      complex * 16 t1048
      complex * 16 t1052
      complex * 16 t1053
      complex * 16 t1055
      complex * 16 t1056
      complex * 16 t1057
      complex * 16 t1059
      complex * 16 t1064
      complex * 16 t1065
      complex * 16 t1067
      complex * 16 t1083
      complex * 16 t1087
      complex * 16 t1089
      complex * 16 t1090
      complex * 16 t1091
      complex * 16 t1092
      complex * 16 t1093
      complex * 16 t1101
      complex * 16 t1105
      complex * 16 t1106
      complex * 16 t1108
      complex * 16 t1115
      complex * 16 t1116
      complex * 16 t1120
      complex * 16 t1121
      complex * 16 t1122
      complex * 16 t1124
      complex * 16 t1126
      complex * 16 t1131
      complex * 16 t1132
      complex * 16 t1136
      complex * 16 t1139
      complex * 16 t1140
      complex * 16 t1142
      complex * 16 t1143
      complex * 16 t1146
      complex * 16 t1148
      complex * 16 t1151
      complex * 16 t1152
      complex * 16 t1154
      complex * 16 t1165
      complex * 16 t1166
      complex * 16 t1167
      complex * 16 t1168
      complex * 16 t1169
      complex * 16 t1174
      complex * 16 t1179
      complex * 16 t1180
      complex * 16 t1185
      complex * 16 t1194
      complex * 16 t1196
      complex * 16 t1201
      complex * 16 t1202
      complex * 16 t1207
      complex * 16 t1208
      complex * 16 t1216
      complex * 16 t1217
      complex * 16 t1218
      complex * 16 t1219
      complex * 16 t1226
      complex * 16 t1227
      complex * 16 t1229
      complex * 16 t1230
      complex * 16 t1232
      complex * 16 t1237
      complex * 16 t1238
      complex * 16 t1239
      complex * 16 t1242
      complex * 16 t1243
      complex * 16 t1245
      complex * 16 t1248
      complex * 16 t1252
      complex * 16 t1253
      complex * 16 t1255
      complex * 16 t1258
      complex * 16 t1259
      complex * 16 t1263
      complex * 16 t1264
      complex * 16 t1266
      complex * 16 t1267
      complex * 16 t1269
      complex * 16 t1270
      complex * 16 t1271
      complex * 16 t1273
      complex * 16 t1276
      complex * 16 t1281
      complex * 16 t1282
      complex * 16 t1283
      complex * 16 t1285
      complex * 16 t1290
      complex * 16 t1293
      complex * 16 t1294
      complex * 16 t1298
      complex * 16 t1299
      complex * 16 t1300
      complex * 16 t1305
      complex * 16 t1308
      complex * 16 t1309
      complex * 16 t1313
      complex * 16 t1314
      complex * 16 t1315
      complex * 16 t1320
      complex * 16 t1321
      complex * 16 t1324
      complex * 16 t1327
      complex * 16 t1333
      complex * 16 t1334
      complex * 16 t1335
      complex * 16 t1338
      complex * 16 t1339
      complex * 16 t1346
      complex * 16 t1347
      complex * 16 t1350
      complex * 16 t1352
      complex * 16 t1361
      complex * 16 t1363
      complex * 16 t1368
      complex * 16 t1369
      complex * 16 t1373
      complex * 16 t1374
      complex * 16 t1377
      complex * 16 t1379
      complex * 16 t1381
      complex * 16 t1382
      complex * 16 t1383
      complex * 16 t1384
      complex * 16 t1387
      complex * 16 t1393
      complex * 16 t1394
      complex * 16 t1395
      complex * 16 t1401
      complex * 16 t1402
      complex * 16 t1403
      complex * 16 t1404
      complex * 16 t1413
      complex * 16 t1414
      complex * 16 t1417
      complex * 16 t1419
      complex * 16 t1422
      complex * 16 t1423
      complex * 16 t1426
      complex * 16 t1427
      complex * 16 t1428
      complex * 16 t1431
      complex * 16 t1432
      complex * 16 t1436
      complex * 16 t1440
      complex * 16 t1441
      complex * 16 t1442
      complex * 16 t1453
      complex * 16 t1474
      complex * 16 t1478
      complex * 16 t1513
      complex * 16 t1514
      complex * 16 t1516
      complex * 16 t1518
      complex * 16 t1519
      complex * 16 t1521
      complex * 16 t1527
      complex * 16 t1537
      complex * 16 t1538
      complex * 16 t1539
      complex * 16 t1540
      complex * 16 t1542
      complex * 16 t1543
      complex * 16 t1546
      complex * 16 t1552
      complex * 16 t1561
      complex * 16 t1564
      complex * 16 t1571
      complex * 16 t1576
      complex * 16 t1577
      complex * 16 t1578
      complex * 16 t1583
      complex * 16 t1584
      complex * 16 t1586
      complex * 16 t1593
      complex * 16 t1595
      complex * 16 t1596
      complex * 16 t1597
      complex * 16 t1598
      complex * 16 t1603
      complex * 16 t1604
      complex * 16 t1608
      complex * 16 t1609
      complex * 16 t1615
      complex * 16 t1616
      complex * 16 t1617
      complex * 16 t1621
      complex * 16 t1622
      complex * 16 t1624
      complex * 16 t1639
      complex * 16 t1640
      complex * 16 t1644
      complex * 16 t1646
      complex * 16 t1648
      complex * 16 t1651
      complex * 16 t1656
      complex * 16 t1658
      complex * 16 t1660
      complex * 16 t1661
      complex * 16 t1662
      complex * 16 t1678
      complex * 16 t1680
      complex * 16 t1681
      complex * 16 t1682
      complex * 16 t1688
      complex * 16 t1692
      complex * 16 t1693
      complex * 16 t1695
      complex * 16 t1705
      complex * 16 t1711
      complex * 16 t1714
      complex * 16 t1717
      complex * 16 t1723
      complex * 16 t1727
      complex * 16 t1731
      complex * 16 t1734
      complex * 16 t1735
      complex * 16 t1743
      complex * 16 t1744
      complex * 16 t1745
      complex * 16 t1768
      complex * 16 t1773
      complex * 16 t1778
      complex * 16 t1780
      complex * 16 t1783
      complex * 16 t1787
      complex * 16 t1793
      complex * 16 t1795
      complex * 16 t1798
      complex * 16 t1800
      complex * 16 t1808
      complex * 16 t1827
      complex * 16 t1829
      complex * 16 t1839
      complex * 16 t1840
      complex * 16 t1850
      complex * 16 t1854
      complex * 16 t1857
      complex * 16 t1867
      complex * 16 t1869
      complex * 16 t1872
      complex * 16 t1876
      complex * 16 t1884
      complex * 16 t1912
      complex * 16 t1914
      complex * 16 t1915
      complex * 16 t1932
      complex * 16 t1933
      complex * 16 t1935
      complex * 16 t1939
      complex * 16 t1943
      complex * 16 t1949
      complex * 16 t1971
      complex * 16 t1979
      complex * 16 t1980
      complex * 16 t1983
      complex * 16 t1986
      complex * 16 t1988
      complex * 16 t1989
      complex * 16 t1992
      complex * 16 t1995
      complex * 16 t1996
      complex * 16 t2002
      complex * 16 t2004
      complex * 16 t2006
      complex * 16 t2012
      complex * 16 t2013
      complex * 16 t2014
      complex * 16 t2021
      complex * 16 t2024
      complex * 16 t2034
      complex * 16 t2037
      complex * 16 t2039
      complex * 16 t2050
      complex * 16 t2053
      complex * 16 t2058
      complex * 16 t2063
      complex * 16 t2069
      complex * 16 t2086
      complex * 16 t2093
      complex * 16 t2098
      complex * 16 t2100
      complex * 16 t2101
      complex * 16 t2104
      complex * 16 t2110
      complex * 16 t2112
      complex * 16 t2113
      complex * 16 t2114
      complex * 16 t2116
      complex * 16 t2119
      complex * 16 t2122
      complex * 16 t2126
      complex * 16 t2130
      complex * 16 t2131
      complex * 16 t2136
      complex * 16 t2141
      complex * 16 t2143
      complex * 16 t2147
      complex * 16 t2151
      complex * 16 t2164
      complex * 16 t2167
      complex * 16 t2169
      complex * 16 t2173
      complex * 16 t2175
      complex * 16 t2179
      complex * 16 t2187
      complex * 16 t2189
      complex * 16 t2195
      complex * 16 t2198
      complex * 16 t2202
      complex * 16 t2204
      complex * 16 t2225
      complex * 16 t2228
      complex * 16 t2229
      complex * 16 t2232
      complex * 16 t2233
      complex * 16 t2234
      complex * 16 t2238
      complex * 16 t2255
      complex * 16 t2271
      complex * 16 t2282
      complex * 16 t2285
      complex * 16 t2286
      complex * 16 t2288
      complex * 16 t2291
      complex * 16 t2305
      complex * 16 t2307
      complex * 16 t2314
      complex * 16 t2317
      complex * 16 t2325
      complex * 16 t2327
      complex * 16 t2329
      complex * 16 t2331
      complex * 16 t2333
      complex * 16 t2336
      complex * 16 t2345
      complex * 16 t2348
      complex * 16 t2355
      complex * 16 t2357
      complex * 16 t2394
      complex * 16 t2402
      complex * 16 t2414
      complex * 16 t2417
      complex * 16 t2418
      complex * 16 t2419
      complex * 16 t2426
      complex * 16 t2432
      complex * 16 t2466
      complex * 16 t2469
      complex * 16 t2471
      complex * 16 t2475
      complex * 16 t2486
      complex * 16 t2494
      complex * 16 t2501
      complex * 16 t2525
      complex * 16 t2533
      complex * 16 t2548
      complex * 16 t2564
      complex * 16 t2567
      complex * 16 t2570
      complex * 16 t2573
      complex * 16 t2576
      complex * 16 t2582
      complex * 16 t2590
      complex * 16 t2595
      complex * 16 t2611
      complex * 16 t2618
      complex * 16 t2620
      complex * 16 t2624
      complex * 16 t2627
      complex * 16 t2629
      complex * 16 t2632
      complex * 16 t2639
      complex * 16 t2653
      complex * 16 t2655
      complex * 16 t2657
      complex * 16 t2671
      complex * 16 t2674
      complex * 16 t2682
      complex * 16 t2683
      complex * 16 t2685
      complex * 16 t2699
      complex * 16 t2700
      complex * 16 t2702
      complex * 16 t2703
      complex * 16 t2709
      complex * 16 t2716
      complex * 16 t2736
      complex * 16 t2760
      complex * 16 t2765
      complex * 16 t2766
      complex * 16 t2776
      complex * 16 t2780
      complex * 16 t2782
      complex * 16 t2785
      complex * 16 t2794
      complex * 16 t2796
      complex * 16 t2805
      complex * 16 t2806
      complex * 16 t2807
      complex * 16 t2808
      complex * 16 t2815
      complex * 16 t2816
      complex * 16 t2822
      complex * 16 t2824
      complex * 16 t2838
      complex * 16 t2839
      complex * 16 t2843
      complex * 16 t2848
      complex * 16 t2849
      complex * 16 t2854
      complex * 16 t2861
      complex * 16 t2866
      complex * 16 t2867
      complex * 16 t2882
      complex * 16 t2883
      complex * 16 t2889
      complex * 16 t2893
      complex * 16 t2897
      complex * 16 t2907
      complex * 16 t2914
      complex * 16 t2937
      complex * 16 t2950
      complex * 16 t2953
      complex * 16 t2956
      complex * 16 t2959
      complex * 16 t2960
      complex * 16 t2963
      complex * 16 t2966
      complex * 16 t2970
      complex * 16 t2973
      complex * 16 t2982
      complex * 16 t2985
      complex * 16 t3009
      complex * 16 t3012
      complex * 16 t3020
      complex * 16 t3023
      complex * 16 t3031
      complex * 16 t3032
      complex * 16 t3034
      complex * 16 t3057
      complex * 16 t3058
      complex * 16 t3060
      complex * 16 t3062
      complex * 16 t3063
      complex * 16 t3067
      complex * 16 t3068
      complex * 16 t3070
      complex * 16 t3072
      complex * 16 t3082
      complex * 16 t3083
      complex * 16 t3084
      complex * 16 t3090
      complex * 16 t3094
      complex * 16 t3100
      complex * 16 t3103
      complex * 16 t3115
      complex * 16 t3117
      complex * 16 t3120
      complex * 16 t3126
      complex * 16 t3129
      complex * 16 t3135
      complex * 16 t3160
      complex * 16 t3163
      complex * 16 t3165
      complex * 16 t3166
      complex * 16 t3167
      complex * 16 t3209
      complex * 16 t3217
      complex * 16 t3223
      complex * 16 t3227
      complex * 16 t3229
      complex * 16 t3238
      complex * 16 t3239
      complex * 16 t3243
      complex * 16 t3248
      complex * 16 t3249
      complex * 16 t3250
      complex * 16 t3253
      complex * 16 t3270
      complex * 16 t3271
      complex * 16 t3272
      complex * 16 t3276
      complex * 16 t3289
      complex * 16 t3295
      complex * 16 t3302
      complex * 16 t3313
      complex * 16 t3321
      complex * 16 t3333
      complex * 16 t3337
      complex * 16 t3343
      complex * 16 t3349
      complex * 16 t3355
      complex * 16 t3379
      complex * 16 t3380
      complex * 16 t3382
      complex * 16 t3385
      complex * 16 t3388
      complex * 16 t3391
      complex * 16 t3395
      complex * 16 t3408
      complex * 16 t3413
      complex * 16 t3417
      complex * 16 t3428
      complex * 16 t3445
      complex * 16 t3467
      complex * 16 t3479
      complex * 16 t3491
      complex * 16 t3494
      complex * 16 t3502
      complex * 16 t3516
      complex * 16 t3517
      complex * 16 t3524
      complex * 16 t3527
      complex * 16 t3556
      complex * 16 t3562
      complex * 16 t3568
      complex * 16 t3575
      complex * 16 t3576
      complex * 16 t3580
      complex * 16 t3583
      complex * 16 t3588
      complex * 16 t3594
      complex * 16 t3596
      complex * 16 t3597
      complex * 16 t3603
      complex * 16 t3607
      complex * 16 t3610
      complex * 16 t3614
      complex * 16 t3625
      complex * 16 t3634
      complex * 16 t3638
      complex * 16 t3641
      complex * 16 t3644
      complex * 16 t3648
      complex * 16 t3660
      complex * 16 t3661
      complex * 16 t3662
      complex * 16 t3668
      complex * 16 t3669
      complex * 16 t3678
      complex * 16 t3709
      complex * 16 t3712
      complex * 16 t3713
      complex * 16 t3714
      complex * 16 t3715
      complex * 16 t3719
      complex * 16 t3723
      complex * 16 t3724
      complex * 16 t3732
      complex * 16 t3733
      complex * 16 t3734
      complex * 16 t3738
      complex * 16 t3748
      complex * 16 t3749
      complex * 16 t3754
      complex * 16 t3755
      complex * 16 t3758
      complex * 16 t3763
      complex * 16 t3764
      complex * 16 t3802
      complex * 16 t3809
      complex * 16 t3813
      complex * 16 t3820
      complex * 16 t3823
      complex * 16 t3851
      complex * 16 t3872
      complex * 16 t3875
      complex * 16 t3888
      complex * 16 t3889
      complex * 16 t3890
      complex * 16 t3891
      complex * 16 t3892
      complex * 16 t3894
      complex * 16 t3914
      complex * 16 t3915
      complex * 16 t3918
      complex * 16 t3946
      complex * 16 t3951
      complex * 16 t3957
      complex * 16 t3962
      complex * 16 t3972
      complex * 16 t3980
      complex * 16 t3998
      complex * 16 t4004
      complex * 16 t4012
      complex * 16 t4018
      complex * 16 t4043
      complex * 16 t4068
      complex * 16 t4070
      complex * 16 t4080
      complex * 16 t4107
      complex * 16 t4123
      complex * 16 t4140
      complex * 16 t4146
      complex * 16 t4154
      complex * 16 t4158
      complex * 16 t4163
      complex * 16 t4169
      complex * 16 t4175
      complex * 16 t4181
      complex * 16 t4190
      complex * 16 t4193
      complex * 16 t4198
      complex * 16 t4212
      complex * 16 t4220
      complex * 16 t4245
      complex * 16 t4247
      complex * 16 t4250
      complex * 16 t4252
      complex * 16 t4298
      complex * 16 t4308
      complex * 16 t4315
      complex * 16 t4320
      complex * 16 t4325
      complex * 16 t4339
      complex * 16 t4354
      complex * 16 t4356
      complex * 16 t4369
      complex * 16 t4373
      complex * 16 t4379
      complex * 16 t4405
      complex * 16 t4428
      complex * 16 t4432
      complex * 16 t4447
      complex * 16 t4471
      complex * 16 t4472
      complex * 16 t4555
      complex * 16 t4570
      complex * 16 t4580
      complex * 16 t4604
      complex * 16 t4621
      complex * 16 t4631
      complex * 16 t4639
      complex * 16 t4680
      complex * 16 t4712
      complex * 16 t4747
      complex * 16 t4792
      complex * 16 t4828
      complex * 16 t4865
      complex * 16 t4867
      complex * 16 t4884
      complex * 16 t4902
      complex * 16 t4909
      complex * 16 t4910
      complex * 16 t4961
      complex * 16 t4980
      complex * 16 t4986
      complex * 16 t5042
      complex * 16 t5048
      complex * 16 s1
      complex * 16 s2
      complex * 16 s3
      complex * 16 s4
      complex * 16 s5
      complex * 16 s6
      complex * 16 s7
      complex * 16 s8
      complex * 16 s9
      complex * 16 s10
      t1 = dlmu2*gp_mu3_q_mu1
      t2 = 4*t1
      t3 = dlmu2*gp_mu3_l_mu1
      t5 = dlmu3*gp_l_mu2_mu1
      t7 = dlmu3*dqmu2
      t8 = dlmu3*dpmu2
      t12 = dpmu1*dlmu2
      t18 = gp_mu3_q_mu1*dqmu2
      t20 = dqmu1*dmu2mu3
      t21 = t20*gp_l
      t23 = dlmu2*dqmu1
      t26 = dqmu1*gp_mu3_l_mu2
      t28 = dqmu1*dlmu3
      t34 = dqmu2*dmu1mu3
      t36 = 8*t34*gp_q
      t37 = dqmu3*dqmu2
      t38 = t37*gp_mu1
      t39 = 8*t38
      t44 = 4*t26
      t45 = dqmu1*dqmu2
      t46 = dqmu1*dpmu2
      t47 = dpmu1*dqmu2
      t54 = dqmu1*dqmu3
      t56 = 16*t54*dppmu2
      t57 = t46*dqmu3
      t58 = 16*t57
      t59 = dlmu2*dqmu3
      t60 = 16*t59
      t61 = 16*dpmu3
      t62 = 64*dqmu3
      t63 = 16*dppmu3
      t64 = 16*dlmu3
      t69 = t37*dpmu1
      t70 = 16*t69
      t71 = dlmu1-dppmu1
      t72 = 16*t71
      t73 = t72*dqmu3
      t79 = 16*t45*dqmu3*gp_l
      t83 = -pq1+q1q2
      t84 = 16*t83
      t85 = t84*dpmu3
      t86 = 16*pq2
      t87 = 32*pq1
      t88 = 8*q1s
      t89 = 8*q2s
      t92 = t84*dlmu3
      t93 = 16*ppq1
      t94 = 16*pq1
      t96 = (-t93+t94+t88)*dppmu3
      t99 = t87+t88
      t100 = t99*dpmu3
      t101 = 16*q1q2
      t102 = 16*q1s
      t103 = 48*pq1
      t106 = t99*dlmu3
      t111 = -t101-t88
      t115 = t111*dlmu3
      t120 = t87+t88-t89-t101-t86
      t122 = t103+t102
      t124 = t120*dlmu3
      t131 = dqmu3+dlmu3+dpmu3
      t132 = 8*t131
      t133 = t132*dqmu1
      t134 = dppmu3*dpmu1
      t136 = 8*dlmu1
      t137 = 16*dppmu1
      t138 = t136-t137
      t142 = dppmu1*dppmu3
      t143 = 8*t142
      t146 = q2s*dpmu3
      t147 = 8*t146
      t148 = dqmu3*q2s
      t149 = 8*t148
      t150 = t88+t101+t86+t94
      t156 = 32*pq2
      t157 = t156+t101+t89
      t158 = t157*dpmu3
      t159 = t157*dqmu3
      t160 = 32*q1q2
      t161 = 48*pq2
      t162 = 16*q2s
      t165 = 16*ppq2
      t174 = t150*dlmu1+(-16*pq2-16*q1q2)*dppmu1
      t177 = t89+t156+t102+t160+t87
      t187 = (-t165+t101+t89+t86)*dppmu3
      t191 = 8*ppq2
      t192 = 12*q1s
      t193 = 24*pq1
      t194 = 4*q2s
      t195 = -t86-t191-t192-t193-t101-t194
      t203 = -ppq2-ppp-ppq1
      t204 = 16*t203
      t206 = 16*ppp
      t209 = -t165+t88+t94+t89+t101-t206-t93+t86
      t215 = -t94+t89-t88+t206+t93
      t217 = -t89+t206-t156-t101-t88-t87
      t218 = t217*dpmu1
      t219 = t206-t86-t101-t102+t93-t87
      t220 = t219*dlmu1
      t222 = -t219*dppmu1
      t233 = -t204*dlmu1
      t243 = dqmu1*dppmu3
      t244 = 16*t243
      t245 = dlmu3+dpmu3+dqmu3+dppmu3
      t246 = 16*t245
      t247 = t246*dpmu1
      t248 = dlmu1*dpmu3
      t249 = 16*t248
      t250 = dlmu1*dqmu3
      t251 = 16*t250
      t252 = dlmu1*dlmu3
      t253 = 16*t252
      t254 = dlmu1*dppmu3
      t255 = 16*t254
      t256 = t88+t101+t94+t89+t86
      t257 = t256*dmu1mu3
      t260 = t101+t89+t86
      t261 = t260*dpmu3
      t262 = t260*dqmu3
      t263 = t260*dlmu3
      t264 = -t260
      t265 = t264*dppmu3
      t268 = -t246
      t269 = t268*dqmu1
      t271 = -48*t131
      t272 = t271*dpmu1
      t273 = -t256
      t274 = t273*dmu1mu3
      t275 = 16*t142
      t278 = t268*dlmu2
      t279 = t268*dqmu2
      t282 = 32*dpmu3
      t283 = 32*dqmu3
      t284 = 32*dlmu3
      t285 = -t282-t283-t63-t284
      t287 = dppmu1*dpmu3
      t288 = 16*t287
      t289 = dppmu1*dqmu3
      t290 = 16*t289
      t292 = 16*dppmu1*dlmu3
      t301 = 16*t131
      t303 = 16*t134
      t306 = -t94-t88
      t307 = t306*dpmu3
      t308 = t306*dqmu3
      t309 = t306*dlmu3
      t310 = -t306
      t311 = t310*dppmu3
      t312 = t307+t308+t309+t311
      t315 = -t72
      t316 = t315*dpmu3
      t317 = t315*dqmu3
      t318 = t315*dlmu3
      t321 = t285*dqmu2
      t322 = t273*dmu2mu3
      t325 = t268*dpmu1
      t332 = t273*dlmu1
      t333 = t256*dppmu1
      t338 = 16*dqmu3
      t339 = 8*dppmu3
      t342 = dpmu2*dppmu3
      t348 = 8*q1q2
      t349 = t348+t193+t192
      t357 = 24*pq2
      t358 = 4*q1s
      t359 = 8*ppq1
      t360 = 12*q2s
      t361 = t357+t101+t358+t359+t94+t360
      t363 = -t361
      t369 = t101+t357+t360
      t378 = -8*q1s-8*pq1-8*q1q2+8*ppq1
      t380 = -t378
      t386 = -q1q2-pq2-q2s+ppq2
      t387 = 8*t386
      s2 = ((t85+(-t86-t87-t88-t89)*dqmu3+t92+t96)*dppmu2+(t100+(t101+t8
     #9+t102+t103+t86)*dqmu3+t106+(-t102+t93-t101-t87)*dppmu3)*dpmu2+(t1
     #11*dpmu3+(t89+t94+t86)*dqmu3+t115+16*ppq1*dppmu3)*dlmu2+(t120*dpmu
     #3+t122*dqmu3+t124+(-t102+t86+t93+t89-t87)*dppmu3)*dqmu2)*gp_mu1+(t
     #133+8*t134+t138*dpmu3+t138*dqmu3+t138*dlmu3+t143)*gp_l_mu2_q+((-t1
     #47-t149+t150*dlmu3+16*ppq2*dppmu3)*dqmu1+(t158+t159+(t160+t161+t16
     #2+t94+t88)*dlmu3+(t165-t162-t101-t156)*dppmu3)*dpmu1+t174*dpmu3+t1
     #74*dqmu3+(t177*dlmu1+(-t156-t89-t160-t94-t88)*dppmu1)*dlmu3+(-t89-
     #t88+t165-t101-t86-t94)*dppmu3*dlmu1+t187*dppmu1)*gp_mu2
      s1 = s2+(t195*dqmu1+t195*dpmu1+t195*dlmu1-t195*dppmu1)*gp_mu3_l_mu
     #2+((t204*dqmu1+(t101+t162+t156+t94-t165-t206)*dpmu1+t209*dlmu1-t20
     #9*dppmu1)*dppmu2+(t215*dqmu1+t218+t220+t222)*dpmu2+((t93+t101+t88+
     #t89+t86+t206+t165+t94)*dlmu2+t215*dqmu2)*dqmu1+((-t89+t206+t88+t16
     #5-t86)*dpmu1+t233+t204*dppmu1)*dlmu2+t217*dqmu2*dpmu1+(t220+t222)*
     #dqmu2)*gp_mu3+((t244+t247-t249-t251-t253+t255+t257)*dppmu2+(t261+t
     #262+t263+t265)*dmu1mu2+(t269+t272-t255+t274+t275)*dpmu2+(t278+t279
     #)*dqmu1+(t285*dpmu1+t288+t290+t292-t255+t274)*dlmu2+t271*dqmu2*dpm
     #u1+(-t255+t274+t275)*dqmu2)*gp_q
      t395 = s1+((t301*dqmu1+t303)*dppmu2+t312*dmu1mu2+(t285*dqmu1+t272+
     #t316+t317+t318+t275-t255)*dpmu2+(t278+t321+t322)*dqmu1+(t325+t275-
     #t255)*dlmu2+(t321+t322)*dpmu1+(t316+t317+t318+t275-t255)*dqmu2+(t3
     #32+t333)*dmu2mu3)*gp_l+((-t338-t64-t61+t339)*dppmu2+8*t342+t132*dl
     #mu2+t132*dqmu2)*gp_l_q_mu1+(t349*dpmu3+t349*dqmu3+t349*dlmu3-t349*
     #dppmu3)*gp_l_mu2_mu1+(t361*dppmu2+t363*dpmu2+t363*dlmu2+t363*dqmu2
     #)*gp_mu3_q_mu1+(t369*dpmu3+t369*dqmu3+t369*dlmu3-t369*dppmu3)*gp_m
     #u2_q_mu1+(t378*dppmu2+t380*dpmu2+t380*dlmu2+t380*dqmu2)*gp_mu3_l_m
     #u1+(t387*dqmu1+t387*dpmu1+t387*dlmu1-t387*dppmu1)*gp_mu3_mu2_q
      t398 = -t59-t7
      t399 = 16*t398
      t403 = t399*dqmu1-16*t37*dlmu1
      t408 = dmu2mu3*dpmu1
      t410 = 8*t408*gp_q
      t411 = dpmu1*gp_mu3_mu2_q
      t412 = 4*t411
      t413 = t47+t46
      t414 = 8*t413
      t415 = t414*gp_mu3
      t416 = dqmu1*dpmu3
      t417 = dpmu1*dqmu3
      t418 = -t416-t417
      t419 = 8*t418
      t424 = 16*t28
      t425 = dpmu1*dlmu3
      t426 = 16*t425
      t427 = 32*dlmu1
      t428 = t137-t427
      t429 = t428*dlmu3
      t432 = -t428
      t433 = t432*dlmu3
      t434 = t424+t426+t249+t251+t433-t255
      t436 = -t63+t61+t338+t284
      t437 = t436*dlmu2
      t438 = 16*t7
      t441 = t436*dpmu1
      t442 = t432*dpmu3
      t443 = t432*dqmu3
      t444 = 48*dlmu1
      t445 = 32*dppmu1
      t446 = t444-t445
      t447 = t446*dlmu3
      t448 = 32*t254
      t451 = t425*dqmu2
      t452 = 16*t451
      t457 = -t436
      t460 = t428*dpmu3
      t461 = t428*dqmu3
      t462 = -t446
      t463 = t462*dlmu3
      t469 = 48*dlmu3
      t470 = 32*dppmu3
      t471 = t282+t469+t283-t470
      t473 = t436*dqmu2
      t479 = 48*dppmu1
      t480 = 64*dlmu1
      t483 = 48*t254
      t484 = 32*t142
      t495 = -dlmu3-dqmu3
      t496 = 16*t495
      t498 = -t495
      t499 = 8*t498
      t504 = -t358-t194-t348
      t508 = -t499
      t510 = -t496
      t512 = 8*t250
      t513 = 8*t252
      t518 = 4*dppmu1
      t520 = 4*dqmu1
      t521 = 8*dpmu1
      t522 = -t518+4*dlmu1+t520-t521
      t525 = 8*t45
      t526 = 16*t47
      t527 = 8*t71
      t528 = t527*dqmu2
      t531 = t499*dmu1mu2
      t532 = dmu1mu3*dpmu2
      t533 = 8*t532
      t534 = 8*t408
      t535 = t531+t533-t534
      t538 = 6*q1s
      t539 = 6*q2s
      t540 = 12*q1q2
      t543 = 8*dlmu3
      t544 = 8*dqmu3
      t545 = 4*dppmu3
      t546 = 4*dpmu3
      t547 = -t543-t544+t545+t546
      t555 = 4*dlmu2
      t556 = 8*dpmu2
      t557 = -4*dqmu2+4*dppmu2-t555+t556
      t560 = (t496*dpmu2+t499*dlmu2+(-8*dlmu3-8*dqmu3+8*dppmu3+8*dpmu3)*
     #dqmu2+t504*dmu2mu3)*gp_mu1+(t508*dqmu1+t510*dpmu1-t512-t513-t504*d
     #mu1mu3)*gp_mu2+t522*gp_mu3_l_mu2+(t504*dmu1mu2+t525-t526+t528)*gp_
     #mu3+t535*gp_q+t535*gp_l+(t538+t539+t540)*gp_mu3_mu2_mu1+t547*gp_l_
     #mu2_mu1-t547*gp_mu2_q_mu1-t522*gp_mu3_mu2_q+t557*gp_mu3_l_mu1+t557
     #*gp_mu3_q_mu1
      t563 = -t28-t250
      t564 = 16*t563
      t565 = t564*dppmu2
      t567 = -t564*dpmu2
      t568 = 48*dqmu3
      t569 = t284+t568-t63+t61
      t571 = 48*t7
      t574 = 16*t417
      t577 = 48*t250
      t593 = t310*dlmu3
      t595 = t310*dpmu3
      t596 = t310*dqmu3
      t597 = t102+t87
      t599 = t306*dppmu3
      t605 = dqmu1+dlmu1-dppmu1+dpmu1
      t606 = 8*t605
      t608 = -t606
      t610 = -dqmu2-dlmu2
      t611 = 8*t610
      t613 = -dpmu1+dppmu1-dlmu1
      t614 = 8*t613
      t615 = t614*dlmu2
      t616 = 8*t47
      t617 = -t527
      t618 = t617*dqmu2
      t619 = t606*dppmu2+t608*dpmu2+t611*dqmu1+t615-t616+t618
      t621 = 8*t28
      t622 = 8*t425
      t623 = 8*t248
      t624 = 8*dppmu1
      t625 = 16*dlmu1
      t626 = -t624+t625
      t628 = 8*t254
      t632 = 16*q2s+16*q1q2+16*pq2
      t634 = dppmu3*q2s
      t635 = 8*t634
      t636 = t147+t149+t632*dlmu3-t635
      t640 = q2s*dppmu1
      t641 = 8*t640
      t642 = t632*dlmu1-t641
      t645 = 24*q2s
      t648 = -t632
      t655 = 8*t142*q2s
      t658 = q2s*dqmu1
      t659 = 8*t658
      t660 = q2s*dpmu1
      t661 = 8*t660
      t662 = t101+t94+t162+t86+t88
      t666 = -t662
      t667 = t666*dlmu1
      t671 = q2s*dqmu2
      t672 = 8*t671
      t682 = 8*t660*dqmu2
      t687 = t426+t318+t255
      t689 = 32*t425
      t696 = 64*dlmu3
      t701 = 32*t451
      t707 = 32*t131
      t708 = t707*dqmu1
      t709 = t707*dpmu1
      t710 = t72*dpmu3
      t713 = 48*dpmu3
      t714 = 80*dlmu3
      t715 = -t568-t713-t714+t63
      t718 = -t480+t445
      t719 = t718*dpmu3
      t720 = t718*dqmu3
      t723 = (t479-80*dlmu1)*dlmu3
      t726 = 64*dpmu3
      t727 = -t714-t62-t726+t63
      t729 = t715*dqmu2
      t733 = t462*dpmu3
      t734 = t462*dqmu3
      t742 = 8*dpmu3
      t743 = -t64+t339-t742-t544
      t745 = -t743
      t747 = 24*dlmu3
      t761 = 4*dppmu2*gp_mu3_q_mu1
      t765 = 4*dppmu3*gp_mu2_q_mu1
      t766 = dppmu2*dpmu1
      t767 = dppmu1*dqmu2
      t768 = dppmu1*dpmu2
      t775 = 8*t134*gp_mu2
      t777 = 4*dppmu1*gp_mu3_mu2_q
      t778 = dqmu2*dppmu3
      t784 = 16*t416
      t785 = -t338+t63-t64-t282
      t789 = -t785
      t793 = 32*t71
      t794 = t793*dpmu3
      t795 = t72*dlmu3
      t798 = dqmu2*dpmu3
      t799 = dlmu2*dpmu3
      t800 = t798+t799
      t801 = 16*t800
      t803 = t789*dpmu1
      t808 = t710*dqmu2
      t809 = (-t784+t785*dpmu1+t316)*dppmu2+(t789*dqmu1+(t713-t470+t283+
     #t284)*dpmu1+t794+t73+t795+t275-t255)*dpmu2+t801*dqmu1+(t803+t710)*
     #dlmu2+t789*dqmu2*dpmu1+t808
      t815 = -t282+t470-t568-t284
      t817 = -t283+t63-t64-t61
      t819 = -t793
      t820 = t819*dqmu3
      t823 = -t815
      t825 = -t817
      t826 = t825*dpmu1
      t827 = t793*dqmu3
      t831 = 48*dppmu3
      t838 = t823*dqmu2
      t840 = 48*t71
      t842 = t793*dlmu3
      t851 = t825*dqmu1+t574+t73
      t858 = t825*dqmu2
      t871 = t54*gp_mu2
      t874 = dqmu3*dmu1mu2
      t876 = -16*t874-16*t34-16*t20
      t878 = t45*gp_mu3
      s2 = (t2+4*t3+4*t5+(8*t7+8*t8)*gp_mu1+8*t12*gp_mu3)*Dij_1345(1,2)+
     #(-12*t18-8*t21+8*t23*gp_mu3-12*t26+8*t28*gp_mu2+(-4*ppq1+4*pq1)*gp
     #_mu3_mu2_mu1-t36+t39)*Dij_2345(1,1)+(4*t18+t44+(8*t45+8*t46+8*t47)
     #*gp_mu3)*Dij_1234(1,2)+((-t56+t58+(t60+(t61+t62-t63+t64)*dqmu2)*dq
     #mu1+t70+t73*dqmu2)*gp_q+t79)*Eij_12345(4,10)
      s1 = s2+t395*Eij_12345(2,4)+(t403*gp_q-t79)*Eij_12345(4,9)+(t410-t
     #412+t415+t419*gp_mu2)*Dij_1235(2,4)+(((-t424-t426-t249-t251+t429+t
     #255)*dppmu2+t434*dpmu2+(t437+t438)*dqmu1+(t441+t442+t443+t447+t275
     #-t448)*dlmu2+t452+(t249+t251+t433-t255)*dqmu2)*gp_q+((t457*dqmu1+t
     #457*dpmu1+t460+t461+t463+t448-t275)*dppmu2+(t436*dqmu1+t441+t442+t
     #443+t447+t275-t448)*dpmu2+(t471*dlmu2+t473)*dqmu1+(t471*dpmu1+t446
     #*dpmu3+t446*dqmu3+(-t479+t480)*dlmu3-t483+t484)*dlmu2+t473*dpmu1+(
     #t442+t443+t447+t275-t448)*dqmu2)*gp_l)*Eij_12345(4,16)+t560*Dij_12
     #45(1,2)
      s3 = s1+((t565+t567+(t569*dlmu2+t571)*dqmu1+(t574+t443)*dlmu2+t452
     #+(t249+t577+t433-t255)*dqmu2)*gp_q+(-t56+t58+(32*t59+t569*dqmu2)*d
     #qmu1+t70+t443*dqmu2)*gp_l)*Eij_12345(4,28)
      s5 = s3
      s8 = (t309*dppmu2+t593*dpmu2+(t595+t596+t597*dlmu3+t599)*dlmu2+t59
     #3*dqmu2)*gp_mu1+t619*gp_mu3_l_q+(t621+t622+t623+t512+t626*dlmu3-t6
     #28)*gp_l_mu2_q+(t636*dqmu1+t636*dpmu1+t642*dpmu3+t642*dqmu3+((t645
     #+t160+t156)*dlmu1+t648*dppmu1)*dlmu3+t648*dppmu3*dlmu1+t655)*gp_mu
     #2
      s9 = s8+((t659+t661+t662*dlmu1-t641)*dppmu2+(-t659-t661+t667+t641)
     #*dpmu2+(t666*dlmu2-t672)*dqmu1+(t666*dpmu1+(-t160-t156-t87-t645-t1
     #02)*dlmu1+t662*dppmu1)*dlmu2-t682+(t667+t641)*dqmu2)*gp_mu3
      s7 = s9+(t687*dppmu2+(-t424-t689-t249-t251-t253)*dpmu2+((-t284-t33
     #8-t61)*dlmu2-t438)*dqmu1+((-t283-t696-t282+t63)*dpmu1-t249-t251-t2
     #53+t275-t255)*dlmu2-t701+(-16*t248-16*t250-16*t252)*dqmu2)*gp_q+((
     #t708+t709+t710+t73-t275+t448)*dppmu2+(t715*dqmu1+t715*dpmu1+t719+t
     #720+t723+t255)*dpmu2+(t727*dlmu2+t729)*dqmu1+(t727*dpmu1+t733+t734
     #+t463+t275-t255)*dlmu2+t729*dpmu1+(t719+t720+t723+t255)*dqmu2)*gp_
     #l+(t743*dppmu2+t745*dpmu2+(t747+t61-t63+t338)*dlmu2+t745*dqmu2)*gp
     #_l_q_mu1
      s8 = Eij_12345(3,16)
      s6 = s7*s8
      s4 = s5+s6
      s2 = s4+(4*t203*gp_mu3_mu2_mu1-4*dppmu1*gp_mu3_l_mu2+t761+4*dppmu2
     #*gp_mu3_l_mu1-t765+(8*t766-8*t767-8*t768)*gp_mu3+4*dppmu3*gp_l_mu2
     #_mu1-t775+t777+(8*t342+8*t778)*gp_mu1)*D0_1245+(t809*gp_q+t809*gp_
     #l)*Eij_12345(4,14)
      t885 = s2+(((t815*dqmu1+t817*dpmu1+t316+t820+t318+t255-t275)*dppmu
     #2+(t823*dqmu1+t826+t710+t827+t795+t275-t255)*dpmu2+(t823*dlmu2+(t6
     #2+t469+t713-t831)*dqmu2)*dqmu1+(t826+t710+t827+t795+t275-t255)*dlm
     #u2+t838*dpmu1+(t794+t840*dqmu3+t842+t484-t448)*dqmu2)*gp_q+((t817*
     #dqmu1-t574+t317)*dppmu2+t851*dpmu2+(t825*dlmu2+t838)*dqmu1+(t574+t
     #73)*dlmu2+t858*dpmu1+(t710+t827+t795+t275-t255)*dqmu2)*gp_l)*Eij_1
     #2345(4,15)-2*gp_mu3_mu2_mu1*C0_135-2*gp_mu3_mu2_mu1*C0_345+(-16*t8
     #71-16*t38+t876*gp_q-16*t878)*Eij_12345(4,37)-2*C0_124*gp_mu3_mu2_m
     #u1
      t886 = dppmu2*dpmu3
      t888 = 16*t886*dpmu1
      t889 = 32*t416
      t894 = t47*dpmu3
      t896 = t12*dpmu3
      t897 = 16*t896
      t901 = -16*t418
      t903 = 16*t894
      t913 = t888+(-t784+(-t338+t63-t713-t64)*dpmu1+t316)*dpmu2-t897-t90
     #3
      t919 = t252*dqmu2
      t921 = t28*dlmu2
      t923 = -16*t919-16*t250*dlmu2-16*t921
      t929 = pq1*dlmu3
      t931 = pq1*dpmu3
      t932 = t931*dlmu2
      t936 = dpmu1*dpmu2
      t938 = 8*t936*gp_mu3_l_q
      t939 = pq2*dlmu3
      t940 = 16*t939
      t941 = t940+t147
      t943 = dlmu1*pq2
      t948 = -pq1-pq2
      t949 = 16*t948
      t950 = t949*dlmu1
      t951 = -t661+t950
      t953 = t949*dpmu1
      t957 = dpmu2*dpmu3
      t962 = (-t338-t63-t64-t713)*dpmu1
      t968 = 8*t957*gp_l_q_mu1
      t972 = pq1*dmu2mu3
      t973 = dpmu2*dqmu3
      t977 = 8*t416
      t981 = pq1*dmu1mu3
      t985 = 8*t766
      t988 = dpmu1-dqmu1
      t989 = 8*t988
      t993 = 12*t411
      t994 = gp_mu3_q_mu1*dpmu2
      t995 = 12*t994
      t996 = dpmu3*dmu1mu2
      t998 = 8*t532-8*t996-8*t408
      t1000 = dpmu3*gp_mu2_q_mu1
      t1001 = 12*t1000
      t1007 = t248+t425
      t1008 = 16*t1007
      t1022 = 8*t8
      t1023 = dpmu3+dqmu3
      t1024 = 8*t1023
      t1026 = 8*t7
      t1031 = t957*gp_mu1
      t1032 = dpmu1*dpmu3
      t1033 = t1032*gp_mu2
      t1034 = gp_mu3*dpmu1
      t1035 = t1034*dpmu2
      t1040 = q1s*dpmu3
      t1041 = 8*t1040
      t1042 = pq1*dqmu3
      t1043 = 16*t1042
      t1044 = t1041+t1043
      t1046 = pq1*dqmu2
      t1048 = 16*t1046*dpmu3
      t1052 = 8*t1032*gp_l_mu2_q
      t1053 = pq2*dpmu3
      t1055 = 16*t1053*dqmu1
      t1056 = pq2*dqmu3
      t1057 = q1q2*dpmu3
      t1059 = 16*t1056+16*t1057
      t1064 = t111*dpmu1
      t1065 = t949*dqmu1+t1064
      t1067 = t949*dqmu2
      t1083 = q1s*dlmu3
      t1087 = 8*t1083*dpmu2
      t1089 = -16*pq1-16*q1s
      t1090 = t1089*dqmu3
      t1091 = 16*t1083
      t1092 = dppmu3*q1s
      t1093 = 8*t1092
      t1101 = 8*dlmu1-8*dqmu1
      t1105 = 16*t45
      t1106 = -t521+t624-t625
      t1108 = 8*t767
      t1115 = 16*t148
      t1116 = -t162-t160-t86
      t1120 = q1q2*dlmu3
      t1121 = 16*t1120
      t1122 = -t149-t1121
      t1124 = q1q2*dlmu1
      t1126 = 16*t1124*dpmu3
      t1131 = q1q2*dppmu1
      t1132 = 16*t1131
      t1136 = 16*t1124*dppmu3
      t1139 = t111*dlmu1
      t1140 = -t659+t1139
      t1142 = -t111
      t1143 = t1142*dlmu1
      t1146 = t86+t94+t162+t102+t160
      t1148 = 16*t671
      t1151 = t1142*dpmu1
      t1152 = t160+t102
      t1154 = t111*dppmu1
      t1165 = 32*t28
      t1166 = 32*t248
      t1167 = 32*t250
      t1168 = t444-t137
      t1169 = t1168*dlmu3
      t1174 = 32*t7
      t1179 = t1168*dpmu3
      t1180 = -t718
      t1185 = 32*t252
      t1194 = t469+t282+t62
      t1196 = 32*t417
      t1201 = 96*dlmu3
      t1202 = 96*dqmu3
      t1207 = t12*dqmu3
      t1208 = 32*t1207
      t1216 = dppmu2*dqmu3
      t1217 = 8*t1216
      t1218 = 8*t973
      t1219 = t339-t338-t64-t742
      t1226 = -dlmu3+dppmu3-dpmu3-dqmu3
      t1227 = 8*t1226
      t1229 = -t1227
      t1230 = t1229*dpmu2
      t1232 = t1229*dqmu2
      t1237 = 16*t8
      t1238 = -t63-t284-t338-t61
      t1239 = t1238*dlmu2
      t1242 = dmu1mu2*gp_mu3_l_q
      t1243 = 8*t1242
      t1245 = 8*dmu1mu3*q2s
      t1248 = dlmu1*gp_mu3_l_mu2
      t1252 = dlmu1*dpmu2
      t1253 = 16*t1252
      t1255 = 16*dpmu1
      t1258 = dlmu1*dqmu2
      t1259 = 16*t1258
      t1263 = dlmu2*dmu1mu3
      t1264 = dlmu1*dmu2mu3
      t1266 = 16*dlmu3*dmu1mu2+16*t1263+16*t1264
      t1267 = t1266*gp_q
      t1269 = 16*t20
      t1270 = 16*t408
      t1271 = 16*t532
      t1273 = 16*t34
      t1276 = dlmu3*gp_mu2_q_mu1
      t1281 = dmu2mu3*gp_l_q_mu1
      t1282 = 8*t1281
      t1283 = dlmu1*gp_mu3_mu2_q
      t1285 = (-t1237+t1239-t438)*gp_mu1-t1243+(-t251-t424-t255-t249+t12
     #45-t426)*gp_mu2+16*t1248+(-8*q2s*dmu1mu2-t1253-32*t23+(-t1255-t444
     #)*dlmu2-t1259)*gp_mu3+t1267+(t268*dmu1mu2-t1269-t1270-t1271+16*t12
     #64-t1273)*gp_l+8*t1276+16*t1+8*t5+8*t3+t1282-8*t1283
      t1290 = 8*t1083
      t1293 = -t1044
      t1294 = t1293*dlmu2
      t1298 = 16*t929
      t1299 = pq1*dppmu3
      t1300 = 16*t1299
      t1305 = 16*t46
      t1308 = -dqmu3-dppmu3
      t1309 = 8*t1308
      t1313 = 32*t1056
      t1314 = pq2*dppmu3
      t1315 = 16*t1314
      t1320 = q1q2*dppmu3
      t1321 = 16*t1320
      t1324 = q1s*dlmu1
      t1327 = pq1*dlmu1
      t1333 = 24*pq1*dqmu1
      t1334 = q1s*dpmu1
      t1335 = 12*t1334
      t1338 = q1q2*dpmu1
      t1339 = pq2*dqmu1
      t1346 = q1s*dppmu1
      t1347 = 8*t1346
      t1350 = dlmu2*pq1
      t1352 = -t94+t156
      t1361 = pq1*dppmu1
      t1363 = -16*t1327+16*t1361
      t1368 = -dppmu3-dlmu3
      t1369 = 16*t1368
      t1373 = -t101-t88-t94
      t1374 = t1373*dmu1mu3
      t1377 = pq2*dmu2mu3
      t1379 = -16*t799+16*t1377
      t1381 = -t282-t63-t338
      t1382 = t1381*dqmu2
      t1383 = -t86-t89
      t1384 = t1383*dmu2mu3
      t1387 = pq2*dqmu2
      t1393 = 16*t54
      t1394 = q1s*dmu1mu3
      t1395 = 8*t1394
      t1401 = 16*t798
      t1402 = dmu2mu3*q1s
      t1403 = 8*t1402
      t1404 = t1401-t1403
      t1413 = 24*t1042
      t1414 = 12*t1040
      t1417 = -t94-t101-t359-t358
      t1419 = 24*t1387
      t1422 = 16*t1057
      t1423 = 24*t1056
      t1426 = q1s*dpmu2
      t1427 = 12*t1426
      t1428 = 24*t1046
      t1431 = 24*t1339
      t1432 = t194-t348+t191+t86
      s1 = (((t103+t101+t88)*dqmu3-t1290-t1093)*dpmu2+t1294+((16*pq1-16*
     #q1q2+16*q1s)*dpmu3+t1043-t1298-t1300)*dqmu2)*gp_mu1+(t1305-t616-t9
     #85)*gp_mu3_l_q+(t784+t1309*dpmu1)*gp_l_mu2_q+((t261-t1313+t1298-t1
     #315)*dqmu1+((t89+t161)*dqmu3+t1290-t1321)*dpmu1+8*t1324*dpmu3+16*t
     #1327*dqmu3)*gp_mu2+(-t1333-t1335)*gp_mu3_l_mu2+((16*t1338+16*t1339
     #)*dppmu2+((-t94-t101-t89+t93+t206-t86)*dqmu1-8*t1324+t1347)*dpmu2+
     #(16*t1350+t1352*dqmu2)*dqmu1+8*t1334*dlmu2+(t206-t89-t161-t102-t94
     #)*dqmu2*dpmu1+t1363*dqmu2)*gp_mu3+(t1059*dmu1mu2+(t1369*dqmu1+(-t4
     #70-t64-t283-t713)*dpmu1+t1374)*dpmu2+t1379*dqmu1-t897+(t1382+t1384
     #)*dpmu1-16*t1387*dmu1mu3)*gp_q
      t1436 = s1+(t1293*dmu1mu2+(t1393-t574+t1395)*dpmu2+(-16*t972-16*t7
     #98)*dqmu1+t1404*dpmu1+16*t981*dqmu2)*gp_l+(-8*t973+8*t798)*gp_l_q_
     #mu1+(t1413+t1414)*gp_l_mu2_mu1+(t1417*dpmu2-t1419)*gp_mu3_q_mu1+(t
     #1422+t1423)*gp_mu2_q_mu1+(t1427+t1428)*gp_mu3_l_mu1+(-t1431+t1432*
     #dpmu1)*gp_mu3_mu2_q
      t1440 = 16*t252*dppmu2
      t1441 = t252*dpmu2
      t1442 = 16*t1441
      t1453 = t1168*dqmu3
      t1474 = t470-t282-t568-t696
      t1478 = -t284-t61-t283+t63
      t1513 = -t12-t1252
      t1514 = 8*t1513
      t1516 = dlmu3*q2s
      t1518 = 8*t1516*dpmu1
      t1519 = dlmu1*q2s
      t1521 = 8*t1519*dpmu3
      t1527 = 8*t1519*dpmu2
      t1537 = t416*dlmu2
      t1538 = 16*t1537
      t1539 = -t64-t282-t63-t338
      t1540 = t1539*dpmu1
      t1542 = t798*dlmu1
      t1543 = 16*t1542
      t1546 = t799+t8
      s2 = ((-t888+(t889+(t64+t283+t713-t63)*dpmu1+t710)*dpmu2+32*t894+t
     #897)*gp_q+(t901*dpmu2+t903)*gp_l)*Eij_12345(4,24)+(t913*gp_q+t913*
     #gp_l)*Eij_12345(4,19)+(t923*gp_q+t403*gp_l)*Eij_12345(4,20)+((16*t
     #929*dpmu2+16*t932)*gp_mu1-t938+(t941*dpmu1+16*t943*dpmu3)*gp_mu2+(
     #t951*dpmu2+t953*dlmu2)*gp_mu3+16*t957*dlmu1*gp_q+((-t784+t962)*dpm
     #u2-t903)*gp_l+t968)*Eij_12345(3,6)
      s1 = s2+((8*t972+8*t798-8*t973)*gp_mu1+(t977+(-8*dpmu3+8*dqmu3+8*d
     #ppmu3)*dpmu1-8*t981)*gp_mu2+(-t985+8*pq1*dmu1mu2+t989*dpmu2-t616)*
     #gp_mu3+t993+t995+t998*gp_q-t1001-12*pq1*gp_mu3_mu2_mu1)*Dij_1235(1
     #,1)+((t1008*dpmu2+t897)*gp_q+(-t888+(t784+(t713-t63+t284+t338)*dpm
     #u1+t442)*dpmu2+t903+32*t896)*gp_l)*Eij_12345(4,25)+(t1022+t1024*dl
     #mu2+t1026)*gp_mu1*Dij_1345(2,4)+(-16*t1031-16*t1033-16*t1035)*Eij_
     #12345(4,36)+((t1044*dpmu2+t1048)*gp_mu1+t938+t1052+(t1055+t1059*dp
     #mu1)*gp_mu2+(t1065*dpmu2+t1067*dpmu1)*gp_mu3+(t962*dpmu2-t903-t897
     #)*gp_q)*Eij_12345(3,5)
      s3 = s1+(t765-t777+4*ppq1*gp_mu3_mu2_mu1-t761+t775-8*t766*gp_mu3)*
     #D0_1235
      s5 = s3
      s8 = (8*t1083*dppmu2-t1087+(-t1041+t1090-t1091+t1093)*dlmu2+t1089*
     #dlmu3*dqmu2)*gp_mu1+(t1101*dppmu2-t1101*dpmu2+t1105+t1106*dlmu2-t1
     #108+t616)*gp_mu3_l_q+(-t424-t622-t623-t251-t626*dlmu3+t628)*gp_l_m
     #u2_q+((-t147-t1115+t1116*dlmu3+t635)*dqmu1+t1122*dpmu1-t1126+(t111
     #6*dlmu1+t641)*dqmu3+(-32*t1124+t1132)*dlmu3+t1136)*gp_mu2
      s7 = s8+(t1140*dppmu2+(t659+t1143)*dpmu2+(t1146*dlmu2+t1148)*dqmu1
     #+(t1151+t1152*dlmu1+t1154)*dlmu2+t682+(t1146*dlmu1-t641)*dqmu2)*gp
     #_mu3+((-16*t254-16*t248-16*t252-16*t425)*dppmu2+(t1165+t689+t1166+
     #t1167+t1169)*dpmu2+((t696+t63+t282+t283)*dlmu2+t1174)*dqmu1+((t696
     #+t282+t568)*dpmu1+t1179+t1167+t1180*dlmu3+t255-t275)*dlmu2+48*t451
     #+(t1166+t1167+t1185+t255)*dqmu2)*gp_q+(((-t283-t63-t61)*dqmu1-t574
     #+t251)*dppmu2+(t1194*dqmu1+t1196+t443)*dpmu2+((t284+t282+t63+t62)*
     #dlmu2+(t726+t1201+t1202)*dqmu2)*dqmu1+t1208+t1194*dqmu2*dpmu1+(t44
     #2+t1180*dqmu3+t433+t255-t275)*dqmu2)*gp_l+(t1217-t1218-t60+t1219*d
     #qmu2)*gp_l_q_mu1
      s8 = Eij_12345(3,20)
      s6 = s7*s8
      s4 = s5+s6
      s2 = s4+(t1227*dppmu2+t1230+t1229*dlmu2+t1232)*gp_mu1*Dij_1345(2,3
     #)+t1285*Eij_12345(3,23)
      s3 = s2+t1436*Eij_12345(2,5)+((-t1440+t1442+32*t921+(t426+t249+t11
     #67+t1169-t255)*dlmu2+32*t919)*gp_q+(t565+t567+((-t63+t469+t61+t283
     #)*dlmu2+t1174)*dqmu1+(t574+t1453)*dlmu2+t452+(t249+t1167+t1169-t25
     #5)*dqmu2)*gp_l)*Eij_12345(4,31)
      s4 = s3+(t2+(8*t12+8*t1258+8*t1252)*gp_mu3+4*t1248)*Dij_1234(1,3)
      s5 = s4
      s7 = (((t1165+t426+t249+t1167+t433-t255)*dppmu2+(-t1165-t426-t249-
     #t1167+t429+t255)*dpmu2+(t1474*dlmu2-t571)*dqmu1+(t1478*dpmu1+t460+
     #t720+t463-t275+t448)*dlmu2-t701+(-t1166-t577+t718*dlmu3+t448)*dqmu
     #2)*gp_q+((-t1478*dqmu1+t574+t443)*dppmu2+(t1478*dqmu1-t574+t461)*d
     #pmu2+((-t469-t62+t470-t282)*dlmu2+t1474*dqmu2)*dqmu1+(-t1196+t734)
     #*dlmu2+t1478*dqmu2*dpmu1+(t460+t720+t463-t275+t448)*dqmu2)*gp_l)*E
     #ij_12345(4,34)
      s8 = (16*t1350*dlmu3*gp_mu1+t1514*gp_mu3_l_q+(t1518+t1521+16*t943*
     #dlmu3)*gp_mu2+(-t1527+t951*dlmu2)*gp_mu3+(16*t1441+16*t248*dlmu2)*
     #gp_q+((-t255-t253-t1166-t251-t689-t424)*dpmu2-t1538+t1540*dlmu2-t4
     #52-t1543)*gp_l+8*t1546*gp_l_q_mu1)*Eij_12345(3,11)
      s6 = s7+s8
      t1552 = s5+s6
      t1561 = t1219*dmu1mu2
      t1564 = -t545+t546+t543+t544
      t1571 = dlmu3*dppmu2
      t1576 = dlmu1*dppmu2
      t1577 = 16*t1576
      t1578 = 16*t23
      t1583 = dppmu2*dmu1mu3
      t1584 = 16*t1583
      t1586 = 32*t1263
      t1593 = dmu2mu3*gp_mu1
      t1595 = dmu1mu2*gp_mu3
      t1596 = 8*t1595
      t1597 = 4*gp_mu3_mu2_mu1
      t1598 = dmu1mu3*gp_mu2
      t1603 = -t1008
      t1604 = t1603*dppmu2
      t1608 = -t63+t282+t284+t283
      t1609 = t1608*dpmu1
      t1615 = -t901
      t1616 = t1615*dppmu2
      t1617 = t1608*dqmu1
      t1621 = 32*t800*dqmu1
      t1622 = t1608*dqmu2
      t1624 = t442*dqmu2
      t1639 = 8*t878
      t1640 = dmu1mu3*gp_l_mu2_q
      t1644 = dqmu2*gp_mu3_l_mu1
      t1646 = dqmu3*gp_l_mu2_mu1
      t1648 = 4*t1276-4*t1281+(-12*dqmu2+t555)*gp_mu3_q_mu1+t44+t39+(4*p
     #pp-4*q1q2)*gp_mu3_mu2_mu1+4*t1283+t1639-4*t1640-10*gp_mu3_l_mu2_q_
     #mu1-4*t1242+4*t1644+4*t1646-t36
      t1651 = 8*t798
      t1656 = 4*t1000
      t1658 = 8*t996*gp_q
      t1660 = 8*t996*gp_l
      t1661 = dpmu3*gp_l_mu2_mu1
      t1662 = 4*t1661
      t1678 = dqmu3*q1s
      t1680 = 8*t1678*dpmu2
      t1681 = t1678*dlmu2
      t1682 = 8*t1681
      t1688 = dqmu1*dppmu2
      t1692 = -dppmu3+dqmu3-dpmu3
      t1693 = 8*t1692
      t1695 = 8*t417
      t1705 = 16*t1338*dqmu3
      t1711 = q1q2*dqmu1
      t1714 = t101-t88
      t1717 = q1s*dlmu2
      t1723 = -t1714
      t1727 = -8*t1324+8*t1346
      t1731 = q1q2*dqmu3
      t1734 = -dpmu3-dppmu3-dqmu3
      t1735 = 16*t1734
      t1743 = 16*t1207
      t1744 = -t63-t64-t61-t283
      t1745 = t1744*dqmu2
      t1768 = (t1680-t1682+(-t1041+t597*dqmu3-t1290-t1093)*dqmu2)*gp_mu1
     #+(8*t45-8*t1688-8*t46+8*t47)*gp_mu3_l_q+(t1693*dqmu1+t1695)*gp_l_m
     #u2_q+((-t1422+(t89+t156)*dqmu3+t1290-t1321)*dqmu1+8*t1324*dqmu3+t1
     #705)*gp_mu2-12*q1s*dqmu1*gp_mu3_l_mu2+(16*t1711*dppmu2+t1714*dqmu1
     #*dpmu2+(8*t1717+(t206-t89-t156-t88+t93)*dqmu2)*dqmu1+t1723*dqmu2*d
     #pmu1+t1727*dqmu2)*gp_mu3+(16*t1731*dmu1mu2+(t1735*dqmu1-t574)*dpmu
     #2+(-t60+(-t64-t470-t61-t338)*dqmu2+t1384)*dqmu1-t1743+t1745*dpmu1+
     #t1374*dqmu2)*gp_q+(-8*t1678*dmu1mu2-t58+t1404*dqmu1+8*t1394*dqmu2)
     #*gp_l+12*t1678*gp_l_mu2_mu1+t1417*dqmu2*gp_mu3_q_mu1+16*t1731*gp_m
     #u2_q_mu1+12*q1s*dqmu2*gp_mu3_l_mu1+t1432*dqmu1*gp_mu3_mu2_q
      t1773 = 16*t799
      t1778 = 16*t766
      t1780 = 32*dpmu1
      t1783 = 16*t12
      t1787 = 16*t996+16*t532+16*t408
      t1793 = t252*gp_mu2
      t1795 = dlmu2*dlmu3
      t1798 = -t1266
      t1800 = dlmu1*dlmu2
      t1808 = t149+(t88+t87-t101)*dlmu3
      t1827 = 8*dppmu1+8*dpmu1
      t1829 = 8*dqmu1
      t1839 = 16*t146
      t1840 = -t101-t86+t89
      t1850 = t1840*dlmu1-t641
      t1854 = 64*pq2
      t1857 = t160+t89+t156
      t1867 = t101+t191+t86+t194
      t1869 = t1867*dpmu1
      t1872 = -t1867
      t1876 = 16*t660
      t1884 = (t86-t89+t101-t93+t94-t206+t88)*dlmu1
      t1912 = -t1369
      t1914 = dmu2mu3*q2s
      t1915 = 8*t1914
      t1932 = -dpmu3-dppmu3
      t1933 = 16*t1932
      t1935 = t306*dmu1mu3
      t1939 = t282+t283+t469
      t1943 = t310*dmu1mu3
      t1949 = t260*dmu2mu3
      t1971 = -t543+t63-t544
      t1979 = 8*t1057
      t1980 = 8*t1731
      t1983 = 8*t1320
      t1986 = q2s*dppmu2
      t1988 = q2s*dpmu2
      t1989 = 12*t1988
      t1992 = 12*t671
      t1995 = 12*t146
      t1996 = 12*t148
      t2002 = t348-t358-t94-t359
      t2004 = -t2002
      t2006 = 8*pq1
      t2012 = 12*t658
      t2013 = 12*t660
      t2014 = 8*pq2
      s2 = (t1808*dppmu2+(-t149-t122*dlmu3)*dpmu2+(t85+(-t162-t87-t88-t8
     #6)*dqmu3+32*t83*dlmu3+t96)*dlmu2+(t147+(t86-t103+t101-t102+t162)*d
     #lmu3-t635)*dqmu2)*gp_mu1+(t1827*dppmu2+(t1829+t136-t137)*dpmu2+t52
     #5-t1827*dlmu2+t138*dqmu2)*gp_mu3_l_q+(-t621-t623-t512+t318)*gp_l_m
     #u2_q
      s3 = s2+((t1839+t1115+t1840*dlmu3)*dqmu1+((-t160-t161-t645)*dlmu3+
     #16*t634)*dpmu1+t1850*dpmu3+t1850*dqmu3+((-t162-64*q1q2-t1854-t87-t
     #102)*dlmu1+t1857*dppmu1)*dlmu3-16*t386*dppmu3*dlmu1-t655)*gp_mu2
      s1 = s3+(t1867*dqmu1+t1869+(t192+t193+t160+t156+t89+t165)*dlmu1+t1
     #872*dppmu1)*gp_mu3_l_mu2+((-t1876+(-16*q1q2+16*ppq1+16*ppp+16*ppq2
     #-16*q2s-16*pq2)*dlmu1+t641)*dppmu2+(-16*t658+t1884+t641)*dpmu2+((-
     #16*q1q2-16*ppp-16*pq2-16*ppq2-16*q2s-16*ppq1)*dlmu2-t1148)*dqmu1+(
     #(16*q2s+16*pq1-16*ppq2+16*pq2-16*ppp)*dpmu1+32*t203*dlmu1+(-t88+t2
     #06-t94+t93+t165)*dppmu1)*dlmu2+(t1884+t641)*dqmu2)*gp_mu3+((-t426+
     #t253-t1245-t255)*dppmu2+(-t147-t149+t648*dlmu3+t635)*dmu1mu2+(t689
     #+t318+t448+t1245)*dpmu2+(t1912*dlmu2-t1915)*dqmu1+((t61+t469+t63+t
     #338)*dpmu1-t292+t448+t662*dmu1mu3)*dlmu2+(-t1915+t1174)*dpmu1+(t31
     #8+t448+t1245)*dqmu2+(-8*t1519+8*t640)*dmu2mu3)*gp_q
      s2 = s1+((-t424+t1933*dpmu1+t1935)*dppmu2+t593*dmu1mu2+(t1939*dqmu
     #1+(t696+t283+t713+t63)*dpmu1+t442+t443+t433+t448+t1943-t275)*dpmu2
     #+(-t1238*dlmu2+t1939*dqmu2+t1949)*dqmu1+((t284+t338+t470+t282)*dpm
     #u1+t448+t1943-t275)*dlmu2+((t469+t713+t63+t283)*dqmu2+t1949)*dpmu1
     #+(t442+t443+t433+t448+t1943-t275)*dqmu2+((t88+t94+t160+t162+t156)*
     #dlmu1+t264*dppmu1)*dmu2mu3)*gp_l+((t64-t742-t339)*dppmu2+t1971*dpm
     #u2+(-t64+t339-t544)*dlmu2+(-t544+t63-t64)*dqmu2)*gp_l_q_mu1+(-t197
     #9-t1980+(-t192-t101-t193)*dlmu3+t1983)*gp_l_mu2_mu1
      t2021 = s2+(-12*t1986+t1989+(t359+t645+t357+t358+t101+t94)*dlmu2+t
     #1992)*gp_mu3_q_mu1+(-t1995-t1996+(-t101-t645-t357)*dlmu3+12*t634)*
     #gp_mu2_q_mu1+(t2002*dppmu2+t2004*dpmu2+(t2006+t93-t358-t101)*dlmu2
     #+t2004*dqmu2)*gp_mu3_l_mu1+(t2012+t2013+(t348+t2014+20*q2s-t191)*d
     #lmu1-12*t640)*gp_mu3_mu2_q
      t2024 = -t23-t1258
      t2034 = 16*t1056
      t2037 = 16*t932
      t2039 = 16*pq1-16*pq2
      t2050 = dpmu1*pq1
      t2053 = dpmu1*pq2
      t2058 = 16*t1361
      t2063 = -t2039
      t2069 = pq2*dmu1mu3
      t2086 = dpmu2*pq1
      t2093 = ((t100+t2034-t1298-t1300)*dpmu2-t2037+t2039*dpmu3*dqmu2)*g
     #p_mu1+(-t1055+(t158-t2034+t1298-t1315)*dpmu1+16*t1327*dpmu3)*gp_mu
     #2-24*t2050*gp_mu3_l_mu2+(16*t2053*dppmu2+(16*t1339+t218-16*t1327+t
     #2058)*dpmu2+16*t2050*dlmu2+t2063*dqmu2*dpmu1)*gp_mu3+(16*t1053*dmu
     #1mu2-16*t2069*dpmu2+16*t2053*dmu2mu3)*gp_q+(-16*t931*dmu1mu2+16*t9
     #81*dpmu2-16*t2050*dmu2mu3)*gp_l+24*t931*gp_l_mu2_mu1-24*pq2*dpmu2*
     #gp_mu3_q_mu1+24*t2086*gp_mu3_l_mu1+24*t1053*gp_mu2_q_mu1-24*t2053*
     #gp_mu3_mu2_q
      t2098 = -t132
      t2100 = t348+t2014+t194+t358+t2006
      t2101 = t2100*dmu2mu3
      t2104 = t1227*dpmu1
      t2110 = -4*t605
      t2112 = t2100*dmu1mu2
      t2113 = dppmu1-dlmu1-dqmu1
      t2114 = 8*t2113
      t2116 = 8*t12
      t2119 = t2098*dmu1mu2
      t2122 = 12*pq1
      t2126 = -t545+t742+t543+t544
      t2130 = dlmu2+dpmu2+dqmu2-dppmu2
      t2131 = 4*t2130
      t2136 = (8*dppmu2*dppmu3+t1230+t2098*dlmu2+t1232+t2101)*gp_mu1+(t1
     #33+t2104+t623+t512+t513-t2100*dmu1mu3-t143)*gp_mu2+t2110*gp_mu3_l_
     #mu2+(-t985+t2112+t2114*dpmu2-t525+t2116+t618)*gp_mu3+t2119*gp_q+t2
     #119*gp_l+(-t539-t2122-12*pq2-t540-t538)*gp_mu3_mu2_mu1+t2126*gp_l_
     #mu2_mu1-t2126*gp_mu2_q_mu1+t2131*gp_mu3_l_mu1-t2110*gp_mu3_mu2_q+t
     #2131*gp_mu3_q_mu1
      t2141 = -t1089
      t2143 = t1041+t2141*dqmu3+t1290-t1093
      t2147 = 24*q1s
      t2151 = t1089*dppmu3
      t2164 = 24*dqmu3
      t2167 = t742+t338+t543-t339
      t2169 = t527*dpmu3
      t2173 = t86+t160+t89
      t2175 = 48*q1q2
      t2179 = -t2173
      t2187 = 16*t1124-16*t1131
      t2189 = t2173*dlmu1
      t2195 = 16*t1131*dppmu3
      t2198 = t102+t94+t89+t160+t86
      t2202 = -t2198
      t2204 = t1142*dppmu1
      t2225 = (t568+t282+t284)*dpmu1
      t2228 = 80*dqmu3
      t2229 = -t726-t696-t2228+t63
      t2232 = (-t2228-t469-t713+t63)*dpmu1
      t2233 = t819*dpmu3
      t2234 = t819*dlmu3
      t2238 = 80*dpmu3
      t2255 = -t284-t62-t282+t63
      t2271 = 8*t59
      s2 = ((t499*dpmu2-t1219*dqmu2)*gp_mu1+t508*dpmu1*gp_mu2+t1561*gp_q
     #+t1561*gp_l+t1564*gp_l_mu2_mu1-t1564*gp_mu2_q_mu1)*Dij_1245(2,6)+(
     #(-16*t1571+t1237+t437+t438)*gp_mu1+t434*gp_mu2+(-t1577+t1253+t1578
     #+(-t137+t1255+t427)*dlmu2+t1259)*gp_mu3+t1267+(-t1584+t436*dmu1mu2
     #+t1271+t1269+t1586+t1270+t1273+t432*dmu2mu3)*gp_l)*Eij_12345(4,45)
     #+(8*t1593+t1596-t1597-8*t1598)*Dij_1235(2,7)+((t1604+(t1165+t689+t
     #1166+t1167+t433-t255)*dpmu2+32*t1537+(t1609+t442)*dlmu2+32*t1542+t
     #701)*gp_q+(t1616+(t1617+t1196+t443)*dpmu2+t1621+t1208+t1622*dpmu1+
     #t1624)*gp_l)*Eij_12345(4,35)
      s1 = s2+t1648*D0_2345+((t508*dpmu2-t1651)*gp_mu1+t499*dpmu1*gp_mu2
     #+t1656+t1658+t1660-t1662)*Dij_1245(2,4)+(4*t1598-4*t1593-4*t1595+6
     #*gp_mu3_mu2_mu1)*C0_245+2*gp_mu3_mu2_mu1*C0_234-t1514*gp_mu3*Dij_1
     #234(2,5)
      s2 = s1+t1768*Eij_12345(2,2)+((-16*t886+t789*dpmu2+t1401+t1773)*gp
     #_mu1+(t784+t803+t710)*gp_mu2+(-t1778+(16*dqmu1+t625+t1780-t137)*dp
     #mu2+t526+t1783)*gp_mu3+t1787*gp_q+t1787*gp_l)*Eij_12345(4,42)+(-16
     #*t1793-16*t1795*gp_mu1+t1798*gp_l-16*t1800*gp_mu3)*Eij_12345(4,38)
     #+t2021*Eij_12345(2,10)
      s3 = s2+(16*t2024*gp_mu3+t564*gp_mu2+t876*gp_l+t1798*gp_q+t399*gp_
     #mu1)*Eij_12345(4,43)+t2093*Eij_12345(2,1)
      s4 = s3+t2136*Dij_1245(1,3)
      s5 = s4
      s9 = ((-t1041+t1090-t1290+t1093)*dppmu2+t2143*dpmu2+t2143*dlmu2+(t
     #2141*dpmu3+(t2147+t87)*dqmu3+t2141*dlmu3+t2151)*dqmu2)*gp_mu1+(t60
     #8*dppmu2+t606*dpmu2-t611*dqmu1-t614*dlmu2+t616+t528)*gp_mu3_l_q+((
     #t61+t64+t2164-t63)*dqmu1+t2167*dpmu1+t2169+t73+t527*dlmu3-t628+t14
     #3)*gp_l_mu2_q+((t2173*dpmu3+(t156+t162+t2175)*dqmu3+t2173*dlmu3+t2
     #179*dppmu3)*dqmu1+(t1422+t2173*dqmu3+t1121-t1321)*dpmu1+t2187*dpmu
     #3+(t2189+t2179*dppmu1)*dqmu3+t2187*dlmu3-t1136+t2195)*gp_mu2
      s10 = s9+((t2198*dqmu1+t1151+t1143+t1154)*dppmu2+(t2202*dqmu1+t106
     #4+t1139+t2204)*dpmu2+(t2202*dlmu2+(-t162-t2147-t87-t2175-t156)*dqm
     #u2)*dqmu1+(t1064+t1139+t2204)*dlmu2+t2202*dqmu2*dpmu1+(t2202*dlmu1
     #+t2198*dppmu1)*dqmu2)*gp_mu3
      s8 = s10+((-t285*dqmu1+t2225+t710+t795+t255-t275)*dppmu2+(t2229*dq
     #mu1+t2232+t2233+t820+t2234)*dpmu2+(t2229*dlmu2+(-t2238-t714-t1202+
     #t63)*dqmu2)*dqmu1+(t2232+t2233+t820+t2234)*dlmu2+(-112*dqmu3+t470-
     #t2238-t714)*dqmu2*dpmu1+(t2233+t820+t2234-t255+t275)*dqmu2)*gp_q+(
     #(t244+t317)*dppmu2+(t2255*dqmu1-t1196+t317)*dpmu2+((-t283-t61-t64)
     #*dlmu2+(-t696-t1202-t726+t470)*dqmu2)*dqmu1-t1743+t2255*dqmu2*dpmu
     #1+(t316+t820+t318)*dqmu2)*gp_l+(-t1217+t1218+t2271+dqmu2*t2167)*gp
     #_l_q_mu1
      s9 = Eij_12345(3,15)
      s7 = s8*s9
      s8 = (-16*t1598-16*t1593-16*t1595)*Eij_12345(4,46)
      s6 = s7+s8
      t2282 = s5+s6
      t2285 = 16*t1800*dlmu3*gp_q
      t2286 = 16*t921
      t2288 = (-t137+t480)*dlmu3
      t2291 = 16*t919
      t2305 = t358+t2006
      t2307 = ppp*dpmu1
      t2314 = t310*dppmu1
      t2317 = ppp*dqmu2
      t2325 = 4*q1s*ppq2
      t2327 = 8*q1q2*ppq1
      t2329 = 4*q2s*ppq1
      t2331 = 8*pq1*ppq2
      t2333 = 8*pq2*ppq1
      t2336 = t348+t2014+t194
      t2345 = -t2336
      t2348 = -t2305
      t2355 = (t306*dpmu2*dppmu3+t306*dqmu2*dppmu3)*gp_mu1+t265*dpmu1*gp
     #_mu2+(t2305*dppmu1+8*t2307)*gp_mu3_l_mu2+(t260*dppmu2*dpmu1+(16*t2
     #307+t2314)*dpmu2+16*t2317*dpmu1+t310*dqmu2*dppmu1)*gp_mu3+(t2325-t
     #2327-t2329+t2331-t2333)*gp_mu3_mu2_mu1+t2336*dppmu1*gp_mu3_mu2_q+(
     #t2336*dppmu2+8*t2317+8*ppp*dpmu2)*gp_mu3_q_mu1+t2345*dppmu3*gp_mu2
     #_q_mu1+t2348*dppmu3*gp_l_mu2_mu1+4*ppp*gp_mu3_l_mu2_q_mu1+t2348*dp
     #pmu2*gp_mu3_l_mu1
      t2357 = dpmu3*gp_q
      t2394 = t99*dqmu3
      t2402 = 8*t1576
      t2414 = t157*dlmu3
      t2417 = -t1122
      t2418 = t2417*dpmu1
      t2419 = t157*dlmu1
      t2426 = 12*t1324
      t2432 = t1714*dlmu1
      t2466 = t1744*dqmu1
      t2469 = 16*t1308
      t2471 = t264*dmu2mu3
      t2475 = 8*t1263*q1s
      t2486 = 12*t1083
      t2494 = 12*t1717
      s1 = (t1087+(-t1041+t2394-t1091-t1093)*dlmu2+(t87-t101+t102)*dlmu3
     #*dqmu2)*gp_mu1+(-t2402+(-8*dlmu1-8*dqmu1)*dpmu2-t525+t2116+t616+(8
     #*dppmu1+8*dlmu1)*dqmu2)*gp_mu3_l_q+(-8*t254+8*t28+8*t425+8*t250-8*
     #t248)*gp_l_mu2_q+((-t147-t1115+t2414-t635)*dqmu1+t2418-t1126+t2419
     #*dqmu3+16*t1324*dlmu3-t1136)*gp_mu2+(t1872*dqmu1-t2426)*gp_mu3_l_m
     #u2+((t659+16*t1124)*dppmu2+(t659+t2432)*dpmu2+((16*ppq2+16*ppp+16*
     #q1q2+16*ppq1)*dlmu2+t1148)*dqmu1+(t1723*dpmu1+t1347)*dlmu2-t682+(-
     #t88+t93-t156-t101-t89+t206)*dlmu1*dqmu2)*gp_mu3+(t2417*dmu1mu2+(-1
     #6*t254-16*t425-16*t248)*dpmu2+(t1369*dlmu2+t1915)*dqmu1+(t1238*dpm
     #u1+t1374)*dlmu2-t701+(-t249-t253-t1245-t448)*dqmu2+t1383*dlmu1*dmu
     #2mu3)*gp_q
      t2501 = s1+(-8*t1083*dmu1mu2+(t2466-t574-t251)*dpmu2+(t2469*dlmu2+
     #t1745+t2471)*dqmu1+t2475+t279*dpmu1+(t249+t1935)*dqmu2-8*t1324*dmu
     #2mu3)*gp_l+(t1218+t2271+t1693*dqmu2)*gp_l_q_mu1+(t2486+t1980)*gp_l
     #_mu2_mu1+(t1417*dlmu2-t1992)*gp_mu3_q_mu1+(t1121+t1996)*gp_mu2_q_m
     #u1+(t2494+t2002*dqmu2)*gp_mu3_l_mu1+(-t2012+t1432*dlmu1)*gp_mu3_mu
     #2_q
      t2525 = dqmu1*gp_mu3_mu2_q
      t2533 = -t284-t338-t282+t63
      t2548 = 32*t799
      t2564 = 16*t972
      t2567 = 16*t2069
      t2570 = gp_mu3_l_mu2*dpmu1
      t2573 = 32*dqmu1
      t2576 = 32*t12
      t2582 = dpmu2*gp_mu3_l_mu1
      t2590 = (t1539*dpmu2-t1773-t1401+t2564)*gp_mu1+(-t784+t1540+t2567+
     #t249)*gp_mu2+16*t2570+(t949*dmu1mu2+(-t625-t1780-t2573)*dpmu2-t526
     #-t2576)*gp_mu3+16*t994+8*t1661-8*t411+8*t2582+(16*t532+16*t996)*gp
     #_q+8*t1000+16*t408*gp_l
      t2595 = 8*t287
      t2611 = q1s*gp_mu1
      t2618 = 8*t252*gp_l_mu2_q
      t2620 = 8*t1516*dqmu1
      t2624 = 8*t1519*dqmu3
      t2627 = dlmu2*q2s
      t2629 = 8*t2627*dqmu1
      t2632 = 8*t1519*dqmu2
      t2639 = t564*dpmu2
      t2653 = t2639-t452-t1743-t1543-t1538
      t2655 = t45*dpmu3
      t2657 = -16*t2655-16*t57-16*t69
      t2671 = t111*dqmu1
      t2674 = t111*dqmu2
      t2682 = -t801*dqmu1
      t2683 = -t63-t282-t283-t64
      t2685 = t2683*dqmu2*dpmu1
      t2699 = 8*t1688
      t2700 = 8*t46
      t2702 = 16*dqmu2
      t2703 = -8*dlmu2-t2702
      t2709 = t617*dqmu3
      t2716 = -t2187
      t2736 = t284+t282+t62
      t2760 = 8*t37*gp_l_q_mu1
      t2765 = 16*t1040
      t2766 = -t597
      t2776 = t624-t1255-t136-t1829
      t2780 = -t338-t61-t543+t339
      t2782 = t617*dpmu3
      t2785 = -t1857
      t2794 = pq2*dppmu1
      t2796 = -16*t943+16*t2794
      t2805 = -t949
      t2806 = t2805*dlmu2
      t2807 = -t948
      t2808 = 32*t2807
      t2815 = t2805*dlmu1
      t2816 = t949*dppmu1
      t2822 = t63+t726+t284+t283
      t2824 = 96*dpmu3
      t2838 = t901*dppmu2
      t2839 = t64+t283+t282
      t2843 = 32*t2655
      t2848 = -t798-t973
      t2849 = 8*t2848
      s2 = (t2285+(-t1440+t1442+t2286+(t426+t249+t251+t2288-t255)*dlmu2+
     #t2291)*gp_l)*Eij_12345(4,13)+t2355*E0_12345-16*t936*t2357*Eij_1234
     #5(4,5)+((t1442+t1008*dlmu2)*gp_q+(t1604+(t424+t689+t1166+t251+t116
     #9-t255)*dpmu2+t1538+((-t63+t469+t338+t282)*dpmu1+t1179)*dlmu2+t154
     #3+t452)*gp_l)*Eij_12345(4,30)
      s3 = s2+(-12*t1-8*t1263*gp_q+(-4*ppq2+4*pq2)*gp_mu3_mu2_mu1-8*t126
     #4*gp_l-12*t1248+(-8*t1258+8*t23+8*t1800)*gp_mu3+8*t1793+8*t59*gp_m
     #u1)*Dij_2345(1,2)+t2501*Eij_12345(2,8)
      s1 = s3+(-8*t1800*gp_mu3_l_q+8*t252*q2s*gp_mu2-8*t1519*dlmu2*gp_mu
     #3+t2285+(-t1442-t2286+(-16*t252-16*t248-16*t250-16*t254-16*t425)*d
     #lmu2-t2291)*gp_l+8*t1795*gp_l_q_mu1)*Eij_12345(3,3)+(t1639+8*t20*g
     #p_q-8*t871-4*t2525)*Dij_1235(2,2)+((t1008*dppmu2+(-t424-t689-t1166
     #-t251+t429+t255)*dpmu2-t1538+(t2533*dpmu1+t460)*dlmu2-t452-t1543)*
     #gp_q+((t784-t2533*dpmu1+t442)*dppmu2+(t2533*dqmu1+(-t713+t470-t283
     #-t696)*dpmu1+t719+t461+t463-t275+t448)*dpmu2+(-t1401-t2548)*dqmu1+
     #((-t469+t470-t283-t726)*dpmu1+t733)*dlmu2+t2533*dqmu2*dpmu1+t460*d
     #qmu2)*gp_l)*Eij_12345(4,33)
      s2 = s1+t2590*Eij_12345(3,21)+((t977+(t61+t544)*dpmu1-t2595)*gp_mu
     #2+((-t1829+t624-t1255)*dpmu2-t616)*gp_mu3-t410+t412)*Dij_1235(2,5)
     #+8*t1595*Dij_1234(2,7)+t619*gp_mu3*Dij_2345(2,3)
      s3 = s2+(8*t1795*t2611+(-8*t1258+8*t1800-8*t23)*gp_mu3_l_q+t2618+(
     #t2620+16*t1124*dlmu3+t2624)*gp_mu2+(-t2629+t1139*dlmu2-t2632)*gp_m
     #u3+(-t1442+(-t426-t249-t255-t1185)*dlmu2)*gp_q+(t2639+(t1744*dlmu2
     #-t1174)*dqmu1-t1743-t452+(-t255-t1167-t249-t253)*dqmu2)*gp_l-8*t39
     #8*gp_l_q_mu1)*Eij_12345(3,12)+(t2653*gp_q+t2657*gp_l)*Eij_12345(4,
     #26)
      s4 = s3+((t1680+t1044*dqmu2)*gp_mu1+t414*gp_mu3_l_q-t419*gp_l_mu2_
     #q+(t1059*dqmu1+t1705)*gp_mu2+(t2671*dpmu2+t1067*dqmu1+t2674*dpmu1)
     #*gp_mu3+((t1539*dqmu1-t1196)*dpmu2+t2682-t1743+t2685)*gp_q)*Eij_12
     #345(3,8)
      s5 = s4
      s9 = (8*t1678*dppmu2-t1680-t1682+(-t1041+(-t2147-t94)*dqmu3-t1290+
     #t1093)*dqmu2)*gp_mu1+(t2699-t2700+t2703*dqmu1-t616+t618)*gp_mu3_l_
     #q+((-t543+t339-t2164-t742)*dqmu1-t1695+t2709)*gp_l_mu2_q+((-t1422+
     #(-t2175-t89-t86)*dqmu3-t1121+t1321)*dqmu1-t1705+t2716*dqmu3)*gp_mu
     #2
      s8 = s9+(t2671*dppmu2+t1142*dqmu1*dpmu2+(t1142*dlmu2+(t89+t2147+t8
     #6+t94+t2175)*dqmu2)*dqmu1+t1142*dqmu2*dpmu1+(t1143+t1154)*dqmu2)*g
     #p_mu3+((t269-t574)*dppmu2+(t2736*dqmu1+t1196+t73)*dpmu2+(t2736*dlm
     #u2+(t726+t63+t1202+t696)*dqmu2)*dqmu1+(t1196+t73)*dlmu2+(t284+t282
     #+t2228)*dqmu2*dpmu1+(t710+t73+t795+t255-t275)*dqmu2)*gp_q+(t56+t58
     #+(t64+t568+t61)*dqmu2*dqmu1+t70)*gp_l-t2760
      s9 = Eij_12345(3,10)
      s7 = s8*s9
      s10 = (t1044*dppmu2+(-t2765+t2766*dqmu3-t1290+t1093)*dpmu2+t1294+(
     #t2766*dpmu3-32*t1042-t1298+t1300)*dqmu2)*gp_mu1+(t985+t2776*dpmu2-
     #t616-t2116)*gp_mu3_l_q+(-t784+t2780*dpmu1+t2782)*gp_l_mu2_q+((t278
     #5*dpmu3-t1313-t940+t1315)*dqmu1+(-32*t1057+t2785*dqmu3-t1121+t1321
     #)*dpmu1+t2716*dpmu3+t2796*dqmu3)*gp_mu2
      s9 = s10+(t1065*dppmu2+(t177*dqmu1+t1152*dpmu1+t1143+t1154)*dpmu2+
     #(t2806+t2808*dqmu2)*dqmu1+t1151*dlmu2+t177*dqmu2*dpmu1+(t2815+t281
     #6)*dqmu2)*gp_mu3+(t1540*dppmu2+(t2822*dqmu1+(t696+t1202+t2824)*dpm
     #u1+t794+t795-t275+t255)*dpmu2+t1621+((t726+t568+t284)*dpmu1+t710)*
     #dlmu2+(t2824+t469+t62)*dqmu2*dpmu1)*gp_q+(t2838+(t2839*dqmu1+t1196
     #)*dpmu2+t2843+t2839*dqmu2*dpmu1)*gp_l+t2849*gp_l_q_mu1
      s10 = Eij_12345(3,18)
      s8 = s9*s10
      s6 = s7+s8
      t2854 = s5+s6
      t2861 = t103+t88
      t2866 = (t94-t101+t88)*dpmu3
      t2867 = 32*t929
      t2882 = -t941*dqmu1
      t2883 = t101+t161+t89
      t2889 = 16*t943*dqmu3
      t2893 = 16*t943*dppmu3
      t2897 = 24*t1327
      t2907 = pq2*dlmu2
      t2914 = t2063*dlmu1
      t2937 = 16*t981
      t2950 = 24*t929
      t2953 = 24*t2907
      t2956 = 24*t939
      t2959 = t2002*dpmu2
      t2960 = 24*t1350
      t2963 = 24*t943
      s1 = ((t149+t2861*dlmu3)*dpmu2+(t2866+t2034-t2867-t1300)*dlmu2+(-t
     #147+t2039*dlmu3)*dqmu2)*gp_mu1+((t625+t624)*dpmu2-t616-t1783)*gp_m
     #u3_l_q+(-t622+t249)*gp_l_mu2_q+(t2882+(-t149+t2883*dlmu3-t635)*dpm
     #u1+t2189*dpmu3-t2889+32*t1327*dlmu3-t2893)*gp_mu2+(t1872*dpmu1-t28
     #97)*gp_mu3_l_mu2+((t661+16*t943)*dppmu2+(t659+(t93+t206-t94-t86-t1
     #60-t89)*dlmu1)*dpmu2+16*t2907*dqmu1+((t165-t86+t206-t94+t101-t88)*
     #dpmu1+t2058)*dlmu2+t682+t2914*dqmu2)*gp_mu3+(t941*dmu1mu2+(-t426+t
     #249-t253-t255-t1245)*dpmu2+(16*t425-16*t2069)*dlmu2+16*t1377*dlmu1
     #+8*t1914*dpmu1+t1543)*gp_q
      t2966 = s1+(-16*t929*dmu1mu2+(-t784+(-t284-t470-t338-t713)*dpmu1+t
     #251+t1935)*dpmu2+t2682+(t2469*dpmu1+t2937)*dlmu2+(t1382+t2471)*dpm
     #u1-t1543-16*t1327*dmu2mu3)*gp_l+(8*t1368*dpmu2+t1651+t1773)*gp_l_q
     #_mu1+(t1979+t2950)*gp_l_mu2_mu1+(-t2953-t1989)*gp_mu3_q_mu1+(t1995
     #+t2956)*gp_mu2_q_mu1+(t2959+t2960)*gp_mu3_l_mu1+(-t2963-t2013)*gp_
     #mu3_mu2_q
      t2970 = t1603*dpmu2-t897
      t2973 = t1615*dpmu2-t903
      t2982 = dlmu2-dqmu2
      t2985 = 8*t1258
      t3009 = (-16*t252-16*t248-16*t254-16*t250)*dqmu2
      t3012 = -t61-t568-t63-t64
      t3020 = t2348*dmu2mu3
      t3023 = 8*t1734
      t3031 = -dppmu1+dpmu1+dqmu1
      t3032 = 8*t3031
      t3034 = t2348*dmu1mu2
      t3057 = -t1226
      t3058 = 16*t3057
      t3060 = -t3058
      t3062 = t3060*dlmu2
      t3063 = t3060*dqmu2
      t3067 = t3060*dpmu1
      t3068 = t3060*dqmu1+t3067+t316+t317+t318+t255-t275
      t3070 = 16*t605
      t3072 = -t3070
      t3082 = 16*t1263
      t3083 = t315*dmu2mu3
      t3084 = t1584+t3060*dmu1mu2-t1271-t1269-t3082-t1270-t1273+t3083
      t3090 = t283-t63+t64+t282
      t3094 = -t568-t726-t284+t470
      t3100 = 48*t798
      t3103 = -t3090
      t3115 = 32*t798
      t3117 = (-t1773-t3115)*dqmu1
      t3120 = t316*dqmu2
      t3126 = -t1024
      t3129 = 8*t289
      t3135 = dppmu1*dmu2mu3
      t3160 = 4*t2570
      t3163 = (8*t936+8*t47)*gp_mu3
      t3165 = 4*dqmu2+4*dpmu2
      t3166 = t3165*gp_mu3_q_mu1
      t3167 = 2*gp_mu3_l_mu2_q_mu1
      t3209 = (q1s+2*pq1)*q2s
      t3217 = 8*dppmu1*dlmu2
      t3223 = t264*dlmu3
      t3227 = t310*dlmu1
      t3229 = t264*dlmu1
      t3238 = 8*ppp
      t3239 = -t358-t2006+t359+t3238
      t3243 = t260*dlmu1
      t3248 = ppp+ppq1
      t3249 = 16*t3248
      t3250 = t3249*dlmu1
      t3253 = ppq2+ppp
      t3270 = -dppmu3-dlmu3-dpmu3
      t3271 = 16*t3270
      t3272 = t3271*dpmu1
      t3276 = 16*t134*dlmu2
      t3289 = 4*ppq2
      t3295 = 8*dlmu2*dppmu3
      t3302 = t88+t193
      t3313 = 24*q1q2
      t3321 = t358+t2006-t359
      s1 = (t106*dpmu2+(t309+t599)*dlmu2+t124*dqmu2+4*t3209*dmu2mu3)*gp_
     #mu1+(-t2402+t2345*dmu1mu2+t1827*dpmu2+t616+t3217+t1108)*gp_mu3_l_q
     #+(-t622+t2345*dmu1mu3-t628)*gp_l_mu2_q+(t3223*dqmu1+(-t147+t2414-t
     #635)*dpmu1+t3227*dlmu3+t3229*dppmu3-4*t3209*dmu1mu3)*gp_mu2+(t2336
     #*dqmu1+(-t191-t86-t194-t348)*dpmu1+t3239*dlmu1)*gp_mu3_l_mu2+((t32
     #43+t661)*dppmu2+4*t3209*dmu1mu2+(t661+t3250)*dpmu2+(16*t3253*dpmu1
     #+t2314)*dlmu2+t3250*dqmu2)*gp_mu3+(t263*dmu1mu2+(-16*t425-16*t254)
     #*dpmu2-t452-16*t778*dlmu1)*gp_q+(t309*dmu1mu2+(t3272+t1935)*dpmu2-
     #t3276+(t1933*dqmu2+t2471)*dpmu1+t1935*dqmu2)*gp_l
      t3333 = s1-2*q2s*(t2006+3*q1s+2*ppq1)*gp_mu3_mu2_mu1+(-t539-t540-t
     #86-t3289)*gp_mu3_l_mu2_q_mu1+(8*t3270*dpmu2-t3295+8*t1932*dqmu2+t2
     #345*dmu2mu3)*gp_l_q_mu1+(t2336*dqmu3+t3302*dlmu3)*gp_l_mu2_mu1+(4*
     #t1986+8*t3253*dlmu2+4*t1988)*gp_mu3_q_mu1+(-4*t146+(t360+t3313+t15
     #6)*dlmu3-4*t634)*gp_mu2_q_mu1+((-t358-t94-t359)*dpmu2+t3321*dlmu2+
     #(-t359-t94-t358+t194+t348+t2014)*dqmu2)*gp_mu3_l_mu1+(4*t660+4*t64
     #0+8*ppq2*dlmu1)*gp_mu3_mu2_q
      t3337 = -8*t2024
      t3343 = 16*t1046*dlmu3
      t3349 = 8*t1007
      t3355 = t949*dlmu2
      t3379 = t1227*dpmu2
      t3380 = t1227*dqmu2
      t3382 = (t3379+t3380)*gp_mu1
      t3385 = t1229*dmu1mu2
      t3388 = 4*t3057
      t3391 = -t3388*gp_l_mu2_mu1
      s2 = -8*t45*gp_mu3*Dij_2345(2,1)+t2966*Eij_12345(2,6)+(t2970*gp_q+
     #t2973*gp_l)*Eij_12345(4,23)+((8*t1681+8*t7*q1s)*gp_mu1+(8*t2982*dq
     #mu1+t2985)*gp_mu3_l_q-8*t563*gp_l_mu2_q+(t2417*dqmu1+16*t1124*dqmu
     #3)*gp_mu2+((t111*dlmu2-t672)*dqmu1+t1139*dqmu2)*gp_mu3+(t2639+(t12
     #39-t438)*dqmu1+(-16*t250-16*t417)*dlmu2-t452+t3009)*gp_q+(-t58+t30
     #12*dqmu2*dqmu1-t70)*gp_l+t2760)*Eij_12345(3,9)
      s3 = s2+((t1218-t1651+t3020)*gp_mu1+(t3023*dqmu1+(8*dpmu3-8*dppmu3
     #)*dpmu1+t2305*dmu1mu3+t143)*gp_mu2+(t3032*dppmu2+t3034-t989*dpmu2+
     #t525)*gp_mu3+(t1024*dmu1mu2-t533-8*t34)*gp_q+(-8*dqmu2-12*dpmu2)*g
     #p_mu3_q_mu1+(-t520-t521)*gp_mu3_mu2_q+(t2122+t538)*gp_mu3_mu2_mu1+
     #(t544+12*dpmu3)*gp_mu2_q_mu1)*Dij_1235(1,3)+((t3058*dppmu2+t3060*d
     #pmu2+t3062+t3063)*gp_mu1+t3068*gp_mu2+(t3070*dppmu2+t3072*dpmu2+16
     #*t610*dqmu1+16*t613*dlmu2-t526+t315*dqmu2)*gp_mu3+t3084*gp_q+t3084
     #*gp_l)*Eij_12345(4,39)
      s1 = s3+(((t889+t3090*dpmu1+t710)*dppmu2+(t3094*dqmu1+(-t713+t470-
     #t62-t284)*dpmu1+t2233+t820+t318-t275+t255)*dpmu2+(-t3100-t2548)*dq
     #mu1+(t3103*dpmu1+t316)*dlmu2+t3094*dqmu2*dpmu1+t2233*dqmu2)*gp_q+(
     #t2838+(t3103*dqmu1-t1196+t317)*dpmu2+t3117-t1743+t3103*dqmu2*dpmu1
     #+t3120)*gp_l)*Eij_12345(4,32)+((t3126*dqmu1+t3126*dpmu1+t2595+t312
     #9)*gp_mu2+(t3032*dpmu2+t525+t616-t1108)*gp_mu3+(8*t20+8*t408-8*t31
     #35)*gp_q-4*t3031*gp_mu3_mu2_q)*Dij_1235(2,3)+(t1024*dppmu2+t2780*d
     #pmu2+t3126*dlmu2+t2780*dqmu2)*gp_mu1*Dij_1345(2,5)
      s2 = s1+t2970*gp_l*Eij_12345(4,18)+(8*t1035-8*t1033)*Dij_1235(2,1)
     #+(t3160+t3163+t3166+t3167)*D0_1234+(8*t37*t2611+8*t45*gp_mu3_l_q+8
     #*t54*gp_l_mu2_q+16*t1711*dqmu3*gp_mu2+t2674*dqmu1*gp_mu3+(-t58+(-t
     #60+t1745)*dqmu1-t70)*gp_q)*Eij_12345(3,2)
      s3 = s2+(((t338+t742)*dqmu1+t1695-t3129)*gp_mu2+(-t2700-t1105-t616
     #+t1108)*gp_mu3+(4*dpmu1-t518+t1829)*gp_mu3_mu2_q+(-t534-t1269+8*t3
     #135)*gp_q)*Dij_1235(2,6)+t3333*Eij_12345(1,3)
      t3395 = s3+t3337*gp_mu3*Dij_1234(2,6)+((t1087+t1044*dlmu2+t3343)*g
     #p_mu1+(t1101*dpmu2-t616+t2116)*gp_mu3_l_q+t3349*gp_l_mu2_q+(t941*d
     #qmu1+t2418+t1126+t2889)*gp_mu2+(t1140*dpmu2+t3355*dqmu1+t1064*dlmu
     #2-t682+t950*dqmu2)*gp_mu3+((-t689-t255-t253-t1166)*dpmu2+((-t63-t2
     #84-t282-t338)*dpmu1-t249)*dlmu2-t452)*gp_q+((t2683*dqmu1-t1196)*dp
     #mu2-t2843+t2685)*gp_l-t2849*gp_l_q_mu1)*Eij_12345(3,17)+(t3382+t12
     #29*dpmu1*gp_mu2+t3385*gp_q+t3385*gp_l+t3388*gp_mu2_q_mu1+t3391)*Di
     #j_1245(2,3)
      t3408 = -t1787
      t3413 = t61+t64+t568-t63
      t3417 = -t3413
      t3428 = t3417*dqmu2
      t3445 = -t1442+t1603*dlmu2
      t3467 = t2805*dqmu1
      t3479 = (t338+t64+t713)*dpmu1
      t3491 = 16*t973
      t3494 = 8*t1640
      t3502 = 48*dqmu2
      t3516 = 8*dqmu3*gp_mu2_q_mu1
      t3517 = (-t3491-t60+t1745+t1403)*gp_mu1+t1243+t3494+(t2466+16*q1q2
     #*dmu1mu3+t251-t574)*gp_mu2+16*t26+(t111*dmu1mu2-t1305+(-32*dlmu2-t
     #3502)*dqmu1-t1259-t526)*gp_mu3+(t3271*dmu1mu2-t3082-t1271-t1270)*g
     #p_q+8*t1646-8*t2525+16*t18+16*t21+8*t1644+t3516
      t3524 = t306*dmu2mu3
      t3527 = t264*dmu1mu3
      t3556 = -t2130
      t3562 = (t268*dppmu2+t707*dpmu2+t707*dlmu2+t707*dqmu2+t3524)*gp_mu
     #1-t3494+(t708+t709-t288-t290-t292+t448+t3527-t275)*gp_mu2+t3072*gp
     #_mu3_l_mu2+((-t625-t1255-t2573)*dppmu2+t256*dmu1mu2+(t1780+48*dqmu
     #1+t427-t137)*dpmu2+(t3502+64*dlmu2)*dqmu1+(48*dpmu1+t444-t445)*dlm
     #u2+32*t47+t432*dqmu2)*gp_mu3+(t1584+16*dppmu3*dmu1mu2+t1270+t3083)
     #*gp_q+(t1584+t301*dmu1mu2+t1271+t1273+t3083)*gp_l-t1282+t1227*gp_l
     #_mu2_mu1+t1227*gp_mu2_q_mu1+16*t3556*gp_mu3_q_mu1+t606*gp_mu3_mu2_
     #q+8*t3556*gp_mu3_l_mu1
      t3568 = (t1024*dpmu2+t1024*dqmu2)*gp_mu1
      t3575 = -t99
      t3576 = t3575*dpmu3
      t3580 = t3575*dppmu3
      t3583 = t100+t1043+t1298-t1300
      t3588 = -t339+t544+t543+t61
      t3594 = t161+t162+t160
      t3596 = -t157
      t3597 = t3596*dppmu3
      t3603 = -t2796
      t3607 = 16*t2794*dppmu3
      t3610 = t87+t156+t89+t88+t101
      t3614 = -t3610
      t3625 = t2805*dppmu1
      t3634 = t282+t64+t338
      t3638 = -t3634
      t3641 = (-t2824-t696-t62+t470)*dpmu1
      t3644 = -t726-t283+t63-t284
      t3648 = t3644*dqmu2*dpmu1
      t3660 = 8*t886
      t3661 = t3588*dpmu2
      t3662 = 8*t799
      t3668 = -t1168
      t3669 = t3668*dlmu3
      t3678 = t63-t61-t338-t469
      t3709 = t617*dlmu3
      t3712 = t264*dpmu3
      t3713 = t264*dqmu3
      t3714 = t260*dppmu3
      t3715 = t3712+t3713+t3223+t3714
      t3719 = t3229+t260*dppmu1
      t3723 = t3243*dppmu3
      t3724 = t265*dppmu1
      t3732 = t256*dpmu1
      t3733 = t256*dlmu1
      t3734 = t273*dppmu1
      t3738 = t256*dqmu2
      t3748 = -t301
      t3749 = t3748*dqmu1
      t3754 = t469+t713-t470+t568
      t3755 = t3754*dpmu1
      t3758 = t1608*dlmu2
      t3763 = t3754*dqmu2
      t3764 = t3763*dpmu1
      t3802 = 4*t498
      t3809 = 8*t1252
      t3813 = -t2402+t3809+8*t23-t1106*dlmu2+t2985
      s2 = (-8*t1031+8*t1033)*Dij_1245(2,1)+t415*Dij_1234(2,4)+(t1603*gp
     #_mu2-16*t1546*gp_mu1+16*t1513*gp_mu3+t3408*gp_l)*Eij_12345(4,41)+(
     #((t3413*dqmu1+t574+t73)*dppmu2+(t3417*dqmu1-t574+t317)*dpmu2+(t341
     #7*dlmu2+(-t469-t1202-t713+t831)*dqmu2)*dqmu1+(-t574+t317)*dlmu2+t3
     #428*dpmu1+(t316-t840*dqmu3+t318-t275+t255)*dqmu2)*gp_q+(t56-t58+(-
     #t60+t3428)*dqmu1-t70+t317*dqmu2)*gp_l)*Eij_12345(4,21)
      s1 = s2+(t3445*gp_q+t2653*gp_l)*Eij_12345(4,29)+((16*t931*dppmu2+(
     #-t2861*dpmu3-t1043-t1298+t1300)*dpmu2-t2037-t1048)*gp_mu1-t1052+(-
     #t1055+(-t2883*dpmu3-t2034-t940+t1315)*dpmu1+t2796*dpmu3)*gp_mu2+(t
     #953*dppmu2+(t3467+(t103+t101+t89+t161+t88)*dpmu1+t2815+t2816)*dpmu
     #2+t2805*dpmu1*dlmu2+t2805*dqmu2*dpmu1)*gp_mu3+((t3479+t316)*dpmu2+
     #t897+t903)*gp_q+(t888+(t784+t3479)*dpmu2+t903)*gp_l-t968)*Eij_1234
     #5(3,7)+t3517*Eij_12345(3,22)+t3562*Eij_12345(3,24)+(t3568+t3166+4*
     #t1023*gp_l_mu2_mu1+t3163+t3165*gp_mu3_l_mu1)*Dij_1345(1,1)
      s4 = s1
      s7 = ((t3576-t1043-t1298+t1300)*dppmu2+(t122*dpmu3+t2394+t106+t358
     #0)*dpmu2+t3583*dlmu2+t3583*dqmu2)*gp_mu1+(t977+t3588*dpmu1+t2169)*
     #gp_l_mu2_q+((t158+t2034+t940-t1315)*dqmu1+(t3594*dpmu3+t159+t2414+
     #t3597)*dpmu1+(t2419+t3596*dppmu1)*dpmu3+t3603*dqmu3+t3603*dlmu3-t2
     #893+t3607)*gp_mu2
      s6 = s7+((t3467+t3610*dpmu1+t2815+t2816)*dppmu2+(t3614*dqmu1+(-t10
     #3-t162-t102-t160-t161)*dpmu1+t3614*dlmu1+t3610*dppmu1)*dpmu2+(t335
     #5+t1067)*dqmu1+(t3614*dpmu1+t950+t3625)*dlmu2+t3614*dqmu2*dpmu1+(t
     #950+t3625)*dqmu2)*gp_mu3+((t3634*dpmu1+t316)*dppmu2+(t3638*dqmu1+t
     #3641+t275-t255)*dpmu2+t2682+t3644*dpmu1*dlmu2+t3648)*gp_q+((t303+t
     #316)*dppmu2+(t3644*dqmu1+t3641+t2233+t317+t318)*dpmu2+t3117+t3638*
     #dpmu1*dlmu2+t3648+t3120)*gp_l+(-t3660+t3661+t1651+t3662)*gp_l_q_mu
     #1
      s7 = Eij_12345(3,14)
      s5 = s6*s7
      s3 = s4+s5
      s4 = s3+((t1440-t1442-t2286+(-t426-t249-t251+t3669+t255)*dlmu2-t22
     #91)*gp_q+((t424+t426+t249+t251+t1169-t255)*dppmu2+(-t424-t426-t249
     #-t251+t3669+t255)*dpmu2+(t3678*dlmu2-t438)*dqmu1+(t3678*dpmu1+t366
     #8*dpmu3+t3668*dqmu3+(-96*dlmu1+t479)*dlmu3-t275+t483)*dlmu2-t452+(
     #-t249-t251+t3669+t255)*dqmu2)*gp_l)*Eij_12345(4,22)
      s5 = s4
      s7 = (t923*gp_l-t2285)*Eij_12345(4,12)
      s10 = ((t595+t596+t593+t599)*dppmu2+t312*dpmu2+t312*dlmu2+t312*dqm
     #u2)*gp_mu1+(t1227*dqmu1+t2104+t2782+t2709+t3709-t143+t628)*gp_l_mu
     #2_q+(t3715*dqmu1+t3715*dpmu1+t3719*dpmu3+t3719*dqmu3+t3719*dlmu3+t
     #3723+t3724)*gp_mu2
      s9 = s10+((t273*dqmu1+t273*dpmu1+t332+t333)*dppmu2+(t256*dqmu1+t37
     #32+t3733+t3734)*dpmu2+(t256*dlmu2+t3738)*dqmu1+(t3732+t3733+t3734)
     #*dlmu2+t3738*dpmu1+(t3733+t3734)*dqmu2)*gp_mu3+((t3749-t1608*dpmu1
     #+t275-t255)*dppmu2+(t1617+t3755+t710+t73+t795)*dpmu2+(t3758+t1622)
     #*dqmu1+(t3755+t710+t73+t795)*dlmu2+t3764+(t710+t73+t795)*dqmu2)*gp
     #_q+((t3749+t3748*dpmu1+t275-t255)*dppmu2+(t3754*dqmu1+t3755+t794+t
     #827+t842+t275-t255)*dpmu2+(t3758+t3763)*dqmu1+(t1609+t710+t73+t795
     #)*dlmu2+t3764+(t794+t827+t842+t275-t255)*dqmu2)*gp_l+(t1229*dppmu2
     #+t3379+t1227*dlmu2+t3380)*gp_l_q_mu1
      s10 = Eij_12345(3,4)
      s8 = s9*s10
      s6 = s7+s8
      s2 = s5+s6
      t3820 = s2+8*t1034*dpmu2*Dij_1234(2,1)+8*t1800*gp_mu3*Dij_1234(2,3
     #)+(t508*dqmu2*gp_mu1+t531*gp_q+t531*gp_l+t3802*gp_mu2_q_mu1-t3802*
     #gp_l_mu2_mu1)*Dij_1245(2,2)+t3813*gp_mu3*Dij_2345(2,6)+8*t1593*Dij
     #_1345(2,7)
      t3823 = t1735*dqmu2
      t3851 = 32*t20
      t3872 = (t1735*dpmu2+t60+t3823+(-t89-t160+t87+t88-t86)*dmu2mu3)*gp
     #_mu1+32*t1242+32*t1640+(t424+t3272-t249-t251-t255+(t160+t156+t89-t
     #88-t94)*dmu1mu3)*gp_mu2+(-40*dqmu1+t1255)*gp_mu3_l_mu2+(t1577+(t20
     #6-t160-t156-t89+t93)*dmu1mu2+(-t2573-t1255+t625)*dpmu2+32*t2982*dq
     #mu1+t1259-t2576-t526)*gp_mu3+(t1369*dmu1mu2+t1271-t1586-t1270+t127
     #3+32*t1264)*gp_q+(16*t1023*dmu1mu2-t1271-t3851-t1273+t1270)*gp_l+3
     #2*t1281+(-t742-t338-t339)*gp_l_mu2_mu1+(t2702-40*dlmu2+16*dpmu2)*g
     #p_mu3_q_mu1+(-8*dppmu2-t556-t2702)*gp_mu3_l_mu1+(-t339-t64-t742)*g
     #p_mu2_q_mu1+(-t359-t94+t160+t194+t191-t358+t86)*gp_mu3_mu2_mu1+36*
     #gp_mu3_l_mu2_q_mu1+(t521-40*dlmu1+t624)*gp_mu3_mu2_q
      t3875 = 8*t1571
      t3888 = pq2*q1s
      t3889 = 8*t3888
      t3890 = pq1*q1q2
      t3891 = 16*t3890
      t3892 = pq1*q2s
      t3894 = t3889-t3891-8*t3892
      t3914 = t260*dqmu1
      t3915 = 32*ppp
      t3918 = t306*dlmu1
      t3946 = t193+t192
      t3951 = pq2*dppmu2
      t3957 = t360+t3313+t357
      t3962 = pq1*dppmu2
      s1 = ((t595+t262+t309+t3580)*dpmu2+t307*dlmu2+((t94+t88-t89-t101-t
     #86)*dpmu3-t1300)*dqmu2+t3894*dmu2mu3)*gp_mu1+(t3712*dqmu1+(t261+t3
     #713+t593+t3597)*dpmu1+t3227*dpmu3-t3894*dmu1mu3)*gp_mu2+((-t193-t1
     #92+t3238)*dpmu1+8*t1361)*gp_mu3_l_mu2+(t157*dpmu1*dppmu2+t3894*dmu
     #1mu2+(t3914+(t3915-t86-t94-t88-t89-t101)*dpmu1+t3918+t99*dppmu1)*d
     #pmu2+t310*dpmu1*dlmu2+(-t94-t88+t89+t101+t86+t206)*dqmu2*dpmu1+16*
     #t1361*dqmu2)*gp_mu3+(t261*dmu1mu2+t3527*dpmu2+t1949*dpmu1)*gp_q+(t
     #307*dmu1mu2+t1943*dpmu2+t3524*dpmu1)*gp_l
      t3972 = s1+(-12*t3888-t2333+12*t3892+24*t3890+t2331)*gp_mu3_mu2_mu
     #1+(t3946*dpmu3-8*t1299)*gp_l_mu2_mu1+(8*t3951+(t3238-t357-t360-t33
     #13)*dpmu2)*gp_mu3_q_mu1+(t3957*dpmu3-8*t1314)*gp_mu2_q_mu1+(-8*t39
     #62+t3946*dpmu2)*gp_mu3_l_mu1+(-t3957*dpmu1+8*t2794)*gp_mu3_mu2_q
      t3980 = dpmu3*gp_l
      t3998 = -t2776
      t4004 = -t162-t156-t101
      t4012 = 16*t2794
      t4018 = t156+t87+t101+t162+t88
      t4043 = t469+t726+t283
      t4068 = q1s*q2s
      t4070 = t3891-t3889-4*t4068
      t4080 = 8*t243
      t4107 = t101+t86+t89+t93+t206
      t4123 = 16*t778*dqmu1
      t4140 = 4*ppq1
      t4146 = t87+t192
      t4154 = t3238-t359-t94-t358-t89-t101-t357
      t4158 = t89+t101+t357
      t4163 = q1s*dppmu2
      t4169 = t191-t348-t2014-t194
      s1 = ((-t1041+t2394-t1093)*dpmu2+t308*dlmu2+(-t1041+t2394+t309+t21
     #51)*dqmu2+t4070*dmu2mu3)*gp_mu1+((-8*dpmu1-8*dqmu1)*dppmu2+t3034-8
     #*t936-t616)*gp_mu3_l_q+(-t4080+t3023*dpmu1+t2348*dmu1mu3)*gp_l_mu2
     #_q+((t3713+t265)*dqmu1+(-t1422+t159-t1321)*dpmu1+t3227*dqmu3-t4070
     #*dmu1mu3)*gp_mu2+(8*t3248*dqmu1+4*t1334+4*t1346)*gp_mu3_l_mu2+((t3
     #914+16*t1338)*dppmu2+t4070*dmu1mu2+(t3249*dqmu1+t1151+t1347)*dpmu2
     #+t4107*dqmu2*dqmu1+(t88+t206-t101-t156-t89)*dqmu2*dpmu1+(t3918+t21
     #41*dppmu1)*dqmu2)*gp_mu3+(t262*dmu1mu2+(-t244+t1735*dpmu1+t1935)*d
     #pmu2-t4123+(t3823+t2471)*dpmu1+t274*dqmu2)*gp_q+(t308*dmu1mu2-16*t
     #973*dpmu1+t903+t1943*dqmu2)*gp_l
      t4175 = s1+(-32*t3890-t2327+6*t4068+16*t3888+t2325)*gp_mu3_mu2_mu1
     #+(-t94-t538-t4140)*gp_mu3_l_mu2_q_mu1+(-t1218+t1651+t3020)*gp_l_q_
     #mu1+(-4*t1040+t4146*dqmu3-4*t1092)*gp_l_mu2_mu1+(8*q1q2*dppmu2+t29
     #59+t2305*dlmu2+t4154*dqmu2)*gp_mu3_q_mu1+(-t1979+t4158*dqmu3+t2305
     #*dlmu3-t1983)*gp_mu2_q_mu1+(-4*t4163-4*t1426+t4146*dqmu2)*gp_mu3_l
     #_mu1+(t4169*dqmu1+t1869+t2305*dlmu1+8*t1131)*gp_mu3_mu2_q
      t4181 = 4*pq2
      t4190 = -t101-t86-t89-t88-t87
      t4193 = (-t93+t87+t102)*dppmu3
      t4198 = -t94+t89+t101+t86
      t4212 = 8*t245
      t4220 = -t156-t89-t101-t88-t94
      t4245 = -t4107
      t4247 = -t3915+t94+t86
      t4250 = (-t93-t206+t88+t94)*dlmu1
      t4252 = (-t102-t87+t93+t206)*dppmu1
      t4298 = -t3302
      t4308 = -t4154
      t4315 = -t4158
      t4320 = -t4169
      t4325 = -t3321
      s2 = (t599*dppmu2+(-16*t931+t4190*dqmu3-t1298+t4193)*dpmu2+(t595+t
     #596+t593+t311)*dlmu2+(t4198*dpmu3+t3575*dqmu3+t4198*dlmu3+t4193)*d
     #qmu2)*gp_mu1+((8*dqmu1+8*dpmu1+8*dlmu1)*dppmu2+t2112-8*t768-t3217-
     #t1108)*gp_mu3_l_q+(t4080+t4212*dpmu1+t628+t2100*dmu1mu3)*gp_l_mu2_
     #q+((t261+t262+t263+t3714)*dqmu1+(-16*t1053-t2034+t4220*dlmu3+(t162
     #+t156-t165+t160)*dppmu3)*dpmu1+t3918*dpmu3+t3918*dqmu3+t3918*dlmu3
     #+t3723+t3724)*gp_mu2
      s1 = s2+((-t348-t2014-t194-t359-t3238)*dqmu1+(t191+t86+t194+t348+t
     #193+t88-t3238)*dpmu1-t3239*dlmu1+t3239*dppmu1)*gp_mu3_l_mu2+((t264
     #*dqmu1+(-t162-t160-t156+t165+t206)*dpmu1+t3229+t333)*dppmu2+(t4245
     #*dqmu1+t4247*dpmu1+t4250+t4252)*dpmu2+t4245*dqmu2*dqmu1+((-t88-t94
     #-t165-t206)*dpmu1+t306*dppmu1)*dlmu2+t4247*dqmu2*dpmu1+(t4250+t425
     #2)*dqmu2)*gp_mu3+((t3712+t3713+t3223)*dmu1mu2+(t244+t247+t255+t257
     #)*dpmu2+t4123+t246*dqmu2*dpmu1+(t255+t257)*dqmu2)*gp_q+((t595+t596
     #+t593)*dmu1mu2+t247*dpmu2+t3276+(16*t778+t256*dmu2mu3)*dpmu1)*gp_l
      t4339 = s1+(4*t3892+t2329+8*t3890-t2331-4*t3888-t2325+t2327+t2333)
     #*gp_mu3_mu2_mu1+(t539+t540+t86+t3289+t94+t538+t4140)*gp_mu3_l_mu2_
     #q_mu1+(t4212*dpmu2+t3295+8*t778+t2101)*gp_l_q_mu1+(t4298*dpmu3+(-t
     #348-t2014-t194-t87-t192)*dqmu3+t4298*dlmu3+t3321*dppmu3)*gp_l_mu2_
     #mu1+((-t2014-t348-t194+t191+t3238)*dppmu2+t4308*dpmu2+(-t358-t2006
     #-t191-t3238)*dlmu2+t4308*dqmu2)*gp_mu3_q_mu1+(t4315*dpmu3+t4315*dq
     #mu3+(-t360-t3313-t156-t358-t2006)*dlmu3+t4320*dppmu3)*gp_mu2_q_mu1
     #+(t3321*dppmu2+t4325*dpmu2+t4325*dlmu2+(t359-t94-t88-t194-t348-t20
     #14)*dqmu2)*gp_mu3_l_mu1+(t4320*dqmu1+t4320*dpmu1+(-t358-t2006-t191
     #)*dlmu1+t4169*dppmu1)*gp_mu3_mu2_q
      t4354 = -8*dpmu3+8*dlmu3+8*dqmu3+8*dppmu3
      t4356 = 8*t948
      t4369 = -t2114*dpmu2
      t4373 = -t998
      t4379 = (t4354*dpmu2+t3662-t1401+t4356*dmu2mu3)*gp_mu1+(-t977-t435
     #4*dpmu1-t623-t4356*dmu1mu3)*gp_mu2+12*t2570+(t985+t4356*dmu1mu2+t4
     #369-t2116+t526)*gp_mu3-12*t1661-t993+t4373*gp_q+t1001-12*t2582+12*
     #t2807*gp_mu3_mu2_mu1-t995+t4373*gp_l
      t4405 = t64-t63+t568+t282
      t4428 = 8*t254*q2s
      t4432 = 8*t1519*dppmu2
      t4447 = t282+t696+t283
      s2 = t3872*Eij_12345(2,11)+(t3875-t1022+t743*dlmu2-t1026)*gp_mu1*D
     #ij_1345(2,6)+t3972*Eij_12345(1,1)-16*t1800*dlmu3*gp_l*Eij_12345(4,
     #3)
      s3 = s2-16*t936*t3980*Eij_12345(4,6)-16*t45*dqmu3*gp_q*Eij_12345(4
     #,2)
      s5 = s3
      s8 = (16*t929*dppmu2+t3575*dlmu3*dpmu2+(t3576-t1043-t2867+t1300)*d
     #lmu2-t3343)*gp_mu1+(-t985+t3998*dpmu2+t616+t2116)*gp_mu3_l_q-t3349
     #*gp_l_mu2_q+(t2882+(-t1839-t149+t4004*dlmu3+t635)*dpmu1+(t4004*dlm
     #u1+t641)*dpmu3-t2889+(-32*t943+t4012)*dlmu3+t2893)*gp_mu2
      s7 = s8+(t951*dppmu2+(t659+t1876+t4018*dlmu1-t641)*dpmu2+t2806*dqm
     #u1+(t4018*dpmu1+t2808*dlmu1+t2816)*dlmu2+t682+t2815*dqmu2)*gp_mu3+
     #(16*t248*dppmu2+(t689+t318+t255)*dpmu2+((t282+t284+t338)*dpmu1+t31
     #6)*dlmu2+t452)*gp_q+((-t784+t1381*dpmu1+t249)*dppmu2+(t4043*dqmu1+
     #(t2824+t1201+t62)*dpmu1+t1180*dpmu3+t443+t433-t275+t255)*dpmu2+t16
     #21+t2822*dpmu1*dlmu2+t4043*dqmu2*dpmu1+t1624)*gp_l+(t3660+(t339-t6
     #4-t61-t544)*dpmu2-t1651-t1773)*gp_l_q_mu1
      s8 = Eij_12345(3,19)
      s6 = s7*s8
      s4 = s5+s6
      s1 = s4+t4175*Eij_12345(1,2)+t3445*gp_l*Eij_12345(4,11)
      s2 = s1+(t3163+t3160+(t4181+4*q1q2+2*q2s)*gp_mu3_mu2_mu1+t3166+t31
     #67)*D0_1345+t4339*Eij_12345(1,4)+t2657*gp_q*Eij_12345(4,8)+((t3661
     #+t1651)*gp_mu1-t3588*dpmu1*gp_mu2-t1656-t1660-t1658+t1662)*Dij_124
     #5(2,5)
      s3 = s2+t4379*Dij_1245(1,1)+(4*t1402*gp_mu1+(-t1309*dqmu1-t1695-4*
     #t1394)*gp_mu2+(-t2699+4*q1s*dmu1mu2-t525+t616)*gp_mu3-4*t988*gp_mu
     #3_mu2_q+8*t18-6*q1s*gp_mu3_mu2_mu1+(8*t408+8*t34-8*t874)*gp_q-t351
     #6)*Dij_1235(1,2)
      s4 = s3+((t1616+(t4405*dqmu1+t1196+t73)*dpmu2+(t3100+t1773)*dqmu1+
     #t1743+t4405*dqmu2*dpmu1+t808)*gp_q-t2657*gp_l)*Eij_12345(4,27)
      t4471 = s4+(t309*dlmu2*gp_mu1+t3813*gp_mu3_l_q-t2618+(-t2620-t1518
     #-t1521-t2624+((-t86-t645-t101)*dlmu1+t641)*dlmu3+t4428)*gp_mu2+(-t
     #4432+t1527+t2629+(t661+(t101+t86+t94+t645+t88)*dlmu1-t641)*dlmu2+t
     #2632)*gp_mu3+(t1440+t687*dlmu2)*gp_q+((-16*t425-16*t28-16*t248-16*
     #t250-16*t254)*dppmu2+(t1165+t689+t1166+t1167+t2288)*dpmu2+(t4447*d
     #lmu2+t1174)*dqmu1+(t4447*dpmu1+t1179+t1453+t1169-t275+t448)*dlmu2+
     #t701+(t1166+t1167+t2288)*dqmu2)*gp_l+(t3875-t1022+(-t544-t747+t339
     #-t742)*dlmu2-t1026)*gp_l_q_mu1)*Eij_12345(3,13)+(16*t936*t3980+16*
     #t936*t2357)*Eij_12345(4,7)
      t4472 = 16*t1216
      t4555 = -dpmu3+dlmu3-dppmu3
      t4570 = (t1808*dlmu2-8*t7*q2s)*gp_mu1+(-t3809+t1827*dlmu2-t2985)*g
     #p_mu3_l_q+t2618+(-t2620+t1518-t1521-t2624+t1857*dlmu1*dlmu3-t4428)
     #*gp_mu2+t1872*dlmu1*gp_mu3_l_mu2+(t4432+t1527+t2629+(-t661+t233)*d
     #lmu2+t2632)*gp_mu3+(8*t1516*dmu1mu2+t1442+(-t426-t1245-t255)*dlmu2
     #+8*t1914*dlmu1+t2291)*gp_q+((-16*t28-16*t248-16*t254-16*t250-16*t2
     #52-16*t425)*dpmu2+(-16*t1795-16*t7)*dqmu1+(t3272-t255+t1935)*dlmu2
     #-t452+t3009+t3229*dmu2mu3)*gp_l+(t1022+8*t4555*dlmu2+t1026)*gp_l_q
     #_mu1+8*t1120*gp_l_mu2_mu1-12*t2627*gp_mu3_q_mu1+12*t1516*gp_mu2_q_
     #mu1+t2002*dlmu2*gp_mu3_l_mu1-12*t1519*gp_mu3_mu2_q
      t4580 = -t2131
      t4604 = 64*pq1
      t4621 = t87-t86
      t4631 = -t1971
      t4639 = -16*q1q2+16*pq2
      t4680 = -t93+t160+t89+t86+t103-t206+t88
      t4712 = (t713+t283+t63+t284)*dpmu1
      t4747 = -t2469
      s2 = ((t2866+t2034-t1298-t1300)*dppmu2+((-t102-t4604)*dpmu3+(-t101
     #-t103-t88-t89-t156)*dqmu3-t1298+(t88+t103+t101-t93)*dppmu3)*dpmu2+
     #((16*q1q2+16*pq1)*dpmu3+t2039*dqmu3+t2867)*dlmu2+((-t102-t103+t156
     #+t89+t101)*dpmu3-t1043+t940+t4621*dppmu3)*dqmu2)*gp_mu1+(-t1778+16
     #*t2113*dpmu2+t1783+t526)*gp_mu3_l_q+(-t784+t4631*dpmu1+(-t625+24*d
     #ppmu1)*dpmu3)*gp_l_mu2_q
      s3 = s2+((t4639*dpmu3+t1313+t2063*dlmu3)*dqmu1+((-t160-t1854-t162)
     #*dpmu3-t2034+(-t88-t101-t161-t89-t87)*dlmu3+(t89+t161-t165)*dppmu3
     #)*dpmu1+((-t86-t88-t160-t89-t87)*dlmu1+t2173*dppmu1)*dpmu3+(t2914-
     #t4012)*dqmu3+(-32*t1327+t2058)*dlmu3+t2805*dppmu3*dlmu1-t3607)*gp_
     #mu2
      s1 = s3+(t1333+(t194+t191+t101+t192+t86+t103)*dpmu1+t2897-24*t1361
     #)*gp_mu3_l_mu2+(((t165-t94+t206-t161-t89-t88)*dpmu1+t950+t3625)*dp
     #pmu2+((-16*pq2-16*ppp+16*q1q2-16*ppq1+16*pq1)*dqmu1+(-t3915+t1854+
     #t160+t102+t4604+t162)*dpmu1+t4680*dlmu1-t4680*dppmu1)*dpmu2+(t3355
     #-t1352*dqmu2)*dqmu1+(-16*pq1-16*ppp+16*pq2-16*q1q2-16*ppq2)*dpmu1*
     #dlmu2+(t103+t102+t86-t206)*dqmu2*dpmu1+(t4621*dlmu1-t4621*dppmu1)*
     #dqmu2)*gp_mu3+((16*t425+16*t248-16*t2069)*dppmu2+(t3596*dpmu3-t203
     #4-t940+t1315)*dmu1mu2+(t1912*dqmu1+t4712+t316+t795-t4220*dmu1mu3)*
     #dpmu2-t1379*dqmu1+(-16*t4555*dpmu1+t2567-t288)*dlmu2+((t282+t338)*
     #dqmu2-16*t1377)*dpmu1+(t316+t2567)*dqmu2+t2796*dmu2mu3)*gp_q
      t4792 = s1+((-t784+(-16*dqmu3+16*dpmu3)*dpmu1+t2937)*dppmu2+t3583*
     #dmu1mu2+(-16*t1692*dqmu1+t4712+t317-t2937)*dpmu2+(t1773+t2564+t311
     #5)*dqmu1+(t4747*dpmu1-t2937)*dlmu2+((t338+t470+t61)*dqmu2-t4190*dm
     #u2mu3)*dpmu1+(t710-t2937)*dqmu2-t1363*dmu2mu3)*gp_l+(24*t886+t4631
     #*dpmu2-t1401-t1773)*gp_l_q_mu1+((-t103-t192-t348)*dpmu3-t1413-t295
     #0+24*t1299)*gp_l_mu2_mu1+(-24*t3951+(t101+t94+t358+t360+t359+t161)
     #*dpmu2+t2953+t1419)*gp_mu3_q_mu1+((-t360-t101-t161)*dpmu3-t1423-t2
     #956+24*t1314)*gp_mu2_q_mu1+(24*t3962+(-t88-t87+t359-t348)*dpmu2-t2
     #960-t1428)*gp_mu3_l_mu1+(t1431+(t156+t89-t191+t348)*dpmu1+t2963-24
     #*t2794)*gp_mu3_mu2_q
      t4828 = -t101-t156-t89-t102-t94
      t4865 = -16*q1q2+16*q1s
      t4867 = -t4865
      t4884 = -t206+t101+t94-t93+t156+t2147+t89
      t4902 = -t1373*dmu1mu3
      t4909 = -t1383
      t4910 = t4909*dmu2mu3
      t4961 = -t1417
      t4980 = -t1432
      s2 = ((-t1041+t2394-t1290-t1093)*dppmu2+((-t101-t2147-t103)*dqmu3+
     #16*t1092)*dpmu2+(t2765+(t88-t94)*dqmu3+t1091)*dlmu2+(t85+(-80*pq1-
     #32*q1s)*dqmu3+t92+(t94-t93+t2147)*dppmu3)*dqmu2)*gp_mu1+(t3998*dpp
     #mu2+t4369-t525-t2116-t526+t618)*gp_mu3_l_q+((-t338-t543+t339)*dqmu
     #1+t1971*dpmu1+t2169-t138*dqmu3-t143+t628)*gp_l_mu2_q
      s3 = s2+((-t4639*dpmu3-t1313+t4828*dlmu3+t187)*dqmu1+(-t3594*dqmu3
     #+t115+32*t1320)*dpmu1+(t2432-t1132)*dpmu3+(t4828*dlmu1+t157*dppmu1
     #)*dqmu3+(-16*t1324+t1347)*dlmu3+t1142*dppmu3*dlmu1-t2195)*gp_mu2
      s4 = s3+((t101+t194+t193+t191+t86+t2147)*dqmu1+t1335+t2426-12*t134
     #6)*gp_mu3_l_mu2
      s5 = s4
      s7 = (((t93-t89+t206-t101+t165-t86)*dqmu1+(-t160+t88)*dpmu1+t1139+
     #t2204)*dppmu2+((16*pq2+16*q1s+16*pq1-16*q1q2-16*ppp-16*ppq1)*dqmu1
     #+t4865*dlmu1+t4867*dppmu1)*dpmu2+((-16*pq1-16*q1s-16*ppp-16*q1q2-1
     #6*ppq2-16*ppq1)*dlmu2+(-t3915+t156+t94-32*ppq1+t102)*dqmu2)*dqmu1+
     #t4867*dpmu1*dlmu2+(t162+t94+t161-t206+t160)*dqmu2*dpmu1+(t4884*dlm
     #u1-t4884*dppmu1)*dqmu2)*gp_mu3
      s8 = ((-t244+t325+t251+t1374)*dppmu2+(-t1422+t2179*dqmu3-t1121+t13
     #21)*dmu1mu2+((t282+t64+t283+t470)*dqmu1+(t62+t63+t713+t284)*dpmu1+
     #t710+t255+t4902-t275)*dpmu2+(-t1744*dlmu2+(t283+t282+t831+t284)*dq
     #mu2+t4910)*dqmu1+(t2225-t290+t4902)*dlmu2+((t470+t726+t2228+t696)*
     #dqmu2+t4910)*dpmu1+(t710+t795+t255+(t87+t160+t89+t102+t86)*dmu1mu3
     #-t275)*dqmu2+(t4909*dlmu1+t1383*dppmu1)*dmu2mu3)*gp_q
      s6 = s7+s8
      s1 = s5+s6
      t4986 = s1+((-t1393+t574+t1395)*dppmu2+t2143*dmu1mu2+(-t3012*dqmu1
     #+t1196+t73-t1395)*dpmu2+(t4747*dlmu2+(t64+t283+t470)*dqmu2+(t101+t
     #94+t102+t86+t89)*dmu2mu3)*dqmu1-t2475+(t510*dqmu2+t1403)*dpmu1+(t3
     #16-t1395)*dqmu2-t1727*dmu2mu3)*gp_l+(t4472-t2271+t1309*dqmu2)*gp_l
     #_q_mu1+(-t1414+(-t348-t193-t2147)*dqmu3-t2486+12*t1092)*gp_l_mu2_m
     #u1+(t1417*dppmu2+t4961*dpmu2+t4961*dlmu2+(t87+t93+t88+t160+t360+t3
     #57)*dqmu2)*gp_mu3_q_mu1+(-t1422+(-t357-t160-t360)*dqmu3-t1121+t132
     #1)*gp_mu2_q_mu1+(12*t4163-t1427-t2494+(-20*q1s-t348+t359-t2006)*dq
     #mu2)*gp_mu3_l_mu1+((t101-t2014+t194-t165)*dqmu1+t4980*dpmu1+t4980*
     #dlmu1+t1432*dppmu1)*gp_mu3_mu2_q
      t5042 = (t3058*dqmu1+t3058*dpmu1+t710+t73+t795-t255+t275)*dppmu2+t
     #3068*dpmu2+(t3062+t3063)*dqmu1+(t3067+t316+t317+t318+t255-t275)*dl
     #mu2+t3063*dpmu1+(t316+t317+t318+t255-t275)*dqmu2
      s1 = ((-t4472+t3491+t60+t858)*gp_mu1+t851*gp_mu2+(-16*t1688+t1305+
     #(32*dqmu2+16*dlmu2)*dqmu1+t526+t72*dqmu2)*gp_mu3+(-t1584+t825*dmu1
     #mu2+t1271+t3851+t3082+t1270+32*t34+t72*dmu2mu3)*gp_q-t876*gp_l)*Ei
     #j_12345(4,44)-8*t1800*gp_mu3*Dij_2345(2,2)+t3568*Dij_1345(2,1)+t29
     #73*gp_q*Eij_12345(4,17)+8*t1795*gp_mu1*Dij_1345(2,2)+(16*t2086*dpm
     #u3*gp_mu1+16*t2053*dpmu3*gp_mu2+t953*dpmu2*gp_mu3)*Eij_12345(3,1)+
     #t4570*Eij_12345(2,3)+(t3160+4*t994+(16*t936+t616)*gp_mu3)*Dij_1234
     #(1,1)+(t3382+t4580*gp_mu3_l_mu1+(8*t766-8*t47-8*t12-8*t936)*gp_mu3
     #+t4580*gp_mu3_q_mu1+t3391)*Dij_1345(1,3)
      t5048 = s1+(t1596-t1597)*Dij_1245(2,7)-t3337*gp_mu3*Dij_2345(2,4)+
     #(-t2699+t2700-t2703*dqmu1+t616+t528)*gp_mu3*Dij_2345(2,5)+t4792*Ei
     #j_12345(2,7)+t4986*Eij_12345(2,9)+8*t45*gp_mu3*Dij_1234(2,2)-8*t15
     #95*Dij_2345(2,7)+((-8*t973-8*t37-8*t59+8*t1216)*gp_mu1+(-t621-t622
     #+t3709)*gp_mu2+12*t605*gp_mu3_l_mu2+(t2699-t2700-t1578+t615+t616+t
     #528)*gp_mu3+(8*t532-8*t1583+8*t34+8*t1263)*gp_q+(8*t20+t534+t527*d
     #mu2mu3)*gp_l+12*t2130*gp_mu3_q_mu1+(t4140+t3289+t3238-4*pq1-t4181)
     #*gp_mu3_mu2_mu1)*Dij_2345(1,3)+(t1615*gp_mu2-16*t413*gp_mu3+t3408*
     #gp_q+16*t2848*gp_mu1)*Eij_12345(4,40)+(t5042*gp_q+t5042*gp_l)*Eij_
     #12345(4,4)
      virt_pent = t885+t1552+t2282+t2854+t3395+t3820+t4471+t5048
      end
