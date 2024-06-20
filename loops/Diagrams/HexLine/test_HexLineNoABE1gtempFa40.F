c************** Calling the Fa functions*************************
       subroutine test_FaFunctionhexlineNoAbe1g40(p1mup2,p1mup3,p1mup4
     -   ,p1mup6,p2mup2,p2mup3,p2mup4,p2mup6,p3mup2,p3mup3,p3mup4,p3m
     -   up6,p4mup2,p4mup3,p4mup4,p4mup6,p5mup2,p5mup3,p5mup4,p5mup6,
     -   p6mup2,p6mup3,p6mup4,p6mup6,mup2mup3,mup2mup4,mup2mup6,mup3m
     -   up4,mup3mup6,mup4mup6,Fa40)
       IMPLICIT NONE
      Complex*16   p1mup2, p1mup3, p1mup4, p1mup6, p2mup2, p2mup3, 
     -          p2mup4, p2mup6, p3mup2, p3mup3, p3mup4, p3mup6, 
     -          p4mup2, p4mup3, p4mup4, p4mup6, p5mup2, p5mup3, 
     -          p5mup4, p5mup6, p6mup2, p6mup3, p6mup4, p6mup6
       Complex*16   mup2mup3, mup2mup4, mup2mup6, mup3mup4, mup3mup6, 
     -          mup4mup6
        common/test_FhlFunctionsNoAbe1g/F
       COMMON/Invariants/p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,s12,s23,s34,s45
     -   ,s56,s16,s123,s234,s345
       Real*8   p1sq 
       Real*8   p2sq 
       Real*8   p3sq 
       Real*8   p4sq 
       Real*8   p5sq 
       Real*8   p6sq 
       Real*8   s12, s16 
       Real*8   s23 
       Real*8   s34 
       Real*8   s45 
       Real*8   s56 
       Real*8   s123,s234,s345
      Complex*16 Fa(587),F(6395)
      Real*8 P(7129) 
        COMMON/test_PFunctionsNoAbe1g/P
       Complex*16 Fa40(547:587)
       COMMON/test_FaFunctionsNoAbe1g/Fa
       Fa40(547) = p1mup2*(p3mup3*p3mup4*F(6098)+p1mup3*p3mup4*F(6100)
     -   +p2mup3*p3mup4*F(6101)+p4mup3*p4mup4*F(6104)+p3mup4*p4mup3*F
     -   (6106)+p3mup4*p5mup3*F(6108)+p4mup4*p5mup3*F(6110)+p1mup4*p5
     -   mup3*F(6112)+p2mup4*p5mup3*F(6113)-mup3mup4*F(6116))+p2mup2*
     -   (p3mup3*p3mup4*F(6099)+p1mup3*p3mup4*F(6102)+p2mup3*p3mup4*F
     -   (6103)+p4mup3*p4mup4*F(6105)+p3mup4*p4mup3*F(6107)+p3mup4*p5
     -   mup3*F(6109)+p4mup4*p5mup3*F(6111)+p1mup4*p5mup3*F(6114)+p2m
     -   up4*p5mup3*F(6115)-mup3mup4*F(6117))
       Fa40(548) = p3mup2*p3mup4*p5mup3*F(6124)+p2mup4*p3mup2*p5mup3*F
     -   (6125)+p3mup2*p4mup4*p5mup3*F(6126)-p1mup3*p1mup4*p5mup2*F(6
     -   127)-p1mup3*p2mup4*p5mup2*F(6128)-p1mup4*p2mup3*p5mup2*F(612
     -   9)-p1mup3*p3mup4*p5mup2*F(6130)-p1mup4*p3mup3*p5mup2*F(6131)
     -   -p1mup3*p4mup4*p5mup2*F(6132)-p2mup3*p2mup4*p5mup2*F(6133)-p
     -   2mup3*p3mup4*p5mup2*F(6134)-p2mup4*p3mup3*p5mup2*F(6135)-p2m
     -   up3*p4mup4*p5mup2*F(6136)-p3mup3*p3mup4*p5mup2*F(6137)-p3mup
     -   3*p4mup4*p5mup2*F(6138)+p2mup3*p3mup2*p4mup4*F(6139)+p3mup2*
     -   p3mup3*p3mup4*F(6140)+p2mup3*p2mup4*p3mup2*F(6141)+p2mup3*p3
     -   mup2*p3mup4*F(6142)+p2mup4*p3mup2*p3mup3*F(6143)+p3mup2*p3mu
     -   p3*p4mup4*F(6144)+p1mup3*p1mup4*p3mup2*F(6145)+p2mup3*p4mup2
     -   *p4mup4*F(6146)+p3mup3*p4mup2*p4mup4*F(6147)+p2mup3*p3mup4*p
     -   4mup2*F(6148)-p2mup4*p3mup2*p4mup3*F(6149)+p2mup2*(p2mup3*p2
     -   mup4*F(6118)+p2mup4*p3mup3*F(6119)+p1mup3*p1mup4*F(6122)+p2m
     -   up4*p4mup3*F(6151))+p3mup3*p3mup4*p4mup2*F(6153)+p1mup2*(p1m
     -   up3*p2mup4*F(6120)+p1mup4*p2mup3*F(6121)+p1mup4*p3mup3*F(612
     -   3)+p1mup4*p4mup3*F(6154))+p1mup3*p1mup4*p4mup2*F(6155)-p3mup
     -   2*p3mup4*p4mup3*F(6156)+p1mup4*p4mup3*p5mup2*F(6157)-p1mup4*
     -   p5mup2*p5mup3*F(6158)+p2mup4*p4mup3*p5mup2*F(6159)-p2mup4*p5
     -   mup2*p5mup3*F(6160)+p3mup4*p4mup3*p5mup2*F(6161)-p3mup4*p5mu
     -   p2*p5mup3*F(6162)+p4mup3*p4mup4*p5mup2*F(6163)-p4mup4*p5mup2
     -   *p5mup3*F(6164)-mup3mup4*p3mup2*F(6165)+mup2mup4*p2mup3*F(61
     -   66)+mup2mup3*p2mup4*F(6167)+mup2mup4*p3mup3*F(6168)+mup2mup3
     -   *p3mup4*F(6169)+mup2mup3*p4mup4*F(6170)+mup2mup4*p1mup3*F(61
     -   71)+mup2mup3*p1mup4*F(6172)+p3mup2*p4mup3*p4mup4*F(6173)+p3m
     -   up4*p4mup2*p4mup3*F(6175)+p4mup2*p4mup3*p4mup4*F(6176)+p4mup
     -   2*p4mup4*p5mup3*F(6177)+2*(p2mup4*p3mup3*p4mup2*F(6150)+p2mu
     -   p3*p2mup4*p4mup2*F(6152)+p2mup4*p4mup2*p4mup3*F(6174)+p2mup4
     -   *p4mup2*p5mup3*F(6178))+p3mup4*p4mup2*p5mup3*F(6179)+mup3mup
     -   4*p4mup2*F(6180)+mup2mup4*p4mup3*F(6181)-mup3mup4*p5mup2*F(6
     -   182)+mup2mup4*p5mup3*F(6183)
       Fa40(549) = p1mup2*(p2mup3*p2mup4*F(6184)+p2mup4*p3mup3*F(6187)
     -   +p1mup3*p1mup4*F(6189)+p2mup4*p4mup3*F(6196))+p1mup3*(p2mup2
     -   *p2mup4*F(6185)+p2mup4*p3mup2*F(6191)+p3mup2*p3mup4*F(6193)+
     -   p3mup2*p4mup4*F(6195)+2*p2mup4*p4mup2*F(6197)+p3mup4*p4mup2*
     -   F(6200)+p4mup2*p4mup4*F(6202))+p1mup4*(p2mup2*(p2mup3*F(6186
     -   )+p3mup3*F(6188)+p4mup3*F(6198))+p3mup2*(p5mup3*F(6190)+p2mu
     -   p3*F(6192)+p3mup3*F(6194)+p4mup3*F(6203))+p4mup2*(p2mup3*F(6
     -   199)+p3mup3*F(6201)+p4mup3*F(6204)+p5mup3*F(6205)))
       Fa40(550) = p5mup4*(p3mup2*p4mup3*F(6210)+p3mup3*p4mup2*F(6211)
     -   -p3mup2*p3mup3*F(6212)+p2mup2*(p3mup3*F(6206)+p2mup3*F(6207)
     -   +p4mup3*F(6208)+p5mup3*F(6213))+p3mup3*p5mup2*F(6215)+p2mup3
     -   *(p4mup2*F(6209)+p5mup2*F(6214)-p3mup2*F(6216))-p4mup2*p4mup
     -   3*F(6217)+p4mup3*p5mup2*F(6218)+p3mup2*p5mup3*F(6219)+p4mup2
     -   *p5mup3*F(6220)+p5mup2*p5mup3*F(6221)+mup2mup3*F(6222))
       Fa40(551) = p2mup2*p4mup4*(p2mup3*F(6231)+p3mup3*F(6232))+p1mup
     -   2*(p1mup3*p5mup4*F(6223)+p2mup3*p5mup4*F(6224)+p3mup3*p5mup4
     -   *F(6226)+p1mup3*p4mup4*F(6227)+p2mup3*p4mup4*F(6228)+p3mup3*
     -   p4mup4*F(6229)+p4mup3*p5mup4*F(6233)+p5mup3*p5mup4*F(6236))+
     -   p1mup3*(p2mup2*p5mup4*F(6225)+p2mup2*p4mup4*F(6230)+p4mup2*p
     -   5mup4*F(6234)-p3mup2*p5mup4*F(6235)+p5mup2*p5mup4*F(6237))
       Fa40(552) = p1mup2*(p3mup3*p3mup4*F(6238)+p1mup3*p3mup4*F(6240)
     -   +p2mup3*p3mup4*F(6241)-p3mup4*p5mup3*F(6244)-p4mup4*p5mup3*F
     -   (6246)-p1mup4*p5mup3*F(6248)-p2mup4*p5mup3*F(6249)+p4mup3*p4
     -   mup4*F(6252)+2*p3mup4*p4mup3*F(6254)+mup3mup4*F(6256))+p2mup
     -   2*(p3mup3*p3mup4*F(6239)+p1mup3*p3mup4*F(6242)+p2mup3*p3mup4
     -   *F(6243)-p3mup4*p5mup3*F(6245)-p4mup4*p5mup3*F(6247)-p1mup4*
     -   p5mup3*F(6250)-p2mup4*p5mup3*F(6251)+p4mup3*p4mup4*F(6253)+2
     -   *p3mup4*p4mup3*F(6255)+mup3mup4*F(6257))
       Fa40(553) = p2mup2*(p2mup3*p2mup4*F(6258)+p2mup4*p3mup3*F(6259)
     -   +p1mup3*p1mup4*F(6262)+p2mup4*p4mup3*F(6264))-p3mup2*p3mup4*
     -   p5mup3*F(6265)+p2mup3*p4mup2*p4mup4*F(6266)+p3mup3*p4mup2*p4
     -   mup4*F(6267)-p2mup3*p3mup2*p4mup4*F(6268)+p2mup3*p3mup4*p4mu
     -   p2*F(6269)+p2mup4*p3mup3*p4mup2*F(6270)-p2mup4*p3mup2*p5mup3
     -   *F(6271)-p3mup2*p3mup3*p3mup4*F(6272)-p3mup2*p4mup4*p5mup3*F
     -   (6273)-p2mup3*p2mup4*p3mup2*F(6274)+p2mup3*p2mup4*p4mup2*F(6
     -   275)-p2mup3*p3mup2*p3mup4*F(6276)-p2mup4*p3mup2*p3mup3*F(627
     -   7)-p3mup2*p3mup3*p4mup4*F(6278)+p3mup3*p3mup4*p4mup2*F(6279)
     -   -p1mup3*p1mup4*p3mup2*F(6280)+p1mup2*(p1mup3*p2mup4*F(6260)+
     -   p1mup4*p2mup3*F(6261)+p1mup4*p3mup3*F(6263)+p1mup4*p4mup3*F(
     -   6281))+p1mup3*p1mup4*p4mup2*F(6282)+p1mup3*p1mup4*p5mup2*F(6
     -   283)+p1mup3*p2mup4*p5mup2*F(6284)+p1mup4*p2mup3*p5mup2*F(628
     -   5)+p1mup3*p3mup4*p5mup2*F(6286)+p1mup4*p3mup3*p5mup2*F(6287)
     -   +p1mup3*p4mup4*p5mup2*F(6288)+p2mup3*p2mup4*p5mup2*F(6289)+p
     -   2mup3*p3mup4*p5mup2*F(6290)+p2mup4*p3mup3*p5mup2*F(6291)+p2m
     -   up3*p4mup4*p5mup2*F(6292)+p3mup3*p3mup4*p5mup2*F(6293)+p3mup
     -   3*p4mup4*p5mup2*F(6294)+p4mup2*p4mup4*p5mup3*F(6295)+p2mup4*
     -   p4mup2*p5mup3*F(6296)+p3mup4*p4mup2*p5mup3*F(6297)-p1mup4*p5
     -   mup2*p5mup3*F(6298)-p2mup4*p5mup2*p5mup3*F(6299)-p3mup4*p5mu
     -   p2*p5mup3*F(6300)-p4mup4*p5mup2*p5mup3*F(6301)+p1mup4*p4mup3
     -   *p5mup2*F(6302)+p2mup4*p3mup2*p4mup3*F(6303)+p2mup4*p4mup2*p
     -   4mup3*F(6304)+p2mup4*p4mup3*p5mup2*F(6305)+p3mup2*p4mup3*p4m
     -   up4*F(6306)+p4mup2*p4mup3*p4mup4*F(6309)+2*(p3mup2*p3mup4*p4
     -   mup3*F(6307)+p3mup4*p4mup2*p4mup3*F(6308)+p3mup4*p4mup3*p5mu
     -   p2*F(6310))+p4mup3*p4mup4*p5mup2*F(6311)+mup3mup4*p3mup2*F(6
     -   312)-mup2mup4*p2mup3*F(6313)+mup2mup3*p2mup4*F(6314)-mup2mup
     -   4*p3mup3*F(6315)+mup2mup3*p3mup4*F(6316)+mup2mup3*p4mup4*F(6
     -   317)-mup2mup4*p1mup3*F(6318)+mup2mup3*p1mup4*F(6319)+mup3mup
     -   4*p4mup2*F(6320)+mup2mup4*p4mup3*F(6321)-mup2mup4*p5mup3*F(6
     -   322)+mup3mup4*p5mup2*F(6323)
       Fa40(554) = p1mup2*(p2mup3*p2mup4*F(6324)+p2mup4*p3mup3*F(6327)
     -   +p1mup3*p1mup4*F(6329)+p2mup4*p4mup3*F(6330))+p1mup3*(p2mup2
     -   *p2mup4*F(6325)-p2mup4*p3mup2*F(6331)+p2mup4*p4mup2*F(6333)-
     -   p3mup2*p3mup4*F(6336)-p3mup2*p4mup4*F(6338)+p3mup4*p4mup2*F(
     -   6339)+p4mup2*p4mup4*F(6342))+p1mup4*(p2mup2*p2mup3*F(6326)+p
     -   2mup2*p3mup3*F(6328)-p2mup3*p3mup2*F(6332)+p2mup2*p4mup3*F(6
     -   334)+p2mup3*p4mup2*F(6335)-p3mup2*p3mup3*F(6337)+p3mup3*p4mu
     -   p2*F(6340)-p3mup2*p5mup3*F(6341)+p4mup2*p5mup3*F(6343)+p3mup
     -   2*p4mup3*F(6344)+p4mup2*p4mup3*F(6345))
       Fa40(555) = p1mup3*p5mup4*F(6347)+p2mup3*p5mup4*F(6348)+p1mup3*
     -   p4mup4*F(6349)+p2mup3*p4mup4*F(6350)+p2mup3*p2mup4*F(6355)+p
     -   1mup3*p1mup4*F(6357)+p1mup3*p2mup4*F(6358)+p1mup4*p2mup3*F(6
     -   359)+p3mup3*(p5mup4*F(6346)+p4mup4*F(6351)+p2mup4*F(6356)+p1
     -   mup4*F(6360))+p4mup3*p5mup4*F(6361)+p3mup4*(2*p4mup3*F(5717)
     -   +p3mup3*F(6352)+p1mup3*F(6353)+p2mup3*F(6354)-p5mup3*F(6362)
     -   )-p4mup4*p5mup3*F(6363)-p1mup4*p5mup3*F(6364)-p2mup4*p5mup3*
     -   F(6365)+p2mup4*p4mup3*F(6366)+p4mup3*p4mup4*F(6367)+p1mup4*p
     -   4mup3*F(6368)+p5mup3*p5mup4*F(6369)+mup3mup4*F(6370)
       Fa40(556) = p1mup2*(p4mup4*F(6371)+p5mup4*F(6373)+p3mup4*F(6375
     -   )+p1mup4*F(6377))+p2mup2*(p4mup4*F(6372)+p5mup4*F(6374)+p3mu
     -   p4*F(6376)+p1mup4*F(6379))+p3mup2*(p4mup4*F(6381)+p5mup4*F(6
     -   382)+p3mup4*F(6383)+p1mup4*F(6384))-p1mup4*p5mup2*F(6387)+p2
     -   mup4*(2*p4mup2*F(5027)+p1mup2*F(6378)+p2mup2*F(6380)+p3mup2*
     -   F(6385)-p5mup2*F(6388))-p5mup2*(p3mup4*F(6389)+p4mup4*F(6390
     -   ))+p4mup2*(p5mup4*F(6386)+p4mup4*F(6391)+p3mup4*F(6392)+p1mu
     -   p4*F(6393))+p5mup2*p5mup4*F(6394)+mup2mup4*F(6395)
       Fa40(557) = mup4mup6*(p1mup3+p2mup3)*p4mup2*F(268)+Fa(170)+16*(
     -   Fa(171)+Fa(172)-F(246)*Fa(173)-Fa(174)-Fa(175)+F(257)*Fa(176
     -   ))+Fa(177)+Fa(178)
       Fa40(558) = 8*(mup4mup6*p2mup2*p2mup3*F(267)+F(267)*Fa(181))+16
     -   *(Fa(182)+Fa(183)+Fa(184))+Fa(185)+Fa(186)
       Fa40(559) = Fa(213)+Fa(214)-8*s45*Fa(215)+Fa(216)-4*Fa(217)
       Fa40(560) = Fa(220)+Fa(221)+Fa(222)+Fa(223)+Fa(224)+8*Fa(225)-F
     -   a(226)-16*(Fa(227)+Fa(228))
       Fa40(561) = -Fa(230)+s45*Fa(231)+Fa(233)+Fa(234)-8*(Fa(232)-Fa(
     -   235)-2*(Fa(236)-Fa(237)-Fa(238)))
       Fa40(562) = Fa(258)-4*Fa(259)-8*(Fa(260)+Fa(261))+Fa(262)
       Fa40(563) = -Fa(263)-4*Fa(264)-8*(Fa(265)-Fa(266))+Fa(267)
       Fa40(564) = Fa(288)-4*Fa(290)+Fa(292)+8*(Fa(289)+Fa(291)+Fa(293
     -   )+Fa(294))
       Fa40(565) = Fa(300)+Fa(301)+Fa(302)+Fa(303)+Fa(304)+Fa(305)-16*
     -   (Fa(306)+Fa(307)+Fa(308)+Fa(309)+Fa(310)+Fa(311)+Fa(312)+Fa(
     -   313))
       Fa40(566) = Fa(319)+8*(Fa(315)+Fa(316)+Fa(317)+Fa(318)+Fa(320))
     -   +16*(Fa(321)+Fa(322)+Fa(323))+Fa(324)+Fa(325)
       Fa40(567) = s12*Fa(327)+Fa(328)+8*Fa(329)+Fa(330)+Fa(331)-Fa(33
     -   2)+16*(Fa(333)+Fa(334)+Fa(335))+Fa(336)+Fa(337)
       Fa40(568) = Fa(339)+Fa(340)+Fa(341)+Fa(342)+Fa(343)
       Fa40(569) = Fa(346)+Fa(347)+4*(Fa(344)+Fa(348))+8*(Fa(345)+Fa(3
     -   49)+Fa(350))
       Fa40(570) = Fa(352)+Fa(353)+Fa(354)+Fa(355)+Fa(356)+8*(Fa(357)+
     -   Fa(358))
       Fa40(571) = Fa(396)+Fa(398)+Fa(399)+Fa(400)+Fa(401)+8*(Fa(397)+
     -   Fa(402))
       Fa40(572) = Fa(403)+8*(Fa(404)+Fa(405)+Fa(406))+Fa(407)+Fa(408)
       Fa40(573) = Fa(409)+Fa(410)+8*(Fa(411)-Fa(412))+Fa(413)
       Fa40(574) = Fa(418)+8*(Fa(417)+Fa(419)+Fa(420))+Fa(421)+Fa(422)
     -   +Fa(428)+Fa(429)+Fa(430)-16*(Fa(423)-Fa(424)-Fa(425)-Fa(426)
     -   -Fa(427)-Fa(431))
       Fa40(575) = Fa(438)+8*(Fa(433)+Fa(434)+Fa(435)+Fa(436)-Fa(437)+
     -   Fa(444))+Fa(445)+Fa(446)+16*(Fa(439)+Fa(440)+Fa(441)+Fa(442)
     -   +Fa(443)+Fa(447))
       Fa40(576) = Fa(450)+16*(Fa(452)-Fa(453)+Fa(455)+Fa(456)+Fa(457)
     -   +Fa(458)+Fa(459))+Fa(460)+8*(Fa(449)+Fa(451)-Fa(454)+Fa(461)
     -   )+Fa(462)+Fa(463)
       Fa40(577) = Fa(465)+Fa(466)+8*(Fa(467)+Fa(468)+Fa(469))
       Fa40(578) = Fa(470)+8*Fa(471)-4*(Fa(472)+Fa(474)+2*(Fa(473)+Fa(
     -   475)+Fa(476)))
       Fa40(579) = Fa(478)-4*Fa(480)+Fa(482)+8*(Fa(479)+Fa(481)-Fa(483
     -   )-Fa(484))
       Fa40(580) = Fa(489)-Fa(493)-2*(Fa(490)-Fa(491)-Fa(492)-Fa(494))
       Fa40(581) = Fa(498)+8*(Fa(499)-Fa(500)+Fa(501))+Fa(502)+Fa(503)
       Fa40(582) = Fa(512)+Fa(513)+8*(Fa(514)+Fa(515))+Fa(516)
       Fa40(583) = Fa(517)+Fa(518)-8*(Fa(519)-Fa(520))+Fa(521)
       Fa40(584) = Fa(530)+Fa(531)+Fa(532)+Fa(533)+Fa(534)
       Fa40(585) = -Fa(540)+2*(Fa(535)-Fa(536)+Fa(537)+Fa(538)+Fa(539)
     -   +Fa(541))-Fa(542)+Fa(543)
       Fa40(586) = Fa(545)+Fa(546)+Fa(547)+Fa(548)+Fa(549)
       Fa40(587) = Fa(550)+Fa(551)+Fa(552)+Fa(553)+Fa(554)
       Return
       End
