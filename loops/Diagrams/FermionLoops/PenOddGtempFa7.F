c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c***********************************************************************************
c************** Calling the Fa functions*************************
       subroutine FaPenOddG7(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1
     -   ,p2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3m
     -   up5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,
     -   p5mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,m
     -   up2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa7)
       IMPLICIT NONE
      Complex*16   p1mup1, p1mup2, p1mup3, p1mup4, p1mup5, p2mup1, 
     -          p2mup2, p2mup3, p2mup4, p2mup5, p3mup1, p3mup2, 
     -          p3mup3, p3mup4, p3mup5, p4mup1, p4mup2, p4mup3, 
     -          p4mup4, p4mup5, p5mup1, p5mup2, p5mup3, p5mup4, p5mup5
       Complex*16   mup1mup2, mup1mup3, mup1mup4, mup1mup5, mup2mup3, 
     -          mup2mup4, mup2mup5, mup3mup4, mup3mup5, mup4mup5
        common/FPenOddG/F
      Complex*16 EpsF(126)
       common/EPSFPenOddG/EPSF
       COMMON/InvariantsPenOddG/p1sq,p2sq,p3sq,p4sq,p5sq,s12,s23,s34,s
     -   45,s234,m1,m1sq,m3,m3sq,m4,m4sq
       Real*8   p1sq 
        Real*8   p2sq 
       Real*8   p3sq 
       Real*8   p4sq 
       Real*8   p5sq 
       Real*8   s12 
       Real*8   s23 
       Real*8   s34 
       Real*8   s45 
       Real*8   s234,m1,m1sq,m3,m3sq,m4,m4sq
      Complex*16 Fa(4390),F(7459)
      Real*8 P(1989) 
        COMMON/PPenOddG/P
       Complex*16 Fa7(1753:2044)
       COMMON/FaPenOddG/Fa
       Fa7(1753) = mup4mup5*F(5322)-2*(p1mup5*(2*(p2mup4*F(359)+p3mup4
     -   *F(888))+p4mup4*F(7441))-p1mup4*(2*(p2mup5*F(359)+p3mup5*F(8
     -   88))+p4mup5*F(7441)))+4*((p3mup5*p4mup4-p3mup4*p4mup5)*F(744
     -   2)+(p2mup5*p4mup4-p2mup4*p4mup5)*F(7448))
       Fa7(1754) = -2*(p1mup5*p3mup4*F(352)-p1mup4*p3mup5*F(352))+mup4
     -   mup5*F(5177)+4*((p3mup5*p4mup4-p3mup4*p4mup5)*F(7444)-(p2mup
     -   5*p3mup4-p2mup4*p3mup5)*F(7448))
       Fa7(1755) = p2mup5*F(5221)+p3mup5*F(5222)+(2*p1mup5*F(276)+p4mu
     -   p5*F(352))*P(1)
       Fa7(1756) = p1mup5*F(5196)+p3mup5*F(5197)+2*p2mup5*F(5198)+p4mu
     -   p5*F(5199)
       Fa7(1757) = -2*p3mup5*F(5189)+p1mup5*F(5200)+p2mup5*F(5201)+p4m
     -   up5*F(5202)
       Fa7(1758) = p3mup5*F(5174)+p2mup5*F(5225)+(p1mup5*F(352)+2*p4mu
     -   p5*F(7444))*P(1)
       Fa7(1759) = p2mup1*F(5175)+p3mup1*F(5194)-p4mup1*F(7441)*P(1)
       Fa7(1760) = p4mup1*F(5177)+p1mup1*F(5178)+p3mup1*F(5195)
       Fa7(1761) = p4mup1*F(368)+p1mup1*F(5192)+p2mup1*F(5193)
       Fa7(1762) = p3mup1*F(374)+p2mup1*F(5180)+p1mup1*F(7441)*P(1)
       Fa7(1763) = p1mup4*F(5183)+p3mup4*F(5184)+2*p2mup4*F(5185)+p4mu
     -   p4*F(5186)
       Fa7(1764) = -(p1mup4*F(5187))+p2mup4*F(5188)+2*p3mup4*F(5189)+p
     -   4mup4*F(5190)
       Fa7(1765) = 2*(p2mup5*p4mup4-p2mup4*p4mup5)*F(352)-4*(p1mup5*p2
     -   mup4*F(276)-p1mup4*p2mup5*F(276)-(p2mup5*p3mup4-p2mup4*p3mup
     -   5)*F(888))+mup4mup5*F(5194)
       Fa7(1766) = 2*p1mup2*F(5324)+p2mup2*F(5325)+p3mup2*F(5326)+p4mu
     -   p2*F(5327)
       Fa7(1767) = p1mup2*F(5328)+p3mup2*F(5329)+2*p2mup2*F(5330)+p4mu
     -   p2*F(5331)
       Fa7(1768) = p1mup2*F(5332)+p2mup2*F(5333)+2*p3mup2*F(5334)+p4mu
     -   p2*F(5335)
       Fa7(1769) = p1mup2*F(5336)+p2mup2*F(5337)+p3mup2*F(5338)+2*p4mu
     -   p2*F(5339)
       Fa7(1770) = p2mup4*F(5181)+p3mup4*F(5182)-(2*p1mup4*F(276)+p4mu
     -   p4*F(352))*P(1)
       Fa7(1771) = p3mup4*F(5140)+p2mup4*F(5191)-(p1mup4*F(352)+2*p4mu
     -   p4*F(7444))*P(1)
       Fa7(1772) = -2*(p3mup5*p4mup4-p3mup4*p4mup5)*F(352)+4*(p1mup5*p
     -   3mup4*F(276)-p1mup4*p3mup5*F(276)+(p2mup5*p3mup4-p2mup4*p3mu
     -   p5)*F(359))+mup4mup5*F(5123)
       Fa7(1773) = mup4mup5*F(5340)-2*(p1mup5*(2*(p2mup4*F(359)+p3mup4
     -   *F(888))+p4mup4*F(7441))-p1mup4*(2*(p2mup5*F(359)+p3mup5*F(8
     -   88))+p4mup5*F(7441)))+4*((p3mup5*p4mup4-p3mup4*p4mup5)*F(744
     -   2)+(p2mup5*p4mup4-p2mup4*p4mup5)*F(7448))
       Fa7(1774) = -2*(p1mup5*p3mup4*F(352)-p1mup4*p3mup5*F(352))+mup4
     -   mup5*F(5125)+4*((p3mup5*p4mup4-p3mup4*p4mup5)*F(7444)-(p2mup
     -   5*p3mup4-p2mup4*p3mup5)*F(7448))
       Fa7(1775) = p2mup5*F(5169)+p3mup5*F(5170)+(2*p1mup5*F(276)+p4mu
     -   p5*F(352))*P(1)
       Fa7(1776) = p1mup5*F(5145)+p3mup5*F(5146)+2*p2mup5*F(5147)+p4mu
     -   p5*F(5148)
       Fa7(1777) = 2*p3mup5*F(5137)+p4mup5*F(5138)+p1mup5*F(5149)+p2mu
     -   p5*F(5150)
       Fa7(1778) = p2mup5*F(5173)+p3mup5*F(5174)+(p1mup5*F(352)+2*p4mu
     -   p5*F(7444))*P(1)
       Fa7(1779) = p2mup1*F(5123)+p3mup1*F(5143)-p4mup1*F(7441)*P(1)
       Fa7(1780) = p4mup1*F(5125)+p1mup1*F(5126)+p3mup1*F(5144)
       Fa7(1781) = p4mup1*F(368)+p1mup1*F(5141)+p2mup1*F(5142)
       Fa7(1782) = p3mup1*F(374)+p2mup1*F(5128)+p1mup1*F(7441)*P(1)
       Fa7(1783) = p1mup4*F(5131)+p3mup4*F(5132)+2*p2mup4*F(5133)+p4mu
     -   p4*F(5134)
       Fa7(1784) = p1mup4*F(5135)-p2mup4*F(5136)+2*p3mup4*F(5137)+p4mu
     -   p4*F(5138)
       Fa7(1785) = 2*(p2mup5*p4mup4-p2mup4*p4mup5)*F(352)-4*(p1mup5*p2
     -   mup4*F(276)-p1mup4*p2mup5*F(276)-(p2mup5*p3mup4-p2mup4*p3mup
     -   5)*F(888))+mup4mup5*F(5143)
       Fa7(1786) = 2*p1mup2*F(5342)+p2mup2*F(5343)+p3mup2*F(5344)+p4mu
     -   p2*F(5345)
       Fa7(1787) = p1mup2*F(5346)+p3mup2*F(5347)+2*p2mup2*F(5348)+p4mu
     -   p2*F(5349)
       Fa7(1788) = p3mup2*F(5352)+M1*M4*(p1mup2*F(5350)+p2mup2*F(5351)
     -   +p4mup2*F(5353))
       Fa7(1789) = p1mup2*F(5354)+p2mup2*F(5355)+p3mup2*F(5356)+2*p4mu
     -   p2*F(5357)
       Fa7(1790) = p2mup4*F(5129)+p3mup4*F(5130)-(2*p1mup4*F(276)+p4mu
     -   p4*F(352))*P(1)
       Fa7(1791) = p2mup4*F(5139)+p3mup4*F(5140)-(p1mup4*F(352)+2*p4mu
     -   p4*F(7444))*P(1)
       Fa7(1792) = 8*p2mup4*p2mup5*F(4991)+mup4mup5*F(5360)+p1mup5*(8*
     -   p2mup4*F(4984)+p3mup4*F(5361)+4*p4mup4*F(5362))+p2mup5*p3mup
     -   4*F(5364)+p3mup4*p4mup5*F(5365)+2*p3mup5*p4mup4*F(5366)+4*(p
     -   1mup4*(p3mup5*F(4983)+p4mup5*F(4986)+2*(p1mup5*F(4982)+p2mup
     -   5*F(4990)))+p2mup4*p3mup5*F(4992)+p3mup4*p3mup5*F(4996)+p2mu
     -   p5*p4mup4*F(5016)+p2mup4*p4mup5*F(5363)+p4mup4*p4mup5*F(5367
     -   ))
       Fa7(1793) = 8*p4mup4*p4mup5*F(5011)+4*p2mup4*(p4mup5*F(5004)+2*
     -   (p3mup5*F(5009)+p2mup5*F(5013)))+mup4mup5*F(5368)+p1mup4*p2m
     -   up5*F(5369)+p1mup4*p3mup5*F(5370)+p1mup4*p4mup5*F(5371)+p3mu
     -   p4*p4mup5*F(5372)+p3mup5*p4mup4*F(5373)+p1mup5*(4*(p2mup4*F(
     -   5000)+p1mup4*F(5001))+2*p4mup4*F(5008)+p3mup4*F(5374))+p3mup
     -   4*p3mup5*F(5375)+p2mup5*p4mup4*F(5376)+p2mup5*p3mup4*F(5377)
       Fa7(1794) = -4*p2mup4*p3mup5*F(5026)+mup4mup5*F(5378)+p1mup4*(-
     -   2*p3mup5*F(5018)+4*(p2mup5*F(4989)+p1mup5*F(5380)+p4mup5*F(5
     -   381)))+p3mup4*p4mup5*F(5383)+p3mup5*p4mup4*F(5384)+p1mup5*(p
     -   3mup4*F(5379)+4*(p2mup4*F(5017)+p4mup4*F(5385)))+8*(p2mup4*p
     -   4mup5*F(5382)+p3mup4*p3mup5*F(5386)+p4mup4*p4mup5*F(5387)+p2
     -   mup5*p4mup4*F(5388)+p2mup4*p2mup5*F(5389))+p2mup5*p3mup4*F(5
     -   390)
       Fa7(1795) = -2*p1mup5*F(5035)+p2mup5*F(5391)+p3mup5*F(5392)+p4m
     -   up5*F(5393)
       Fa7(1796) = -2*p2mup5*F(5039)+p1mup5*F(5394)+p4mup5*F(5395)+p3m
     -   up5*F(5396)
       Fa7(1797) = 2*p3mup5*F(5043)+p1mup5*F(5397)+p2mup5*F(5398)+p4mu
     -   p5*F(5399)
       Fa7(1798) = -2*p4mup5*F(5047)+p1mup5*F(5400)+p2mup5*F(5401)+p3m
     -   up5*F(5402)
       Fa7(1799) = p2mup1*F(4978)+p3mup1*F(5053)+p4mup1*F(5057)
       Fa7(1800) = p1mup1*F(5048)+p3mup1*F(5054)+p4mup1*F(5056)
       Fa7(1801) = p1mup1*F(5049)+p2mup1*F(5051)+p4mup1*F(5058)
       Fa7(1802) = p1mup1*F(5050)+p2mup1*F(5052)+p3mup1*F(5055)
       Fa7(1803) = p1mup4*F(5403)+p2mup4*F(5404)+p3mup4*F(5405)+p4mup4
     -   *F(5406)
       Fa7(1804) = p1mup4*F(5407)+p2mup4*F(5408)+p3mup4*F(5409)+p4mup4
     -   *F(5410)
       Fa7(1805) = -2*(4*p3mup4*p3mup5*F(4563)+p2mup4*p4mup5*F(5070)+p
     -   2mup5*p4mup4*F(5071)+2*(p3mup4*p4mup5*F(4560)+p3mup5*p4mup4*
     -   F(4561)+p4mup4*p4mup5*F(4564)+p1mup4*(p4mup5*F(4558)+2*(p1mu
     -   p5*F(4556)+p3mup5*F(4562))+p2mup5*F(5068))+p1mup5*(2*p3mup4*
     -   F(4557)+p4mup4*F(4559)-p2mup4*F(5069))+p2mup4*p2mup5*F(5072)
     -   -p2mup4*p3mup5*F(5073)+p2mup5*p3mup4*F(5074)))+mup4mup5*F(54
     -   11)
       Fa7(1806) = 8*(p3mup4*p4mup5*F(4585)+p4mup4*p4mup5*F(4587)+p3mu
     -   p5*p4mup4*F(4589)+p2mup4*p2mup5*F(4590)+p3mup4*p3mup5*F(4590
     -   ))+2*p1mup4*(2*(p3mup5*F(4561)+p1mup5*F(4581)+p4mup5*F(4582)
     -   )+p2mup5*F(5076))+p1mup5*(4*(p3mup4*F(4580)+p4mup4*F(4586))-
     -   2*p2mup4*F(5077))-4*(p2mup4*p4mup5*F(5078)-p2mup5*p4mup4*F(5
     -   079)+p2mup4*p3mup5*F(5080)-p2mup5*p3mup4*F(5081))+mup4mup5*F
     -   (5412)
       Fa7(1807) = p1mup2*F(5414)+p2mup2*F(5415)+p3mup2*F(5416)+p4mup2
     -   *F(5417)
       Fa7(1808) = p1mup2*F(5418)+p3mup2*F(5419)+p4mup2*F(5420)+2*p2mu
     -   p2*F(5421)
       Fa7(1809) = p1mup2*F(5422)+p2mup2*F(5423)+p4mup2*F(5424)+p3mup2
     -   *F(5425)
       Fa7(1810) = p1mup2*F(5427)+p2mup2*F(5428)+p3mup2*F(5429)+p4mup2
     -   *F(5430)
       Fa7(1811) = p2mup4*F(5433)+p3mup4*F(5434)+p1mup4*F(5435)+p4mup4
     -   *F(5436)
       Fa7(1812) = p2mup4*F(5440)+p1mup4*F(5441)+p3mup4*F(5442)-2*p4mu
     -   p4*F(5443)
       Fa7(1813) = mup4mup5*F(5057)-4*(p1mup5*p3mup4*F(4528)-p1mup4*p3
     -   mup5*F(4663)+p3mup4*p4mup5*F(4668)+p3mup4*p3mup5*F(4670)+p3m
     -   up5*p4mup4*F(4672)+p1mup4*p2mup5*F(5116)-p2mup5*p4mup4*F(512
     -   2)-p1mup5*p2mup4*F(5447)-p2mup4*p4mup5*F(5448)-p2mup4*p2mup5
     -   *F(5449))+2*p2mup4*p3mup5*F(5450)+p2mup5*p3mup4*F(5451)
       Fa7(1814) = 2*p3mup5*p4mup2*F(5457)+p3mup2*p4mup5*F(5459)+4*(p1
     -   mup2*(p3mup5*F(5452)+p4mup5*F(5453))+p1mup5*(p4mup2*F(5447)+
     -   p3mup2*F(5454))+p2mup5*p3mup2*F(5455)+p2mup5*p4mup2*F(5456)+
     -   p3mup2*p3mup5*F(5458)+p4mup2*p4mup5*F(5460)+p2mup2*p3mup5*F(
     -   5461)-p2mup2*p4mup5*F(5462))
       Fa7(1815) = p1mup2*(2*(p1mup5*F(5463)+p2mup5*F(5464)+p3mup5*F(5
     -   465))+p4mup5*F(5466))+p1mup5*(2*(p2mup2*F(5467)+p3mup2*F(547
     -   2))-p4mup2*F(5473))+4*(p2mup2*p2mup5*F(5468)+p3mup2*p3mup5*F
     -   (5469)+p4mup2*p4mup5*F(5471)+(p2mup5*p3mup2+p2mup2*p3mup5)*F
     -   (5474))-2*(p3mup5*p4mup2*F(5470)+p3mup2*p4mup5*F(5470)+(p2mu
     -   p5*p4mup2+p2mup2*p4mup5)*F(5475))
       Fa7(1816) = 4*(p1mup2*p1mup5*F(5476)+p1mup2*p2mup5*F(5477)+p1mu
     -   p2*p4mup5*F(5479))+p2mup5*p3mup2*F(5482)+p2mup2*(4*p1mup5*F(
     -   5480)+8*(p4mup5*F(5389)+p2mup5*F(5481))+p3mup5*F(5482))-2*(p
     -   1mup2*p3mup5*F(5478)+p1mup5*p3mup2*F(5483))-8*p3mup2*p3mup5*
     -   F(5484)+p3mup2*p4mup5*F(5485)+p4mup2*(8*(p4mup5*F(5382)+p2mu
     -   p5*F(5389))+p3mup5*F(5485)+4*p1mup5*F(5486))
       Fa7(1817) = p2mup5*(2*(p3mup2*F(5472)+p1mup2*F(5488))+p4mup2*F(
     -   5493))+p2mup2*(2*(p3mup5*F(5465)+p1mup5*F(5489)+p2mup5*F(549
     -   2))+p4mup5*F(5495))+2*(-(p1mup2*p4mup5*F(5490))+p1mup5*p4mup
     -   2*F(5491)+p3mup5*p4mup2*F(5494)+p4mup2*p4mup5*F(5497)-p3mup2
     -   *p4mup5*F(5498)+2*(p1mup2*p1mup5*F(5487)+p1mup2*p3mup5*F(549
     -   6)+p3mup2*p3mup5*F(5499)+p1mup5*p3mup2*F(5500)))
       Fa7(1818) = p2mup2*(4*(p4mup5*F(5479)+p1mup5*F(5502)+p2mup5*F(5
     -   505))+p3mup5*F(5506))+2*p2mup5*p3mup2*F(5507)+4*(p2mup5*p4mu
     -   p2*F(5486)+p1mup5*p3mup2*F(5504)+p3mup2*p3mup5*F(5508)-p1mup
     -   2*(p2mup5*F(5501)-p3mup5*F(5503)+2*(p1mup5*F(4982)+p4mup5*F(
     -   5509)))+p3mup2*p4mup5*F(5510)+p3mup5*p4mup2*F(5511))+8*(p4mu
     -   p2*p4mup5*F(5512)+p1mup5*p4mup2*F(5513))
       Fa7(1819) = p1mup2*(p2mup5*F(5515)+2*(p1mup5*F(5514)+p3mup5*F(5
     -   516)+p4mup5*F(5517)))+p1mup5*p2mup2*F(5519)+2*(p2mup2*p3mup5
     -   *F(5518)+p1mup5*p3mup2*F(5520)+p1mup5*p4mup2*F(5524)+p2mup2*
     -   p4mup5*F(5525)+p2mup5*(p3mup2*F(5518)+2*p2mup2*F(5521)+p4mup
     -   2*F(5525)))+4*(p3mup2*p3mup5*F(5522)+p4mup2*p4mup5*F(5523)+(
     -   p3mup5*p4mup2+p3mup2*p4mup5)*F(5526))
       Fa7(1820) = p3mup5*(p4mup2*F(1622)-2*(p2mup2*F(1335)+p1mup2*F(1
     -   646)))+p3mup2*(p4mup5*F(352)+2*(p1mup5*F(276)+p2mup5*F(359)+
     -   p3mup5*F(1843)))
       Fa7(1821) = p1mup2*(2*(p1mup5*F(668)+p2mup5*F(670))+p4mup5*F(10
     -   33))-p1mup5*p4mup2*F(1034)+4*p2mup2*p2mup5*F(1036)+2*(p1mup5
     -   *p2mup2*F(676)+(p2mup5*p4mup2+p2mup2*p4mup5)*F(1037)+p3mup5*
     -   (p4mup2*F(1046)+2*p2mup2*F(1060)+p1mup2*F(1864))+p3mup2*(p4m
     -   up5*F(1046)+2*(p3mup5*F(304)+p2mup5*F(1060))+p1mup5*F(1873))
     -   +2*p4mup2*p4mup5*F(7454))
       Fa7(1822) = p1mup5*p3mup2*F(352)+p1mup2*p3mup5*F(352)+2*(2*p3mu
     -   p2*p3mup5*F(7442)+p3mup5*p4mup2*F(7444)+p3mup2*p4mup5*F(7444
     -   )+(p2mup5*p3mup2+p2mup2*p3mup5)*F(7448))
       Fa7(1823) = p2mup5*p4mup2*F(1024)+p2mup2*p4mup5*F(1024)+4*(p3mu
     -   p5*p4mup2*F(1047)+p3mup2*p3mup5*F(1059))+2*(p1mup2*(2*(p1mup
     -   5*F(130)+p3mup5*F(278))+p4mup5*F(673)+p2mup5*F(962))+p2mup2*
     -   p2mup5*F(1022)+p3mup2*p4mup5*F(1026)+p4mup2*p4mup5*F(1027)+p
     -   1mup5*(p2mup2*F(964)+2*(p4mup2*F(338)+p3mup2*F(1028)))+p2mup
     -   2*p3mup5*F(1864)+p2mup5*p3mup2*F(1873))
       Fa7(1824) = -(p2mup5*p3mup2*F(7441))+p3mup5*(p2mup2*F(1626)+2*(
     -   p1mup2*F(1012)+p3mup2*F(7442)+p4mup2*F(7444)))
       Fa7(1825) = p3mup5*(2*(p1mup2*F(1325)+p2mup2*F(1623))+p4mup2*F(
     -   1645))+p3mup2*(p4mup5*F(441)+2*(p1mup5*F(276)+p2mup5*F(359)+
     -   p3mup5*F(1843)))
       Fa7(1826) = p1mup2*p4mup5*F(1070)+p1mup5*p4mup2*F(1071)-4*(p3mu
     -   p2*p3mup5*F(304)+p2mup2*p2mup5*F(1036)+p2mup5*p3mup2*F(1060)
     -   +p2mup2*p3mup5*F(1060)-p4mup2*p4mup5*F(1076))-2*(p1mup2*p1mu
     -   p5*F(668)+p1mup2*p2mup5*F(670)+p1mup5*p2mup2*F(676)-p2mup5*p
     -   4mup2*F(1073)-p2mup2*p4mup5*F(1073)+p1mup2*p3mup5*F(1864)+p1
     -   mup5*p3mup2*F(1873)-(p3mup5*p4mup2+p3mup2*p4mup5)*F(5527))
       Fa7(1827) = 2*(p2mup5*p3mup2*F(434)+p2mup2*p3mup5*F(434)+p3mup5
     -   *p4mup2*F(436)+p3mup2*p4mup5*F(436))+p1mup5*p3mup2*F(438)+p1
     -   mup2*p3mup5*F(438)+4*p3mup2*p3mup5*F(439)
       Fa7(1828) = 4*(p3mup5*p4mup2*F(1047)+p3mup2*p3mup5*F(1059))-p2m
     -   up2*p4mup5*F(1062)+p2mup5*p4mup2*F(1063)+2*(p2mup2*p2mup5*F(
     -   1022)+p1mup5*(p2mup2*F(964)+2*(p4mup2*F(338)+p3mup2*F(1028))
     -   )+p1mup2*(2*(p1mup5*F(130)+p3mup5*F(278))+p2mup5*F(962)-p4mu
     -   p5*F(1064))-p3mup2*p4mup5*F(1066)-p4mup2*p4mup5*F(1068)+p2mu
     -   p2*p3mup5*F(1864)+p2mup5*p3mup2*F(1873))
       Fa7(1829) = -(p2mup5*p3mup2*F(1095))+p3mup5*(p2mup2*F(1653)+2*(
     -   p4mup2*F(436)+p3mup2*F(439)+p1mup2*F(1654)))
       Fa7(1830) = p3mup5*(p4mup2*F(1645)-2*(p2mup2*F(1335)+p1mup2*F(1
     -   646)))+p3mup2*(-(p4mup5*F(429))+2*(p1mup5*F(276)+p2mup5*F(35
     -   9)+p3mup5*F(1843)))
       Fa7(1831) = 4*p2mup2*p2mup5*F(1036)-p1mup5*p4mup2*F(1071)+p1mup
     -   2*(2*(p1mup5*F(668)+p2mup5*F(670))+p4mup5*F(1099))+2*(p1mup5
     -   *p2mup2*F(676)+p2mup5*p4mup2*F(1100)+p2mup2*p4mup5*F(1100)+2
     -   *(p3mup2*p3mup5*F(304)+p2mup5*p3mup2*F(1060)+p2mup2*p3mup5*F
     -   (1060)+p4mup2*p4mup5*F(1101))+p1mup2*p3mup5*F(1864)+p1mup5*p
     -   3mup2*F(1873)+(p3mup5*p4mup2+p3mup2*p4mup5)*F(5528))
       Fa7(1832) = 4*(p3mup5*p4mup2*F(1047)+p3mup2*p3mup5*F(1059))+p2m
     -   up5*p4mup2*F(1063)+p2mup2*p4mup5*F(1096)+2*(p2mup2*p2mup5*F(
     -   1022)+p1mup5*(p2mup2*F(964)+2*(p4mup2*F(338)+p3mup2*F(1028))
     -   )+p1mup2*(2*(p1mup5*F(130)+p3mup5*F(278))+p2mup5*F(962)-p4mu
     -   p5*F(1064))+p3mup2*p4mup5*F(1097)+p4mup2*p4mup5*F(1098)+p2mu
     -   p2*p3mup5*F(1864)+p2mup5*p3mup2*F(1873))
       Fa7(1833) = p3mup5*(p4mup2*F(1622)+2*(p1mup2*F(1325)+p2mup2*F(1
     -   623)))+p3mup2*(p4mup5*F(352)+2*(p1mup5*F(276)+p2mup5*F(359)+
     -   p3mup5*F(1843)))
       Fa7(1834) = p3mup5*p4mup2*F(5530)+2*p3mup2*p4mup5*F(5531)+4*(p1
     -   mup5*p4mup2*F(5117)-p1mup2*(p3mup5*F(5452)+p4mup5*F(5453))-p
     -   1mup5*p3mup2*F(5454)-p2mup5*p3mup2*F(5455)-p3mup2*p3mup5*F(5
     -   458)-p4mup2*p4mup5*F(5460)+p2mup2*p4mup5*F(5462)+p2mup5*p4mu
     -   p2*F(5529)+p2mup2*p3mup5*F(5532))
       Fa7(1835) = p1mup5*(4*(p2mup2*F(5467)+p3mup2*F(5472))-2*p4mup2*
     -   F(5473))+8*(p2mup2*p2mup5*F(5468)+p3mup2*p3mup5*F(5469)+p4mu
     -   p2*p4mup5*F(5471)+p2mup5*p3mup2*F(5474)+p2mup2*p3mup5*F(5474
     -   ))+p1mup2*(4*(p1mup5*F(5463)+p2mup5*F(5464)+p3mup5*F(5465))+
     -   p4mup5*F(5533))+p3mup5*p4mup2*F(5534)+p3mup2*p4mup5*F(5534)+
     -   (p2mup5*p4mup2+p2mup2*p4mup5)*F(5535)
       Fa7(1836) = p1mup2*p3mup5*F(5478)+p1mup5*p3mup2*F(5483)+4*p3mup
     -   2*p3mup5*F(5484)+2*(p1mup2*p1mup5*F(5536)+p1mup2*p2mup5*F(55
     -   37)+p1mup2*p4mup5*F(5538)+p2mup5*p3mup2*F(5541)+p2mup2*(p1mu
     -   p5*F(5539)+2*(p4mup5*F(5030)+p2mup5*F(5540))+p3mup5*F(5541))
     -   +p3mup2*p4mup5*F(5542)+p4mup2*(2*(p4mup5*F(5022)+p2mup5*F(50
     -   30))+p3mup5*F(5542)+p1mup5*F(5543)))
       Fa7(1837) = 2*p2mup5*(2*(p3mup2*F(5472)+p1mup2*F(5488))+p4mup2*
     -   F(5493))+p2mup2*(4*(p3mup5*F(5465)+p1mup5*F(5489)+p2mup5*F(5
     -   492))+p4mup5*F(5544))+4*(-(p1mup2*p4mup5*F(5490))+p1mup5*p4m
     -   up2*F(5491)+p3mup5*p4mup2*F(5494)+p4mup2*p4mup5*F(5497)+2*(p
     -   1mup2*p1mup5*F(5487)+p1mup2*p3mup5*F(5496)+p3mup2*p3mup5*F(5
     -   499)+p1mup5*p3mup2*F(5500))+p3mup2*p4mup5*F(5545))
       Fa7(1838) = -2*(p2mup5*p3mup2*F(5507)-p2mup2*p3mup5*F(5548))-4*
     -   (p3mup2*p3mup5*F(5508)-p1mup2*(p2mup5*F(5501)-p3mup5*F(5503)
     -   +2*(p1mup5*F(4982)+p4mup5*F(5509)))+p3mup2*p4mup5*F(5510)-p2
     -   mup2*p4mup5*F(5538)-p2mup5*p4mup2*F(5543)-p2mup2*p2mup5*F(55
     -   47)-p3mup5*p4mup2*F(5549))+8*p4mup2*p4mup5*F(5550)+p1mup5*(-
     -   4*(p3mup2*F(5504)-p2mup2*F(5546))+8*p4mup2*F(5551))
       Fa7(1839) = p1mup2*(8*(p1mup4*F(5552)+p2mup4*F(5553))+p3mup4*F(
     -   5554)+4*p4mup4*F(5556))+p3mup4*p4mup2*F(5560)-8*p2mup2*p2mup
     -   4*F(5561)+p2mup2*p3mup4*F(5563)+p3mup2*p3mup4*F(5565)-2*p3mu
     -   p2*p4mup4*F(5566)-4*(p2mup4*p3mup2*F(5558)+p2mup4*p4mup2*F(5
     -   559)+p1mup4*(p3mup2*F(5555)+p4mup2*F(5557)+2*p2mup2*F(5562))
     -   +p2mup2*p4mup4*F(5564)-p4mup2*p4mup4*F(5567))
       Fa7(1840) = p1mup2*(p2mup4*F(5569)+8*(p1mup4*F(5568)+p3mup4*F(5
     -   571))+p4mup4*F(5572))+2*(p2mup4*p4mup2*F(5575)+4*p3mup2*p3mu
     -   p4*F(5577)+p2mup2*p4mup4*F(5581)+2*(p2mup2*p2mup4*F(5574)+p3
     -   mup4*p4mup2*F(5576)+p1mup4*(p2mup2*F(5570)+p4mup2*F(5573)+2*
     -   p3mup2*F(5578))+p2mup2*p3mup4*F(5579)+p2mup4*p3mup2*F(5580)+
     -   p3mup2*p4mup4*F(5582)+p4mup2*p4mup4*F(5583)))
       Fa7(1841) = p1mup2*(4*(p1mup4*F(5584)+p2mup4*F(5585))+p3mup4*F(
     -   5588)+2*p4mup4*F(5590))+p1mup4*p4mup2*F(5591)+p2mup4*p4mup2*
     -   F(5593)+8*(p2mup4*p3mup2*F(5592)+p3mup2*p3mup4*F(5594)+p4mup
     -   2*p4mup4*F(5595))+p3mup4*p4mup2*F(5597)+4*(p1mup4*p3mup2*F(5
     -   589)+p3mup2*p4mup4*F(5596)+p2mup2*(p1mup4*F(5587)+2*(p2mup4*
     -   F(5586)+p3mup4*F(5598))+p4mup4*F(5599)))
       Fa7(1842) = p3mup2*p3mup4*F(5606)+p2mup2*p3mup4*F(5607)-2*p2mup
     -   4*p3mup2*F(5608)+p1mup2*(4*p2mup4*F(5601)+p3mup4*F(5603)+8*(
     -   p1mup4*F(5600)-p4mup4*F(5609)))+4*(p2mup2*p2mup4*F(5605)+p1m
     -   up4*(p2mup2*F(5602)-p3mup2*F(5604)-2*p4mup2*F(5609))-p2mup2*
     -   p4mup4*F(5610)-p2mup4*p4mup2*F(5611)+p3mup2*p4mup4*F(5612))+
     -   p3mup4*p4mup2*F(5613)-8*p4mup2*p4mup4*F(5614)
       Fa7(1843) = 2*p1mup4*p3mup2*F(5621)+p1mup2*(p3mup4*F(5617)+4*(p
     -   1mup4*F(5615)+p2mup4*F(5616)+p4mup4*F(5624)))+p3mup4*p4mup2*
     -   F(5627)+8*(p3mup2*p3mup4*F(5620)-p2mup4*p4mup2*F(5626)-p4mup
     -   2*p4mup4*F(5628))+4*(p2mup4*p3mup2*F(5623)-p1mup4*p4mup2*F(5
     -   625)+p3mup2*p4mup4*F(5629))+p2mup2*(4*p1mup4*F(5619)+p3mup4*
     -   F(5622)+8*(p2mup4*F(5618)+p4mup4*F(5630)))
       Fa7(1844) = p1mup2*(p2mup4*F(5632)+2*(p1mup4*F(5631)+p3mup4*F(5
     -   633)+p4mup4*F(5639)))+p2mup2*(4*p2mup4*F(5634)+p1mup4*F(5636
     -   )+2*(p3mup4*F(5638)+p4mup4*F(5644)))+2*(p4mup2*(p1mup4*F(564
     -   0)+p2mup4*F(5641)+2*(p3mup4*F(5642)+p4mup4*F(5643)))+p3mup2*
     -   (p1mup4*F(5637)+p2mup4*F(5638)+2*(p3mup4*F(5635)+p4mup4*F(56
     -   45))))
       Fa7(1845) = -2*mup2mup4*(p1mup5*F(5454)+p2mup5*F(5455)+p3mup5*F
     -   (5648))+mup2mup4*p4mup5*F(5649)+mup2mup5*(2*(p1mup4*F(5646)+
     -   p2mup4*F(5647)+p3mup4*F(5650))+p4mup4*F(5651))+mup4mup5*(-2*
     -   (p1mup2*F(5452)-p2mup2*F(5532)-p3mup2*F(5652))-p4mup2*F(5653
     -   ))
       Fa7(1846) = mup4mup5*p1mup2*F(5654)-mup2mup5*p1mup4*F(5654)+mup
     -   2mup4*p1mup5*F(5654)+2*(mup4mup5*p2mup2*F(5655)-mup2mup5*p2m
     -   up4*F(5655)+mup2mup4*p2mup5*F(5655)+mup4mup5*p3mup2*F(5656)-
     -   mup2mup5*p3mup4*F(5656)+mup2mup4*p3mup5*F(5656)+(mup4mup5*p4
     -   mup2-mup2mup5*p4mup4+mup2mup4*p4mup5)*F(5657))
       Fa7(1847) = mup2mup5*(4*(p1mup4*F(5116)-p4mup4*F(5122)+p2mup4*F
     -   (5658))+p3mup4*F(5659))+mup2mup4*(4*(p1mup5*F(5117)+p4mup5*F
     -   (5118)+p2mup5*F(5529))+p3mup5*F(5660))+2*mup4mup5*(p3mup2*F(
     -   5661)-2*(p1mup2*F(5453)-p2mup2*F(5462)+p4mup2*F(5662)))
       Fa7(1848) = mup2mup5*(2*(p1mup4*F(4663)-p4mup4*F(4672)+p3mup4*F
     -   (4961))+p2mup4*F(5664))-mup2mup4*(2*(p1mup5*F(4528)+p3mup5*F
     -   (4529)+p4mup5*F(4668))+p2mup5*F(5665))+mup4mup5*(p2mup2*F(56
     -   66)+2*(p1mup2*F(5663)+p3mup2*F(5667)+p4mup2*F(5668)))
       Fa7(1849) = p3mup4*(p4mup2*F(1622)+2*(p1mup2*F(1325)+p2mup2*F(1
     -   623)))+p3mup2*(p4mup4*F(352)+2*(p1mup4*F(276)+p2mup4*F(359)+
     -   p3mup4*F(1843)))
       Fa7(1850) = p2mup4*p4mup2*F(1024)+p2mup2*p4mup4*F(1024)+4*(p3mu
     -   p4*p4mup2*F(1047)+p3mup2*p3mup4*F(1059))+2*(p1mup2*(2*(p1mup
     -   4*F(130)+p3mup4*F(278))+p4mup4*F(673)+p2mup4*F(962))+p2mup2*
     -   p2mup4*F(1022)+p3mup2*p4mup4*F(1026)+p4mup2*p4mup4*F(1027)+p
     -   1mup4*(p2mup2*F(964)+2*(p4mup2*F(338)+p3mup2*F(1028)))+p2mup
     -   2*p3mup4*F(1864)+p2mup4*p3mup2*F(1873))
       Fa7(1851) = p1mup2*(2*(p1mup4*F(668)+p2mup4*F(670))+p4mup4*F(10
     -   33))-p1mup4*p4mup2*F(1034)+4*p2mup2*p2mup4*F(1036)+2*(p1mup4
     -   *p2mup2*F(676)+(p2mup4*p4mup2+p2mup2*p4mup4)*F(1037)+p3mup4*
     -   (p4mup2*F(1046)+2*p2mup2*F(1060)+p1mup2*F(1864))+p3mup2*(p4m
     -   up4*F(1046)+2*(p3mup4*F(304)+p2mup4*F(1060))+p1mup4*F(1873))
     -   +2*p4mup2*p4mup4*F(7454))
       Fa7(1852) = -(p2mup4*p3mup2*F(7441))+p3mup4*(p2mup2*F(1626)+2*(
     -   p1mup2*F(1012)+p3mup2*F(7442)+p4mup2*F(7444)))
       Fa7(1853) = p1mup4*p3mup2*F(352)+p1mup2*p3mup4*F(352)+2*(2*p3mu
     -   p2*p3mup4*F(7442)+p3mup4*p4mup2*F(7444)+p3mup2*p4mup4*F(7444
     -   )+(p2mup4*p3mup2+p2mup2*p3mup4)*F(7448))
       Fa7(1854) = -(mup2mup5*(p4mup4*F(352)+2*(p1mup4*F(276)+p2mup4*F
     -   (359)+p3mup4*F(888))))+mup2mup4*(p4mup5*F(352)+2*(p1mup5*F(2
     -   76)+p2mup5*F(359)+p3mup5*F(888)))+mup4mup5*(p4mup2*F(1622)+2
     -   *(p1mup2*F(1325)+p2mup2*F(1623)+p3mup2*F(2717)))
       Fa7(1855) = mup4mup5*p3mup2+mup2mup5*p3mup4-mup2mup4*p3mup5
       Fa7(1856) = p3mup4*(p4mup2*F(1645)-2*(p2mup2*F(1335)+p1mup2*F(1
     -   646)))+p3mup2*(p4mup4*F(430)+2*(p1mup4*F(349)+p2mup4*F(1050)
     -   +p3mup4*F(5669)))
       Fa7(1857) = p2mup2*p4mup4*F(1658)+p2mup4*p4mup2*F(1659)+2*(p1mu
     -   p4*(2*(p3mup2*F(278)+p4mup2*F(333)+p1mup2*F(1627))+p2mup2*F(
     -   1631))+p2mup2*p2mup4*F(1633)+p3mup2*p4mup4*F(1660)+p4mup2*p4
     -   mup4*F(1661)+p1mup2*(p2mup4*F(1630)+p4mup4*F(1657)+2*p3mup4*
     -   F(5670))+p2mup2*p3mup4*F(5671)+p2mup4*p3mup2*F(5673))+4*(p3m
     -   up2*p3mup4*F(5672)+p3mup4*p4mup2*F(5674))
       Fa7(1858) = p1mup2*p4mup4*F(1647)+4*(p2mup2*p3mup4*F(5676)+p3mu
     -   p2*p3mup4*F(5677))+2*(p1mup2*p1mup4*F(1635)+p1mup4*p2mup2*F(
     -   1641)+p2mup4*(p1mup2*F(1636)+2*(p3mup2*F(890)+p2mup2*F(1642)
     -   ))+p2mup2*p4mup4*F(1649)+p1mup2*p3mup4*F(5675)+p1mup4*p3mup2
     -   *F(5678)+p3mup2*p4mup4*F(5679))+p4mup2*(4*p4mup4*F(306)+p1mu
     -   p4*F(1648)+2*(p2mup4*F(1650)+p3mup4*F(5680)))
       Fa7(1859) = -(p2mup4*p3mup2*F(1095))+p3mup4*(p2mup2*F(1653)+2*(
     -   p4mup2*F(436)+p3mup2*F(439)+p1mup2*F(1654)))
       Fa7(1860) = 2*(p2mup4*p3mup2*F(434)+p2mup2*p3mup4*F(434)+p3mup4
     -   *p4mup2*F(436)+p3mup2*p4mup4*F(436))+p1mup4*p3mup2*F(438)+p1
     -   mup2*p3mup4*F(438)+4*p3mup2*p3mup4*F(439)
       Fa7(1861) = mup2mup4*(p4mup5*F(429)-2*(p1mup5*F(276)+p2mup5*F(3
     -   59)+p3mup5*F(888)))+mup2mup5*(p4mup4*F(430)+2*(p1mup4*F(349)
     -   +p2mup4*F(1050)+p3mup4*F(5257)))+mup4mup5*(-(p4mup2*F(1645))
     -   +2*(p2mup2*F(1335)+p1mup2*F(1646)+p3mup2*F(5681)))
       Fa7(1862) = p3mup4*(2*(p1mup2*F(1325)+p2mup2*F(1623))+p4mup2*F(
     -   1645))+p3mup2*(p4mup4*F(430)+2*(p1mup4*F(349)+p2mup4*F(1050)
     -   +p3mup4*F(5669)))
       Fa7(1863) = p1mup2*p4mup4*F(1647)+4*(p2mup2*p3mup4*F(5676)+p3mu
     -   p2*p3mup4*F(5677))+2*(p1mup2*p1mup4*F(1635)+p1mup4*p2mup2*F(
     -   1641)+p2mup4*(p1mup2*F(1636)+2*(p3mup2*F(890)+p2mup2*F(1642)
     -   ))+p2mup2*p4mup4*F(1649)+p1mup2*p3mup4*F(5675)+p1mup4*p3mup2
     -   *F(5678)+p3mup2*p4mup4*F(5679))+p4mup2*(4*p4mup4*F(306)-p1mu
     -   p4*F(1663)-2*(p2mup4*F(1664)+p3mup4*F(5682)))
       Fa7(1864) = mup2mup4*(p4mup5*F(441)+2*(p1mup5*F(276)+p2mup5*F(3
     -   59)+p3mup5*F(888)))+mup4mup5*(p4mup2*F(1645)+2*(p1mup2*F(132
     -   5)+p2mup2*F(1623)+p3mup2*F(2717)))-mup2mup5*(p4mup4*F(430)+2
     -   *(p1mup4*F(349)+p2mup4*F(1050)+p3mup4*F(5257)))
       Fa7(1865) = p3mup4*(p4mup2*F(1622)-2*(p2mup2*F(1335)+p1mup2*F(1
     -   646)))+p3mup2*(p4mup4*F(352)+2*(p1mup4*F(276)+p2mup4*F(359)+
     -   p3mup4*F(1843)))
       Fa7(1866) = mup2mup5*(p4mup4*F(352)+2*(p1mup4*F(276)+p2mup4*F(3
     -   59)+p3mup4*F(888)))-mup2mup4*(p4mup5*F(352)+2*(p1mup5*F(276)
     -   +p2mup5*F(359)+p3mup5*F(888)))+mup4mup5*(-(p4mup2*F(1622))+2
     -   *(p2mup2*F(1335)+p1mup2*F(1646)+p3mup2*F(5681)))
       Fa7(1867) = 8*p2mup2*p2mup4*F(5561)+2*p3mup2*p4mup4*F(5566)+4*(
     -   p2mup4*p3mup2*F(5558)+p2mup4*p4mup2*F(5559)+p1mup4*(p3mup2*F
     -   (5555)+p4mup2*F(5557)+2*p2mup2*F(5562))+p2mup2*p4mup4*F(5564
     -   )-p4mup2*p4mup4*F(5567))+p1mup2*(-8*(p1mup4*F(5552)+p2mup4*F
     -   (5553))-4*p4mup4*F(5556)+p3mup4*F(5683))+p3mup4*p4mup2*F(568
     -   4)+p2mup2*p3mup4*F(5685)+p3mup2*p3mup4*F(5686)
       Fa7(1868) = p2mup4*p4mup2*F(5575)+4*p3mup2*p3mup4*F(5577)+p2mup
     -   2*p4mup4*F(5581)+2*(p2mup2*p2mup4*F(5574)+p3mup4*p4mup2*F(55
     -   76)+p1mup4*(p2mup2*F(5570)+p4mup2*F(5573)+2*p3mup2*F(5578))+
     -   p2mup2*p3mup4*F(5579)+p2mup4*p3mup2*F(5580)+p3mup2*p4mup4*F(
     -   5582)+p4mup2*p4mup4*F(5583)+p1mup2*(2*(p1mup4*F(5568)+p3mup4
     -   *F(5571))+p2mup4*F(5687)+p4mup4*F(5688)))
       Fa7(1869) = -8*(p2mup4*p3mup2*F(5592)+p3mup2*p3mup4*F(5594)+p4m
     -   up2*p4mup4*F(5595))-2*p1mup2*(p4mup4*F(5590)+2*(p1mup4*F(558
     -   4)+p2mup4*F(5585)+p3mup4*F(5689)))+p1mup4*p4mup2*F(5690)-4*(
     -   p1mup4*p3mup2*F(5589)+p3mup2*p4mup4*F(5596)+p2mup2*(p1mup4*F
     -   (5587)+2*(p2mup4*F(5586)+p3mup4*F(5598))+p4mup4*F(5599))-p2m
     -   up4*p4mup2*F(5691)-p3mup4*p4mup2*F(5692))
       Fa7(1870) = 2*p2mup4*p3mup2*F(5608)+p1mup4*(-4*(p2mup2*F(5602)-
     -   p3mup2*F(5604))+8*p4mup2*F(5609))+8*p4mup2*p4mup4*F(5614)+p1
     -   mup2*(-4*p2mup4*F(5601)-8*(p1mup4*F(5600)-p4mup4*F(5609))+p3
     -   mup4*F(5693))+p3mup2*p3mup4*F(5694)+p2mup2*p3mup4*F(5695)-4*
     -   (p2mup2*p2mup4*F(5605)-p2mup2*p4mup4*F(5610)-p2mup4*p4mup2*F
     -   (5611)-p3mup2*p4mup4*F(5696))+p3mup4*p4mup2*F(5697)
       Fa7(1871) = -2*p1mup4*p3mup2*F(5621)-8*(p3mup2*p3mup4*F(5620)-p
     -   2mup4*p4mup2*F(5626)-p4mup2*p4mup4*F(5628))-4*(p2mup4*p3mup2
     -   *F(5623)-p1mup4*p4mup2*F(5625)+p3mup2*p4mup4*F(5629))+p2mup2
     -   *(-8*(p2mup4*F(5618)+p4mup4*F(5630))+4*p1mup4*F(5699)+p3mup4
     -   *F(5700))+p1mup2*(p3mup4*F(5698)-4*(p1mup4*F(5615)+p2mup4*F(
     -   5616)-p4mup4*F(5701)))+p3mup4*p4mup2*F(5702)
       Fa7(1872) = 2*p2mup2*(4*p2mup4*F(5634)+p1mup4*F(5636)+2*(p3mup4
     -   *F(5638)+p4mup4*F(5644)))+4*(p4mup2*(p1mup4*F(5640)+p2mup4*F
     -   (5641)+2*(p3mup4*F(5642)+p4mup4*F(5643)))+p3mup2*(p1mup4*F(5
     -   637)+p2mup4*F(5638)+2*(p3mup4*F(5635)+p4mup4*F(5645))))+p1mu
     -   p2*(2*p2mup4*F(5632)+4*(p1mup4*F(5631)+p3mup4*F(5633))+p4mup
     -   4*F(5703))
       Fa7(1873) = mup2mup4*(4*(p1mup5*F(5454)+p2mup5*F(5455)+p3mup5*F
     -   (5648))+p4mup5*F(5704))-2*(mup2mup5*(2*(p1mup4*F(5646)+p2mup
     -   4*F(5647)+p3mup4*F(5650))+p4mup4*F(5651))-mup4mup5*(p4mup2*F
     -   (5653)+2*(p1mup2*F(5452)+p2mup2*F(5461)+p3mup2*F(5705))))
       Fa7(1874) = mup2mup5*(-4*(p1mup4*F(5116)-p4mup4*F(5122)+p2mup4*
     -   F(5658))+p3mup4*F(5706))+2*(mup4mup5*(-(p3mup2*F(5661))+2*(p
     -   1mup2*F(5453)-p2mup2*F(5462)+p4mup2*F(5662)))+mup2mup4*(2*(p
     -   1mup5*F(5447)+p4mup5*F(5448)+p2mup5*F(5456))+p3mup5*F(5707))
     -   )
       Fa7(1875) = -4*p3mup1*(p1mup5*F(5454)+p2mup5*F(5455)+p3mup5*F(5
     -   648))+2*p3mup1*p4mup5*F(5649)
       Fa7(1876) = 2*p3mup1*(p1mup5*F(5654)+2*(p2mup5*F(5655)+p3mup5*F
     -   (5656)+p4mup5*F(5657)))
       Fa7(1877) = -4*p1mup5*(p3mup1*F(4528)-p2mup1*F(5447))-2*(p3mup1
     -   *(2*(p3mup5*F(4529)+p4mup5*F(4668))+p2mup5*F(5665))-p2mup1*(
     -   2*(p4mup5*F(5448)+p2mup5*F(5456))+p3mup5*F(5707)))
       Fa7(1878) = p2mup1*(p1mup5*F(5654)+2*(p2mup5*F(5655)+p3mup5*F(5
     -   656)+p4mup5*F(5657)))
       Fa7(1879) = 4*(p1mup4*p3mup1*F(4663)-p3mup1*p4mup4*F(4672)+p3mu
     -   p1*p3mup4*F(4961)-p1mup4*p2mup1*F(5116)+p2mup1*p4mup4*F(5122
     -   )-p2mup1*p2mup4*F(5658))+2*p2mup4*p3mup1*F(5664)+p2mup1*p3mu
     -   p4*F(5706)
       Fa7(1880) = p3mup1*(p1mup4*F(5654)+2*(p2mup4*F(5655)+p3mup4*F(5
     -   656)+p4mup4*F(5657)))
       Fa7(1881) = 2*p2mup1*(p1mup4*F(5654)+2*(p2mup4*F(5655)+p3mup4*F
     -   (5656)+p4mup4*F(5657)))
       Fa7(1882) = 2*p3mup1*(2*(p1mup4*F(5646)+p2mup4*F(5647)+p3mup4*F
     -   (5650))+p4mup4*F(5651))
       Fa7(1883) = mup1mup2*(4*(p1mup5*F(5117)+p4mup5*F(5118)+p2mup5*F
     -   (5529))+p3mup5*F(5660))
       Fa7(1884) = 2*mup1mup2*(2*(p1mup5*F(4528)+p3mup5*F(4529)+p4mup5
     -   *F(4668))+p2mup5*F(5665))
       Fa7(1885) = p2mup1*(2*(p1mup4*F(5646)+p2mup4*F(5647)+p3mup4*F(5
     -   650))+p4mup4*F(5651))
       Fa7(1886) = p2mup1*(4*(p1mup5*F(5454)+p2mup5*F(5455)+p3mup5*F(5
     -   648))+p4mup5*F(5704))
       Fa7(1887) = -4*mup1mup2*(p1mup5*F(5454)+p2mup5*F(5455)+p3mup5*F
     -   (5648))+2*mup1mup2*p4mup5*F(5649)
       Fa7(1888) = mup1mup2*(p1mup5*F(5654)+2*(p2mup5*F(5655)+p3mup5*F
     -   (5656)+p4mup5*F(5657)))
       Fa7(1889) = p3mup1*(p4mup5*F(352)+2*(p1mup5*F(276)+p2mup5*F(359
     -   )+p3mup5*F(888)))
       Fa7(1890) = p3mup1*(p1mup5*F(352)+2*(p3mup5*F(7442)+p4mup5*F(74
     -   44)+p2mup5*F(7448)))
       Fa7(1891) = p2mup5*p3mup1-p2mup1*p3mup5
       Fa7(1892) = -(p2mup4*p3mup1)+p2mup1*p3mup4
       Fa7(1893) = p3mup1*(p1mup4*F(352)+2*(p3mup4*F(7442)+p4mup4*F(74
     -   44)+p2mup4*F(7448)))
       Fa7(1894) = p3mup1*(p4mup4*F(352)+2*(p1mup4*F(276)+p2mup4*F(359
     -   )+p3mup4*F(888)))
       Fa7(1895) = p2mup1*(p4mup4*F(352)+2*(p1mup4*F(276)+p2mup4*F(359
     -   )+p3mup4*F(888)))
       Fa7(1896) = p2mup1*(p4mup5*F(352)+2*(p1mup5*F(276)+p2mup5*F(359
     -   )+p3mup5*F(888)))
       Fa7(1897) = mup1mup2*(p4mup5*F(352)+2*(p1mup5*F(276)+p2mup5*F(3
     -   59)+p3mup5*F(888)))
       Fa7(1898) = p3mup1*(-(p4mup5*F(429))+2*(p1mup5*F(276)+p2mup5*F(
     -   359)+p3mup5*F(888)))
       Fa7(1899) = p3mup1*(p1mup5*F(438)+2*(p2mup5*F(434)+p4mup5*F(436
     -   )+p3mup5*F(439)))
       Fa7(1900) = -(p2mup5*p3mup1)+p2mup1*p3mup5
       Fa7(1901) = p2mup4*p3mup1-p2mup1*p3mup4
       Fa7(1902) = p3mup1*(p1mup4*F(438)+2*(p2mup4*F(434)+p4mup4*F(436
     -   )+p3mup4*F(439)))
       Fa7(1903) = p3mup1*(p4mup4*F(430)+2*(p1mup4*F(349)+p2mup4*F(105
     -   0)+p3mup4*F(5257)))
       Fa7(1904) = p2mup1*(p4mup4*F(430)+2*(p1mup4*F(349)+p2mup4*F(105
     -   0)+p3mup4*F(5257)))
       Fa7(1905) = p2mup1*(p4mup5*F(441)+2*(p1mup5*F(276)+p2mup5*F(359
     -   )+p3mup5*F(888)))
       Fa7(1906) = mup1mup2*(-(p4mup5*F(429))+2*(p1mup5*F(276)+p2mup5*
     -   F(359)+p3mup5*F(888)))
       Fa7(1907) = p3mup1*(p4mup5*F(441)+2*(p1mup5*F(276)+p2mup5*F(359
     -   )+p3mup5*F(888)))
       Fa7(1908) = p2mup1*(-(p4mup5*F(429))+2*(p1mup5*F(276)+p2mup5*F(
     -   359)+p3mup5*F(888)))
       Fa7(1909) = mup1mup2*(p4mup5*F(441)+2*(p1mup5*F(276)+p2mup5*F(3
     -   59)+p3mup5*F(888)))
       Fa7(1910) = p3mup1*(4*(p1mup5*F(5454)+p2mup5*F(5455)+p3mup5*F(5
     -   648))+p4mup5*F(5704))
       Fa7(1911) = p3mup1*(p1mup5*F(5654)+2*(p2mup5*F(5655)+p3mup5*F(5
     -   656)+p4mup5*F(5657)))
       Fa7(1912) = 4*p1mup5*(p3mup1*F(4528)+p2mup1*F(5117))+p2mup1*(4*
     -   (p4mup5*F(5118)+p2mup5*F(5529))+p3mup5*F(5660))+2*p3mup1*(2*
     -   (p3mup5*F(4529)+p4mup5*F(4668))+p2mup5*F(5665))
       Fa7(1913) = 2*p2mup1*(p1mup5*F(5654)+2*(p2mup5*F(5655)+p3mup5*F
     -   (5656)+p4mup5*F(5657)))
       Fa7(1914) = -4*(p1mup4*p3mup1*F(4663)-p3mup1*p4mup4*F(4672)+p3m
     -   up1*p3mup4*F(4961)-p1mup4*p2mup1*F(5116)+p2mup1*p4mup4*F(512
     -   2)-p2mup1*p2mup4*F(5658))+p2mup1*p3mup4*F(5659)-2*p2mup4*p3m
     -   up1*F(5664)
       Fa7(1915) = 2*p3mup1*(p1mup4*F(5654)+2*(p2mup4*F(5655)+p3mup4*F
     -   (5656)+p4mup4*F(5657)))
       Fa7(1916) = p2mup1*(p1mup4*F(5654)+2*(p2mup4*F(5655)+p3mup4*F(5
     -   656)+p4mup4*F(5657)))
       Fa7(1917) = p3mup1*(2*(p1mup4*F(5646)+p2mup4*F(5647)+p3mup4*F(5
     -   650))+p4mup4*F(5651))
       Fa7(1918) = 2*mup1mup2*(2*(p1mup5*F(5447)+p4mup5*F(5448)+p2mup5
     -   *F(5456))+p3mup5*F(5707))
       Fa7(1919) = mup1mup2*(2*(p1mup5*F(4528)+p3mup5*F(4529)+p4mup5*F
     -   (4668))+p2mup5*F(5665))
       Fa7(1920) = 2*p2mup1*(2*(p1mup4*F(5646)+p2mup4*F(5647)+p3mup4*F
     -   (5650))+p4mup4*F(5651))
       Fa7(1921) = -4*p2mup1*(p1mup5*F(5454)+p2mup5*F(5455)+p3mup5*F(5
     -   648))+2*p2mup1*p4mup5*F(5649)
       Fa7(1922) = mup1mup2*(4*(p1mup5*F(5454)+p2mup5*F(5455)+p3mup5*F
     -   (5648))+p4mup5*F(5704))
       Fa7(1923) = 2*mup1mup2*(p1mup5*F(5654)+2*(p2mup5*F(5655)+p3mup5
     -   *F(5656)+p4mup5*F(5657)))
       Fa7(1924) = mup1mup5*F(5048)+p4mup1*(4*(p1mup5*F(5117)+p4mup5*F
     -   (5118)+p2mup5*F(5529))+p3mup5*F(5660))
       Fa7(1925) = mup1mup5*F(5051)-2*p4mup1*(p1mup5*F(5654)+2*(p2mup5
     -   *F(5655)+p3mup5*F(5656)+p4mup5*F(5657)))
       Fa7(1926) = p1mup1*(4*(p1mup4*F(5116)-p4mup4*F(5122)+p2mup4*F(5
     -   658))+p3mup4*F(5659))
       Fa7(1927) = p1mup1*(2*(p1mup4*F(4663)-p4mup4*F(4672)+p3mup4*F(4
     -   961))+p2mup4*F(5664))
       Fa7(1928) = -2*(p2mup1*p3mup2*F(5661)-p2mup2*p3mup1*F(5666))-4*
     -   (p2mup1*p2mup2*F(5462)-p1mup2*(p2mup1*F(5453)+p3mup1*F(5663)
     -   )-p3mup1*p3mup2*F(5667))
       Fa7(1929) = p3mup1*(p1mup2*F(5654)+2*(p2mup2*F(5655)+p3mup2*F(5
     -   656)))+p1mup1*(-2*(p1mup2*F(5453)-p2mup2*F(5462))+p3mup2*F(5
     -   661))
       Fa7(1930) = p1mup2*p2mup1*F(5654)+p1mup1*p2mup2*F(5666)+2*(p2mu
     -   p1*p2mup2*F(5655)+p2mup1*p3mup2*F(5656)+p1mup1*p1mup2*F(5663
     -   )+p1mup1*p3mup2*F(5667))
       Fa7(1931) = p4mup1*(4*(p1mup4*F(5116)-p4mup4*F(5122)+p2mup4*F(5
     -   658))+p3mup4*F(5659))
       Fa7(1932) = mup1mup4*(4*(p1mup5*F(5454)+p2mup5*F(5455)+p3mup5*F
     -   (5648))+p4mup5*F(5704))
       Fa7(1933) = 2*mup1mup4*(2*(p1mup5*F(5447)+p4mup5*F(5448)+p2mup5
     -   *F(5456))+p3mup5*F(5707))
       Fa7(1934) = p3mup1*(p1mup2*F(5452)-p2mup2*F(5532)-p3mup2*F(5652
     -   ))
       Fa7(1935) = p2mup1*(p1mup2*F(5452)+p2mup2*F(5461)+p3mup2*F(5705
     -   ))
       Fa7(1936) = p4mup1*(2*(p1mup4*F(4663)-p4mup4*F(4672)+p3mup4*F(4
     -   961))+p2mup4*F(5664))
       Fa7(1937) = 2*mup1mup2*(2*(p1mup4*F(5646)+p2mup4*F(5647)+p3mup4
     -   *F(5650))+p4mup4*F(5651))+4*mup1mup4*(-(p1mup2*F(5452))+p2mu
     -   p2*F(5532)+p3mup2*F(5652))
       Fa7(1938) = mup1mup2*(4*(p1mup4*F(5116)-p4mup4*F(5122)+p2mup4*F
     -   (5658))+p3mup4*F(5659))+mup1mup4*(-4*(p1mup2*F(5453)-p2mup2*
     -   F(5462))+2*p3mup2*F(5661))-2*mup2mup4*p3mup1*F(5708)
       Fa7(1939) = mup1mup2*(2*(p1mup4*F(4663)-p4mup4*F(4672)+p3mup4*F
     -   (4961))+p2mup4*F(5664))+mup1mup4*(p2mup2*F(5666)+2*(p1mup2*F
     -   (5663)+p3mup2*F(5667)))-mup2mup4*p2mup1*F(5708)
       Fa7(1940) = p1mup1*(2*(p1mup4*F(5646)+p2mup4*F(5647)+p3mup4*F(5
     -   650))+p4mup4*F(5651))+p4mup1*(p1mup4*F(5654)+2*(p2mup4*F(565
     -   5)+p3mup4*F(5656)+p4mup4*F(5657)))
       Fa7(1941) = 2*mup1mup4*(p1mup5*F(5654)+2*(p2mup5*F(5655)+p3mup5
     -   *F(5656)+p4mup5*F(5657)))
       Fa7(1942) = mup1mup5*F(5049)+2*p4mup1*(2*(p1mup5*F(4528)+p3mup5
     -   *F(4529)+p4mup5*F(4668))+p2mup5*F(5665))
       Fa7(1943) = mup1mup4*(2*(p1mup5*F(4528)+p3mup5*F(4529)+p4mup5*F
     -   (4668))+p2mup5*F(5665))
       Fa7(1944) = p1mup1*(-(p1mup2*F(5452))+p2mup2*F(5532)+p3mup2*F(5
     -   652))
       Fa7(1945) = -2*(mup1mup4*(p1mup2*F(5654)+2*(p2mup2*F(5655)+p3mu
     -   p2*F(5656)))-mup1mup2*(p1mup4*F(5654)+2*(p2mup4*F(5655)+p3mu
     -   p4*F(5656)+p4mup4*F(5657)))+mup2mup4*p1mup1*F(5708))
       Fa7(1946) = mup1mup5*F(5126)+2*p3mup5*p4mup1*F(7441)
       Fa7(1947) = 2*p1mup5*p4mup1*F(352)+mup1mup5*F(5142)+4*p4mup1*(p
     -   3mup5*F(7442)+p4mup5*F(7444)+p2mup5*F(7448))
       Fa7(1948) = p2mup1*p3mup2*F(7441)-p3mup1*(p2mup2*F(1626)+2*(p1m
     -   up2*F(1012)+p3mup2*F(7442)))
       Fa7(1949) = p1mup2*p3mup1*F(352)+p1mup1*p3mup2*F(7441)+2*p3mup1
     -   *(p3mup2*F(7442)+p2mup2*F(7448))
       Fa7(1950) = mup1mup4*(p4mup5*F(352)+2*(p1mup5*F(276)+p2mup5*F(3
     -   59)+p3mup5*F(888)))
       Fa7(1951) = p3mup1*(p1mup2*F(1325)+p2mup2*F(1623)+p3mup2*F(2717
     -   ))
       Fa7(1952) = p2mup1*(p2mup2*F(1335)+p1mup2*F(1646)+p3mup2*F(5681
     -   ))
       Fa7(1953) = mup1mup2*(p4mup4*F(352)+2*(p1mup4*F(276)+p2mup4*F(3
     -   59)+p3mup4*F(888)))-2*mup1mup4*(p1mup2*F(1325)+p2mup2*F(1623
     -   )+p3mup2*F(2717))
       Fa7(1954) = -(mup2mup4*p3mup1)+mup1mup4*p3mup2+mup1mup2*p3mup4
       Fa7(1955) = p1mup1*(p4mup4*F(352)+2*(p1mup4*F(276)+p2mup4*F(359
     -   )+p3mup4*F(888)))+p4mup1*(p1mup4*F(352)+2*(p3mup4*F(7442)+p4
     -   mup4*F(7444)+p2mup4*F(7448)))
       Fa7(1956) = mup1mup5*F(5141)-2*p2mup5*p4mup1*F(7441)
       Fa7(1957) = p1mup1*(p1mup2*F(1325)+p2mup2*F(1623)+p3mup2*F(2717
     -   ))
       Fa7(1958) = mup1mup5*F(5178)+2*p3mup5*p4mup1*F(7441)
       Fa7(1959) = 2*p1mup5*p4mup1*F(352)+mup1mup5*F(5193)+4*p4mup1*(p
     -   3mup5*F(7442)+p4mup5*F(7444)+p2mup5*F(7448))
       Fa7(1960) = mup1mup5*F(5192)-2*p2mup5*p4mup1*F(7441)
       Fa7(1961) = -2*p3mup5*p4mup1*F(1095)+mup1mup5*F(5229)
       Fa7(1962) = 2*p4mup1*(p1mup5*F(438)+2*(p2mup5*F(434)+p4mup5*F(4
     -   36)+p3mup5*F(439)))-mup1mup5*F(5230)
       Fa7(1963) = -(p2mup1*p3mup2*F(1095))+p3mup1*(p2mup2*F(1653)+2*(
     -   p3mup2*F(439)+p1mup2*F(1654)))
       Fa7(1964) = p1mup2*p3mup1*F(438)+2*(p2mup2*p3mup1*F(434)+p3mup1
     -   *p3mup2*F(439))+p1mup1*p3mup2*F(1095)
       Fa7(1965) = mup1mup4*(p4mup5*F(441)+2*(p1mup5*F(276)+p2mup5*F(3
     -   59)+p3mup5*F(888)))
       Fa7(1966) = p3mup1*(p2mup2*F(1335)+p1mup2*F(1646)+p3mup2*F(5681
     -   ))
       Fa7(1967) = p2mup1*(p1mup2*F(1325)+p2mup2*F(1623)+p3mup2*F(2717
     -   ))
       Fa7(1968) = mup1mup2*(p4mup4*F(430)+2*(p1mup4*F(349)+p2mup4*F(1
     -   050)+p3mup4*F(5257)))+2*mup1mup4*(p2mup2*F(1335)+p1mup2*F(16
     -   46)+p3mup2*F(5681))
       Fa7(1969) = p4mup1*(p1mup4*F(438)+2*(p2mup4*F(434)+p4mup4*F(436
     -   )+p3mup4*F(439)))+p1mup1*(p4mup4*F(430)+2*(p1mup4*F(349)+p2m
     -   up4*F(1050)+p3mup4*F(5257)))
       Fa7(1970) = 2*p2mup5*p4mup1*F(1095)+mup1mup5*F(5246)
       Fa7(1971) = p1mup1*(p2mup2*F(1335)+p1mup2*F(1646)+p3mup2*F(5681
     -   ))
       Fa7(1972) = 2*p3mup5*p4mup1*F(1095)+mup1mup5*F(5226)
       Fa7(1973) = 2*p4mup1*(p1mup5*F(438)+2*(p2mup5*F(434)+p4mup5*F(4
     -   36)+p3mup5*F(439)))+mup1mup5*F(5247)
       Fa7(1974) = mup1mup4*(-(p4mup5*F(429))+2*(p1mup5*F(276)+p2mup5*
     -   F(359)+p3mup5*F(888)))
       Fa7(1975) = -2*mup1mup4*(p1mup2*F(1325)+p2mup2*F(1623)+p3mup2*F
     -   (2717))+mup1mup2*(p4mup4*F(430)+2*(p1mup4*F(349)+p2mup4*F(10
     -   50)+p3mup4*F(5257)))
       Fa7(1976) = -2*p2mup5*p4mup1*F(1095)+mup1mup5*F(5248)
       Fa7(1977) = mup1mup5*F(5175)-2*p3mup5*p4mup1*F(7441)
       Fa7(1978) = 2*p1mup5*p4mup1*F(352)-mup1mup5*F(5195)+4*p4mup1*(p
     -   3mup5*F(7442)+p4mup5*F(7444)+p2mup5*F(7448))
       Fa7(1979) = p2mup2*p3mup1*F(1626)-p2mup1*p3mup2*F(7441)+2*(p1mu
     -   p2*p3mup1*F(1012)+p3mup1*p3mup2*F(7442))
       Fa7(1980) = mup1mup2*(p4mup4*F(352)+2*(p1mup4*F(276)+p2mup4*F(3
     -   59)+p3mup4*F(888)))+2*mup1mup4*(p2mup2*F(1335)+p1mup2*F(1646
     -   )+p3mup2*F(5681))
       Fa7(1981) = mup1mup5*F(5194)+2*p2mup5*p4mup1*F(7441)
       Fa7(1982) = mup1mup5*F(5123)-2*p3mup5*p4mup1*F(7441)
       Fa7(1983) = 2*p1mup5*p4mup1*F(352)-mup1mup5*F(5144)+4*p4mup1*(p
     -   3mup5*F(7442)+p4mup5*F(7444)+p2mup5*F(7448))
       Fa7(1984) = mup1mup5*F(5143)+2*p2mup5*p4mup1*F(7441)
       Fa7(1985) = mup1mup5*F(4978)+2*p4mup1*(2*(p1mup5*F(5447)+p4mup5
     -   *F(5448)+p2mup5*F(5456))+p3mup5*F(5707))
       Fa7(1986) = mup1mup5*F(5054)+2*p4mup1*(p1mup5*F(5654)+2*(p2mup5
     -   *F(5655)+p3mup5*F(5656)+p4mup5*F(5657)))
       Fa7(1987) = p1mup1*(-4*(p1mup4*F(5116)-p4mup4*F(5122)+p2mup4*F(
     -   5658))+p3mup4*F(5706))
       Fa7(1988) = 2*p1mup1*(2*(p1mup4*F(4663)-p4mup4*F(4672)+p3mup4*F
     -   (4961))+p2mup4*F(5664))
       Fa7(1989) = -(p2mup1*p3mup2*F(5661))+p2mup2*p3mup1*F(5666)-2*(p
     -   2mup1*p2mup2*F(5462)-p1mup2*(p2mup1*F(5453)+p3mup1*F(5663))-
     -   p3mup1*p3mup2*F(5667))
       Fa7(1990) = p4mup1*(-4*(p1mup4*F(5116)-p4mup4*F(5122)+p2mup4*F(
     -   5658))+p3mup4*F(5706))
       Fa7(1991) = -4*mup1mup4*(p1mup5*F(5454)+p2mup5*F(5455)+p3mup5*F
     -   (5648))+2*mup1mup4*p4mup5*F(5649)
       Fa7(1992) = mup1mup4*(4*(p1mup5*F(5117)+p4mup5*F(5118)+p2mup5*F
     -   (5529))+p3mup5*F(5660))
       Fa7(1993) = p3mup1*(p1mup2*F(5452)+p2mup2*F(5461)+p3mup2*F(5705
     -   ))
       Fa7(1994) = p2mup1*(-(p1mup2*F(5452))+p2mup2*F(5532)+p3mup2*F(5
     -   652))
       Fa7(1995) = 2*p4mup1*(2*(p1mup4*F(4663)-p4mup4*F(4672)+p3mup4*F
     -   (4961))+p2mup4*F(5664))
       Fa7(1996) = -2*mup1mup2*(2*(p1mup4*F(5646)+p2mup4*F(5647)+p3mup
     -   4*F(5650))+p4mup4*F(5651))+4*mup1mup4*(p1mup2*F(5452)+p2mup2
     -   *F(5461)+p3mup2*F(5705))
       Fa7(1997) = mup1mup4*(4*(p1mup2*F(5453)-p2mup2*F(5462))-2*p3mup
     -   2*F(5661))+mup1mup2*(-4*(p1mup4*F(5116)-p4mup4*F(5122)+p2mup
     -   4*F(5658))+p3mup4*F(5706))+2*mup2mup4*p3mup1*F(5708)
       Fa7(1998) = mup1mup4*(p1mup5*F(5654)+2*(p2mup5*F(5655)+p3mup5*F
     -   (5656)+p4mup5*F(5657)))
       Fa7(1999) = mup1mup5*F(5053)-2*p4mup1*(2*(p1mup5*F(4528)+p3mup5
     -   *F(4529)+p4mup5*F(4668))+p2mup5*F(5665))
       Fa7(2000) = 2*mup1mup4*(2*(p1mup5*F(4528)+p3mup5*F(4529)+p4mup5
     -   *F(4668))+p2mup5*F(5665))
       Fa7(2001) = p1mup1*(p1mup2*F(5452)+p2mup2*F(5461)+p3mup2*F(5705
     -   ))
       Fa7(2002) = 2*(mup1mup4*(p1mup2*F(5654)+2*(p2mup2*F(5655)+p3mup
     -   2*F(5656)))-mup1mup2*(p1mup4*F(5654)+2*(p2mup4*F(5655)+p3mup
     -   4*F(5656)+p4mup4*F(5657)))+mup2mup4*p1mup1*F(5708))
       Fa7(2003) = -4*(p1mup2*p4mup1*F(5453)-p2mup2*p4mup1*F(5462))-2*
     -   (p3mup1*p4mup2*F(5653)-p3mup2*p4mup1*F(5661))
       Fa7(2004) = p4mup1*(p2mup2*F(5666)+2*(p1mup2*F(5663)+p3mup2*F(5
     -   667)))
       Fa7(2005) = p4mup1*(p1mup2*F(5654)+2*(p2mup2*F(5655)+p3mup2*F(5
     -   656)))
       Fa7(2006) = p3mup1*p4mup2*F(1622)-p3mup2*p4mup1*F(7441)
       Fa7(2007) = p4mup1*(p2mup2*F(1626)+2*(p1mup2*F(1012)+p3mup2*F(7
     -   442)))
       Fa7(2008) = p4mup1*(p1mup2*F(352)+2*(p3mup2*F(7442)+p2mup2*F(74
     -   48)))
       Fa7(2009) = p3mup2*p4mup1*F(1095)-p3mup1*p4mup2*F(1645)
       Fa7(2010) = p4mup1*(p2mup2*F(1653)+2*(p3mup2*F(439)+p1mup2*F(16
     -   54)))
       Fa7(2011) = p4mup1*(p1mup2*F(438)+2*(p2mup2*F(434)+p3mup2*F(439
     -   )))
       Fa7(2012) = -(p3mup2*p4mup1*F(1095))+p3mup1*p4mup2*F(1645)
       Fa7(2013) = -(p3mup1*p4mup2*F(1622))+p3mup2*p4mup1*F(7441)
       Fa7(2014) = 4*(p1mup2*p4mup1*F(5453)-p2mup2*p4mup1*F(5462))+2*(
     -   p3mup1*p4mup2*F(5653)-p3mup2*p4mup1*F(5661))
       Fa7(2015) = 2*p4mup1*(p2mup2*F(5666)+2*(p1mup2*F(5663)+p3mup2*F
     -   (5667)))
       Fa7(2016) = 2*p4mup1*(p1mup2*F(5654)+2*(p2mup2*F(5655)+p3mup2*F
     -   (5656)))
       Fa7(2017) = p4mup2*(p2mup1*F(5657)+p1mup1*F(5668))
       Fa7(2018) = 2*p4mup2*(p2mup1*F(5653)-2*p4mup1*F(5668))
       Fa7(2019) = mup1mup4*p4mup2*F(5653)-mup2mup4*p4mup1*F(5708)
       Fa7(2020) = p4mup2*(p1mup1*F(5653)+2*p4mup1*F(5657))
       Fa7(2021) = -2*p2mup1*p4mup2*F(5653)+4*p4mup1*p4mup2*F(5668)
       Fa7(2022) = 2*p4mup2*(p1mup1*F(5653)+2*p4mup1*F(5657))
       Fa7(2023) = p1mup1*(p2mup5*F(5455)+p3mup5*F(5648))
       Fa7(2024) = p1mup1*(4*p2mup5*F(5529)+p3mup5*F(5660))
       Fa7(2025) = 2*p1mup1*(2*p3mup5*F(4529)+p2mup5*F(5665))
       Fa7(2026) = mup1mup5*(2*(p1mup4*F(5646)+p2mup4*F(5647)+p3mup4*F
     -   (5650))+p4mup4*F(5651))
       Fa7(2027) = mup1mup5*(-4*(p1mup4*F(5116)-p4mup4*F(5122)+p2mup4*
     -   F(5658))+p3mup4*F(5706))
       Fa7(2028) = 2*mup1mup5*(2*(p1mup2*F(5453)-p2mup2*F(5462))-p3mup
     -   2*F(5661))
       Fa7(2029) = 2*(mup1mup5*(p2mup2*F(5666)+2*(p1mup2*F(5663)+p3mup
     -   2*F(5667)))-mup2mup5*p2mup1*F(5708))
       Fa7(2030) = mup1mup5*(p1mup4*F(5654)+2*(p2mup4*F(5655)+p3mup4*F
     -   (5656)+p4mup4*F(5657)))
       Fa7(2031) = 2*mup1mup5*(2*(p1mup4*F(4663)-p4mup4*F(4672)+p3mup4
     -   *F(4961))+p2mup4*F(5664))
       Fa7(2032) = mup1mup5*(p1mup2*F(5452)+p2mup2*F(5461)+p3mup2*F(57
     -   05))
       Fa7(2033) = mup1mup5*(p1mup2*F(5654)+2*(p2mup2*F(5655)+p3mup2*F
     -   (5656)))+mup2mup5*p1mup1*F(5708)
       Fa7(2034) = p1mup1*(p2mup5*F(359)+p3mup5*F(888))
       Fa7(2035) = mup1mup5*(p4mup4*F(352)+2*(p1mup4*F(276)+p2mup4*F(3
     -   59)+p3mup4*F(888)))
       Fa7(2036) = mup1mup5*(p2mup2*F(1335)+p1mup2*F(1646)+p3mup2*F(56
     -   81))
       Fa7(2037) = mup1mup5*(p4mup4*F(430)+2*(p1mup4*F(349)+p2mup4*F(1
     -   050)+p3mup4*F(5257)))
       Fa7(2038) = mup1mup5*(p1mup2*F(1325)+p2mup2*F(1623)+p3mup2*F(27
     -   17))
       Fa7(2039) = 2*p1mup1*(2*p2mup5*F(5456)+p3mup5*F(5707))
       Fa7(2040) = p1mup1*(2*p3mup5*F(4529)+p2mup5*F(5665))
       Fa7(2041) = 2*mup1mup5*(2*(p1mup4*F(5646)+p2mup4*F(5647)+p3mup4
     -   *F(5650))+p4mup4*F(5651))
       Fa7(2042) = mup1mup5*(4*(p1mup4*F(5116)-p4mup4*F(5122)+p2mup4*F
     -   (5658))+p3mup4*F(5659))
       Fa7(2043) = 2*mup1mup5*(-2*(p1mup2*F(5453)-p2mup2*F(5462))+p3mu
     -   p2*F(5661))
       Fa7(2044) = mup1mup5*(p2mup2*F(5666)+2*(p1mup2*F(5663)+p3mup2*F
     -   (5667)))-mup2mup5*p2mup1*F(5708)
       Return
       End
