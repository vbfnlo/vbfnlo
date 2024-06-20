c************** Calling the Fa functions*************************
       subroutine test_FaFunctionhexlineAbe36(p1mup2,p1mup3,p1mup4,p1m
     -   up5,p2mup2,p2mup3,p2mup4,p2mup5,p3mup2,p3mup3,p3mup4,p3mup5,
     -   p4mup2,p4mup3,p4mup4,p4mup5,p5mup2,p5mup3,p5mup4,p5mup5,p6mu
     -   p2,p6mup3,p6mup4,p6mup5,mup2mup3,mup2mup4,mup2mup5,mup3mup4,
     -   mup3mup5,mup4mup5,Fa36)
       IMPLICIT NONE
      Complex*16   p1mup2, p1mup3, p1mup4, p1mup5, p2mup2, p2mup3, 
     -          p2mup4, p2mup5, p3mup2, p3mup3, p3mup4, p3mup5, 
     -          p4mup2, p4mup3, p4mup4, p4mup5, p5mup2, p5mup3, 
     -          p5mup4, p5mup5, p6mup2, p6mup3, p6mup4, p6mup5
       Complex*16   mup2mup3, mup2mup4, mup2mup5, mup3mup4, mup3mup5, 
     -          mup4mup5
        common/test_FhlFunctionsAbe/F
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
      Complex*16 Fa(445),F(5930)
      Real*8 P(11674) 
        COMMON/test_PFunctionsAbe/P
       Complex*16 Fa36(386:396)
       COMMON/test_FaFunctionsAbe/Fa
       Fa36(386) = p1mup2*(p1mup5*p3mup3*F(4450)-p3mup3*p6mup5*F(4460)
     -   +p3mup3*p3mup5*F(4487)-p3mup3*p4mup5*F(4820)-p1mup5*p2mup3*F
     -   (5242)-p2mup3*p6mup5*F(5245)+p1mup3*(p1mup5*F(4445)+p3mup5*F
     -   (4451)-p6mup5*F(4455)+p2mup5*F(5241)+p4mup5*F(5249))-p2mup3*
     -   p2mup5*F(5251)+p2mup3*p3mup5*F(5252)-p2mup5*p3mup3*F(5253)+p
     -   2mup3*p4mup5*F(5255))+p2mup2*(p1mup3*(p2mup5*F(5243)+p1mup5*
     -   F(5246)+p3mup5*F(5247)-p6mup5*F(5248)+p4mup5*F(5250))-p2mup3
     -   *(p2mup5*F(4331)+p1mup5*F(5244)-p3mup5*F(5254))+p1mup5*p3mup
     -   3*F(5256))
       Fa36(387) = p2mup2*(p2mup3*(-(p4mup5*F(5218))+p6mup5*F(5257))+p
     -   3mup3*(p4mup5*F(4821)+p6mup5*F(5258)+p2mup5*F(5259)-p3mup5*F
     -   (5260)))
       Fa36(388) = -(mup3mup5*p1mup2*F(5262))+mup3mup5*p2mup2*F(5263)+
     -   mup2mup5*p2mup3*F(5264)-mup3mup5*p6mup2*F(5265)-mup2mup3*p6m
     -   up5*F(5266)-mup2mup5*p4mup3*F(5270)+mup2mup5*p6mup3*F(5271)+
     -   p3mup3*(-(mup2mup5*F(5269))+2*p6mup5*(p6mup2*F(4529)-p4mup2*
     -   F(5273)+p3mup2*F(5275)))+p1mup3*(-(mup2mup5*F(5261))+2*p6mup
     -   5*(p6mup2*F(4538)-p4mup2*F(5274)+p3mup2*F(5279)))+2*(p1mup2*
     -   p6mup3*p6mup5*F(4540)-p6mup2*p6mup3*p6mup5*F(4553)+p4mup2*p6
     -   mup3*p6mup5*F(4835)-p2mup2*p4mup3*p6mup5*F(5267)-p1mup2*p4mu
     -   p3*p6mup5*F(5268)-p2mup3*p4mup2*p6mup5*F(5272)-p3mup2*p4mup3
     -   *p6mup5*F(5273)+p2mup2*p6mup3*p6mup5*F(5276)+p2mup3*p3mup2*p
     -   6mup5*F(5277)+p2mup3*p6mup2*p6mup5*F(5278)+p4mup2*p4mup3*p6m
     -   up5*F(5280)+p4mup3*p6mup2*p6mup5*F(5281)+p3mup2*p6mup3*p6mup
     -   5*F(5282))
       Fa36(389) = mup3mup5*(p3mup2*F(5283)+p4mup2*F(5284))+2*p6mup3*(
     -   -(p3mup5*p6mup2*F(4551))+p2mup2*p4mup5*F(5231)+p2mup2*p2mup5
     -   *F(5286)+p2mup2*p3mup5*F(5288)+p1mup2*(p3mup5*F(4527)-p1mup5
     -   *F(4539)+p2mup5*F(5285)+p4mup5*F(5289))+p1mup5*(-(p6mup2*F(4
     -   546))-p2mup2*F(5287)+p3mup2*F(5290)+p4mup2*F(5291))+p2mup5*p
     -   3mup2*F(5292)+p2mup5*p4mup2*F(5293)-p2mup5*p6mup2*F(5294)+p3
     -   mup2*p3mup5*F(5295)+p3mup5*p4mup2*F(5296)+p3mup2*p4mup5*F(52
     -   96)+p4mup2*p4mup5*F(5297)-p4mup5*p6mup2*F(5298))
       Fa36(390) = mup2mup3*(-(p1mup5*F(5299))+p2mup5*F(5300)+p3mup5*F
     -   (5301)+p4mup5*F(5302))+2*(p1mup2*p4mup3*(p3mup5*F(4820)-p4mu
     -   p5*F(4828)+p2mup5*F(5303)-p1mup5*F(5304))+p2mup2*p4mup3*(p2m
     -   up5*F(4705)+p3mup5*F(4821)-p4mup5*F(4845)-p1mup5*F(5305)))
       Fa36(391) = p1mup3*p3mup5*p6mup2*F(4526)+p1mup3*p1mup5*p6mup2*F
     -   (4537)-p1mup3*p1mup5*p4mup2*F(4848)+p3mup2*p3mup5*p4mup3*F(4
     -   860)+p3mup5*p4mup2*p4mup3*F(4862)+p3mup2*p4mup3*p4mup5*F(486
     -   2)-p3mup5*p4mup3*p6mup2*F(4886)-p4mup2*p4mup3*p4mup5*F(4888)
     -   -p4mup3*p4mup5*p6mup2*F(4889)+p2mup5*p3mup2*p4mup3*F(5236)-p
     -   1mup3*p2mup5*p3mup2*F(5306)-p1mup3*p2mup5*p4mup2*F(5308)+p1m
     -   up3*p2mup5*p6mup2*F(5310)-p1mup3*p3mup5*p4mup2*F(5312)-p1mup
     -   3*p3mup2*p4mup5*F(5312)+p1mup3*p1mup5*p3mup2*F(5313)+p1mup3*
     -   p3mup2*p3mup5*F(5314)+p1mup3*p4mup2*p4mup5*F(5319)+p1mup3*p4
     -   mup5*p6mup2*F(5320)+p2mup3*(-(p3mup2*p4mup5*F(4519))+p1mup5*
     -   p3mup2*F(5307)+p1mup5*p4mup2*F(5309)-p1mup5*p6mup2*F(5311)+p
     -   2mup5*(-(p4mup2*F(4386))+p3mup2*F(5315)+p6mup2*F(5316))+p3mu
     -   p5*(-(p4mup2*F(4519))+p3mup2*F(5317)+p6mup2*F(5318))+p4mup2*
     -   p4mup5*F(5321)+p4mup5*p6mup2*F(5322))+p1mup5*p3mup2*p4mup3*F
     -   (5324)+p1mup5*p4mup2*p4mup3*F(5325)-p1mup5*p4mup3*p6mup2*F(5
     -   326)+p2mup5*p4mup2*p4mup3*F(5329)-p2mup5*p4mup3*p6mup2*F(533
     -   0)+p3mup3*(p3mup2*p4mup5*F(4860)+p4mup2*p4mup5*F(4862)+p2mup
     -   5*p4mup2*F(5236)+p1mup5*(-(p6mup2*F(4545))+p3mup2*F(5323)+p4
     -   mup2*F(5324))+p2mup5*p3mup2*F(5327)-p6mup2*(p4mup5*F(4886)+p
     -   2mup5*F(5328))+p3mup5*(-(p6mup2*F(4495))+p4mup2*F(4860)+p3mu
     -   p2*F(5331)))
       Fa36(392) = -(mup2mup5*F(5332))+2*(p1mup2*p1mup5*F(5333)-p1mup2
     -   *p6mup5*F(5335)+p1mup2*p3mup5*F(5338)+p1mup2*p4mup5*F(5339)+
     -   p3mup2*(-3*(p4mup5*F(4903)-p3mup5*F(5133))-p6mup5*F(5337)+p1
     -   mup5*F(5341))+p4mup2*(3*(p4mup5*F(4635)-p3mup5*F(4903))-p6mu
     -   p5*F(5340)+p1mup5*F(5342))-p2mup2*(3*(p4mup5*F(4612)+p3mup5*
     -   F(5109))+p6mup5*F(5336)-p1mup5*F(5343))-p1mup5*p6mup2*F(5344
     -   )+p2mup5*(-3*(p4mup2*F(4612)-p2mup2*F(4623)+p3mup2*F(5109))+
     -   p1mup2*F(5334)+p6mup2*F(5345))+p6mup2*p6mup5*F(5346)+p3mup5*
     -   p6mup2*F(5347)+p4mup5*p6mup2*F(5348))
       Fa36(393) = p1mup5*F(5349)+p2mup5*F(5350)+p6mup5*F(5351)+p3mup5
     -   *F(5352)+p4mup5*F(5353)
       Fa36(394) = mup3mup5*F(5354)-2*(-3*(p3mup3*p4mup5*F(4626)-p4mup
     -   3*p4mup5*F(4637)+p2mup3*p4mup5*F(4903))-p1mup5*p2mup3*F(5356
     -   )+p1mup3*p1mup5*F(5357)-p1mup3*p6mup5*F(5359)-p2mup3*p6mup5*
     -   F(5360)-p3mup3*p6mup5*F(5361)+p1mup3*p4mup5*F(5362)-p4mup3*p
     -   6mup5*F(5363)+p2mup5*(-3*(p4mup3*F(4607)+p2mup3*F(5109)+p3mu
     -   p3*F(5112))+p1mup3*F(5355)+p6mup3*F(5364))+p1mup5*p4mup3*F(5
     -   365)+p1mup5*p3mup3*F(5366)+p6mup3*p6mup5*F(5367)-p1mup5*p6mu
     -   p3*F(5368)+p3mup5*(-3*(p4mup3*F(4626)-p2mup3*F(5133)-p3mup3*
     -   F(5134))+p1mup3*F(5358)+p6mup3*F(5369))+p4mup5*p6mup3*F(5370
     -   ))
       Fa36(395) = mup4mup5*p1mup3*F(5371)+mup3mup4*p1mup5*F(5372)+mup
     -   3mup5*p1mup4*F(5373)-mup3mup4*p2mup5*F(5375)-mup3mup4*p3mup5
     -   *F(5377)+mup3mup5*p3mup4*F(5378)-mup3mup4*p4mup5*F(5379)+mup
     -   3mup5*p6mup4*F(5380)+mup3mup5*p4mup4*F(5381)-mup3mup4*p6mup5
     -   *F(5382)+mup4mup5*p3mup3*F(5383)+mup4mup5*p4mup3*F(5384)+mup
     -   4mup5*p6mup3*F(5385)+p2mup4*(mup3mup5*F(5376)-2*(p1mup5*p4mu
     -   p3*F(4633)+p2mup5*p4mup3*F(4635)+p4mup3*p4mup5*F(4637)+p4mup
     -   5*p6mup3*F(4638)+p1mup3*p2mup5*F(4641)+p2mup5*p6mup3*F(4642)
     -   +p1mup3*p1mup5*F(4927)-p1mup5*p6mup3*F(4933)+p3mup3*(-(p4mup
     -   5*F(4626))-p2mup5*F(4903)+p1mup5*F(4930)+p6mup5*F(4934))-p6m
     -   up3*p6mup5*F(4935)+p3mup5*(-(p4mup3*F(4626))+p1mup3*F(4928)+
     -   p6mup3*F(4932)+p3mup3*F(5134))-p1mup3*p4mup5*F(5389)+p1mup3*
     -   p6mup5*F(5397)+p4mup3*p6mup5*F(5415)))-2*(p1mup3*p2mup5*p4mu
     -   p4*F(5386)-p1mup4*p2mup5*p3mup3*F(5387)-p2mup5*p3mup3*p3mup4
     -   *F(5388)+p2mup5*p3mup4*p4mup3*F(5392)-p1mup3*p1mup4*p2mup5*F
     -   (5393)-p1mup3*p2mup5*p3mup4*F(5394)-p1mup3*p2mup5*p6mup4*F(5
     -   400)-p2mup5*p3mup4*p6mup3*F(5403)+p1mup4*p2mup5*p4mup3*F(540
     -   8)+p1mup4*p2mup5*p6mup3*F(5409)+p2mup5*p3mup3*p4mup4*F(5412)
     -   -p2mup5*p3mup3*p6mup4*F(5413)+p2mup5*p4mup3*p4mup4*F(5414)-p
     -   2mup5*p4mup4*p6mup3*F(5416)-p2mup5*p4mup3*p6mup4*F(5416)-p2m
     -   up5*p6mup3*p6mup4*F(5417))+p2mup3*(-(mup4mup5*F(5374))+2*(-(
     -   p2mup4*(p3mup5*F(5390)-p4mup5*F(5391)+p1mup5*F(5396)+p6mup5*
     -   F(5398)))+p3mup4*(-(p3mup5*F(4626))+p1mup5*F(4633)+p4mup5*F(
     -   4637)+p2mup5*F(5133)-p6mup5*F(5399))+p1mup4*(-(p3mup5*F(5389
     -   ))+p2mup5*F(5395)+p1mup5*F(5404)+p4mup5*F(5406)+p6mup5*F(540
     -   7))+p4mup4*(p4mup5*F(4636)-p2mup5*F(4903)-p3mup5*F(4904)+p1m
     -   up5*F(5405)+p6mup5*F(5418))+p6mup4*(p2mup5*F(5401)+p3mup5*F(
     -   5402)+p4mup5*F(5410)-p1mup5*F(5411)+p6mup5*F(5419))))
       Fa36(396) = p1mup2*(p1mup3*(p1mup4*F(5420)+p2mup4*F(5421)+p3mup
     -   4*F(5422)+p4mup4*F(5423))+p1mup4*p2mup3*F(5424)+p2mup3*p2mup
     -   4*F(5425)+p2mup3*p3mup4*F(5426)+p2mup3*p4mup4*F(5427)-p1mup4
     -   *p3mup3*F(5428)-p2mup4*p3mup3*F(5429)-p3mup3*p3mup4*F(5430)-
     -   p3mup3*p4mup4*F(5431))+p2mup2*(p1mup3*(p1mup4*F(5432)+p2mup4
     -   *F(5433)+p3mup4*F(5434)+p4mup4*F(5435))+p1mup4*p2mup3*F(5436
     -   )+p2mup3*p2mup4*F(5437)+p2mup3*p3mup4*F(5438)+p2mup3*p4mup4*
     -   F(5439)-p1mup4*p3mup3*F(5440)-p2mup4*p3mup3*F(5441)-p3mup3*p
     -   3mup4*F(5442)-p3mup3*p4mup4*F(5443))
       Return
       End
