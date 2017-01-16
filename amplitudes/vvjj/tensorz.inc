c Dieters version from Apr. 12, 2005, extended for real emission contributions
c
c  precalculated lepton tensors for ZZ --> l lbar l' lbar' 
c	always write ze (zu) for the Z which goes to e+e- (mu+mu-)
c
c  lepton spinors and Z polarization vectors

      double complex lep(6),lem(6),lup(6),lum(6)
      common /lspinorZ/ lep,lem,lup,lum
      
      integer ie,iu
      
      double complex zep,zem,zup,zum,ze,zu,ae,au,xxea,xxua,
     1		     xxez,xxuz,zet,zut,aet,aut,ea
      double precision qe, qu,qa

      common /zemspinor/ ze(6), zu(6), ae(6), au(6), xxea, xxua ,xxez, 
     2  		 xxuz ,zet(6),zut(6),aet(6),aut(6), qe(0:4), qu(0:4)
      common /vaMom/ ea(6),qa(0:4)


c  leptonic tensors
      double complex aazz,azzz,zazz,zzzz,wwzz5,wwzz6, 
     3           wwhh5, wwhh6, zzhh,
     3		 aaee,azee,zaee,zzee,aauu,azuu,zauu,zzuu,
     3		 CCee, CCuu, CCee6, CCuu6,
     3		 azz,zzztens
      logical lzs
      common /tenszz/ aazz(0:3,0:3,3),azzz(0:3,0:3,3),
     3  	      zazz(0:3,0:3,3),zzzz(0:3,0:3,3),
     3  	      wwzz5(0:3,0:3,3),wwzz6(0:3,0:3,3),
     3                wwhh5(0:3,0:3,3), wwhh6(0:3,0:3,3),
     3                zzhh(0:3,0:3,3),
     3		      aaee(0:3,0:3,2,3),azee(0:3,0:3,2,3),
     3		      zaee(0:3,0:3,2,3),zzee(0:3,0:3,2,3),
     3		      aauu(0:3,0:3,2,3),azuu(0:3,0:3,2,3),
     3		      zauu(0:3,0:3,2,3),zzuu(0:3,0:3,2,3),
     3  	      CCee(0:3,0:3,2,3),  CCuu(0:3,0:3,2,3),
     3  	      CCee6(0:3,0:3,2,3), CCuu6(0:3,0:3,2,3),
     3		      azz(0:5),zzztens(0:5),
     8		      lzs(8)


cc  lzs remembers whether or not an fsign combination needs to be calculated 
cc     lzs(1) = .false.      fsign(1)=fsign(2) = -1 has been calcd and stored
cc     lzs(2) = .false.      fsign(1)=fsign(2) = +1 has been calcd and stored
cc     lzs(3) = .false.      fsign(3)=fsign(4) = -1 has been calcd and stored
cc     lzs(4) = .false.      fsign(3)=fsign(4) = +1 has been calcd and stored
cc lzs(1,+1) to lzs(4,+1) are for gsign=+1 (outgoing gluon) in case of real emission
cc for incoming gluons use lzs(1,-1) to lzs(4,-1):
cc     lzs(5) = .false.      fsign(1) = -1 fsign(2) = +1 has been calcd ...
cc     lzs(6) = .false.      fsign(1) = +1 fsign(2) = -1 has been calcd ...
cc     lzs(7) = .false.      fsign(3) = -1 fsign(4) = +1 has been calcd ...
cc     lzs(8) = .false.      fsign(3) = +1 fsign(4) = -1 has been calcd ...


c for testing purposes
      double complex tampz,tamprz
      integer ipz,iprz
c      integer is1,is3
      common /zztest/  tampz(600,-1:1,-1:1,6), ipz
      common /zzrtest/ tamprz(3000,-1:1,-1:1,-1:1,6), iprz               
