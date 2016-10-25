
c  photon polarization vectors
      double complex sp2au, sp2ae
      common /sp2aapol/ sp2au(1:6), sp2ae(1:6)

c  leptonic tensors
      double complex sp2wwaa5, sp2wwaa6, sp2aaaa, sp2azaa, sp2zaaa,
     1               sp2zzaa


      common /tenssp2aa/  sp2wwaa5(0:3,0:3,3), sp2wwaa6(0:3,0:3,3),
     1        sp2aaaa(0:3,0:3,3), sp2azaa(0:3,0:3,3), sp2zaaa(0:3,0:3,3),
     2        sp2zzaa(0:3,0:3,3)


      double complex sp2wwww5, sp2wwww6, sp2aaww, sp2azww, sp2zaww,
     1               sp2zzww

      common /tenssp2ww/  sp2wwww5(0:3,0:3,3), sp2wwww6(0:3,0:3,3),
     1        sp2aaww(0:3,0:3,3), sp2azww(0:3,0:3,3), sp2zaww(0:3,0:3,3),
     2        sp2zzww(0:3,0:3,3)


      double complex sp2wwzz5, sp2wwzz6, sp2aazz, sp2azzz, sp2zazz,
     1               sp2zzzz

      common /tenssp2zz/  sp2wwzz5(0:3,0:3,3), sp2wwzz6(0:3,0:3,3),
     1        sp2aazz(0:3,0:3,3), sp2azzz(0:3,0:3,3), sp2zazz(0:3,0:3,3),
     2        sp2zzzz(0:3,0:3,3)


