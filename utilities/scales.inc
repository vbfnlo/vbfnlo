c
c alfas, scales etc. Notation for the arrays is
c
c  variable(i,L):   i=1 upper line,  i=2 lower line
c                   L=1 for use with main amplitude^2
c                   L=2,3 for use with upper/lower line finite
c                         Born^2 terms in 3 parton run 
c
      real*8 mufsq(2,25), mursq(2,25), als(3,25), xir, xif
      real*8 muf_user, mur_user
      integer id_mur, id_muf
      common /cscales/ mufsq, mursq, als, xir, xif, muf_user, mur_user, id_mur, id_muf

