      real*8 alphaDip,alphaDipFF
      parameter (alphaDip = 1d0)
      parameter (alphaDipFF = 1d0 )
      real*8 GaugeAccuracy

#ifdef WITH_QUAD
      parameter (GaugeAccuracy = 1d-2)
#else
      parameter (GaugeAccuracy = 1d-1)
#endif
c      integer NF
c      parameter (NF = 5)
