c----------------------------------------------------------------
c
c  Define common block containing all Kaluza-Klein coupling 
c  constants, masses, widths and tower height
c
c----------------------------------------------------------------

      integer  kkhz,kkha,kkhw,kkhh,kkh,kkw,kkz,kka
      parameter(kkz=20,kka=20,kkw=20,kkh=20)
      double precision w1mass_k(0:kkw),z1mass_k(0:kkz),z1width_k(0:kkz),w1width_k(0:kkw),
     &                 a1mass_k(0:kka),a1width_k(0:kka),h1mass_k(0:kkh),h1width_k(0:kkh),
     &                 gwwww_k(0:kkw,0:kkw,0:kkw,0:kkw),gwwzz_k(0:kkw,0:kkw,0:kkz,0:kkz),
     &                 gwwaa_k(0:kkw,0:kka,0:kkz,0:kka),gwwaz_k(0:kkw,0:kkw,0:kka,0:kkz),
     &                 gwwz_k(0:kkw,0:kkw,0:kkz),gwwa_k(0:kkw,0:kkw,0:kka),
     &                 gwwh_k(0:kkw,0:kkw,0:kkh),gzzh_k(0:kkz,0:kkz,0:kkh)
      logical hswitchww,hswitchwz,hswitchzz,hswitch
      logical hswitchwww,hswitchwwz,hswitchzzw, hswitchzzz
	logical mhlmswitch

      common /COUPL_KK/ w1mass_k,z1mass_k,z1width_k,w1width_k,
     &                 a1mass_k,a1width_k,h1mass_k,h1width_k,
     &                 gwwww_k,gwwzz_k,
     &                 gwwaa_k,gwwaz_k,
     &                 gwwz_k,gwwa_k,
     &                 gwwh_k,gzzh_k,
     &                 kkhz,kkha,kkhw,kkhh,
     &                 hswitchww,hswitchwz,hswitchzz,hswitch,
     &               hswitchwww,hswitchwwz,hswitchzzw,hswitchzzz,
     &               mhlmswitch
     
	
