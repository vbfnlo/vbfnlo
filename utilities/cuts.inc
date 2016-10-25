      double precision ptmin_def, ymax_def, rll_min, rll_max, ymax_part
      double precision jetgenktp
      integer min_jets, njet_min
      double precision rjjmin, rjl_min, ymax_l, ptlmin,ptlmax, etajjmin, ptmiss_min
      double precision dely_jlep,dely_jveto
      double precision mdijmin,delmll,mdijmax
      double precision ptmin_def_default(max_jets)
      double precision ptlmin_def_default(max_v), ptlmax_def_default(max_v)

      logical ysign, lrapidgap,grapidgap,lveto
      double precision ptmin_veto,ymax_veto
      double precision rgg_min, rgg_max, rjg_min, rlg_min, ymax_g
      double precision ptgmin,ptgmax, dely_jphot
      double precision isol_phopart, eff_phopart,ptinmin
      double precision mllmin,mllmax,mlAmin, mlAmax, MTmax,MTmin,mtwa_min,mza_min
      logical mllosonly

      common/jetdef/rjjmin,ymax_part,njet_min
      common/jetalgo/jetgenktp

      double precision eta_central, ptmin_tag_1, ptmin_tag_2, ptmin_central, semilep_v_mass_range
      integer semilep_fatjet, def_tagjet, semilep_reconst_v_had
      logical hard_central, semileptonic, decayPartonsCombined(max_kin), semilepAlwaysVBFCuts

      integer tag1,tag2
      logical vbfprocess, lneutrinoCut

      common/tempNeutrinoCut/ lneutrinoCut
      
      common/basiccuts/ptmin_def(max_jets), ymax_def, ptmiss_min,mtwa_min,mza_min, min_jets


      common/photoncuts/isol_phopart, eff_phopart

      common/vbfcuts/rgg_min,rgg_max,rjg_min,rlg_min,ymax_g,ptgmin(max_v),dely_jphot,
     &               rll_min,rll_max,rjl_min,ymax_l,ptlmax(max_v),ptlmin(max_v),etajjmin,
     &               mdijmin,mdijmax,ymax_veto,ptmin_veto,delmll,dely_jlep,
     &               dely_jveto,grapidgap,lveto,ysign,lrapidgap,ptgmax(max_v)

      common/tribosoncuts/mllmin,mllmax,mlAmin,mlAmax,MTmax,MTmin,ptinmin,
     &               mllosonly

      common /semilepcuts/ eta_central, ptmin_tag_1, ptmin_tag_2, ptmin_central,
     &                     semilep_v_mass_range,
     &                     def_tagjet, semilep_reconst_v_had, semilep_fatjet, 
     &                     hard_central, semileptonic, decayPartonsCombined,
     &                     semilepAlwaysVBFCuts

      common /tagjetdef/tag1,tag2,vbfprocess

