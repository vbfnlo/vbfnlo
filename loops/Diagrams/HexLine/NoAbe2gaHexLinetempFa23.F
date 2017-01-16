c************** Calling the Fa functions*************************
       subroutine NoAbe2gaHexLineFa23(p1mup2,p1mup3,p1mup5,p1mup6,p2mu
     -   p2,p2mup3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,p
     -   4mup3,p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mup
     -   3,p6mup5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup6
     -   ,mup5mup6,Fa23)
       IMPLICIT NONE
      Complex*16   p1mup2, p1mup3, p1mup5, p1mup6, p2mup2, p2mup3, 
     -          p2mup5, p2mup6, p3mup2, p3mup3, p3mup5, p3mup6, 
     -          p4mup2, p4mup3, p4mup5, p4mup6, p5mup2, p5mup3, 
     -          p5mup5, p5mup6, p6mup2, p6mup3, p6mup5, p6mup6
       Complex*16   mup2mup3, mup2mup5, mup2mup6, mup3mup5, mup3mup6, 
     -          mup5mup6
        common/NoAbe2gaHexLineFhlFunctions/F
       COMMON/NoAbe2gaHexLineInvariants/p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,
     -   s12,s23,s34,s45,s56,s16,s123,s234,s345
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
      Complex*16 Fa(872),F(9706)
      Real*8 P(9687) 
        COMMON/NoAbe2gaHexLinePFunctions/P
       Complex*16 Fa23(529:552)
       COMMON/NoAbe2gaHexLineFaFunctions/Fa
       Fa23(529) = p6mup3*(p5mup5*p5mup6*F(2841)+p1mup6*p4mup5*F(2847)
     -   -p2mup6*p4mup5*F(2848)-p3mup6*p4mup5*F(2851)-4*(p1mup6*p2mup
     -   5*F(2844)+p2mup5*p2mup6*F(2845)+p1mup6*p3mup5*F(2849)+p2mup5
     -   *p3mup6*F(2853))-2*(p3mup5*p4mup6*F(2850)+p2mup5*p4mup6*F(28
     -   52)-p1mup5*(p4mup6*F(2846)+2*(p1mup6*F(2842)+p2mup6*F(2843)+
     -   p3mup6*F(2854)))))
       Fa23(530) = -(p2mup6*p5mup5*p6mup3*F(2857))-p2mup3*p5mup5*p6mup
     -   6*F(2860)-p4mup5*p5mup6*p6mup3*F(2868)+p4mup6*p5mup5*p6mup3*
     -   F(2869)-p2mup3*p4mup5*p6mup6*F(2874)-2*(p1mup6*p5mup5*p6mup3
     -   *F(2856)-p1mup5*p2mup3*p6mup6*F(2858)+p2mup3*p2mup5*p6mup6*F
     -   (2862)-p3mup6*p5mup5*p6mup3*F(2867)+p2mup3*p3mup5*p6mup6*F(2
     -   870)-p4mup3*p4mup6*p6mup5*F(2872)-p2mup3*p4mup6*p6mup5*F(287
     -   5))+p4mup5*p4mup6*p6mup3*F(2876)+4*(p1mup5*p5mup6*p6mup3*F(2
     -   855)+p1mup6*p2mup3*p6mup5*F(2859)-p2mup5*p5mup6*p6mup3*F(286
     -   1)+p2mup3*p2mup6*p6mup5*F(2863)+p2mup3*p5mup6*p6mup5*F(2864)
     -   +p2mup6*p3mup5*p6mup3*F(2865)-p3mup5*p5mup6*p6mup3*F(2866)+p
     -   2mup3*p3mup6*p6mup5*F(2871)+p5mup3*p5mup6*p6mup5*F(2873)-p1m
     -   up3*p2mup6*p6mup5*F(2877)+p3mup5*p3mup6*p6mup3*F(2878)+p1mup
     -   3*p1mup6*p6mup5*F(2879)+p3mup3*p3mup6*p6mup5*F(2880))
       Fa23(531) = p6mup5*(p2mup6*p4mup3*F(2885)+2*(p4mup6*p5mup3*F(28
     -   91)-p3mup3*p4mup6*F(2892)+p1mup3*(p4mup6*F(2886)+2*(p5mup6*F
     -   (2882)+p3mup6*F(2893))))+4*(p2mup6*p5mup3*F(2884)+p3mup6*p4m
     -   up3*F(2887)+p3mup3*p5mup6*F(2888)-p3mup6*p5mup3*F(2889)+p4mu
     -   p3*p5mup6*F(2890)+p2mup6*p3mup3*F(2894)+p1mup6*(-(p4mup3*F(2
     -   881))+p5mup3*F(2883)+p3mup3*F(2895))))
       Fa23(532) = p6mup3*p6mup5*(-(p4mup6*F(2898))+2*(p1mup6*F(2896)-
     -   p2mup6*F(2897)-p5mup6*F(2899)+p3mup6*F(2900)))
       Fa23(533) = p6mup6*(p2mup3*p6mup5*F(2901)+p2mup5*p5mup3*F(2902)
     -   )
       Fa23(534) = p6mup6*(p4mup3*(p5mup5*F(2903)+p2mup5*F(2909)+8*(p1
     -   mup5*F(2905)-p3mup5*F(2911))+4*p4mup5*F(2917))+4*(p5mup3*(p5
     -   mup5*F(2904)+2*(p1mup5*F(2907)-p3mup5*F(2913))-p4mup5*F(2916
     -   ))+p1mup3*(p5mup5*F(2906)+2*(p1mup5*F(2908)+p2mup5*F(2910)-p
     -   3mup5*F(2914))+p4mup5*F(2918))+p3mup3*(-(p5mup5*F(2912))+2*(
     -   p1mup5*F(2915)+p3mup5*F(2919)+p2mup5*F(2920))+p4mup5*F(2921)
     -   )))
       Fa23(535) = p6mup3*p6mup6*(p5mup5*F(2922)+2*(p1mup5*F(2923)+p3m
     -   up5*F(2924)+p2mup5*F(2925))+p4mup5*F(2926))
       Fa23(536) = p6mup5*p6mup6*(p4mup3*F(2927)+p5mup3*F(2928)+8*(p1m
     -   up3*F(2929)-p3mup3*F(2930)-p6mup3*F(2931)))
       Fa23(537) = mup3mup6*(p2mup2*F(2934)+p5mup2*F(2935)+p4mup2*F(29
     -   38))+p1mup6*p2mup3*p5mup2*F(2942)+p2mup3*p5mup2*p6mup6*F(294
     -   4)+p1mup6*p3mup2*p5mup3*F(2947)+p1mup6*p3mup3*p5mup2*F(2948)
     -   +p3mup3*p5mup2*p6mup6*F(2949)+p1mup6*p5mup2*p5mup3*F(2950)+p
     -   5mup2*p5mup3*p5mup6*F(2951)+p2mup6*p4mup2*p5mup3*F(2956)+p2m
     -   up6*p4mup3*p5mup2*F(2957)+p1mup6*p4mup2*p5mup3*F(2958)+p1mup
     -   6*p4mup3*p5mup2*F(2959)+p4mup2*p5mup3*p6mup6*F(2960)+p5mup2*
     -   p5mup3*p6mup6*F(2963)+p2mup3*p2mup6*p5mup2*F(2965)+mup2mup3*
     -   p6mup6*F(2969)+p2mup6*p3mup2*p5mup3*F(2970)+p3mup2*p5mup3*p5
     -   mup6*F(2971)+p3mup2*p5mup3*p6mup6*F(2972)+p3mup3*p5mup2*p5mu
     -   p6*F(2973)+p1mup6*p3mup2*p3mup3*F(2976)+mup3mup6*p3mup2*F(29
     -   78)+p2mup6*p3mup3*p5mup2*F(2985)+p2mup3*p5mup2*p5mup6*F(2986
     -   )+p2mup6*p5mup2*p5mup3*F(2987)+p2mup2*(p5mup3*(p5mup6*F(2939
     -   )+p1mup6*F(2940)+p6mup6*F(2941)+p2mup6*F(2943))+p2mup3*p5mup
     -   6*F(2946)+p1mup6*(p2mup3*F(2945)+4*p4mup3*F(2953)+p3mup3*F(2
     -   962))+p2mup3*p2mup6*F(2964)+p2mup6*p3mup3*F(2979)-4*(p4mup3*
     -   p5mup6*F(2954)-p4mup3*p6mup6*F(2955)-p3mup3*p6mup6*F(2980)+p
     -   2mup6*p4mup3*F(2989))+p3mup3*p5mup6*F(2991)+p2mup3*p6mup6*F(
     -   2992))+p1mup6*p2mup3*p3mup2*F(2993)+p3mup2*p3mup3*p5mup6*F(2
     -   996)+p2mup3*p3mup2*p5mup6*F(3000)+p2mup6*p3mup2*p3mup3*F(300
     -   1)+8*(mup2mup3*p2mup6*F(2967)+mup2mup3*p1mup6*F(2968)+mup2mu
     -   p3*p5mup6*F(3004))+p2mup3*p2mup6*p3mup2*F(3005)+4*(mup2mup6*
     -   (p2mup3*F(2932)+p5mup3*F(2933)+p4mup3*F(2936)+p3mup3*F(2937)
     -   )+p4mup3*p5mup2*p5mup6*F(2952)+p4mup3*p5mup2*p6mup6*F(2961)+
     -   p4mup2*p5mup3*p5mup6*F(2966)+p2mup6*p3mup2*p4mup3*F(2974)+p2
     -   mup6*p3mup3*p4mup2*F(2975)+p3mup2*p4mup3*p5mup6*F(2977)+p1mu
     -   p6*p3mup3*p4mup2*F(2981)+p3mup2*p4mup3*p6mup6*F(2982)+p4mup2
     -   *p4mup3*p5mup6*F(2983)+p4mup2*p4mup3*p6mup6*F(2984)+p1mup6*p
     -   2mup3*p4mup2*F(2988)+p2mup3*p2mup6*p4mup2*F(2990)+p1mup6*p4m
     -   up2*p4mup3*F(2994)+p1mup6*p3mup2*p4mup3*F(2995)-p2mup6*p4mup
     -   2*p4mup3*F(2997)+p2mup3*p4mup2*p5mup6*F(2998)+p2mup3*p4mup2*
     -   p6mup6*F(2999)+p3mup2*p3mup3*p6mup6*F(3002)+p3mup3*p4mup2*p5
     -   mup6*F(3003)+p2mup3*p3mup2*p6mup6*F(3006)+p3mup3*p4mup2*p6mu
     -   p6*F(3007))
       Fa23(538) = mup3mup6*p1mup2*F(3008)+p1mup3*(4*mup2mup6*F(3009)+
     -   p1mup6*p2mup2*F(3010)+p2mup2*p5mup6*F(3011))+p1mup3*(p2mup2*
     -   p2mup6*F(3016)+p2mup2*p6mup6*F(3023)+p1mup6*(p5mup2*F(3013)+
     -   p3mup2*F(3022)-4*p4mup2*F(3026))+p5mup2*(p6mup6*F(3015)+p2mu
     -   p6*F(3018)+p5mup6*F(3029))+p2mup6*p3mup2*F(3033)+p3mup2*p5mu
     -   p6*F(3034)-4*(p2mup6*p4mup2*F(3027)-p4mup2*p5mup6*F(3032)+p3
     -   mup2*p6mup6*F(3035)+p4mup2*p6mup6*F(3038)))+p1mup2*(p2mup3*p
     -   5mup6*F(3021)+p3mup3*p5mup6*F(3028)+p5mup3*(p6mup6*F(3014)+p
     -   2mup6*F(3017)+p5mup6*F(3031))+p2mup3*p2mup6*F(3036)+p2mup6*p
     -   3mup3*F(3037)+p1mup6*(p5mup3*F(3012)+p2mup3*F(3019)+p3mup3*F
     -   (3020)+4*p4mup3*F(3024)+p1mup3*F(3039))+p1mup3*p5mup6*F(3040
     -   )+p1mup3*p6mup6*F(3041)+p2mup3*p6mup6*F(3042)+p1mup3*p2mup6*
     -   F(3044)+4*(p2mup6*p4mup3*F(3025)-p4mup3*p5mup6*F(3030)+p4mup
     -   3*p6mup6*F(3043)+p3mup3*p6mup6*F(3045)))
       Fa23(539) = -4*p1mup2*p4mup6*p5mup3*F(3049)+p2mup2*(p6mup3*p6mu
     -   p6*F(3048)-4*(p5mup6*p6mup3*F(3047)+p4mup6*p5mup3*F(3050)+p4
     -   mup6*p6mup3*F(3051)))+p1mup2*p6mup3*(p6mup6*F(3053)-4*(p5mup
     -   6*F(3052)-p4mup6*F(3054)))
       Fa23(540) = p4mup3*p6mup2*(-4*(p1mup6*F(3055)-p5mup6*F(3056))+p
     -   6mup6*F(3057))
       Fa23(541) = p3mup6*p4mup3*p5mup2*F(3063)+p4mup3*p4mup6*p5mup2*F
     -   (3068)-4*(p2mup2*p3mup6*p5mup3*F(3060)+p2mup2*p3mup6*p6mup3*
     -   F(3061)+p1mup3*p2mup2*(p3mup6*F(3059)+p4mup6*F(3062))+p2mup2
     -   *p3mup6*p4mup3*F(3064)+p2mup2*p2mup3*p3mup6*F(3066)-p2mup2*p
     -   4mup3*p4mup6*F(3067)-p2mup2*p2mup3*p4mup6*F(3070)-p3mup2*p4m
     -   up3*p4mup6*F(3073)-p3mup6*p4mup2*p4mup3*F(3074)+p4mup2*p4mup
     -   3*p4mup6*F(3075)+p3mup2*p3mup6*p4mup3*F(3076)+p1mup2*(p4mup6
     -   *(p4mup3*F(3069)+p2mup3*F(3072)-p1mup3*F(3078))+p3mup6*(p5mu
     -   p3*F(3058)+p2mup3*F(3065)+p4mup3*F(3071)-p6mup3*F(3077)-p1mu
     -   p3*F(3079))))
       Fa23(542) = -(p5mup3*p6mup2*(p2mup6*F(3081)+p3mup6*F(3083)+p5mu
     -   p6*F(3084)))+p2mup6*p4mup3*p6mup2*F(3085)+p1mup3*p3mup6*p5mu
     -   p2*F(3086)-p4mup6*p5mup3*p6mup2*F(3087)+p1mup3*p4mup6*p5mup2
     -   *F(3088)-p4mup6*p5mup2*p5mup3*F(3089)-p3mup3*p3mup6*p5mup2*F
     -   (3090)-p3mup3*p5mup6*p6mup2*F(3091)+p3mup6*p5mup2*p5mup3*F(3
     -   092)-p1mup2*p3mup3*p3mup6*F(3093)+p3mup6*p4mup3*p6mup2*F(309
     -   4)+p3mup6*p4mup2*p5mup3*F(3095)-p4mup3*p4mup6*p6mup2*F(3096)
     -   -p4mup2*p4mup6*p5mup3*F(3097)-p2mup3*p3mup6*p5mup2*F(3098)-p
     -   2mup3*p5mup6*p6mup2*F(3099)-p2mup3*p4mup6*p5mup2*F(3100)-p3m
     -   up3*p4mup6*p5mup2*F(3101)-p1mup2*p3mup3*p4mup6*F(3102)-p2mup
     -   3*p4mup6*p6mup2*F(3104)+p2mup2*(p2mup6*p6mup3*F(3082)+p3mup3
     -   *p4mup6*F(3103)-p3mup3*p3mup6*F(3105))-p3mup3*p3mup6*p6mup2*
     -   F(3106)-4*(mup2mup3*p3mup6*F(3107)+mup2mup3*p4mup6*F(3108))-
     -   p2mup3*p2mup6*p6mup2*F(3109)-p2mup3*p3mup6*p6mup2*F(3110)-p2
     -   mup6*p3mup3*p6mup2*F(3111)-p3mup3*p4mup6*p6mup2*F(3112)+p1mu
     -   p6*p6mup3*(p2mup2*F(3080)-p1mup2*F(3113))+p1mup2*p2mup6*p6mu
     -   p3*F(3114)+p1mup3*p3mup6*p4mup2*F(3115)+p3mup2*p3mup6*p5mup3
     -   *F(3116)+p1mup3*p4mup2*p4mup6*F(3117)+p3mup2*p4mup6*p5mup3*F
     -   (3118)+p2mup3*p3mup6*p4mup2*F(3119)+p2mup3*p4mup2*p4mup6*F(3
     -   120)+p3mup3*p3mup6*p4mup2*F(3121)+p3mup3*p4mup2*p4mup6*F(312
     -   2)+p1mup3*p3mup2*p3mup6*F(3123)+p1mup3*p3mup2*p4mup6*F(3124)
     -   +p2mup3*p3mup2*p3mup6*F(3125)+p2mup3*p3mup2*p4mup6*F(3126)+p
     -   3mup2*p3mup3*p3mup6*F(3127)+p3mup2*p3mup3*p4mup6*F(3128)
       Fa23(543) = mup2mup6*p6mup3*F(3130)+p2mup6*p5mup2*p6mup3*F(3132
     -   )+p3mup6*p5mup2*p6mup3*F(3133)+p5mup2*p5mup6*p6mup3*F(3134)+
     -   p1mup3*p5mup6*p6mup2*F(3135)+p5mup6*p6mup2*p6mup3*F(3136)-p1
     -   mup6*p6mup3*(p5mup2*F(3131)+p3mup2*F(3137)+p4mup2*F(3138))-p
     -   2mup6*p4mup2*p6mup3*F(3139)+p4mup6*p5mup2*p6mup3*F(3140)-p4m
     -   up2*p5mup6*p6mup3*F(3141)+p1mup3*p3mup6*p6mup2*F(3142)+p3mup
     -   6*p6mup2*p6mup3*F(3143)+p2mup6*p3mup2*p6mup3*F(3144)-p3mup2*
     -   p5mup6*p6mup3*F(3145)+p1mup3*p2mup6*p6mup2*F(3146)+p2mup6*p6
     -   mup2*p6mup3*F(3147)+p1mup3*p4mup6*p6mup2*F(3148)+p4mup6*p6mu
     -   p2*p6mup3*F(3149)+p3mup6*p4mup2*p6mup3*F(3150)-p4mup2*p4mup6
     -   *p6mup3*F(3151)+p3mup2*p3mup6*p6mup3*F(3152)+p3mup2*p4mup6*p
     -   6mup3*F(3153)
       Fa23(544) = p1mup6*p6mup2*(p5mup3*F(3154)+p1mup3*F(3155)-4*(p6m
     -   up3*F(3156)-p3mup3*F(3157)))
       Fa23(545) = p6mup3*p6mup6*(p5mup2*F(3158)+p3mup2*F(3159)+p4mup2
     -   *F(3160))
       Fa23(546) = p6mup2*p6mup6*(p5mup3*F(3161)+p2mup3*F(3162)+4*p6mu
     -   p3*F(9616)*P(102)+(p3mup3*F(9624)+p1mup3*F(9626))*P(5539))
       Fa23(547) = p6mup3*F(3166)+2*(p2mup3*F(3163)+p5mup3*F(3164)+p1m
     -   up3*F(3165)+p4mup3*F(3167)+p3mup3*F(3168))
       Fa23(548) = mup5mup6*(p2mup2*F(3169)+p5mup2*F(3170)+p3mup2*F(31
     -   71)+p4mup2*F(3172))
       Fa23(549) = p4mup2*p5mup5*p6mup6*F(3178)+p2mup2*p2mup6*p5mup5*F
     -   (3181)+p2mup2*p2mup6*p6mup5*F(3182)+p2mup5*p4mup2*p6mup6*F(3
     -   185)+p2mup6*p4mup2*p5mup5*F(3186)-8*(p4mup2*p6mup5*p6mup6*F(
     -   3179)+p3mup2*p6mup5*p6mup6*F(3194))+p4mup2*p5mup5*p5mup6*F(3
     -   195)+p2mup5*p5mup2*p5mup6*F(3202)+p2mup2*p2mup5*p2mup6*F(320
     -   3)+p1mup5*(p4mup2*p6mup6*F(3175)+p2mup2*p2mup6*F(3180)+p1mup
     -   6*(p2mup2*F(3188)-4*(p4mup2*F(3174)+p3mup2*F(3189)))+p2mup2*
     -   p5mup6*F(3199)+4*(p2mup6*p4mup2*F(3183)-p3mup2*p6mup6*F(3190
     -   )-p4mup2*p5mup6*F(3200)+p5mup2*p5mup6*F(3201)+p2mup6*p3mup2*
     -   F(3206)+p3mup2*p5mup6*F(3207)))+p1mup6*(p4mup2*p5mup5*F(3176
     -   )+p2mup5*p4mup2*F(3184)+p3mup2*p5mup5*F(3191)+4*(p4mup2*p6mu
     -   p5*F(3177)-p3mup2*p6mup5*F(3192)+p2mup5*p3mup2*F(3208)))+p2m
     -   up5*p3mup2*p6mup6*F(3209)+p2mup6*p3mup2*p5mup5*F(3210)+p3mup
     -   2*p5mup6*p6mup5*F(3213)+4*(p2mup6*p4mup2*p6mup5*F(3187)-p3mu
     -   p2*p5mup5*p6mup6*F(3193)+p4mup2*p5mup6*p6mup5*F(3196)-p5mup2
     -   *p5mup5*p5mup6*F(3197)+p5mup2*p5mup6*p6mup5*F(3198)-p2mup5*p
     -   2mup6*p4mup2*F(3204)-p2mup5*p4mup2*p5mup6*F(3205)+p2mup6*p3m
     -   up2*p6mup5*F(3211)+p3mup2*p5mup5*p5mup6*F(3212)-p2mup5*p3mup
     -   2*p5mup6*F(3214)-p2mup5*p2mup6*p3mup2*F(3215))
       Fa23(550) = 2*mup5mup6*p6mup2*F(3219)+mup2mup6*(p1mup5*F(3216)+
     -   p2mup5*F(3217)+p5mup5*F(3218)+p6mup5*F(3220))+p1mup6*p2mup2*
     -   p5mup5*F(3222)+p2mup2*p5mup5*p5mup6*F(3223)+p2mup2*p5mup5*p6
     -   mup6*F(3224)+p1mup6*p5mup2*p6mup5*F(3233)-8*(p2mup2*p6mup5*p
     -   6mup6*F(3226)-p5mup2*p6mup5*p6mup6*F(3236))+p1mup6*p2mup2*p2
     -   mup5*F(3237)+p1mup6*p2mup2*p6mup5*F(3238)+p1mup5*(p2mup2*p6m
     -   up6*F(3221)+p5mup2*(p6mup6*F(3230)+4*(p1mup6*F(3229)+p2mup6*
     -   F(3240))))+4*(p2mup2*p5mup6*p6mup5*F(3225)-p2mup6*p5mup2*p5m
     -   up5*F(3227)+p2mup6*p5mup2*p6mup5*F(3228)+p1mup6*p2mup5*p5mup
     -   2*F(3231)+p1mup6*p5mup2*p5mup5*F(3232)+p2mup5*p5mup2*p6mup6*
     -   F(3234)-p5mup2*p5mup5*p6mup6*F(3235)+p2mup2*p2mup5*p6mup6*F(
     -   3239)-p2mup5*p2mup6*p5mup2*F(3241))+p2mup2*p2mup5*p5mup6*F(3
     -   242)
       Fa23(551) = mup2mup6*p4mup5*F(3243)+p2mup2*p4mup5*(p5mup6*F(325
     -   0)-4*p4mup6*F(3255)+p6mup6*F(3258))+p1mup2*(p1mup5*p6mup6*F(
     -   3245)+p5mup5*p5mup6*F(3246)+p5mup5*p6mup6*F(3247)-8*p6mup5*p
     -   6mup6*F(3249)+p2mup5*p6mup6*F(3251)+p2mup6*p5mup5*F(3252)+p4
     -   mup5*p6mup6*F(3257)+4*(p5mup6*p6mup5*F(3248)+p4mup5*p5mup6*F
     -   (3256)+p4mup5*p4mup6*F(3259))+p1mup6*(p5mup5*F(3244)+p1mup5*
     -   F(3253)+p6mup5*F(3254)+p2mup5*F(3260))+p1mup5*p2mup6*F(3261)
     -   +p2mup6*p6mup5*F(3262)+p1mup5*p5mup6*F(3263)+p2mup5*p5mup6*F
     -   (3264)+p2mup5*p2mup6*F(3265))
       Fa23(552) = mup2mup6*p3mup5*F(3269)+mup2mup5*(4*(p1mup6*F(3266)
     -   -p2mup6*F(3267))+p6mup6*F(3268)+p5mup6*F(3270))-16*p2mup5*p5
     -   mup6*p6mup2*F(3273)+p2mup2*p3mup5*p6mup6*F(3280)+p3mup5*p4mu
     -   p2*p6mup6*F(3283)+p1mup2*p2mup6*p4mup5*F(3285)+p1mup2*p3mup5
     -   *p6mup6*F(3286)+p4mup2*p4mup5*p6mup6*F(3287)+p4mup5*p5mup2*p
     -   6mup6*F(3289)+p1mup6*(p2mup2*(p4mup5*F(3271)+4*p3mup5*F(3297
     -   ))+4*(p2mup5*p6mup2*F(3272)-p4mup2*p4mup5*F(3277)-p4mup5*p5m
     -   up2*F(3279)-p3mup2*p4mup5*F(3281)+p3mup5*p4mup2*F(3282)-p1mu
     -   p2*p4mup5*F(3284)+p3mup5*p5mup2*F(3291)-p3mup2*p3mup5*F(3294
     -   )+p1mup2*p3mup5*F(3298)))-8*p2mup5*p6mup2*p6mup6*F(3299)+p3m
     -   up2*p4mup5*p6mup6*F(3302)+p1mup2*p3mup5*p5mup6*F(3307)+4*(p2
     -   mup2*p2mup6*p4mup5*F(3274)+p3mup5*p5mup2*p6mup6*F(3275)-p2mu
     -   p6*p4mup2*p4mup5*F(3276)+p2mup6*p4mup5*p5mup2*F(3278)+p4mup5
     -   *p5mup2*p5mup6*F(3288)-p4mup2*p4mup5*p5mup6*F(3290)-p2mup6*p
     -   3mup5*p4mup2*F(3292)-p2mup2*p3mup5*p5mup6*F(3293)+p3mup2*p3m
     -   up5*p6mup6*F(3295)-p2mup2*p2mup6*p3mup5*F(3296)+p2mup6*p3mup
     -   2*p4mup5*F(3300)-p3mup2*p4mup5*p5mup6*F(3301)-p3mup5*p4mup2*
     -   p5mup6*F(3303)-p2mup6*p3mup5*p5mup2*F(3304)-p3mup2*p3mup5*p5
     -   mup6*F(3305)-p3mup5*p5mup2*p5mup6*F(3306)-p2mup6*p3mup2*p3mu
     -   p5*F(3308))+p1mup2*p2mup6*p3mup5*F(3309)
       Return
       End
