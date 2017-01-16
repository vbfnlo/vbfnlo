***********************************************************************
c  author: Michael Kubocz
c  date     : 14.08.06    
c  modified : 20.07.08
***********************************************************************
      real*8 q_dot_q
      dimension q_dot_q(4,4)
      complex*16 q_dot_ep
      dimension q_dot_ep(4,4)
      complex*16 ep_dot_ep
      dimension ep_dot_ep(4,4)
      real*8 ep
      dimension ep(4,4)
      real*8 det,Xi(4,4)

      COMMON/auxiliary_calculation_GG/q_dot_ep,ep_dot_ep,ep,q_dot_q,
     &                                det,Xi

      complex*16 E0t(12),E0b(12),E1it(4,12),E1ib(4,12),
     &           DDEijtop(4,46,12),DDEijbottom(4,46,12)

      COMMON/form_factors_E0_GG/E0t,E0b,E1it,E1ib,
     &       DDEijtop,DDEijbottom



