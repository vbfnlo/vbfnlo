c----------------------------
c  user functions
c----------------------------

      double precision function MZPf(MW,MZ,MWP)
      implicit none
      double precision MW, MZ, MWP, x, t
      double precision xf, tf      
      x = xf(MW,MWP)
      t = tf(MW,MZ,MWP)
      MZPf = (MWP*sqrt(2d0 + x**2 + t**2*x**2 + 
     -      sqrt(4d0 + x**4 - 2*t**2*x**4 + t**4*x**4)))
     -   /sqrt(2d0 + x**2 + sqrt(4d0 + x**4))
      return
      end

      double precision function MFPf(MF,MD,MW,MWP)
      implicit none
      double precision MF, MD, MW, MWP, EpL, EpR
      double precision EpLf, EpRf
      EpL = EpLf(MW,MWP)
      EpR = EpRf(MF,MD,MW,MWP)
      MFPf = (sqrt(1d0 + EpL**2 + EpR**2 + 
     -      sqrt(1d0 + 2d0*EpL**2 + EpL**4 + 2d0*EpR**2 - 
     -        2d0*EpL**2*EpR**2 + EpR**4))*MD)/sqrt(2d0)
      return
      end


      double precision function xf(MW,MWP)
      implicit none
      double precision MW, MWP
      xf =  ((1d0 - sqrt(1d0 + MW**4/MWP**4 - (6d0*MW**2)/MWP**2) + 
     -     MW**2/MWP**2)*MWP)/(2d0*MW)
      return
      end


      double precision function EpLf(MW,MWP)
      implicit none
      double precision MW, MWP, x, xf
      x = xf(MW,MWP)
      EpLf = (sqrt(2d0)*x)/sqrt(2d0 - x**2 + sqrt(4d0 + x**4))
      return
      end

      double precision function tf(MW,MZ,MWP)
      implicit none
      double precision MW, MZ, MWP, x, xf
      x = xf(MW,MWP)
      tf =   sqrt(-(((-1d0 + mW**2/mZ**2)*mZ**2*
     -        ((2d0*mW**4*(x**2 + x**4))/mZ**4 + 
     -          x**2*(2d0 + x**2 - sqrt(4d0 + x**4)) + 
     -          (mW**2*(-x**6 + 
     -               2d0*x**2*(-4d0 + sqrt(4d0 + x**4)) + 
     -               x**4*(-4d0 + sqrt(4d0 + x**4)) + 
     -               2d0*(-2d0 + sqrt(4d0 + x**4))))/mZ**2))/
     -      (mW**2*x**2*(x**2 + 
     -          (mW**4*(1d0 + x**2)**2)/mZ**4 - 
     -          (mW**2*(2d0 + 3d0*x**2 + x**4))/mZ**2))))/sqrt(2d0)
      return
      end

      double precision function EpRf(MF,MD,MW,MWP)
      implicit none
      double precision MF, MD, MW, MWP, EpL, EpLf
      EpL = EpLf(MW,MWP)
      EpRf = sqrt((MF**2*(MF**2 - MD**2 - EpL**2*MD**2))/
     -   (MD**2*(MF**2 - EpL**2*MD**2)))
      return
      end

      double precision function gf(ee,MW,MZ,MWP)
      implicit none
      double precision ee, MW, MZ, MWP, x, xf, t, tf
      x = xf(MW,MWP)
      t = tf(MW,MZ,MWP)
      gf = ee*sqrt(1d0 + 1d0/t**2 + x**2)
      return
      end

      double precision function gtf(ee,MW,MZ,MWP)
      implicit none
      double precision ee, MW, MZ, MWP, x, xf, t, tf
      x = xf(MW,MWP)
      t = tf(MW,MZ,MWP)
      gtf = ee*sqrt(1d0 + 1d0/x**2 + 1d0/(t**2*x**2))
      return
      end

      double precision function gpf(ee,MW,MZ,MWP)
      implicit none
      double precision ee, MW, MZ, MWP, x, xf, t, tf
      x = xf(MW,MWP)
      t = tf(MW,MZ,MWP)
      gpf = ee*sqrt(1d0 + t**2 + t**2*x**2)
      return
      end

      double precision function vW0f(MW,MWP)
      implicit none
      double precision MW, MWP, x, xf
      x = xf(MW,MWP)
      vW0f = (2d0 - x**2 + sqrt(4d0 + x**4))/
     -  (sqrt(2d0)*sqrt(4d0 + x**4 + 2d0*sqrt(4d0 + x**4) - 
     -      x**2*sqrt(4d0 + x**4)))
      return
      end
      
      double precision function vW1f(MW,MWP)
      implicit none
      double precision MW, MWP, x, xf
      x = xf(MW,MWP)
      vW1f = 1d0/
     -  sqrt(1d0 + (2d0 - x**2 + Sqrt(4d0 + x**4))**2/
     -     (4d0*x**2))
      return
      end

      double precision function vWP0f(MW,MWP)
      implicit none
      double precision MW, MWP, x, xf
      x = xf(MW,MWP)
      vWP0f =  -(-2d0 + x**2 + sqrt(4d0 + x**4))/
     -  (2d0*x*sqrt(1d0 + 
     -      (-2d0 + x**2 + sqrt(4d0 + x**4))**2/
     -       (4d0*x**2)))
      return
      end
      
      double precision function vWP1f(MW,MWP)
      implicit none
      double precision MW, MWP, x, xf
      x = xf(MW,MWP)
      vWP1f =   1d0/
     -  sqrt(1d0 + (-2d0 + x**2 + sqrt(4d0 + x**4))**
     -      2/(4d0*x**2))
      return
      end


      double precision function vZ0f(MW,MZ,MWP)
      implicit none
      double precision MW, MZ, MWP, x, xf, t, tf
      x = xf(MW,MWP)
      t = tf(MW,MZ,MWP)
      vZ0f =  ((-1d0 + t**2)*x**2 + 
     -    sqrt(4d0 + (-1d0 + t**2)**2*x**4))/
     -  (t*sqrt(4d0 + (-2d0 + (-1d0 + t**2)*x**2 + 
     -          sqrt(4d0 + (-1d0 + t**2)**2*x**4))**2/
     -       (t**2*x**2) + 
     -      ((-1d0 + t**2)*x**2 + 
     -          sqrt(4d0 + (-1d0 + t**2)**2*x**4))**2/
     -       t**2))
      return
      end
		
      
      double precision function vZ1f(MW,MZ,MWP)
      implicit none
      double precision MW, MZ, MWP, x, xf, t, tf
      x = xf(MW,MWP)
      t = tf(MW,MZ,MWP)
      vZ1f =   -((-2d0 + (-1d0 + t**2)*x**2 + 
     -      sqrt(4d0 + (-1d0 + t**2)**2*x**4))/
     -    (sqrt(2d0)*sqrt(4d0 + (-1d0 + t**2)**2*x**6 - 
     -        2d0*sqrt(4d0 + (-1d0 + t**2)**2*x**4) + 
     -        (-1d0 + t**2)*x**4*
     -         (-1d0 + t**2 + 
     -           sqrt(4d0 + (-1d0 + t**2)**2*x**4)) + 
     -        x**2*(4d0 - 
     -           sqrt(4d0 + (-1d0 + t**2)**2*x**4) + 
     -           t**2*sqrt(4d0 + (-1d0 + t**2)**2*x**4)))
     -      ))
      return
      end

      double precision function vZ2f(MW,MZ,MWP)
      implicit none
      double precision MW, MZ, MWP, x, xf, t, tf
      x = xf(MW,MWP)
      t = tf(MW,MZ,MWP)
      vZ2f = -2d0/
     -  sqrt(4d0 + (-2d0 + (-1d0 + t**2)*x**2 + 
     -        sqrt(4d0 + (-1d0 + t**2)**2*x**4))**2/
     -     (t**2*x**2) + 
     -    ((-1d0 + t**2)*x**2 + 
     -        sqrt(4d0 + (-1d0 + t**2)**2*x**4))**2/t**2)  
      return
      end
     
      double precision function vZP0f(MW,MZ,MWP)
      implicit none
      double precision MW, MZ, MWP, x, xf, t, tf
      x = xf(MW,MWP)
      t = tf(MW,MZ,MWP)
      vZP0f =   ((-1d0 + t**2)*x**2 - 
     -    sqrt(4d0 + (-1d0 + t**2)**2*x**4))/
     -  (t*sqrt(4d0 + (-((-1d0 + t**2)*x**2) + 
     -          sqrt(4d0 + (-1d0 + t**2)**2*x**4))**2/
     -       t**2 + (2d0 - (-1d0 + t**2)*x**2 + 
     -          sqrt(4d0 + (-1d0 + t**2)**2*x**4))**2/
     -       (t**2*x**2)))
      return
      end
		
		
      double precision function vZP1f(MW,MZ,MWP)
      implicit none
      double precision MW, MZ, MWP, x, xf, t, tf
      x = xf(MW,MWP)
      t = tf(MW,MZ,MWP)
      vZP1f =   (2d0 - (-1d0 + t**2)*x**2 + 
     -    sqrt(4d0 + (-1d0 + t**2)**2*x**4))/
     -  (t*x*sqrt(4d0 + 
     -      (-((-1d0 + t**2)*x**2) + 
     -          sqrt(4d0 + (-1d0 + t**2)**2*x**4))**2/
     -       t**2 + (2d0 - (-1d0 + t**2)*x**2 + 
     -          sqrt(4d0+ (-1d0 + t**2)**2*x**4))**2/
     -       (t**2*x**2)))
      return
      end

      double precision function vZP2f(MW,MZ,MWP)
      implicit none
      double precision MW, MZ, MWP, x, xf, t, tf
      x = xf(MW,MWP)
      t = tf(MW,MZ,MWP)
      vZP2f =  -2d0/
     -  sqrt(4d0 + (-((-1d0 + t**2)*x**2) + 
     -        sqrt(4d0 + (-1d0 + t**2)**2*x**4))**2/
     -     t**2 + (2d0 - (-1d0 + t**2)*x**2 + 
     -        sqrt(4d0 + (-1d0 + t**2)**2*x**4))**2/
     -     (t**2*x**2))
      return
      end
      
      double precision function vL0f(MW,MWP)
      implicit none
      double precision MW, MWP, EpL, EpLf
      EpL = EpLf(MW,MWP)
      vL0f =  -(1d0 - EpL**2 + sqrt((1d0 + EpL**2)**2))/
     -  (2d0*EpL*sqrt(1d0 + 
     -      (1d0 - EpL**2 + sqrt((1d0 + EpL**2)**2))**2/
     -       (4d0*EpL**2)))
      return 
      end
		
      double precision function vL1f(MW,MWP)
      implicit none
      double precision MW, MWP, EpL, EpLf
      EpL = EpLf(MW,MWP)
      vL1f =   1d0/sqrt(1d0 + (1d0 - EpL**2 + 
     -        sqrt((1d0 + EpL**2)**2))**2/(4d0*EpL**2))
      return
      end
		
      double precision function vLt0f(MW,MWP,mt,MD)
      implicit none
      double precision MW, MWP, mt, MD, EpL, EpLf, EpR, EpRf
      EpL = EpLf(MW,MWP)
      EpR = EpRf(mt,MD,MW,MWP)
      vLt0f =  -(1d0 - EpL**2 + EpR**2 + 
     -     sqrt(-4d0*EpL**2*EpR**2 + 
     -       (1d0 + EpL**2 + EpR**2)**2))/
     -  (2d0*EpL*sqrt(1d0 + 
     -      (1d0 - EpL**2 + EpR**2 + 
     -          sqrt(-4d0*EpL**2*EpR**2 + 
     -            (1d0 + EpL**2 + EpR**2)**2))**2/
     -       (4d0*EpL**2)))
      return
      end
		
      double precision function vLt1f(MW,MWP,mt,MD)
      implicit none
      double precision MW, MWP, mt, MD, EpL, EpLf, EpR, EpRf
      EpL = EpLf(MW,MWP)
      EpR = EpRf(mt,MD,MW,MWP)
      vLt1f =  1d0/sqrt(1d0 + (1d0 - EpL**2 + EpR**2 + 
     -        sqrt(-4d0*EpL**2*EpR**2 + 
     -          (1d0 + EpL**2 + EpR**2)**2))**2/
     -     (4d0*EpL**2)) 
      return
      end
		
      double precision function vRt1f(MW,MWP,mt,MD)
      implicit none
      double precision MW, MWP, mt, MD, EpL, EpLf, EpR, EpRf
      EpL = EpLf(MW,MWP)
      EpR = EpRf(mt,MD,MW,MWP)
      vRt1f = -(-1d0 - EpL**2 + EpR**2 + 
     -     sqrt(-4d0*EpL**2*EpR**2 + 
     -       (1d0 + EpL**2 + EpR**2)**2))/
     -  (2d0*EpR*sqrt(1d0 + 
     -      (-1d0 - EpL**2 + EpR**2 + 
     -          sqrt(-4d0*EpL**2*EpR**2 + 
     -            (1d0 + EpL**2 + EpR**2)**2))**2/
     -       (4d0*EpR**2)))
      return
      end
		
      double precision function vRt2f(MW,MWP,mt,MD)
      implicit none
      double precision MW, MWP, mt, MD, EpL, EpLf, EpR, EpRf
      EpL = EpLf(MW,MWP)
      EpR = EpRf(mt,MD,MW,MWP)
      vRt2f = 1d0/sqrt(1d0 + (-1d0 - EpL**2 + EpR**2 + 
     -        sqrt(-4d0*EpL**2*EpR**2 + 
     -          (1d0 + EpL**2 + EpR**2)**2))**2/
     -     (4d0*EpR**2))
      return
      end
		
      double precision function vLP0f(MW,MWP)
      implicit none
      double precision MW, MWP, EpL, EpLf
      EpL = EpLf(MW,MWP)
      vLP0f = (1d0 - EpL**2 - sqrt((1d0 + EpL**2)**2))/
     -  (2d0*EpL*sqrt(1d0 + 
     -      (-1d0 + EpL**2 + sqrt((1d0 + EpL**2)**2))**2/
     -       (4d0*EpL**2)))
      return
      end
		
      double precision function vLP1f(MW,MWP)
      implicit none
      double precision MW, MWP, EpL, EpLf
      EpL = EpLf(MW,MWP)
      vLP1f = -(1d0/
     -    sqrt(1d0 + (-1d0 + EpL**2 + 
     -          sqrt((1d0 + EpL**2)**2))**2/(4d0*EpL**2)
     -      ))
      return
      end
			
      double precision function vLTP0f(MW,MWP,mt,MD)
      implicit none
      double precision MW, MWP, mt, MD, EpL, EpLf, EpR, EpRf
      EpL = EpLf(MW,MWP)
      EpR = EpRf(mt,MD,MW,MWP)
      vLTP0f = (1d0 - EpL**2 + EpR**2 - 
     -    sqrt(-4d0*EpL**2*EpR**2 + 
     -      (1d0 + EpL**2 + EpR**2)**2))/
     -  (2d0*EpL*sqrt(1d0 + 
     -      (-1d0 + EpL**2 - EpR**2 + 
     -          sqrt(-4d0*EpL**2*EpR**2 + 
     -            (1d0 + EpL**2 + EpR**2)**2))**2/
     -       (4d0*EpL**2)))
      return
      end
		
      double precision function vLTP1f(MW,MWP,mt,MD)
      implicit none
      double precision MW, MWP, mt, MD, EpL, EpLf, EpR, EpRf
      EpL = EpLf(MW,MWP)
      EpR = EpRf(mt,MD,MW,MWP)
      vLTP1f = -(1d0/
     -    sqrt(1d0 + (-1d0 + EpL**2 - EpR**2 + 
     -          sqrt(-4d0*EpL**2*EpR**2 + 
     -            (1d0 + EpL**2 + EpR**2)**2))**2/
     -       (4d0*EpL**2)))
      return
      end
			
      double precision function vRTP1f(MW,MWP,mt,MD)
      implicit none
      double precision MW, MWP, mt, MD, EpL, EpLf, EpR, EpRf
      EpL = EpLf(MW,MWP)
      EpR = EpRf(mt,MD,MW,MWP)
      vRTP1f = (-1d0 - EpL**2 + EpR**2 - 
     -    sqrt(-4d0*EpL**2*EpR**2 + 
     -      (1d0 + EpL**2 + EpR**2)**2))/
     -  (2d0*EpR*sqrt(1d0 + 
     -      (1d0 + EpL**2 - EpR**2 + 
     -          sqrt(-4d0*EpL**2*EpR**2 + 
     -            (1d0 + EpL**2 + EpR**2)**2))**2/
     -       (4d0*EpR**2)))
      return
      end		
		
      double precision function vRTP2f(MW,MWP,mt,MD)
      implicit none
      double precision MW, MWP, mt, MD, EpL, EpLf, EpR, EpRf
      EpL = EpLf(MW,MWP)
      EpR = EpRf(mt,MD,MW,MWP)
      vRTP2f = -(1d0/
     -    sqrt(1d0 + (1d0 + EpL**2 - EpR**2 + 
     -          sqrt(-4d0*EpL**2*EpR**2 + 
     -            (1d0 + EpL**2 + EpR**2)**2))**2/
     -       (4d0*EpR**2)))
      return
      end		
		
      
c----------------------------
c  user functions
c----------------------------
