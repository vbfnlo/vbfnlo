      function theta(x)
      implicit none
      double precision theta,x
      if (x.gt.0d0) then
         theta = 1d0
      else
         theta = 0d0
      endif
      end

      function sgn(x)
      implicit none
      double precision sgn,x
      if (x.gt.0d0) then
         sgn = +1d0
      else
         sgn = -1d0
      endif
      end

c     imaginary part of a/b  
c     IF a/b > 0 then returns 0
c     IF a/b < 0 then
c        if a < 0 then return +1
c        else return -1
      function im_part(a,b)
      implicit none
      double precision a,b,im_part
      if ((a/b).gt.0d0) then
         im_part = 0d0
      else
         if (a.lt.0d0) then
            im_part = +1d0
         else
            im_part = -1d0
         endif
      endif
      
      end

      function eta3(im_a,im_b,im_ab)
      implicit none
      double precision eta3,im_a,im_b,im_ab
      double precision theta
      external theta
      eta3 = theta(-im_a)*theta(-im_b)*theta(im_ab)
     #     - theta(im_a)*theta(im_b)*theta(-im_ab)
      end


***********************************************************************
        FUNCTION ETA(C1,C2)                                            
***********************************************************************
*       COMPLEX ETA-FUNKTION                                           
*---------------------------------------------------------------------*
*       8.06.90    ANSGAR DENNER                                       
***********************************************************************
        IMPLICIT   LOGICAL(A-Z)                                        
        COMPLEX*16 ETA,C1,C2                                           
        REAL*8     PI,IM1,IM2,IM12                                     
                                                                       
        PI     = 4D0*DATAN(1D0)                                        
        IM1    = DIMAG(C1)                                             
        IM2    = DIMAG(C2)                                             
        IM12   = DIMAG(C1*C2)                                          
                                                                       
        IF(IM1.LT.0D0.AND.IM2.LT.0D0.AND.IM12.GT.0D0) THEN             
            ETA = DCMPLX(0D0,2D0*PI)                                   
        ELSE IF (IM1.GT.0D0.AND.IM2.GT.0D0.AND.IM12.LT.0D0) THEN       
            ETA = DCMPLX(0D0,-2D0*PI)                                  
        ELSE                                                           
            ETA = DCMPLX(0D0)                                          
        END IF                                                         
        END                                                            

**********************************************************************
        FUNCTION V_ETA(C1,C2)                                            
***********************************************************************
*       COMPLEX V_ETA-FUNKTION                                           
*---------------------------------------------------------------------*
*       8.06.90    ANSGAR DENNER                                       
***********************************************************************
        implicit none                                       
        COMPLEX*16 V_ETA,C1,C2                                           
        REAL*8     PI,IM1,IM2,IM12                                     
                                                                       
        PI     = 4D0*DATAN(1D0)                                        
        IM1    = DIMAG(C1)                                             
        IM2    = DIMAG(C2)                                             
        IM12   = DIMAG(C1*C2)                                          
                                                                       
        IF(IM1.LT.0D0.AND.IM2.LT.0D0.AND.IM12.GT.0D0) THEN             
            V_ETA = DCMPLX(0D0,2D0*PI)                                   
        ELSE IF (IM1.GT.0D0.AND.IM2.GT.0D0.AND.IM12.LT.0D0) THEN       
            V_ETA = DCMPLX(0D0,-2D0*PI)                                  
        ELSE                                                           
            V_ETA = DCMPLX(0D0)                                          
        END IF                                                         
        END                                                            

***********************************************************************
        FUNCTION ETAS(Y,R,RS)                                            
***********************************************************************
*       MODIFIED ETA-FUNKTION                                           
*---------------------------------------------------------------------*
*       18.1.94   SD                                       
***********************************************************************
        implicit none
        COMPLEX*16 V_ETA,ETAS,Y,R,RS
        REAL*8     PI,IMY,IMRS
                                                                       
        PI     = 4D0*DATAN(1D0)                                        

	IF( DIMAG(R).NE.0D0 ) THEN
	    ETAS = V_ETA(Y,R)
	ELSE	    
	    IF( DREAL(R).GT.0D0 ) THEN
		ETAS = DCMPLX(0D0,0D0)
	    ELSE
	 	IMY  = DIMAG(Y)
		IMRS = DIMAG(RS)
		ETAS = 2D0*DCMPLX(0D0,PI)*(
     *			(1D0+SIGN(1D0,-IMY))*(1D0+SIGN(1D0,-IMRS))-
     *			(1D0+SIGN(1D0, IMY))*(1D0+SIGN(1D0, IMRS))
     *					  )/4D0
	    ENDIF
	ENDIF
        END


***********************************************************************
        FUNCTION SQE(A,B,C)                                            
***********************************************************************
*       SOLUTION OF QUADRATIC EQUATION				      *
*---------------------------------------------------------------------*
*       13.1.92  SD						      *
***********************************************************************
        implicit none
        COMPLEX*16 A,B,C,SQE,X1,X2

	X1=(-B+SQRT(B**2-4D0*A*C))/2D0/A
	X2=(-B-SQRT(B**2-4D0*A*C))/2D0/A

	IF (ABS(X1).GT.ABS(X2)) THEN
	   SQE=X1
	ELSE
	   SQE=X2
	ENDIF

        END         


      function delta(x,y)
      implicit none
      real * 8 delta
      integer x,y
      if (x.eq.y) then
         delta = 1d0
      else
         delta = 0d0
      endif
      end

