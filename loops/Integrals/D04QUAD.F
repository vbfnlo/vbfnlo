
************************************************************************
c        FUNCTION D04(P1,P2,P3,P4,P12,P23,M1,M2,M3,M4)
        FUNCTION D04QUAD1(P1t,P2t,P3t,P4t,P12t,P23t,M1t,M2t,M3t,M4t)
************************************************************************
*  SCALAR 4-POINT FUNCTION WITH AT LEAST ONE MASS ZERO                 *
*  P1,P2,P3,P4 = SQUARED EXTERNAL MOMENTA			       *
*  P12 = (p1+p2)**2,  P23 = (p2+p3)**2				       *
*----------------------------------------------------------------------*
*  2.1.92  SD	         					       *
*  rearrangement to exploit massless external momenta   14.3.01  DZ    *
*  Modified: Michael Kubocz                                            *
*  Interception of NANs e.g. caused by Log(0) etc. (see below)         *
************************************************************************
c        IMPLICIT REAL*16 (A-Z)
        implicit none
	REAL*16 M(4),P(4,4),K(4,4)
        real*16 pi,eps,eps1
        real*16 im1,im2
        real*16 m1,m2,m3,m4
        real*16 m1t,m2t,m3t,m4t
        real*16 m02,m12,m22,m32,m42
        real*16 mm0,mm1,mm2,mm3,mm4 
        real*16 p1,p2,p3,p4,p12,p23
        real*16 p1t,p2t,p3t,p4t,p12t,p23t
        real*16 q0,q1,q2,q3,q4,q00,q12,q23
	COMPLEX*32 A1,A2,A3,A4,SWAP
	COMPLEX*32 SS(4), XX(2), X(2,4),RS(4,4)
	COMPLEX*32 S0(4),XX0(2),X0(2,4), R(4,4),G(2)
	COMPLEX*32 AA,BB,CC,DD,IEPS,H,HH,L1,L2,L3,L4
	INTEGER I,J, i1,i2,i3,i4
        complex*32 D1,D2,D3,D4
        COMPLEX*32 D04QUAD1,QUADCDILOG,ETA_QUAD,SQE_QUAD,ETAS_QUAD
        EXTERNAL QUADCDILOG,ETA_QUAD,SQE_QUAD,ETAS_QUAD

        D1=CMPLX(0q0,0q0,16)
        D2=CMPLX(0q0,0q0,16)
        D3=CMPLX(0q0,0q0,16)
        D4=CMPLX(0q0,0q0,16)
        eps1=1d-7

        if(abs(P1t).le.eps1) then
           P1=0q0
        else
           P1=P1t
        endif
        if(abs(P2t).le.eps1) then
           P2=0q0
        else
           P2=P2t
        endif
        if(abs(P3t).le.eps1) then
           P3=0q0
        else
           P3=P3t
        endif
        if(abs(P4t).le.eps1) then
           P4=0q0
        else
           P4=P4t
        endif
        if(abs(P12t).le.eps1) then
           P12=0q0
        else
           P12=P12t
        endif
        if(abs(P23t).le.eps1) then
           P23=0q0
        else
           P23=P23t
        endif
        if(abs(M1t).le.eps1) then
           M1=0q0
        else
           M1=M1t
        endif
        if(abs(M2t).le.eps1) then
           M2=0q0
        else
           M2=M2t
        endif
        if(abs(M3t).le.eps1) then
           M3=0q0
        else
           M3=M3t
        endif
        if(abs(M4t).le.eps1) then
           M4=0q0
        else
           M4=M4t
        endif

        MM1=M1
        MM2=M2
        MM3=M3
        MM4=M4
        M12=M1*M1
        M22=M2*M2
        M32=M3*M3
        M42=M4*M4
        Q1=P1
        Q2=P2
        Q3=P3
	Q4=P4
        Q12=P12
        Q23=P23

C	IS AT LEAST ONE MASS ZERO ???
	IF (MM1*MM2*MM3*MM4.NE.0q0) GOTO 130
C	PERMUTATE UNTIL MM3=0D0
	GOTO 20
10	CONTINUE
	MM0=MM1
	MM1=MM2
	MM2=MM3
	MM3=MM4
	MM4=MM0
	M02=M12
	M12=M22
	M22=M32
	M32=M42
	M42=M02
	Q00=Q12
	Q12=Q23
	Q23=Q00
	Q0=Q1
	Q1=Q2
	Q2=Q3
	Q3=Q4
	Q4=Q0
20	IF (MM3.NE.0D0) GOTO 10
C	ONLY MM3 IS ZERO
	IF (MM1*MM2*MM4.NE.0D0) GOTO 30
C	ONLY MM3 AND MM4 ARE ZERO ==> 3->2, 4->3...
	IF ((MM1*MM2.NE.0D0).AND.(MM4.EQ.0D0)) GOTO 10
C	ONLY MM2 AND MM3 ARE ZERO
	IF ((MM1*MM4.NE.0D0).AND.(MM2.EQ.0D0)) GOTO 40
	WRITE(*,*)'CASE OF THIS SPECIAL D0-FUNCTION NOT IMPLEMENTED!'
	STOP

C	****** NO MASS EQUAL TO ZERO ******
130	CONTINUE
	EPS=1q-38
	IEPS=CMPLX(0q0,EPS,16)
c check for massless external momentum: excellent candidate for p13,
c leading to r13 >=1 and real.
c$$$        if (q1*q2*q3*q4.eq.0d0) then !org
        if (abs(q1*q2*q3*q4).le.eps1) then
c$$$           if (q2.eq.0d0) then
           if (abs(q2).le.eps1) then
              I1 = 2
              I2 = 3
              I3 = 1
              I4 = 4
c$$$           elseif (q1.eq.0d0) then
           elseif (abs(q1).le.eps1) then
              I1 = 1
              I2 = 3
              I3 = 2
              I4 = 4
c$$$           elseif (q3.eq.0d0) then
           elseif (abs(q3).le.eps1) then
              I1 = 2
              I2 = 4
              I3 = 1
              I4 = 3
           else
              I1 = 1
              I2 = 4
              I3 = 2
              I4 = 3
           endif
           M(i1)=MM1
           M(i2)=MM2
           M(i3)=MM3
           M(i4)=MM4
           P(i1,i2)=Q1
           P(i3,i2)=Q2
           P(i3,i4)=Q3
           P(i1,i4)=Q4
           P(i1,i3)=Q12
           P(i3,i1)=Q12
           P(i2,i4)=Q23
           P(i4,i2)=Q23
	ELSEIF( ABS((MM1**2+MM3**2-Q12)/MM1/MM3).LT.2q0 ) THEN
C	R13 WOULD BE NOT REAL. -> PERMUTATION! -> R(2,4) IS NOT REAL.
	   M(1)=MM2
	   M(2)=MM3
	   M(3)=MM4
	   M(4)=MM1
	   P(1,2)=Q2
	   P(1,3)=Q23
	   P(1,4)=Q1
	   P(2,3)=Q3
	   P(2,4)=Q12
	   P(3,4)=Q4
	ELSE
C	R(1,3) IS REAL.
	   M(1)=MM1
	   M(2)=MM2
	   M(3)=MM3
	   M(4)=MM4
	   P(1,2)=Q1
	   P(1,3)=Q12
	   P(1,4)=Q4
	   P(2,3)=Q2
	   P(2,4)=Q23
	   P(3,4)=Q3
	ENDIF

	DO 11 J=2,4
	DO 11 I=1,J-1
	K(I,J)=(M(I)**2+M(J)**2-P(I,J))/M(I)/M(J)
	R(I,J) =SQE_QUAD(CMPLX(1q0,0q0,16),CMPLX(-K(I,J),0q0,16),
     *	            CMPLX(1q0,0q0,16))
        IF( K(I,J).LT.-2q0 ) THEN
c        IF( DIMAG(R(I,J)).EQ.0D0 ) THEN
	   RS(I,J)=SQE_QUAD(CMPLX(1q0,0q0,16),CMPLX(-K(I,J),EPS,16),
     *	               CMPLX(1q0,0q0,16))
	ELSE
	   RS(I,J)=R(I,J)
	ENDIF
11	CONTINUE


	SS(1)=RS(1,2)
	SS(2)=RS(2,3)
	SS(3)=RS(3,4)
	SS(4)=RS(1,4)
	S0(1)=R(1,2)
	S0(2)=R(2,3)
	S0(3)=R(3,4)
	S0(4)=R(1,4)
	AA=K(3,4)/R(2,4)+R(1,3)*K(1,2)-K(1,4)*R(1,3)/R(2,4)-K(2,3)
	BB=(R(2,4)-1q0/R(2,4))*(R(1,3)-1q0/R(1,3))
     *		+K(1,2)*K(3,4)-K(1,4)*K(2,3)
	CC=K(1,2)/R(1,3)+R(2,4)*K(3,4)-K(1,4)*R(2,4)/R(1,3)-K(2,3)
	DD=K(2,3)-R(1,3)*K(1,2)-R(2,4)*K(3,4)+R(1,3)*R(2,4)*K(1,4)
	XX(1)=SQE_QUAD(AA,BB,CC+IEPS*DD)
	XX(2)=(CC+IEPS*DD)/AA/XX(1)
	XX0(1)=SQE_QUAD(AA,BB,CC)
	XX0(2)=CC/AA/XX0(1)
c	IF (ABS(DREAL(XX0(1)-XX(2))).LT.ABS(DREAL(XX0(1)-XX(1)))) THEN
	IF (ABS(XX0(1)-XX(2)).LT.ABS(XX0(1)-XX(1))) THEN
	  SWAP  =XX0(1)
	  XX0(1)=XX0(2)
	  XX0(2)=SWAP
	ENDIF

	DO 12 I=1,2
	G(I)  =SIGN( 1Q0,REAL(AA*(XX(I)-XX(3-I))) )
	 X(I,1)= XX(I)/R(2,4)
	X0(I,1)=XX0(I)/R(2,4)
	 X(I,2)= XX(I)/R(2,4)*R(1,3)
	X0(I,2)=XX0(I)/R(2,4)*R(1,3)
	 X(I,3)= XX(I)*R(1,3)
	X0(I,3)=XX0(I)*R(1,3)
	 X(I,4)= XX(I)
	X0(I,4)=XX0(I)
12	CONTINUE

	D04QUAD1 = CMPLX(0q0,0q0,16)
	DO 13 I=1,2
	DO 13 J=1,4
	A1 = 1q0+X0(I,J)*S0(J) + ABS(1q0+X0(I,J)*S0(J))*IEPS*
     *				  SIGN(1q0,IMAG(X(I,J)*SS(J)))
	A2 = 1q0+X0(I,J)/S0(J) + ABS(1q0+X0(I,J)/S0(J))*IEPS*
     *				  SIGN(1q0,IMAG(X(I,J)/SS(J)))

c  org. code:
c$$$           D04 = D04 + (-1D0)**(I+J)*(
c$$$     &          LI2(A1)+ETA(-X(I,J),SS(J))*LOG(A1)
c$$$     &          +LI2(A2)+ETA(-X(I,J),1D0/SS(J))*LOG(A2))

        if(abs(ETA_QUAD(-X(I,J),SS(J))).ne.0q0) then
           D1=ETA_QUAD(-X(I,J),SS(J))*LOG(A1)
        else
           D1=CMPLX(0q0,0q0,16)
        endif
        if(abs(ETA_QUAD(-X(I,J),1q0/SS(J))).ne.0q0) then
           D2=ETA_QUAD(-X(I,J),1q0/SS(J))*LOG(A2)
        else
           D2=CMPLX(0q0,0q0,16)
        endif
        D04QUAD1=D04QUAD1+(-1q0)**(I+J)*(QUADCDILOG(A1)+QUADCDILOG(A2)+D1+D2)
c   The enquiry avoids occurrence of NANs causing by LOG(A1) for A1=0 
c   and LOG(A2) for A2=0. At that points also ETA_QUAD(-X(I,J),1D0/SS(J)) 
c   or ETA_QUAD(-X(I,J),SS(J) are 0. (Michael Kubocz)

13	CONTINUE


c        print*,'DIMAG(R(1,3))',DIMAG(R(1,3))
c	IF( DIMAG(R(1,3)).EQ.0D0 ) THEN !org (makes troubles in squark pentagons)
	IF( abs(IMAG(R(1,3))).le.eps ) THEN
	DO 14 I=1,2
	   A1 = (K(1,3)-2q0*R(1,3))/XX0(I)
     *		      -R(1,3)*K(1,4)+K(3,4)
     	   A2 = ((K(2,4)-2q0*R(2,4))*R(1,3)*XX0(I)
     *		      -R(2,4)*K(3,4)+K(2,3))/DD
	   A3 = (K(1,3)-2q0*R(1,3))*R(2,4)/XX0(I)
     *		      -R(1,3)*K(1,2)+K(2,3)
	   A4 = ((K(2,4)-2q0*R(2,4))*XX0(I)
     *		      -R(2,4)*K(1,4)+K(1,2))/DD
	   L1 = LOG( A1-ABS(A1)*IEPS )
     	   L2 = LOG( A2+ABS(A2)*IEPS*G(I)*SIGN(1q0,REAL(R(1,3))
     *				        	  *IMAG(RS(2,4))) )
	   L3 = LOG( A3-ABS(A3)*IEPS )
	   L4 = LOG( A4+ABS(A4)*IEPS*G(I)*SIGN(1q0,IMAG(RS(2,4))) )

c org. code:
c$$$	   D04 = D04 
c$$$     &         + (3D0-2D0*I)*(
c$$$     *		       ETAS_QUAD( -XX(I),R(1,3),RS(1,3) )
c$$$     *		          *( LOG(R(1,3)*XX(I)) + L1 + L2 )
c$$$     *		     + ETAS_QUAD( -XX(I),1D0/R(2,4),1D0/RS(2,4) )
c$$$     *		          *( LOG(XX(I)/R(2,4)) + L3 + L4 )
c$$$     *		     - ( ETAS_QUAD( -XX(I),R(1,3)/R(2,4),RS(1,3)/RS(2,4) )
c$$$     *		       + ETA( RS(1,3),1D0/RS(2,4) )                  )
c$$$     *		        *( LOG(XX(I)*R(1,3)/R(2,4)) + L3 + L2 )
c$$$     *	  	     + ETA( RS(1,3),1D0/RS(2,4) )
c$$$     *		       *ETAS_QUAD(-XX(I),-R(1,3)/R(2,4),-RS(1,3)/RS(2,4))   )

            if(abs(ETAS_QUAD(-XX(I),R(1,3),RS(1,3))).ne.0q0) then
               D1=ETAS_QUAD(-XX(I),R(1,3),RS(1,3))*(LOG(R(1,3)*XX(I))+L1+L2)
            else
               D1=CMPLX(0q0,0q0,16)
            endif
            if(abs(ETAS_QUAD(-XX(I),1q0/R(2,4),1q0/RS(2,4))).ne.0q0) then
               D2=ETAS_QUAD(-XX(I),1q0/R(2,4),1q0/RS(2,4))*(LOG(XX(I)/R(2,4))
     &           +L3+L4)
            else
               D2=CMPLX(0q0,0q0,16)
            endif
            if((abs(ETAS_QUAD(-XX(I),R(1,3)/R(2,4),RS(1,3)/RS(2,4))).ne.0q0)
     &          .or.(abs(ETA_QUAD(RS(1,3),1q0/RS(2,4))).ne.0q0)) then 
               D3=-(ETAS_QUAD(-XX(I),R(1,3)/R(2,4),RS(1,3)/RS(2,4))
     &            +ETA_QUAD(RS(1,3),1q0/RS(2,4)))*(LOG(XX(I)*R(1,3)/R(2,4))+L3+L2)
            else
               D3=CMPLX(0q0,0q0,16)
            endif
            D4=ETA_QUAD(RS(1,3),1q0/RS(2,4))*ETAS_QUAD(-XX(I),-R(1,3)/R(2,4),-RS(1,3)/RS(2,4))
            D04QUAD1=D04QUAD1+(3q0-2q0*I)*(D1+D2+D3+D4)
c   The enquiry avoids occurrence of NANs causing by LOG(0). At that points 
c   also ETA(...) and ETAS_QUAD(...) are 0. (Michael Kubocz)


14	CONTINUE
	ELSE
	DO 15 I=1,2
	   L1 = LOG( R(2,4)/XX0(I)+XX0(I)/R(2,4)+K(1,2)
     *		     -XX0(I)/R(2,4)*EPS*BB*G(I) )
	   L2 = LOG( R(1,3)*XX0(I)+1Q0/XX0(I)/R(1,3)+K(3,4)
     *		     -XX0(I)*R(1,3)*EPS*BB*G(I) )
	   L3 = LOG( R(1,3)/R(2,4)*XX0(I)+R(2,4)/XX0(I)/R(1,3)+K(2,3)
     *		     -XX0(I)*R(1,3)/R(2,4)*EPS*BB*G(I) )

	   D04QUAD1 = D04QUAD1 
     &          + (3q0-2q0*I)*(
     *		     ETA_QUAD(-XX(I),1q0/R(2,4))
     *		      *( LOG(XX(I)/R(2,4)) + L1 )
     *		    +ETA_QUAD(-XX(I),R(1,3))
     *		      *( LOG(R(1,3)*XX(I)) + L2 )
     *		    -( ETA_QUAD(-XX(I),R(1,3)/R(2,4))
     *		      +ETA_QUAD(R(1,3),1Q0/R(2,4)) )
     *		      *( LOG(XX(I)*R(1,3)/R(2,4)) + L3 )
     *	  	    +ETA_QUAD(R(1,3),1Q0/R(2,4))
     *		      *ETA_QUAD(-XX(I),-R(1,3)/R(2,4))
     *		      *(1q0-G(I)*SIGN(1q0,REAL(BB)))	    )

                  
15	CONTINUE
	ENDIF

	D04QUAD1 = D04QUAD1/M(1)/M(2)/M(3)/M(4)/AA/(XX(1)-XX(2))
	RETURN
 
C--->	***************** SPEZIELL ( --> T.SACK-PROMOTION )
C	D1=Q12-M12
C	D2=Q2 -M22
C	D3=Q3 -M42
C	IF ((D1*D2.LE.0D0).OR.(D2*D3.LE.0D0)) THEN
C	   WRITE(*,*) 'THE CASE OF DIFFERENT SIGNS OF THE D1,D2,D3'
C	   WRITE(*,*) 'IN D04(...) IS NOT IMPLEMENTED !!!'
C	   STOP
C	ENDIF
C	NM1=ABS(MM1/D1)
C	NM2=ABS(MM2/D2)
C	NM3=ABS(MM4/D3)
C	NP1=Q2/D2**2+Q12/D1**2+(Q1-Q2-Q12)/D1/D2
C	NP2=Q2/D2**2+ Q3/D3**2+(Q23-Q2-Q3)/D2/D3
C	NP3=Q3/D3**2+Q12/D1**2+(Q4-Q3-Q12)/D1/D3
C	D04=C04(NP1,NP2,NP3,NM1,NM2,NM3)/D1/D2/D3

C	*************** ALLGEMEIN


C	****** ONLY MM3 IS ZERO ******
30	CONTINUE
	EPS=1q-47
	IEPS=CMPLX(0Q0,EPS,16)
	M(1)=MM1
	M(2)=MM2
	M(3)=10Q0
	M(4)=MM4
	P(1,2)=Q1
	P(1,3)=Q12
	P(1,4)=Q4
	P(2,3)=Q2
	P(2,4)=Q23
	P(3,4)=Q3
	DO 1 J=2,4
	DO 1 I=1,J-1
	K(I,J)=(M(I)**2+M(J)**2-P(I,J))/M(I)/M(J)
	IF (I.EQ.3) K(I,J)=K(I,J)-M(I)/M(J)
	IF (J.EQ.3) K(I,J)=K(I,J)-M(J)/M(I)
	R(I,J) =SQE_QUAD(CMPLX(1Q0,0Q0,16),CMPLX(-K(I,J),0Q0,16),
     *	            CMPLX(1Q0,0Q0,16))
	IF( IMAG(R(I,J)).EQ.0Q0 ) THEN
	   RS(I,J)=SQE_QUAD(CMPLX(1Q0,0Q0,16),CMPLX(-K(I,J),EPS,16),
     *	               CMPLX(1Q0,0Q0,16))
	ELSE
	   RS(I,J)=R(I,J)
	ENDIF
1	CONTINUE
	SS(1)=RS(1,2)
	SS(2)=RS(2,3)
	SS(3)=RS(3,4)
	SS(4)=RS(1,4)
	AA=K(3,4)/R(2,4)-K(2,3)
	BB=K(1,3)*(1Q0/R(2,4)-R(2,4))+K(1,2)*K(3,4)-K(1,4)*K(2,3)
	CC=K(1,2)*K(1,3)-K(1,3)*K(1,4)*R(2,4)+R(2,4)*K(3,4)-K(2,3)
	DD=K(2,3)-R(2,4)*K(3,4)
	XX(1)=SQE_QUAD(AA,BB,CC+IEPS*DD)
	XX(2)=(CC+IEPS*DD)/AA/XX(1)
	DO 2 I=1,2
	X(I,1)=XX(I)/R(2,4)
	X(I,2)=XX(I)/R(2,4)*R(1,3)
	X(I,3)=XX(I)*R(1,3)
	X(I,4)=XX(I)
2	CONTINUE
	D04QUAD1 = CMPLX(0Q0,0Q0,16)
	DO 3 I=1,2
	D04QUAD1 = D04QUAD1 + (2Q0*I-3Q0)*(
     *		quadcdilog(1Q0+SS(4)*X(I,4))
     *	       -quadcdilog(1Q0+SS(1)*X(I,1))
     *	       +quadcdilog(1Q0+X(I,4)/SS(4))
     *	       -quadcdilog(1Q0+X(I,1)/SS(1))
     *	       +ETA_QUAD(-X(I,4),SS(4))*LOG(1Q0+SS(4)*X(I,4))
     *	       -ETA_QUAD(-X(I,1),SS(1))*LOG(1Q0+SS(1)*X(I,1))
     *	       +ETA_QUAD(-X(I,4),1Q0/SS(4))*LOG(1Q0+X(I,4)/SS(4))
     *	       -ETA_QUAD(-X(I,1),1Q0/SS(1))*LOG(1Q0+X(I,1)/SS(1))
     *	       -quadcdilog(1Q0+X(I,4)*(K(3,4)-IEPS)/(K(1,3)-IEPS))
     *	       +quadcdilog(1Q0+X(I,1)*(K(2,3)-IEPS)/(K(1,3)-IEPS))
     *	       -ETA_QUAD(-X(I,4),(K(3,4)-IEPS)/(K(1,3)-IEPS))
     *	           *LOG(1Q0+X(I,4)*(K(3,4)-IEPS)/(K(1,3)-IEPS))
     *	       +ETA_QUAD(-X(I,1),(K(2,3)-IEPS)/(K(1,3)-IEPS))
     *	           *LOG(1Q0+X(I,1)*(K(2,3)-IEPS)/(K(1,3)-IEPS))   )
	IF (IMAG(R(2,4)).NE.0D0) THEN
	   H=ETA_QUAD(-1Q0/XX(I),R(2,4))
	ELSE
	   H=CMPLX(0Q0,0Q0,16)
	   IF (REAL(R(2,4)).LT.0D0) THEN
	      HH=-1Q0/XX(I)
	      IM1=IMAG(HH)
	      IM2=IMAG(RS(2,4))
              pi = 4.Q0*atan(1.Q0)
	      IF ((IM1.GT.0Q0).AND.(IM2.GT.0Q0)) THEN
	         H=-CMPLX(0Q0,2Q0*PI,16)
	      ENDIF
	      IF ((IM1.LT.0Q0).AND.(IM2.LT.0Q0)) THEN
	         H=+CMPLX(0Q0,2Q0*PI,16)
	      ENDIF
	   ENDIF
	ENDIF
	D04QUAD1 = D04QUAD1 + (2Q0*I-3Q0)*
     *	          H*( LOG( (K(1,2)-R(2,4)*K(1,4)
     *			  +XX(I)*(1Q0/R(2,4)-R(2,4)))/DD )
     *		     +LOG(K(1,3)-IEPS) )
3	CONTINUE
	D04QUAD1 = D04QUAD1/M(1)/M(2)/M(3)/M(4)/AA/(XX(1)-XX(2))
	RETURN

C	****** ONLY MM2 AND MM3 ARE ZERO ******
40	CONTINUE
	EPS=1D-17
	IEPS=CMPLX(0Q0,EPS,16)

	M(1)=MM1
	M(2)=10Q0
	M(3)=10Q0
	M(4)=MM4
	P(1,2)=Q1
	P(1,3)=Q12
	P(1,4)=Q4
	P(2,3)=Q2
	P(2,4)=Q23
	P(3,4)=Q3
	DO 4 J=2,4
	DO 4 I=1,J-1
	K(I,J)=(M(I)**2+M(J)**2-P(I,J))/M(I)/M(J)
	IF (I.EQ.2) K(I,J)=K(I,J)-M(I)/M(J)
	IF (J.EQ.2) K(I,J)=K(I,J)-M(J)/M(I)
	IF (I.EQ.3) K(I,J)=K(I,J)-M(I)/M(J)
	IF (J.EQ.3) K(I,J)=K(I,J)-M(J)/M(I)
	R(I,J) =SQE_QUAD(CMPLX(1Q0,0Q0,16),CMPLX(-K(I,J),0Q0,16),
     *	            CMPLX(1Q0,0Q0,16))
	IF( IMAG(R(I,J)).EQ.0Q0 ) THEN
	   RS(I,J)=SQE_QUAD(CMPLX(1Q0,0Q0,16),CMPLX(-K(I,J),EPS,16),
     *	               CMPLX(1Q0,0Q0,16))
	ELSE
	   RS(I,J)=R(I,J)
	ENDIF
4	CONTINUE
	SS(1)=RS(1,2)
	SS(2)=RS(2,3)
	SS(3)=RS(3,4)
	SS(4)=RS(1,4)
	AA=K(2,4)*K(3,4)-K(2,3)
	BB=K(1,3)*K(2,4)+K(1,2)*K(3,4)-K(1,4)*K(2,3)
	CC=K(1,2)*K(1,3)-K(2,3)
	DD=K(2,3)
	XX(1)=SQE_QUAD(AA,BB,CC+IEPS*DD)
	XX(2)=(CC+IEPS*DD)/AA/XX(1)
	DO 5 I=1,2
	X(I,1)=XX(I)/R(2,4)
	X(I,2)=XX(I)/R(2,4)*R(1,3)
	X(I,3)=XX(I)*R(1,3)
	X(I,4)=XX(I)
5	CONTINUE
	D04QUAD1 = CMPLX(0Q0,0Q0,16)
	DO 6 I=1,2
	D04QUAD1 = D04QUAD1 + (2Q0*I-3Q0)*(
     *		quadcdilog(1Q0+SS(4)*X(I,4))
     *	       +quadcdilog(1Q0+X(I,4)/SS(4))
     *	       +ETA_QUAD(-X(I,4),SS(4))*LOG(1Q0+SS(4)*X(I,4))
     *	       +ETA_QUAD(-X(I,4),1Q0/SS(4))*LOG(1Q0+X(I,4)/SS(4))
     *	       -quadcdilog(1Q0+XX(I)*(K(3,4)-IEPS)/(K(1,3)-IEPS))
     *	       -quadcdilog(1Q0+XX(I)*(K(2,4)-IEPS)/(K(1,2)-IEPS))
     *	       -ETA_QUAD(-XX(I),(K(3,4)-IEPS)/(K(1,3)-IEPS))
     *	           *LOG(1Q0+XX(I)*(K(3,4)-IEPS)/(K(1,3)-IEPS))
     *	       -ETA_QUAD(-XX(I),(K(2,4)-IEPS)/(K(1,2)-IEPS))
     *	           *LOG(1Q0+XX(I)*(K(2,4)-IEPS)/(K(1,2)-IEPS))
     *	       +LOG(-XX(I))*( LOG(K(1,2)-IEPS)
     *			     +LOG(K(1,3)-IEPS)-LOG(K(2,3)-IEPS) ) )
6	CONTINUE
	D04QUAD1 = D04QUAD1/M(1)/M(2)/M(3)/M(4)/AA/(XX(1)-XX(2))

	RETURN


	END
