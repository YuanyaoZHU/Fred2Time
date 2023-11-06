
C/************************************************************************/
C/*                                                                      */
C/* ALGORITHM 412 - COMPUTE THE CONVOLUTION OF RETADATION FUNCTION       */
C/*                                                                      */
C/* Yingyi LIU 17 Dec 2013.                                              */
C/* With modifications suggested by Yingyi LIU,                          */
C/* Research Institute of Applied Mechanics, Kyshu-U, Japan,             */
C/* 8 December 2013.                                                     */
C/*                                                                      */
C/************************************************************************/
  
      SUBROUTINE CONVRETAR1(L,TAG,TSP,NSTP,RFK,DX,FTAR,RK)  ! 迟滞函数时间步细分为二
                                                  ! 此程序计算结果觼ELiu(2012) 一致,特别是振幅一致
      IMPLICIT NONE 
      INTEGER I,J,K,L,NSTP,TAG
      REAL*8 TSP,RFK(6,6,0:2*NSTP),DX(6,0:NSTP),FTAR(6),RK(6,4,0:NSTP)
      REAL*8 RTF,VCT
      
      
      IF (TAG.EQ.1) THEN     
        DO 100 I=1,6
         FTAR(I)=0.D0
        DO 100 J=1,6
         DO K=0,2*L
          IF (MOD(K,2).EQ.0) THEN
            VCT= DX(J,K/2)
          ELSE
            VCT=( DX(J,(K-1)/2) + DX(J,(K+1)/2) )/2.D0            
          ENDIF    
          RTF= RFK(I,J,2*L-K)
          FTAR(I)= FTAR(I)-RTF*VCT*TSP/2.D0
         ENDDO

100     CONTINUE   
      
      ELSEIF (TAG.EQ.2) THEN    
        DO 200 I=1,6
         FTAR(I)=0.D0          
        DO 200 J=1,6
         DO K=0,2*L
          IF (MOD(K,2).EQ.0) THEN
            VCT= DX(J,K/2) + RK(J,1,K/2)/2.D0
          ELSE
            VCT=( DX(J,(K-1)/2)   + DX(J,(K+1)/2)   )/2.D0+ 
     1          ( RK(J,1,(K-1)/2) + RK(J,1,(K+1)/2) )/4.D0           
          ENDIF    
          RTF= RFK(I,J,2*L+1-K)
          FTAR(I)= FTAR(I)-RTF*VCT*TSP/2.D0
         ENDDO
200     CONTINUE       
      
      ELSEIF (TAG.EQ.3) THEN    
        DO 300 I=1,6
         FTAR(I)=0.D0          
        DO 300 J=1,6
         DO K=0,2*L
          IF (MOD(K,2).EQ.0) THEN
            VCT= DX(J,K/2) + RK(J,2,K/2)/2.D0
          ELSE
            VCT=( DX(J,(K-1)/2)   + DX(J,(K+1)/2)   )/2.D0+ 
     1          ( RK(J,2,(K-1)/2) + RK(J,2,(K+1)/2) )/4.D0            
          ENDIF    
          RTF= RFK(I,J,2*L+1-K)
          FTAR(I)= FTAR(I)-RTF*VCT*TSP/2.D0
         ENDDO
300     CONTINUE         
        
      ELSEIF (TAG.EQ.4) THEN    
        DO 400 I=1,6
         FTAR(I)=0.D0          
        DO 400 J=1,6
         DO K=0,2*L
          IF (MOD(K,2).EQ.0) THEN
            VCT= DX(J,K/2) + RK(J,3,K/2)/2.D0
          ELSE
            VCT=( DX(J,(K-1)/2) + RK(J,3,(K-1)/2)/2.D0+ 
     1            DX(J,(K+1)/2) + RK(J,3,(K+1)/2)/2.D0 )/2.D0 
          ENDIF    
          RTF= RFK(I,J,2*L+2-K)
          FTAR(I)= FTAR(I)-RTF*VCT*TSP/2.D0
         ENDDO         
400     CONTINUE 
        
      ENDIF
      
      RETURN 
      END
      
C/************************************************************************/
C/*                                                                      */
C/* ALGORITHM 412 - COMPUTE THE CONVOLUTION OF RETADATION FUNCTION       */
C/*                                                                      */
C/* Yingyi LIU 17 Dec 2013.                                              */
C/* With modifications suggested by Yingyi LIU,                          */
C/* Research Institute of Applied Mechanics, Kyshu-U, Japan,             */
C/* 8 December 2013.                                                     */
C/*                                                                      */
C/************************************************************************/
  
      SUBROUTINE CONVRETAR(L,TAG,TSP,NSTP,RFK,DX,FTAR,RK)  ! 迟滞函数时间步细分为二

      IMPLICIT NONE
      INTEGER I,J,K,L,NSTP,TAG
      REAL*8 TSP,RFK(6,6,0:2*NSTP),DX(6,0:NSTP),FTAR(6),RK(6,4,0:NSTP)
      REAL*8 RTF,VCT
      
      
      IF (TAG.EQ.1) THEN     
        DO 100 I=1,6
         FTAR(I)=0.D0
        DO 100 J=1,6
         DO K=0,L             
          VCT=DX(J,K)
          RTF=RFK(I,J,2*L-2*K)
          FTAR(I)=FTAR(I)-RTF*VCT*TSP
         ENDDO

100     CONTINUE   
      
      ELSEIF (TAG.EQ.2) THEN    
        DO 200 I=1,6
         FTAR(I)=0.D0          
        DO 200 J=1,6
         DO K=0,L  
          VCT=DX(J,K)+RK(J,1,K)/2.D0
          RTF=RFK(I,J,2*L+1-2*K)      
          FTAR(I)=FTAR(I)-RTF*VCT*TSP       
         ENDDO
200     CONTINUE       
      
      ELSEIF (TAG.EQ.3) THEN    
        DO 300 I=1,6
         FTAR(I)=0.D0
        DO 300 J=1,6
         DO K=0,L
          VCT=DX(J,K)+RK(J,2,K)/2.D0
          RTF=RFK(I,J,2*L+1-2*K)
          FTAR(I)=FTAR(I)-RTF*VCT*TSP  ! should be TSP/2.d0
         ENDDO
300     CONTINUE
        
      ELSEIF (TAG.EQ.4) THEN
        DO 400 I=1,6
         FTAR(I)=0.D0
        DO 400 J=1,6
         DO K=0,L  
          VCT=DX(J,K)+RK(J,3,K)/2.D0   ! 
          RTF=RFK(I,J,2*L+2-2*K)
          FTAR(I)=FTAR(I)-RTF*VCT*TSP       
         ENDDO
400     CONTINUE 
        
      ENDIF
      
      RETURN 
      END


!C/************************************************************************/
!C/*                                                                      */
!C/* ALGORITHM 412 - Calculation of survival function of Kij              */
!C/*                                                                      */
!C/* Yingyi LIU 16 Dec 2013.                                              */
!C/* With modifications suggested by Yingyi LIU,                          */
!C/* Research Institute of Applied Mechanics, Kyshu-U, Japan,             */
!C/* 8 December 2013.                                                     */
!C/*                                                                      */
!C/************************************************************************/
      SUBROUTINE SUVQ_KIJ(WLB,BIJUB,T,MODE,VAL)

      IMPLICIT NONE
      INTEGER MODE
      REAL*8 WLB,BIJUB,T,VAL
      REAL*8 POW,CK,PI,CI,SI
      REAL*8 A1,A2,A3,B1,B2,B3,A,B,C,G
      REAL*8 WSQ,WSQ2,WSQ3,TSQ,TSQ2,TSQ3,TA
      
      PI=4.D0*DATAN(1.D0)
      IF (MODE.EQ.3.OR.MODE.EQ.5) THEN
       POW=7
      ELSE
       POW=3
      ENDIF
 
      CK=BIJUB*WLB**POW
      
      WSQ=WLB*WLB
      WSQ2=WSQ*WSQ
      WSQ3=WSQ*WSQ2
      TSQ=T*T
      TSQ2=TSQ*TSQ
      TSQ3=TSQ*TSQ2
      TA=WLB*T
      
      IF (MODE.EQ.3.OR.MODE.EQ.5) THEN
      
       IF (DABS(TA).LT.1.E-8) THEN
           
        A1=1.D0/3.D0/WSQ3
        A=A1*DCOS(TA)
        VAL=A*CK/PI
        
       ELSE

        A1=1.D0/3.D0/WSQ3/TSQ3
        A2=1.D0/60.D0/WSQ2/TSQ2
        A3=1.D0/360.D0/WSQ/TSQ
        B1=TA*A1/5.D0
        B2=TA*A2/3.D0
        B3=TA*A3 
       
        A=(A1-A2+A3)*DCOS(TA)
        B=(-B1+B2-B3)*DSIN(TA)
        CALL CISIA(TA,CI,SI)
        G=1.D0/360.D0*CI
       
        C=A+B+G
        IF (DABS(C).LT.1.E-20) THEN
         VAL=0.D0
        ELSE
         VAL=C*CK/PI*TSQ3
        ENDIF
        
       ENDIF
       
      ELSE
           
       A=DCOS(TA)/WSQ
       B=DSIN(TA)/WLB*T   
       CALL CISIA(TA,CI,SI)
       G=TSQ*CI
              
       VAL=(A-B+G)*CK/PI
       
      ENDIF
      
      
      RETURN 
      END
      

C/************************************************************************/
C/*                                                                      */
C/* ALGORITHM 353 - FILON QUADRATURE                                     */
C/*                                                                      */
C/* Yingyi LIU 25 Nov 2013.                                              */
C/* With modifications suggested by Yingyi LIU,                          */
C/* Research Institute of Applied Mechanics, Kyshu-U, Japan,             */
C/* 8 December 2013.                                                     */
C/*                                                                      */
C/************************************************************************/

C/* FILON evaluates the integrals                                        */

C/*      INF                         INF                                 */
C/* C = S F(X) cos (T X) dX,     S = S F(X) sin (T X) dX                 */
C/*      0                            0                                  */

C/* using the Filon quadrature algorithm.                                 */
C/* The user may request an evaluation of C only, S only, or both C and S. */

C/FILON: PROCEDURE (FX,LEN,T,LB,UB,OPT,EPS,VAL);
C/*        FX  = name of a function to evaluate F(X), supplied by the user, */
C/*              or values in discrete vector form.  
C/*        LEN = length of the vector FX.                                   */
C/*        T   = parameter in the formula.                                  */
C/*        LB  = lower bound of the integral.                               */
C/*        UB  = upper bound of the integral.                               */
C/*        OPT = 1 if the cosine integral is to be computed                 */
C/*            = 0 if the sine integral is to be computed.                  */
C/*        EPS = approximation error for the algorithm.                     */
C/*        VAL = result of the filon integration.                           */

      SUBROUTINE FILON_KIJ(FX,LEN,T,LB,UB,OPT,VAL)
      
      INTEGER,INTENT(IN):: LEN,OPT
      REAL*8,INTENT(IN):: FX(0:LEN),LB,UB,T
      REAL*8,INTENT(OUT):: VAL      
      REAL*8 H,TP,TSQ,S1SQ,S1,C1,P
      REAL*8 A,B,B1,B2,B3,B4,B5,G
      REAL*8 X(0:LEN),C0,C2N,C2N_1,PI
      INTEGER I,J,N,TAG
      
      DATA PI/3.141592653589793D0/
C
C  COMPUTE THE APPROXIMATION OF ERROR.
C  IF FX REPRESENTS DAMPING VECTOR, THE FOLLOWING TYPE IS CALCULATED..
C              INF 
C  LIJ(T)=2/PI*S BIJ(W)/W*SIN(W*T) dW               OPT=0
C               0      
C  IF FX REPRESENTS ADDED MASS VECTOR, THE FOLLOWING TYPE IS CALCULATED..
C              INF 
C  LIJ(T)=2/PI*S (AIJ(W)-AIJ(INF))*COS(W*T) dW      OPT=1    
C               0            
      
      IF (MOD(LEN,2).NE.0) THEN
        PRINT*, "SIZE OF INPUT VECTOR MUST BE EVEN"  
      ENDIF
      
      N=LEN/2
      H=(UB-LB)/LEN
      TP=T*H
      TSQ=TP*TP
      
      X(0)=LB
      X(LEN)=UB
      DO 50 I=1,LEN-1
       X(I)=LB+I*H
50    CONTINUE   
  
      
      IF (TP.LT.0.1D0) THEN
C
C  74 IS THE POWER SERIES FOR SMALL T, 75 IS THE CLOSED FORM USED WITH
C  LARGER VALUES OF T.
C  THE POWER SERIES (WITH 'TN' = TP**N)
C  A = (2/45)*T3 - (2/315)*T5 + (2/4725)*T7
C  B = (2/3) + (2/15)*T2 - (4/105)*T4 + (2/567)*T6 - (4/22275)*T8
C  G = (4/3) - (2/15)*T2 + (1/210)*T4 - (1/11340)*T6
C  THE NEXT TERM IN G IS TOO SMALL.  IT IS (1/997920)*T8.
C
74    A = TP * TSQ * ( 1.0E+00 - TSQ * ( 1.0E+00 - TSQ / 15.0E+00 )
     &  / 7.0E+00 ) / 22.5E+00
      B1 = 2.0E+00 / 3.0E+00      
      B2 = B1 * TSQ * 0.2E+00
      B3 = B2 * TSQ * 2.0E+00 / 7.0E+00
      B4 = B3 * TSQ / 10.8E+00
      B5 = B4 * TSQ * 14.0E+00 / 275.0E+00
      B = B1 + B2 - B3 + B4 - B5
      G = 2.0E+00 * B1 - B2 + B3 / 8.0E+00 - B4 / 40.0E+00
      

C
C  G = 2*B1 - B2 + B3/8 - B4/40 + 5*B5/896.  IF YOU WANT THE T8
C  TERM INCLUDED IN G.
C
      ELSE
C
C  CLOSED FORM OF THE COEFFICIENTS, WHERE AGAIN 'TN' MEANS TP**N.
C  A = 1/TP + COS(TP)*SIN(TP)/T2 - 2*(SIN(TP))**2/T3
C  B = 2*((1+(COS(TP))**2)/T2 - 2*SIN(TP)*COS(TP)/T3).
C  G = 4*( SIN(TP)/T3 - COS(TP)/T2 ).
C      
      S1 = DSIN ( TP )
      C1 = DCOS ( TP ) 
      P = S1 * C1
      S1SQ = S1 * S1
      
      A = ((( -2.0E+00 * S1SQ / TP ) + P ) / TP + 1.0E+00 ) / TP
      B = 2.0E+00 * ( ( -2.0E+00 * P / TP ) + 2.0E+00 - S1SQ ) / TSQ
      G = 4.0E+00 * ( S1 / TP - C1 ) / TSQ        
      

      ENDIF      
      
C
C  COMPUTE THE VALUE OF INTEGRAL.
C  FOR COSINE TRANSFORM...
C  C0 = F(X2N)*SIN(X2N*T)-F(0)*SIN(X0*T)
C  C2N = -0.5*(F(X2N)*SIN(X2N*T)-F(0)*SIN(X0*T))
C        +THE SUM FROM I=0,...,N of [ F(2I)*COS(X2I*T) ]
C  C2N_1 = THE SUM FROM I=1,...,N of [ F(2I-1)*COS(X2I_1*T) ]
C
C  FOR SINE TRANSFORM...
C  C0 = F(0)*SIN(X0*T)-F(X2N)*SIN(X2N*T)
C  C2N = -0.5*(F(X2N)*SIN(X2N*T)-F(0)*SIN(X0*T))
C        +THE SUM FROM I=0,...,N of [ F(2I)*SIN(X2I*T) ]
C  C2N_1 = THE SUM FROM I=1,...,N of [ F(2I-1)*SIN(X2I_1*T) ]    
C
      IF (OPT.EQ.1) THEN
        
        C0=FX(LEN)*DSIN(X(LEN)*T)-FX(0)*DSIN(X(0)*T)
          
        C2N=-0.5D0*(FX(0)*DCOS(LB*T)+FX(LEN)*DCOS(UB*T))
        DO 100 I=0,N
          C2N=C2N+FX(2*I)*DCOS(X(2*I)*T)
100     CONTINUE  
        
        C2N_1=0.D0
        DO 150 I=1,N
          C2N_1=C2N_1+FX(2*I-1)*DCOS(X(2*I-1)*T)
150     CONTINUE  
        
        VAL=H*(A*C0+B*C2N+G*C2N_1)
          
      ELSE
          
        C0=FX(0)*DCOS(X(0)*T)-FX(LEN)*DCOS(X(LEN)*T)
          
        C2N=-0.5D0*(FX(0)*DSIN(LB*T)+FX(LEN)*DSIN(UB*T))
        DO 200 I=0,N
          C2N=C2N+FX(2*I)*DSIN(X(2*I)*T)
200     CONTINUE  
        
        C2N_1=0.D0
        DO 250 I=1,N
          C2N_1=C2N_1+FX(2*I-1)*DSIN(X(2*I-1)*T)
250     CONTINUE  
        
        VAL=H*(A*C0+B*C2N+G*C2N_1)
        
      ENDIF

          
      RETURN
      END      
      
C/************************************************************************/
C/*                                                                      */
C/* ALGORITHM 353 - FILON QUADRATURE                                     */
C/*                                                                      */
C/* Yingyi LIU 25 Nov 2013.                                              */
C/* With modifications suggested by Yingyi LIU,                          */
C/* Research Institute of Applied Mechanics, Kyshu-U, Japan,             */
C/* 8 December 2013.                                                     */
C/*                                                                      */
C/************************************************************************/

C/* FILON evaluates the integrals                                        */

C/*      INF                         INF                                 */
C/* C = S F(X) cos (T X) dX,     S = S F(X) sin (T X) dX                 */
C/*      0                            0                                  */

C/* using the Filon quadrature algorithm.                                 */
C/* The user may request an evaluation of C only, S only, or both C and S. */

C/FILON: PROCEDURE (FX,LEN,T,LB,UB,OPT,EPS,VAL);
C/*        FX  = name of a function to evaluate F(X), supplied by the user, */
C/*              or values in discrete vector form.  
C/*        LEN = length of the vector FX.                                   */
C/*        T   = parameter in the formula.                                  */
C/*        LB  = lower bound of the integral.                               */
C/*        UB  = upper bound of the integral.                               */
C/*        OPT = 1 if the cosine integral is to be computed                 */
C/*            = 0 if the sine integral is to be computed.                  */
C/*        EPS = approximation error for the algorithm.                     */
C/*        VAL = result of the filon integration.                           */

      SUBROUTINE FILON_LIJ(FX,LEN,T,LB,UB,OPT,EPS,VAL)
      
      INTEGER,INTENT(IN):: LEN,OPT
      REAL*8,INTENT(IN):: FX(0:LEN),LB,UB,T
      REAL*8,INTENT(OUT):: EPS,VAL      
      REAL*8 H,TP,TSQ,S1SQ,S1,C1,P,TA,FAC
      REAL*8 A,B,B1,B2,B3,B4,B5,G,CI,SI
      REAL*8 X(0:LEN),C0,C2N,C2N_1,PI
      INTEGER I,J,N
      
      DATA PI/3.141592653589793D0/
C
C  COMPUTE THE APPROXIMATION OF ERROR.
C  IF FX REPRESENTS DAMPING VECTOR, THE FOLLOWING TYPE IS CALCULATED..
C              INF 
C  LIJ(T)=2/PI*S BIJ(W)/W*SIN(W*T) dW               OPT=0
C               0      
C  IF FX REPRESENTS ADDED MASS VECTOR, THE FOLLOWING TYPE IS CALCULATED..
C              INF 
C  LIJ(T)=2/PI*S (AIJ(W)-AIJ(INF))*COS(W*T) dW      OPT=1    
C               0            
      
      IF (MOD(LEN,2).NE.0) THEN
        PRINT*, "SIZE OF INPUT VECTOR MUST BE EVEN"  
      ENDIF
      
      N=LEN/2
      H=(UB-LB)/LEN
      TP=T*H
      TSQ=TP*TP
      
      X(0)=LB
      X(LEN)=UB
      DO 50 I=1,LEN
       X(I)=LB+I*H
50    CONTINUE               
      
      IF (TP.LT.0.1D0) THEN
C
C  74 IS THE POWER SERIES FOR SMALL T, 75 IS THE CLOSED FORM USED WITH
C  LARGER VALUES OF T.
C  THE POWER SERIES (WITH 'TN' = TP**N)
C  A = (2/45)*T3 - (2/315)*T5 + (2/4725)*T7
C  B = (2/3) + (2/15)*T2 - (4/105)*T4 + (2/567)*T6 - (4/22275)*T8
C  G = (4/3) - (2/15)*T2 + (1/210)*T4 - (1/11340)*T6
C  THE NEXT TERM IN G IS TOO SMALL.  IT IS (1/997920)*T8.
C
74    A = TP * TSQ * ( 1.0E+00 - TSQ * ( 1.0E+00 - TSQ / 15.0E+00 )
     &  / 7.0E+00 ) / 22.5E+00
      B1 = 2.0E+00 / 3.0E+00      
      B2 = B1 * TSQ * 0.2E+00
      B3 = B2 * TSQ * 2.0E+00 / 7.0E+00
      B4 = B3 * TSQ / 10.8E+00
      B5 = B4 * TSQ * 14.0E+00 / 275.0E+00
      B = B1 + B2 - B3 + B4 - B5
      G = 2.0E+00 * B1 - B2 + B3 / 8.0E+00 - B4 / 40.0E+00
C
C  G = 2*B1 - B2 + B3/8 - B4/40 + 5*B5/896.  IF YOU WANT THE T8
C  TERM INCLUDED IN G.
C
      ELSE
C
C  CLOSED FORM OF THE COEFFICIENTS, WHERE AGAIN 'TN' MEANS TP**N.
C  A = 1/TP + COS(TP)*SIN(TP)/T2 - 2*(SIN(TP))**2/T3
C  B = 2*((1+(COS(TP))**2)/T2 - 2*SIN(TP)*COS(TP)/T3).
C  G = 4*( SIN(TP)/T3 - COS(TP)/T2 ).
C      
      S1 = DSIN ( TP )
      C1 = DCOS ( TP )      
      P = S1 * C1
      S1SQ = S1 * S1
      
      A = ((( -2.0E+00 * S1SQ / TP ) + P ) / TP + 1.0E+00 ) / TP
      B = 2.0E+00 * ( ( -2.0E+00 * P / TP ) + 2.0E+00 - S1SQ ) / TSQ
      G = 4.0E+00 * ( S1 / TP - C1 ) / TSQ         
      
      ENDIF
C
C  COMPUTE THE VALUE OF INTEGRAL.
C  FOR COSINE TRANSFORM...
C  C0 = F(X2N)*SIN(X2N*T)-F(0)*SIN(X0*T)
C  C2N = -0.5*(F(X2N)*SIN(X2N*T)-F(0)*SIN(X0*T))
C        +THE SUM FROM I=0,...,N of [ F(2I)*COS(X2I*T) ]
C  C2N_1 = THE SUM FROM I=1,...,N of [ F(2I-1)*COS(X2I_1*T) ]
C
C  FOR SINE TRANSFORM...
C  C0 = F(0)*SIN(X0*T)-F(X2N)*SIN(X2N*T)
C  C2N = -0.5*(F(X2N)*SIN(X2N*T)-F(0)*SIN(X0*T))
C        +THE SUM FROM I=0,...,N of [ F(2I)*SIN(X2I*T) ]
C  C2N_1 = THE SUM FROM I=1,...,N of [ F(2I-1)*SIN(X2I_1*T) ]    
C
      IF (OPT.EQ.1) THEN
        
        C0=FX(LEN)*DSIN(X(LEN)*T)-FX(0)*DSIN(X(0)*T)
          
        C2N=-0.5D0*(FX(0)*DCOS(LB*T)+FX(LEN)*DCOS(UB*T))
        DO 100 I=0,N
          C2N=C2N+FX(2*I)*DCOS(X(2*I)*T)
100     CONTINUE  
        
        C2N_1=0.D0
        DO 150 I=1,N
          C2N_1=C2N_1+FX(2*I-1)*DCOS(X(2*I-1)*T)
150     CONTINUE  
        
        VAL=H*(A*C0+B*C2N+G*C2N_1)
        
      ELSE
          
        C0=FX(0)*DCOS(X(0)*T)-FX(LEN)*DCOS(X(LEN)*T)
          
        C2N=-0.5D0*(FX(0)*DSIN(LB*T)+FX(LEN)*DSIN(UB*T))
        DO 200 I=0,N
          C2N=C2N+FX(2*I)*DSIN(X(2*I)*T)
200     CONTINUE  
        
        C2N_1=0.D0
        DO 250 I=1,N
          C2N_1=C2N_1+FX(2*I-1)*DSIN(X(2*I-1)*T)
250     CONTINUE  
        
        VAL=H*(A*C0+B*C2N+G*C2N_1)
        
      ENDIF
C
C  COMPUTE THE APPROXIMATION OF ERROR.
C  IF FX REPRESENTS DAMPING VECTOR...
C                                                UB 
C  W(T)=-4/(PI**2*UB)*(COS(UB*T)+UB*T*SI(UB*T))*S F(X) dX      
C                                                0      
C  IF FX REPRESENTS ADDED MASS VECTOR...
C                                       UB 
C  W(T)=2/PI*(COS(UB*T)+UB*T*SI(UB*T))*S (F(X)-F(INF)) dX      
C                                       0 
      IF (OPT.EQ.1) THEN  
          
        EPS=0.D0  
        DO 350 I=1,N-1
         EPS=EPS+(FX(I+1)-FX(I))*(X(I+1)-X(I))  
350     CONTINUE 
        TA=UB*T
        CALL CISIA(TA,CI,SI)
        FAC= 2.D0/PI*(DCOS(TA)+TA*SI)
        EPS= FAC*EPS 
        
      ELSE
          
        EPS=0.D0  
        DO 300 I=1,N-1
         EPS=EPS+(FX(I+1)*X(I+1)-FX(I)*X(I))*(X(I+1)-X(I))  
300     CONTINUE 
        TA=UB*T
        CALL CISIA(TA,CI,SI)
        FAC= -4.D0/(PI**2*UB)*(DCOS(TA)+TA*SI)
        EPS= FAC*EPS  
          
      ENDIF          
          
      RETURN
      END       

C
C       ========================================================
C       Purpose: This program computes the cosine and sine 
C                integrals using subroutine CISIA
C       Input :  x  --- Argument of Ci(x) and Si(x)
C       Output:  CI --- Ci(x)
C                SI --- Si(x)
C       Example:
C                    x         Ci(x)          Si(x)
C                 ------------------------------------
C                   0.0     - ?            .00000000
C                   5.0     -.19002975     1.54993124
C                  10.0     -.04545643     1.65834759
C                  20.0      .04441982     1.54824170
C                  30.0     -.03303242     1.56675654
C                  40.0      .01902001     1.58698512
C       ========================================================
C
        SUBROUTINE CISIA(X,CI,SI)
C
C       =============================================
C       Purpose: Compute cosine and sine integrals
C                Si(x) and Ci(x)  ( x ?0 )
C       Input :  x  --- Argument of Ci(x) and Si(x)
C       Output:  CI --- Ci(x)
C                SI --- Si(x)
C       =============================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION BJ(101)
        P2=1.570796326794897D0
        EL=.5772156649015329D0
        EPS=1.0D-15
        X2=X*X
        IF (X.EQ.0.0D0) THEN
           CI=-1.0D+300
           SI=0.0D0
        ELSE IF (X.LE.16.0D0) THEN
           XR=-.25D0*X2
           CI=EL+DLOG(X)+XR
           DO 10 K=2,40
              XR=-.5D0*XR*(K-1)/(K*K*(2*K-1))*X2
              CI=CI+XR
              IF (DABS(XR).LT.DABS(CI)*EPS) GO TO 15
10         CONTINUE
15         XR=X
           SI=X
           DO 20 K=1,40
              XR=-.5D0*XR*(2*K-1)/K/(4*K*K+4*K+1)*X2
              SI=SI+XR
              IF (DABS(XR).LT.DABS(SI)*EPS) RETURN
20         CONTINUE
        ELSE IF (X.LE.32.0D0) THEN
           M=INT(47.2+.82*X)
           XA1=0.0D0
           XA0=1.0D-100
           DO 25 K=M,1,-1
              XA=4.0D0*K*XA0/X-XA1
              BJ(K)=XA
              XA1=XA0
25            XA0=XA
           XS=BJ(1)
           DO 30 K=3,M,2
30            XS=XS+2.0D0*BJ(K)
           BJ(1)=BJ(1)/XS
           DO 35 K=2,M
35            BJ(K)=BJ(K)/XS
           XR=1.0D0
           XG1=BJ(1)
           DO 40 K=2,M
              XR=.25D0*XR*(2.0*K-3.0)**2/((K-1.0)*(2.0*K-1.0)**2)*X
40            XG1=XG1+BJ(K)*XR
           XR=1.0D0
           XG2=BJ(1)
           DO 45 K=2,M
              XR=.25D0*XR*(2.0*K-5.0)**2/((K-1.0)*(2.0*K-3.0)**2)*X
45            XG2=XG2+BJ(K)*XR
           XCS=DCOS(X/2.0D0)
           XSS=DSIN(X/2.0D0)
           CI=EL+DLOG(X)-X*XSS*XG1+2*XCS*XG2-2*XCS*XCS
           SI=X*XCS*XG1+2*XSS*XG2-DSIN(X)
        ELSE
           XR=1.0D0
           XF=1.0D0
           DO 50 K=1,9
              XR=-2.0D0*XR*K*(2*K-1)/X2
50            XF=XF+XR
           XR=1.0D0/X
           XG=XR
           DO 55 K=1,8
              XR=-2.0D0*XR*(2*K+1)*K/X2
55            XG=XG+XR
           CI=XF*DSIN(X)/X-XG*DCOS(X)/X
           SI=P2-XF*DCOS(X)/X-XG*DSIN(X)/X
        ENDIF
        RETURN
        END