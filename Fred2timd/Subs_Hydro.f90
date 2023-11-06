     
! ***************************************************************
! *                                                             *
! *  Read the data of gravity center, mass matrix and restoring *
! *  matrix.                                                    *
! *                                                             *
! ***************************************************************
! 
     SUBROUTINE MATCRS

	 USE Platform_mod
     IMPLICIT NONE  

	 INTEGER I,J
     REAL*8 IXX,IYY,IZZ
!
! =======================================================
        
      READ(20,*) 
      READ(20,*)  XCM0(1),XCM0(2),XCM0(3)
      READ(20,*) 
      DO I=1,6
       READ(20,120) (MATX(I,J), J=1, 6)
      ENDDO
      
       READ(20,*) 
      DO I=1,6
       READ(20,120) (CRS(I,J), J=1, 6)
      ENDDO
      READ(20,*) 
      DO I=1,6
       READ(20,120) (LRS(I,J), J=1, 6)
      ENDDO
      READ(20,*)
      DO I=1,6
       READ(20,120) (ZDAMP(I,J),J=1, 6)
      ENDDO
      BLNR=0.D0
      KSPR=0.D0
           
      BLNR(1,1)= 2.23E+05
      BLNR(2,2)= 2.25E+05
      BLNR(3,3)= 3.59E+05
      BLNR(4,4)= 1.61E+10
      BLNR(5,5)= 1.56E+10
      BLNR(6,6)= 1.77E+10
      BLNR(1,5)=-3.72E+07
      BLNR(5,1)=-1.67E+06
      BLNR(2,4)= 3.67E+07
      BLNR(4,2)= 1.68E+06
      BLNR(3,5)=-1.09E+08
      BLNR(6,2)= 5.21E+03
      BLNR(6,4)=-2.09E+08 
        
      BYC= 48851720.9806301D0  !8.371183757269338e+007
  
120	  FORMAT(6(2x,E12.5))          
      RETURN        
      END      
      

   !-------------------------------------------------------------------------------------------------
   ! Buffer function for external force in the beginning of motion...
   !-------------------------------------------------------------------------------------------------       

      REAL*8 FUNCTION BUFR(T,TP)
      
      IMPLICIT NONE
      REAL*8 T,TP,TM,PI
      
      TM=2.D0*TP
      PI=4.D0*DATAN(1.D0)
      
      IF (T.LT.TM) THEN
        BUFR=0.5D0*(1.D0-DCOS(PI*T/TM))  
      ELSE
        BUFR=1.D0
      ENDIF
      
      RETURN
      END    
 

      
! ==================================================================================
!   SOLVER for linear algibric EQUATIONS 
!
      SUBROUTINE SOLVER(IP,A,B,NA,NB,N,NRHS,NSYS)
	  IMPLICIT    NONE

! 
	  INTEGER  N1,N,LR,IBT,IS,ITBP,ITB,LRI,NA,NB,NR,NRHS,IP,NSYS
        REAL*8 A(NA,NA,NSYS),B(NA,NB,NSYS),ALRLR,ALRLRI,BLRNR,BT,ASUM(NA+10)   
!
!   MATRIX FACTORISATION
!
  
      N1=N-1
      IF(N1.NE.0)   GOTO 1
      IF(DABS(A(N,N,IP)).LT.1.E-6)  WRITE(6,9999)N,N
      GOTO 2
    1 DO 1000 LR=1,N1  
      ALRLR=A(LR,LR,IP)
      IF(DABS(ALRLR).LT.1.E-6)  WRITE(6,9999)LR,LR
      IS=LR+1
      DO 10 ITBP=IS,N
   10 A(ITBP,LR,IP)=A(ITBP,LR,IP)/ALRLR
      DO 100 LRI=IS,N
      ALRLRI=A(LR,LRI,IP)
      DO 100 ITB=IS,N
  100 A(ITB,LRI,IP)=A(ITB,LRI,IP)-A(ITB,LR,IP)*ALRLRI
!
!   RHS FACTORISATION
!
      DO 1000 NR=1,NRHS
      BLRNR=B(LR,NR,IP)
      DO 1000 ITB=IS,N
 1000 B(ITB,NR,IP)=B(ITB,NR,IP)-A(ITB,LR,IP)*BLRNR
!
!   BACK SUBSTITUTION
!
    2 CONTINUE
      DO 2000 NR=1,NRHS
      B(N,NR,IP)=B(N,NR,IP)/A(N,N,IP)
      IF(N1.EQ.0)GOTO 2000
      BT=B(N,NR,IP)
      DO 20 ITB=1,N
   20 ASUM(ITB)=A(ITB,N,IP)*BT
      DO 200 IBT=N1,1,-1
      B(IBT,NR,IP)=(B(IBT,NR,IP)-ASUM(IBT))/A(IBT,IBT,IP)
      BT=B(IBT,NR,IP)
      DO 200 ITB=1,IBT
  200 ASUM(ITB)=ASUM(ITB)+A(ITB,IBT,IP)*BT
 2000 CONTINUE
 9999 FORMAT(3X,' ** WARNING **',2I6,'TH DIAG. COEF. ZERO IN SOLVE')
      RETURN
      END            
      
      
      
! 
!
!
!  *************************************************
!  *    The subroutine computes wave number        *
!  *  by wave frequency and water depth.           *
!  *    For infinite water depth, give negative H. *
!  *************************************************

     SUBROUTINE WAVECK(SIGMA,H,WK)
     IMPLICIT  NONE

	 REAL*8 SIGMA,H,WK,B,G,A,Y,C
!
!  H: WATER DEPTH    SIGMA: WAVE FREQUENCY
!  C: WAVE CELERITY  Wk: wave number
!
     DATA G/9.807D0/
!
     IF( SIGMA .LE. 0.0D0 )  THEN
	    Print *,' IN WAVECK:  W=',SIGMA
		STOP
	
     ELSE 

       IF( H .GT. 0.0D0) THEN
           B = G * H
           Y = SIGMA * SIGMA *H/G
  12       A=1.0D0/(1.0D0+Y*(0.66667D0+Y*(0.35550D0+     &
            Y*(0.16084D0+Y*(0.063201D0+Y*(0.02174D0+     &
            Y*(0.00654D0+Y*(0.00171D0+Y*(0.00039D0+Y*0.00011D0)))))))))
  22       C=SQRT(B/(Y+A))
           WK=SIGMA/C
           ELSE IF( H .LE. 0.0D0) THEN
           WK=SIGMA**2/G
       END IF
       RETURN
!
	 END IF

      RETURN
      END


!
!  *******************************************************
!  *                                                     *
!  *    The subroutine computes the real and imaginary   *
!  *  roots of dispersion equation by wave frequency     *
!  *  W and water depth H.                               *
!  *                                                     *
!  *******************************************************
!
!
      SUBROUTINE DISPERS(WAVENO,MNO,W,H)
!
!CC   EVALUATION OF THE ROOTS OF THE FOLLOWING EQUATIONS
!CC   BY NEWTON-RAPHSON METHOD,RESULTS ARE GIVEN IN ARRAY WAVENO
!CC   MNO ARE THE NUMBER OF ROOTS REQUIRED, FIRST ROOT IS FROM EQN. (I)
!CC   THE REST ARE THE FIRST +VE (MNO-1) ROOTS OF (II)
!CC   I) W*W/G = K TANH( KH )
!CC   II) -W*W/G = M TAN( MH )
!

     IMPLICIT NONE
	 INTEGER,INTENT(IN)::MNO
	 REAL*8, INTENT(IN)::W,H
	 REAL*8, INTENT(OUT)::WAVENO(800)
!
	 INTEGER I,M,MM,IFLAG
     REAL*8 UPLIM,LOLIM,G,PI
     REAL*8 WWH,FUN,DFUN,TRIAL,EXX,EXXOR,CHECK


     DATA G,PI/9.807d0,3.141592653589793d0/
!
	 IF(MNO .GT. 799) THEN
	   Print *,' MNO=',MNO,' To enlarge WAVENO(*)'
	   STOP
	 END IF

!
      DO 1 I=1,9
    1 WAVENO(I)=0.D0
      WWH=W*W*H
!
!CC   CALCULATION OF WAVE NUMBER (ROOT OF EQN. (I))
!
      TRIAL=WWH/G
      M=0
      IF (TRIAL.GT.1.D1) GO TO 20
      EXX=0.D0
      IFLAG=0
   10 FUN=G*TRIAL - WWH/DTANH(TRIAL)
      DFUN=G + WWH/DSINH(TRIAL)/DSINH(TRIAL)
      TRIAL=TRIAL - FUN/DFUN
      EXXOR=DABS(TRIAL - EXX)
      IF (EXXOR.LE.1.0D-10) GO TO 20
      EXX=TRIAL
      GO TO 10
   20 MM=M + 1
      WAVENO(MM)=TRIAL/H
      CHECK=DABS(W*W/G - WAVENO(MM)*DTANH(TRIAL))
      IF (CHECK.GT.1.0D-5) GO TO 999
      IF (MNO.LE.1) RETURN
!
!CC   CALCULATION OF FIRST +VE (MNO-1) ROOTS OF EQN. (II)
!
      M=1
      IFLAG=1
      EXX=0.D0
      IF (WWH.LE.2.25D2) GO TO 120
      GO TO 110
  100 M=MM
      EXX=0.D0
      IF (MM.EQ.MNO) GO TO 9999
      IF (IFLAG.EQ.2) GO TO 120
  110 TRIAL=(DBLE(FLOAT(M)) - 1.D-1)*PI
      GO TO 140
  120 IFLAG=2
      TRIAL=(DBLE(FLOAT(M)) - 5.D-1)*PI + 1.D-1
  140 IF(IFLAG .EQ. 1)GO TO 160
      IF(IFLAG .EQ. 2)GO TO 170
  150 TRIAL=TRIAL - FUN/DFUN
      EXXOR=DABS(TRIAL - EXX)
      IF (EXXOR.LE.1.D-10) GO TO 180
      EXX=TRIAL
      IF(IFLAG .EQ. 2)GO TO 170
  160 FUN=G*TRIAL + WWH/DTAN(TRIAL)
      DFUN=G - WWH/DSIN(TRIAL)/DSIN(TRIAL)
      GO TO 150
  170 FUN=WWH/(G*TRIAL) + DTAN(TRIAL)
      DFUN=-WWH/(G*TRIAL*TRIAL) + 1.D0/DCOS(TRIAL)/DCOS(TRIAL)
      GO TO 150
  180 UPLIM=(DBLE(FLOAT(M)) + 5.D-1)*PI
      LOLIM=UPLIM - PI
      IF ((TRIAL.GT.UPLIM).OR.(TRIAL.LT.LOLIM)) GO TO 190
      MM=M + 1
      WAVENO(MM)=TRIAL/H
      CHECK=DABS(W*W/G + WAVENO(MM)*DTAN(TRIAL))
      IF (CHECK.GT.1.0D-5) GO TO 999
      GO TO 100
  190 IF (IFLAG.EQ.1) GO TO 120
      WRITE(6,200)M
  200 FORMAT('  **ERROR**  OCCURS AT M',I3)
      STOP
  999 WRITE(6,1000)CHECK
      GO TO 190
 1000 FORMAT('  CHECK=',D11.4)
 9999 CONTINUE

      RETURN
      END
      
      
      