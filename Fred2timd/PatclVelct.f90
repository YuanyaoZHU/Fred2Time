!    ----------------------------------------------------------------------------------------------------
!       For derivatives of incident wave potential in regular wave, i.e. incident wave orbital velocity
!    ----------------------------------------------------------------------------------------------------
!
        SUBROUTINE  DINP_RGU(XLC,T,DPOT)
        
	    USE Const_mod
	    USE WaveDyn_mod        
        
	    IMPLICIT    NONE
	  
	    REAL*8,INTENT(IN):: XLC(3),T
        REAL*8,INTENT(OUT)::  DPOT(3)

	    REAL*8 Z,WKX,C2CH,C2SH,S2CH,S2SH,C1,S1,THKDH,THKDV
        
        Z=XLC(3)
        !IF (Z.GT.0.D0) WRITE(*,*) 'Z>0, IN DINP', 'Z=',Z
   
        WKX=WK*(XLC(1)*DCOS(BETA)+XLC(2)*DSIN(BETA))
        C1=DCOS(WKX-W1*T)
        S1=DSIN(WKX-W1*T)

	    IF (H .LT. 0.0d0) THEN
          THKDH= DEXP(WK*Z)
          THKDV= DEXP(WK*Z)
        ELSE
          THKDH= C2SH(WK,Z,H)
          THKDV= S2SH(WK,Z,H)
	    ENDIF

        DPOT(1)=    W1*HS*THKDH*C1*DCOS(BETA)
        DPOT(2)=    W1*HS*THKDH*C1*DSIN(BETA)
        DPOT(3)=    W1*HS*THKDV*S1
        
        RETURN
        END SUBROUTINE DINP_RGU
        
!    -------------------------------------------------------------------------------------------------------
!       For derivatives of incident wave potentialin irregular waves, i.e. incident wave orbital velocity
!    -------------------------------------------------------------------------------------------------------
!
        SUBROUTINE  DINP_IRGU(XLC,T,DPOT)
        
	    USE Const_mod
	    USE WaveDyn_mod        
        USE mt95_mod      
        
	    IMPLICIT    NONE
	  
	    REAL*8,INTENT(IN):: XLC(3),T
        REAL*8,INTENT(OUT):: DPOT(3)
        
        INTEGER M,I,J
        REAL*8 WMIN,WMAX,DW,SW,EPSL,KX,X,Y,Z
        REAL*8 C2SH,S2SH,C1,S1,THKDH,THKDV
        REAL*8,ALLOCATABLE::W(:),WW(:),WI(:),KI(:),AI(:),RA(:)
        REAL*8,EXTERNAL:: WaveNumber
      
        X=XLC(1)
        Y=XLC(2)
        Z=XLC(3)
        
        M=200
        ALLOCATE(W(1:M+1),WW(1:M),WI(1:M),KI(1:M),AI(1:M),RA(M))      
      
        WMIN=0.0D0
        WMAX=3.D0
        DW=(WMAX-WMIN)/M
      
        DO I=0,M
          W(I+1)=WMIN+I*DW
        ENDDO
!      
        call genrand_init()
        call genrand_real1(RA)
      
        DO I=1,M     
          WW(I)=(W(I)+W(I+1))/2.D0
          WI(I)=W(I)+DW*RA(I)
          KI(I)= WaveNumber ( WI(I), G, H )
          CALL JONSWAP_goda( WW(I), Hs, Tp , Sw )
          AI(I)=DSQRT(2.D0*SW*DW)
        ENDDO

        call genrand_real1(RA)       

        DPOT=0.D0
       
        DO 100 I=1,M          
        
          KX=KI(I)*(X*DCOS(BETA)+Y*DSIN(BETA))
          EPSL=2.D0*PI*RA(I)
          C1=DCOS(KX-WI(I)*T+EPSL)
          S1=DSIN(KX-WI(I)*T+EPSL)
          
	      IF (H .LT. 0.0d0) THEN
            THKDH= DEXP(KI(I)*Z)
            THKDV= DEXP(KI(I)*Z)
          ELSE
            THKDH= C2SH(KI(I),Z,H)
            THKDV= S2SH(KI(I),Z,H)
          ENDIF
          
          DPOT(1)= DPOT(1) + WI(I)*AI(I)*THKDH*C1*DCOS(BETA)
          DPOT(2)= DPOT(2) + WI(I)*AI(I)*THKDH*C1*DSIN(BETA)
          DPOT(3)= DPOT(3) + WI(I)*AI(I)*THKDV*S1
          
100     CONTINUE
        
      
        DEALLOCATE(W,WW,WI,KI,AI,RA)
        
        RETURN
        END SUBROUTINE DINP_IRGU
        
!    ----------------------------------------------------------------------------------------
!       For caculation of cosh(k(z+h))/cosh(kh)
!    ----------------------------------------------------------------------------------------
!
        REAL*8 FUNCTION  C2CH(K,Z,H)       
        
	    IMPLICIT    NONE
	  
	    REAL*8,INTENT(IN):: K,Z,H
 
        IF (K*(Z+H).GT.89.4d0) THEN
          C2CH=DEXP(K*Z)
        ELSE
          C2CH=DCOSH(K*(Z+H))/DCOSH(K*H)  
        ENDIF
        
        
        RETURN
        END        
        
!    ----------------------------------------------------------------------------------------
!       For caculation of sinh(k(z+h))/cosh(kh)
!    ----------------------------------------------------------------------------------------
!
        REAL*8 FUNCTION  S2CH(K,Z,H)
        
	    IMPLICIT    NONE
	  
	    REAL*8,INTENT(IN):: K,Z,H
 
        IF (K*(Z+H).GT.89.4d0) THEN
          S2CH=DEXP(K*Z)
        ELSE
          S2CH=DSINH(K*(Z+H))/DCOSH(K*H)
        ENDIF
        
        
        RETURN
        END      
        
!    ----------------------------------------------------------------------------------------
!       For caculation of cosh(k(z+h))/sinh(kh)
!    ----------------------------------------------------------------------------------------
!
        REAL*8 FUNCTION  C2SH(K,Z,H)
        
	    IMPLICIT    NONE
	  
	    REAL*8,INTENT(IN):: K,Z,H
 
        IF (K*(Z+H).GT.89.4d0) THEN
          C2SH=DEXP(K*Z)
        ELSE
          C2SH=DCOSH(K*(Z+H))/DSINH(K*H)  
        ENDIF
        
        
        RETURN
        END        
        
!    ----------------------------------------------------------------------------------------
!       For caculation of sinh(k(z+h))/sinh(kh)
!    ----------------------------------------------------------------------------------------
!
        REAL*8 FUNCTION  S2SH(K,Z,H)
        
	    IMPLICIT    NONE
	  
	    REAL*8,INTENT(IN):: K,Z,H

        IF (K*(Z+H).GT.89.4d0) THEN
          S2SH=DEXP(K*Z)
        ELSE
          S2SH=DSINH(K*(Z+H))/DSINH(K*H)
        ENDIF
        
        
        RETURN
        END
        
 !=======================================================================
      FUNCTION WaveNumber ( Omega, g, depth )


         ! This FUNCTION solves the finite depth dispersion relationship:
         !
         !                   k*tanh(k*h)=(Omega^2)/g
         !
         ! for k, the wavenumber (WaveNumber) given the frequency, Omega,
         ! gravitational constant, g, and water depth, h, as inputs.  A
         ! high order initial guess is used in conjunction with a quadratic
         ! Newton's method for the solution with seven significant digits
         ! accuracy using only one iteration pass.  The method is due to
         ! Professor J.N. Newman of M.I.T. as found in routine EIGVAL of
         ! the SWIM-MOTION-LINES (SML) software package in source file
         ! Solve.f of the SWIM module.


      IMPLICIT                        NONE


         ! Passed Variables:

      REAL*8, INTENT(IN )      :: g                                               ! Gravitational acceleration (m/s^2)
      REAL*8, INTENT(IN )      :: depth                                           ! Water depth (meters)
      REAL*8, INTENT(IN )      :: Omega                                           ! Wave frequency (rad/s)
      REAL*8                   :: WaveNumber                                      ! This function = wavenumber, k (1/m)


         ! Local Variables:
         
      REAL*8                   :: h                                               ! Positive water depth (meters)
      REAL*8                   :: A                                               ! A temporary variable used in the solution.
      REAL*8                   :: B                                               ! A temporary variable used in the solution.
      REAL*8                   :: C                                               ! A temporary variable used in the solution.
      REAL*8                   :: C2                                              ! A temporary variable used in the solution.
      REAL*8                   :: CC                                              ! A temporary variable used in the solution.
      REAL*8                   :: E2                                              ! A temporary variable used in the solution.
      REAL*8                   :: X0                                              ! A temporary variable used in the solution.
      
      
      IF ( depth < 0.0 )  THEN
          
         h = DABS(depth) 
         
      ENDIF   

         ! Compute the wavenumber, unless Omega is zero, in which case, return
         !   zero:

      IF ( Omega == 0.0 )  THEN  ! When .TRUE., the formulation below is ill-conditioned; thus, the known value of zero is returned.


         WaveNumber = 0.0


      ELSE                       ! Omega > 0.0; solve for the wavenumber as usual.


         C  = Omega*Omega*h/g
         CC = C*C


         ! Find X0:

         IF ( C <= 2.0 )  THEN

            X0 = SQRT(C)*( 1.0 + C*( 0.169 + (0.031*C) ) )

         ELSE

            E2 = EXP(-2.0*C)

            X0 = C*( 1.0 + ( E2*( 2.0 - (12.0*E2) ) ) )

         ENDIF


         ! Find the WaveNumber:

         IF ( C <= 4.8 )  THEN

            C2 = CC - X0*X0
            A  = 1.0/( C - C2 )
            B  = A*( ( 0.5*LOG( ( X0 + C )/( X0 - C ) ) ) - X0 )

            WaveNumber = ( X0 - ( B*C2*( 1.0 + (A*B*C*X0) ) ) )/h

         ELSE

            WaveNumber = X0/h

         ENDIF


      ENDIF



      RETURN
      END FUNCTION WaveNumber       
        