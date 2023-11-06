!/************************************************************************/
!/*                                                                      */
!/* ALGORITHM 412 - Calculation of irrgular wave height                  */
!/*                                                                      */
!/* Yingyi LIU 16 Dec 2013.                                              */
!/* With modifications suggested by Yingyi LIU,                          */
!/* Research Institute of Applied Mechanics, Kyshu-U, Japan,             */
!/* 8 December 2013.                                                     */
!/*                                                                      */
!/************************************************************************/
  
      SUBROUTINE IWHC(XLC,TSP,NSTP,H,HS,TP,BETA,ETA)  
      
	  USE Const_mod    
      USE mt95_mod
      
      IMPLICIT NONE
      
      INTEGER M,I,J,L,NSTP
      REAL*8 WMIN,WMAX,DW,SW,EPSL,KX,X,Y
      REAL*8 H,HS,TP,BETA,TSP,T,XLC(3),ETA(0:NSTP)
      REAL*8,ALLOCATABLE::W(:),WW(:),WI(:),KI(:),AI(:),RA(:)
      REAL*8,EXTERNAL:: WaveNumber
      
      X=XLC(1)
      Y=XLC(2)
      
      M=96
      ALLOCATE(W(1:M+1),WW(1:M),WI(1:M),KI(1:M),AI(1:M),RA(M))      
      
      WMIN=0.1D0
      WMAX=2.98D0
      DW=(WMAX-WMIN)/M
      
      DO I=0,M
       W(I+1)=WMIN+I*DW 
      ENDDO 
!      
      call genrand_init()
      call genrand_real1(RA)
!     call random_seed()
      DO I=1,M
!       CALL RANDOM_NUMBER(RA(I))
      ENDDO
      
      DO I=1,M     
        WW(I)=(W(I)+W(I+1))/2.D0
        WI(I)=W(I)+DW*RA(I)
        KI(I)= WaveNumber ( WI(I), G, H )
!       CALL PM(HS,TP,WW(I),PI,SW)        
        CALL JONSWAP_Goda( WW(I), Hs, Tp , Sw )
        AI(I)=DSQRT(2.D0*SW*DW)        
      ENDDO
      
!     call sleepqq(10000)

      call genrand_real1(RA)
!     call random_seed()
      DO I=1,M
!       CALL RANDOM_NUMBER(RA(I))
      ENDDO 
       
      DO 100 L=0,NSTP
        T=TSP*L
        ETA(L)=0.D0
       
      DO 100 I=1,M
          
        KX=KI(I)*(X*DCOS(BETA)+Y*DSIN(BETA))
        EPSL=2.D0*PI*RA(I)        
        ETA(L)=ETA(L)+AI(I)*DCOS(KX-WI(I)*T+EPSL)

100   CONTINUE
      
      DO L=0,NSTP
        T=TSP*L             
        WRITE(62,*)  T, ETA(L)
      ENDDO
      
      DO I=1,M
        CALL JONSWAP_Goda( WW(I), Hs, Tp , Sw )
        WRITE(75,*)  WW(I), Sw
      ENDDO     
      
      DEALLOCATE(W,WW,WI,KI,AI,RA)
      
      RETURN 
      END SUBROUTINE IWHC
!============================================================      
subroutine init_random_seed()
   integer :: i, n, clock
   integer, dimension(:), allocatable :: seed
   
   call random_seed(size = n)
   allocate(seed(n))
   
   call system_clock(count=clock)
   
   seed = clock + 37 * (/ (i - 1, i = 1, n) /)
   call random_seed(put = seed)
   
   deallocate(seed)
end subroutine init_random_seed       
       
!===============================================================
       SUBROUTINE RGWAVEXFC(W,DW,I,T,TFOR,APH)

 	   USE WaveDyn_mod
       IMPLICIT   NONE  

	   INTEGER I,J,NUM
       REAL*8 W,DW,APH,AMP,T,TFOR,EXF(2)
          
       
       NUM = INT((W-AW(0))/DW)
        
!        WRITE(*,*) 'NUM', NUM       
!        WRITE(*,*) 'AW',AW(NUM),AW(NUM+1)
        
       DO J=1,2
        EXF(J) =((EFOR(I,J,NUM+1)-EFOR(I,J,NUM))/(AW(NUM+1)-AW(NUM)))*(W-AW(NUM))+EFOR(I,J,NUM)
       ENDDO
        
       AMP=DSQRT(EXF(1)**2+EXF(2)**2)
       APH=DATAN2(EXF(2),EXF(1))
        
       TFOR=AMP*DCOS(W*T-APH)
        
        
       RETURN        
       END               
!/************************************************************************/
!/*                                                                      */
!/* ALGORITHM 412 - Calculation of Pierson-Moskowitz spectrum            */
!/*                                                                      */
!/* Yingyi LIU 26 Nov 2013.                                              */
!/* With modifications suggested by Yingyi LIU,                          */
!/* Research Institute of Applied Mechanics, Kyshu-U, Japan,             */
!/* 8 December 2013.                                                     */
!/*                                                                      */
!/************************************************************************/
  
      SUBROUTINE PM(HS,TP,W,PI,SW)  

      IMPLICIT NONE
      REAL*8 HS,TP,W,PI,SW,PI2,CON
      
      PI2=2.D0*PI
      CON=W*TP/PI2
      
      SW=1.D0/PI2*5.D0/16.D0*HS**2*TP*CON**(-5.D0)*DEXP(-1.25D0*CON**(-4.D0))      
      
      RETURN 
      END
      
!=======================================================================      
        SUBROUTINE JONSWAP_Goda ( Omega, Hs, Ts , Sptr )


         ! This FUNCTION computes the JOint North Sea WAve Project
         ! (JONSWAP) representation of the one-sided power spectral density
         ! or wave spectrum given the frequency, Omega, peak shape
         ! parameter, Gamma, significant wave height, Hs, and significant spectral
         ! period, Ts, as inputs.  The value of Gamma is set to 3.30.
         !
         ! There are several different versions of the JONSWAP spectrum
         ! formula.  This version is based on the one documented in the
         ! Goda(1999).



      IMPLICIT NONE


         ! Passed Variables:

      REAL(8)              :: Gamma     ! Peak shape parameter (-)
      REAL(8), INTENT(IN ) :: Hs     ! Significant wave height (meters)
      REAL(8)              :: Sptr   ! This function = JONSWAP wave 
                                        !  spectrum,S (m^2/(rad/s))
      REAL(8), INTENT(IN ) :: Omega     ! Wave frequency (rad/s)
      REAL(8), INTENT(IN ) :: Ts        ! Peak spectral period (sec)


         ! Local Variables:

      REAL(8) :: Alpha   ! Exponent on Gamma used in the spectral
                         ! formulation (-)
      REAL(8) :: C, BJ       ! Normalising factor used in the spectral
                         ! formulation (-)
      REAL(8) :: f       ! Wave frequency (Hz)
      REAL(8) :: Tp,fp      ! Peak spectral frequency (Hz)
      REAL(8) :: Tpf_4   ! (Tp*f)^(-4)
      REAL(8) :: Sigma     ! Scaling factor used in the spectral 
                         ! formulation (-)
      REAL(8) :: Inv2Pi  ! 0.5D0 /PI
      REAL(8) :: Pi  ! 0.5D0 /PI     
      
      ! the peak shape parameter of the incident wave spectrum conditioned 
      ! on significant wave height and peak spectral period (-)      
            
      
      ! Compute the default peak shape parameter of the incident wave spectrum,
      !   conditioned on significant wave height and peak spectral period:
      
      PI = 4.D0 * DATAN( 1.D0 )
      Inv2Pi = 0.5D0 /PI

      Tp=Ts/(1.D0-0.132D0*(Gamma+0.2D0)**(-0.559D0))
      
!      WRITE(*,'(3(A8,F12.4))') 'Ts',Ts,'Tp',Tp,'W1',2.D0*PI/TP
!      PAUSE  
      
      Gamma = 3.30D0

         ! Compute the JONSWAP wave spectrum, unless Omega is zero, in which case,
         !   return zero:

      IF ( Omega == 0.0 )  THEN  ! When .TRUE., the formulation below is
          !  ill-conditioned; thus, the known value of zero is returned.


         Sptr  = 0.0


      ELSE     ! Omega > 0.0; forumulate the JONSWAP spectrum.


         ! Compute the wave frequency and peak spectral frequency in Hz:

         f        = Inv2Pi*Omega
         fp       = 1/Tp
         Tpf_4  = (Tp*f)**(-4.0)


         ! Compute the normalising factor:

         C   = 1.094D0 - 0.01915D0*DLOG(Gamma)
         BJ  = 0.06238D0/(0.230D0+0.0336D0*Gamma-0.185D0/(1.9D0+Gamma))*C

         ! Compute Alpha:

         IF ( f <= fp )  THEN
            Sigma = 0.07D0
         ELSE
            Sigma = 0.09D0
         ENDIF

         Alpha    = DEXP(  -0.5D0*( ( (f/fp) - 1.0 )/Sigma )**2.0 )


         ! Compute the wave spectrum:

         Sptr  = ( Inv2Pi*BJ*Hs*Hs*Tpf_4/f )*DEXP( ( -1.25*Tpf_4 ) )*( Gamma**Alpha )


      ENDIF


      RETURN
        END SUBROUTINE JONSWAP_Goda
        
        
        
        
!=======================================================================      
        SUBROUTINE JONSWAP_IEC ( Omega, Hs, Tp , Sptr )


         ! This FUNCTION computes the JOint North Sea WAve Project
         ! (JONSWAP) representation of the one-sided power spectral density
         ! or wave spectrum given the frequency, Omega, peak shape
         ! parameter, Gamma, significant wave height, Hs, and peak spectral
         ! period, Tp, as inputs.  If the value of Gamma is 1.0, the
         ! Pierson-Moskowitz wave spectrum is returned.
         !
         ! There are several different versions of the JONSWAP spectrum
         ! formula.  This version is based on the one documented in the
         ! IEC61400-3 wind turbine design standard for offshore wind
         ! turbines.



      IMPLICIT NONE


         ! Passed Variables:

      REAL(8)              :: Gamma     ! Peak shape parameter (-)
      REAL(8), INTENT(IN ) :: Hs     ! Significant wave height (meters)
      REAL(8)              :: Sptr   ! This function = JONSWAP wave 
                                        !  spectrum,S (m^2/(rad/s))
      REAL(8), INTENT(IN ) :: Omega     ! Wave frequency (rad/s)
      REAL(8), INTENT(IN ) :: Tp        ! Peak spectral period (sec)


         ! Local Variables:

      REAL(8) :: Alpha   ! Exponent on Gamma used in the spectral
                         ! formulation (-)
      REAL(8) :: C, BJ       ! Normalising factor used in the spectral
                         ! formulation (-)
      REAL(8) :: f       ! Wave frequency (Hz)
      REAL(8) :: fp      ! Peak spectral frequency (Hz)
      REAL(8) :: Tpf_4   ! (Tp*f)^(-4)
      REAL(8) :: Sigma     ! Scaling factor used in the spectral formulation (-)
      REAL(8) :: TpOvrSqrtHs     ! Tp over square root of Hs
      REAL(8) :: Inv2Pi  ! 0.5D0 /PI
      REAL(8) :: Pi  ! 0.5D0 /PI
      
      
      ! the peak shape parameter of the incident wave spectrum conditioned 
      ! on significant wave height and peak spectral period (-)      
            
      
      ! Compute the default peak shape parameter of the incident wave spectrum,
      !   conditioned on significant wave height and peak spectral period:
      
      PI = 4.D0 * DATAN( 1.D0 )
      Inv2Pi = 0.5D0 /PI
      TpOvrSqrtHs = Tp/DSQRT(Hs)
      

         ! Compute the JONSWAP wave spectrum, unless Omega is zero, in which case,
         !   return zero:

      IF ( Omega == 0.0 )  THEN  ! When .TRUE., the formulation below is
          !  ill-conditioned; thus, the known value of zero is returned.


         Sptr  = 0.0


      ELSE     ! Omega > 0.0; forumulate the JONSWAP spectrum.


         ! Compute the wave frequency and peak spectral frequency in Hz:

         f        = Inv2Pi*Omega
         fp       = 1/Tp
         Tpf_4  = (Tp*f)**(-4.0)

         ! Compute wave peak shape default value:

         IF ( TpOvrSqrtHs <= 3.6D0 )  THEN          
           Gamma = 5.0D0          
         ELSEIF ( TpOvrSqrtHs >= 5.0D0 )  THEN          
           Gamma = 1.0D0          
         ELSE          
           Gamma = DEXP( 5.75D0 - 1.15D0*TpOvrSqrtHs )
         ENDIF
         
         ! Compute the normalising factor:

         C   = 1.D0 - 0.287D0*DLOG(Gamma)
         BJ  = 5.D0/16.D0*C

         ! Compute Alpha:

         IF ( f <= fp )  THEN
            Sigma = 0.07D0
         ELSE
            Sigma = 0.09D0
         ENDIF

         Alpha    = DEXP(  -0.5D0*( ( (f/fp) - 1.0 )/Sigma )**2.0 )


         ! Compute the wave spectrum:

         Sptr  = ( Inv2Pi*BJ*Hs*Hs*Tpf_4/f )*DEXP( ( -1.25*Tpf_4 ) )*( Gamma**Alpha )


      ENDIF


      RETURN
      END SUBROUTINE JONSWAP_IEC         
        
        
        