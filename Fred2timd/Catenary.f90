      
      SUBROUTINE SmllRotTrans( Theta1, Theta2, Theta3, TransMat )
      
      


      ! This routine computes the 3x3 transformation matrix, TransMat,
      !   to a coordinate system x (with orthogonal axes x1, x2, x3)
      !   resulting from three rotations (Theta1, Theta2, Theta3) about the
      !   orthogonal axes (X1, X2, X3) of coordinate system X.  All angles
      !   are assummed to be small, as such, the order of rotations does
      !   not matter and Euler angles do not need to be used.  This routine
      !   is used to compute the transformation matrix (TransMat) between
      !   undeflected (X) and deflected (x) coordinate systems.  In matrix
      !   form:
      !      {x1}   [TransMat(Theta1, ] {X1}
      !      {x2} = [         Theta2, ]*{X2}
      !      {x3}   [         Theta3 )] {X3}

      ! The transformation matrix, TransMat, is the closest orthonormal
      !   matrix to the nonorthonormal, but skew-symmetric, Bernoulli-Euler
      !   matrix:
      !          [   1.0    Theta3 -Theta2 ]
      !      A = [ -Theta3   1.0    Theta1 ]
      !          [  Theta2 -Theta1   1.0   ]
      !
      !   In the Frobenius Norm sense, the closest orthornormal matrix is:
      !      TransMat = U*V^T,
      !
      !   where the columns of U contain the eigenvectors of A*A^T and the
      !   columns of V contain the eigenvectors of A^T*A (^T = transpose).
      !   This result comes directly from the Singular Value Decomposition
      !   (SVD) of A = U*S*V^T where S is a diagonal matrix containing the
      !   singular values of A, which are SQRT( eigenvalues of A*A^T ) =
      !   SQRT( eigenvalues of A^T*A ).

      ! The algebraic form of the transformation matrix, as implemented
      !   below, was derived symbolically by J. Jonkman by computing U*V^T
      !   by hand with verification in Mathematica.

      ! NOTE: This routine is an exact copy of SUBROUTINE SmllRotTrans() given in
      !       FAST source file FAST.f90.


       IMPLICIT                        NONE


      ! Passed Variables:

       REAL(8), INTENT(IN )      :: Theta1    ! The small rotation about X1, (rad).
       REAL(8), INTENT(IN )      :: Theta2                                          ! The small rotation about X2, (rad).
       REAL(8), INTENT(IN )      :: Theta3                                          ! The small rotation about X3, (rad).
       REAL(8), INTENT(OUT)      :: TransMat (3,3)                                  ! The resulting transformation matrix from X to x, (-).


       ! Local Variables:

       REAL(8)                   :: ComDenom                                        ! = ( Theta1^2 + Theta2^2 + Theta3^2 )*SQRT( 1.0 + Theta1^2 + Theta2^2 + Theta3^2 )
       REAL(8), PARAMETER        :: LrgAngle  = 0.4                                 ! Threshold for when a small angle becomes large (about 23deg).  This comes from: COS(SmllAngle) ~ 1/SQRT( 1 + SmllAngle^2 ) and SIN(SmllAngle) ~ SmllAngle/SQRT( 1 + SmllAngle^2 ) results in ~5% error when SmllAngle = 0.4rad.
       REAL(8)                   :: Theta11                                         ! = Theta1^2
       REAL(8)                   :: Theta12S                                        ! = Theta1*Theta2*[ SQRT( 1.0 + Theta1^2 + Theta2^2 + Theta3^2 ) - 1.0 ]
       REAL(8)                   :: Theta13S                                        ! = Theta1*Theta3*[ SQRT( 1.0 + Theta1^2 + Theta2^2 + Theta3^2 ) - 1.0 ]
       REAL(8)                   :: Theta22                                         ! = Theta2^2
       REAL(8)                   :: Theta23S                                        ! = Theta2*Theta3*[ SQRT( 1.0 + Theta1^2 + Theta2^2 + Theta3^2 ) - 1.0 ]
       REAL(8)                   :: Theta33                                         ! = Theta3^2
       REAL(8)                   :: SqrdSum                                         ! = Theta1^2 + Theta2^2 + Theta3^2
       REAL(8)                   :: SQRT1SqrdSum                                    ! = SQRT( 1.0 + Theta1^2 + Theta2^2 + Theta3^2 )

       LOGICAL,    SAVE             :: FrstWarn  = .TRUE.                              ! When .TRUE., indicates that we're on the first warning.

      ! Display a warning message if at least one angle gets too large in
      !   magnitude:

       IF ( ( ( DABS(Theta1) > LrgAngle ) .OR. ( DABS(Theta2) > LrgAngle ) .OR. ( DABS(Theta3) > LrgAngle ) ) .AND. FrstWarn )  THEN
           
!          WRITE(*,*)Theta1,'  ',Theta2,'    ',Theta3
!
       Print*, ' WARNING: '
       Print*, ' Small angle assumption violated in SUBROUTINE SmllRotTrans()'
!
       FrstWarn = .FALSE.   ! Don't enter here again!
!
       ENDIF

!      ! Compute some intermediate results:
!
       Theta11      = Theta1*Theta1
       Theta22      = Theta2*Theta2
       Theta33      = Theta3*Theta3
!
       SqrdSum      = Theta11 + Theta22 + Theta33
       SQRT1SqrdSum = DSQRT( 1.0 + SqrdSum )
       ComDenom     = SqrdSum*SQRT1SqrdSum
!
       Theta12S     = Theta1*Theta2*( SQRT1SqrdSum - 1.0 )
       Theta13S     = Theta1*Theta3*( SQRT1SqrdSum - 1.0 )
       Theta23S     = Theta2*Theta3*( SQRT1SqrdSum - 1.0 )
!
!
       ! Define the transformation matrix:
!
       IF ( ComDenom == 0.0 )  THEN  ! All angles are zero and matrix is ill-conditioned (the matrix is derived assuming that the angles are not zero); return identity
!
       TransMat(1,:) = (/ 1.0, 0.0, 0.0 /)
       TransMat(2,:) = (/ 0.0, 1.0, 0.0 /)
       TransMat(3,:) = (/ 0.0, 0.0, 1.0 /)
!
       ELSE                          ! At least one angle is nonzero
!
       TransMat(1,1) = ( Theta11*SQRT1SqrdSum + Theta22              + Theta33              )/ComDenom
       TransMat(2,2) = ( Theta11              + Theta22*SQRT1SqrdSum + Theta33              )/ComDenom
       TransMat(3,3) = ( Theta11              + Theta22              + Theta33*SQRT1SqrdSum )/ComDenom
       TransMat(1,2) = (  Theta3*SqrdSum + Theta12S )/ComDenom
       TransMat(2,1) = ( -Theta3*SqrdSum + Theta12S )/ComDenom
       TransMat(1,3) = ( -Theta2*SqrdSum + Theta13S )/ComDenom
       TransMat(3,1) = (  Theta2*SqrdSum + Theta13S )/ComDenom
       TransMat(2,3) = (  Theta1*SqrdSum + Theta23S )/ComDenom
       TransMat(3,2) = ( -Theta1*SqrdSum + Theta23S )/ComDenom
!
       ENDIF
!
!
       RETURN
      END SUBROUTINE SmllRotTrans
      
      
! ====================================================================   
      SUBROUTINE LargRotTrans( Theta1, Theta2, Theta3, TransMat )


      ! This routine computes the 3x3 transformation matrix, TransMat,
      !   to a coordinate system x (with orthogonal axes x1, x2, x3)
      !   resulting from three rotations (Theta1, Theta2, Theta3) about the
      !   orthogonal axes (X1, X2, X3) of coordinate system X.  All angles
      !   are assummed to be small, as such, the order of rotations does
      !   not matter and Euler angles do not need to be used.  This routine
      !   is used to compute the transformation matrix (TransMat) between
      !   undeflected (X) and deflected (x) coordinate systems.  In matrix
      !   form:
      !      {x1}   [TransMat(Theta1, ] {X1-XG1}
      !      {x2} = [         Theta2, ]*{X2-XG2}
      !      {x3}   [         Theta3 )] {X3-XG3}
      
      !   and the inverse transform:
      !      {X1}   [XG1]+[TransMat(Theta1, ]^T {x1}
      !      {X2} = [XG2]+[         Theta2, ]^T*{x2}
      !      {X3}   [XG3]+[         Theta3 )]^T {x3}      

      ! The transformation matrix, TransMat, is the closest orthonormal
      !   matrix to the nonorthonormal, but skew-symmetric, Bernoulli-Euler
      !   matrix:
      !          [   c5c6   c4s6+s4s5c6   s4s6-c4s5c6 ]
      !      A = [  -c5s6   c4c6-s4s5s6   s4c6+c4s5s6 ]
      !          [    s5       -s4c5          c4c5    ]
      ! where ci=cos(Thetai), si=sin(Thetai),ti=tan(Thetai)
       
      ! NOTE: This routine is a MODIFICATION of SUBROUTINE SmllRotTrans() given in
      !       FAST source file FAST.f90.


       IMPLICIT                        NONE

      ! Passed Variables:

       REAL(8), INTENT(IN )      :: Theta1    ! The large rotation about X1, (rad).
       REAL(8), INTENT(IN )      :: Theta2                                          ! The large rotation about X2, (rad).
       REAL(8), INTENT(IN )      :: Theta3                                          ! The large rotation about X3, (rad).
       REAL(8), INTENT(OUT)      :: TransMat (3,3)                                  ! The resulting transformation matrix from X to x, (-).


       ! Local Variables:

       REAL(8)                   :: C4,C5,C6,S4,S5,S6                               

!      ! Compute some intermediate results:
!
       C4 = DCOS(Theta1)
       C5 = DCOS(Theta2)
       C6 = DCOS(Theta3)       
       
       S4 = DSIN(Theta1)
       S5 = DSIN(Theta2)
       S6 = DSIN(Theta3)   
!
!
       ! Define the transformation matrix:
!
       TransMat(1,1) = C5*C6
       TransMat(2,2) = C4*C6-S4*S5*S6
       TransMat(3,3) = C4*C5
       TransMat(1,2) = C4*S6+S4*S5*C6
       TransMat(2,1) = -C5*S6
       TransMat(1,3) = S4*S6-C4*S5*C6
       TransMat(3,1) = S5
       TransMat(2,3) = S4*C6+C4*S5*S6
       TransMat(3,2) = -S4*C5

!
       RETURN
      END SUBROUTINE LargRotTrans
      
! ====================================================================   
      SUBROUTINE CatenarySystem( StageType, NumLines, X0, X, F_Lines ,LFairTe )
      USE omp_lib
      USE CABLE
      USE Const_mod
       IMPLICIT                        NONE

  ! Passed Variables:
      
   INTEGER(4), INTENT(IN )   :: NumLines                                        ! Number of mooring lines (-)    
   LOGICAL, INTENT(IN )      :: StageType                                       ! When .FALSE., indicates that we're on the first step.  
   REAL(8), INTENT(IN )      :: X0       (6)                                    ! The initial 3 components of the translational displacement of the mass center of platform in inertial frame (in m  )
   REAL(8), INTENT(IN )      :: X        (6)                                    !  The 3 components of the translational displacement    (in m  )        of the platform reference and the 3 components of the rotational displacement        (in rad  )        of the platform relative to the inertial frame
   REAL(8), INTENT(OUT)      :: F_Lines  (6)                                    ! Total load contribution from all mooring lines (N, N-m)
   REAL(8)                   :: F_ILine  (6)                                    ! Load contribution from EVERY mooring lines (N, N-m)
      
   REAL(8)                   :: COSPhi                                          ! Cosine of the angle between the xi-axis of the inertia frame and the X-axis of the local coordinate system of the current mooring line (-)
   REAL(8)                   :: SINPhi                                          ! Sine   of the angle between the xi-axis of the inertia frame and the X-axis of the local coordinate system of the current mooring line (-)
   REAL(8)                   :: TransMat (3,3)                                  ! Transformation matrix from the inertial frame to the tower base / platform coordinate system (-)
!jmj End of proposed change.  v6.02b-jmj  15-Nov-2006.
     ! Local Variables:
   REAL(8)                   :: XF                                              ! Horizontal distance between anchor and fairlead of the current mooring line (meters)
   REAL(8)                   :: ZF                                              ! Vertical   distance between anchor and fairlead of the current mooring line (meters)
   REAL(8)                   :: LFairxi                                         ! xi-coordinate of the current fairlead in the inertial frame coordinate system (meters)
   REAL(8)                   :: LFairxiRef                                      ! xi-coordinate of the current fairlead relative to the platform reference point in the inertial frame coordinate system (meters)
   REAL(8)                   :: LFairxiTe                                       ! xi-component of the effective tension at the fairlead of the current mooring line (N)
   REAL(8)                   :: LFairyi                                         ! yi-coordinate of the current fairlead in the inertial frame coordinate system (meters)
   REAL(8)                   :: LFairyiRef                                      ! yi-coordinate of the current fairlead relative to the platform reference point in the inertial frame coordinate system (meters)
   REAL(8)                   :: LFairyiTe                                       ! yi-component of the effective tension at the fairlead of the current mooring line (N)
   REAL(8)                   :: LFairzi                                         ! zi-coordinate of the current fairlead in the inertial frame coordinate system (meters)
   REAL(8)                   :: LFairziRef                                      ! zi-coordinate of the current fairlead relative to the platform reference point in the inertial frame coordinate system (meters)
   REAL(8)                   :: LFairziTe                                       ! zi-component of the effective tension at the fairlead of the current mooring line (N)

   REAL(8)                   :: LAnchxi   (NumLines)                         ! xi-coordinate of each anchor   in the inertial frame        coordinate system (meters)
   REAL(8)                   :: LAnchyi   (NumLines)                         ! yi-coordinate of each anchor   in the inertial frame        coordinate system (meters)
   REAL(8)                   :: LAnchzi   (NumLines)                         ! zi-coordinate of each anchor   in the inertial frame        coordinate system (meters)
   REAL(8)                   :: LFairxt   (NumLines)                         ! xt-coordinate of each fairlead in the tower base / platform coordinate system (meters)
   REAL(8)                   :: LFairyt   (NumLines)                         ! yt-coordinate of each fairlead in the tower base / platform coordinate system (meters)
   REAL(8)                   :: LFairzt   (NumLines)                         ! zt-coordinate of each fairlead in the tower base / platform coordinate system (meters)
   REAL(8)                   :: LFairHTe  (NumLines)
   REAL(8)                   :: LFairVTe  (NumLines) 
   REAL(8)                   :: LFairTe  (NumLines)    
   REAL(8)                   :: LAnchHTe  (NumLines)
   REAL(8)                   :: LAnchVTe  (NumLines)   
   REAL(8)                   :: LAnchTe  (NumLines) 
   
   REAL(8)                   :: LMassDen  (NumLines)                         ! Mass density of each mooring line (kg/m)
   REAL(8)                   :: LSeabedCD (NumLines)                         ! Coefficient of seabed static friction drag of each mooring line (a negative value indicates no seabed) (-)
   REAL(8)                   :: LUnstrLen (NumLines)                         ! Unstretched length of each mooring line (meters)
   REAL(8)                   :: LEAStff   (NumLines)                         ! Extensional stiffness of line (N)   
   REAL(8)                   :: LFldWght  (NumLines)                         ! Weight of each mooring line in fluid per unit length (N/m)  
   REAL(8)                   :: LDiam     (NumLines)                         ! Effective diameter of each mooring line for calculation of the line buoyancy (meters)
   

   REAL(8)                   :: XF2                                          ! = XF*XF
   REAL(8)                   :: ZF2                                          ! = ZF*ZF 

   REAL(8)                   :: WtrDens      
   REAL(8)                   :: Gravity  
   REAL(8)                   :: Lamda0                                       ! Catenary parameter used to generate the initial guesses of the horizontal and vertical tensions at the fairlead for the Newton-Raphson iteration (-)
   !REAL(8)                   :: PI 
   !REAL(8)                   :: PiOvr4    
  
   REAL(8)                   :: AngLines(NumLines)  
   
   INTEGER(4)                :: I,NumThread,J                                  ! Index for counting iterations or looping through line nodes (-)
   REAL(8)                   :: F_LINE_THREAD(6,Numlines)                    ! Cont the F_ILINE in multi-Thread, zhu.
  
   ! Initilization of the parameters:
     
     WtrDens = RHOW
     Gravity =  9.807d0
     !PiOvr4 =  DATAN(1.d0)
     !PI = 4.D0 * PiOvr4
     
      ! First initialize the total load contribution from all mooring lines to
      ! zero:

! the following is for 3 mooring lines
     DO I = 1,NumLines
         AngLines(I) = PI/180.0D0*LAngle(I)
     END DO 
     

! the following is for 6 mooring lines:
      !
      !AngLines(1) = PI/3.D0 
      !AngLines(2) = PI/3.D0*2.D0   
      !AngLines(3) = PI/6.D0*5.D0 
      !AngLines(4) = PI/6.D0*7.D0   
      !AngLines(5) = PI/6.D0*9.D0   
      !AngLines(6) = PI/6.D0*11.D0 
      
! the following is for 9 mooring lines:

!      AngLines(1) = 0.D0      
!      AngLines(2) = PI/6.D0 
!      AngLines(3) = PI/6.D0*3.D0 
!      AngLines(4) = PI/6.D0*4.D0 
!      AngLines(5) = PI/6.D0*5.D0      
!      AngLines(6) = PI/6.D0*7.D0    
!      AngLines(7) = PI/6.D0*8.D0   
!      AngLines(8) = PI/6.D0*9.D0       
!      AngLines(9) = PI/6.D0*11.D0 
      
      F_Lines(:) = 0.0
      
     DO I = 1,NumLines ! Loop through all mooring lines     
       LDiam(I) =      LineDiam                             !放到输入文件中，从外部输入
       LMassDen(I) =   LineMassDen                             !放到输入文件中，从外部输入
       LUnstrLen(I) =  LineUnstrLen                           !放到输入文件中，从外部输入
       LSeabedCD(I) =  LineSeabedCD                            !放到输入文件中，从外部输入
       LEAStff(I) = LineEAStff                           !放到输入文件中，从外部输入
       
      ! Compute the weight of each mooring line in fluid per unit length based on
      !   their mass density and effective diameter, water density, and gravity:
      ! NOTE: The buoyancy is calculated assuming that the entire length of the
      !       mooring line is submerged in the water.
      ! 计算湿重，以及锚链两端的初始坐标
       LFldWght(I) = ( LMassDen(I) - WtrDens*Pi*LDiam(I)*LDiam(I)/4 )*Gravity 
     
       LFairxt(I) =   RFair *DCOS( AngLines(I) ) - X0(1)
       LFairyt(I) =   RFair *DSIN( AngLines(I) ) - X0(2)
       LFairzt(I) = - DepthFair - X0(3)
       
       LAnchxi(I) =   RAnch *DCOS( AngLines(I) ) 
       LAnchyi(I) =   RAnch *DSIN( AngLines(I) ) 
       LAnchzi(I) = - DepthAnch 
       
     ENDDO   ! I - All mooring lines
     
     
     IF ( StageType == .FALSE.) THEN
         
         !WRITE(*,*)'X =  ', X(4),'    ',X(5),'    ',X(6)
     
      CALL SmllRotTrans (  X(4), X(5), X(6), TransMat ) 
      
      DO I = 1,NumLines ! Loop through all mooring lines     
          
      ! First initialize the load contribution from EVERY mooring lines to
      !   zero:

      F_ILine(:) = 0.0
      
      ! Transform the fairlead location from the initial platform to the inertial
      !    frame coordinate system:
      ! NOTE: TransMat0^T = TransMat0^-1 where ^T = matrix transpose and ^-1 =
      !       matrix inverse. 
      
         LFairxiRef = TransMat(1,1)*LFairxt(I) + TransMat(2,1)*LFairyt(I) + TransMat(3,1)*LFairzt(I)
         LFairyiRef = TransMat(1,2)*LFairxt(I) + TransMat(2,2)*LFairyt(I) + TransMat(3,2)*LFairzt(I)
         LFairziRef = TransMat(1,3)*LFairxt(I) + TransMat(2,3)*LFairyt(I) + TransMat(3,3)*LFairzt(I)

         LFairxi    = X(1) + LFairxiRef
         LFairyi    = X(2) + LFairyiRef
         LFairzi    = X(3) + LFairziRef
         WRITE(96,*) 'line number: ',I,LFairxi,LFairyi,LFairzi

      ! Transform the fairlead location from the inertial frame coordinate system
      !   to the local coordinate system of the current line (this coordinate
      !   system lies at the current anchor, Z being vertical, and X directed from
      !   the current anchor to the current fairlead):

         XF      = DSQRT( ( LFairxi - LAnchxi(I) )**2.0D0 + ( LFairyi - LAnchyi(I) )**2.0D0 )
         ZF      =          LFairzi - LAnchzi(I)

         XF2     = XF*XF
         ZF2     = ZF*ZF
         
         IF ( XF == 0.0D0 )  THEN  ! .TRUE. if the current mooring line is exactly vertical; thus, the solution below is ill-conditioned because the orientation is undefined; so set it such that the tensions and nodal positions are only vertical
            COSPhi  = 0.0D0
            SINPhi  = 0.0D0
            Print *, 'WARNING: Current mooring line is exact vertical.'
         ELSE                    ! The current mooring line must not be vertical; use simple trigonometry
            COSPhi  =       ( LFairxi - LAnchxi(I) )/XF
            SINPhi  =       ( LFairyi - LAnchyi(I) )/XF
         ENDIF

      ! Generate the initial guess values for the horizontal and vertical tensions
      !   at the fairlead in the Newton-Raphson iteration for the catenary mooring
      !   line solution.  Use starting values documented in: Peyrot, Alain H. and
      !   Goulois, A. M., "Analysis Of Cable Structures," Computers & Structures,
      !   Vol. 10, 1979, pp. 805-813:

         IF     ( XF           == 0.0D0             )  THEN ! .TRUE. if the current mooring line is exactly vertical
            Lamda0 = 1.0E+06
         ELSEIF ( LUnstrLen(I) <= DSQRT( XF2 + ZF2 ) )  THEN ! .TRUE. if the current mooring line is taut
            Lamda0 = 0.2D0
         ELSE                                               ! The current mooring line must be slack and not vertical
            Lamda0 = DSQRT( 3.0D0*( ( LUnstrLen(I)*LUnstrLen(I) - ZF2 )/XF2 - 1.0D0 ) )
         ENDIF

         LFairHTe(I) = DABS( 0.5D0*LFldWght(I)*  XF/     Lamda0                   )
         LFairVTe(I) =       0.5D0*LFldWght(I)*( ZF/DTANH(Lamda0) + LUnstrLen(I) )
         
         LFairxiTe  = LFairHTe(I)*COSPhi
         LFairyiTe  = LFairHTe(I)*SINPhi
         LFairziTe  = LFairVTe(I)

         F_ILine(1) = F_ILine(1) -            LFairxiTe
         F_ILine(2) = F_ILine(2) -            LFairyiTe
         F_ILine(3) = F_ILine(3) -            LFairziTe
         F_ILine(4) = F_ILine(4) - LFairyiRef*LFairziTe + LFairziRef*LFairyiTe
         F_ILine(5) = F_ILine(5) - LFairziRef*LFairxiTe + LFairxiRef*LFairziTe
         F_ILine(6) = F_ILine(6) - LFairxiRef*LFairyiTe + LFairyiRef*LFairxiTe
         
         F_Lines(1) = F_Lines(1) + F_ILine(1)
         F_Lines(2) = F_Lines(2) + F_ILine(2)
         F_Lines(3) = F_Lines(3) + F_ILine(3)         
         F_Lines(4) = F_Lines(4) + F_ILine(4)
         F_Lines(5) = F_Lines(5) + F_ILine(5)         
         F_Lines(6) = F_Lines(6) + F_ILine(6)
         
     ENDDO             ! I - All mooring lines
               
   ELSE              ! Standard quasi-static.
     
      ! Get the transformation matrix, TransMat, from the inertial frame to the
      !   tower base / platform coordinate system:
    
      CALL LargRotTrans ( X(4), X(5), X(6), TransMat )
       
      !DO I = 1,NumLines ! Loop through all mooring lines
    !Parallel calculation zhu
      !WRITE(88,*)!'Output\ZHUWMOOR1.txt'
      !WRITE(89,*)!'Output\ZHUWMOOR2.txt'
      !WRITE(90,*)!'Output\ZHUWMOOR3.txt'
      !$omp parallel num_threads(Numlines) &
      !$omp private(NumThread,LFairxiRef,LFairyiRef,LFairziRef,LFairxi,LFairyi,LFairzi,XF,ZF,COSPhi,SINPhi) &
      !$omp private(LFairxiTe,LFairyiTe,LFairziTe,F_ILine,I)
      NumThread = omp_get_thread_num()  !ZHU
      I = NumThread+1                   !ZHU
      
      
          
      ! First initialize the total load contribution from all mooring lines to
      !   zero:

      F_ILine(:) = 0.0

      ! Transform the fairlead location from the platform to the inertial
      !   frame coordinate system:
      ! NOTE: TransMat^T = TransMat^-1 where ^T = matrix transpose and ^-1 =
      !       matrix inverse.

         LFairxiRef = TransMat(1,1)*LFairxt(I) + TransMat(2,1)*LFairyt(I) + TransMat(3,1)*LFairzt(I)
         LFairyiRef = TransMat(1,2)*LFairxt(I) + TransMat(2,2)*LFairyt(I) + TransMat(3,2)*LFairzt(I)
         LFairziRef = TransMat(1,3)*LFairxt(I) + TransMat(2,3)*LFairyt(I) + TransMat(3,3)*LFairzt(I)
         
         !WRITE(*,*)'LFair    ',LFairxiRef,'    ',LFairyiRef,'    ',LFairziRef

         LFairxi    = X(1) + LFairxiRef
         LFairyi    = X(2) + LFairyiRef
         LFairzi    = X(3) + LFairziRef

      ! Transform the fairlead location from the inertial frame coordinate system
      !   to the local coordinate system of the current line (this coordinate
      !   system lies at the current anchor, Z being vertical, and X directed from
      !   current anchor to the current fairlead).  Also, compute the orientation
      !   of this local coordinate system:

         XF         = DSQRT( ( LFairxi - LAnchxi(I) )**2.0D0 + ( LFairyi - LAnchyi(I) )**2.0D0 )
         ZF         =          LFairzi - LAnchzi(I)
         IF ( XF == 0.0D0 )  THEN  ! .TRUE. if the current mooring line is exactly vertical; thus, the solution below is ill-conditioned because the orientation is undefined; so set it such that the tensions and nodal positions are only vertical
            COSPhi  = 0.0D0
            SINPhi  = 0.0D0
            Print *, 'WARNING: Current mooring line is exactly vertical.'
         ELSE                    ! The current mooring line must not be vertical; use simple trigonometry
            COSPhi  =       ( LFairxi - LAnchxi(I) )/XF
            SINPhi  =       ( LFairyi - LAnchyi(I) )/XF
         ENDIF

      ! Solve the analytical, static equilibrium equations for a catenary (or
      !   taut) mooring line with seabed interaction in order to find the
      !   horizontal and vertical tensions at the fairlead in the local coordinate
      !   system of the current line:
      ! NOTE: The values for the horizontal and vertical tensions at the fairlead
      !       from the previous time step are used as the initial guess values at
      !       at this time step (because the LAnchHTe(:) and LAnchVTe(:) arrays
      !       are stored in a module and thus their values are saved from CALL to
      !       CALL).

!       WRITE(53,*)  LEAStff (I)    , LFldWght (I),LSeabedCD(I)
        !WRITE(*,*) XF,'    ',ZF 
        
        !WRITE(87+I,*)'I=    ',I    
        !WRITE(I+87,1000)I,      XF,         ZF,         LUnstrLen(I),   LEAStff(I), LFldWght(I),    LSeabedCD(I)
         
        CALL CatenarySingLine1 ( XF           , ZF          , LUnstrLen(I) , LEAStff (I) , LFldWght (I) ,  160, &
                                 LSeabedCD(I) , LFairHTe(I) , LFairVTe (I) , LAnchHTe(I) , LAnchVTe (I), I)
        WRITE(I+87,*)   LFairHTe(I),'    ',LFairVTe(I)
        
1000            FORMAT(I5,5x,   F20.2,5x,   F20.2,5x,   F20.2,5x,       F20.2,5x,   F20.2,5x,       F20.2,5x,       F20.2,5x,       F20.2,5x    )
      ! Transform the positions of each node on the current line from the local
      !   coordinate system of the current line to the inertial frame coordinate
      !   system:

!         DO J = 1,LineNodes ! Loop through all nodes per line where the line position and tension can be output
!            LNodesPi(I,J,1) = LAnchxi(I) + LNodesX(J)*COSPhi
!            LNodesPi(I,J,2) = LAnchyi(I) + LNodesX(J)*SINPhi
!            LNodesPi(I,J,3) = LAnchzi(I) + LNodesZ(J)
!         ENDDO              ! J - All nodes per line where the line position and tension can be output


      ! Compute the portion of the mooring system load on the platform associated
      !   with the current line:

         LFairxiTe  = LFairHTe(I)*COSPhi
         LFairyiTe  = LFairHTe(I)*SINPhi
         LFairziTe  = LFairVTe(I)

         F_ILine(1) = F_ILine(1) -            LFairxiTe
         F_ILine(2) = F_ILine(2) -            LFairyiTe
         F_ILine(3) = F_ILine(3) -            LFairziTe
         F_ILine(4) = F_ILine(4) - LFairyiRef*LFairziTe + LFairziRef*LFairyiTe
         F_ILine(5) = F_ILine(5) - LFairziRef*LFairxiTe + LFairxiRef*LFairziTe
         F_ILine(6) = F_ILine(6) - LFairxiRef*LFairyiTe + LFairyiRef*LFairxiTe
       
         
         DO J = 1,6
         F_LINE_THREAD(J,I) = F_ILine(J)
         END DO
         !F_Lines(1) = F_Lines(1) + F_ILine(1)
         !F_Lines(2) = F_Lines(2) + F_ILine(2)
         !F_Lines(3) = F_Lines(3) + F_ILine(3) 
         !F_Lines(4) = F_Lines(4) + F_ILine(4)
         !F_Lines(5) = F_Lines(5) + F_ILine(5) 
         !F_Lines(6) = F_Lines(6) + F_ILine(6)
         
      !ENDDO             ! I - All mooring lines
    !$omp end parallel
         DO I =1,NumLines
           DO J = 1,6               
             F_Lines(J) = F_LINE_THREAD(J,I)+F_Lines(J)
           END DO
        END DO
      
   ENDIF
   
    DO I = 1,NumLines    
       LFairTe(I)=DSQRT(LFairHTe(I)**2.D0+LFairVTe(I)**2.D0) 
    ENDDO
    
   
    
  END SUBROUTINE CatenarySystem   
 ! ====================================================================   
      SUBROUTINE CatenarySingLine( XF, ZF, L,  EA, W , &
                                   CB, HF, VF, HA, VA )                                   


         ! This routine solves the analytical, static equilibrium equations
         ! for a catenary (or taut) mooring line with seabed interaction.
         ! Stretching of the line is accounted for, but bending stiffness
         ! is not.  Given the mooring line properties and the fairlead
         ! position relative to the anchor, this routine finds the line
         ! configuration and tensions.  Since the analytical solution
         ! involves two nonlinear equations (XF and  ZF) in two unknowns
         ! (HF and VF), a Newton-Raphson iteration scheme is implemented in
         ! order to solve for the solution.  The values of HF and VF that
         ! are passed into this routine are used as the initial guess in
         ! the iteration.  The Newton-Raphson iteration is only accurate in
         ! double precision, so all of the input/output arguments are
         ! converteds to/from double precision from/to default precision.


      IMPLICIT                        NONE

         ! Local Variables:

      REAL(8), INTENT(IN )      :: CB                                              ! Coefficient of seabed static friction (a negative value indicates no seabed) (-)
      REAL(8)                   :: CBOvrEA                                         ! = CB/EA
      REAL(8), INTENT(IN )      :: EA                                              ! Extensional stiffness of line (N)
      REAL(8), INTENT(OUT)      :: HA                                              ! Effective horizontal tension in line at the anchor   (N)
      REAL(8), INTENT(OUT)      :: HF                                              ! Effective horizontal tension in line at the fairlead (N)
      REAL(8)                   :: HFOvrW                                          ! = HF/W
      REAL(8)                   :: HFOvrWCB                                        ! = HF/W/CB      
      REAL(8), INTENT(IN )      :: L                                               ! Unstretched length of line (meters)
      REAL(8)                   :: Lamda0                                          ! Catenary parameter used to generate the initial guesses of the horizontal and vertical tensions at the fairlead for the Newton-Raphson iteration (-)
      REAL(8)                   :: LMinVFOvrW                                      ! = L - VF/W
      REAL(8)                   :: LOvrEA                                          ! = L/EA
      REAL(8)                   :: Tol                                             ! Convergence tolerance within Newton-Raphson iteration specified as a fraction of tension (-)
      REAL(8), INTENT(OUT)      :: VA                                              ! Effective vertical   tension in line at the anchor   (N)
      REAL(8), INTENT(OUT)      :: VF                                              ! Effective vertical   tension in line at the fairlead (N)
      REAL(8), INTENT(IN )      :: W                                               ! Weight of line in fluid per unit length (N/m)
      REAL(8)                   :: WL                                              ! Total weight of line in fluid (N): W*L
      REAL(8), INTENT(IN )      :: XF                                              ! Horizontal distance between anchor and fairlead (meters)
      REAL(8)                   :: XF2                                             ! = XF*XF
      REAL(8), INTENT(IN )      :: ZF                                              ! Vertical   distance between anchor and fairlead (meters)
      REAL(8)                   :: ZF2                                             ! = ZF*ZF
      REAL(8)                   :: TempXF                                          ! New XF in the interation
      REAL(8)                   :: TempZF                                          ! New ZF in the interation      
      REAL(8)                   :: TempHF                                          ! New HF in the interation
      REAL(8)                   :: TempVF                                          ! New VF in the interation      
      REAL(8)                   :: DXF                                          
      REAL(8)                   :: DZF       
      REAL(8),EXTERNAL          :: DSQR1X
      REAL(8),EXTERNAL          :: INVSINH    
      
      INTEGER(4)                :: I                                               ! Index for counting iterations or looping through line nodes (-)
      INTEGER(4)                :: J                                               ! Index for counting iterations or looping through line nodes (-)

      LOGICAL              :: IfRestbed  = .FALSE.                              ! When .TRUE., indicates that we're on the first warning.
  


    ! Initialize some commonly used terms that don't depend on the iteration:
   
      WL      =          W  *L
      LOvrEA  =          L  /EA
      CBOvrEA =          CB /EA
      Tol     =          1.E-3  

      ! Generate the initial guess values for the horizontal and vertical tensions
      !   at the fairlead in the Newton-Raphson iteration for the catenary mooring
      !   line solution.  Use starting values documented in: Peyrot, Alain H. and
      !   Goulois, A. M., "Analysis Of Cable Structures," Computers & Structures,
      !   Vol. 10, 1979, pp. 805-813:

      XF2 = XF*XF
      ZF2 = ZF*ZF

      IF ( DABS( XF ) <= Tol )  THEN ! .TRUE. if the current mooring line is taut
          Lamda0 = 1.E+6    
      ELSEIF ( L <= DSQRT( XF2 + ZF2 ) )  THEN ! .TRUE. if the current mooring line is taut
          Lamda0 = 0.2D0
      ELSE                                ! The current mooring line must be slack and not vertical
          Lamda0 = DSQRT( 3.0D0*( ( L*L - ZF2 )/XF2 - 1.0D0 ) )
      ENDIF

      HF  = MAX( DABS( 0.5D0*W*  XF/     Lamda0      ), Tol )   ! As above, set the lower limit of the guess value of HF to the tolerance
      VF  =            0.5D0*W*( ZF/DTANH(Lamda0) + L )
       WRITE(53,*) HF,VF 
      TempXF = XF + 1.D0
      TempZF = ZF + 2.D0
      
      DXF = TempXF - XF
      DZF = TempZF - ZF  
      
      TempHF = HF * ( 1.D0 - DXF /XF ) + 1.D0
      TempVF = VF * ( 1.D0 - DZF /ZF ) + 1.D0      
 
    ! Start iteration for the Equation without seabed interaction:
    
      J = 0
    
      DO WHILE ( DABS( TempHF-HF ) > Tol .OR. DABS( TempVF-VF ) > Tol )
          
    ! Initialize some commonly used terms that depend on HF and VF:

        HFOvrW             =      TempHF/W
      
        TempXF  =  HFOvrW *( InvSinh( TempVF, TempHF ) - InvSinh( TempVF - WL, TempHF ) ) + TempHF * LOvrEA
        TempZF  =  HFOvrW * Dsqr1x( TempVF, TempHF, WL ) + LOvrEA * ( TempVF - 0.5D0*WL )     
      
        HF = TempHF
        VF = TempVF
        
        DXF = TempXF - XF
        DZF = TempZF - ZF      

        TempHF = HF * ( 1.D0 - DXF /XF )
        TempVF = VF * ( 1.D0 - DZF /ZF )
        
        J = J + 1
!       WRITE(53,*) TempXF,TempZF          
      ENDDO
      
     ! Fairlead tensions:
     
        HF = TempHF
        VF = TempVF
        
     ! Anchor tensions:

        HA = HF
        VA = VF - WL
         
    ! Check if the line rests on the seabed     
        
       IF ( TempVF >= WL ) THEN
           
         IfRestbed = .FALSE.
         
       ELSE
           
         IfRestbed = .TRUE.
         
       ENDIF
         PRINT*, 'DODO'      
    ! Check if the line rests on the seabed,  solve another catenary equation        
       
     IF ( IfRestbed ) THEN
           
       HF  = MAX( DABS( 0.5D0*W*  XF/     Lamda0      ), Tol )   ! As above, set the lower limit of the guess value of HF to the tolerance
       VF  =            0.5D0*W*( ZF/DTANH(Lamda0) + L )

       TempXF = XF + 1.D0
       TempZF = ZF + 1.D0
      
       DXF = TempXF - XF
       DZF = TempZF - ZF  
      
       TempHF = HF * ( 1.D0 - DXF /XF ) + 1.D0
       TempVF = VF * ( 1.D0 - DZF /ZF ) + 1.D0               
      
    ! Start iteration for the Equation with seabed interaction:
    
       J = 0

       DO WHILE ( DABS( TempHF-HF ) > Tol .OR. DABS( TempVF-VF ) > Tol )
          
    ! Initialize some commonly used terms that depend on HF and VF:
    
        LMinVFOvrW         = L  - TempVF/W
        HFOvrW             =      TempHF/W
        HFOvrWCB           =      TempHF/W/CB        

        TempXF  =  LMinVFOvrW + HFOvrW * InvSinh( TempVF, TempHF ) + TempHF * LOvrEA  &
                              + 0.5D0* CBOvrEA *W *( -LMinVFOvrW**2 + ( LMinVFOvrW-HFOvrWCB )*MAX( LMinVFOvrW-HFOvrWCB, 0.D0) )
        
        TempZF  =  HFOvrW * Dsqr1x( TempVF, TempHF, WL ) + LOvrEA * ( TempVF - 0.5D0*WL )     
      
        HF = TempHF
        VF = TempVF
        
        DXF = TempXF - XF
        DZF = TempZF - ZF      

        TempHF = HF * ( 1.D0 - DXF /XF )
        TempVF = VF * ( 1.D0 - DZF /ZF )
        
        J = J + 1
        
      ENDDO
      
     ! Fairlead tensions:
     
        HF = TempHF
        VF = TempVF
        
     ! Anchor tensions:

         HA = MAX( HF - CB *W *LMinVFOvrW, 0.D0)
         VA = 0.D0
        
    ENDIF         


      RETURN
      END SUBROUTINE CatenarySingLine
       
 ! ========================================================
!  Inverse sinh function    

      REAL*8 FUNCTION InvSinh(avar,bvar)
      
      IMPLICIT NONE
      REAL*8 x,avar,bvar
      
      x=avar/bvar
      InvSinh=DLOG(x+DSQRT(1.d0+x*x))
      
      RETURN
      END        

!  Difference of two square functions   

      REAL*8 FUNCTION Dsqr1x(avar,bvar,cvar)
      
      IMPLICIT NONE
      REAL*8 x1,x2,avar,bvar,cvar
      
      x1= avar/bvar
      x2=(avar-cvar)/bvar      
      Dsqr1x=DSQRT(1.d0+x1*x1)-DSQRT(1.d0+x2*x2)
      
      RETURN
      END        
 ! ====================================================================   
      SUBROUTINE CatenarySingLine1( XF, ZF, L,  EA, W , N, &
                                   CB, HF, VF, HA, VA, NoThread )                                   


         ! This routine solves the analytical, static equilibrium equations
         ! for a catenary (or taut) mooring line with seabed interaction.
         ! Stretching of the line is accounted for, but bending stiffness
         ! is not.  Given the mooring line properties and the fairlead
         ! position relative to the anchor, this routine finds the line
         ! configuration and tensions.  Since the analytical solution
         ! involves two nonlinear equations (XF and  ZF) in two unknowns
         ! (HF and VF), a Newton-Raphson iteration scheme is implemented in
         ! order to solve for the solution.  The values of HF and VF that
         ! are passed into this routine are used as the initial guess in
         ! the iteration.  The Newton-Raphson iteration is only accurate in
         ! double precision, so all of the input/output arguments are
         ! converteds to/from double precision from/to default precision.


      IMPLICIT                        NONE


         ! Passed Variables:

      INTEGER(4), INTENT(IN   )    :: N                                               ! Number of nodes where the line position and tension can be output (-)

      REAL(8), INTENT(IN   )    :: CB                                            ! Coefficient of seabed static friction drag (a negative value indicates no seabed) (-)
      REAL(8), INTENT(IN   )    :: EA                                            ! Extensional stiffness of line (N)
      REAL(8), INTENT(  OUT)    :: HA                                            ! Effective horizontal tension in line at the anchor   (N)
      REAL(8), INTENT(  OUT)    :: HF                                            ! Effective horizontal tension in line at the fairlead (N)
      REAL(8), INTENT(IN   )    :: L                                             ! Unstretched length of line (meters)
      REAL(8)  s      (N)                                    ! Unstretched arc distance along line from anchor to each node where the line position and tension can be output (meters)
      REAL(8)  Te     (N)                                    ! Effective line tensions at each node (N)
      REAL(8)  Tol                                           ! Convergence tolerance within Newton-Raphson iteration specified as a fraction of tension (-)
      REAL(8), INTENT(  OUT)    :: VA                                            ! Effective vertical   tension in line at the anchor   (N)
      REAL(8), INTENT(  OUT)    :: VF                                            ! Effective vertical   tension in line at the fairlead (N)
      REAL(8), INTENT(IN   )    :: W                                             ! Weight of line in fluid per unit length (N/m)
      REAL(8)  X      (N)                                    ! Horizontal locations of each line node relative to the anchor (meters)
      REAL(8)  XF                                            ! Horizontal distance between anchor and fairlead (meters)
      REAL(8)  Z      (N)                                    ! Vertical   locations of each line node relative to the anchor (meters)
      REAL(8)  ZF                                            ! Vertical   distance between anchor and fairlead (meters)
      INTEGER, INTENT(IN)     :: NoThread


         ! Local Variables:

      REAL(8)                   :: CBOvrEA                                         ! = CB/EA
      REAL(8)                   :: DET                                             ! Determinant of the Jacobian matrix (m^2/N^2)
      REAL(8)                   :: dHF                                             ! Increment in HF predicted by Newton-Raphson (N)
      REAL(8)                   :: dVF                                             ! Increment in VF predicted by Newton-Raphson (N)
      REAL(8)                   :: dXFdHF                                          ! Partial derivative of the calculated horizontal distance with respect to the horizontal fairlead tension (m/N): dXF(HF,VF)/dHF
      REAL(8)                   :: dXFdVF                                          ! Partial derivative of the calculated horizontal distance with respect to the vertical   fairlead tension (m/N): dXF(HF,VF)/dVF
      REAL(8)                   :: dZFdHF                                          ! Partial derivative of the calculated vertical   distance with respect to the horizontal fairlead tension (m/N): dZF(HF,VF)/dHF
      REAL(8)                   :: dZFdVF                                          ! Partial derivative of the calculated vertical   distance with respect to the vertical   fairlead tension (m/N): dZF(HF,VF)/dVF
      REAL(8)                   :: EXF                                             ! Error function between calculated and known horizontal distance (meters): XF(HF,VF) - XF
      REAL(8)                   :: EZF                                             ! Error function between calculated and known vertical   distance (meters): ZF(HF,VF) - ZF
      REAL(8)                   :: HFOvrW                                          ! = HF/W
      REAL(8)                   :: HFOvrWEA                                        ! = HF/WEA
      REAL(8)                   :: Lamda0                                          ! Catenary parameter used to generate the initial guesses of the horizontal and vertical tensions at the fairlead for the Newton-Raphson iteration (-)
      REAL(8)                   :: LMax                                            ! Maximum stretched length of the line with seabed interaction beyond which the line would have to double-back on itself; here the line forms an "L" between the anchor and fairlead (i.e. it is horizontal along the seabed from the anchor, then vertical to the fairlead) (meters)
      REAL(8)                   :: LMinVFOvrW                                      ! = L - VF/W
      REAL(8)                   :: LOvrEA                                          ! = L/EA
      REAL(8)                   :: sOvrEA                                          ! = s(I)/EA
      REAL(8)                   :: SQRT1VFOvrHF2                                   ! = SQRT( 1.0D0 + VFOvrHF2      )
      REAL(8)                   :: SQRT1VFMinWLOvrHF2                              ! = SQRT( 1.0D0 + VFMinWLOvrHF2 )
      REAL(8)                   :: SQRT1VFMinWLsOvrHF2                             ! = SQRT( 1.0D0 + VFMinWLsOvrHF*VFMinWLsOvrHF )
      REAL(8)                   :: VFMinWL                                         ! = VF - WL
      REAL(8)                   :: VFMinWLOvrHF                                    ! = VFMinWL/HF
      REAL(8)                   :: VFMinWLOvrHF2                                   ! = VFMinWLOvrHF*VFMinWLOvrHF
      REAL(8)                   :: VFMinWLs                                        ! = VFMinWL + Ws
      REAL(8)                   :: VFMinWLsOvrHF                                   ! = VFMinWLs/HF
      REAL(8)                   :: VFOvrHF                                         ! = VF/HF
      REAL(8)                   :: VFOvrHF2                                        ! = VFOvrHF*VFOvrHF
      REAL(8)                   :: VFOvrWEA                                        ! = VF/WEA
      REAL(8)                   :: WEA                                             ! = W*EA
      REAL(8)                   :: WL                                              ! Total weight of line in fluid (N): W*L
      REAL(8)                   :: Ws                                              ! = W*s(I)
      REAL(8)                   :: XF2                                             ! = XF*XF
      REAL(8)                   :: ZF2                                             ! = ZF*ZF

      INTEGER(4)                   :: I                                               ! Index for counting iterations or looping through line nodes (-)
      INTEGER(4)                   :: MaxIter                                         ! Maximum number of Newton-Raphson iterations possible before giving up (-)

      LOGICAL                      :: FirstIter                                       ! Flag to determine whether or not this is the first time through the Newton-Raphson interation (flag)



         ! Abort when there is no solution or when the only possible solution is
         !   illogical:

      IF (    Tol <= 0.0D0 )  THEN   ! .TRUE. when the convergence tolerance is specified incorrectly

!         Print*, ' Convergence tolerance must be greater than zero in routine Catenary().'  


      ELSEIF ( XF <  0.0D0 )  THEN   ! .TRUE. only when the local coordinate system is not computed correctly

         Print*,  ' The horizontal distance between an anchor and its'// &
                      ' fairlead must not be less than zero in routine Catenary().' 


      ELSEIF ( ZF <  0.0D0 )  THEN   ! .TRUE. if the fairlead has passed below its anchor

         Print*, ' A fairlead has passed below its anchor.' 


      ELSEIF ( L  <= 0.0D0 )  THEN   ! .TRUE. when the unstretched line length is specified incorrectly

         Print*, ' Unstretched length of line must be greater than zero in routine Catenary().' 


      ELSEIF ( EA <= 0.0D0 )  THEN   ! .TRUE. when the unstretched line length is specified incorrectly

         Print*, ' Extensional stiffness of line must be greater than zero in routine Catenary().' 


      ELSEIF ( W  == 0.0D0 )  THEN   ! .TRUE. when the weight of the line in fluid is zero so that catenary solution is ill-conditioned

         Print*, ' The weight of the line in fluid must not be zero. '// &
                      ' Routine Catenary() cannot solve quasi-static mooring line solution.' 


      ELSEIF ( W  >  0.0D0 )  THEN   ! .TRUE. when the line will sink in fluid

         LMax      = XF - EA/W + DSQRT( (EA/W)*(EA/W) + 2.0D0*ZF*EA/W )  ! Compute the maximum stretched length of the line with seabed interaction beyond which the line would have to double-back on itself; here the line forms an "L" between the anchor and fairlead (i.e. it is horizontal along the seabed from the anchor, then vertical to the fairlead)
        
         IF ( ( L  >=  LMax   ) .AND. ( CB >= 0.0D0 ) )  &  ! .TRUE. if the line is as long or longer than its maximum possible value with seabed interaction
            Print*, ' Unstretched mooring line length too large. '// &
                         ' Routine Catenary() cannot solve quasi-static mooring line solution.' 
!         write(*,*) 'LMax',LMax,L   

      ENDIF



         ! Initialize some commonly used terms that don't depend on the iteration:
      Tol     =          1.E-6
      WL      =          W  *L
      WEA     =          W  *EA
      LOvrEA  =          L  /EA
      CBOvrEA =          CB /EA
      MaxIter = INT(1.0D0/Tol)   ! Smaller tolerances may take more iterations, so choose a maximum inversely proportional to the tolerance



         ! To avoid an ill-conditioned situation, ensure that the initial guess for
         !   HF is not less than or equal to zero.  Similarly, avoid the problems
         !   associated with having exactly vertical (so that HF is zero) or exactly
         !   horizontal (so that VF is zero) lines by setting the minimum values
         !   equal to the tolerance.  This prevents us from needing to implement
         !   the known limiting solutions for vertical or horizontal lines (and thus
         !   complicating this routine):

      HF = MAX( HF, Tol )
      XF = MAX( XF, Tol )
      ZF = MAX( ZF, Tol )



         ! Solve the analytical, static equilibrium equations for a catenary (or
         !   taut) mooring line with seabed interaction:

         ! Begin Newton-Raphson iteration:

      I         = 1        ! Initialize iteration counter
      FirstIter = .TRUE.   ! Initialize iteration flag

      DO


         ! Initialize some commonly used terms that depend on HF and VF:

         VFMinWL            = VF - WL
         LMinVFOvrW         = L  - VF/W
         HFOvrW             =      HF/W
         HFOvrWEA           =      HF/WEA
         VFOvrWEA           =      VF/WEA
         VFOvrHF            =      VF/HF
         VFMinWLOvrHF       = VFMinWL/HF
         VFOvrHF2           = VFOvrHF     *VFOvrHF
         VFMinWLOvrHF2      = VFMinWLOvrHF*VFMinWLOvrHF
         SQRT1VFOvrHF2      = DSQRT( 1.0D0 + VFOvrHF2      )
         SQRT1VFMinWLOvrHF2 = DSQRT( 1.0D0 + VFMinWLOvrHF2 )


         ! Compute the error functions (to be zeroed) and the Jacobian matrix
         !   (these depend on the anticipated configuration of the mooring line):

         IF ( ( CB <  0.0D0 ) .OR. ( W  <  0.0D0 ) .OR. ( VFMinWL >  0.0D0 ) )  THEN   ! .TRUE. when no portion of the line      rests on the seabed

            EXF    = (   DLOG( VFOvrHF      +               SQRT1VFOvrHF2      )                                       &
                       - DLOG( VFMinWLOvrHF +               SQRT1VFMinWLOvrHF2 )                                         )*HFOvrW &
                   + LOvrEA*  HF                         - XF
            EZF    = (                                     SQRT1VFOvrHF2                                              &
                       -                                   SQRT1VFMinWLOvrHF2                                           )*HFOvrW &
                   + LOvrEA*( VF - 0.5D0*WL )         - ZF

            dXFdHF = (   DLOG( VFOvrHF      +               SQRT1VFOvrHF2      )                                       &
                       - DLOG( VFMinWLOvrHF +               SQRT1VFMinWLOvrHF2 )                                         )/     W &
                   - (      ( VFOvrHF      + VFOvrHF2     /SQRT1VFOvrHF2      )/( VFOvrHF      + SQRT1VFOvrHF2      ) &
                       -    ( VFMinWLOvrHF + VFMinWLOvrHF2/SQRT1VFMinWLOvrHF2 )/( VFMinWLOvrHF + SQRT1VFMinWLOvrHF2 )   )/     W &
                   + LOvrEA
            dXFdVF = (      ( 1.0D0     + VFOvrHF      /SQRT1VFOvrHF2      )/( VFOvrHF      + SQRT1VFOvrHF2      ) &
                       -    ( 1.0D0     + VFMinWLOvrHF /SQRT1VFMinWLOvrHF2 )/( VFMinWLOvrHF + SQRT1VFMinWLOvrHF2 )   )/     W
            dZFdHF = (                                     SQRT1VFOvrHF2                                              &
                       -                                   SQRT1VFMinWLOvrHF2                                           )/     W &
                   - (                       VFOvrHF2     /SQRT1VFOvrHF2                                              &
                       -                     VFMinWLOvrHF2/SQRT1VFMinWLOvrHF2                                           )/     W
            dZFdVF = (                       VFOvrHF      /SQRT1VFOvrHF2                                              &
                       -                     VFMinWLOvrHF /SQRT1VFMinWLOvrHF2                                           )/     W &
                   + LOvrEA


         ELSEIF (                                           -CB*VFMinWL <  HF         )  THEN   ! .TRUE. when a  portion of the line      rests on the seabed and the anchor tension is nonzero

            EXF    =     DLOG( VFOvrHF      +               SQRT1VFOvrHF2      )                                          *HFOvrW &
                   - 0.5D0*CBOvrEA*W*  LMinVFOvrW*LMinVFOvrW                                                                  &
                   + LOvrEA*  HF           + LMinVFOvrW  - XF
            EZF    = (                                     SQRT1VFOvrHF2                                   - 1.0D0   )*HFOvrW &
                   + 0.5D0*VF*VFOvrWEA                - ZF

            dXFdHF =     DLOG( VFOvrHF      +               SQRT1VFOvrHF2      )                                          /     W &
                   - (      ( VFOvrHF      + VFOvrHF2     /SQRT1VFOvrHF2      )/( VFOvrHF      + SQRT1VFOvrHF2        ) )/     W &
                   + LOvrEA
            dXFdVF = (      ( 1.0D0     + VFOvrHF      /SQRT1VFOvrHF2      )/( VFOvrHF      + SQRT1VFOvrHF2        ) )/     W &
                   + CBOvrEA*LMinVFOvrW - 1.0D0/W
            dZFdHF = (                                     SQRT1VFOvrHF2                                   - 1.0D0 &
                       -                     VFOvrHF2     /SQRT1VFOvrHF2                                                )/     W
            dZFdVF = (                       VFOvrHF      /SQRT1VFOvrHF2                                                )/     W &
                   + VFOvrWEA


         ELSE                                                ! 0.0D0 <  HF  <= -CB*VFMinWL   !             A  portion of the line must rest  on the seabed and the anchor tension is    zero

            EXF    =     DLOG( VFOvrHF      +               SQRT1VFOvrHF2      )                                          *HFOvrW &
                   - 0.5D0*CBOvrEA*W*( LMinVFOvrW*LMinVFOvrW - ( LMinVFOvrW - HFOvrW/CB )*( LMinVFOvrW - HFOvrW/CB ) )        &
                   + LOvrEA*  HF           + LMinVFOvrW  - XF
            EZF    = (                                     SQRT1VFOvrHF2                                   - 1.0D0   )*HFOvrW &
                   + 0.5D0*VF*VFOvrWEA                - ZF

            dXFdHF =     DLOG( VFOvrHF      +               SQRT1VFOvrHF2      )                                          /     W &
                   - (      ( VFOvrHF      + VFOvrHF2     /SQRT1VFOvrHF2      )/( VFOvrHF      + SQRT1VFOvrHF2      )   )/     W &
                   + LOvrEA - ( LMinVFOvrW - HFOvrW/CB )/EA
            dXFdVF = (      ( 1.0D0     + VFOvrHF      /SQRT1VFOvrHF2      )/( VFOvrHF      + SQRT1VFOvrHF2      )   )/     W &
                   + HFOvrWEA           - 1.0D0/W
            dZFdHF = (                                     SQRT1VFOvrHF2                                   - 1.0D0 &
                       -                     VFOvrHF2     /SQRT1VFOvrHF2                                                )/     W
            dZFdVF = (                       VFOvrHF      /SQRT1VFOvrHF2                                                )/     W &
                   + VFOvrWEA


         ENDIF


         ! Compute the determinant of the Jacobian matrix and the incremental
         !   tensions predicted by Newton-Raphson:

         DET = dXFdHF*dZFdVF - dXFdVF*dZFdHF

         dHF = ( -dZFdVF*EXF + dXFdVF*EZF )/DET    ! This is the incremental change in horizontal tension at the fairlead as predicted by Newton-Raphson
         dVF = (  dZFdHF*EXF - dXFdHF*EZF )/DET    ! This is the incremental change in vertical   tension at the fairlead as predicted by Newton-Raphson

         dHF = dHF*( 1.0D0 - Tol*I )            ! Reduce dHF by factor (between 1 at I = 1 and 0 at I = MaxIter) that reduces linearly with iteration count to ensure that we converge on a solution even in the case were we obtain a nonconvergent cycle about the correct solution (this happens, for example, if we jump to quickly between a taut and slack catenary)
         dVF = dVF*( 1.0D0 - Tol*I )            ! Reduce dHF by factor (between 1 at I = 1 and 0 at I = MaxIter) that reduces linearly with iteration count to ensure that we converge on a solution even in the case were we obtain a nonconvergent cycle about the correct solution (this happens, for example, if we jump to quickly between a taut and slack catenary)

         dHF = MAX( dHF, ( Tol - 1.0D0 )*HF )   ! To avoid an ill-conditioned situation, make sure HF does not go less than or equal to zero by having a lower limit of Tol*HF [NOTE: the value of dHF = ( Tol - 1.0D0 )*HF comes from: HF = HF + dHF = Tol*HF when dHF = ( Tol - 1.0D0 )*HF]


         ! Check if we have converged on a solution, or restart the iteration, or
         !   Abort if we cannot find a solution:

         IF ( ( DABS(dHF) <= DABS(Tol*HF) ) .AND. ( DABS(dVF) <= DABS(Tol*VF) ) )  THEN ! .TRUE. if we have converged; stop iterating! [The converge tolerance, Tol, is a fraction of tension]

            EXIT


         ELSEIF ( ( I == MaxIter )        .AND. (       FirstIter         ) )  THEN ! .TRUE. if we've iterated MaxIter-times for the first time;

         ! Perhaps we failed to converge because our initial guess was too far off.
         !   (This could happen, for example, while linearizing a model via large
         !   pertubations in the DOFs.)  Instead, use starting values documented in:
         !   Peyrot, Alain H. and Goulois, A. M., "Analysis Of Cable Structures,"
         !   Computers & Structures, Vol. 10, 1979, pp. 805-813:
         ! NOTE: We don't need to check if the current mooring line is exactly
         !       vertical (i.e., we don't need to check if XF == 0.0), because XF is
         !       limited by the tolerance above.

            XF2 = XF*XF
            ZF2 = ZF*ZF

            IF ( L <= DSQRT( XF2 + ZF2 ) )  THEN ! .TRUE. if the current mooring line is taut
               Lamda0 = 0.2D0
            ELSE                                ! The current mooring line must be slack and not vertical
               Lamda0 = DSQRT( 3.0D0*( ( L*L - ZF2 )/XF2 - 1.0D0 ) )
            ENDIF

            HF  = MAX( DABS( 0.5D0*W*  XF/     Lamda0      ), Tol )   ! As above, set the lower limit of the guess value of HF to the tolerance
            VF  =           0.5D0*W*( ZF/DTANH(Lamda0) + L )


         ! Restart Newton-Raphson iteration:

            I         = 0
            FirstIter = .FALSE.
            dHF       = 0.0D0
            dVF       = 0.0D0


         ELSEIF ( ( I == MaxIter )        .AND. ( .NOT. FirstIter         ) )  THEN ! .TRUE. if we've iterated as much as we can take without finding a solution; Abort

            Print*,  ' Iteration not convergent. '// &
                         ' Routine Catenary() cannot solve quasi-static mooring line solution.'  


         ENDIF


         ! Increment fairlead tensions and iteration counter so we can try again:

         HF = HF + dHF
         VF = VF + dVF

         I  = I  + 1


      ENDDO



         ! We have found a solution for the tensions at the fairlead!

         ! Now compute the tensions at the anchor and the line position and tension
         !   at each node (again, these depend on the configuration of the mooring
         !   line):

      IF ( ( CB <  0.0D0 ) .OR. ( W  <  0.0D0 ) .OR. ( VFMinWL >  0.0D0 ) )  THEN   ! .TRUE. when no portion of the line      rests on the seabed

         ! Anchor tensions:

         HA = HF
         VA = VFMinWL


         ! Line position and tension at each node:

         DO I = 1,N  ! Loop through all nodes where the line position and tension are to be computed

!BJJ START of proposed change v7.00.01a-bjj
!bjj: line is too long
!rm            IF ( ( s(I) <  0.0D0 ) .OR. ( s(I) >  L ) )  &
!rm               CALL ProgAbort ( ' All line nodes must be located between the anchor and fairlead (inclusive) in routine Catenary().' )
             IF ( ( s(I) <  0.0D0 ) .OR. ( s(I) >  L ) )  THEN
!               Print*,  ' All line nodes must be located between the anchor ' &
!                              //'and fairlead (inclusive) in routine Catenary().' 
             END IF
!bjj end of proposed change            

            Ws                  = W       *s(I)                                  ! Initialize
            VFMinWLs            = VFMinWL + Ws                                   ! some commonly
            VFMinWLsOvrHF       = VFMinWLs/HF                                    ! used terms
            sOvrEA              = s(I)    /EA                                    ! that depend
            SQRT1VFMinWLsOvrHF2 = DSQRT( 1.0D0 + VFMinWLsOvrHF*VFMinWLsOvrHF ) ! on s(I)

            X (I)    = (   DLOG( VFMinWLsOvrHF + SQRT1VFMinWLsOvrHF2 ) &
                         - DLOG( VFMinWLOvrHF  + SQRT1VFMinWLOvrHF2  )   )*HFOvrW                     &
                     + sOvrEA*  HF
            Z (I)    = (                        SQRT1VFMinWLsOvrHF2   &
                         -                      SQRT1VFMinWLOvrHF2      )*HFOvrW                     &
                     + sOvrEA*(         VFMinWL + 0.5D0*Ws    )
            Te(I)    = DSQRT( HF*HF +    VFMinWLs*VFMinWLs )

         ENDDO       ! I - All nodes where the line position and tension are to be computed


      ELSEIF (                                           -CB*VFMinWL <  HF         )  THEN   ! .TRUE. when a  portion of the line      rests on the seabed and the anchor tension is nonzero

         ! Anchor tensions:

         HA = HF + CB*VFMinWL
         VA = 0.0D0


         ! Line position and tension at each node:

         DO I = 1,N  ! Loop through all nodes where the line position and tension are to be computed

!BJJ START of proposed change v7.00.01a-bjj
!bjj: line is too long
!rm            IF ( ( s(I) <  0.0D0 ) .OR. ( s(I) >  L ) )  &
!rm               CALL ProgAbort ( ' All line nodes must be located between the anchor and fairlead (inclusive) in routine Catenary().' )
            IF ( ( s(I) <  0.0D0 ) .OR. ( s(I) >  L ) )  THEN
!               Print*, ' All line nodes must be located between the anchor ' &
!                              //'and fairlead (inclusive) in routine Catenary().' 
            END IF                              
!bjj end of proposed change

            Ws                  = W       *s(I)                                  ! Initialize
            VFMinWLs            = VFMinWL + Ws                                   ! some commonly
            VFMinWLsOvrHF       = VFMinWLs/HF                                    ! used terms
            sOvrEA              = s(I)    /EA                                    ! that depend
            SQRT1VFMinWLsOvrHF2 = DSQRT( 1.0D0 + VFMinWLsOvrHF*VFMinWLsOvrHF ) ! on s(I)

            IF (     s(I) <= LMinVFOvrW             )  THEN ! .TRUE. if this node rests on the seabed and the tension is nonzero

               X (I) = s(I)                                                                          &
!jmj Start of proposed change.  v6.02c-jmj  02-Feb-2007.
!jmj Bug fix: the s(I) in this line should not be here:
!remove6.02c                     + sOvrEA*( HF + CB*VFMinWL + 0.5D0*Ws*CB ) + s(I)
                     + sOvrEA*( HF + CB*VFMinWL + 0.5D0*Ws*CB )
!jmj End of proposed change.  v6.02c-jmj  02-Feb-2007.
               Z (I) = 0.0D0
               Te(I) =       HF    + CB*VFMinWLs

            ELSE           ! LMinVFOvrW < s <= L            !           This node must be above the seabed

               X (I) =     DLOG( VFMinWLsOvrHF + SQRT1VFMinWLsOvrHF2 )    *HFOvrW                     &
                     + sOvrEA*  HF + LMinVFOvrW                    - 0.5D0*CB*VFMinWL*VFMinWL/WEA
               Z (I) = ( - 1.0D0           + SQRT1VFMinWLsOvrHF2     )*HFOvrW                     &
                     + sOvrEA*(         VFMinWL + 0.5D0*Ws    ) + 0.5D0*   VFMinWL*VFMinWL/WEA
               Te(I) = DSQRT( HF*HF +    VFMinWLs*VFMinWLs )

            ENDIF

         ENDDO       ! I - All nodes where the line position and tension are to be computed


      ELSE                                                ! 0.0D0 <  HF  <= -CB*VFMinWL   !             A  portion of the line must rest  on the seabed and the anchor tension is    zero

         ! Anchor tensions:

         HA = 0.0D0
         VA = 0.0D0


         ! Line position and tension at each node:

         DO I = 1,N  ! Loop through all nodes where the line position and tension are to be computed

!BJJ START of proposed change v7.00.01a-bjj
!bjj: line is too long
!RM            IF ( ( s(I) <  0.0D0 ) .OR. ( s(I) >  L ) )  &
!RM               CALL ProgAbort ( ' All line nodes must be located between the anchor and fairlead (inclusive) in routine Catenary().' )
            IF ( ( s(I) <  0.0D0 ) .OR. ( s(I) >  L ) )  THEN
!               Print*,  ' All line nodes must be located between the anchor ' &
!                              //'and fairlead (inclusive) in routine Catenary().' 
            END IF
!BJJ END OF PROPOSED CHANGE

            Ws                  = W       *s(I)                                  ! Initialize
            VFMinWLs            = VFMinWL + Ws                                   ! some commonly
            VFMinWLsOvrHF       = VFMinWLs/HF                                    ! used terms
            sOvrEA              = s(I)    /EA                                    ! that depend
            SQRT1VFMinWLsOvrHF2 = DSQRT( 1.0D0 + VFMinWLsOvrHF*VFMinWLsOvrHF ) ! on s(I)

            IF (     s(I) <= LMinVFOvrW - HFOvrW/CB )  THEN ! .TRUE. if this node rests on the seabed and the tension is    zero

               X (I) = s(I)
               Z (I) = 0.0D0
               Te(I) = 0.0D0

            ELSEIF ( s(I) <= LMinVFOvrW             )  THEN ! .TRUE. if this node rests on the seabed and the tension is nonzero

               X (I) = s(I)                     - ( LMinVFOvrW - 0.5D0*HFOvrW/CB )*HF/EA          &
                     + sOvrEA*( HF + CB*VFMinWL + 0.5D0*Ws*CB ) + 0.5D0*CB*VFMinWL*VFMinWL/WEA
               Z (I) = 0.0D0
               Te(I) =       HF    + CB*VFMinWLs

            ELSE           ! LMinVFOvrW < s <= L            !           This node must be above the seabed

               X (I) =     DLOG( VFMinWLsOvrHF + SQRT1VFMinWLsOvrHF2 )    *HFOvrW                     &
                     + sOvrEA*  HF + LMinVFOvrW - ( LMinVFOvrW - 0.5D0*HFOvrW/CB )*HF/EA
               Z (I) = ( - 1.0D0           + SQRT1VFMinWLsOvrHF2     )*HFOvrW                     &
                     + sOvrEA*(         VFMinWL + 0.5D0*Ws    ) + 0.5D0*   VFMinWL*VFMinWL/WEA
               Te(I) = DSQRT( HF*HF +    VFMinWLs*VFMinWLs )

            ENDIF

         ENDDO       ! I - All nodes where the line position and tension are to be computed


      ENDIF
      
     ! WRITE(87+NoThread,1000) NoThread,     HF,         VF
      
1000 FORMAT(I5,5x,  F20.2,5x,   F20.2)      


      RETURN
      END SUBROUTINE CatenarySingLine1        