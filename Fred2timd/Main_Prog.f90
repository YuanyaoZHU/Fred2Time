
   !-------------------------------------------------------------------------------------------------
   ! MAIN PROGRAM FOR UNSTEADY DYNAMIC COMPUTATION OF MULTI-DUCTED WIND TURBINE SYSTEM...
   !-------------------------------------------------------------------------------------------------
   
      PROGRAM LENSDYN
      
      USE Const_mod
	  USE STKMesh_mod
 	  USE WaveDyn_mod
	  USE Platform_mod
      USE Blade_mod
      USE LGKT4K
      USE CABLE

 
      IMPLICIT NONE
      
         
      INTEGER NREC,NFILE,NSTP,NUMLINES
      INTEGER I,J,K,L,M,N,LEN,IOSTAT,ERR,IRR,IFWKO
 
      REAL*8 WingAng(3),URL(3),Xrl(3),XLC(3),BUFR,PltfSdLnth
      REAL*8 RMAS(6,6),TRMAS(6,6),ZEROS(6),FTAL(6),Force(6)
      
      REAL*8 SIMT,TSP,T,DISP(6),VLCT(6),ERS,ERC
      REAL*8 WMIN,WMAX,SW,DW,SVF,ARST,DLTF,ILTF,PHAW
      REAL*8 TBGN,TEND
      
      
       REAL(8) :: Roll    ! The large rotation about X1, (rad).
       REAL(8) :: Pitch    ! The large rotation about X2, (rad).
       REAL(8) :: Yaw    ! The large rotation about X3, (rad).     
       REAL(8) :: YawAng    ! wind turbine yaw angle, (rad).
       REAL(8) :: TiltAng    ! wind turbine tilt angle, (rad).
       REAL(8) :: ConeAng    ! Cone angle       

       REAL(8) :: VSurge    ! The large rotation about X1, (rad).
       REAL(8) :: VSway    ! The large rotation about X2, (rad).
       REAL(8) :: VHeave    ! The large rotation about X3, (rad).       
       REAL(8) :: VRoll    ! The large rotation about X1, (rad).
       REAL(8) :: VPitch    ! The large rotation about X2, (rad).
       REAL(8) :: VYaw    ! The large rotation about X3, (rad).       
  
       REAL(8) :: Xpc(3)    ! Platfrom mass center in Coordinate System 1 
       REAL(8) :: Xyc(3)    ! Wind turbine yaw center in Coordinate System 2. 
       REAL(8) :: Xrc(3)    ! Wind turbine hub/rotor center in Coordinate System 4. 
       REAL(8) :: ZHub
       REAL(8) :: VHub
       REAL(8) :: WindDrt
       REAL(8) :: TransMat(3,3)
       REAL(8) :: InitialPosition(6)
       REAL(8) :: InitialVerticalForce
       
       REAL(8) :: LDAMP_ZHU(6) !ZHU:linear damping
       REAL(8) :: CDAMP_ZHU(6) !ZHU:coefficiency of linear damping
      
   !-------------------------------------------------------------------------------------------------
   ! Input data files for aerodynamics computation...
   !-------------------------------------------------------------------------------------------------  
   
      OPEN(0, FILE='Input\InputControl.txt',    STATUS='OLD')     
      OPEN(1, FILE='Input\BladeInput.txt',    STATUS='OLD')  
      OPEN(2, FILE='Output\CheckAero.txt',    STATUS='UNKNOWN')      
      OPEN(3, FILE='Output\CoefPower.txt',    STATUS='UNKNOWN')   
      OPEN(4, FILE='Output\CoefTrust.txt',    STATUS='UNKNOWN')      
      OPEN(5, FILE='Output\InflwWind.txt',    STATUS='UNKNOWN')         
      OPEN(6, FILE='Output\RotSpeed.txt',    STATUS='UNKNOWN')          
      OPEN(7, FILE='Output\RotAccer.txt',    STATUS='UNKNOWN')   
      OPEN(8, FILE='Output\Thrust.txt',    STATUS='UNKNOWN')  
      OPEN(9, FILE='Output\Torque.txt',    STATUS='UNKNOWN')        
      OPEN(10, FILE='Input\STKBDMS.txt',     STATUS='OLD') 
      
      OPEN(11, FILE='Input\CYLINDER1.txt',    STATUS='OLD') 
      OPEN(12, FILE='Input\CYLINDER2.txt',    STATUS='OLD') 
      OPEN(13, FILE='Input\DU21_A17.txt',    STATUS='OLD')
      OPEN(14, FILE='Input\DU25_A17.txt',    STATUS='OLD') 
      OPEN(15, FILE='Input\DU30_A17.txt',    STATUS='OLD') 
      OPEN(16, FILE='Input\DU35_A17.txt',    STATUS='OLD')
      OPEN(17, FILE='Input\DU40_A17.txt',    STATUS='OLD') 
      OPEN(18, FILE='Input\NACA64_A17.txt',    STATUS='OLD')    
      OPEN(19, FILE='Output\a0.txt',    STATUS='unknown')
      OPEN(20, FILE='Input\DATMASS.txt',     STATUS='OLD')       
      
   !-------------------------------------------------------------------------------------------------
   ! Input data files for integral dynamic computation...
   !-------------------------------------------------------------------------------------------------   

      OPEN(21, FILE='Input\OAMASS1.txt',    STATUS='OLD')  
      OPEN(22, FILE='Input\OAMASS2.txt',    STATUS='OLD') 
      OPEN(23, FILE='Input\OAMASS3.txt',    STATUS='OLD')       
      OPEN(24, FILE='Input\OAMASS4.txt',    STATUS='OLD')       
      OPEN(25, FILE='Input\OAMASS5.txt',    STATUS='OLD') 
      OPEN(26, FILE='Input\OAMASS6.txt',    STATUS='OLD')    
      
      OPEN(31, FILE='Input\ODAMPING1.txt',    STATUS='OLD')  
      OPEN(32, FILE='Input\ODAMPING2.txt',    STATUS='OLD') 
      OPEN(33, FILE='Input\ODAMPING3.txt',    STATUS='OLD')       
      OPEN(34, FILE='Input\ODAMPING4.txt',    STATUS='OLD')       
      OPEN(35, FILE='Input\ODAMPING5.txt',    STATUS='OLD') 
      OPEN(36, FILE='Input\ODAMPING6.txt',    STATUS='OLD')
      
      OPEN(41,FILE='Input\OEXFOR1.txt', STATUS='OLD') 
      OPEN(42,FILE='Input\OEXFOR2.txt', STATUS='OLD')  
      OPEN(43,FILE='Input\OEXFOR3.txt', STATUS='OLD') 
      OPEN(44,FILE='Input\OEXFOR4.txt', STATUS='OLD') 
      OPEN(45,FILE='Input\OEXFOR5.txt', STATUS='OLD') 
      OPEN(46,FILE='Input\OEXFOR6.txt', STATUS='OLD') 
      OPEN(47,FILE='Input\CURVE.txt', STATUS='OLD') !CSQ
      !OPEN(48,FILE='Input\TORQUECSQ.txt', STATUS='OLD') !CSQ
      
      OPEN(51, FILE='Output\ORTAR1.txt',    STATUS='UNKNOWN')  
      OPEN(52, FILE='Output\ORTAR2.txt',    STATUS='UNKNOWN')        
      OPEN(53, FILE='Output\ORTAR3.txt',    STATUS='UNKNOWN')        
      OPEN(54, FILE='Output\ORTAR4.txt',    STATUS='UNKNOWN')  
      OPEN(55, FILE='Output\ORTAR5.txt',    STATUS='UNKNOWN')        
      OPEN(56, FILE='Output\ORTAR6.txt',    STATUS='UNKNOWN')        
      
      OPEN(60, FILE='Output\OHLTF.txt',    STATUS='UNKNOWN')  
      OPEN(61, FILE='Output\OFWAVE.txt',    STATUS='UNKNOWN')
      OPEN(62, FILE='Output\OWVHT.txt',    STATUS='UNKNOWN')
      OPEN(63, FILE='Output\OAINF.txt',    STATUS='UNKNOWN')  
      OPEN(64, FILE='Output\ODISP.txt',    STATUS='UNKNOWN')
      OPEN(65, FILE='Output\OVLCT.txt',    STATUS='UNKNOWN')       
      OPEN(66, FILE='Output\OACLR.txt',    STATUS='UNKNOWN')      
      
      OPEN(70, FILE='Output\CheckHydro.txt',    STATUS='UNKNOWN')   
      OPEN(71, FILE='Output\RUNGE4.txt',    STATUS='UNKNOWN')       
      OPEN(72, FILE='Output\GKIJ.txt',    STATUS='UNKNOWN')
      OPEN(73, FILE='Output\CATENARY.txt',    STATUS='UNKNOWN')
      OPEN(74, FILE='Output\MATRX.txt',    STATUS='UNKNOWN')      
      OPEN(75, FILE='Output\JSPTR.txt',    STATUS='UNKNOWN')  
      OPEN(76, FILE='Output\PowerOutput.txt',    STATUS='UNKNOWN')
      
      OPEN(80, FILE='Output\ErrorOutput.txt',    STATUS='UNKNOWN')
      OPEN(81, FILE='Output\ZHUHSFORCE.txt',    STATUS='UNKNOWN')
      OPEN(82, FILE='Output\ZHUMOORFORCE.txt',    STATUS='UNKNOWN')
      OPEN(83, FILE='Output\ZHUMOORFORCE0.txt',    STATUS='UNKNOWN')
      OPEN(84, FILE='Output\ZHUWAVEFORCE.txt',    STATUS='UNKNOWN')
      OPEN(85, FILE='Output\ZHUREATFORCE.txt',    STATUS='UNKNOWN')
      OPEN(86, FILE='Output\ZHUDAMPINGFORCE.txt',    STATUS='UNKNOWN')
      OPEN(87, FILE='Output\ZHUWINDFORCE.txt',    STATUS='UNKNOWN')
      OPEN(88, FILE='Output\ZHUWMOOR1.txt',    STATUS='UNKNOWN')
      OPEN(89, FILE='Output\ZHUWMOOR2.txt',    STATUS='UNKNOWN')
      OPEN(90, FILE='Output\ZHUWMOOR3.txt',    STATUS='UNKNOWN')
      OPEN(91, FILE='Output\ZHULDAMPING.txt',    STATUS='UNKNOWN')
      OPEN(92, FILE='Output\ZHUOFFSET.txt',    STATUS='UNKNOWN')
      OPEN(93, FILE='Output\ZHUMASS.txt',    STATUS='UNKNOWN')
      OPEN(94, FILE='Output\MORISON.txt',    STATUS='UNKNOWN')
      OPEN(95, FILE='Output\TESTMORISON.txt',    STATUS='UNKNOWN')
      OPEN(96, FILE='Output\fairxyz.txt',    STATUS='UNKNOWN')
      CALL CPU_TIME (TBGN)    
    
   !-------------------------------------------------------------------------------------------------
   ! Read the sample file and determine length of the signal...
   !------------------------------------------------------------------------------------------------- 
      PRINT*
      PRINT*,'Computation starts...'


      NREC=0
      DO WHILE(.NOT.EOF(41))
	   READ(41,*,IOSTAT=ERR)  
	   IF(ERR/=0) EXIT
	   NREC=NREC+1
      ENDDO

      LEN=NREC-1
      ALLOCATE(AW(0:LEN),AMAS(6,6,0:LEN),BDMP(6,6,0:LEN))
      ALLOCATE(EFOR(6,2,0:LEN),HAM(6,6,0:LEN),HBM(6,6,0:LEN))

      DO NFILE=21,26       
      REWIND(NFILE)
      ENDDO 
      
      DO NFILE=31,36  
      REWIND(NFILE)
      ENDDO   
      
      DO NFILE=41,46  
      REWIND(NFILE)
      ENDDO    
      
   !-------------------------------------------------------------------------------------------------
   ! Read the data of added mass and damping coefficients...
   !-------------------------------------------------------------------------------------------------      
      PRINT*
      PRINT*,'Read the data of added mass and damping coefficients...'
      
      DO 100 I=1,6
      DO 100 K=0,LEN
       READ(20+I,*) WK, AW(K),(AMAS(I,J,K),J=1,6)
100   CONTINUE
      
      DO 200 I=1,6
      DO 200 K=0,LEN
       READ(30+I,*) WK, AW(K),(BDMP(I,J,K),J=1,6) 
200   CONTINUE      
      
      WMIN=AW(0)
      WMAX=AW(LEN)
      DW=(WMAX-WMIN)/LEN
     
      DO 300 I=1,6
      DO 300 J=1,6
      DO 300 K=0,LEN
       HAM(I,J,K)=AW(K)*(AMAS(I,J,K)-AMAS(I,J,LEN))     
300   CONTINUE  
      
      DO 400 I=1,6
      DO 400 J=1,6
      DO 400 K=0,LEN
       HBM(I,J,K)=BDMP(I,J,K)
400   CONTINUE 
      
      DO 450 I=1,6
      DO 450 K=0,LEN
       READ(40+I,*) WK, AW(K),(EFOR(I,J,K),J=1,2) 
450   CONTINUE
      
   !-------------------------------------------------------------------------------------------------
   ! Initialization for time domain computation...
   !-------------------------------------------------------------------------------------------------      
      PRINT*
      PRINT*,'Initialize for time domain computation...'
      
      

      READ(0,*) 
      READ(0,*) 
      READ(0,*) RHOW, RHOA
      READ(0,*) 
      READ(0,*) IFDUCT, IFWAKE
      READ(0,*) 
      READ(0,*) IFVISC, IFCOUP
      READ(0,*) 
      READ(0,*) H, HS, IRR
      READ(0,*) 
      READ(0,*) IFWKO, WK
      READ(0,*) 
      READ(0,*) BETA, WindDrt, VHub ,ZHub
      READ(0,*) 
      READ(0,*) MASSRNA, xRNA
      READ(0,*) 
      READ(0,*) SIMT, NSTP
      READ(0,*) 
      READ(0,*) InitialPosition(1),InitialPosition(2),InitialPosition(3),InitialPosition(4),InitialPosition(5),InitialPosition(6)
      READ(0,*) 
      READ(0,*) InitialVerticalForce
      READ(0,*) 
      READ(0,*) CD1,CD2,CD3,CD4,CD5
      READ(0,*) 
      READ(0,*) CDX(1),CDX(2),CDX(3),CDX(4),CDX(5),CDX(6)
      READ(0,*)
      READ(0,*) NUMLINES, RFair, RAnch, DepthAnch, DepthFair
      READ(0,*)
      READ(0,*) LineDiam,LineMassDen,LineUnstrLen,LineSeabedCD,LineEAStff
      ALLOCATE(LAngle(NUMLINES))
      READ(0,*)
      READ(0,*), (LAngle(I),I=1,NUMLINES)
      
      BETA=BETA/180.d0*pi
      WindDrt=WindDrt/180.d0*pi

      
      IF (IFWKO .EQ. 0)  THEN
        IF (H .LE. 0.0D0) THEN
          W1=DSQRT(G*WK)
        ELSE
          W1=DSQRT(G*WK*DTANH(WK*H))
        END IF
        WL=2.D0*PI/WK
        TP=2.D0*PI/W1
      ELSEIF (IFWKO .EQ. 1)  THEN
        W1=WK
        IF(H .LE. 0.0D0) THEN
          WK=W1*W1/G
        ELSE
          CALL WAVECK(W1,H,WK)			  ! Compute wave number
        END IF
        WL=2.D0*PI/WK
        TP=2.D0*PI/W1
      ELSEIF (IFWKO .EQ. 2)  THEN
        TP=WK
        W1=2.D0*PI/TP
        IF(H .LE. 0.0D0) THEN
          WK=W1*W1/G
        ELSE
          CALL WAVECK(W1,H,WK)			  ! Compute wave number
        END IF
        WL=2.D0*PI/WK
      END IF
 
      TSP=SIMT/NSTP
     
      ALLOCATE(WDSP(3,0:NSTP+1),WTSP(3,0:NSTP+1),ACROT(3,0:NSTP+1))
      ALLOCATE(FAERO(3,0:NSTP+1),TAERO(3,0:NSTP+1),MCPWR(3,0:NSTP+1),PITH(3,0:NSTP+1))
      ALLOCATE(COT(3,0:NSTP+1),COP(3,0:NSTP+1))
      
      
      ALLOCATE(RFK(6,6,0:2*NSTP),LTF(6,0:NSTP),ETA(0:NSTP),RK(6,4,0:NSTP),RK1(6,4,0:NSTP))
      ALLOCATE(X(6,0:NSTP),DX(6,0:NSTP),XCM(6,0:NSTP))
      ALLOCATE(FWAVE(6,0:NSTP),FWIND(6,0:NSTP))
      
      CALL MATCRS
      
      READ(10,*)
      READ(10,*)  NSTRE,NPART
        
      ALLOCATE( RELM(NSTRE,3),EPN(NSTRE),DELM(NPART,3),IDL(NPART),NELM(NPART),DCYL(NPART),SCYL(NPART),LCYL(NPART) )          
      ALLOCATE( RSTR(NPART,3),REND(NPART,3),CM(NSTRE),CD(NSTRE) )   
      
      CALL mMESHDA
      
      CALL INITIAL(WDSP(:,0),WTSP(:,0),PITH(:,0)) 
           
   !-------------------------------------------------------------------------------------------------
   ! Compute the retardation function Kij(t) and output...
   !-------------------------------------------------------------------------------------------------   
      PRINT*
      PRINT*,'Computing the retardation function...'  
      
      DO 500 L=0,NSTP*2
       T=TSP/2.D0*L
      DO 500 I=1,6    
      DO 500 J=1,6 

       CALL FILON_KIJ(HBM(I,J,:),LEN,T,WMIN,WMAX,1,RFK(I,J,L)) 
       CALL SUVQ_KIJ(WMAX,HBM(I,J,LEN),T,I,SVF)
      
       RFK(I,J,L)=2.D0/PI*RFK(I,J,L)+SVF

500   CONTINUE

      DO 600 I=1,6      
      DO 600 L=0,NSTP*2
       T=TSP/2.D0*L
       WRITE(50+I,1010) T,(RFK(I,J,L),J=1,6)
600   CONTINUE 
      
      DO 650 L=0,NSTP*2
       T=TSP/2.D0*L
!       WRITE(72,1010) T,RFK(4,4,L),RFK(5,5,L),RFK(6,6,L)
650   CONTINUE 
      
   !-------------------------------------------------------------------------------------------------
   ! Compute the infinite added mass...
   !-------------------------------------------------------------------------------------------------    
      PRINT*
      PRINT*,'Computing the infinite added mass...'       
      
      WRITE(63,142) 
      
      DO 750 I=1,6
      DO 700 J=1,6
          
       CALL FILON_KIJ(RFK(I,J,:),NSTP,WMAX,0.D0,SIMT,2,ARST)
       AINF(I,J)=AMAS(I,J,LEN)+ARST/WMAX
      
700   CONTINUE

       WRITE(63,'(6E14.6)') (AINF(I,J),J=1,6)      
750   CONTINUE
      
   !-------------------------------------------------------------------------------------------------
   ! Compute the linear transfer function and wave excitation force in time series...
   !-------------------------------------------------------------------------------------------------        
      PRINT*  
      PRINT*,'Computing the linear transfer function...'    
      DO 800 L=0,NSTP
       T=TSP*L 
      DO 800 I=1,6    

       CALL FILON_KIJ(EFOR(I,1,:),LEN,T,WMIN,WMAX,1,DLTF) 
       CALL FILON_KIJ(EFOR(I,2,:),LEN,T,WMIN,WMAX,2,ILTF) 
      
       LTF(I,L)=(DLTF-ILTF)/PI
       WRITE(60,1010) T,DLTF,ILTF
800   CONTINUE          
  
      
      PRINT*          
      PRINT*,'Reading the wave elevation data in time series...' 
      
      XLC(1)=0.D0
      XLC(2)=0.D0
      XLC(3)=0.D0
       
      IF (IRR.EQ.1) THEN
       CALL IWHC(XLC,TSP,NSTP,H,HS,TP,BETA,ETA)       
       PRINT*      
       PRINT*,'Computing the wave excitation force in time series...'    
       DO 900 I=1,6
       DO 900 L=0,NSTP
        T=TSP*L
        FWAVE(I,L)=0.D0
       DO 900 J=0,L
        FWAVE(I,L)=FWAVE(I,L)+LTF(I,L-J)*ETA(J)*TSP
900    CONTINUE
      ELSEIF (IRR.EQ.2) THEN
          DO L=0,NSTP
              READ(62,*) T, ETA(L)
          ENDDO
       DO 920 I=1,6
       DO 920 L=0,NSTP
        T=TSP*L
        FWAVE(I,L)=0.D0
       DO 920 J=0,L
        FWAVE(I,L)=FWAVE(I,L)+LTF(I,L-J)*ETA(J)*TSP
920    CONTINUE
             
      ELSE          
          
       DO 910 I=1,6
       DO 910 L=0,NSTP 
        T=TSP*L
        CALL RGWAVEXFC(W1,DW,I,T,FWAVE(I,L),PHAW)
        FWAVE(I,L)=(HS/2.D0)*FWAVE(I,L)
910    CONTINUE
      
      ENDIF
      
      DO 950 L=0,NSTP
        T=TSP*L
        WRITE(61,1010) T,(BUFR(T,TP)*FWAVE(I,L)/(HS/2.D0),I=1,6)
950   CONTINUE     
       
            
   !-------------------------------------------------------------------------------------------------
   ! Solve the motion equation using the 4th-order Runge-Kutta scheme...
   !-------------------------------------------------------------------------------------------------         
      PRINT*      
      PRINT*,'Solving the motion equation starts...' 
      
      
      X=0.D0    !ZHU
      X(1,0) = InitialPosition(1)!*3.1415926/180
      X(2,0) = InitialPosition(2)
      X(3,0) = InitialPosition(3)
      X(4,0) = InitialPosition(4)*3.1415926/180
      X(5,0) = InitialPosition(5)*3.1415926/180
      X(6,0) = InitialPosition(6)*3.1415926/180
      
      DX=0.D0   !ZHU
      FGRVT=0.D0
      FGRVT(3)= -79101100.D0
      FLINE0=0.D0
      FLINE0(3)= 0!-1607000.D0!3.34941E+05 * 40    !0.962826E+06
      FBOYC(3) = 80708100.D0
      
      FRTAR=0.D0      
      RMAS=MATX+AINF 
      ZEROS=0.D0
      DO I=1,    6
        WRITE(93,102) I, (MATX(I,J),  J=1,6)
      END DO
      
      
      WRITE(63,143)
      DO   I=1,   6
       WRITE(63,102) I, (RMAS(I,J),  J=1, 6) 
      END  DO
      
 	  CALL INVERSE(RMAS,TRMAS,6,6) 
      
      WRITE(63,144)
      DO   I=1,   6
       WRITE(63,102) I, (TRMAS(I,J),  J=1, 6) 
      END  DO   
        
      WingAng= 0.d0
      PltfSdLnth= 90.d0 !useless value
      !ZHub= 90.07d0
 
     !---------------------------------------------------
     !  Start of the recycle
     !--------------------------------------------------
      DO 1000 L=0,NSTP-1  !zhu: Start of the recycle

       T=TSP*L 
       WRITE(*,*) '   Current time is :', T, 'sec'
       WRITE(74,*) 'Time step: ',T   
       
   !-------------------------------------------------------------------------------------------------
   ! Compute the angles for wind turbine...
   !-------------------------------------------------------------------------------------------------          

       TiltAng=0.d0
       ConeAng=0.d0
       YawAng=WindDrt-X(6,L)
       WingAng(:)=WingAng(:)+WTSP(:,L)*TSP
       
       DO J=1,6
        XCM(J,L)=X(J,L)+XCM0(J)
       ENDDO       
      
       Roll=  XCM(4,L)
       Pitch= XCM(5,L)
       Yaw=   XCM(6,L)
       
       VSurge= DX(1,L)
       VSway=  DX(2,L)
       VHeave= DX(3,L)          
       VRoll=  DX(4,L)
       VPitch= DX(5,L)
       VYaw=   DX(6,L)
       
   !-------------------------------------------------------------------------------------------------
   ! Compute the unsteady aerodynamic force in time series...
   !-------------------------------------------------------------------------------------------------       
!      PRINT*
!      PRINT*,'Computing the unsteady aerodynamic force...' 
   
       FWIND= 0.D0

       WRITE(2,*)       
       WRITE(2,*) 'T=',T
      
!---- Calculate the first heading wind turbine
       
       Xpc(1)=XCM(1,L)
       Xpc(2)=XCM(2,L)
       Xpc(3)=XCM(3,L)
 
       Xyc(1)= 0
       Xyc(2)= 0
       Xyc(3)=  ZHub-XCM0(3)
       
       Xrc(1)=-1.75d0
       Xrc(2)=  0.d0
       Xrc(3)=  0.d0
       
       Ifheading=.true.
       
       CALL RotorCenterLocation( YawAng, TiltAng,  Roll, Pitch, Yaw, Xpc, Xyc, Xrc, URL )  ! Transform the coordinate of rotor hub center from the 4th system to the 1st system
 
       CALL RelativeVelocity( YawAng, TiltAng, WingAng(3), ConeAng, Roll, Pitch, Yaw, VSurge, VSway, VHeave, VRoll, VPitch, VYaw, Xpc, Xyc, Xrc, URL(3), VHub, WindDrt, NBEM, NBL, RND, COT(1,L), RTB, PltfSdLnth, URL, IFWAKE, IFCOUP, Ifheading, VRltv )
 
       WRITE(2,*)
       WRITE(2,*) 'The wind turbine'
       WRITE(2,*)
       
       CALL AERO_STEP(VRltv,WTSP(1,L),PITH(1,L),FAERO(1,L),TAERO(1,L),MCPWR(1,L),COT(1,L),COP(1,L))
        
       !CALL CONTROL(WTSP(1,L),TAERO(1,L),PITH(1,L),TSP,WTSP(1,L+1),ACROT(1,L),PITH(1,L+1))

       CAll CalForceTorque( WindDrt, FAERO(1,L), TAERO(1,L), URL, Xpc, Force )
       
       FWIND(1,L)= FWIND(1,L)+Force(1)
       FWIND(2,L)= FWIND(2,L)+Force(2)
       FWIND(3,L)= FWIND(3,L)+Force(3)
       FWIND(4,L)= FWIND(4,L)+Force(4)       
       FWIND(5,L)= FWIND(5,L)+Force(5)       
       FWIND(6,L)= FWIND(6,L)+Force(6)
       

!       
!!      PRINT*
!       PRINT*,'Output the aerodynamic simulation results...'      
      
       WRITE(3,1010) T,(COP(1,L))
       WRITE(4,1010) T,(COT(1,L))
       WRITE(5,1010) T,(WDSP(1,L))  
       WRITE(6,1010) T,(WTSP(1,L))     
       WRITE(7,1010) T,(ACROT(1,L))   
       WRITE(8,1010) T,(FAERO(1,L))   
       WRITE(9,1010) T,(TAERO(1,L))     
       WRITE(76,1010) T,(MCPWR(1,L))  
       !WRITE(3,1010) T,(COP(J,L),J=1,3)
       !WRITE(4,1010) T,(COT(J,L),J=1,3)
       !WRITE(5,1010) T,(WDSP(J,L),J=1,3)  
       !WRITE(6,1010) T,(WTSP(J,L),J=1,3)     
       !WRITE(7,1010) T,(ACROT(J,L),J=1,3)   
       !WRITE(8,1010) T,(FAERO(J,L),J=1,3)   
       !WRITE(9,1010) T,(TAERO(J,L),J=1,3)     
       !WRITE(76,1010) T,(MCPWR(J,L),J=1,3)   
       
       
       !--------------ZHU--------------------------------
       
       DO J=1,6
           FORCE_ZHU(J,1)=0
       END DO
       
       STIFF=CRS+LRS
              
       MOTION = LungKuta4K(T,TSP,X(:,L),DX(:,L),TRMAS,ZDAMP,STIFF,FORCE_ZHU,L,XCM0,NUMLINES,NSTP,IRR,InitialVerticalForce)
       
       DO J=1,6
       X(J,L+1) = MOTION(J,1)
       DX(J,L+1)= MOTION(J,2)
       END DO
!--------------------------------------------------------------
!    Delete start
!---------------------------------------------------------------

!-------------------------------------------------------------------------
!    Delete end
!-------------------------------------------------------------------------
       
      WRITE(64,1010) T, (X(I,L),I=1,6)
      WRITE(65,1010) T, (DX(I,L),I=1,6)
      
1000  CONTINUE 
      
      
      CALL CPU_TIME (TEND) 
	  WRITE(6,1600)  TEND-TBGN
	  WRITE(63,1600)  TEND-TBGN
      
     
1600  FORMAT(//, 'Time of operation was',F12.3, 'seconds')      
142   FORMAT(/,10x,'  MATRIX of Added MASS') 
143   FORMAT(/,10x,'  ORIGINAL  MATRIX of MASS')   
144   FORMAT(/,10x,'  INVERSE  MATRIX of MASS')          
102   FORMAT('M',I1,'J',3x,6(2x,E10.3))       
1010  FORMAT(F10.2,1x,6E14.6)
1011  FORMAT(F20.2,5x,F20.2,5x,F20.2,5x,F20.2,5x,F20.2,5x,F0.2)
1012  FORMAT(F20.2,5x,F20.2)    
      
      END PROGRAM LENSDYN
      
