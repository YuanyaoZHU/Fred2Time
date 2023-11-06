 
   !-------------------------------------------------------------------------------------------------
   ! Aerodynamic calculation for the rotor at every time step...
   !-------------------------------------------------------------------------------------------------            
        
      SUBROUTINE AERO_STEP(VWND,WROT,PTAG,THRUST,MTB,POWER,CT,CP)
      
      USE Const_mod     
      USE Blade_mod  
      IMPLICIT NONE   
      
      !INTEGER IEL,IBL,J
      INTEGER*4 :: line=21, column=20
      REAL*8 :: THRUST !CSQ
      CHARACTER(LEN=20) :: filename !CSQ filename的字符长度定义
      REAL*8  VWND(3),WROT,PTAG,FTB,MTB,CT,CP,power
      REAL*8  VRELWND(3)
      REAL*8  Vlct(3), Roll, Pitch, Yaw, VSurge, VSway, VHeave, VRoll, VPitch, VYaw, Rcg, Rh
      !REAL*8  VRt(3), YawAng, TiltAng, WingAng, ConeAng, Xpc(3), Xyc(3), Xrc(3), ZHub, VHub, WindDrt, relm(17), PltfSdLnth, URL(3)
      
      WROT=0.D0 
      POWER=0.D0
      MTB=0.D0
      CT=0.D0
      CP=0.D0
      VRELWND(:)=0.D0
   
    ! 推力曲线调用
      
    filename = 'Input\CURVE.txt'
    ! OPEN(47,FILE='Input\CURVE.txt') !CSQ
    ! CALL THRUST(liju, FILENAMe,VREL, mm, nn)  !  mm是文件的行数 nn是文件的列数，也就是说调用几个变量
    CALL PlatformVelocity(Roll, Pitch, Yaw, VSurge, VSway, VHeave, VRoll, VPitch, VYaw, Rcg, Rh, Vlct)
    !CALL RelativeVelocity( YawAng, TiltAng, WingAng, ConeAng, Roll, Pitch, Yaw, VSurge, VSway, VHeave, VRoll, VPitch, VYaw, Xpc, Xyc, Xrc, ZHub, VHub, WindDrt, NBEM, NBL, relm, CT, RTB, PltfSdLnth, URL, IFWAKE, IFCOUP, Ifheading, VRltv )
    
    VRELWND(1)=VWND(1)-Vlct(1)
    VRELWND(2)=VWND(2)-Vlct(2)
    VRELWND(3)=VWND(3)-Vlct(3)
    !VRELWND(:)=VRltv(1,1)
    
    !write(2,*) "The Value of VWND is",VWND(1),VWND(2),VWND(3)
    !write(2,*) "The Value of VRELWND is",VRELWND(1),VRELWND(2),VRELWND(3)
    
    CALL CURVE(THRUST, MTB, filename, VRELWND(1), line, column, WROT, CT, CP)    ! line=215 column=2
    !CT=THRUST/(0.5*1.205*VRELWND*VRELWND*PI*RTB*RTB)
    !CP=POWER/(0.5*1.205*VRELWND*VRELWND*VRELWND*PI*RTB*RTB)
    POWER=CP*0.5*1.205*PI*RTB*RTB*VRELWND(1)**3
    
      WRITE(2,*) 'Calculation results:'  
      WRITE(2,*) 'VWND=',VWND(1),'m/s'
      WRITE(2,*) 'WROT=',WROT,'rpm'
      WRITE(2,*) 'VRELWND=',VRELWND,'m/s'
      WRITE(2,*) 'POWER=',POWER,'Watt'
      WRITE(2,*) 'CT=',CT
      WRITE(2,*) 'CP=',CP
      WRITE(2,*) 'MTB=',MTB,'N*m'
      WRITE(2,'(A10,F16.6,A5)') 'Thrust=',THRUST,'N'
      
      END SUBROUTINE AERO_STEP
    
    subroutine CURVE(THRUST,MTB,filename,VREL,line,column,WROT,CT,CP) 
    Implicit None
     
    character(LEN=20) :: filename
    integer*4 :: line, column
    real*8, dimension(0:line-1,1:6) :: P
    !real*8 :: result
    real*8 :: eps, VREL, THRUST,MTB,WROT,CP,CT
    integer*4 :: i, j,  counter

      eps= 0.05         
    
      open(unit=47,file= filename)
      do i=0,20
     ! do i=0,214     
          read(47,*) (P(i,j),j=1,6)
          !write(*,*) ,(P(i,j),j=1,6)
!          line = i
      enddo
      close(47)
      
      do i=0,19
      !do i=0,212
          !write(*,*) abs(p(i,1)-VREL)
          if(abs(p(i,1)-VREL)<eps) exit
      enddo
      counter = i

      !write(*,*) "The Value of counter is",counter
      
       THRUST = P(counter,2)
       !write(*,*) "The Value of TRUST is",THRUST
       MTB = P(counter,3)
       !write(*,*) "The Value of MTB is",MTB
       CT = P(counter,4)
       !write(*,*) "The Value of CT is",CT
       CP = P(counter,5)
       !write(*,*) "The Value of CP is",CP
       WROT = P(counter,6)
       !write(*,*) "The Value of WROT is",WROT
       
    END SUBROUTINE CURVE  