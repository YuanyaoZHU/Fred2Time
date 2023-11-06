 
   !-------------------------------------------------------------------------------------------------
   ! Aerodynamic calculation for the rotor at every time step...
   !-------------------------------------------------------------------------------------------------            
        
      SUBROUTINE AERO_STEP(VWND,WROT,PTAG,THRUST,MTB,POWER,CT,CP)
      
      USE Const_mod     
      USE Blade_mod  
      IMPLICIT NONE   
      
      !INTEGER IEL,IBL,J
      INTEGER*4 :: line=215, column=2
      REAL*8 :: THRUST !CSQ
      CHARACTER(LEN=20) :: filename !CSQ filename的字符长度定义
      REAL*8  VWND(3),WROT,PTAG,FTB,MTB,CT,CP,power
      REAL*8  VRELWND
      REAL*8  Vlct(3),Roll, Pitch, Yaw, VSurge, VSway, VHeave, VRoll, VPitch, VYaw, Rcg, Rh
      
      WROT=0.D0
      POWER=0.D0
      MTB=0.D0
      CT=0.D0
      CP=0.D0
      VRELWND=VWND(1) 
   
    ! 推力曲线调用
      
    filename = 'Input\CURVE.txt'
    ! OPEN(47,FILE='Input\CURVE.txt') !CSQ
    ! call THRUST(liju, FILENAMe,VREL, mm, nn)  !  mm是文件的行数 nn是文件的列数，也就是说调用几个变量
    CALL PlatformVelocity(Roll, Pitch, Yaw, VSurge, VSway, VHeave, VRoll, VPitch, VYaw, Rcg, Rh, Vlct)
    
    VRELWND=VRELWND-Vlct(1)
    
    CALL CURVE(THRUST,MTB,filename,VRELWND,line,column)    ! line=215 column=2
    CT=THRUST/(0.5*1.205*VRELWND*VRELWND*PI*RTB*RTB)
    CP=POWER/(0.5*1.205*VRELWND*VRELWND*VRELWND*PI*RTB*RTB)
    POWER=WROT*MTB
    
      WRITE(2,*) 'Calculation results:'  
      WRITE(2,*) 'VWND=',VWND(1),'WROT=',WROT
      !,'VRELWND=',VRELWND(1,1)     
      
      WRITE(2,*)  
     ! WRITE(2,*) '   CT=',CT
      WRITE(2,'(A10,F16.6,A5)') 'Thrust=',THRUST,'N'    
    
      END SUBROUTINE AERO_STEP
    
    subroutine CURVE(THRUST,MTB,filename,VREL,line,column) 
    Implicit None
     
    character(LEN=20) :: filename
    integer*4 :: line, column
    real*8, dimension(0:line-1,1:2) :: P
    !real*8 :: result
    real*8 :: eps, VREL, THRUST
    integer*4 :: i, j,  counter

      eps= 0.05         
      
      open(unit=47,file= filename)
      do i=0,214
          read(47,*) (P(i,j),j=1,2)
          !write(*,*) i,(P(i,j),j=1,2)
!          line = i
      enddo
      close(47)
      
      do i=0,212
          !write(*,*) abs(p(i,1)-VREL)
          if(abs(p(i,1)-VREL)<eps) exit
      enddo
      counter = i
       THRUST = P(counter,2)
       !write(*,*) THRUST
    END SUBROUTINE CURVE  