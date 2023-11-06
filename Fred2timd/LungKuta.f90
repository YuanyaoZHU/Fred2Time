
    MODULE LGKT4K
    USE  WaveDyn_mod
    USE  Platform_mod
    IMPLICIT NONE
    
    
    CONTAINS
    
    
    
    FUNCTION LungKuta4K(t,dt,A0,B0,MASSinv,DAMPING,RESTORING,FORCE_ZHU,L,XCM0,NUMLINES,NSTP,IRR,InitialVerticalForce)    
      IMPLICIT NONE
      
      REAL(8),INTENT(IN)::t,dt
      REAL(8),INTENT(IN)::A0(6,1),B0(6,1) ! 位移 速度
      REAL(8),INTENT(IN)::MASSinv(6,6),DAMPING(6,6),RESTORING(6,6),FORCE_ZHU(6,1)
      INTEGER,INTENT(IN)::L,NUMLINES,NSTP,IRR
      REAL(8),INTENT(IN)::XCM0(6)
      REAL(8),INTENT(IN)::InitialVerticalForce
      REAL(8),DIMENSION(6,2):: LungKuta4K
      
      
      REAL(8),DIMENSION(6,6)::F2inv,F1,F0
      REAL(8),DIMENSION(6,1)::F,Fgrav,Ffloat,Frna,Fwvex
      
      REAL(8),DIMENSION(6,1)::K1A,K1B,K2A,K2B,K3A,K3B,K4A,K4B,Fmoor,Fretar,Fmorison,Fdamping_2kc,Fmorison2
      REAL(8)::t1,t2,t3,t4,MassRna,g
      REAL(8),DIMENSION(6,1)::A01,A02,A03,A04,B01,B02,B03,B04,Y1,Y2
      REAL(8),DIMENSION(6)::A,B
      INTEGER :: n
      
      F2inv = MASSinv !逆质量矩阵
      F1 = DAMPING !线性阻尼
      F0 = RESTORING !回复矩阵
      Fgrav = 0.D0
      
      
      !Fgrav(3,1) = -MATX(1,1)*9.8066
      Ffloat(3,1) = InitialVerticalForce !80708100.D0+1547628.79
      !MassRna = 630291
      g=9.8066
      Frna(5,1)= MassRna*g*xRNA
      
      
      t1 = t
      t2 = t+dt/2
      t3 = t+dt/2
      t4 = t+dt
!------------------------      
      A01 = A0 
      B01 = B0
      DO N=1,6
        A(N) = A01(N,1)
        B(N) = B01(N,1)
      END DO
      
      !WRITE(*,*)'A01 = ',A01(4,1),'    ',A01(5,1),'    ',A01(6,1)
      !Fmoor = Fmoc(NUMLINES,L,XCM0,A01)
      Fmoor = Fmoc(NUMLINES,L,XCM0,A)  !计算锚链力的函数
      WRITE(82,1011)t,Fmoor(1,1),Fmoor(2,1),Fmoor(3,1),Fmoor(4,1),Fmoor(5,1),Fmoor(6,1) !Output\ZHUMOORFORCE.txt
      WRITE(92,1011)t,A01(1,1),A01(2,1),A01(3,1),A01(4,1),A01(5,1),A01(6,1) !Output\ZHUOFFSET.txt
            
      Fretar = Frec(L,1,dt,NSTP) !计算迟滞函数
      WRITE(85,1011)t,Fretar(1,1),Fretar(2,1),Fretar(3,1),Fretar(4,1),Fretar(5,1),Fretar(6,1) !Output\ZHUREATFORCE.txt
      
      Fmorison = Fmoric(L,T,IRR,A0,B01)
      Fmorison2 = Fmoric2(L,T,IRR,A0,B01)  !算fast时加上的计算轴向morison力的，目前没有用上
      WRITE(94,1011)t,Fmorison(1,1),Fmorison(2,1),Fmorison(3,1),Fmorison(4,1),Fmorison(5,1),Fmorison(6,1)
      WRITE(95,1011)t,Fmorison2(1,1),Fmorison2(2,1),Fmorison2(3,1),Fmorison2(4,1),Fmorison2(5,1),Fmorison2(6,1)
      Fwvex = Fwvec(L,1)
      
      Fdamping_2kc = Fdamping_2k(B)
      
      
      
      F = FORCE_ZHU+Fmoor+Fretar+Fgrav+Ffloat+Frna+Fwvex+Fdamping_2kc !目前已经屏蔽风载荷+Fmorison
      
      
      
      K1A = FF(t1,A01,B01)  !速度
      K1B = GG(t1,A01,B01,F2inv,F1,F0,F) !加速度 GG = Massmatrixinv * (F - DAMPING * velocity - Restoring * displacement)
      
      DO N = 1,6
      RK(N,1,L)=K1B(N,1)*dt 
      END DO
      
      
!-------------------------------------------      
      A02 = A0+dt/2*K1A
      B02 = B0+dt/2*K1B
      
      DO N=1,6
        A(N) = A02(N,1)
        B(N) = B02(N,1)
      END DO
      
      !WRITE(*,*)'A02 = ',A02(4,1),'    ',A02(5,1),'    ',A02(6,1)
      
      !Fmoor = Fmoc(NUMLINES,L,XCM0,A)
      WRITE(82,1011)t,Fmoor(1,1),Fmoor(2,1),Fmoor(3,1),Fmoor(4,1),Fmoor(5,1),Fmoor(6,1) !Output\ZHUMOORFORCE.txt
      WRITE(92,1011)t,A02(1,1),A02(2,1),A02(3,1),A02(4,1),A02(5,1),A02(6,1) !Output\ZHUOFFSET.txt
            
      Fretar = Frec(L,2,dt,NSTP)
      WRITE(85,1011)t,Fretar(1,1),Fretar(2,1),Fretar(3,1),Fretar(4,1),Fretar(5,1),Fretar(6,1) !Output\ZHUREATFORCE.txt
      
      Fmorison = Fmoric(L,T,IRR,A0,B02)
      
      Fwvex = Fwvec(L,2)
      
      Fdamping_2kc = Fdamping_2k(B)
      
      F = Force_ZHU+Fmoor+Fretar+Fgrav+Ffloat+Frna+Fwvex+Fdamping_2kc!+Fmorison
      
      K2A = FF(t2,A02,B02)
      K2B = GG(t2,A02,B02,F2inv,F1,F0,F)
      
      DO N = 1,6
      RK(N,2,L)=K2B(N,1)*dt
      END DO
!-------------------------------------------
      A03 = A0+dt/2*K2A
      B03 = B0+dt/2*K2B
      
      !WRITE(*,*)'A03 = ',A03(4,1),'    ',A03(5,1),'    ',A03(6,1)
      
      DO N=1,6
        A(N) = A03(N,1)
        B(N) = B03(N,1)
      END DO      
      
      !Fmoor = Fmoc(NUMLINES,L,XCM0,A)
      WRITE(82,1011)t,Fmoor(1,1),Fmoor(2,1),Fmoor(3,1),Fmoor(4,1),Fmoor(5,1),Fmoor(6,1) !Output\ZHUMOORFORCE.txt
      WRITE(92,1011)t,A03(1,1),A03(2,1),A03(3,1),A03(4,1),A03(5,1),A03(6,1) !Output\ZHUOFFSET.txt
      
      Fretar = Frec(L,3,dt,NSTP)
      WRITE(85,1011)t,Fretar(1,1),Fretar(2,1),Fretar(3,1),Fretar(4,1),Fretar(5,1),Fretar(6,1) !Output\ZHUREATFORCE.txt
      
      Fmorison = Fmoric(L,T,IRR,A0,B03)
      
      Fwvex = Fwvec(L,2)
      
      Fdamping_2kc = Fdamping_2k(B)
      
      !WRITE(*,*)Fdamping_2kc(1,1),'    ',Fdamping_2kc(2,1),'    ',Fdamping_2kc(3,1),'    ',Fdamping_2kc(4,1),'    ',Fdamping_2kc(5,1),'    ',Fdamping_2kc(6,1)
      
      F = Force_ZHU+Fmoor+Fretar+Fgrav+Ffloat+Frna+Fwvex+Fdamping_2kc!+Fmorison
      
      K3A = FF(t3,A03,B03)
      K3B = GG(t3,A03,B03,F2inv,F1,F0,F)
      
      DO N = 1,6
      RK(N,3,L)=K3B(N,1)*dt
      END DO
!-------------------------------------------
      A04 = A0+dt*K3A
      B04 = B0+dt*K3B
      
      !WRITE(*,*)'A04 = ',A04(4,1),'    ',A04(5,1),'    ',A04(6,1)
      
      DO N=1,6
        A(N) = A04(N,1)
        B(N) = B04(N,1)
      END DO
      
      !Fmoor = Fmoc(NUMLINES,L,XCM0,A)
      WRITE(82,1011)t,Fmoor(1,1),Fmoor(2,1),Fmoor(3,1),Fmoor(4,1),Fmoor(5,1),Fmoor(6,1) !Output\ZHUMOORFORCE.txt
      WRITE(92,1011)t,A04(1,1),A04(2,1),A04(3,1),A04(4,1),A04(5,1),A04(6,1) !Output\ZHUOFFSET.txt
      
      Fretar = Frec(L,4,dt,NSTP)
      WRITE(85,1011)t,Fretar(1,1),Fretar(2,1),Fretar(3,1),Fretar(4,1),Fretar(5,1),Fretar(6,1) !Output\ZHUREATFORCE.txt
      
      Fmorison = Fmoric(L,T,IRR,A0,B04)
      
      Fwvex = Fwvec(L,1)
      
      Fdamping_2kc = Fdamping_2k(B)
      
      F = Force_ZHU+Fmoor+Fretar+Fgrav+Ffloat+Frna+Fwvex+Fdamping_2kc!+Fmorison
      
      
      K4A = FF(t4,A04,B04)
      K4B = GG(t4,A04,B04,F2inv,F1,F0,F)
      
      DO N = 1,6
      RK(N,4,L)=K4B(N,1)*dt
      END DO
!--------------------------------------------      
      Y1 = A0+dt/6*(K1A+2*K2A+2*K3A+K4A)
      Y2 = B0+dt/6*(K1B+2*K2B+2*K3B+K4B)
      
      do n=1,6
          Lungkuta4K(n,1)=Y1(n,1)
          Lungkuta4K(n,2)=Y2(n,1)
      end do
      
1011  FORMAT(F20.3,5x,F20.2,5x,F20.2,5x,F20.2,5x,F20.2,5x,F20.2,5x,F0.2)
      
      RETURN
    END
!------------------------------------------------------------
    FUNCTION FF(t,A,B)
      IMPLICIT NONE
      
      REAL(8),INTENT(IN)::t
      REAL(8),INTENT(IN)::A(6,1),B(6,1)
      REAL(8),DIMENSION(6,1):: FF
      
      FF = B
      RETURN
    END FUNCTION FF
    
!-------------------------------------------------------------
    
    FUNCTION GG(t,A,B,F2inv,F1,F0,F)
      IMPLICIT NONE
      
      REAL(8),INTENT(IN)::t
      REAL(8),INTENT(IN)::A(6,1),B(6,1)
      
      REAL(8),INTENT(IN)::F2inv(6,6),F1(6,6),F0(6,6),F(6,1)
      REAL(8),DIMENSION(6,1)::GG
     
      GG = matmul(F2inv,F-(matmul(F1,B)+matmul(F0,A)))
      
      RETURN
    END FUNCTION GG
!---------------------------------------------------------------------------------    
    
    FUNCTION Fmoc(NUMLINES,L,XCM0,A)
      IMPLICIT NONE
      
      INTEGER,INTENT(IN)::NUMLINES,L
      REAL(8),INTENT(IN),DIMENSION(6)::XCM0
      REAL(8),INTENT(IN),DIMENSION(6)::A
      REAL(8),DIMENSION(6,1)::Fmoc
      
      REAL(8),DIMENSION(6)::XCM,FLINE,FILINE
      INTEGER::J
      
      !WRITE(*,*) 'A = ',A(4),'     ',A(5),'    ',A(6)
      DO J=1,6
        XCM(J)=A(J)+XCM0(J)
      END DO
            
      IF (L.EQ.0) THEN           
       CALL CATENARYSYSTEM( .FALSE., NUMLINES, XCM0, XCM, FLINE ,FILINE)
      ELSE
       CALL CATENARYSYSTEM( .TRUE., NUMLINES, XCM0, XCM, FLINE ,FILINE)
      ENDIF
      
      DO J=1,6
        Fmoc(J,1) = FLINE(J)
      END DO
      
      RETURN
      
    END FUNCTION Fmoc
    
!-----------------------------------------------
    FUNCTION Frec(L,TAG,dt,NSTP) !retardation calculation
      IMPLICIT NONE
      
      INTEGER,INTENT(IN)::L,TAG,NSTP
      REAL(8),INTENT(IN)::dt
      REAL(8),DIMENSION(6)::Fretar
      REAL(8),DIMENSION(6,1)::Frec
      
      INTEGER::N
            
      !WRITE(85,"('Iner',5X,F10.2)")RFK(1,1,L) !85, FILE='Output\ZHUREATFORCE.txt'
      
      CALL CONVRETAR(L,TAG,dt,NSTP,RFK,DX,Fretar,RK)
    
      DO N=1,6
        Frec(N,1)=Fretar(N)
      END DO
    
    END FUNCTION Frec
    
!------------------------------------------------------------------    
    FUNCTION Fmoric(L,T,IRR,ROTATION,VE)
      IMPLICIT NONE
      INTEGER,INTENT(IN)::L,IRR
      REAL(8),INTENT(IN)::T,ROTATION(6,1),VE(6,1)
      REAL(8),DIMENSION(6,1)::Fmoric
      
      REAL(8)::ROLL,PITCH,YAW,VELOCITY(6)
      INTEGER::N
      
      ROLL = ROTATION(4,1)
      PITCH = ROTATION(5,1)
      YAW =  ROTATION(6,1)
      
      FDAMP = 0.D0
      
      DO N = 1,6
        VELOCITY(N)= VE(N,1)   
      END DO
      
      CALL ViscousForce(L, T, IRR, Roll, Pitch, Yaw, VELOCITY)
      
      DO N = 1,6
        Fmoric(N,1)= FDAMP(N)
      END DO
      
      RETURN    
    
    END FUNCTION Fmoric
    
    !-------------------------------
    FUNCTION Fwvec(L,TAG)
    IMPLICIT NONE
    
    INTEGER(4), INTENT(IN)::L,TAG
    REAL(8),DIMENSION(6,1)::Fwvec
    INTEGER(4) N
    
    IF (TAG .EQ. 1) THEN
      DO N = 1,6 
        Fwvec(N,1) = Fwave(N,L)
      END DO        
    ELSEIF (TAG .EQ. 2) THEN
      DO N = 1,6
        Fwvec(N,1) = (Fwave(N,L)+Fwave(N,L+1))/2.D0
      END DO        
    ENDIF
    
    END FUNCTION Fwvec
    
    !-----------------------------
    FUNCTION Fdamping_2k(B)
    IMPLICIT NONE
    
    REAL(8),INTENT(IN)::B(6)
    
    REAL(8):: Fdamping_2k(6,1)
    INTEGER(4)::N,FUHAO
    
    DO N=1,6
        IF (B(N).LT.0) THEN
            FUHAO = 1
        ELSE
            FUHAO = -1
        END IF
      Fdamping_2k(N,1) = -CDX(N)*B(N)**2*FUHAO
      !WRITE(*,*) Fdamping_2k(N,1),' CDX=',CDX(N),'  B=',B(N)
    END DO
    
    RETURN
    END FUNCTION Fdamping_2k
    
    FUNCTION Fmoric2(L,T,IRR,ROTATION,VE)
      IMPLICIT NONE
      INTEGER,INTENT(IN)::L,IRR
      REAL(8),INTENT(IN)::T,ROTATION(6,1),VE(6,1)
      REAL(8),DIMENSION(6,1)::Fmoric2
      
      REAL(8)::ROLL,PITCH,YAW,VELOCITY(6)
      INTEGER::N
      
      ROLL = ROTATION(4,1)
      PITCH = ROTATION(5,1)
      YAW =  ROTATION(6,1)
      
      FDAMP = 0.D0
      
      DO N = 1,6
        VELOCITY(N)= VE(N,1)   
      END DO
      
      CALL ViscousForce(L, T, IRR, Roll, Pitch, Yaw, VELOCITY)
      
      DO N = 1,6
        Fmoric2(N,1)= FDAMP2(N)
      END DO
      
      RETURN    
    
    END FUNCTION Fmoric2    
    

      
    END MODULE LGKT4K  