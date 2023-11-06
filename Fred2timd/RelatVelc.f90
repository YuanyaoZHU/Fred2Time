! ====================================================================   
      SUBROUTINE PlatformVelocity( Roll, Pitch, Yaw, VSurge, VSway, VHeave, VRoll, VPitch, VYaw, Rcg, Rh, Vlct )


       IMPLICIT                        NONE

      ! Passed Variables:

       REAL(8), INTENT(IN )      :: Roll    ! The large rotation about X1, (rad).
       REAL(8), INTENT(IN )      :: Pitch    ! The large rotation about X2, (rad).
       REAL(8), INTENT(IN )      :: Yaw    ! The large rotation about X3, (rad).
       REAL(8), INTENT(IN )      :: VSurge    ! The large rotation about X1, (m/s).
       REAL(8), INTENT(IN )      :: VSway    ! The large rotation about X2, (m/s).
       REAL(8), INTENT(IN )      :: VHeave    ! The large rotation about X3, (m/s).         
       REAL(8), INTENT(IN )      :: VRoll    ! The large rotation about X1, (rad).
       REAL(8), INTENT(IN )      :: VPitch    ! The large rotation about X2, (rad).
       REAL(8), INTENT(IN )      :: VYaw    ! The large rotation about X3, (rad).
       REAL(8), INTENT(IN )      :: Rcg 
       REAL(8), INTENT(IN )      :: Rh
       
       REAL(8), INTENT(OUT)      ::  Vlct(3)         ! The resulting transformation matrix from X to x, (-).
       
       ! Local Variables:

       REAL(8)                   :: C4,C5,C6,S4,S5,S6
       
!      ! Compute some intermediate results:
!
       C4 = DCOS(Roll)
       C5 = DCOS(Pitch)
       C6 = DCOS(Yaw)
       
       S4 = DSIN(Roll)
       S5 = DSIN(Pitch)
       S6 = DSIN(Yaw)
!
       ! Define the transformation matrix:
!      
       Vlct(1)= VSurge + Rcg*VPitch*C5 - Rh*VYaw*S6
       Vlct(2)= VSway  - Rcg*VRoll*C4  + Rh*VYaw*C6       
       Vlct(3)= VHeave - Rcg*VRoll*S4  - Rcg*VPitch*S5
       !write(2,*) 'Vlct=' ,Vlct(1),Vlct(2),Vlct(3)

!
      RETURN
      END SUBROUTINE PlatformVelocity
      
! ====================================================================   
      SUBROUTINE WindProfile( Z, ZHub, VHub, Vz )


       IMPLICIT                        NONE

      ! Passed Variables:
      
       REAL(8), INTENT(IN )      :: Z
       REAL(8), INTENT(IN )      :: ZHub
       REAL(8), INTENT(IN )      :: VHub
       
       REAL(8), INTENT(OUT)      :: Vz         
       
       ! Local Variables:

       REAL(8)                   :: aph
       
!      ! Compute some intermediate results:
!
       aph=0.2D0
!
       ! Define the transformation matrix:
!      
       Vz=VHub*(Z/ZHub)**aph

!
      RETURN
      END SUBROUTINE WindProfile
      
! ====================================================================   
      SUBROUTINE CalWakeCenter( WindDrt, x1, y1, x2, y2, x, y )


       IMPLICIT                        NONE

      ! Passed Variables:
      
       REAL(8), INTENT(IN )      :: WindDrt
       REAL(8), INTENT(IN )      :: x1
       REAL(8), INTENT(IN )      :: x2
       REAL(8), INTENT(IN )      :: y1
       REAL(8), INTENT(IN )      :: y2
       
       REAL(8), INTENT(OUT)      :: x
       REAL(8), INTENT(OUT)      :: y 
       
       ! Local Variables:

       REAL(8)                   :: k
       
!      ! Compute some intermediate results:
!
         k=DTAN(WindDrt)
!
       ! Define the transformation matrix:
!      
         x=(k*(y2-y1)+k**2*x1+x2)/(k**2+1.d0)
         y=y1+(k**2*(y2-y1)+k*(x2-x1))/(k**2+1.d0)
!
      RETURN
      END SUBROUTINE CalWakeCenter
      
! ====================================================================   
      SUBROUTINE RelativeVelocity( YawAng, TiltAng, WingAng, ConeAng, Roll, Pitch, Yaw, VSurge, VSway, VHeave, VRoll, VPitch, VYaw, Xpc, Xyc, Xrc, ZHub, VHub, WindDrt, NBEM, NBL, relm, CT, RTB, PltfSdLnth, URL, IFWAKE, IFCOUP, Ifheading, VRltv)

       IMPLICIT                        NONE

      ! Passed Variables:
      
       REAL(8), INTENT(IN )      :: Roll    ! The large rotation about X1, (rad).
       REAL(8), INTENT(IN )      :: Pitch    ! The large rotation about X2, (rad).
       REAL(8), INTENT(IN )      :: Yaw    ! The large rotation about X3, (rad).
       REAL(8), INTENT(IN )      :: YawAng    ! wind turbine yaw angle, (rad).
       REAL(8), INTENT(IN )      :: TiltAng    ! wind turbine tilt angle, (rad).
       REAL(8), INTENT(IN )      :: WingAng    ! Wing angle
       REAL(8), INTENT(IN )      :: ConeAng    ! Cone angle       
       REAL(8), INTENT(IN )      :: VSurge    ! The large rotation about X1, (rad).
       REAL(8), INTENT(IN )      :: VSway    ! The large rotation about X2, (rad).
       REAL(8), INTENT(IN )      :: VHeave    ! The large rotation about X3, (rad).  
       REAL(8), INTENT(IN )      :: VRoll    ! The large rotation about X1, (rad).
       REAL(8), INTENT(IN )      :: VPitch    ! The large rotation about X2, (rad).
       REAL(8), INTENT(IN )      :: VYaw    ! The large rotation about X3, (rad).       
  
       REAL(8), INTENT(IN )      :: Xpc(3)    ! Platfrom mass center in Coordinate System 1 
       REAL(8), INTENT(IN )      :: Xyc(3)    ! Wind turbine yaw center in Coordinate System 2. 
       REAL(8), INTENT(IN )      :: Xrc(3)    ! Wind turbine hub/rotor center in Coordinate System 4. 
       REAL(8), INTENT(IN )      :: URL(3)  ! UpstRotLoc
       REAL(8), INTENT(IN )      :: ZHub
       REAL(8), INTENT(IN )      :: VHub
       REAL(8), INTENT(IN )      :: WindDrt

       REAL(8), INTENT(IN )      :: relm(NBEM)
       REAL(8), INTENT(IN )      :: CT
       REAL(8), INTENT(IN )      :: RTB
       REAL(8), INTENT(IN )      :: PltfSdLnth
       
       INTEGER(4), INTENT(IN )   :: NBEM, NBL
       LOGICAL,    INTENT(IN)    :: IFWAKE, IFCOUP, Ifheading
       
       REAL(8), INTENT(OUT)      :: VRltv(NBEM,NBL)        
       
       ! Local Variables:
       
       INTEGER                   :: IBL, IEL
       REAL(8)                   :: WingAngP, WakeDfc, Vz, Rh, Rcg, Dist, rr, Rxyz, Wc(3), pi, xyz1(3), VPlfm(3), VWind(3), VRt(3), VRt6(3), TransMat(3,3)
!       
       pi=3.141592653589793d0
       
       ! Define the transformation matrix:
!      
       DO IBL = 1, NBL 
           
         WingAngP= WingAng + (IBL-1)*2.d0/3.d0*PI
         
       DO IEL = 1, NBEM

         CALL BladePointLocation( relm(IEL), YawAng, TiltAng, WingAngP, ConeAng, Roll, Pitch, Yaw, Xpc, Xyc, Xrc, xyz1 ) ! Transform the coordinate of any point on the blade from the 6th system to the 1st system
         
         Rh=DSQRT((xyz1(1)-Xpc(1))**2+(xyz1(2)-Xpc(2))**2)
         Rcg=DSQRT(Rh**2+(xyz1(3)-Xpc(3))**2)         
         CALL PlatformVelocity( Roll, Pitch, Yaw, VSurge, VSway, VHeave, VRoll, VPitch, VYaw, Rcg, Rh, VPlfm )
         IF (.NOT.IFCOUP) VPlfm=0.d0
         
         
         CALL WindProfile( xyz1(3), ZHub, VHub, Vz )
         
         CALL CalWakeCenter( WindDrt, URL(1), URL(2), xyz1(1), xyz1(2), Wc(1), Wc(2) )
         Wc(3)=URL(3)
         Dist=DSQRT((URL(1)-Wc(1))**2+(URL(2)-Wc(2))**2+(URL(3)-Wc(3))**2)
         rr=DSQRT((xyz1(1)-Wc(1))**2+(xyz1(2)-Wc(2))**2+(xyz1(3)-Wc(3))**2)
         Rxyz=DSQRT(RTB**2-(xyz1(3)-Wc(3))**2)
         

         
!         CALL ParkWAKE(Ifheading,IFWAKE,rr,Dist,ZHub,VHub,CT,RTB,WakeDfc)
!         CALL B_PWAKE(Ifheading,IFWAKE,rr,Dist,xyz1(2)-Wc(2),xyz1(3)-Wc(3),ZHub,VHub,CT,RTB,WakeDfc)
         CALL LarsenWAKE(Ifheading,IFWAKE,rr,Dist,xyz1(2)-Wc(2),xyz1(3)-Wc(3),ZHub,VHub,CT,RTB,WakeDfc)
         
         VWind(1)= (Vz-WakeDfc)*DCOS(WindDrt)
         VWind(2)= (Vz-WakeDfc)*DSIN(WindDrt)
         VWind(3)= 0.d0
         
         VRt(1)=VWind(1)-VPlfm(1)
         VRt(2)=VWind(2)-VPlfm(2)
         VRt(3)=VWind(3)-VPlfm(3)
         
         CALL Trans1to6( YawAng, TiltAng, WingAngP, ConeAng, Roll, Pitch, Yaw, TransMat )
         
         VRt6(1)=TransMat(1,1)*VRt(1) + TransMat(1,2)*VRt(2) + TransMat(1,3)*VRt(3)
         VRt6(2)=TransMat(2,1)*VRt(1) + TransMat(2,2)*VRt(2) + TransMat(2,3)*VRt(3)
         VRt6(3)=TransMat(3,1)*VRt(1) + TransMat(3,2)*VRt(2) + TransMat(3,3)*VRt(3)
         
         VRltv(IEL,IBL)=VRt6(1)
         !Write(*,*) "The Value of VRltv is",VRltv

       ENDDO
       ENDDO

       !Write(2,*) "The Value of VRltv is",VRltv
       !Write(2,*) "The Value of VWind is",VWind(1),VWind(2),VWind(3)
       !Write(2,*) "The Value of VRt6 is",VRt6(1),VRt6(2),VRt6(3)
      RETURN
      END SUBROUTINE RelativeVelocity
      
! ====================================================================   
      SUBROUTINE BladePointLocation( relm, YawAng, TiltAng, WingAng, ConeAng, Roll, Pitch, Yaw, Xpc, Xyc, Xrc, xyz1 )
! Transform the coordinate of any point on the blade from the 6th system to the 1st system
      IMPLICIT                        NONE

      ! Passed Variables:
      
       REAL(8), INTENT(IN )      :: relm
       REAL(8), INTENT(IN )      :: Roll    ! The large rotation about X1, (rad).
       REAL(8), INTENT(IN )      :: Pitch    ! The large rotation about X2, (rad).
       REAL(8), INTENT(IN )      :: Yaw    ! The large rotation about X3, (rad).     
       REAL(8), INTENT(IN )      :: YawAng    ! wind turbine yaw angle, (rad).
       REAL(8), INTENT(IN )      :: TiltAng    ! wind turbine tilt angle, (rad).
       REAL(8), INTENT(IN )      :: WingAng    ! Wing angle
       REAL(8), INTENT(IN )      :: ConeAng    ! Cone angle       
       REAL(8), INTENT(IN )      :: Xpc(3)    ! Platfrom mass center in Coordinate System 1 
       REAL(8), INTENT(IN )      :: Xyc(3)    ! Wind turbine yaw center in Coordinate System 2. 
       REAL(8), INTENT(IN )      :: Xrc(3)    ! Wind turbine hub/rotor center in Coordinate System 4. 
       
       REAL(8), INTENT(OUT)      :: xyz1(3)
       
       ! Local Variables:

       REAL(8)                   :: xyz2(3),xyz4(3),xyz5(3),xyz6(3),TransMat(3,3)
       
!      ! Compute some intermediate results:

       xyz6(1)    = 0.d0
       xyz6(2)    = 0.d0
       xyz6(3)    = relm

       CALL Trans5to6( ConeAng, TransMat )
       
       xyz5(1) = TransMat(1,1)*xyz6(1) + TransMat(2,1)*xyz6(2) + TransMat(3,1)*xyz6(3)
       xyz5(2) = TransMat(1,2)*xyz6(1) + TransMat(2,2)*xyz6(2) + TransMat(3,2)*xyz6(3)
       xyz5(3) = TransMat(1,3)*xyz6(1) + TransMat(2,3)*xyz6(2) + TransMat(3,3)*xyz6(3)
       
       CALL Trans4to5( WingAng, TransMat )
       
       xyz4(1) = TransMat(1,1)*xyz5(1) + TransMat(2,1)*xyz5(2) + TransMat(3,1)*xyz5(3)
       xyz4(2) = TransMat(1,2)*xyz5(1) + TransMat(2,2)*xyz5(2) + TransMat(3,2)*xyz5(3)
       xyz4(3) = TransMat(1,3)*xyz5(1) + TransMat(2,3)*xyz5(2) + TransMat(3,3)*xyz5(3)       

       xyz4(1)    = Xrc(1) + xyz4(1)
       xyz4(2)    = Xrc(2) + xyz4(2)
       xyz4(3)    = Xrc(3) + xyz4(3)
         
       CALL Trans3to4( YawAng, TiltAng, TransMat )
       
       xyz2(1) = TransMat(1,1)*xyz4(1) + TransMat(2,1)*xyz4(2) + TransMat(3,1)*xyz4(3)
       xyz2(2) = TransMat(1,2)*xyz4(1) + TransMat(2,2)*xyz4(2) + TransMat(3,2)*xyz4(3)
       xyz2(3) = TransMat(1,3)*xyz4(1) + TransMat(2,3)*xyz4(2) + TransMat(3,3)*xyz4(3)       

       xyz2(1)    = Xyc(1) + xyz2(1)
       xyz2(2)    = Xyc(2) + xyz2(2)
       xyz2(3)    = Xyc(3) + xyz2(3)
       
       CALL Trans1to2( Roll, Pitch, Yaw, TransMat )
       
       xyz1(1) = TransMat(1,1)*xyz2(1) + TransMat(2,1)*xyz2(2) + TransMat(3,1)*xyz2(3)
       xyz1(2) = TransMat(1,2)*xyz2(1) + TransMat(2,2)*xyz2(2) + TransMat(3,2)*xyz2(3)
       xyz1(3) = TransMat(1,3)*xyz2(1) + TransMat(2,3)*xyz2(2) + TransMat(3,3)*xyz2(3)       

       xyz1(1)    = Xpc(1) + xyz1(1)
       xyz1(2)    = Xpc(2) + xyz1(2)
       xyz1(3)    = Xpc(3) + xyz1(3)
!
       
      RETURN
      END SUBROUTINE BladePointLocation
      
! ====================================================================   
      SUBROUTINE RotorCenterLocation( YawAng, TiltAng,  Roll, Pitch, Yaw, Xpc, Xyc, Xrc, xyz1 )
! Transform the coordinate of rotor hub center from the 4th system to the 1st system
      IMPLICIT                        NONE

      ! Passed Variables:
      
       REAL(8), INTENT(IN )      :: Roll    ! The large rotation about X1, (rad).
       REAL(8), INTENT(IN )      :: Pitch    ! The large rotation about X2, (rad).
       REAL(8), INTENT(IN )      :: Yaw    ! The large rotation about X3, (rad).     
       REAL(8), INTENT(IN )      :: YawAng    ! wind turbine yaw angle, (rad).
       REAL(8), INTENT(IN )      :: TiltAng    ! wind turbine tilt angle, (rad).
     
       REAL(8), INTENT(IN )      :: Xpc(3)    ! Platfrom mass center in Coordinate System 1 
       REAL(8), INTENT(IN )      :: Xyc(3)    ! Wind turbine yaw center in Coordinate System 2. 
       REAL(8), INTENT(IN )      :: Xrc(3)    ! Wind turbine hub/rotor center in Coordinate System 4. 
       
       REAL(8), INTENT(OUT)      :: xyz1(3)
       
       ! Local Variables:

       REAL(8)                   :: xyz2(3),xyz4(3),TransMat(3,3)
       
!      ! Compute some intermediate results:

       xyz4(1)    = Xrc(1)
       xyz4(2)    = Xrc(2)
       xyz4(3)    = Xrc(3)
         
       CALL Trans3to4( YawAng, TiltAng, TransMat )
       
       xyz2(1) = TransMat(1,1)*xyz4(1) + TransMat(2,1)*xyz4(2) + TransMat(3,1)*xyz4(3)
       xyz2(2) = TransMat(1,2)*xyz4(1) + TransMat(2,2)*xyz4(2) + TransMat(3,2)*xyz4(3)
       xyz2(3) = TransMat(1,3)*xyz4(1) + TransMat(2,3)*xyz4(2) + TransMat(3,3)*xyz4(3)       

       xyz2(1)    = Xyc(1) + xyz2(1)
       xyz2(2)    = Xyc(2) + xyz2(2)
       xyz2(3)    = Xyc(3) + xyz2(3)
       
       CALL Trans1to2( Roll, Pitch, Yaw, TransMat )
       
       xyz1(1) = TransMat(1,1)*xyz2(1) + TransMat(2,1)*xyz2(2) + TransMat(3,1)*xyz2(3)
       xyz1(2) = TransMat(1,2)*xyz2(1) + TransMat(2,2)*xyz2(2) + TransMat(3,2)*xyz2(3)
       xyz1(3) = TransMat(1,3)*xyz2(1) + TransMat(2,3)*xyz2(2) + TransMat(3,3)*xyz2(3)       

       xyz1(1)    = Xpc(1) + xyz1(1)
       xyz1(2)    = Xpc(2) + xyz1(2)
       xyz1(3)    = Xpc(3) + xyz1(3)      
!
       
      RETURN
      END SUBROUTINE RotorCenterLocation
      
! ====================================================================   
      SUBROUTINE CalForceTorque( WindDrt, TrstForc, AeroTorq, Xrl, Xpc, Force )


       IMPLICIT                        NONE

      ! Passed Variables:
      
       REAL(8), INTENT(IN )      :: WindDrt
       REAL(8), INTENT(IN )      :: TrstForc
       REAL(8), INTENT(IN )      :: AeroTorq
       REAL(8), INTENT(IN )      :: Xrl(3)
       REAL(8), INTENT(IN )      :: Xpc(3)       
       REAL(8), INTENT(OUT)      :: Force(6)
       
       ! Local Variables:

       REAL(8)                   :: r(3),tForce(3)
       
!      ! Compute some intermediate results:
 
         r=Xrl-Xpc
         
         Force(1)=TrstForc*DCOS(WindDrt)
         Force(2)=TrstForc*DSIN(WindDrt)
         Force(3)=0.D0
         
         tForce(1)=AeroTorq*DCOS(WindDrt)
         tForce(2)=AeroTorq*DSIN(WindDrt)
         tForce(3)=0.D0
         
       ! Define the transformation matrix:
!      
         Force(4)=r(2)*Force(3)-r(3)*Force(2)+ tForce(1)
         Force(5)=r(3)*Force(1)-r(1)*Force(3)+ tForce(2)
         Force(6)=r(1)*Force(2)-r(2)*Force(1)+ tForce(3)
         
!
      RETURN
      END SUBROUTINE CalForceTorque
      