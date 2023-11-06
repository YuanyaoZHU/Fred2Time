      
   !-------------------------------------------------------------------------------------------------
   ! Variable-speed generator-torque controller and PI collective blade-pitch controller...
   !-------------------------------------------------------------------------------------------------   
   
      SUBROUTINE CONTROL(WROT,ATORQ,PTAG,TSP,NWROT,DWROT,NPTAG)
      
      USE Const_mod      
      USE Blade_mod   
      IMPLICIT NONE
      
      INTEGER I,J         
      REAL*8  KOPT,COPT,TSROPT,NGEAR,IDRV,EFM
      REAL*8  WROT,DWROT,NWROT,TSP,TGEN,ATORQ,WGEN,PTAG,NPTAG
 
   !-------------------------------------------------------------------------------------------------
   ! Initialization of control sytem parameters...
   !-------------------------------------------------------------------------------------------------        
      NGEAR=97.D0
      EFM=0.95D0
      IDRV=5025500/SHRO**4
      COPT=0.435994D0
      TSROPT=9.02795D0
      KOPT=PI*RHOA*RTB**5*COPT/(2.D0*NGEAR**3*TSROPT**3*EFM)
 
   !-------------------------------------------------------------------------------------------------
   ! Optimum mode torque gain controller...
   !-------------------------------------------------------------------------------------------------        
   
      WGEN=NGEAR*WROT    
      TGEN=KOPT*WGEN**2      
      DWROT=(ATORQ-NGEAR*TGEN)/IDRV
      
      NWROT=WROT+TSP*DWROT
      NPTAG=PTAG
      
    
      WRITE(2,*)   
      WRITE(2,*) 'AFTER CONTROL' 
      WRITE(2,*) 'NWROT=',NWROT       
      WRITE(2,*) 'DWROT=',DWROT      
      WRITE(2,*) 'ATORQ=',ATORQ,'NGEAR*TGEN=',NGEAR*TGEN          
   !-------------------------------------------------------------------------------------------------
   ! Variable-speed generator-torque controller...
   !-------------------------------------------------------------------------------------------------       
          
      
      END SUBROUTINE CONTROL      
      
      
      
   !-------------------------------------------------------------------------------------------------
   ! Multiple-wake modelling using Park's model...
   !-------------------------------------------------------------------------------------------------   
   
      SUBROUTINE ParkWAKE(Ifheading,IFWAKE,rr,X01,ZHub,U0,CT,R0,DV01)
           
      IMPLICIT NONE   
      
      LOGICAL,INTENT(IN)::IFWAKE,Ifheading
      REAL*8  U0,CT,R0,ZHub,X01,DV01,rr
      REAL*8  KCST,RX
   
      KCST=0.5D0/DLOG(ZHub/0.0002)
      RX=R0+KCST*X01
      
      IF (IFWAKE.AND..NOT.Ifheading.AND.rr.LE.RX) THEN

!      DV01=U0*(1.D0-DSQRT(DABS(1.D0-CT)))/(1.D0+KCST*X01/R0)**2
       IF (CT.GT.1.D0) PRINT*,'CT.GT.1.D0'
       DV01=U0*(1.D0-DSQRT(1.D0-CT))/(1.D0+KCST*X01/R0)**2
       
       Write(80,100) RX,R0,X01,CT,KCST,ZHub,DV01

      ELSE
       DV01=0.D0
      ENDIF

100   Format(10E14.6)          
      END SUBROUTINE ParkWAKE
      
      
      
   !-------------------------------------------------------------------------------------------------
   ! Multiple-wake modelling using B_P's model...
   !-------------------------------------------------------------------------------------------------   
   
      SUBROUTINE B_PWAKE(Ifheading,IFWAKE,rr,X01,Y01,Z01,ZHub,U0,CT,R0,DV01)

      IMPLICIT NONE   
      
      LOGICAL,INTENT(IN)::IFWAKE,Ifheading
      REAL*8  U0,CT,R0,X01,Y01,Z01,ZHub,DV01,RR
      REAL*8  KCST,KGRT,D,B,E,H1,H2,RX,RX1
      
       D=2.D0*R0
       KGRT=0.0005D0*D
       KCST=0.5D0/DLOG(ZHub/0.000002)
       B=(1.D0+DSQRT(1.D0-CT))/2.D0/DSQRT(1.D0-CT)
       E=0.25D0*DSQRT(B)
       RX=R0*DSQRT(B+KCST*X01/D)
       RX1=R0+0.04D0*X01

      IF (IFWAKE.AND..NOT.Ifheading.AND.rr.LE.RX1) THEN

       H1=2.D0*(KGRT*X01/D+E)**2
       H2=4.D0*H1
      
       DV01=U0*(1.D0-DSQRT(1.D0-CT/H2))*DEXP( -1.D0/H1*( ((Z01-ZHub)/D)**2+(Y01/D)**2 ) )
         write(80,100)  RX,R0,X01,CT,KCST,ZHub,DV01

      ELSE
       DV01=0.D0
      ENDIF 
100   Format(10E14.6)       
      END SUBROUTINE B_PWAKE
      
      
   !-------------------------------------------------------------------------------------------------
   ! Multiple-wake modelling using Larsen's model...
   !-------------------------------------------------------------------------------------------------   
   
      SUBROUTINE LarsenWAKE(Ifheading,IFWAKE,rr,X01,Y01,Z01,ZHub,U0,CT,R0,DV01)

      IMPLICIT NONE
      
      LOGICAL,INTENT(IN)::IFWAKE,Ifheading
      REAL*8  U0,CT,R0,X01,Y01,Z01,ZHub,DV01,RR
      REAL*8  IA,D,A,RNB,R95,DEFF,X0,C1,RX,PI,RX1
      
      PI=3.141592653589793D0
      
      
       IA=0.06D0
       D=2.D0*R0
       A=PI*R0**2
       RNB=MAX(1.08D0*D,1.08D0*D+21.7D0*D*(IA-0.05D0))
       R95=0.5D0*(RNB+MIN(ZHub,RNB))
       
       DEFF=D*DSQRT((1.D0+DSQRT(1.D0-CT))/(2.D0*DSQRT(1.D0-CT)))
       X0=9.5D0*D/((2.D0*R95/DEFF)**3-1.D0)
       C1=(DEFF/2.0)**2.5D0/DSQRT(105.D0/(2.D0*PI))/(CT*A*X0)**(5.D0/6.D0)
       
       RX=(35.D0/(2.D0*PI))**0.2D0*(3.D0*C1**2)**0.2D0*(CT*A*(X01+X0))**(1.D0/3.D0)
       RX1=R0+0.04D0*X01
       

         
      IF (IFWAKE.AND..NOT.Ifheading.AND.rr.LE.RX1) THEN

       DV01=U0/9.D0*(CT*A/(X01+X0)**2)**(1.D0/3.D0)*( RR**1.5D0/DSQRT(3.D0*C1**2*CT*A*(X01+X0))-(35.D0/(2.D0*PI))**0.3D0/(3.D0*C1**2)**0.2D0 )**2
       
!         write(80,100) DEFF,CT,A,X0,DV01
  
      ELSE
       DV01=0.D0
      ENDIF 
      
100   Format(10E14.6)      
      END SUBROUTINE LarsenWAKE