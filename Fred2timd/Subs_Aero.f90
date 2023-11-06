   !-------------------------------------------------------------------------------------------------
   ! Initialization for the use of unseady bem routine...
   !-------------------------------------------------------------------------------------------------        
   
      SUBROUTINE INITIAL(VWND,WROT,PTAG)
      
      USE Const_mod     
      USE Blade_mod
      
      IMPLICIT NONE    
      
      INTEGER IEL,I,J,ERR
      REAL*8  VWND(3),WROT(3),PTAG(3)
      
   !-------------------------------------------------------------------------------------------------
   ! Read blade element parameters and the airfoil shape data.
   !-------------------------------------------------------------------------------------------------      
   
       !SHRO= 63.D0/11.692D0
      SHRO= 1.D0
       
!       IFDUCT=.FALSE.
       
        READ(1,*)  
        READ(1,*) 
        I=1
       DO WHILE(.NOT.EOF(1))
	  READ(1,*,IOSTAT=ERR) J,RND(I),TWST(I),DRND(I),CHRD(I),NFOIL(I)
	  IF(ERR/=0) EXIT 
	  I=I+1
       ENDDO
       
       DO J=1,I-1
        RND(J)=RND(J)/SHRO
        DRND(J)=DRND(J)/SHRO
        CHRD(J)=CHRD(J)/SHRO 
       ENDDO        
       
       DO I=11,18       
       REWIND(I)
       ENDDO     
       
       READ(11,*)
       DO I=1,NCYLD
       READ(11,'(F8.2,2X,3(2X,F7.4))') (CYLD1(I,J),J=1,4)
       ENDDO       
       
       READ(12,*)
       DO I=1,NCYLD
       READ(12,'(F8.2,2X,3(2X,F7.4))') (CYLD2(I,J),J=1,4)
       ENDDO      
       
       READ(13,*)
       DO I=1,NDU21
       READ(13,'(F8.2,2X,3(2X,F7.4))') (DU21(I,J),J=1,4)
       ENDDO
       
       READ(14,*)
       DO I=1,NDU25
       READ(14,'(F8.2,2X,3(2X,F7.4))') (DU25(I,J),J=1,4)
       ENDDO 
       
       READ(15,*)
       DO I=1,NDU30
       READ(15,'(F8.2,2X,3(2X,F7.4))') (DU30(I,J),J=1,4)
       ENDDO
       
       READ(16,*)
       DO I=1,NDU35
       READ(16,'(F8.2,2X,3(2X,F7.4))') (DU35(I,J),J=1,4)
       ENDDO        
       
       READ(17,*)
       DO I=1,NDU40
       READ(17,'(F8.2,2X,3(2X,F7.4))') (DU40(I,J),J=1,4)
       ENDDO
       
       READ(18,*)
       DO I=1,NNACA64
       READ(18,'(F8.2,2X,3(2X,F7.4))') (NACA64(I,J),J=1,4)
       ENDDO                

       DO I=11,18   
       CLOSE(I)
       ENDDO  
      
   !-------------------------------------------------------------------------------------------------
   ! Initialization of the parameters...
   !-------------------------------------------------------------------------------------------------         

       PTAG(:)= 0.D0       
       VWND(:)=11.4D0    ! Wind speed (m/s)
       !WROT(:)=(15.6D0)/60.D0*SHRO
                         ! Rotation speed of the wind turbine (rad/s)
       RTB=63.D0/SHRO    ! Radius of wind turbine (m)
       RHB=1.5D0/SHRO    ! Radius of hub (m)        
       
      DO IEL=1,NBEM 
        
        XR(IEL)=RND(IEL)/RTB 
     
        SDT(IEL)=NBL*CHRD(IEL)/(2.D0*PI*RND(IEL)) ! Solidity
        
        CALL ADUCT(XR(IEL)/2.D0,A0(IEL))   

 
      ENDDO        
      
      END SUBROUTINE INITIAL
      
   !-------------------------------------------------------------------------------------------------
   ! Function integrand of Cp integral...
   !-------------------------------------------------------------------------------------------------        
   
      REAL*8 FUNCTION CPX(AXF,TGF,TSR,XR,FTL,LTD,IFDUCT,OPTION,A0) 
      
      IMPLICIT NONE
      INTEGER OPTION
      REAL*8  AXF,TGF,TSR,XR,FTL,LTD,A0,LXF
      LOGICAL IFDUCT
!      
!  Note: Option=1, use M.O.L.Hansen's Cp formula
!        Option=2, use P.Jamieson's Cp formula
!        Option=3, use Momentum theory Cp formula
!   
      IF (IFDUCT) THEN

        LXF=(AXF-A0)/(1.D0-A0)
      
      ELSE
           
        LXF=AXF
                  
      ENDIF
       
      IF (OPTION.EQ.1) THEN           

        CPX=FTL*TGF*(1.D0-LXF)*XR**3  
        CPX=CPX*8.D0*TSR**2  
                    
      ELSEIF (OPTION.EQ.2) THEN          
      
        CPX=8.D0*LXF*(1.D0-LXF)*FTL*(LTD*(1.D0-LXF)-TSR*XR*(1.D0+TGF))*TSR*XR**2/(LTD*TSR*XR*(1.D0+TGF)+(1.D0-LXF)) 
        
      ELSE
          
        CPX=4.D0*LXF*(1.D0-LXF)**2*FTL
      
      ENDIF
      
      CPX=CPX/3.D0
      
      RETURN
      END FUNCTION CPX
 
   !-------------------------------------------------------------------------------------------------
   ! Compute a0 by polyfit ...
   !-------------------------------------------------------------------------------------------------        
   
      SUBROUTINE ADUCT(X,A0) 
      
      IMPLICIT NONE
      REAL*8 X,UX0,A0
      
      UX0=5.387188213350107*X**3-0.316704873345186*X**2-0.461512203246055*X+1.341198842000954
    
      A0=1.D0-UX0 !-0.3D0
            
      RETURN
      END SUBROUTINE ADUCT           
      
      
   !-------------------------------------------------------------------------------------------------
   ! Calculate tip-loss, hub-loss and total-loss factors...
   !-------------------------------------------------------------------------------------------------        
   
      SUBROUTINE TIPLS(R,BETA,OPTION,RTB,RHB,FTP,FHB,FTL)        
      
      IMPLICIT NONE
      INTEGER OPTION
      REAL*8 R,BETA,FTP,FHB,FTL,F1,F2,RTB,RHB,PI
!      
!  Note: Option=0, no N-S equation correction
!        Otherwise, use N-S equation correction (Xu and Sankar,2002)
!
      PI=4.D0*DATAN(1.D0)
      
      IF (OPTION.EQ.0) THEN 
          
       F1=0.5D0*(RTB-R)/R/DSIN(BETA)
       F2=0.5D0*(R-RHB)/RHB/DSIN(BETA)
       
       FTP=2.D0/PI*DACOS(DEXP(-F1))
       FHB=2.D0/PI*DACOS(DEXP(-F2))
       
       FTL=FTP*FHB  
 
      ELSE 
      
       IF (R/RTB.GE.0.7D0.AND.R/RTB.LE.1.D0) THEN
           
        F1=0.5D0*(RTB-R)/R/DSIN(BETA)
        F2=0.5D0*(R-RHB)/RHB/DSIN(BETA)
       
        FTP=2.D0/PI*DACOS(DEXP(-F1))
        FTP=(FTP**0.85D0+0.5D0)/2.D0          
        
        FHB=2.D0/PI*DACOS(DEXP(-F2))  
       
        FTL=FTP*FHB 
        
       ELSEIF (R/RTB.LT.0.7D0) THEN  
  
        F1=1.5D0/7.D0/DSIN(BETA)
        F2=0.5D0*(R-RHB)/RHB/DSIN(BETA)  
        
        FTP=2.D0/PI*DACOS(DEXP(-F1))
        FTP=1.D0-R/RTB*(1.D0-FTP)/0.7D0          
        
        FHB=2.D0/PI*DACOS(DEXP(-F2))  
       
        FTL=FTP*FHB
        
       ELSE
           
        PRINT*, ' Error in tip-loss calculation.'          
        PRINT*, ' Local radius larger than turbine radius.'
        
       ENDIF       
      ENDIF
      
      
      IF (OPTION.NE.0.AND.R.EQ.RTB.AND.FTP.GT.0.D0) THEN  
!     Ref. to: (UAE Phase 6, Hand et al. 2001)          
        PRINT*, ' Error in tip-loss calculation.'          
        PRINT*, ' NS correction leading to mistake at tip.'   
        PRINT*, ' Tip-loss factor larger than zero.'
        
      ENDIF                     
      
      RETURN      
      END SUBROUTINE TIPLS
      
   
   !-------------------------------------------------------------------------------------------------
   ! Glauert correction and modified Glauert correction...
   !-------------------------------------------------------------------------------------------------   
   
      SUBROUTINE GLAUERT(CT,FTL,BETA,SDT,XR,LTD,A0,C1,C2,OPTION,IFDUCT,AXF,TGF)         
        
      IMPLICIT NONE
      INTEGER OPTION
      LOGICAL IFDUCT
      REAL*8  CT,FTL,BETA,SDT,XR,LTD,C1,C2,AXF,TGF,A0
      REAL*8  AC,K,P1,P2
!      
!  Note: Option=0, normal Glauert correction  (Spera, 1994)
!        Otherwise, modified Glauert correction  (Buhl, 2004)    

      AC=0.2D0
      TGF=1.D0/(-1.D0+4.D0*FTL*DSIN(BETA)*DCOS(BETA)/SDT/C2)
      
      IF (IFDUCT) THEN
          
        K=4.D0*FTL*DSIN(BETA)**2/SDT/C1
        AXF=(1.D0+A0*K/(1.D0-A0)**2)/(1.D0+K/(1.D0-A0)**2)    
        
        IF (AXF.GT.A0+AC*(1.D0-A0)) THEN
            
         P1=2.D0*(1.D0-A0)
         P2=K*(1.D0-2.D0*AC)
         
         AXF=(P1+P2-DSQRT((P1+P2)**2-P1*(P1+2.D0*P2*A0-K*AC**2*P1)))/P1
        ELSE            
           RETURN 
        ENDIF  
        
      ELSE
          
        K=4.D0*FTL*DSIN(BETA)**2/SDT/C1
        AXF=1.D0/(1.D0+K)            
        
      IF (OPTION.EQ.0) THEN
      
        IF (AXF.GT.AC) THEN
         AXF=(2.D0+K*(1.D0-2.D0*AC)-DSQRT((K*(1.D0-2.D0*AC)+2.D0)**2+4.D0*(K*AC**2-1.D0)))/2.D0
        ELSE            
           RETURN 
        ENDIF   
      
      ELSE
          
        IF (CT.GT.0.96D0*FTL) THEN
            
        AXF=(18.D0*FTL-20.D0-3.D0*DSQRT(CT*(50.D0-36.D0*FTL)+12.D0*FTL*(3.D0*FTL-4.D0)))/(36.D0*FTL-50.D0) 
        
        ELSE
           RETURN
        ENDIF 
        
      ENDIF        
          
      ENDIF
        

100     RETURN
        END SUBROUTINE GLAUERT    
            
      
   !-------------------------------------------------------------------------------------------------
   ! Find the lift and drag coefficients...
   !-------------------------------------------------------------------------------------------------    
      
      SUBROUTINE COFDL(INC,NF,CCL,CCD,INC0L) 
      
      USE Blade_mod      
      IMPLICIT NONE
      
      INTEGER NUM,ID,I,NF
      REAL*8 INC,CCL,CCD,INC0L 
      REAL*8,ALLOCATABLE:: X(:),Y1(:),Y2(:)

      IF (NF.EQ.1) THEN
       NUM=NCYLD
       ALLOCATE(X(NUM),Y1(NUM),Y2(NUM))
       X=CYLD1(:,1)
       Y1=CYLD1(:,2)
       Y2=CYLD1(:,3)
      ELSEIF (NF.EQ.2) THEN
       NUM=NCYLD
       ALLOCATE(X(NUM),Y1(NUM),Y2(NUM))
       X=CYLD2(:,1)
       Y1=CYLD2(:,2)
       Y2=CYLD2(:,3)
      ELSEIF (NF.EQ.3) THEN
       NUM=NDU21 
       ALLOCATE(X(NUM),Y1(NUM),Y2(NUM))
       X=DU21(:,1)
       Y1=DU21(:,2)
       Y2=DU21(:,3)
      ELSEIF (NF.EQ.4) THEN
       NUM=NDU25 
       ALLOCATE(X(NUM),Y1(NUM),Y2(NUM))
       X=DU25(:,1)
       Y1=DU25(:,2)
       Y2=DU25(:,3)
      ELSEIF (NF.EQ.5) THEN
       NUM=NDU30
       ALLOCATE(X(NUM),Y1(NUM),Y2(NUM))
       X=DU30(:,1)
       Y1=DU30(:,2)
       Y2=DU30(:,3)
      ELSEIF (NF.EQ.6) THEN
       NUM=NDU35 
       ALLOCATE(X(NUM),Y1(NUM),Y2(NUM))
       X=DU35(:,1)
       Y1=DU35(:,2)
       Y2=DU35(:,3)      
      ELSEIF (NF.EQ.7) THEN
       NUM=NDU40
       ALLOCATE(X(NUM),Y1(NUM),Y2(NUM))
       X=DU40(:,1)
       Y1=DU40(:,2)
       Y2=DU40(:,3)       
      ELSEIF (NF.EQ.8) THEN
       NUM=NNACA64
       ALLOCATE(X(NUM),Y1(NUM),Y2(NUM))
       X=NACA64(:,1)
       Y1=NACA64(:,2)
       Y2=NACA64(:,3)       
      ELSE
       PRINT*,'Error. No such file found'
      ENDIF
  
      ID=0
      DO I=1,NUM-1
      IF(inc>=X(i).AND.inc<X(I+1)) THEN
      ID=1
      EXIT
      END IF
      END DO

      IF(ID==0) THEN
      PRINT*, ' Error. Incident angle out of range'
      STOP
      ELSE
      ccl=( (y1(i+1)-y1(i))/(x(i+1)-x(i)) )* (inc-x(i)) + y1(i)
      ccd=( (y2(i+1)-y2(i))/(x(i+1)-x(i)) )* (inc-x(i)) + y2(i)      
      ENDIF
      
      IF (NF.GE.3) THEN
       ID=0
       DO I=1,NUM-1
       IF(0.D0>=Y1(i).AND.0.D0<Y1(I+1).AND.X(I)>=-10.D0.AND.X(I+1)<5.D0) THEN
       ID=1
       EXIT
       END IF
       END DO

       IF(ID==0) THEN
       PRINT*, ' Error. Cannot find the zero lift force angle.'
       STOP
       ELSE
       inc0l=-( (x(i+1)-x(i))/(y1(i+1)-y1(i)) )* y1(i) + x(i) 
       ENDIF
      
      ELSE
       inc0l=0.d0      
      ENDIF
                
      RETURN
      END SUBROUTINE COFDL      
      
      
 