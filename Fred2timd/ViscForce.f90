!-------------------------------------------------------------------------------   
!                  Computing wave drag force.
!-------------------------------------------------------------------------------  
! 
      SUBROUTINE ViscousForce(L, T, IRR, Roll, Pitch, Yaw, VCG)
 
	   USE Const_mod
	   USE STKMesh_mod
	   USE WaveDyn_mod
	   USE Platform_mod
       
       IMPLICIT NONE  

	   INTEGER IEL,I,J,K,L,IRR
	   REAL*8 NMAT(3,3),NMAT2(3,3),DMAT(3,6),NDMAT(3,6),RMAT(6,3),XMP(3),AVRT,T,ROLL,PITCH,YAW
       REAL*8 VCG(6),VF(3),VFT(3),VMT(3),VRT(3),PVRT(3),MVCT(6),DINC(3),TransMat(3,3)
       
       FDAMP= 0.D0
!       VCG= DX(:,L)

       DO 10 IEL=1,NSTRE
           
!         CALL Trans1to2( Roll, Pitch, Yaw, TransMat )
       
!         XMP(1) = TransMat(1,1)*RELM(IEL,1) + TransMat(2,1)*RELM(IEL,2) + TransMat(3,1)*RELM(IEL,3) + XCM(1,L)
!         XMP(2) = TransMat(1,2)*RELM(IEL,1) + TransMat(2,2)*RELM(IEL,2) + TransMat(3,2)*RELM(IEL,3) + XCM(2,L)
!         XMP(3) = TransMat(1,3)*RELM(IEL,1) + TransMat(2,3)*RELM(IEL,2) + TransMat(3,3)*RELM(IEL,3) + XCM(3,L) 

         XMP(1)=RELM(IEL,1)+XCM0(1)
         XMP(2)=RELM(IEL,2)+XCM0(2)
         XMP(3)=RELM(IEL,3)+XCM0(3)
         
         IF (IRR.EQ.0) THEN
           CALL DINP_RGU(XMP,T,DINC)
         ELSE
           CALL DINP_IRGU(XMP,T,DINC)
         ENDIF

         DO J=1,3
           VF(J)=DINC(J)
         ENDDO
         
         CALL NMATRX(DELM(EPN(IEL),:),NMAT)  ! DELM为单位向量
         CALL BRMUL(NMAT,VF,3,3,1,VFT)

         CALL DMATRX(RELM(IEL,:),DMAT)       ! RELM为每个morison单元的质心，计算全局-局部转换矩阵
         CALL BRMUL(NMAT,DMAT,3,3,6,NDMAT)
         CALL BRMUL(NDMAT,VCG,3,6,1,VMT)       
         
         VRT=VFT-VMT
         
         
         AVRT=DSQRT(DABS(VRT(1))**2+DABS(VRT(2))**2+DABS(VRT(3))**2)
!        WRITE(70,120) (VCG(I),I=1,6)
         
         DO J=1,3
          PVRT(J)= AVRT*VRT(J)
         ENDDO
         
         CALL RMATRX(RELM(IEL,:),RMAT) 
         CALL BRMUL(RMAT,PVRT,6,3,1,MVCT) 
         
         FDAMP=FDAMP+1.D0/2.D0*RHOW*CD(IEL)*DCYL(EPN(IEL))*LCYL(EPN(IEL))*MVCT
         
         IF ( IEL >= 220 ) THEN
             CALL NMATRX2(DELM(EPN(IEL),:),NMAT2)  ! DELM为单位向量
             CALL BRMUL(NMAT2,VF,3,3,1,VFT)
         
             CALL BRMUL(NMAT2,DMAT,3,3,6,NDMAT)
             CALL BRMUL(NDMAT,VCG,3,6,1,VMT)       
         
             VRT=VFT-VMT         
         
            AVRT=DSQRT(DABS(VRT(1))**2+DABS(VRT(2))**2+DABS(VRT(3))**2)
            DO J=1,3  
                PVRT(J)= AVRT*VRT(J)
            ENDDO
            CALL RMATRX(RELM(IEL,:),RMAT)
            CALL BRMUL(RMAT,PVRT,6,3,1,MVCT)
            FDAMP2 = 1.D0/2.D0*RHOW*CD5*DCYL(EPN(IEL))*LCYL(EPN(IEL))*MVCT*PI
            FDAMP=FDAMP+ FDAMP2
         ENDIF
         
             
         
   
10     CONTINUE
       
120	  FORMAT(12(2x,E12.5))  
      RETURN
      END SUBROUTINE 
      

!  -----------------------------------------------------------------------------   
!        Read the input mesh of stick model
!-------------------------------------------------------------------------------      
     SUBROUTINE mMESHDA

	  USE STKMesh_mod
	  USE WaveDyn_mod
	  USE Platform_mod      
      IMPLICIT   NONE  

	  INTEGER M,IEL,IPT
!
! ========================================================
!
	  Print *,' In mMESHDA'
!
	    READ(10,*)
	    READ(10,*)         
        DO 10 IEL=1, NSTRE
          READ(10,*) M, EPN(IEL),RELM(IEL,1),RELM(IEL,2),RELM(IEL,3),CD(IEL)
10      CONTINUE
	   Print *,' After 10  M=',M
       
	    READ(10,*) 
	    READ(10,*)         
        DO 20 IPT=1, NPART
          READ(10,*) M, DELM(IPT,1),DELM(IPT,2),DELM(IPT,3)
20      CONTINUE
	   Print *,' After 20  M=',M
       
	    READ(10,*) 
	    READ(10,*)         
        DO 30  IPT=1, NPART
	    READ(10, *) M,NELM(IPT),DCYL(IPT),SCYL(IPT),LCYL(IPT)
30      CONTINUE
	   Print *,' After 30  M=',M          
       
	    READ(10,*) 
	    READ(10,*)         
        DO 40  IPT=1, NPART
	    READ(10, *) M, RSTR(IPT,1),RSTR(IPT,2),RSTR(IPT,3),REND(IPT,1),REND(IPT,2),REND(IPT,3)
40      CONTINUE
	   Print *,' After 40  M=',M 

!  Ajust the position vector for each element corresponding to the CG of platform
       
        DO 50 IEL=1,NSTRE
          RELM(IEL,1)=RELM(IEL,1)-XCM0(1)
          RELM(IEL,2)=RELM(IEL,2)-XCM0(2)
          RELM(IEL,3)=RELM(IEL,3)-XCM0(3)
50      CONTINUE

        DO 60 IPT=1,NPART
          RSTR(IPT,1)=RSTR(IPT,1)-XCM0(1)
          RSTR(IPT,2)=RSTR(IPT,2)-XCM0(2)
          RSTR(IPT,3)=RSTR(IPT,3)-XCM0(3)
          
          REND(IPT,1)=REND(IPT,1)-XCM0(1)
          REND(IPT,2)=REND(IPT,2)-XCM0(2)
          REND(IPT,3)=REND(IPT,3)-XCM0(3)
60      CONTINUE  
        
!  Give the inertia and drag coefficient for each element

!        DO 70 IEL=1,NSTRE        
!         IF (IDL(EPN(IEL)).EQ.1) THEN
!             CM(IEL)=1.5D0
!             CD(IEL)=CD1
!         ELSEIF (IDL(EPN(IEL)).EQ.2) THEN
!             CM(IEL)=1.5D0
!             CD(IEL)=CD2
!         ELSEIF (IDL(EPN(IEL)).EQ.3) THEN
!             CM(IEL)=2.2D0
!             CD(IEL)=CD3
!         ELSEIF (IDL(EPN(IEL)).EQ.4) THEN
!             CM(IEL)=1.5D0
!             CD(IEL)=CD4
!         ELSEIF (IDL(EPN(IEL)).EQ.5) THEN
!             CM(IEL)=1.5D0
!             CD(IEL)=CD5
!         ELSE
!             PRINT*,'No such type of element.Please check the input mesh.'
!         ENDIF        
!70      CONTINUE
        
       RETURN
     END SUBROUTINE mMESHDA      
     
     
!--------------------------------------------------------
!      Calculate global-local transformation matrix 
!--------------------------------------------------------
!
        SUBROUTINE DMATRX(LOCT,DMAT)

        IMPLICIT  NONE

	    REAL*8,INTENT(IN):: LOCT(3)
	    REAL*8,INTENT(OUT):: DMAT(3,6) 
	    REAL*8 X,Y,Z
! 
        
        X=LOCT(1)
        Y=LOCT(2)        
        Z=LOCT(3) 
         
        DMAT=0.D0
        DMAT(1,1)= 1.D0
        DMAT(2,2)= 1.D0 
        DMAT(3,3)= 1.D0
        
        DMAT(1,5)=  Z
        DMAT(1,6)= -Y        
        DMAT(2,4)= -Z
        DMAT(2,6)=  X        
        DMAT(3,4)=  Y
        DMAT(3,5)= -X
        
        RETURN
        END         
!       
!---------------------------------------------------------------------
!           Calculate normal matrix of the Morison member
!---------------------------------------------------------------------
!
        SUBROUTINE NMATRX(AXDV,NMAT)

        IMPLICIT  NONE

	    REAL*8,INTENT(IN):: AXDV(3)
	    REAL*8,INTENT(OUT):: NMAT(3,3) 
	    REAL*8 LX,LY,LZ
        ! 单位向量的X,Y,Z
        LX=AXDV(1)  
        LY=AXDV(2)
        LZ=AXDV(3)         
!  
        NMAT(1,1)= 1.D0-LX**2
        NMAT(2,2)= 1.D0-LY**2 
        NMAT(3,3)= 1.D0-LZ**2
        
        NMAT(1,2)= -LX*LY
        NMAT(1,3)= -LX*LZ
        NMAT(2,1)= -LX*LY
        NMAT(2,3)= -LY*LZ    
        NMAT(3,1)= -LZ*LX
        NMAT(3,2)= -LY*LZ  
        
        RETURN
        END 
        
!---------------------------------------------------------------------
!   Calculate moment matrix of the Morison member with respect to CG
!---------------------------------------------------------------------
!
        SUBROUTINE RMATRX(LOCT,RMAT)

        IMPLICIT  NONE

	    REAL*8,INTENT(IN):: LOCT(3)
	    REAL*8,INTENT(OUT):: RMAT(6,3) 
	    REAL*8 X,Y,Z       
        
        X=LOCT(1)
        Y=LOCT(2)        
        Z=LOCT(3)         
!  
        RMAT=0.D0
        RMAT(1,1)= 1.D0
        RMAT(2,2)= 1.D0
        RMAT(3,3)= 1.D0
        
        RMAT(4,2)= -Z
        RMAT(4,3)=  Y        
        RMAT(5,1)=  Z
        RMAT(5,3)= -X        
        RMAT(6,1)= -Y
        RMAT(6,2)=  X
        
!        RMAT(1,2)= -Z
!        RMAT(1,3)=  Y        
!        RMAT(2,1)=  Z
!        RMAT(2,3)= -X        
!        RMAT(3,1)= -Y
!        RMAT(3,2)=  X        
        RETURN
        END             
!       
!---------------------------------------------------------------------
!           Calculate normal matrix of the Morison member(axial direction)
!---------------------------------------------------------------------
!     
        SUBROUTINE NMATRX2(AXDV,NMAT2)

        IMPLICIT  NONE

	    REAL*8,INTENT(IN):: AXDV(3)
	    REAL*8,INTENT(OUT):: NMAT2(3,3) 
	    REAL*8 LX,LY,LZ
        ! 单位向量的X,Y,Z
        LX=AXDV(1)  
        LY=AXDV(2)
        LZ=AXDV(3)         
!  
        NMAT2(1,1)= LX**2
        NMAT2(2,2)= LY**2 
        NMAT2(3,3)= LZ**2
        
        NMAT2(1,2)= LX*LY
        NMAT2(1,3)= LX*LZ
        NMAT2(2,1)= LX*LY
        NMAT2(2,3)= LY*LZ    
        NMAT2(3,1)= LZ*LX
        NMAT2(3,2)= LY*LZ  
        
        RETURN
        END     