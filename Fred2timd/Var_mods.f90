        MODULE CABLE
        REAL(8) :: RFair, RAnch, DepthAnch, DepthFair
        REAL(8) :: LineDiam,LineMassDen,LineUnstrLen,LineSeabedCD,LineEAStff    !锚链密度（干重）;直径
        REAL*8,ALLOCATABLE:: LAngle(:)


        END MODULE CABLE
    
    
    
    !---------------------------------------------------------------------------------------------
!        Data module for declaring variables in stick mesh generator
!---------------------------------------------------------------------------------------------				
!
        MODULE STKMesh_mod
!
	    INTEGER NPART,NSTRD,NSTRE,NBRAC,NCLMN,NCLMN1,NCLMN2,NPOTN,NCONC
!     
	    INTEGER,ALLOCATABLE:: EPN(:),IDL(:),NELM(:)
        REAL*8,ALLOCATABLE:: RELM(:,:),DELM(:,:),RNOD(:,:),DCYL(:),SCYL(:),LCYL(:)
        REAL*8,ALLOCATABLE:: RSTR(:,:),REND(:,:)
!
        END MODULE STKMesh_mod
        
!---------------------------------------------------------------------------------------------
!        Data module for declaring wave relevant variables in time domain simulation
!---------------------------------------------------------------------------------------------	
!
        MODULE WaveDyn_mod
!
        REAL*8 H,HS,BETA
        REAL*8 WK,W1,WL,TP
!
        REAL*8,ALLOCATABLE:: CM(:),CD(:)
        REAL*8,ALLOCATABLE:: AMAS(:,:,:),BDMP(:,:,:),HAM(:,:,:),HBM(:,:,:),RFK(:,:,:)
        REAL*8,ALLOCATABLE:: EFOR(:,:,:),FWAVE(:,:),LTF(:,:),ETA(:) ,RK(:,:,:),AW(:),RK1(:,:,:)
        REAL*8 :: CD1,CD2,CD3,CD4,CD5
        REAL*8 :: CDX(6)
      
         REAL*8 FRTAR(6),FHDRS(6),FDAMP(6),FPLFM(6),FGRVT(6),FLINE(6),FLINE0(6),FILINE(6),AINF(6,6),FBOYC(6),FDAMP2(6)
!
        END MODULE WaveDyn_mod

        
!---------------------------------------------------------------------------------------------
!        Data module for declaring body relevant variables in Morison_stick numerical model
!---------------------------------------------------------------------------------------------	
!
        MODULE Platform_mod
       
        REAL*8 VOL,MASS,XCM0(6),BYC
!
        REAL*8 IB(3,3),MATX(6,6),CRS(6,6),KSPR(6,6),BLNR(6,6),LRS(6,6),ZDAMP(6,6),FORCE_ZHU(6,1),MOTION(6,2),STIFF(6,6)
        REAL*8,ALLOCATABLE:: XCM(:,:),X(:,:),DX(:,:)
        REAL*8 MassRna,xRNA
        REAL*8 CWT,ZWT,CWF,ZWF,CCF,ZCF
!
        END MODULE Platform_mod
        
        
!---------------------------------------------------------------------------------------------
!        Data module for declaring variables in aerodynamic subroutines
!---------------------------------------------------------------------------------------------						
!
        MODULE Blade_mod
!       
        LOGICAL IFDUCT,IFWAKE,Ifheading,IFVISC, IFCOUP
        INTEGER NBEM,NBL,NFOIL(100)
        INTEGER NDU21,NDU25,NDU30,NDU35,NDU40,NNACA64,NCYLD
        PARAMETER (NBEM=17,NBL=3,NDU21=142,NDU25=141,NDU30=143)  
        PARAMETER (NDU35=135,NDU40=136,NNACA64=127,NCYLD=3)

        REAL*8  RTB,RHB,INC0,SHRO
      
        REAL*8  DU21(NDU21,4),DU25(NDU25,4),DU30(NDU30,4),DU35(NDU35,4) 
        REAL*8  DU40(NDU40,4),NACA64(NNACA64,4),CYLD1(NCYLD,4),CYLD2(NCYLD,4)

        REAL*8  RND(NBEM),TWST(NBEM),DRND(NBEM),CHRD(NBEM),SDT(NBEM),XR(NBEM),A0(NBEM),VRltv(NBEM,NBL)
        
        REAL*8,ALLOCATABLE:: WDSP(:,:),WTSP(:,:),PITH(:,:),ACROT(:,:)
        REAL*8,ALLOCATABLE:: FAERO(:,:),TAERO(:,:),MCPWR(:,:),FWIND(:,:),COT(:,:),COP(:,:)
   
        END MODULE Blade_mod
        
!---------------------------------------------------------------------------------------------
!        Data module for declaring variables of general computation
!---------------------------------------------------------------------------------------------						
!
        MODULE Const_mod
!     
        REAL*8 G,PI,RHOW,RHOA,EPS
!
        COMPLEX*16 CI
!
        DATA G,PI,EPS/9.807d0, 3.141592653589793d0, 1.E-8/  
        DATA CI/(0.0D0, 1.0D0)/
!
    END MODULE Const_mod
    
!---------------------------------------------------------------------------------------------
!       Data module for ODE calculation
!---------------------------------------------------------------------------------------------

    
