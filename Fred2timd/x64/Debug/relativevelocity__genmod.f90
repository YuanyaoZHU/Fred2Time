        !COMPILER-GENERATED INTERFACE MODULE: Fri Mar 17 20:19:43 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE RELATIVEVELOCITY__genmod
          INTERFACE 
            SUBROUTINE RELATIVEVELOCITY(YAWANG,TILTANG,WINGANG,CONEANG, &
     &ROLL,PITCH,YAW,VSURGE,VSWAY,VHEAVE,VROLL,VPITCH,VYAW,XPC,XYC,XRC, &
     &ZHUB,VHUB,WINDDRT,NBEM,NBL,RELM,CT,RTB,PLTFSDLNTH,URL,IFWAKE,     &
     &IFCOUP,IFHEADING,VRLTV)
              INTEGER(KIND=4), INTENT(IN) :: NBL
              INTEGER(KIND=4), INTENT(IN) :: NBEM
              REAL(KIND=8), INTENT(IN) :: YAWANG
              REAL(KIND=8), INTENT(IN) :: TILTANG
              REAL(KIND=8), INTENT(IN) :: WINGANG
              REAL(KIND=8), INTENT(IN) :: CONEANG
              REAL(KIND=8), INTENT(IN) :: ROLL
              REAL(KIND=8), INTENT(IN) :: PITCH
              REAL(KIND=8), INTENT(IN) :: YAW
              REAL(KIND=8), INTENT(IN) :: VSURGE
              REAL(KIND=8), INTENT(IN) :: VSWAY
              REAL(KIND=8), INTENT(IN) :: VHEAVE
              REAL(KIND=8), INTENT(IN) :: VROLL
              REAL(KIND=8), INTENT(IN) :: VPITCH
              REAL(KIND=8), INTENT(IN) :: VYAW
              REAL(KIND=8), INTENT(IN) :: XPC(3)
              REAL(KIND=8), INTENT(IN) :: XYC(3)
              REAL(KIND=8), INTENT(IN) :: XRC(3)
              REAL(KIND=8), INTENT(IN) :: ZHUB
              REAL(KIND=8), INTENT(IN) :: VHUB
              REAL(KIND=8), INTENT(IN) :: WINDDRT
              REAL(KIND=8), INTENT(IN) :: RELM(NBEM)
              REAL(KIND=8), INTENT(IN) :: CT
              REAL(KIND=8), INTENT(IN) :: RTB
              REAL(KIND=8), INTENT(IN) :: PLTFSDLNTH
              REAL(KIND=8), INTENT(IN) :: URL(3)
              LOGICAL(KIND=4), INTENT(IN) :: IFWAKE
              LOGICAL(KIND=4), INTENT(IN) :: IFCOUP
              LOGICAL(KIND=4), INTENT(IN) :: IFHEADING
              REAL(KIND=8), INTENT(OUT) :: VRLTV(NBEM,NBL)
            END SUBROUTINE RELATIVEVELOCITY
          END INTERFACE 
        END MODULE RELATIVEVELOCITY__genmod
