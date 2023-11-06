        !COMPILER-GENERATED INTERFACE MODULE: Fri Mar 17 20:19:43 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE PLATFORMVELOCITY__genmod
          INTERFACE 
            SUBROUTINE PLATFORMVELOCITY(ROLL,PITCH,YAW,VSURGE,VSWAY,    &
     &VHEAVE,VROLL,VPITCH,VYAW,RCG,RH,VLCT)
              REAL(KIND=8), INTENT(IN) :: ROLL
              REAL(KIND=8), INTENT(IN) :: PITCH
              REAL(KIND=8), INTENT(IN) :: YAW
              REAL(KIND=8), INTENT(IN) :: VSURGE
              REAL(KIND=8), INTENT(IN) :: VSWAY
              REAL(KIND=8), INTENT(IN) :: VHEAVE
              REAL(KIND=8), INTENT(IN) :: VROLL
              REAL(KIND=8), INTENT(IN) :: VPITCH
              REAL(KIND=8), INTENT(IN) :: VYAW
              REAL(KIND=8), INTENT(IN) :: RCG
              REAL(KIND=8), INTENT(IN) :: RH
              REAL(KIND=8), INTENT(OUT) :: VLCT(3)
            END SUBROUTINE PLATFORMVELOCITY
          END INTERFACE 
        END MODULE PLATFORMVELOCITY__genmod
