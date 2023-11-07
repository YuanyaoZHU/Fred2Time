        !COMPILER-GENERATED INTERFACE MODULE: Tue Nov  7 10:00:44 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE AERO_STEP__genmod
          INTERFACE 
            SUBROUTINE AERO_STEP(VWND,WROT,PTAG,THRUST,MTB,POWER,CT,CP)
              REAL(KIND=8) :: VWND(3)
              REAL(KIND=8) :: WROT
              REAL(KIND=8) :: PTAG
              REAL(KIND=8) :: THRUST
              REAL(KIND=8) :: MTB
              REAL(KIND=8) :: POWER
              REAL(KIND=8) :: CT
              REAL(KIND=8) :: CP
            END SUBROUTINE AERO_STEP
          END INTERFACE 
        END MODULE AERO_STEP__genmod
