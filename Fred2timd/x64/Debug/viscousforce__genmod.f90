        !COMPILER-GENERATED INTERFACE MODULE: Fri Mar 17 20:19:43 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE VISCOUSFORCE__genmod
          INTERFACE 
            SUBROUTINE VISCOUSFORCE(L,T,IRR,ROLL,PITCH,YAW,VCG)
              INTEGER(KIND=4) :: L
              REAL(KIND=8) :: T
              INTEGER(KIND=4) :: IRR
              REAL(KIND=8) :: ROLL
              REAL(KIND=8) :: PITCH
              REAL(KIND=8) :: YAW
              REAL(KIND=8) :: VCG(6)
            END SUBROUTINE VISCOUSFORCE
          END INTERFACE 
        END MODULE VISCOUSFORCE__genmod
