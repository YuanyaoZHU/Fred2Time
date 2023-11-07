        !COMPILER-GENERATED INTERFACE MODULE: Tue Nov  7 10:00:44 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CURVE__genmod
          INTERFACE 
            SUBROUTINE CURVE(THRUST,MTB,FILENAME,VREL,LINE,COLUMN,WROT, &
     &CT,CP)
              INTEGER(KIND=4) :: LINE
              REAL(KIND=8) :: THRUST
              REAL(KIND=8) :: MTB
              CHARACTER(LEN=20) :: FILENAME
              REAL(KIND=8) :: VREL
              INTEGER(KIND=4) :: COLUMN
              REAL(KIND=8) :: WROT
              REAL(KIND=8) :: CT
              REAL(KIND=8) :: CP
            END SUBROUTINE CURVE
          END INTERFACE 
        END MODULE CURVE__genmod
