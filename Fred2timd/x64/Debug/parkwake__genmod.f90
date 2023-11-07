        !COMPILER-GENERATED INTERFACE MODULE: Tue Nov  7 10:00:44 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE PARKWAKE__genmod
          INTERFACE 
            SUBROUTINE PARKWAKE(IFHEADING,IFWAKE,RR,X01,ZHUB,U0,CT,R0,  &
     &DV01)
              LOGICAL(KIND=4), INTENT(IN) :: IFHEADING
              LOGICAL(KIND=4), INTENT(IN) :: IFWAKE
              REAL(KIND=8) :: RR
              REAL(KIND=8) :: X01
              REAL(KIND=8) :: ZHUB
              REAL(KIND=8) :: U0
              REAL(KIND=8) :: CT
              REAL(KIND=8) :: R0
              REAL(KIND=8) :: DV01
            END SUBROUTINE PARKWAKE
          END INTERFACE 
        END MODULE PARKWAKE__genmod
