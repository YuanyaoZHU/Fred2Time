        !COMPILER-GENERATED INTERFACE MODULE: Tue Nov  7 10:32:06 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SMLLROTTRANS__genmod
          INTERFACE 
            SUBROUTINE SMLLROTTRANS(THETA1,THETA2,THETA3,TRANSMAT)
              REAL(KIND=8), INTENT(IN) :: THETA1
              REAL(KIND=8), INTENT(IN) :: THETA2
              REAL(KIND=8), INTENT(IN) :: THETA3
              REAL(KIND=8), INTENT(OUT) :: TRANSMAT(3,3)
            END SUBROUTINE SMLLROTTRANS
          END INTERFACE 
        END MODULE SMLLROTTRANS__genmod
