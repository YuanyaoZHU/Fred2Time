        !COMPILER-GENERATED INTERFACE MODULE: Fri Mar 17 20:21:20 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE LARGROTTRANS__genmod
          INTERFACE 
            SUBROUTINE LARGROTTRANS(THETA1,THETA2,THETA3,TRANSMAT)
              REAL(KIND=8), INTENT(IN) :: THETA1
              REAL(KIND=8), INTENT(IN) :: THETA2
              REAL(KIND=8), INTENT(IN) :: THETA3
              REAL(KIND=8), INTENT(OUT) :: TRANSMAT(3,3)
            END SUBROUTINE LARGROTTRANS
          END INTERFACE 
        END MODULE LARGROTTRANS__genmod
