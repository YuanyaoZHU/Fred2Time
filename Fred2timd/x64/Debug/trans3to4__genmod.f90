        !COMPILER-GENERATED INTERFACE MODULE: Fri Mar 17 20:19:43 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE TRANS3TO4__genmod
          INTERFACE 
            SUBROUTINE TRANS3TO4(YAWANG,TILTANG,TRANSMAT)
              REAL(KIND=8), INTENT(IN) :: YAWANG
              REAL(KIND=8), INTENT(IN) :: TILTANG
              REAL(KIND=8), INTENT(OUT) :: TRANSMAT(3,3)
            END SUBROUTINE TRANS3TO4
          END INTERFACE 
        END MODULE TRANS3TO4__genmod
