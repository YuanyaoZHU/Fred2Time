        !COMPILER-GENERATED INTERFACE MODULE: Fri Mar 17 20:19:43 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE TRANS1TO6__genmod
          INTERFACE 
            SUBROUTINE TRANS1TO6(YAWANG,TILTANG,WINGANG,CONEANG,THETA1, &
     &THETA2,THETA3,TRANSMAT)
              REAL(KIND=8), INTENT(IN) :: YAWANG
              REAL(KIND=8), INTENT(IN) :: TILTANG
              REAL(KIND=8), INTENT(IN) :: WINGANG
              REAL(KIND=8), INTENT(IN) :: CONEANG
              REAL(KIND=8), INTENT(IN) :: THETA1
              REAL(KIND=8), INTENT(IN) :: THETA2
              REAL(KIND=8), INTENT(IN) :: THETA3
              REAL(KIND=8), INTENT(OUT) :: TRANSMAT(3,3)
            END SUBROUTINE TRANS1TO6
          END INTERFACE 
        END MODULE TRANS1TO6__genmod
