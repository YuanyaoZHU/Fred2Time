        !COMPILER-GENERATED INTERFACE MODULE: Fri Mar 17 20:19:43 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE BRMUL__genmod
          INTERFACE 
            SUBROUTINE BRMUL(A,B,M,N,K,C)
              INTEGER(KIND=4), INTENT(IN) :: K
              INTEGER(KIND=4), INTENT(IN) :: N
              INTEGER(KIND=4), INTENT(IN) :: M
              REAL(KIND=8), INTENT(IN) :: A(1:M,1:N)
              REAL(KIND=8), INTENT(IN) :: B(1:N,1:K)
              REAL(KIND=8), INTENT(OUT) :: C(1:M,1:K)
            END SUBROUTINE BRMUL
          END INTERFACE 
        END MODULE BRMUL__genmod
