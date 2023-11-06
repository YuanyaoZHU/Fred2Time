        !COMPILER-GENERATED INTERFACE MODULE: Fri Mar 17 20:19:43 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE FILON_LIJ__genmod
          INTERFACE 
            SUBROUTINE FILON_LIJ(FX,LEN,T,LB,UB,OPT,EPS,VAL)
              INTEGER(KIND=4), INTENT(IN) :: LEN
              REAL(KIND=8), INTENT(IN) :: FX(0:LEN)
              REAL(KIND=8), INTENT(IN) :: T
              REAL(KIND=8), INTENT(IN) :: LB
              REAL(KIND=8), INTENT(IN) :: UB
              INTEGER(KIND=4), INTENT(IN) :: OPT
              REAL(KIND=8), INTENT(OUT) :: EPS
              REAL(KIND=8), INTENT(OUT) :: VAL
            END SUBROUTINE FILON_LIJ
          END INTERFACE 
        END MODULE FILON_LIJ__genmod
