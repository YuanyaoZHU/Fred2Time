        !COMPILER-GENERATED INTERFACE MODULE: Fri May 30 12:00:54 2014
        MODULE FILON_KIJ__genmod
          INTERFACE 
            SUBROUTINE FILON_KIJ(FX,LEN,T,LB,UB,OPT,VAL)
              INTEGER(KIND=4), INTENT(IN) :: LEN
              REAL(KIND=8), INTENT(IN) :: FX(0:LEN)
              REAL(KIND=8), INTENT(IN) :: T
              REAL(KIND=8), INTENT(IN) :: LB
              REAL(KIND=8), INTENT(IN) :: UB
              INTEGER(KIND=4), INTENT(IN) :: OPT
              REAL(KIND=8), INTENT(OUT) :: VAL
            END SUBROUTINE FILON_KIJ
          END INTERFACE 
        END MODULE FILON_KIJ__genmod
