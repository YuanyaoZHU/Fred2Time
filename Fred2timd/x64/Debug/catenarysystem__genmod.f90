        !COMPILER-GENERATED INTERFACE MODULE: Fri Mar 17 20:21:20 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CATENARYSYSTEM__genmod
          INTERFACE 
            SUBROUTINE CATENARYSYSTEM(STAGETYPE,NUMLINES,X0,X,F_LINES,  &
     &LFAIRTE)
              INTEGER(KIND=4), INTENT(IN) :: NUMLINES
              LOGICAL(KIND=4), INTENT(IN) :: STAGETYPE
              REAL(KIND=8), INTENT(IN) :: X0(6)
              REAL(KIND=8), INTENT(IN) :: X(6)
              REAL(KIND=8), INTENT(OUT) :: F_LINES(6)
              REAL(KIND=8) :: LFAIRTE(NUMLINES)
            END SUBROUTINE CATENARYSYSTEM
          END INTERFACE 
        END MODULE CATENARYSYSTEM__genmod
