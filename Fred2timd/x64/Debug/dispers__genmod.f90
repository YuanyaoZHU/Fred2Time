        !COMPILER-GENERATED INTERFACE MODULE: Tue Nov  7 10:00:44 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DISPERS__genmod
          INTERFACE 
            SUBROUTINE DISPERS(WAVENO,MNO,W,H)
              REAL(KIND=8), INTENT(OUT) :: WAVENO(800)
              INTEGER(KIND=4), INTENT(IN) :: MNO
              REAL(KIND=8), INTENT(IN) :: W
              REAL(KIND=8), INTENT(IN) :: H
            END SUBROUTINE DISPERS
          END INTERFACE 
        END MODULE DISPERS__genmod
