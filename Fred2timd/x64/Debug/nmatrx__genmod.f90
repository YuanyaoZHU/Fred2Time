        !COMPILER-GENERATED INTERFACE MODULE: Tue Nov  7 10:00:45 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE NMATRX__genmod
          INTERFACE 
            SUBROUTINE NMATRX(AXDV,NMAT)
              REAL(KIND=8), INTENT(IN) :: AXDV(3)
              REAL(KIND=8), INTENT(OUT) :: NMAT(3,3)
            END SUBROUTINE NMATRX
          END INTERFACE 
        END MODULE NMATRX__genmod
