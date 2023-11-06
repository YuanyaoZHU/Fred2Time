        !COMPILER-GENERATED INTERFACE MODULE: Fri Mar 17 20:19:43 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE INVERSE__genmod
          INTERFACE 
            SUBROUTINE INVERSE(A,Y,N,NP)
              INTEGER(KIND=4) :: NP
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: A(NP,NP)
              REAL(KIND=8) :: Y(NP,NP)
            END SUBROUTINE INVERSE
          END INTERFACE 
        END MODULE INVERSE__genmod
