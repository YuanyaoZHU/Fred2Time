        !COMPILER-GENERATED INTERFACE MODULE: Fri May 30 12:00:56 2014
        MODULE LUDCMP__genmod
          INTERFACE 
            SUBROUTINE LUDCMP(A,N,NP,INDX,D)
              INTEGER(KIND=4) :: NP
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: A(NP,NP)
              INTEGER(KIND=4) :: INDX(N)
              REAL(KIND=8) :: D
            END SUBROUTINE LUDCMP
          END INTERFACE 
        END MODULE LUDCMP__genmod
