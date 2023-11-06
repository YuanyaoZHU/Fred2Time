        !COMPILER-GENERATED INTERFACE MODULE: Fri May 30 12:00:56 2014
        MODULE LUBKSB__genmod
          INTERFACE 
            SUBROUTINE LUBKSB(A,N,NP,INDX,B)
              INTEGER(KIND=4) :: NP
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: A(NP,NP)
              INTEGER(KIND=4) :: INDX(N)
              REAL(KIND=8) :: B(N)
            END SUBROUTINE LUBKSB
          END INTERFACE 
        END MODULE LUBKSB__genmod
