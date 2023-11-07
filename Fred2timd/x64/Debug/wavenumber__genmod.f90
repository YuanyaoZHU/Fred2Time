        !COMPILER-GENERATED INTERFACE MODULE: Tue Nov  7 10:00:44 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE WAVENUMBER__genmod
          INTERFACE 
            FUNCTION WAVENUMBER(OMEGA,G,DEPTH)
              REAL(KIND=8), INTENT(IN) :: OMEGA
              REAL(KIND=8), INTENT(IN) :: G
              REAL(KIND=8), INTENT(IN) :: DEPTH
              REAL(KIND=8) :: WAVENUMBER
            END FUNCTION WAVENUMBER
          END INTERFACE 
        END MODULE WAVENUMBER__genmod
