        !COMPILER-GENERATED INTERFACE MODULE: Fri Mar 17 20:19:43 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ROTORCENTERLOCATION__genmod
          INTERFACE 
            SUBROUTINE ROTORCENTERLOCATION(YAWANG,TILTANG,ROLL,PITCH,YAW&
     &,XPC,XYC,XRC,XYZ1)
              REAL(KIND=8), INTENT(IN) :: YAWANG
              REAL(KIND=8), INTENT(IN) :: TILTANG
              REAL(KIND=8), INTENT(IN) :: ROLL
              REAL(KIND=8), INTENT(IN) :: PITCH
              REAL(KIND=8), INTENT(IN) :: YAW
              REAL(KIND=8), INTENT(IN) :: XPC(3)
              REAL(KIND=8), INTENT(IN) :: XYC(3)
              REAL(KIND=8), INTENT(IN) :: XRC(3)
              REAL(KIND=8), INTENT(OUT) :: XYZ1(3)
            END SUBROUTINE ROTORCENTERLOCATION
          END INTERFACE 
        END MODULE ROTORCENTERLOCATION__genmod
