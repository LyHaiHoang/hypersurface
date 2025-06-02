PROGRAM DistributePointsCircle
  IMPLICIT NONE
  INTEGER, PARAMETER :: N = 12
  REAL(8), PARAMETER :: pi = 3.141592653589793d0
  REAL(8) :: x(N), y(N)
  INTEGER :: i
  REAL(8) :: theta

  PRINT *, "Phan bo ", N, " diem tren duong tron don vi:"

  DO i = 1, N
     theta = 2.0d0 * pi * (i - 1) / N
     x(i) = COS(theta)
     y(i) = SIN(theta)
     PRINT "(I2,2F12.6)", i, x(i), y(i)
  END DO

END PROGRAM DistributePointsCircle
