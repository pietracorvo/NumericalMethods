MODULE CONSTANTS
   IMPLICIT NONE
   INTEGER, PARAMETER   :: p = 8 ! precision parameter
   REAL(kind=p),parameter :: Pi = 3.1415926_p
!  INSERT GRAVITATIONAL CONSTANT IN CORRECT UNITS
   REAL(KIND=p), PARAMETER :: G=8.811358E-10     ! gravitational constant in AU³/Me/D^2
!  INSERT SOLAR MASS IN CORRECT UNITS
   REAL(KIND=p), PARAMETER :: Ms=0.33305E06      ! mass of sun in units of earth mass
   REAL(KIND=p), PARAMETER :: Me=1.0_p           ! mass of earth in units of earth mass
   REAL(KIND=p), PARAMETER :: Re=1.0_p           ! radius in AU of approximate circular orbit of earth in solar system

   ! BONUS
   ! More constants if you want to play around 
   REAL(KIND=p), PARAMETER :: Mm=0.1745_p        ! mass of mars in units of earth mass
   REAL(KIND=p), PARAMETER :: Mmoon=0.012300_p   ! mass of moon in units of earth mass
   REAL(KIND=p), PARAMETER :: Rm=1.52371_p       ! radius of approximate circular orbit of mars in AU
   REAL(KIND=p), PARAMETER :: Ope=365.25         ! orbital period of earth in days
   REAL(KIND=p), PARAMETER :: Opm=686.97         ! orbital period of mars in days
END MODULE CONSTANTS

