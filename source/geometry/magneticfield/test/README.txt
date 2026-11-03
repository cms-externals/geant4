run tests using ./test.sh

to successfuly run tests you have:
1. build Geant4 with -DGEANT4_BUILD_TESTS=ON flag
2. set enviroment variables
    G4BUILD_DIR=/path/to/Geant4/build/
    G4BIN_CMAKE=${G4BUILD_DIR}/BuildProducts/bin 


Stepper/driver combinations for use with testProElectroMagField



$G4BIN_CMAKE/testProElectroMagField  stepper_driver_id 


- Templated steppers for use with simple (later templated) equation of motion

New (October 2020):
108: Templated Cash/Karp      5th-order (i.e. templated version of G4CashKarpRKF45 
345: Templated Dormand/Prince 5th-order (i.e. templated version of G4DormandPrince745 )

- Use of FSAL integration driver:

241: FSAL Integration driver with  G4RKF547FEq1  ( 5th order 'equilibrium' RK of Higham & Hall 
242: FSAL Integration driver with  G4RKF547FEq2  ( similar / stepper 2 of H & H )
243: FSAL Integration driver with  G4RKF547FEq3  ( similar / stepper 3 of H & H )

New (Spring 2020?)
245: FSAL Integration driver using Dormand/Prince  ( G4DormandPrince745 )
