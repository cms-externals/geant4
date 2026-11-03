//
//  Test GUIntegrationDriver
//   * compares with the output of a reference stepper (high accuracy) - ok for small steps
//
//  Based on testStepperFixed.cc
//    which was started from the work of Somnath Banerjee in GSoC 2015
//

#include "G4FieldTrack.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"

#include <iomanip>

using CLHEP::degree;
using CLHEP::eplus;
using CLHEP::meter;
using CLHEP::millimeter;
using CLHEP::second;
using CLHEP::tesla;

using ThreeVector = G4ThreeVector;  //  vecgeom::Vector3D<double> ThreeVector;

// #include "G4VIntegrationDriver.h"
#include "G4QSSDriver.hh"

// #define  COMPARE_TO_G4  1

#include "G4PhysicalConstants.hh"

// #include <system.h>
// #include "G4Types.hh"
#include "G4CachedMagneticField.hh"
#include "G4DormandPrince745.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4TCachedMagneticField.hh"
#include "G4TDormandPrince45.hh"
#include "G4TMagFieldEquation.hh"
#include "G4TUniformMagneticField.hh"
#include "G4UniformMagField.hh"
// #include "G4ExactHelixStepper.hh"
// #include "G4ClassicalRK4.hh"

#include "G4InterpolationDriver.hh"
#include "G4QSSDriverCreator.hh"

// #include "G4.hh"

using namespace std;
// using namespace CLHEP;

constexpr G4double gDistConstField = 0.05 * millimeter;

#ifndef QUAD_FIELD

G4UniformMagField uniformMagField(10. * tesla, 0., 0.);
G4String fieldName("Uniform 10 Tesla");

// Fields for the two drivers - separate caches for number of calls
G4CachedMagneticField simpleField(&uniformMagField, gDistConstField);  // For test/Qss

G4TUniformMagneticField tUniformMagField(10. * tesla, 0., 0.);
using Field_t = G4TCachedMagneticField<G4TUniformMagneticField>;  // <G4UniformMagField>;
Field_t tMagField(&tUniformMagField, gDistConstField);  // For ref/DoPri5

// G4CachedMagneticField  simpleField( &uniformMagField, gDistConstField );
#else
using Field_t = G4TCachedMagneticField<G4QuadrupoleMagField>;

// G4RotationMatrix       rotM( 0.1, 0.2, 0.05 );  // RotationMatrix( phi, theta, psi );
// G4ThreeVector          origin( 5.0 * millimeter, 2.5 * millimeter, 1.25 * millimeter );
// G4QuadrupoleMagField  quadrupoleMagField( 1.*tesla/meter, origin, &rotM );
G4QuadrupoleMagField quadrupoleMagField(10. * tesla / meter);

G4String fieldName("Cached Quadropole field, 1 T/meter, cache=1cm");

G4CachedMagneticField simpleField(&simpleQuadrupoleMagField, gDistConstField);
Field_t tMagField(quadrupoleMagField, 1.0 * cm);
#endif

constexpr int nvar = 6;  // Number of integration variables
// using StepperCK45_t   = G4TCashKarpRKF45<Equation_t, nvar>;

using Equation_t = G4TMagFieldEquation<Field_t>;
using SimpleEquation_t = G4Mag_UsualEqRhs;

using StepperDoPri5_t = G4TDormandPrince45<Equation_t, nvar>;

int main(int /*argc*/, char** /*args*/)
{
  /* -----------------------------SETTINGS-------------------------------- */

  /* Parameters of test
   - Modify values  */

  int no_of_steps = 20;  // No. of Steps for the stepper
  // int stepper_no     = 5;    // Choose stepper no., for refernce see above
  double step_len_mm = 1.;  // 50, 100 // meant as millimeter;  //Step length
  double z_field_in = DBL_MAX;
  double epsTolInp = -1.0;

  double step_len = step_len_mm * millimeter;

  // Set Charge etc.
  double particleCharge = +1.0;  // in e+ units

  // Choice of output coordinates
  int columns[] = {
    0, 1, 1,  // position  x, y, z
    0, 1, 1,  // momentum  x, y, z
    0, 0, 0,  // dydx pos  x, y, z
    0, 0, 0  // dydx mom  x, y, z
  };  // Variables in yOut[] & dydx[] we want to display - 0 for No, 1 for yes

  bool printErr = 1,  // Print the predicted Error
    printRef = 1,  // Print the reference Solution
    printDiff = 1;  // Print the difference
  bool printSep = 0;  // separator  '|'
  int printOk = 3;  // "Ok" or "Bad" - for Vec (val>0) & for each component (val>2)
                    //   Format:   1 Vec-Terse:    "Ok" Or "Bad"
                    //   Format:  >1 Vec-Verbose:  "0 Ok" "1 Bad"
                    //            >3 [i]-Verbose:
  bool printRatio = 0;  //  yDiff / yAverage - per component
  bool printRatioVec = 1;  //                   - per vector  (ie Pos or Vector )
  // bool printMagShift = 1;  //  Print  | p | / | p_0 | - 1.0
  bool printInp = 0;  // print the input values
  bool printInpX = 0;  // print the input values for Ref

  const int epsDigits = 5;  //  Digits - for printing the difference

  const unsigned int wd = 12;  // Width for float/double
  const unsigned int wdf = wd + epsDigits - 3;
  const unsigned int wok = 5;  //  2 for (0/1," ") + 3 for "Bad"/"Ok "

  // Set coordinates here
  double x_pos = 0.,  // pos - position  : input unit = mm
    y_pos = 0., z_pos = 0.;
  double x_mom = 0.0,  // mom - momentum  : input unit = GeV / c
    y_mom = 1.0, z_mom = 1.0;

#ifndef QUAD_FIELD
  double x_field = 0.,  // Uniform Magnetic Field (x,y,z)
    y_field = 0.;  //  Tesla // *tesla ;
  double z_field;
  if (z_field_in < DBL_MAX)
    z_field = z_field_in;
  else
    z_field = -1.0 * tesla;

  uniformMagField.SetFieldValue(G4ThreeVector(x_field, y_field, z_field));
#endif

  /*-------------------------PREPARING test STEPPER-----------------------------*/

  // Create a stepper : --- Not relevant for QSS
  //  auto myStepper = new G4DormandPrince745( gvEquation );

  /*-------------------------PREPARING test DRIVER -----------------------------*/
  const double stepMin = 0.0001 * millimeter;
  const double epsTolDef = 1.0e-5;
  const double errAmplif = 1.0 / epsTolDef;  //  Amplification for printing of error

  // int statisticsVerbosity = 1;

  const double epsTol = (epsTolInp < 0.0) ? epsTolDef : epsTolInp;
  cout << "# Driver parameters:  eps_tol= " << epsTol << "  stepMin= " << stepMin << endl;

  SimpleEquation_t* uEquation = new SimpleEquation_t(&tMagField);  // G4Mag_UsualEqRhs(&tMagField);
  // G4Mag_UsualEqRhs *uEquation = new G4Mag_UsualEqRhs(tMagField);

  G4cout << "# Set up field (name) = " << fieldName << G4endl;

  //  QSS driver -- target
  // -----------------------
  // auto integrDriver = G4QSSDriverCreator::CreateQss2Stepper(uEquation);
  int qssOrder = 2;
  G4MagIntegratorStepper* qssStepper = nullptr;

  G4cout << "# Driver = Qss Driver with order = " << qssOrder << G4endl;
  switch (qssOrder)
  {
    case 2:
      qssStepper = G4QSSDriverCreator::CreateQss2Stepper(uEquation);
      break;

    case 3:
      qssStepper = G4QSSDriverCreator::CreateQss3Stepper(uEquation);
      break;

      //===============================================
    default:
      qssStepper = nullptr;
      G4ExceptionDescription ErrorMsg;
      ErrorMsg << " Incorrect Stepper qssOrder requested. Value was id= " << qssOrder
               << " .  Correct choices are: 2, 3 " << G4endl;
      ErrorMsg << " NO replacement stepper chosen! " << G4endl;
      G4Exception("Application::SetupField", "Runtime Error",
                  FatalErrorInArgument,  //  use JustWarning,
                  " Invalid value of stepper qssOrder");
      break;
  }
  auto testDriver = G4QSSDriverCreator::CreateDriver(qssStepper, stepMin);

  // auto testDriver = new G4InterpolationDriver<TemplatedStepper_t>(stepMin, dopri5_stepper);
  /// Ad interim to check InterpolationDriver 'protocol'

  //  Alternative drivers - for comparison / development
  // ------------------------------------------------------
  Equation_t* pEquation = new Equation_t(&tMagField);

  // Creating the baseline stepper
  // auto exactStepper = new StepperCK45_t(qvEquation2);
  // cout << "#  Reference stepper is: GUTCashKarpRKF45 <tEquationType_t,nvar>(gvEquation2);" <<
  // endl;

  using TemplatedStepper_t = G4TDormandPrince45<Equation_t>;
  TemplatedStepper_t* dopri5_stepper = new G4TDormandPrince45<Equation_t>(pEquation);

  auto baseIntegrDrv = new G4IntegrationDriver<TemplatedStepper_t>(stepMin, dopri5_stepper);
  // new G4InterpolationDriver<TemplatedStepper_t>(stepMin, dopri5_stepper);
  //   else          = new G4IntegrationDriver< ...

  auto referenceStepper =
    new G4DormandPrince745(pEquation);  // A separate stepper instance
                                        // Fix the type and value of 'referenceStepper' - was auto &
                                        // exactStepper = dopri5_stepper;
  G4cout << "# Reference Stepper: DoPri5 = Dormand-Prince 7 stage 4/5 order" << G4endl;

  // Initialising coordinates
  // --------------------------------------------------------
  const double mmTry = CLHEP::millimeter;
  const double ppTry = CLHEP::GeV;

  double yIn[] = {x_pos * mmTry, y_pos * mmTry, z_pos * mmTry,
                  x_mom * ppTry, y_mom * ppTry, z_mom * ppTry};

  const double mmRef = mmTry;
  const double ppRef = ppTry;

  std::cout << "# step_len_mm = " << step_len_mm;
  std::cout << " mmRef= " << mmRef << "   ppRef= " << ppRef << std::endl;

  double yInX[] = {x_pos * mmRef, y_pos * mmRef, z_pos * mmRef,
                   x_mom * ppRef, y_mom * ppRef, z_mom * ppRef};

  double stepLengthRef = step_len_mm * mmRef;

  // Empty buckets for results
  double dydx[8] = {0., 0., 0., 0., 0., 0., 0., 0.},  // 2 extra safety buffer
    dydxRef[8] = {0., 0., 0., 0., 0., 0., 0., 0.}, yout[8] = {0., 0., 0., 0., 0., 0., 0., 0.},
         youtX[8] = {0., 0., 0., 0., 0., 0., 0., 0.},
         //  yerr [8] = {0.,0.,0.,0.,0.,0.,0.,0.},
    yerrX[8] = {0., 0., 0., 0., 0., 0., 0., 0.}, yDiff[8] = {0., 0., 0., 0., 0., 0., 0., 0.},
         yAver[8] = {0., 0., 0., 0., 0., 0., 0., 0.};
  /*-----------------------END PREPARING STEPPER---------------------------*/

  cout.precision(3);

  const G4ThreeVector startPosition(yIn[0], yIn[1], yIn[2]);
  const G4ThreeVector startMomentum(yIn[3], yIn[4], yIn[5]);

  const G4double momentumMagInit = startMomentum.mag();
  cout << "# momentumMagInit = " << momentumMagInit << endl;
  const G4double momentumMagSq = startMomentum.mag2();
  const G4double momentumMag = sqrt(momentumMagSq);  // Use later ..

  const G4double restMass_c2 = electron_mass_c2;  // G4Electron::G4Electron()->GetRestMass();

  const G4ThreeVector startMomentumDir = startMomentum.unit();
  const G4double kineticEnergy = momentumMagSq / (sqrt(momentumMagSq + restMass_c2) + restMass_c2);
  // = sqrt( momentumMagSq + restMass_c2 ) - restMass_c2;
  const G4FieldTrack yStart(startPosition, startMomentumDir,
                            0.0,  // curve_len
                            kineticEnergy, restMass_c2, particleCharge);

  G4ThreeVector chkMomentum = yStart.GetMomentum();
  G4cout << "Field Track: check of momentum vec value p= " << chkMomentum
         << " change from input = " << chkMomentum - startMomentum << G4endl;

  const G4double velocity = momentumMag / (restMass_c2 + kineticEnergy);
  const G4double labTof = 10.0 * ns, properTof = 0.1 * ns;
  // const G4ThreeVector Spin(1.0, 0.0, 0.0);

  // Driver begins at the start !
  const G4FieldTrack yTrackIn(startPosition, startMomentum.unit(),
                              0.0,  // starting S curve len
                              kineticEnergy, restMass_c2, velocity, labTof, properTof,
                              0  // or &Spin
  );

  double total_step = 0;

  testDriver->OnStartTracking();
  baseIntegrDrv->OnStartTracking();

  G4double spin = 0.0;  // ignore the spin
  G4double magneticMoment = 0.0;  // ignore the magnetic moment
  G4ChargeState chargeState(particleCharge,  // The charge can change (dynamic)
                            spin = 0.0, magneticMoment = 0.0);

  pEquation->SetChargeMomentumMass(chargeState,
                                   momentumMag,  // Momentum in Mev/c
                                   restMass_c2);

  cout << " Error amplification factor = " << errAmplif << endl;

  G4String hlineSep =
    "----------------------------------------------------------------------------------------------"
    "----------------------"
    "----------------------------------------------------------------------------------------------"
    "------\n";
  cout << hlineSep;
  cout << setw(4) << "Step" << " OK? " << setw(10) << "StepLen/mm";  //  << " ";
  for (int i = 0; i < 3; i++)
  {
    if (columns[i])
    {
      if (printSep) cout << " | ";  // Separator
      if (printInp) cout << setw(wd - 6) << "xIn/Try[" << setw(1) << i << "]" << " ";
      if (printInpX) cout << setw(wd - 6) << "xIn/Ref[" << setw(1) << i << "]" << " ";
      cout << setw(wd - 4) << "xOut[" << setw(1) << i << "] ";
      if (printRef)
        cout << setw(wd - 4) << "xOut/Ref[" << setw(1) << i << "]" << " ";  // Reference Solution
      if (printDiff) cout << setw(wdf - 3) << "xDiff[" << setw(1) << i << "] ";
      if (printErr) cout << setw(wd - 1) << "f*X_errRef ";  // [" << setw(1) << i << "] ";

      if (printOk > 2)
      {
        if (printOk > 1)
          cout << setw(wok + 1) << "Ok-P/i";
        else
          cout << setw(3) << "Ok?";
      }
      if (printRatio) cout << " " << setw(wd - 2) << "diff-P/P/eps" << " ";
    }
  }

  if (printOk)
  {
    if (printOk > 1)
      cout << setw(wok + 1) << "Ok-P";
    else
      cout << setw(3) << "Ok-3";
  }
  if (printRatioVec) cout << setw(wd) << " ratioPos ";
  // for (int i = 3; i < 6; i++) { // We will use 0 - 2 instead (x,y,z)
  for (int i = 0; i < 3; i++)
  {
    if (columns[i])
    {
      if (printSep) cout << " | ";  // Separator
      if (printInp) cout << setw(wd - 6) << "pIn/Try[" << setw(1) << i << "]" << " ";
      if (printInpX) cout << setw(wd - 6) << "pIn/Ref[" << setw(1) << i << "]" << " ";
      cout << setw(wd - 4) << "pOut[" << setw(1) << i << "] ";
      if (printRef)
        cout << setw(wd - 4) << "P_Out/Ref[" << setw(1) << i << "]" << " ";  // Reference Solution
      if (printErr) cout << setw(wd - 4) << "f*P_errRef[" << setw(1) << i << "]";
      if (printDiff) cout << setw(wdf - 4) << "pDiff[" << setw(1) << i << "]" << " ";
    }
  }

  cout << endl << hlineSep;
  // ---- End Banner ------------------------------------------------------

  // ----------------Now
  // Stepping--------------------------------------------------------------------
  for (int j = 0; j < no_of_steps; j++)
  {
    bool goodAdvance = true;

    // cerr << " -- j = " << j << " step_len = " << step_len << std::endl;

    referenceStepper->RightHandSide(yInX,
                                    dydxRef);  // compute the value of dydx for the exact stepper
    // Can also add , fieldVal ); to get field value --- given G4double fieldVal[3] );

    G4FieldTrack yTrackOut = yTrackIn;
    G4FieldTrack yTrackOutRef = yTrackIn;

    if (j > 0)  // Let's print the initial points!
    {
      total_step += step_len;

      testDriver->OnComputeStep(&yTrackOut);  // Restart integration ...

      G4double chordDist = 0.25 * total_step;
      G4double hDid = testDriver->AdvanceChordLimited(yTrackOut, total_step, epsTol, chordDist);
      // ---                     *******************
      // 1. Necessary for Interpolation drivers, see comment below  .... 2022.10.20 (Recall/JA)

      if (hDid < total_step)
      {
        G4cerr << " AdvanceChordLimited did not finish interval.  hdid=" << hDid << " vs. hStep"
               << total_step << " - fraction= " << 100.0 * hDid / total_step << " % " << G4endl;
      }

      // 2. For Interpolation drivers, AccurateAdvanced only does Interpolation .... 2022.10.20
      // (Recall/JA) goodAdvance = testDriver->AccurateAdvance(yTrackOut, total_step, epsTol);
      //               *****************************   In/Out

      ThreeVector PositionOut = yTrackOut.GetPosition();
      yout[0] = PositionOut.x();
      yout[1] = PositionOut.y();
      yout[2] = PositionOut.z();
      ThreeVector MomentumOut = yTrackOut.GetMomentum();
      yout[3] = MomentumOut.x();
      yout[4] = MomentumOut.y();
      yout[5] = MomentumOut.z();

      // Compare with a high-quality stepper -- expect it works well for this step size (check!)
      //   This builds on its previous step to create the solution !
      referenceStepper->Stepper(yInX, dydxRef, stepLengthRef, youtX, yerrX);

      // User regular driver instead ...
      baseIntegrDrv->OnComputeStep(
        &yTrackIn);  // Argument yTrackIn is NOT needed for RK drivers ...
      baseIntegrDrv->AccurateAdvance(yTrackOutRef, total_step, epsTol);
    }

    cout << setw(4) << j;  // Printing Step number
    if (goodAdvance)
      cout << " o";  // OK
    else
      cout << " x";  // Bad advance;

    // Check the results
    double sumDiffPos2 = 0.0;  // , sumDiffMom2 = 0.0; // Mag^2 |X-Xref|, |P-Pref|
    double maxDiffPos = 0.0, maxDiffMom = 0.0;
    for (int i = 0; i < 3; i++)
    {
      yDiff[i] = yout[i] / mmTry - youtX[i] / mmRef;
      yAver[i] = 0.5 * (yout[i] / mmTry + youtX[i] / mmRef);

      sumDiffPos2 += (yDiff[i] * yDiff[i]);
      maxDiffPos = std::max(maxDiffPos, fabs(yDiff[i]));
    }

    double magDiffPos = sqrt(sumDiffPos2);
    double relDiffPos = (total_step > 0.0) ? magDiffPos / (total_step / mmTry) : 0.0;
    bool goodPos = relDiffPos < epsTol;
    double ratioPos = relDiffPos / epsTol;

    // Check the momentum
    ThreeVector MomentumNew(yout[3], yout[4], yout[5]);
    MomentumNew *= 1.0 / ppTry;
    ThreeVector MomentumRef(youtX[3], youtX[4], youtX[5]);
    MomentumRef *= 1.0 / ppRef;
    double sumMom2new = 0;
    for (int i = 3; i < 6; i++)
    {
      yDiff[i] = yout[i] / ppTry - youtX[i] / ppRef;
      yAver[i] = 0.5 * (yout[i] / ppTry + youtX[i] / ppRef);

      sumMom2new += yout[i] * yout[i];
      // sumDiffMom2 += (yDiff[i] * yDiff[i]);
      maxDiffMom = std::max(maxDiffMom, fabs(yDiff[i]));
    }
    ThreeVector MomentumDiff(yDiff[3], yDiff[4], yDiff[5]);
    MomentumRef *= 1.0 / ppRef;

    double magDiffMom = MomentumDiff.mag();  // sqrt( sumDiffMom2 );
    double relDiffMom = (momentumMagInit > 0.0) ? magDiffMom / momentumMagInit : 0.0;
    bool goodMom = relDiffMom < epsTol;
    double ratioMom = relDiffMom / epsTol;
    double MomMagShift = sqrt(sumMom2new) / momentumMagInit - 1.0;  // MomentumNew.Mag();

    static ThreeVector totErrMomRefV = ThreeVector(0., 0., 0.);
    ThreeVector errMomRefV = ThreeVector(fabs(yerrX[3]), fabs(yerrX[4]), fabs(yerrX[5]));
    totErrMomRefV += errMomRefV;  // .fabs();

    //-> Then print the data
    cout.setf(ios_base::fixed);

    cout.precision(1);
    cout << setw(9) << total_step / mmTry << " ";

    const int digits = 4;
    cout.precision(digits);

    // Report Position
    for (int i = 0; i < 3; i++)
    {
      if (columns[i])
      {
        cout.precision(digits);
        if (printSep) cout << " | ";  // Separator
        if (printInp) cout << setw(wd - 2) << yIn[i] / mmTry;
        if (printInpX) cout << setw(wd - 2) << yInX[i] / mmRef;
        cout << setw(wd) << yout[i] / mmTry;
        if (printRef) cout << setw(wd) << youtX[i] / mmRef;  // Reference Solution
        if (printDiff)
        {
          cout.precision(epsDigits);
          cout << setw(wdf) << yDiff[i];
        }
        if (printErr) cout << setw(wd) << errAmplif * (yerrX[i] / mmTry);
        bool goodI = fabs(yDiff[i]) < epsTol * fabs(yAver[i]);
        if (printOk > 2)
        {
          if (printOk > 3)
            cout << " " << setw(wok) << (goodI ? "1 ok " : "0 bad");
          else
            cout << " " << setw(3) << (goodI ? "ok " : "bad");
        }
        double ratioPi = (yAver[i] != 0.0) ? yDiff[i] / yAver[i] : 0.0;
        if (printRatio)
          cout  // << "ratio="
            << " " << setw(wd - 2) << ratioPi / epsTol << " ";
      }
    }
    if (printOk)
    {
      if (printOk > 1)
        cout << " " << setw(wok) << (goodPos ? "1 Ok " : "0 Bad");
      else
        cout << " " << setw(3) << (goodPos ? "Ok " : "Bad");
    }
    if (printRatioVec)
    {
      cout << " " << setw(wd) << ratioPos << " ";
    }
    cout.unsetf(ios_base::fixed);

    // Report momentum
    for (int i = 3; i < 6; i++)
    {
      if (columns[i])
      {
        if (printSep) cout << " | ";  // Separator
        if (printInp) cout << setw(wd - 1) << yIn[i] / ppTry << " ";
        if (printInpX) cout << setw(wd - 1) << yInX[i] / ppRef << " ";
        cout << setw(wd) << yout[i] / ppTry;
        if (printRef) cout << setw(wd) << youtX[i] / ppRef;  // Reference Solution
        cout.setf(ios_base::scientific);
        if (printDiff)
        {
          cout.precision(epsDigits);
          cout << setw(wdf) << yDiff[i];
        }
        if (printErr)
        {
          cout.precision(epsDigits);
          cout << " " << setw(wd) << errAmplif * yerrX[i] / ppTry;
        }
        cout.unsetf(ios_base::scientific);
        if (printOk || printRatio)
        {
          bool goodI = fabs(yDiff[i]) < epsTol * fabs(yAver[i]);
          // if( printOk>1 ) cout<< " " << setw(3)<< (goodI ? "ok " : "bad" );
          if (printOk > 2)
          {
            if (printOk > 3)
              cout << " " << setw(wok) << (goodI ? "1 ok " : "0 bad");
            else
              cout << " " << setw(3) << (goodI ? "ok " : "bad");
          }
          double ratioMi = (yAver[i] != 0.0) ? yDiff[i] / yAver[i] : 0.0;
          if (printRatio)
          {
            // cout.unsetf (ios_base::scientific);
            cout.setf(ios_base::fixed);
            cout << " " << setw(wd) << ratioMi / epsTol << " ";
            cout.unsetf(ios_base::fixed);
            // cout.setf (ios_base::scientific);
          }
        }
      }
    }

    if (printOk)
    {
      // cout<< " Vec: ";
      cout << " ";
      if (printOk > 1) cout << goodMom;  // << setw(wok)<< (goodMom ? "Ok " : "Bad" );
      cout << " " << setw(3) << (goodMom ? "Ok " : "Bad");
    }
    if (printRatioVec)
    {
      cout.setf(ios_base::fixed);
      cout << " " << setw(wd) << ratioMom << " ";
      cout.unsetf(ios_base::fixed);
    }

    cout << setw(wd) << MomMagShift / epsTol;

    for (int i = 0; i < 6; i++)  // Print auxiliary components
    {
      double unitTry = 1;
      double unitRef = 1;

      if (i >= 3)
      {
        unitTry = ppTry / mmTry;  //  dp / ds
        unitRef = ppRef / mmRef;  //  dp / ds
      }

      if (i == 0)
      {
        cout.setf(ios_base::fixed);
      }
      else if (i == 3)
      {
        cout.unsetf(ios_base::fixed);
        cout.setf(ios_base::scientific);
      }

      if (columns[i + 6])
      {
        // cout << " dy/dx["<<i<<"] = ";
        if (printSep) cout << " | ";  // Separator
        cout << setw(wd) << dydxRef[i] / unitRef;  // Reference Solution
        cout.precision(epsDigits);
        if (0)
        {  //  ( printDiff ){  -- Not meaningful
          cout << setw(wdf) << (dydx[i] / unitTry) - (dydxRef[i] / unitRef);
        }
      }
    }
    cout.unsetf(ios_base::scientific);
    if (j > 0)  // Step 0 did not move -- printed the starting values
    {
      cout.setf(ios_base::fixed);
      cout.precision(2);
      cout << setw(wd) << atan2(yout[1], yout[0]) / degree;
      if (printRef) cout << setw(wd) << atan2(youtX[1], youtX[0]) / degree;

      // Copy yout into yIn
      for (int i = 0; i < 6; i++)
      {
        // yIn[i] = yout[i];
        yInX[i] = youtX[i];
      }
    }

    cout << "\n";
  }

  /*-----------------END-STEPPING------------------*/
  delete uEquation;
  delete pEquation;

  delete testDriver;
  delete baseIntegrDrv;

  // Old stuff
  // myStepper->InformDone();
  delete dopri5_stepper;
  delete referenceStepper;
  // delete referenceDriver;
  // delete tMagneticField;

  cout << "\n\n#-------------End of output-----------------\n";
}
