//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//

#include "G4GenericMessenger.hh"

#include "G4String.hh"
#include "G4ThreeVector.hh"
#include "G4Types.hh"
#include "G4UImanager.hh"

#include <gtest/gtest.h>

class testGenericMessenger
{
  public:

    testGenericMessenger();

  public:

    bool b;
    int i;
    unsigned int ui;
    long l;
    unsigned long ul;
    float f;
    double d;
    double d_u;
    G4complex c;
    G4String s;
    G4String s2;
    G4ThreeVector v3;
    G4ThreeVector v3_u;
    bool m0_called;
    bool m1_called;
    int m1_arg;
    bool m2_called;
    G4String m2_arg;
    bool m3_called;
    std::string m3_arg0;
    int m3_arg1;
    void method0() { m0_called = true; }
    void method1(double a)
    {
      m1_called = true;
      m1_arg = a;
    }
    void method2(const G4String& sIn)
    {
      m2_called = true;
      m2_arg = sIn;
    }
    int method3(const std::string& sIn, int n)
    {
      m3_called = true;
      m3_arg0 = sIn;
      m3_arg1 = n;
      return n;
    }

    void setflag(bool param) { b = param; }

    void setvalueandflag(int value, bool flag)
    {
      i = value;
      b = flag;
    }

  private:

    G4GenericMessenger messenger;
};

testGenericMessenger::testGenericMessenger()
  : b(false),
    i(0),
    ui(0),
    l(0),
    ul(0),
    f(0.0),
    d(0.0),
    d_u(0.0),
    m0_called(false),
    m1_called(false),
    m1_arg(0),
    m2_called(false),
    m3_called(false),
    m3_arg1(0),
    messenger(this, "/mytest/")
{
  messenger.DeclareProperty("bool", b, "bool property");
  messenger.DeclareMethod("flag", &testGenericMessenger::setflag, "set boolean with method");
  messenger.DeclareMethod("valueandflag", &testGenericMessenger::setvalueandflag,
                          "set value and boolean with method");

  messenger.DeclareProperty("int", i, "int property");
  messenger.DeclareProperty("uint", ui, "unsigned int property").SetRange("value>=0");
  messenger.DeclareProperty("long", l, "long property");
  messenger.DeclareProperty("ulong", ul, "unsigned long property");
  messenger.DeclareProperty("float", f, "float property");
  messenger.DeclareProperty("double", d, "double property");
  messenger.DeclareProperty("complex", c, "complex property");
  messenger.DeclareProperty("string", s, "string property");
  messenger.DeclareProperty("string2", s2, "string property with other options")
    .SetParameterName("String", true)
    .SetCandidates("string1 string2 string3")
    .SetDefaultValue("string2");
  messenger.DeclareProperty("3vector", v3, "G4ThreeVector property");

  messenger.DeclareProperty("doubleunit", d_u)
    .SetGuidance("double energy property")
    .SetGuidance("additional guidance")
    .SetDefaultUnit("MeV")
    .SetStates(G4State_PreInit, G4State_Idle);
  messenger.DeclareProperty("3vectorunit", v3_u, "3 vector length property")
    .SetUnitCategory("Length");

  messenger.DeclareMethod("method0", &testGenericMessenger::method0, "Method with 0 arguments");
  messenger.DeclareMethod("method1", &testGenericMessenger::method1, "Method with 1 arguments");
  messenger.DeclareMethod("method1unit", &testGenericMessenger::method1, "Method with 1 arguments")
    .SetUnit("mm");
  messenger.DeclareMethod("method2", &testGenericMessenger::method2,
                          "Method with 1 const ref argument");
  messenger.DeclareMethod("method3", &testGenericMessenger::method3, "Method with 2 arguments");
}

TEST(G4GenericMessenger, BooleanParameter)
{
  testGenericMessenger tst;
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  UImanager->ApplyCommand("/mytest/bool 1");
  EXPECT_TRUE(tst.b);

  UImanager->ApplyCommand("/mytest/bool 0");
  EXPECT_FALSE(tst.b);

  UImanager->ApplyCommand("/mytest/bool F");
  EXPECT_FALSE(tst.b);

  UImanager->ApplyCommand("/mytest/bool T");
  EXPECT_TRUE(tst.b);

  UImanager->ApplyCommand("/mytest/bool false");
  EXPECT_FALSE(tst.b);

  UImanager->ApplyCommand("/mytest/bool true");
  EXPECT_TRUE(tst.b);

  EXPECT_EQ(UImanager->GetCurrentIntValue("/mytest/bool"), 1);

  // Invalid parameter can only be tested for fact that it doesn't alter state
  UImanager->ApplyCommand("/mytest/bool ABC");
  EXPECT_TRUE(tst.b);
}

TEST(G4GenericMessenger, Bugzilla2606)
{
  testGenericMessenger tst;
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  // Bugzilla 2606 regression tests
  // 1. FuncRef1 (single boolean function call)
  UImanager->ApplyCommand("/mytest/flag 1");
  EXPECT_TRUE(tst.b);
  UImanager->ApplyCommand("/mytest/flag 0");
  EXPECT_FALSE(tst.b);
  UImanager->ApplyCommand("/mytest/flag true");
  EXPECT_TRUE(tst.b);
  UImanager->ApplyCommand("/mytest/flag false");
  EXPECT_FALSE(tst.b);
  UImanager->ApplyCommand("/mytest/flag T");
  EXPECT_TRUE(tst.b);
  UImanager->ApplyCommand("/mytest/flag F");
  EXPECT_FALSE(tst.b);
  UImanager->ApplyCommand("/mytest/flag y");
  EXPECT_TRUE(tst.b);
  UImanager->ApplyCommand("/mytest/flag n");
  EXPECT_FALSE(tst.b);
  UImanager->ApplyCommand("/mytest/flag YeS");
  EXPECT_TRUE(tst.b);
  UImanager->ApplyCommand("/mytest/flag nO");
  EXPECT_FALSE(tst.b);

  // 2. FuncRef2 (function with two parameters, either/both can be bool)
  // Reset to known values
  tst.i = 0;
  tst.b = false;

  UImanager->ApplyCommand("/mytest/valueandflag 42 yes");
  EXPECT_EQ(tst.i, 42);
  EXPECT_TRUE(tst.b);
  UImanager->ApplyCommand("/mytest/valueandflag 43 T");
  EXPECT_EQ(tst.i, 43);
  EXPECT_TRUE(tst.b);
  UImanager->ApplyCommand("/mytest/valueandflag 44 1");
  EXPECT_EQ(tst.i, 44);
  EXPECT_TRUE(tst.b);
  UImanager->ApplyCommand("/mytest/valueandflag 52 no");
  EXPECT_EQ(tst.i, 52);
  EXPECT_FALSE(tst.b);
  UImanager->ApplyCommand("/mytest/valueandflag 53 F");
  EXPECT_EQ(tst.i, 53);
  EXPECT_FALSE(tst.b);
  UImanager->ApplyCommand("/mytest/valueandflag 54 0");
  EXPECT_EQ(tst.i, 54);
  EXPECT_FALSE(tst.b);

  // ... and more
}

TEST(G4GenericMessenger, IntegralParameter)
{
  testGenericMessenger tst;
  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  EXPECT_EQ(tst.i, 0);
  EXPECT_EQ(tst.ui, 0);
  EXPECT_EQ(tst.l, 0);
  EXPECT_EQ(tst.ul, 0);

  // Signed
  UImanager->ApplyCommand("/mytest/int 99");
  EXPECT_EQ(tst.i, 99);

  UImanager->ApplyCommand("/mytest/int -99");
  EXPECT_EQ(tst.i, -99);
  EXPECT_EQ(UImanager->GetCurrentIntValue("/mytest/int"), -99);

  UImanager->ApplyCommand("/mytest/int 12334556694584792872726727627127");
  EXPECT_EQ(tst.i, -99);

  // Unsigned
  UImanager->ApplyCommand("/mytest/uint 88");
  EXPECT_EQ(tst.ui, 88);
  EXPECT_EQ(UImanager->GetCurrentIntValue("/mytest/uint"), 88);

  // Long
  UImanager->ApplyCommand("/mytest/long 999");
  EXPECT_EQ(tst.l, 999L);

  // Unsigned Long
  UImanager->ApplyCommand("/mytest/ulong 888");
  EXPECT_EQ(tst.ul, 888L);
}

TEST(G4GenericMessenger, FloatingPointParameter)
{
  testGenericMessenger tst;
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  // Float
  EXPECT_EQ(tst.f, 0.0f);
  UImanager->ApplyCommand("/mytest/float 9.9");
  EXPECT_EQ(tst.f, 9.9f);
  UImanager->ApplyCommand("/mytest/float -9.9E+10");
  EXPECT_EQ(tst.f, -9.9E+10f);

  // Double
  EXPECT_EQ(tst.d, 0.0);
  UImanager->ApplyCommand("/mytest/double 99.99");
  EXPECT_EQ(tst.d, 99.99);
  UImanager->ApplyCommand("/mytest/double -99.99E+10");
  EXPECT_EQ(tst.d, -99.99E+10);
  EXPECT_EQ(UImanager->GetCurrentDoubleValue("/mytest/double"), -99.99E+10);

  // Double plus unit
  EXPECT_EQ(tst.d_u, 0.0);
  UImanager->ApplyCommand("/mytest/doubleunit 77.8 TeV");
  EXPECT_EQ(tst.d_u, 77800000.0);
  EXPECT_EQ(UImanager->GetCurrentDoubleValue("/mytest/doubleunit"), 77800000.0);
}

TEST(G4GenericMessenger, ComplexParameter)
{
  testGenericMessenger tst;
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  EXPECT_EQ(tst.c, G4complex(0.0, 0.0));
  UImanager->ApplyCommand("/mytest/complex (8.8,9.9)");
  EXPECT_EQ(tst.c, G4complex(8.8, 9.9));
}

TEST(G4GenericMessenger, StringParameter)
{
  testGenericMessenger tst;
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  EXPECT_EQ(tst.s, G4String());
  UImanager->ApplyCommand("/mytest/string mystring");
  EXPECT_EQ(tst.s, "mystring");
  UImanager->ApplyCommand("/mytest/string \"another string\"");
  EXPECT_EQ(tst.s, "another string");

  EXPECT_EQ(tst.s2, G4String());
  UImanager->ApplyCommand("/mytest/string2");
  EXPECT_EQ(tst.s2, "string2");
  UImanager->ApplyCommand("/mytest/string2 string3");
  EXPECT_EQ(tst.s2, "string3");
}

TEST(G4GenericMessenger, ThreeVectorParameter)
{
  testGenericMessenger tst;
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  EXPECT_EQ(tst.v3, G4ThreeVector());
  UImanager->ApplyCommand("/mytest/3vector 1. 2. 3.");
  EXPECT_EQ(tst.v3, G4ThreeVector(1.0, 2.0, 3.0));

  EXPECT_EQ(tst.v3_u, G4ThreeVector());
  UImanager->ApplyCommand("/mytest/3vectorunit 4.1 5.1 6.1 cm");
  EXPECT_EQ(tst.v3_u, G4ThreeVector(41., 51., 61));
}

TEST(G4GenericMessenger, FunctionCall)
{
  testGenericMessenger tst;
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  //---Method with 0 arguments----------------------------------------
  EXPECT_EQ(tst.m0_called, false);
  UImanager->ApplyCommand("/mytest/method0");
  EXPECT_EQ(tst.m0_called, true);

  //---Method with 1 arguments----------------------------------------
  EXPECT_EQ(tst.m1_called, false);
  UImanager->ApplyCommand("/mytest/method1 999.");
  EXPECT_EQ(tst.m1_called, true);
  EXPECT_EQ(tst.m1_arg, 999.);

  UImanager->ApplyCommand("/mytest/method1unit 9.0 km");
  EXPECT_EQ(tst.m1_arg, 9.0E+6);

  //---Method with 1 const ref arguments----------------------------------------
  EXPECT_EQ(tst.m2_called, false);
  UImanager->ApplyCommand("/mytest/method2 abcdef");
  EXPECT_EQ(tst.m2_called, true);
  EXPECT_EQ(tst.m2_arg, "abcdef");

  //---Method with 2 arguments----------------------------------------
  EXPECT_EQ(tst.m3_called, false);
  UImanager->ApplyCommand("/mytest/method3 abcdef 999");
  EXPECT_EQ(tst.m3_called, true);
  EXPECT_EQ(tst.m3_arg0, "abcdef");
  EXPECT_EQ(tst.m3_arg1, 999);
}
