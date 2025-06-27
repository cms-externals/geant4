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
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>

#include "Rtypes.h"
#include "TFile.h"
#include "TH1.h"

using namespace std;

string ComposeG4Tag( string );
string ComposeFileNameTest47( string, string, string, string );
string ComposeFileNameGeneric( string, string, string, string ); 

void SQLCompose( int tdid,
                 string beam, string target, string energy, string momentum, string secondary,
	         string model, string g4version,
	         string histoname,
		 string reaction_details,
	         string xlabel, string ylabel, // x-y axis label(s)
	         string experiment,
	         string observable,
	         string status,
	         int mask ) // mask=0 -> energy; mask=1 -> momentum
{

   // this location is on tev.fnal.gov (Wilson cluster) 
   // where several group of tests usually run
   //
   string location = "/g4/g4p/pbs/g4-had-validation/regression-test-files/";
   
   string tname = "test";
   ostringstream os;
   os << tdid;
   tname += os.str();
   
   // need to reformat g4version --> g4tag: geant4-XX.Y-<p/ref>ZZ --> geant4-XX-YY-<patch/ref>-ZZ
   //  
   string g4tag = ComposeG4Tag( g4version );
            
   string fname = location + tname + "/" + g4tag;
      
   if ( !experiment.empty() )
   {
      for( unsigned int i=0; i<experiment.size(); ++i ) experiment[i] = tolower(experiment[i]);
      fname += "/" + experiment + "-histo";
   }
   
   if ( tdid == 47 )
   {
      fname += ComposeFileNameTest47( beam, target, momentum, model );
   }
   else
   {
      fname += ComposeFileNameGeneric( beam, target, momentum, model );
   }
   fname += ".root";
   
   TFile* histofile = TFile::Open( fname.c_str() );
   
   ofstream sqlfile;
   string fsqlname = tname;
//
//   if ( !experiment.empty() )
//   {
//      fsqlname += ("-" + experiment);
//   }
//
   fsqlname += "-histo.sql";
   sqlfile.open( fsqlname.c_str(), ios::ate|ios::app ); 
   
   TH1F* histo = (TH1F*)histofile->Get( histoname.c_str() );
   
   if ( !histo )
   {
      cout << " NULL histogram found by name " << histoname << " in Root file " << fname << endl;
      return;
   }
   
   string bm = beam;
   if ( beam == "piplus" )
   {
      bm = "pi+";
   }
   else if ( beam == "piminus" )
   {
      bm = "pi-";
   }
   else if ( beam == "muminus" )
   {
      bm = "mu-";
   }
   
   string sec = secondary;
   if ( secondary == "piplus" )
   {
      sec = "pi+";
   }
   else if ( secondary == "piminus" )
   {
      sec = "pi-";
   }
   
   sqlfile <<"INSERT INTO histogramvb (tdid,histoname,g4version,status,model,target,reaction,beam,";     
   sqlfile <<"beammomentum,beamenergy,secondary,observable,xdes,ydes,nbins,minx,maxx,binvalues,bincenter,binwidth,errorup,errorlow,mask)";   
   sqlfile << " VALUES (";
   sqlfile << tdid <<",";
   sqlfile << "'"<< tname <<"',";
   sqlfile << "'"<< g4version <<"',";
   sqlfile << "'"<< status <<"',";
   sqlfile << "'"<< model <<"',";
   sqlfile << "'"<< target <<"',";
//   sqlfile << "'"<< histo->GetTitle() <<"',";
   sqlfile << "'"<< "beam on target -> secondary +X" << reaction_details <<"',";
// -->-->   sqlfile << "'"<< beam <<"',";
   sqlfile << "'"<< bm <<"',";
   sqlfile << "'"<< momentum <<"',"; // beammomentum in the DB table
   sqlfile << "'"<< energy <<"',";   // beamenergy in the DB table
   sqlfile << "'"<< sec <<"',";
   sqlfile << "'"<< observable <<"',";
   sqlfile << "'"<< xlabel<<"',";  // xdes in the DB table
   sqlfile << "'"<< ylabel<<"',";  // ydes in the DB table
   sqlfile << histo->GetNbinsX()<<",";
   sqlfile << histo->GetXaxis()->GetXmin()<<",";
   sqlfile << histo->GetXaxis()->GetXmax()<<",";
   sqlfile << "'{";
   sqlfile <<endl;
  
   int nx = histo->GetNbinsX();
   
   for ( int k=1; k<nx; ++k ) 
   {
      sqlfile <<  histo->GetBinContent(k) << "," << endl;	
   }
   sqlfile << histo->GetBinContent(nx) << "}'," << endl;
   
   sqlfile << " '{";
   for ( int k=1; k<nx; ++k ) 
   {
      sqlfile <<  histo->GetBinCenter(k) << "," << endl;	
   }
   sqlfile << histo->GetBinCenter(nx) << "}'," << endl;
  
   sqlfile << " '{";
   for ( int k=1; k<nx; ++k ) 
   {
      sqlfile <<  histo->GetBinWidth(k) << "," << endl;	
   }
   sqlfile << histo->GetBinWidth(nx) << "}'," << endl;

   sqlfile << " '{";
   for ( int k=1; k<nx; ++k ) 
   {
      sqlfile << histo->GetBinErrorUp(k) << "," << endl;	
   }
   sqlfile << histo->GetBinErrorUp(nx) << "}'," << endl;
  
   sqlfile << " '{";
   for ( int k=1; k<nx; ++k ) 
   {
      sqlfile << histo->GetBinErrorLow(k) << "," << endl;	
   }
   sqlfile << histo->GetBinErrorLow(nx) << "}'," << endl;
   sqlfile << mask << ");" << endl;
   
   histofile->Close();
 
   return;

}

string ComposeG4Tag( string g4version )
{

   string g4tag = "geant4-";
   
   size_t itr1 = g4version.find_first_of("-");
   
   size_t itr2 = g4version.find_first_of(".");
   
   if ( (itr2-itr1) < 3 ) g4tag += "0";     
   
   for ( size_t i=itr1+1; i<itr2; ++i ) g4tag += g4version[i];
   g4tag += "-";
      
   size_t itr3 = g4version.find_first_of("-",itr2);

   if ( itr3 == string::npos ) // no patch or ref
   {
      if ( (g4version.length()-itr2) < 3 ) g4tag += "0";
      for ( size_t i=itr2+1; i<g4version.length(); ++i ) g4tag += g4version[i];
      return g4tag;
   }

   if ( (itr3-itr2) < 3 ) g4tag += "0";
   for ( size_t i=itr2+1; i<itr3; ++i ) g4tag += g4version[i];
   
   size_t reltype = g4version.find("ref",itr3);
   
   if ( reltype != string::npos )
   {
      // reference release
      //
      g4tag += "-ref-";
      for ( size_t i=reltype+3; i<g4version.length(); ++i ) g4tag += g4version[i];
      return g4tag;
   }

   reltype = g4version.find("p",itr3);
   
   if ( reltype != string::npos )
   {
      // patch release
      //
      g4tag += "-patch-";
      for ( size_t i=reltype+1; i<g4version.length(); ++i ) g4tag += g4version[i];
      return g4tag;
   }
   
   return g4tag;

}

string ComposeFileNameGeneric( string beam, string target, 
                               string momentum, 
			       string model )
{

   string fname = "";

   if ( beam == "gamma" )
   {
      fname += ( "/" + beam + momentum + target + model );
   }
   else
   {
      fname += ( "/" + beam + target + momentum + model );
   }

   return fname;

}

string ComposeFileNameTest47( string beam, string target, 
                              string momentum, 
			      string model )
{

   string fname = "";

   fname += ( "/" + beam + target + model + momentum );

   return fname;

}

