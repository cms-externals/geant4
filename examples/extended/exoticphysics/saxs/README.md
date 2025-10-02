\page Examplesaxs Example saxs

\author Gianfranco Patern&ograve; - INFN and University Ferrara (Italy) \n
 paterno@fe.infn.it

The example saxs implements the typical setup of a Small Angle X-ray 
Scattering (SAXS) experiment. It is meant to illustrate the 
usage of molecular interference (MI) of Rayleigh (coherent) scattering
of photons inside the matter, which is implemented in the
G4PenelopeRayleighModelMI model.


## Geometry

The setup consists of a phantom/sample under investigation, slits
to collimate the photon beam and a shielded detector to collect
the photons scattered by the phantom (see SAXSDetectorConstruction).  
   
The geometry is scalable through the interactive commands defined 
in the SAXSDetectorConstructionMessenger class. All the significant
quantities, such as the setup (scattering) rotation angle, the position 
and size of all the volumes, as well as the phantom material can be
set via macro commands. 
   
Two macro files come with this example: saxs.in and saxs_slits.in.

-#  In the saxs.in macro, the phantom is a cylinder with a diameter and
a height of 10 mm made of a mixture of 80% fat and 20% water. 
In general, if the argument of <code>/det/setPhantomMaterial</code> command is 2, 
as in this case, the material is a biological tissue ("MedMat") 
defined as a mixture of fat, water, collagen and hydroxyapatite. 
The weight fraction of the mixture components can be set through commands 
```
 /det/setComp0
 /det/setComp1
 /det/setComp2
 /det/setComp3
```
respectively. 
The tissue form factor (including MI) is automatically calculated as a 
weighed sum of the form factors of the basis components.     
In this case, no slits are foreseen and the sensitive detector 
positioned 400 mm downstream of the phantom collects all the photons 
transmitted and scattered by the phantom, which is irradiated 
by a pencil beam with an energy of 20 keV.
   
-# In the saxs_slits.in macro, the phantom is again a cylinder with a
diameter and a height of 10 mm. The phantom is made of a custom
material ("CustomMat") whose density and composition is set through
```
 /det/setCustomMatDensity
 /det/setCustomMatHmassfract,
 /det/setCustomMatNmassfract
 /det/setCustomMatOmassfract
```
commands. In general, a custom material can be defined by 
specifying the mass fraction of H, C, N, O, Na, P, S, Cl, K, and Ca via
commands analogous to those mentioned above. In this case, the material 
composition corresponds to that of ammonium nitrate
(NH<sub>4</sub>NO<sub>3</sub>).
For a custom material, the user can provide the path of the file with
the material form factor (with MI) through the
```
 /det/SetCustomMatFF
```
command. As an example, the file myFF.dat contains the form factor of 
NH<sub>4</sub>NO<sub>3</sub> measured by Harding in 1999. 
In this case the slits upstream and downstream the phantom 
are present. This setup is suitable for both monochromatic and 
polychromatic beams. To speed-up the simulation, a monochromatic 
photon beam was chosen, but a polychromatic beam can be easily defined.    

## Physics

In this example, only electromagnetic processes and decays are considered.

They are defined in a custom physics list that allows the user to
choose among various EM PhysicsList constructors. In particular,
by choosing G4EmPenelopePhysicsMI and setting fUseMIFlag as true, 
it is possible to enable the molecular interference effects. This
is the default configuration.

## Action Initialization

SAXSActionInitialization class instantiates and registers to 
Geant4 kernel all user action classes. While in sequential mode 
the action classes are instatiated just once, by invoking the 
method: SAXSActionInitialization::Build(),
in multi-threading mode the same method is invoked for each thread 
worker and so all user action classes are defined thread-local.

A run action class is instantiated both thread-local and global.
That's why its instance is created also in the method
SAXSActionInitialization::BuildForMaster(), which is 
invoked only in multi-threading mode.

## Primary generator

The primary generator action class employs the G4GeneralParticleSource (GPS)
generator. The primary beam has to be defined via the G4 built-in 
commands of the G4GeneralParticleSource in a input macro file. 
In particular, a photon beam directed toward the phantom must be defined
to test the MI effects. The X-ray beam can be monochromatic or
polychromatic, parallel or divergent.     


## Event and Detector Response

An event consists of the generation of a single particle which is 
transported through the phantom and then to the sensitive detector. 
   
The interactions of the photons inside the phantom, and in particular,
the scattering events, are scored in a dedicated ntuple through the
SAXSSteppingAction class.
   
The hits of the particles on the sensitive detector positioned 
downstream of the phantom (SAXSSensitiveDetectorHit) are recorded
in a dedicated ntuple through the SAXSSensitiveDetector class. 

## Analysis

The analysis tools are used to accumulate statistics.
ntuple are created in SAXSRunAction::SAXSRunAction() 
constructor for the following quantities:
  
Ntuple1 (<code>part</code>) - Particles impinging on the Sensitive Detector (SD):

- energy of the particles
- position of the hits
- momentum of the particles
- time of the hits
- type of impinging particles
- ID number of the impinging particles 
- number of scattering events a primary had before hitting the SD
- event number of the hits
   
Ntuple2 (<code>scatt</code>) - Interactions of photons inside the phantom:

- ID of the process occurred 
      (0-> transportation, 1->Rayleigh, 2->Compton, 3->Photoelectic)
- initial energy of the particles
- scattering angle

The ntuples are saved in the output file in the Root format.
When running in multi-threading mode, the ntuples accumulated 
on threads are automatically merged in a single output file.
   
The default output format is root. Two root scripts come with 
this example to analyze the output file: scattAnalysis.C and 
ADXRD.C. The first can be used to analyze the <code>scatt</code> ntuple, while 
the second can be used for <code>part</code> ntuple.

## How to run

Execute <code>saxs</code> in the 'interactive mode' with visualization:
```
  % ./saxs
```
and type in the commands line by line:
```
  Idle> /control/verbose 2
  Idle> /tracking/verbose 1
  Idle> ...
  Idle> /run/beamOn 10 
  Idle> ...
  Idle> exit
```
or it is possible to run a macro file (test.in is a simple macro where the 
primary beam is defined through the usual GPS commands):
```
  Idle> /control/execute test.in
  Idle> /run/beamOn 10 
  ....
  Idle> exit
```

Execute saxs in the 'batch' mode from macro files (without visualization)
```
  % ./saxs saxs.in [Ncores]   
  % ./saxs saxs_slits.in [Ncores]
```
  <code>Ncores</code> (optional argument) is the number of threads the user wants
  to use in MT mode.
      
## Appendix

The following paragraphs are common to all basic examples

### Visualization

The visualization manager is set via the G4VisExecutive class
in the main() function in saxs.cc.    
The initialisation of the drawing is done via a set of /vis/ commands
in the macro vis.mac. This macro is automatically read from
the main function when the example is used in interactive running mode.

By default, vis.mac opens an OpenGL viewer (<code>/vis/open OGL</code>).
The user can change the initial viewer by commenting out this line
and instead uncommenting one of the other <code>/vis/open</code>
statements, such as
HepRepFile or DAWNFILE (which produce files that can be viewed with the
HepRApp and DAWN viewers, respectively).  Note that one can always
open new viewers at any time from the command line.  For example, if
you already have a view in, say, an OpenGL window with a name
"viewer-0", then
```
  /vis/open DAWNFILE
```
then to get the same view
```
  /vis/viewer/copyView viewer-0
```
or to get the same view *plus* scene-modifications
```
  /vis/viewer/set/all viewer-0
```
then to see the result
```
  /vis/viewer/flush
```

The DAWNFILE, HepRepFile drivers are always available
(since they require no external libraries), but the OGL driver requires
that the Geant4 libraries have been built with the OpenGL option.

From Release 9.6 the vis.mac macro in example B1 has additional commands
that demonstrate additional functionality of the vis system, such as
displaying text, axes, scales, date, logo and shows how to change
viewpoint and style.  Consider copying these to other examples or
your application.  To see even more commands use help or
ls or browse the available UI commands in the Application
Developers Guide, Section 7.1.

Since 11.1, the TSG visualization driver can also produce the "offscrean"
file output in png, jpeg, gl2ps formats without drawing on the screen.
It can be controlled via UI commands provided in '/vis/tsg' which are
demonstrated in the tsg_offscreen.mac macro in example B5.

For more information on visualization, including information on how to
install and run DAWN, OpenGL and HepRApp, see the visualization tutorials,
for example,\n
- <a href="http://geant4.slac.stanford.edu/Presentations/vis/G4OpenGLTutorial/G4OpenGLTutorial.html">
  OpenGL Tutorial </a>
- <a href="http://geant4.slac.stanford.edu/Presentations/vis/G4DAWNTutorial/G4DAWNTutorial.html">
  DAWN Tutorial </a>
- <a href="http://geant4.slac.stanford.edu/Presentations/vis/G4HepRAppTutorial/G4HepRAppTutorial.html">
  HepRApp Tutorial </a>

The tracks are automatically drawn at the end of each event, accumulated
for all events and erased at the beginning of the next run.

### User Interfaces
 
The user command interface is set via the G4UIExecutive class
in the main() function in saxs.cc.

The selection of the user command interface is then done automatically
according to the Geant4 configuration or it can be done explicitly via
the third argument of the G4UIExecutive constructor (see exampleB4a.cc).
