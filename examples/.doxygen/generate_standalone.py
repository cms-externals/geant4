import os
import shutil
import re
import subprocess

# Script to generate standalone Doxygen documentation
# for the examples with non unique class names.

CURDIR = os.getcwd()

DOXYFILE_PATH = '.doxygen/doc'

BACK_PATH2 = "../../../.doxygen/doc"
BACK_PATH3 = "../../../../.doxygen/doc"
BACK_PATH4 = "../../../../../.doxygen/doc"

# Placeholder for ADD_SHARED and ADD_COMMON - if these are meant to be dynamic,
# you'll need to define how they get their values.
ADD_SHARED = ""
ADD_COMMON = ""

## --- NEW: Flag to track if geant4.tag warning has been issued ---
#_geant4_tag_warning_issued = False
# Use a list to hold the flag. This avoids the 'global' keyword issue. ---
# The boolean is now the first (and only) element of this list.
_geant4_tag_warning_state = [False]

def generate(example_category, example_dir):
    """
    Generates Doxygen documentation for a given example.

    Args:
        example_category (str): basic or extended
        example_dir (str): The relative path from the category directory to the example directory.
        level (int): The level of example directory wrt category
    """
    global _geant4_tag_warning_issued # Declare intent to modify the global variable

    print(f"processing {example_dir}")
    original_cwd = os.getcwd()
    try:
        os.chdir(os.path.join('../' + example_category, example_dir))
        example_name = os.path.basename(os.getcwd())
        # print(f"Generating Doxyfile for example: {example_name}")

        # Evaluate relative path to Doxygen template file
        level = example_dir.count('/') + 2
        relative_path_to_doc = ''
        for i in range(0,level):
            relative_path_to_doc += '../'
        relative_path_to_doc += DOXYFILE_PATH

        doxyfile_template_path = os.path.join(CURDIR, "Doxyfile_standalone")
        with open(doxyfile_template_path, 'r') as f:
            doxyfile_content = f.read()

        # Perform substitutions using re.sub for clarity and consistency
        doxyfile_content = re.sub(r"EXAMPLE_NAME", example_name, doxyfile_content)
        doxyfile_content = re.sub(r"BACK_PATH", relative_path_to_doc, doxyfile_content)
        doxyfile_content = re.sub(r"ADD_SHARED", ADD_SHARED, doxyfile_content)
        doxyfile_content = re.sub(r"ADD_COMMON", ADD_COMMON, doxyfile_content)

        with open("Doxyfile", 'w') as f:
            f.write(doxyfile_content)

        # Create a symbolic link
        geant4_tag_src = os.path.join(CURDIR, "geant4.tag")
        geant4_tag_dest = "geant4.tag"
        if os.path.exists(geant4_tag_src):
            if os.path.lexists(geant4_tag_dest): # Check if the symlink already exists
                os.unlink(geant4_tag_dest) # Remove existing symlink
            os.symlink(geant4_tag_src, geant4_tag_dest)
        else:
            # --- MODIFIED: Print warning only once ---
            if not _geant4_tag_warning_state[0]: # Check the first element of the list
                print(f"Warning: geant4.tag not found at {geant4_tag_src}")
                _geant4_tag_warning_state[0] = True # Set the flag to True in the list
            # --- END MODIFIED ---

        # Run doxygen
        output_file_path = os.path.join(CURDIR, f"doxygen_{example_name}.out")
        print(output_file_path)
        with open(output_file_path, 'w') as outfile:
            subprocess.run(["doxygen"], stdout=outfile, stderr=outfile, check=False)

    except FileNotFoundError as e:
        print(f"Error: A required file or directory was not found: {e}")
    except Exception as e:
        print(f"An error occurred during processing {example_dir}: {e}")
    finally:
        # Clean up
        if os.path.exists("Doxyfile"):
            os.remove("Doxyfile")
        if os.path.lexists("geant4.tag"):
            os.remove("geant4.tag")
        os.chdir(original_cwd)

# ... (rest of your script remains the same) ...

# Process examples in second level directory in extended
basic_examples = [
    "B1", "B2", "B3", "B4", "B5"
]

extended_examples = [
    "analysis/AnaEx01", "analysis/AnaEx02", "analysis/AnaEx03", "analysis/B1Con",
    #
    "biasing/B01", "biasing/B02", "biasing/B03", "biasing/GB01", "biasing/GB02",
    "biasing/GB03", "biasing/GB04", "biasing/GB05", "biasing/GB06", "biasing/GB07",
    "biasing/ReverseMC01",
    #
    "common/exCommon",
    #
    "electromagnetic/TestEm0", "electromagnetic/TestEm1", "electromagnetic/TestEm2",
    "electromagnetic/TestEm3", "electromagnetic/TestEm4", "electromagnetic/TestEm5",
    "electromagnetic/TestEm6", "electromagnetic/TestEm7", "electromagnetic/TestEm8",
    "electromagnetic/TestEm9", "electromagnetic/TestEm10", "electromagnetic/TestEm11",
    "electromagnetic/TestEm12", "electromagnetic/TestEm13", "electromagnetic/TestEm14",
    "electromagnetic/TestEm15", "electromagnetic/TestEm16", "electromagnetic/TestEm17",
    "electromagnetic/TestEm18",
    #
    "errorpropagation/errProp",
    #
    "eventgenerator/HepMC/HepMCEx01", "eventgenerator/HepMC/HepMCEx02",
    "eventgenerator/HepMC/MCTruth", "eventgenerator/exgps", "eventgenerator/particleGun",
    "eventgenerator/pythia/decayer6", "eventgenerator/pythia/py8decayer",
    "eventgenerator/userPrimaryGenerator",
    #
    "exoticphysics/channeling/ch0", "exoticphysics/channeling/ch1",
    "exoticphysics/channeling/ch2", "exoticphysics/channeling/ch3",
    "exoticphysics/dmparticle", "exoticphysics/monopole", "exoticphysics/phonon",
    "exoticphysics/saxs", "exoticphysics/ucn",
    #
    "field/BlineTracer", "field/field01", "field/field02", "field/field03", "field/field04",
    "field/field05", "field/field06",
    #
    "g3tog4/clGeometry",
    #    
    "geometry/transforms", "geometry/vecGeomNavigation",
    #
    "hadronic/FissionFragment", "hadronic/FlukaCern/ProcessLevel/CrossSection",
    "hadronic/FlukaCern/ProcessLevel/FinalState",
    "hadronic/Hadr00", "hadronic/Hadr01", "hadronic/Hadr02", "hadronic/Hadr03",
    "hadronic/Hadr04", "hadronic/Hadr05", "hadronic/Hadr06", "hadronic/Hadr07",
    "hadronic/Hadr08", "hadronic/Hadr09", "hadronic/Hadr10", "hadronic/NeutronSource",
    "hadronic/ParticleFluence/Calo", "hadronic/ParticleFluence/ConcentricSpheres",
    "hadronic/ParticleFluence/Layer", "hadronic/ParticleFluence/Sphere",
    #
    "medical/dna/AuNP", "medical/dna/UHDR", "medical/dna/chem1", "medical/dna/chem2",
    "medical/dna/chem3", "medical/dna/chem4", "medical/dna/chem5", "medical/dna/chem6",
    "medical/dna/clustering", "medical/dna/dnadamage1", "medical/dna/dnadamage2",
    "medical/dna/dnaphysics", "medical/dna/icsd", "medical/dna/jetcounter",
    "medical/dna/mfp", "medical/dna/microdosimetry", "medical/dna/microprox",
    "medical/dna/microyz", "medical/dna/neuron", "medical/dna/pdb4dna",
    "medical/dna/range", "medical/dna/slowing", "medical/dna/scavenger",
    "medical/dna/splitting", "medical/dna/spower", "medical/dna/svalue",
    "medical/dna/wholeNuclearDNA", "medical/dna/wvalue",
    #
    "medical/DICOM/DICOM1", "medical/DICOM/DICOM2", "medical/GammaTherapy",
    "medical/electronScattering", "medical/electronScattering2",
    "medical/fanoCavity", "medical/fanoCavity2", "medical/radiobiology",
    #
    "optical/LXe", "optical/OpNovice","optical/OpNovice2", "optical/wls",
    #
    "parallel/MPI/exMPI01", "parallel/MPI/exMPI02", "parallel/MPI/exMPI03",
    "parallel/MPI/exMPI04", "parallel/ThreadSafeScorers",
    #
    "parameterisations/Par01", "parameterisations/Par02", "parameterisations/Par03",
    "parameterisations/Par04", "parameterisations/gflash/gflash1",
    "parameterisations/gflash/gflash2", "parameterisations/gflash/gflash3",
    "parameterisations/gflash/gflasha",
    #
    "persistency/P01", "persistency/P02", "persistency/P03",
    "persistency/gdml/G01", "persistency/gdml/G02", "persistency/gdml/G03",
    "persistency/gdml/G04",
    #
    "physicslists/extensibleFactory", "physicslists/factory", "physicslists/genericPL",
    #
    "polarisation/Pol01",
    #
    "radioactivedecay/Activation", "radioactivedecay/rdecay01", "radioactivedecay/rdecay02",
    #
    "runAndEvent/RE01", "runAndEvent/RE02", "runAndEvent/RE03", "runAndEvent/RE04",
    "runAndEvent/RE05", "runAndEvent/RE06", "runAndEvent/RE07",
    #
    "visualization/movies", "visualization/perspective", "visualization/standalone",
    "visualization/userVisAction", "visualization/vtk"
]

for dir_name in basic_examples:
    generate('basic', dir_name)

for dir_name in extended_examples:
    generate('extended', dir_name)

# Move all outputs to an 'outputs' directory
outputs_dir = os.path.join(CURDIR, "outputs")
if not os.path.exists(outputs_dir):
    os.makedirs(outputs_dir)

for filename in os.listdir(CURDIR):
    if filename.endswith(".out") and filename.startswith("doxygen_"):
        source_path = os.path.join(CURDIR, filename)
        destination_path = os.path.join(outputs_dir, filename)

        if os.path.exists(destination_path):
            os.remove(destination_path)

        shutil.move(source_path, destination_path)

print("Doxygen generation script finished.")
