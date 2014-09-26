__author__ = 'sebastians'

# PyMOL Animation Tutorials -- 03 - Manual Viewpoints
# This tutorial will show you how to refine the camera position with manual viewpoints

# Import all necessary modules
import pymol
from pymol import cmd


def setup_pymol():
    """Sets up PyMOL for making animations."""
    pymol.finish_launching()  # Prevent threading errors
    # Configure global settings
    cmd.set('scene_buttons', 1)
    cmd.set('matrix_mode', 1)
    cmd.set('movie_panel', 1)
    # Configure quality settings
    cmd.mset("1 x500")
    cmd.set('ray_trace_frames', 1)
    cmd.viewport(800, 800)


def initial_representations():
    """Configure the initial representations for the protein and the ligand"""
    cmd.hide('everything', 'all')
    cmd.show('cartoon', 'all')
    cmd.select('ligand', 'resn NFT')
    cmd.deselect()
    cmd.show("sticks", "ligand")

setup_pymol()
cmd.load('../../input-files/Cathepsin.pdb')  # Load the PDB file
initial_representations()

# After loading the file and setting up the initial representations, let's have a look at the viewpoints again.
# The first one we set with cmd.zoom is fine, showing the whole complex
cmd.zoom('Cathepsin', 10)  # Zoom out to get a view on the whole complex
cmd.scene('001', 'store', message='This is the first scene with a view on the complex!')  # Save the first scene (001)

# The second one, however, showed the back of the binding site with NFT.
# In the normal case, we want to look inside the binding pocket!
# If just one animation should be produced, it's a good choice to select the viewpoint manually.
# Run the script as is and move the camera manually to get a good view on the ligand and the binding pocket.
# Type get_view to display the view matrix. Store it in a variable like this in the script:

closeup = '''
    -0.775189936,    0.267432511,   -0.572329581,\
     0.387867898,    0.916590214,   -0.097048827,\
     0.498639554,   -0.297219634,   -0.814257801,\
     0.000021780,   -0.000062047,  -62.138366699,\
    -3.786274910,   25.372997284,    6.908325195,\
    45.995002747,   78.286071777,  -20.000000000 '''

cmd.set_view(closeup)  # Get a close-up of the ligand by using the manually chosen viewpoint
cmd.scene('002', 'store', message='This is the second scene with a close-up on the ligand!')  # Save the second scene


# Running the script, we should now wee our custom viewpoint in scene 2
# Although there are ways to choose such viewpoints in a more automatic fashion, manual selection often yields much
# better results.
