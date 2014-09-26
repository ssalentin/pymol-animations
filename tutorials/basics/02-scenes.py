__author__ = 'sebastians'

# PyMOL Animation Tutorials -- 02 - Scenes
# This tutorial will show you how to use scenes for your movie.

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


setup_pymol()

# After setting up PyMOL, we will load a PDB file an prepare it for a simple use case.
# Let's show the complex between Cathepsin K and one of its inhibitors (PDB ID 1VSN)

cmd.load('../../input-files/Cathepsin.pdb')  # Load the PDB file

# A good starting point for an animation is to hide all representation to start from scratch
# Then, we want to show the protein in cartoon representation and its ligand in ball-and-stick representation
cmd.hide('everything', 'all')  # Hide everything
cmd.show('cartoon', 'all')  # Show protein in cartoon representation
cmd.select('ligand', 'resn NFT')  # Select the ligand (NFT)
cmd.deselect()  # Deselect everything to hide the selection markers in PyMOL
cmd.show("sticks", "ligand")  # Show the ligand in ball-and-stick representation

# Now let's zoom out to show the whole protein and then do a close-up on the ligand.
# In both cases, we will save the camera perspective as a scene.
# We can switch between saved scenes and connect them later. It's just like in a movie storyboard.
cmd.zoom('Cathepsin', 10)  # Zoom out to get a view on the whole complex
cmd.scene('001', 'store', message='This is the first scene with a view on the complex!')  # Save the first scene (001)
cmd.zoom('ligand', 5)  # Get a close-up of the ligand
cmd.scene('002', 'store', message='This is the second scene with a close-up on the ligand!')  # Save the second scene


# Running the script, PyMOL should display you the complex with the chosen representations and two scenes
# The scenes are accessible in the lower left corner. You can click on the buttons to switch between them.
# We can use this feature later to connect the scenes and make a smooth transition for the movie.
# Note that the scene messages are optional, i.e. you can omit the message argument.