__author__ = 'sebastians'

# PyMOL Animation Tutorials -- 04 - Animation
# This tutorial will show you how to connect the scenes and thus animate the movie.

# Import all necessary modules
import pymol
from pymol import cmd

# The custom viewpoint
closeup = '''
    -0.775189936,    0.267432511,   -0.572329581,\
     0.387867898,    0.916590214,   -0.097048827,\
     0.498639554,   -0.297219634,   -0.814257801,\
     0.000021780,   -0.000062047,  -62.138366699,\
    -3.786274910,   25.372997284,    6.908325195,\
    45.995002747,   78.286071777,  -20.000000000 '''


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


def set_up_scenes():
    """Set up the scenes for a global view and a close-up on the binding site"""
    cmd.zoom('Cathepsin', 10)  # Zoom out to get a view on the whole complex
    cmd.scene('001', 'store', message='This is the first scene with a view on the complex!')
    cmd.set_view(closeup)  # Get a close-up of the ligand by using the manually chosen viewpoint
    cmd.scene('002', 'store', message='This is the second scene with a close-up on the ligand!')

setup_pymol()
cmd.load('../../input-files/Cathepsin.pdb')  # Load the PDB file
initial_representations()
set_up_scenes()

# With the scenes ready, we need to assign them to frames, i.e. set how long a scene should appear in the animation
cmd.scene('001', animate=0)  # First, select the scene to assign
cmd.mview('store', 1)  # We assign it to the first frame, i.e. the movie should start showing the first scene
cmd.mview('store', 150)  # Also assign it to frame 200. Scene 001 is now shown for 150 frames
cmd.scene('002', animate=0)  # Now, choose the close-up scene
cmd.mview('store', 250)  # Assign it to frame 250
cmd.mview('store', 400)  # Also assign it to frame 400

# Using scenes, we don't have to bother about transitions (the camera flights) between scenes.
# PyMOL takes care of that. It interpolates a smooth transition between the scenes.
# Let's rewind to frame 1 and see how the animation looks. To get a fast preview, we have to turn off ray tracing.
cmd.rewind()  # Rewind the movie to frame 1
cmd.set('ray_trace_frames', 0)  # Turn ray-tracing for frames off


# By default, PyMOL shows the video in a loop and also takes care of the transition of the last to the first scene.
# This setting is perfect for producing videos for presentations that should loop over and over again.
