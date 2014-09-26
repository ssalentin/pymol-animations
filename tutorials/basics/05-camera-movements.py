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


def scenes_to_frames():
    """Assign scenes to frames"""
    # Scene 001 from frames 1-150
    cmd.scene('001', animate=0)
    cmd.mview('store', 1)
    cmd.mview('store', 150)
    # Scene 002 from frames 250-400
    cmd.scene('002', animate=0)
    cmd.mview('store', 250)
    cmd.mview('store', 400)

setup_pymol()
cmd.load('../../input-files/Cathepsin.pdb')  # Load the PDB file
initial_representations()
set_up_scenes()
scenes_to_frames()

# Apart from moving between scenes, independent camera movements can be added to the animation
# One of them, the zoom (cmd.zoom) was already used
# Two other movements are rotation (cmd.turn) and panning (cmd.move)
# Let's add a rotation to the global view and a panning to the close-up to make them more interesting

cmd.scene('001', animate=0)  # Select the first scene again (we need the camera view as a starting point)
cmd.turn('y', -40)  # Turn the camera 40 degrees left around the y-axis
cmd.mview('store', 80)  # Store the changed view/assign it to frame 60
cmd.turn('y', 40)  # Turn the camera 40 degrees in the other direction around the y-axis
cmd.mview('store', 140)  # Assign the view to frame 120

cmd.scene('002', animate=0)  # Select the second scene
cmd.move('x', 5)  # Move along the x-axis
cmd.mview('store', 320)  # Assign the view to frame 320

# Rewind and turn off ray-tracing
cmd.rewind()  # Rewind the movie to frame 1
cmd.set('ray_trace_frames', 0)  # Turn ray-tracing for frames off

# Running the script, you should now see the additional camera movements additional to the scene transitions
# Note that there are also transitions back to original scene 001 and scene 002 views after the movements
