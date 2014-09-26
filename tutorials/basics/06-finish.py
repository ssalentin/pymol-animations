__author__ = 'sebastians'

# PyMOL Animation Tutorials -- 06 - Finish
# This tutorial will show you how to pretty up your animation and output a movie

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


def additional_camera():
    # Add a rotation to the global view
    cmd.scene('001', animate=0)
    cmd.turn('y', -40)
    cmd.mview('store', 80)
    cmd.turn('y', 40)
    cmd.mview('store', 140)
    # Add panning to the close-up
    cmd.scene('002', animate=0)
    cmd.move('x', 5)
    cmd.mview('store', 320)


setup_pymol()
cmd.load('../../input-files/Cathepsin.pdb')  # Load the PDB file
initial_representations()
set_up_scenes()
scenes_to_frames()
additional_camera()

# Rewind and turn off ray-tracing
cmd.rewind()  # Rewind the movie to frame 1

#cmd.save('/tmp/movie_session.pse')  # It's also a good idea to save the session
cmd.set('ray_trace_frames', 1)  # Turn ray-tracing for frames on this time to get a nice output

#cmd.save('/tmp/movie_session.pse')  # It's also
# Choose a nice-fitting colorset for your molecules and general setup.
# For a single complex, one good choice is a light green for the protein, orange for the ligand and coloring by atoms.
cmd.color('palegreen', 'Cathepsin')  # Light green protein
cmd.color('tv_orange', 'ligand')  # Orange ligand
cmd.util.cnc('ligand')  # Color ligand atoms by type
cmd.set('bg_rgb', 'white')  # A white background is especially suited for presentations

# Render the frames into png images with the next command
# The output folder here would be /tmp and the prefix for the png files 'movie'
# The png files are called movie0001.png, movie0002.png, etc.
#cmd.mpng('/tmp/movie')  # Uncomment this to render the movie when running the script

# Before you can watch the movie, you need to stick together the png files into a movie
# This can be done with free tools like ffmpeg
# To create a movie from the png files in this tutorial, tun the following command (on Linux):
# ffmpeg -f image2 -i movie%04d.png -r 25 -sameq movie.mp4
# This will create a movie file called movie.mp4 with 25 frames per second and the same quality as the images

# Now it's time to get some popcorn! :)

