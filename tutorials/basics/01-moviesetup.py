__author__ = 'sebastians'

# PyMOL Animation Tutorials -- 01 - Movie Setup
# This tutorial will show you how to set up a python script for writing an PyMOL animation.

# To start PyMOL from a python script, import the pymol modules
import pymol
from pymol import cmd

# The next line is essential to prevent any threading errors. Call this function before using any PyMOl methods!
pymol.finish_launching()


# Next, we configure some global settings which are needed for our movies.
cmd.set('scene_buttons', 1)  # This installs scene buttons to switch between scenes
cmd.set('matrix_mode', 1)  # This sets up a mode needed for animations.
cmd.set('movie_panel', 1)  # This setting shows a slider to interactively move between frames

# Finally, we have to think about the length and quality of our animation
cmd.mset("1 x1000")  # How many frames do we need?. 1000 frames is enough for a 40 seconds video with 25 frames/second
cmd.set('ray_trace_frames', 1)  # Should the frames be ray-traced? If yes, rendering will take longer, but look nicer.
cmd.viewport(800, 800)  # Choose the desired resolution for the movie here (in pixels).


# Launching this script with Python should invoke PyMOL and show a black screen with everything set up for a movie.

# Note that in the cmd.mset() command we don't care about the first number, which has to do with the
# ratio of states and frames in the movie. Let's leave it at '1' for the moment.