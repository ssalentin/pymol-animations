__author__ = 'sebastians'

# The function provided here can be used if you get no image output with cmd.png (can be no or a black picture).
# Can be also used if you experience segmentation faults with cmd.ray

from pymol import cmd
import os


def png_workaround(filepath, width=1024, height=768):
    """Workaround for (a) severe bug(s) in PyMOL preventing ray-traced images to be produced in command-line mode.
    Use this function in case neither cmd.ray() or cmd.png() work.
    """
    cmd.set('ray_trace_frames', 1)  # Frames are raytraced before saving an image.
    cmd.viewport(width, height)  # Set resolution
    ### Workaround for raytracing in command-line mode
    cmd.mpng(filepath, 1, 1)  # Use batch png mode with 1 frame only
    cmd.mplay()  # cmd.mpng needs the animation to 'run'
    os.rename("".join([filepath[:-4], '0001.png']), "".join([filepath[:-4], '.png']))  # Remove frame number in filename