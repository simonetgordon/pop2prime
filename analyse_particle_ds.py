from __future__ import print_function
import yt
import ytree
from yt.extensions.astro_analysis.halo_analysis import HaloCatalog
from yt.data_objects.particle_filters import add_particle_filter
import os
import sys
import itertools
import numpy as np
from scipy.spatial.distance import pdist
import numpy as np
import math

# macros
print_basic_info = True

# load dm particle field .h5 data file
root_dir = "/home/sgordon/pop2prime/analysis_scripts"
ds = yt.load(os.path.join(root_dir, sys.argv[1]))


def min_max_magnitude(pos):
    #row_sums = np.sum(pos, axis=1).tolist().sort()
    pos_mag_max = np.max(np.apply_along_axis(pos_magnitude, 1, pos))
    pos_mag_min = np.min(np.apply_along_axis(pos_magnitude, 1, pos))
    return pos_mag_min, pos_mag_max

def pos_magnitude(pos):
    return np.sqrt(pos[0] ** 2 + pos[1] ** 2 + pos[2] ** 2)

def velocity_magnitude(vx, vy, vz):
    return np.sqrt(vx ** 2 + vy ** 2 + vz ** 2)


# print basic attributes
print("===========================================================================================")
print("At z = {}, time = {:.2f}".format(ds.current_redshift, ds.current_time.in_units("Myr")))
print("Number of DM particles: {}".format(len(ds.data[("data", "temperature")])))
print("Max temp: {:.2f}, Min temp: {:.2f}, Mean temp: {:.2f}".format(np.max(ds.data[("data", "temperature")]),
                                                                     np.min(ds.data[("data", "temperature")]),
                                                                     np.mean(ds.data[("data", "temperature")])))
print("Max density: {:.3g}, Min density: {:.3g}, Mean density: {:.3g}".format(np.max(ds.data[("data", "density")]),
                                                                              np.min(ds.data[("data", "density")]),
                                                                              np.mean(ds.data[("data", "density")])))
print("Max velocity: {:.2f}, Min velocity: {:.2f}, Mean velocity: {:.2f}".format(
    np.max(velocity_magnitude(ds.data[("data", "velocity_x")],
                              ds.data[("data", "velocity_y")],
                              ds.data[("data", "velocity z")])).in_units("m/s"),
    np.min(velocity_magnitude(ds.data[("data", "velocity_x")],
                              ds.data[("data", "velocity_y")],
                              ds.data[("data", "velocity z")])).in_units("m/s"),
    np.mean(velocity_magnitude(ds.data[("data", "velocity_x")],
                               ds.data[("data", "velocity_y")],
                               ds.data[("data", "velocity z")])).in_units("m/s")
))
print("Max distance: {:.2f}, Min distance: {:.2f}".format(
    min_max_magnitude(ds.data[("data", "positions")])[0],
    min_max_magnitude(ds.data[("data", "positions")])[1]
))
print("===========================================================================================")
