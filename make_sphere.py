"""
Produces dataset of traced DM particle positions across entire simulation.

Usage: mpirun -np 4 python make_sphere.py DD0558/DD0558 DD0559/DD0559

Requirements: ds_dm_particles.h5 in current directory and cc_512_no_dust_continue path provided in script.
"""
import yt
import ytree
import os
import sys
import numpy as np
import time
from os.path import exists

yt.enable_parallelism()

# Load all DDs in dir and load the MMhalo tree
# Find position and radius of each one and make a sphere.h5 file for each

# macros
radius_kpccm = 1
use_sphere = True

def halo_attributes(arbor, i):
    """
    Return mass, position and radius of halo at snapshot 560-i
    """
    mass = arbor[i]["mass"].to('Msun')
    pos = arbor[i]["position"].to('unitary')
    pos = arbor.arr(pos.d, "unitary")
    rad = arbor[i]["virial_radius"].to('unitary')
    rad = arbor.quan(rad.d, "unitary")
    return mass, pos, rad

# define particle filter
dm_indices = None


def _precious(pfilter, data):
    global dm_indices
    pids = data["nbody", "particle_index"].astype(int)
    if dm_indices is None:
        return np.zeros(pids.size, dtype=bool)
    return np.in1d(pids, dm_indices, assume_unique=True)


yt.add_particle_filter("precious", function=_precious,
                       filtered_type="nbody",
                       requires=["particle_index"])


if __name__ == "__main__":

    # load data
    root_dir = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/cc_512_no_dust_continue"
    # dds = yt.load(os.path.join(root_dir, sys.argv[1])) # something like DD%04d/DD%04d for multiple DDs
    dm_indices = yt.load("ds_dm_particles.h5").data[("data", "particle_id")]
    a = ytree.load(os.path.join(root_dir, 'merger_trees/p2p_nd/p2p_nd.h5'))

    # final version
    # es = yt.load(sys.argv[1])
    # test version - time series
    ts = yt.DatasetSeries([os.path.join(root_dir, sys.argv[1]), os.path.join(root_dir, sys.argv[2])])
    data_dir = root_dir


    storage = {}
    i = 0
    for store, ds in ts.piter(storage=storage):

        fields = [('nbody', 'particle_index'), ('nbody', 'particle_type'), ('gas', 'temperature'), ('gas', 'density'),
                  ("all", "particle_position_x"), ("all", "particle_position_y"), ("all", "particle_position_z"),
                  ("all", "particle_position"), ('nbody', 'particle_mass')]

        # Make sphere from refined region, centred on MMHalo centre
        # region = ds.box(ds.parameters["RefineRegionLeftEdge"],

        # isolate the dm particles you're tracing
        print(dm_indices)
        ds.add_particle_filter("precious") # takes a long time

        if use_sphere:
            print(str(ds) + "_sphere.h5")
            if not exists(str(ds) + "_sphere.h5"):
                print("making sphere data container")
                sp_halo = ds.sphere(halo_attributes(a, i)[1], (radius_kpccm, 'kpc'))
                sphere = sp_halo.save_as_dataset(fields=fields)  # save as dataset
            else:
                sphere = str(ds) + "_sphere.h5"
            print("sphere: ", sphere)
            sphere_ds = yt.load(sphere)

        ad = sphere_ds.all_data()
        dm_mass_all = sphere_ds["precious", "particle_mass"].to('Msun')
        print(len(dm_mass_all))
        store.result = (ds.current_time.in_units("Myr"), )

        print("execution time (s): ", time.perf_counter()/60)
        i += 1
