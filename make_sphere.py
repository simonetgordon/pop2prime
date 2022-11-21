import yt
import ytree
import os
import sys
import numpy as np

yt.enable_parallelism()

# Load all DDs in dir and load the MMhalo tree
# Find position and radius of each one and make a sphere.h5 file for each

# macros
radius_kpccm = 0.1

def halo_attributes(arbor, i):
    """
    Return mass, position and radius of halo at snapshot 560-i
    """
    # arbor = list(arbor[:])
    mass = arbor[i]["mass"].to('Msun')
    pos = arbor[i]["position"].to('unitary')
    pos = arbor.arr(pos.d, "unitary")
    rad = arbor[i]["virial_radius"].to('unitary')
    rad = arbor.quan(rad.d, "unitary")
    return mass, pos, rad

# define particle filter
@yt.particle_filter(requires=["particle_index"], filtered_type="nbody")
def _traced_dm(pfilter, data):
    bool_mask = data[(pfilter.filtered_type, "particle_index")] == np.any(dm_indices)
    print(bool_mask)
    return data[(pfilter.filtered_type, "particle_index")][bool_mask]


if __name__ == "__main__":

    # load data
    root_dir = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/cc_512_no_dust_continue"
    #dds = yt.load(os.path.join(root_dir, sys.argv[1])) # something like DD%04d/DD%04d for multiple DDs
    dm_indices = yt.load("ds_dm_particles.h5").data[("data", "particle_id")]
    a = ytree.load(os.path.join(root_dir, 'merger_trees/p2p_nd/p2p_nd.h5'))

    # final version
    # es = yt.load(sys.argv[1])
    # test version - time series
    ts = yt.DatasetSeries([os.path.join(root_dir, sys.argv[1]), os.path.join(root_dir, sys.argv[2])])
    #fns = ts.data['filename'].astype(str)[::-1]
    data_dir = root_dir

    storage = {}
    i = 0
    for store, ds in ts.piter(storage=storage):

        # Make sphere from refined region, centred on MMHalo centre
        # region = ds.box(ds.parameters["RefineRegionLeftEdge"],
        #                 ds.parameters["RefineRegionRightEdge"])

        sp_halo = ds.sphere(halo_attributes(a, i)[1], (radius_kpccm, 'kpc'))
        fields = [('nbody', 'particle_index'), ('nbody', 'particle_type'), ('gas', 'temperature'), ('gas', 'density'),
                  ("all", "particle_position_x"), ("all", "particle_position_y"), ("all", "particle_position_z"),
                  ("all", "particle_position"), ('nbody', 'particle_mass')]
        fn = sp_halo.save_as_dataset(fields=fields)  # save as dataset
        sphere_ds = yt.load(fn)

        # isolate the dm particles you're following
        yt.add_particle_filter(
            "traced_dm", function=_traced_dm, filtered_type="nbody", requires=["particle_index"]
        )
        sphere_ds.add_particle_filter("traced_dm")
        print("dm filter added")
        ad = sphere_ds.all_data()
        dm_mass_all = ad["nbody", "particle_mass"][ad["nbody", "particle_type"] == 1].to('Msun')
        print(len(dm_mass_all))
        store.result = (ds.current_time.in_units("Myr"), )

        i += 1

        # dm_ids = sp_halo["particle_index"][sp_halo["particle_type"] == 1]
        # print("number of dm particles in halo: ", len(dm_ids))  # 547,338 within 0.1 kpc