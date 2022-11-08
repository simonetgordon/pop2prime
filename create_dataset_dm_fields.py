import numpy as np
import os
import sys
import yt
import ytree

from yt.extensions.p2p import \
    add_p2p_particle_filters

# load data
root_dir = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/cc_512_no_dust_continue"
a = ytree.load(os.path.join(root_dir, 'merger_trees/p2p_nd/p2p_nd.h5'))
ds = yt.load("DD0560_sphere.h5")

output_file = os.path.join("dm", "{ds.basename}.h5")
# if os.path.exists(output_file):

#add_p2p_particle_filters(ds)
region = ds.box(ds.parameters["RefineRegionLeftEdge"],
                ds.parameters["RefineRegionRightEdge"])

fields = [('nbody', 'particle_index'), ('nbody', 'particle_type'), ('gas', 'temperature'), ('gas', 'density'),
          ("all", "particle_position_x"), ("all", "particle_position_y"), ("all", "particle_position_z"),
          ("all", "particle_position")]
data = dict((field, region[field]) for field in fields)
ftypes = dict(field for field in fields)

yt.save_as_dataset(
    ds, filename=output_file,
    data=data, field_types=ftypes)