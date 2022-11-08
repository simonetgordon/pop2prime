import yt
import ytree
from yt.extensions.astro_analysis.halo_analysis import HaloCatalog
import os
from os.path import exists
import matplotlib.pyplot as plt
import random

plot_halos = False
plot_mass_z = False
print_halo_attributes = False
new_param = False
plot_particles = False
radius_kpc = 0.2
select_random_dm_particles = True

# load data
root_dir = "/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/cc_512_no_dust_continue"
a = ytree.load(os.path.join(root_dir, 'merger_trees/p2p_nd/p2p_nd.h5'))
ds = yt.load(os.path.join(root_dir, "DD0560/DD0560"))

# useful functions
def most_massive_halo_tree(a, write_out=False):
    """
    Find most massive halo and return its tree object
    """
    my_tree = a[a["mass"].argmax()]
    prog_mass = my_tree["prog", "mass"]
    prog_redshift = my_tree["prog", "redshift"]
    if write_out:
        print("progenitor masses: ", prog_mass)
        print("progenitor redshifts: ", prog_redshift)
    return my_tree

def halo_attributes(arbor, ds):
    """
    Return mass, position and radius of halo at last snapshot
    """
    mass = arbor[0]["mass"].to('Msun')
    pos = arbor[0]["position"].to('unitary')
    pos = arbor.arr(pos.d, "unitary")
    rad = arbor[0]["virial_radius"].to('unitary')
    rad = arbor.quan(rad.d, "unitary")
    return mass, pos, rad

""" Main """
if __name__ == "__main__":
    # load most massive halo tree object
    if exists("tree_43576459"):
        arbor = ytree.load("tree_43576459/tree_43576459.h5")
    else:
        arbor = ytree.load(most_massive_halo_tree(a).save_tree())

    # print halo attributes
    if print_halo_attributes:
        print("==================================================")
        print("At z = {}, time = {}".format(ds.current_redshift, ds.current_time.in_units("Myr")))
        print("Halo mass:     ", halo_attributes(arbor, ds)[0])
        print("Halo position: ", halo_attributes(arbor, ds)[1])
        print("Halo radius:   ", halo_attributes(arbor, ds)[2])
        print("==================================================")

    #sp_halo.write_out("sphere_0.2kpc.txt") #  fields=[('nbody', 'particle_index')]
    if new_param is False:
        sphere_ds = yt.load("DD0560_sphere.h5")
        ad = sphere_ds.all_data()
        dm_ids = ad[("nbody", "particle_index")]
        print("dm ids: ", dm_ids)
        print("number of dm particles: ", len(dm_ids))
    else: # make sphere
        sp_halo = ds.sphere(halo_attributes(arbor, ds)[1], (radius_kpc, "kpc"))
        fields = [('nbody', 'particle_index'), ('nbody', 'particle_type'), ('gas', 'temperature'), ('gas', 'density'),
                  ("all", "particle_position_x"), ("all", "particle_position_y"), ("all", "particle_position_z"),
                  ("all", "particle_position")]
        fn = sp_halo.save_as_dataset(fields=fields) # save as dataset
        dm_ids = sp_halo["particle_index"][sp_halo["particle_type"] == 1]
        print("dm ids: ", dm_ids)
        print("number of dm particles: ", len(dm_ids))


    """ Plotting """

    # plot halo mass vs z of my_tree
    if plot_mass_z:
        prog_mass = most_massive_halo_tree(a)["prog", "mass"]
        prog_redshift = most_massive_halo_tree(a)["prog", "redshift"]

        plt.ylabel(r'Halo Mass $M_\odot$')
        plt.xlabel(r'Redshift z')
        plt.semilogy(prog_redshift, prog_mass)
        plt.savefig('mass_vs_redshift.png')

        # merger tree plot
        p = ytree.TreePlot(most_massive_halo_tree(a))
        p.min_mass_ratio = 0.1
        p.save('tree_most_massive_halo.png')


    # Make plot of full simulation volume and circle halos
    if plot_halos:
        hc = HaloCatalog(data_ds=ds, finder_method='hop')
        hc.create()

        prj_halos = yt.ProjectionPlot(ds, 'z', 'all_cic', data_source=ds)
        prj_halos.annotate_halos(hc)
        prj_halos.save('ProjectionPlot_halos_circled.png')

    if plot_particles:
        if new_param:
            data_source = sp_halo
        else:
            data_source = ds.sphere(halo_attributes(arbor, ds)[1], (radius_kpc, "kpc"))
        p = yt.ParticleProjectionPlot(ds, "z", color="g", data_source=data_source)
        p.set_unit("particle_mass", "Msun")
        p.save('particle_plot.png')

    if select_random_dm_particles:
        """
        Chooses at random 1% of all dm particles in the sp_halo region.
        """
        random.seed(0)
        if not new_param:
            sp_halo = ds.sphere(halo_attributes(arbor, ds)[1], (radius_kpc, "kpc"))

        # grab all dm postions and ids
        dm_pos_all = sp_halo[("all", "particle_position")][sp_halo["particle_type"] == 1].d
        dm_ids_all = sp_halo["particle_index"][sp_halo["particle_type"] == 1]

        # total number of particles in sphere
        no_particles = int(len(dm_ids_all)*0.01)
        print(no_particles)

        # sample random indices and form a list of dm positions from these ids
        # this assumes that particle_index[3] has particle_position[3]
        dm_indices = random.sample(range(len(dm_ids_all)-1), no_particles)
        dm_pos = [dm_pos_all[i] for i in dm_indices]
        #print("dm positions: ", dm_pos)

        with open("dm_particles.csv", mode="w") as f:
            col_format = "{:<5}" * 2 + "\n"  # 2 left-justfied columns with 5 character width
            col_format = "{:<15} , {1}"
            col_format = "{0}, {1}" + "\n"
            for (dm_indices, dm_pos) in zip(dm_indices, dm_pos):
                f.write(col_format.format(dm_indices, dm_pos))
            print('file created')
