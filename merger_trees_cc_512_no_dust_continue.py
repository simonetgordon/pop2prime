import yt
import ytree

a = ytree.load("/disk12/brs/pop2-prime/firstpop2_L2-Seed3_large/cc_512_no_dust_continue/merger_trees/p2p_nd/p2p_nd.h5")

# Find most massive halo                                                                                                                                                                                   
mass_array = a["mass"]

# Index of most massive halo
index = mass_array.argmax()

# Most massive halo tree called my_tree
my_tree = a[index]

p = ytree.TreePlot(my_tree)
p.min_mass_ratio = 0.1
p.save('tree_most_massive_halo.png')
