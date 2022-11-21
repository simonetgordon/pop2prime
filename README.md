# pop2prime
Post-processing of pop2prime simulation data to investigate BH growth. 

The `ds_dm_particles.h5` file is a dataset containing all the relevant fields of the chosen dark matter particles for \
calculating the Bondi-Hoyle accretion rate. The field dict is:

```bash
"particle_id": dm_indices,
"position": dm_pos,
"mass": dm_pos_mass,
"density": dm_pos_den,
"temperature": dm_pos_temp,
"velocity_x": dm_pos_velx,
"velocity_y": dm_pos_vely,
"velocity z": dm_pos_velz,
```

The dataset is queried like: 
```python
import yt
dm_ds = yt.load("ds_dm_particles.h5")
dm_indices = dm_ds.data[("data", "particle_id")]
```

Once the particles ids have been gathered from a certain snapshot (e.g the final), the particles are traced back in \
time. The mechanism for doing so is making a sphere centred on the halo and then querying the dm particles from that \
instead of the entire dataset (this is time-intensive).