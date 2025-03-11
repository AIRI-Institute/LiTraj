

### Download dataset

```python
from litraj.data import download_dataset
download_dataset('dataset_name', 'folder', unzip = True, remove_zip = True)
```



### Load dataset

```python
from litraj.data import download_dataset
index = load_data('dataset_name', 'folder')
```

#### Structure of the nebDFT2k dataset
    nebDFT2k/
    ├── nebDFT2k_index.csv              # Table with material_id, edge_id, chemsys, _split, .., columns
    ├── edge-id1_init.xyz               # Initial (BVSE-NEB optimized) trajectory file, edge_id1 = mp-id1_source_target_offsetx_offsety_offsetz
    ├── edge-id1_relaxed.xyz            # Final (DFT-NEB optimized) trajectory file 
    ├── ...
    └── nebDFT2k_centroids.xyz          # File with centroid supercells



#### Structure of the MPLiTrj dataset
    MPLiTraj/
    ├── MPLiTrj_train.xyz               # Training file
    ├── MPLiTrj_val.xyz                 # Validation file
    └── MPLiTrj_test.xyz                # Test file 



#### Structure of the MPLiTrj_subsample dataset
    MPLiTraj_subsample/
    ├── MPLiTrj_subsample_train.xyz      # Training file
    ├── MPLiTrj_subsample_val.xyz        # Validation file
    └── MPLiTrj_subsample_test.xyz       # Test file 

#### Structure of the MPLiTrj_raw dataset
    MPLiTraj_raw/
    ├── edge-id1.neb/                    # Folder with the raw files
    │   ├── band_optim_step_0.traj       # DFT-NEB trajectory optimization step 0
    │   ├── ...
    │   ├── band_optim_step_n.traj       # DFT-NEB trajectory optimization final n-th step (can be up to 100)
    │   ├── source.xyz                   # Source vacancy optimization trajectory
    │   ├── target.xyz                   # Target vacancy optimization trajectory
    │   ├── traj_init.xyz                # Initial (BVSE-NEB optimized) trajectory
    ├── ...
    └── edge-idn.neb/                    # Folder with the raw files

#### Structure of the BVEL13k dataset
    BVEL13k/
    ├── BVEL13k_index.csv                # .csv with material_id, chemsys, _split, E_1D, E_2D, E_3D columns
    ├── BVEL13k_train.xyz                # Training .xyz file 
    ├── BVEL13k_val.xyz                  # Validation .xyz file 
    └── BVEL13k_test.xyz                 # Test .xyz file 

#### Structure of the nebBVSE122k dataset
    nebBVSE122k/
    ├── nebBVSE122k_index.csv            # .csv with material_id, chemsys, _split, em columns
    ├── nebBVSE122k_train.xyz            # Training .xyz file 
    ├── nebBVSE122k_val.xyz              # Validation .xyz file 
    └── nebBVSE122k_test.xyz             # Test .xyz file 

