# Anatomically and physiologically realistic numerical breast phantoms for USCT virtual imaging studies

This software project aims at generating stochastic numerical breast phantoms 
with anatomically realistic tissue structures and physiologically realistic
acoustic properties for use in ultrasound computed tomography (USCT) virtual imaging studies.

This software builds on top of tools for the _Virtual Imaging Clinical Trial for Regulatory Evaluation_ ((VICTRE)[https://github.com/DIDSR/VICTRE]) project
to generatee the numerical breast and tumor phantoms.

The output of the software consists of three-dimensional high-resolution maps of the breast acoustic properties, including speed-of-sound, density, and acoustic attenuation.


## The dependecies

numpy, scipy, h5py, hdf5storage

## Using the software

### 1. Preparation of tissue structure data(tissue label map)

An example of VICTRE phantom raw data is shared on Harvard dataverse. https://doi.org/10.7910/DVN/1KJK4G 
The example data for testing our code can be downloaded by:
```
wget https://dataverse.harvard.edu/api/access/datafile/4950808 -O ./data/Phantom_set/p_324402160.mhd
wget https://dataverse.harvard.edu/api/access/datafile/4950809 -O ./data/Phantom_set/p_324402160.raw.gz
```
The downloaded data will be saved in `./data/Phantom\_set`. 

Alternatively, users can create their own structure phantoms by the use of VICTRE software.
1. VICTRE (BreastPhantom)[https://github.com/DIDSR/breastPhantom]: a c++ opensource software for the generation of anthropomorphic breast phantoms
2. VICTRE (BreastMass)[https://github.com/DIDSR/breastMass]: a c++ opensource software for the generation of three-dimensional breast lesions

Where the folder `./data/Phantom\_cfg` contains the VICTRE config file templates for each breast type.

Two data are required for the executation of the next step:

- p_{phantom_id}.raw.gz: A compressed 3D anatomical data.

- p_{phantom_id}.mhd: A header file for this phantom that contains the physical location of the voxel with index [0, 0, 0], the total number of voxels in the x, y, and z, and the voxel size.

### 2. Assignment of acoustic properties

In this step
1. We read the 3d tissue labels map 
2. We remove labels corresponding to tissues that are invisible in USCT imaging
3. We assign stochastic acoustic properties maps

Before executing this step, make sure all paths for data loading are set correctly

To excute this step
```sh
python run_assign_properties.py -t <phantom_type> -mass <tumor> -rseed <seed_number>
python3 run_assign_properties.py -phantom_id <the digit identifier of the breast phantom> -raw_data_path <the data folder contains tissue structure data> -target_slice <the central slice of the generated phantom> -thickness <the thickness of the phantom to be generated> -output_path <output path>
```
where

- `phantom_id` is the digit identifier of the breast phantom used to generate the anatomical structures which is included in raw data file name.
- `raw_data_path` is the folder with raw anatomical data and header file
- `target_slice` is the target slice to be extracted in 3D anatomical data.
- `thickness` is the thickness of acoustic phantom to be generated (z-axis).
- `output_path` is the folder for saving the output data.

An example script is given in file `./run_assign_properties.sh`


## Data formats

### Tissue label maps

- `*.mhd`:            the header file includes phantom size
- `p_*raw.gz`:        phantom anatomical data 

```python
    Labels = {
        'Water':         0,
        'Fat':           1,
        'Skin':          2,
        'Glandular':    29,
        'Nipple':       33,
        'Muscle':       40,
        'Ligament':     88,
        'TDLU':         95,
        'Duct':        125,
        'Artery':      150,
        'Tumor':       200,
        'Vein':        225,
    }
```

Of these labels, only `Water`, `Fat`, `Skin` `Glandular`, `Ligament`, and `Tumor` are preserved.
 
### Acoustic properties maps  (2D)

   	speed of sound map:      a mat file with data type float32,unit mm/\mus\\
    density map:             a mat file with data type float32,unit kg/mm^3\\
    attenuation map:         a mat file with data type float32,unit dB/mm/Mhz^y\\
    label map:               a mat file with data type unit8 \\
<!---  
## Files description

Code file:

**acoustic_phantom/config.py**:  The configuration file include rawdata directory and acoustic properties range.\
**acoustic_phantom/power_est.py**: Functions for attenuation exponent estimation.\
**acoustic_phantom/utils.py**: Functions for perprossing steps.\
**main.py**: The execute file for reading rawdata and generating acoustic properties map.

Datafile:

**Loc_1024_makeCartCircle.DAT**: Locations (2D coordinates measured in mm) of the1024 transducers arranged in a circular array of radius 110mm.  The datadimension is 2*1024 (transducer 0 x-coordinate, transducer 0 y-coordinate,transducer 1 x-coordinate, transducer 1 y-coordinate, . . . )  and data typeisfloat32.\
**Phantom_source500.DAT**:  The shape of the input pulse used to generate the sono-gram.  It consists of 500 time samples at a frequency of 30MHz.  Data isstored asfloat32.\
**ws_kspace_tmp.ini:** The template configuration file for wavesolver\

  
### Attenuation exponent estimation

the attenuation degree is computed by Beer's law:

 P\_1 = p\_0 \*exp(-\alpha\_1\*length\_1-\alpha\_2\*length\_2)
 
 P\_2 = p\_0 \*exp(-\alpha\*length)  using constant exponent 
 
 alpha = alpha\_0\**f**^b;\
  alpha\_1 = alpha\_1\**f*^b\_1\
   alpha\_2 = alpha\_2\**f*^b\_2\
   
   Given list of frequency *f* based on source specturm.
   
   Find b which can minimaze the mismatch \sum|| p\_1 - p\_2 ||^2
   
  **Regression results in two cases**
  ![fig1](https://gitlab.engr.illinois.edu/anastasio-lab/usct/fdaphantom-acoustic-properties-assignment/-/blob/master/docs/2.png)
  ![fig2](https://gitlab.engr.illinois.edu/anastasio-lab/usct/fdaphantom-acoustic-properties-assignment/-/blob/master/docs/1.png)
 -->
    


