# Copyright (c) 2021,  University of Illinois Urbana-Champaign
# & Washington University in St Louis.
#
#
# This file is part of the usct-breast-phantom library. For more information and
# source code
# availability see https://github.com/comp-imaging-sci/usct-breast-phantom.
#
# usct-breast-phantom is free software; you can redistribute it and/or modify it
# under the
# terms of the GNU General Public License (as published by the Free
# Software Foundation) version 2.0 dated June 1991.






import argparse
import struct
from usct_vit import *
import h5py
import hdf5storage
import scipy

'''
This code assign acoustic properties to the NBP (2D slice, 3D slab, or full phantom)  
'''

def input_check():
    # check whether raw data file exist
    assert(os.path.exists(raw_data_path)) # The raw data file doesn't exist
    rawgzFile = os.path.join(raw_data_path, 'p_'+phantom_id+'.raw.gz')
    headerFile = os.path.join(raw_data_path, 'p_'+phantom_id+'.mhd')
    assert(os.path.isfile(rawgzFile)) #The raw data file doesn't exist
    assert(os.path.isfile(headerFile)) #The raw data file doesn't exist
    newfolder = os.path.join(output_path,phantom_id)
    if not os.path.exists(newfolder): #creat folder for saving data
        print('create a new folder',newfolder)
        os.makedirs(newfolder)
    
        
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-phantom_id', type=str, help="Digit identifier of breast phantom")
    parser.add_argument('-raw_data_path', type=str, help="Location of VICTRE raw data")
    parser.add_argument('-target_slice', type=int, default=-1, help="z-index of slice or mid point of slab to be extracted")
    parser.add_argument('-thickness', type=int, default=0, help="Thickness of the slab")
    parser.add_argument('-output_path', type=str, help="Output path")
    parser.add_argument('-resolution', type=float, default=0.1, help="Voxel size (mm)")


    args = parser.parse_args()
    phantom_id = args.phantom_id
    raw_data_path = args.raw_data_path
    target_slice = args.target_slice
    thickness = args.thickness
    output_path = args.output_path
    voxel_size  = args.resolution
    
    # -----------------------------------
    # 1. Read 3D phantom label data
    # -------------------------------------
    # Target slice and the thinkness
    input_check(); #check the correctness of the path
    volume = GetVolume(raw_data_path, phantom_id, target_slice, thickness)
    # -------------------------------
    # 2. Removel extral labels and extract the slice contain tumor
    # -------------------------------
    volume = Labelprocessing3d(volume)
    print(volume.shape)
    # overwrite the volume
    # -------------------------------
    # 3. Assign acoustic properties
    # -------------------------------
    map_sos = np.zeros((volume.shape),'float32')
    map_density = np.zeros((volume.shape),'float32')
    map_atten = np.zeros((volume.shape),'float32')
    for key in SOS:
        print (key)
        indices = np.where(volume==Labels[key])
        map_sos[indices] = SetPropValue(SOS[key], key)
        map_density[indices] = SetPropValue(Density[key], key)
        map_atten[indices] = SetPropValue(Atten[key], key)
    # -------------------------------
    # 4. Downsampling (from 0.05mm to 0.1mm)
    # -------------------------------
    downsampling_factor = 0.05/voxel_size
    map_sos = scipy.ndimage.zoom(map_sos, downsampling_factor)
    map_density = scipy.ndimage.zoom(map_density,downsampling_factor)
    map_atten = scipy.ndimage.zoom(map_atten, downsampling_factor)
    volume = scipy.ndimage.zoom(volume,downsampling_factor, mode='nearest')
    map_sos, map_density = AddTexture3D(map_sos, map_density, volume)
    map_sos = map_sos.astype('float32')
    map_atten = map_atten.astype('float32')
    map_density = map_density.astype('float32')

    # ---------------------------------
    # 7. save all the data
    # ---------------------------------

    newfolder = os.path.join(output_path,phantom_id)
    data = {'sos':map_sos}
    out_name = os.path.join(newfolder,'sos_'+phantom_id+'_z'+str(target_slice)+'.mat')
    if os.path.exists(out_name):
        os.remove(out_name)
    hdf5storage.write(data, filename=out_name, matlab_compatible=True)
    data= {'aa': map_atten}
    out_name = os.path.join(newfolder,'aa_'+phantom_id+'_z'+str(target_slice)+'.mat')
    if os.path.exists(out_name):
        os.remove(out_name)
    hdf5storage.write(data, filename=out_name, matlab_compatible=True)
    data = {'dd': map_density}
    out_name = os.path.join(newfolder,'density_'+phantom_id+'_z'+str(target_slice)+'.mat')
    if os.path.exists(out_name):
       	os.remove(out_name)
    hdf5storage.write(data, filename=out_name, matlab_compatible=True)
    data = {'label': volume}
    out_name = os.path.join(newfolder,'label_'+phantom_id+'_z'+str(target_slice)+'.mat')
    if os.path.exists(out_name):
       	os.remove(out_name)
    hdf5storage.write(data, filename=out_name, matlab_compatible=True)
