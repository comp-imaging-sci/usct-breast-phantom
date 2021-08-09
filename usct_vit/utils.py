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



from .config import *
import numpy.fft as fft
import math
import scipy.stats as stats
import gzip
import numpy as np
import scipy.io as sio
import os
import h5py
import hdf5storage


def GetVolume(_path, phantom_id, zz, thickness ):
    '''
    read the raw data of fda phantom
    fdaphantom files description:
    *.mhd: the header file includes phantom size
    p_*raw.gz: phantom label data without tumor
    
    Input:
    _path: directory that contains raw phantom data
    phantim_id: the seed number of phantom data.
    zz: the target slice for extraction
    thickness: the thickness of 3d phantom
    Output: a (2*thickness+1) layers label phantom
    '''
    headerFile = os.path.join(_path, 'p_'+phantom_id+'.mhd');
    #rawFile = os.path.join(unzip_out, 'p_'+seed+'.raw');
    rawgzFile = os.path.join(_path, 'p_'+phantom_id+'.raw.gz');
    cmdout=os.popen('cat '+headerFile+' |grep DimSize').read()
    cmdout=cmdout.split()
    xDim = int(cmdout[2]);
    yDim = int(cmdout[3]);
    zDim = int(cmdout[4]);
    print ("VICTRE dims: ", xDim,yDim,zDim)

    #xDim =100
    fid = gzip.open(rawgzFile,'rb')
    _buffer = fid.read(xDim*yDim*zDim)
    fid.close()
    volume = np.fromstring(_buffer, 'uint8')
    volume = np.reshape(volume, (xDim, yDim, zDim), order='F')
    #check target slice range
    if zz==-1 or thickness==-1:
        print ('Generate whole 3d volume data')
        return volume

    print ('Generate 2D slice or 3D slab')
    assert(zz>=0 and zz<zDim) # need a reasonable target slice number
    lb = zz-thickness;
    ub = zz+thickness;
    assert(lb>=0) # out of the bounds
    assert(ub<zDim) # out of the bounds
    
    lb = max(lb-2, 0);
    ub = min(ub+2, zDim-1)
    
    return volume[lb:ub+1,:,:]



def Labelprocessing3d(volume):
    '''
    Remove extra labels
    Input,
    volume: 3d label data
    
    Output: cleaned 3d label data
    '''
    #keep the slice above and beblow the tumor slice for the next step
    volume[np.where(volume==Labels['TDLU'])] = Labels['Glandular']
    volume[np.where(volume==Labels['Duct'])] = Labels['Glandular']
    volume[np.where(volume==Labels['Nipple'])] = Labels['Skin']
    #volume = RemoveLabel(volume, Labels['Ligament'])
    volume = RemoveLabel(volume, Labels['Artery'])
    volume = RemoveLabel(volume, Labels['Vein'])
    # print np.unique(volume)
    #check if extra labels are removed clearly
    labelset = np.unique(volume)
    itern = 0
    while 1:
        # until all extral labels are removed
        volume = RemoveLabel(volume, Labels['Artery'])
        volume = RemoveLabel(volume, Labels['Artery'])
        volume = RemoveLabel(volume, Labels['Vein'])
        volume = RemoveLabel(volume, Labels['Vein'])
        ar = np.sum(volume[2:-2,:,:]==Labels['Artery'])
        ve = np.sum(volume[2:-2,:,:]==Labels['Vein'])
#        print (ar)
        if ar==0 and ve==0:
            break
    print (np.unique(volume[2:-2,:,:]))
    return volume[2:-2,:,:]


def RemoveLabel(img, label):
    '''
    Replace the extra label by its neigbors
    first, counting the labels frequency in adjacent 18 voxels
    18 voxels include 5 voxels at above slice, 5 voxels at the below slice, 8 voxels at the current slice
    then, replace the extra label by the highest frequency label
    NOTICE: this method can't guarantee the extra label can be removed clearly
    this function will be executed several times utils all extra labels are removed.  
    '''
    z,x,y = img.shape
    z_list,x_list,y_list=np.where(img==label)
    new_label_link = [] #list for saving tuples (3d positon,replaced label)
    for idx in range(len(x_list)):
        ii = z_list[idx]
        jj = x_list[idx]
        kk = y_list[idx]
        if ii==0 or ii==z-1:
            continue
        # z = 0
        neighbor = []#list for saving the labels frequency of neighbors 
        if jj!=x and img[ii,jj+1,kk]!=label:
            neighbor.append(img[ii,jj+1,kk])
        if jj!=0 and img[ii,jj-1,kk]!=label:
            neighbor.append(img[ii,jj-1,kk])
        if kk!=0 and img[ii,jj,kk-1]!=label:
            neighbor.append(img[ii,jj,kk-1])
        if kk!=y and img[ii,jj,kk+1]!=label:
            neighbor.append(img[ii,jj,kk+1])
        if (kk!=y and jj!=x) and img[ii,jj+1,kk+1]!=label:
            neighbor.append(img[ii,jj+1,kk+1])
        if (kk!=y and jj!=0) and img[ii,jj-1,kk+1]!=label:
            neighbor.append(img[ii,jj-1,kk+1])
        if (kk!=0 and jj!=x) and img[ii,jj+1,kk-1]!=label:
            neighbor.append(img[ii,jj+1,kk-1])
        if (kk!=0 and jj!=0) and img[ii,jj-1,kk-1]!=label:
            neighbor.append(img[ii,jj-1,kk-1])
        # z = 1
        if img[ii+1,jj,kk]!=label:
            neighbor.append(img[ii+1,jj,kk])
        if jj!=x and img[ii+1,jj+1,kk]!=label:
            neighbor.append(img[ii+1,jj+1,kk])
        if jj!=0 and img[ii+1,jj-1,kk]!=label:
            neighbor.append(img[ii+1,jj-1,kk])
        if kk!=0 and img[ii+1,jj,kk-1]!=label:
            neighbor.append(img[ii+1,jj,kk-1])
        if kk!=y and img[ii+1,jj,kk+1]!=label:
            neighbor.append(img[ii+1,jj,kk+1])
        # z = -1
        if img[ii-1,jj,kk]!=label:
            neighbor.append(img[ii-1,jj,kk])
        if jj!=x and img[ii-1,jj+1,kk]!=label:
            neighbor.append(img[ii-1,jj+1,kk])
        if jj!=0 and img[ii-1,jj-1,kk]!=label:
            neighbor.append(img[ii-1,jj-1,kk])
        if kk!=0 and img[ii-1,jj,kk-1]!=label:
            neighbor.append(img[ii-1,jj,kk-1])
        if kk!=y and img[ii-1,jj,kk+1]!=label:
            neighbor.append(img[ii-1,jj,kk+1])
        if len(neighbor)>0:
            newlabel = max(set(neighbor), key = neighbor.count) # get the highest frequency one
            new_label_link.append(((ii,jj,kk), newlabel))  
            #saving the link between position and label
            #img[ii,jj,kk] = newlabel;
    
    #replace labels
    for item in new_label_link:
        pos, newlabel = item
        ii,jj,kk = pos
        img[ii,jj,kk] = newlabel
    return img

def SetPropValue(Prop, tissue):
    '''
    Assign acoustic properties to label data
    Acoustic value is sampling from turncated gaussian distribution
    mean values, std values are defined in config.py file
    '''
    if type(Prop)==int or type(Prop)==float:
        return Prop
    mu = float(Prop['mean'])
    sigma = float(Prop['sd'])
    lw = float(Prop['min'])
    up = float(Prop['max'])
    X = stats.truncnorm((lw-mu)/sigma, (up-mu)/sigma, loc=mu, scale=sigma)
    val = X.rvs(1)
    return float(val)







def sampler2D(b, kappa, h):
    '''
    generate texture signal by gauss spectral function
    more details in paper
    Franceschini E, Mensah S, Amy D, Lefebvre JP. A 2-D anatomic breast ductal computer phantom for ultrasonic imaging. IEEE Trans Ultrason Ferroelectr Freq Control. 2006;53(7):1281-1288. doi:10.1109/tuffc.2006.1665076
    b: input whitenoise
    kappa: correlation length [mm]
    h: pixel size[mm]
    '''
    bhat = fft.fftn(b);
    nx,ny = b.shape;
    Lx = nx*h;
    Ly = ny*h;
    ky = np.concatenate((np.linspace(0,int(ny/2)-1,int(ny/2)),np.linspace(int(-ny/2),-1,int(ny/2))), axis=None)
    kx = np.concatenate((np.linspace(0,int(nx/2)-1,int(nx/2)),np.linspace(int(-nx/2),-1,int(nx/2))), axis=None)
    ky = (2*math.pi/float(Ly))*ky
    kx = (2*math.pi/float(Lx))*kx#[np.linspace(0,nx/2-1,nx/2), np.linspace(nx/2,-1,nx/2)]
    kkx, kky = np.meshgrid(ky,kx)
    #d = np.power(((np.power(kkx,2) + np.power(kky,2)) + 1./(kappa*kappa)),alpha_over_2);
    d2 = np.exp(-kappa*kappa*(np.power(kkx,2) + np.power(kky,2))/8)
    m = np.real(fft.ifftn(bhat*d2))
    return m

def sampler3D(b, kappa, h):
    '''
    generate texture signal by gauss spectral function
    more details in paper
    Franceschini E, Mensah S, Amy D, Lefebvre JP. A 2-D anatomic breast ductal computer phantom for ultrasonic imaging. IEEE Trans Ultrason Ferroelectr Freq Control. 2006;53(7):1281-1288. doi:10.1109/tuffc.2006.1665076
    b: input whitenoise
    kappa: correlation length [mm]
    h: pixel size[mm]
    '''
    bhat = fft.fftn(b);
    nx,ny,nz = b.shape;
    Lx = nx*h;
    Ly = ny*h;
    Lz = nz*h;
    ky = np.concatenate((np.linspace(0,int(ny/2)-1,int(ny/2)),np.linspace(int(-ny/2),-1,int(ny/2))), axis=None)
    kx = np.concatenate((np.linspace(0,int(nx/2)-1,int(nx/2)),np.linspace(int(-nx/2),-1,int(nx/2))), axis=None)
    kz = np.concatenate((np.linspace(0,int(nz/2)-1,int(nz/2)),np.linspace(int(-nz/2),-1,int(nz/2))), axis=None)
    ky = (2*math.pi/float(Ly))*ky
    kz = (2*math.pi/float(Lz))*kz
    kx = (2*math.pi/float(Lx))*kx#[np.linspace(0,nx/2-1,nx/2), np.linspace(nx/2,-1,nx/2)]
    kkx, kky, kkz = np.meshgrid(ky,kx,kz)
    #d = np.power(((np.power(kkx,2) + np.power(kky,2)) + 1./(kappa*kappa)),alpha_over_2);
    d2 = np.exp(-kappa*kappa*(np.power(kkx,2) + np.power(kky,2) + np.power(kkz,2))/8)
    m = np.real(fft.ifftn(bhat*d2))
    return m


def AddTexture3D(sos, density, label):
    '''
    add texture to sos and density map
    Input:
    sos: sos map
    density: density map
    label: label map
    '''
  
    vshape = sos.shape
    print (vshape)
    #gland_dens_rn = np.random.normal(0,1,vshape)
    mean = 0; sigma =1; lw = mean-0.9*sigma; up = mean+0.9*sigma;
    X = stats.truncnorm((lw-mean)/sigma, (up-mean)/sigma, loc=mean, scale=sigma)
    #b3 = np.random.normal(0,1,(2560,2560))
    #b4 = np.random.normal(0,1,(2560,2560))
    indices_fat = np.where(label==Labels['Fat'])
    indices_gland = np.where(label==Labels['Glandular'])

    #fat_sos_rn = X.rvs(vshape)
    #fat_dens_rn = X.rvs(vshape)
    kappa = 0.21; h = 0.1 #mm

    rn = np.random.normal(0,1,vshape)
    if len(rn.shape)==2:
        text = sampler2D(rn, kappa, h)
    else:
        text =sampler3D(rn, kappa, h) # for gland sos
    text = text*1451*0.02
    sos[indices_gland] = sos[indices_gland]+text[indices_gland] #gland sos
    

    rn = np.random.normal(0,1,vshape)
    if len(rn.shape)==2:
        text = sampler2D(rn, kappa, h)
    else:
        text =sampler3D(rn, kappa, h)  # for gland dens
    text = text*999*0.02
    density[indices_gland] = density[indices_gland]+text[indices_gland] #gland dens

    rn = X.rvs(vshape)
    if len(rn.shape)==2:
        text = sampler2D(rn, kappa, h)
    else:    
        text =sampler3D(rn, kappa, h)  # for fat sos
    text = text*1420*0.02
    sos[indices_fat] = sos[indices_fat]+text[indices_fat] # fat sos

    rn = X.rvs(vshape)
    if len(rn.shape)==2:
        text = sampler2D(rn, kappa, h)
    else:
        text =sampler3D(rn, kappa, h)  # for fat dens
    # texture = texture;
    # adjust std to 2%:
    # mapping current disribution to truncated guassian distribution
    # sigma1 = 1440*0.02; sigma2 = 911*0.02
    text = text*915*0.02
    density[indices_fat] = density[indices_fat]+text[indices_fat] #fat dens
    return sos, density

