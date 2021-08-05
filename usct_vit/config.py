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


'''
Config file for acoustic properties assigment
'''

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

'''
Property                   |  Uint
Speed of Sound             |  [m/s]
Density                    |  [kg/m^3]
Attenuation cofficient     |  [Np/m/Mhz^b]  alpha_0 here
Attenuatio  power          |
'''


SOS = {
    'Water'   :1500.0,
    'Fat'     :{'mean':1440.2, 'sd':20.9, 'min':1412, 'max': 1485},
    'Skin'    :{'mean':1555, 'sd':10, 'min':1530, 'max':1580},
    'Glandular' :{'mean':1520, 'sd': 10, 'min':1505, 'max':1550},
    'Tumor': 1548.0,
    'Ligament':{'mean':1440.0, 'sd':10, 'min':1410, 'max':1470}
}

Density ={
    'Water'   :994,
    'Fat'     :{'mean':911, 'sd':53, 'min':812, 'max': 961},
    'Skin'    :{'mean':1109, 'sd': 14, 'min':1100, 'max':1125},
    'Glandular' :{'mean':1041, 'sd': 45.3, 'min':990, 'max':1092},
    'Tumor':945,
    'Ligament':1142	
}

Atten = {
    'Water':0.025328436023,
    'Fat':4.3578,
    'Skin':21.158,
    'Glandular':8.635,
    'Tumor':31,
    'Ligament':14.506,
}
