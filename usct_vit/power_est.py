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



import math
import numpy as np

def b_estimate(perc, rd=0.5):
    '''
    estimate the power exponent by 1d phantom given the percentage of fatty tissue
    perc: the percentage of fat
    rd: radius of phantom
    '''
    #print perc,rd
    f_list = np.linspace(0.1,2.3,23)
    y0_fat = 1.08;
    y0_gland = 1.5;

    a0_fat = 4.3578; #[Np/m/MHz]
    a0_gland =8.635;
    #rd = 0.05; #[m]
    rs = []
    #print f_list
    for i in range(len(f_list)):
        f0 = f_list[i]
        alpha_fat = a0_fat*(f0**y0_fat)
        alpha_gland = a0_gland*(f0**y0_gland)
        ratio = rd*alpha_fat*perc*2+rd*alpha_gland*(1-perc)*2;
        ratio = np.exp(-ratio)
        rs.append(ratio)

    #print rs
    b0 = 1.2
    rs = np.array(rs); f_list = np.array(f_list)
    while 1:
        m = rd*2*(a0_fat*perc+a0_gland*(1-perc));
        f1 =sum(np.power(np.exp(-np.power(f_list, b0)*m) - rs,2));
        g = sum((np.exp(-np.power(f_list, b0)*m) - rs)*np.exp(-m*np.power(f_list,b0))\
                *np.log(f_list)*np.power(f_list, b0)*(-m))
        # g = sum((np.power(f_list,b0)*m - rs)*np.log(f_list)*(np.power(f_list, b0)*m));
        step = 1;
        b1 = b0-step*g;
        while b1>1.5 or b1<1.08:
            step = step/2.
            b1 = b0-step*g
        f2 = sum(np.power(np.exp(-np.power(f_list, b1)*m) - rs,2));
        while f1<f2:
            step = step/2.
            #print step
            b1 = b0-step*g;
            f2 = sum(np.power(np.exp(-np.power(f_list, b1)*m) - rs,2));
            #f2 = sum(np.power(np.power(f_list,b1)*m - rs,2));
        #print f1,f2
        if f1-f2<0.000000000001 or abs(b1-b0)<0.000000001:
            # print f1-f2
            return b1
        #print b0,b1
        b0 = b1


    return 1.5

#print b_estimate2(0.5)
