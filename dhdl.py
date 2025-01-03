"""Extract dH/dl values from .xvg files into pickle files"""
from alchemlyb.parsing.gmx import extract_dHdl, extract_u_nk
import os.path
from os.path import join
import pandas as pd
from alchemlyb.parsing import gmx
from alchemlyb.preprocessing import (slicing, statistical_inefficiency,
                                     equilibrium_detection)
from alchemlyb.estimators import TI
from alchemlyb.estimators import MBAR
import argparse
import pickle
import os
p = ['vdw', 'coul']
for x in p:
    if not os.path.isdir(x): os.mkdir(x)
    lambdanr=20
    for i in range(0, lambdanr):
        data_folder = os.path.join('lambda_%s' % i, 'Production_MD')
        xvg = os.path.join(data_folder, 'prod.xvg')
        dhdl = extract_dHdl(xvg, T=298.15)
        print(dhdl)
        vdw = dhdl[x].values.tolist()
        vdw = [i*0.59 for i in vdw]
        file = open('%s/dhdl_%s.pickle'%(x,i), 'wb')
        pickle.dump(vdw, file)
        file.close()
