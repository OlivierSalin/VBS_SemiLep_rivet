#!/usr/bin/env python3
import re
import pandas as pd
import os
import numpy as np
import itertools
import glob

from optparse import OptionParser

def diboson_br(process: str, decay: str):
    """
    process: 'WpZ','WmZ','WpWp','WpWm','WmWm','ZZ','Zy','Wpy' (charge doesn’t matter for BR)
    decay:   'llll','lllv','lvlv','lvqq','llqq','vvqq','lly','lvy'
    lepton_set: 'emu' or 'all' (all = e, mu, tau)
    returns: branching ratio as float
    """
    # PDG 2024/25 values
    W_e, W_mu, W_tau = 0.1071, 0.1063, 0.1138      # W → eν, μν, τν
    W_had = 1.0 - 3.0*0.1086                       # from universal fit → 0.6742
    Z_e, Z_mu, Z_tau = 0.03363, 0.03366, 0.03369   # Z → ee, μμ, ττ
    Z_had, Z_inv = 0.6991, 0.2000                  # Z → hadrons, Z → νν (sum)


    W_l = W_e + W_mu + W_tau
    Z_ll = Z_e + Z_mu + Z_tau


    # Map final-state tags to products of single-boson BRs
    if decay == "llll":   # ZZ
        return Z_ll * Z_ll
    if decay == "lllv":   # WZ
        return Z_ll * W_l
    if decay == "lvlv":   # WW
        return W_l * W_l
    if decay == "lvqq":   # WZ or WW semileptonic
        if "Z" in process:
            return W_l * Z_had
        else:
            return W_l * W_had
    if decay == "llqq":   # WZ or ZZ semileptonic with dilepton
        if "Z" in process and process.count("Z") == 1:
            return Z_ll * W_had
        else:
            return Z_ll * Z_had
    if decay == "vvqq":   # ZZ → νν qq
        return Z_inv * Z_had
    if decay == "lly":    # Zγ with Z → ℓℓ
        return Z_ll
    if decay == "lvy":    # Wγ with W → ℓν
        return W_l
    raise ValueError("Unknown decay tag")

valid_combinations = {
    1: [("WmZ", "lvqq"), ("WpZ", "lvqq"), ("WmWm", "lvqq"), ("WpWm", "lvqq"), ("WpWp", "lvqq")],
    2: [("WpZ", "llqq"), ("WmZ", "llqq"), ("ZZ", "llqq")],
    0: [("WpZ", "vvqq"), ("WmZ", "vvqq"), ("ZZ", "vvqq")]
}

valid_combi_aQGC = [f"{proc}_{dec}" for pairs in valid_combinations.values() for proc, dec in pairs]

Branching_ratios = {}
for key in valid_combinations:
    for proc, dec in valid_combinations[key]:
        br = diboson_br(proc, dec)
        Branching_ratios[f"{proc}_{dec}"] = br
print(Branching_ratios)
        
