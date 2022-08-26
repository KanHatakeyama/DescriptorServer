from .GaussianWrapper import parse_alpha, parse_homo_lumo, parse_energy, parse_dipole
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdDistGeom import EmbedMolecule
from rdkit.Chem import Descriptors

conversion_const = 0.148185


def calc_n(alpha, volume):
    conv_alpha = alpha*conversion_const
    # volume=volume*conversion_const
    const_K = 4*np.pi/3*conv_alpha/volume
    eps = (1+2*const_K)/(1-const_K)
    n = eps**0.5
    if type(n) == type(1j):
        return np.nan
    return n


def calc_molecule(sm):
    m = Chem.MolFromSmiles(sm)
    try:
        mw = Descriptors.MolWt(m)
    except:
        mw = np.nan
    try:
        m2 = Chem.AddHs(m)
        params = Chem.rdDistGeom.ETKDGv2()
        EmbedMolecule(m2, params)
        AllChem.MMFFOptimizeMolecule(m2)
        vol = AllChem.ComputeMolVolume(m2)
    except:
        vol = np.nan
    return vol, mw


def parse_result(dat):
    res_dict = {}

    gauss_log_str = dat["log"]
    alpha_dict = parse_alpha(gauss_log_str)
    sm = dat["SMILES"]

    vol, mw = calc_molecule(sm)

    try:
        alpha = list(alpha_dict.values())[0]
    except:
        alpha = np.nan

    try:
        res_dict["DFT_energy"] = parse_energy(gauss_log_str)
    except:
        pass

    try:
        res_dict["dipoleX"], res_dict["dipoleY"], res_dict["dipoleZ"], res_dict["dipoleTot"] = parse_dipole(
            gauss_log_str)
    except:
        res_dict["dipoleX"], res_dict["dipoleY"], res_dict["dipoleZ"], res_dict["dipoleTot"] = np.nan, np.nan, np.nan, np.nan

    try:
        res_dict["HOMO"], res_dict["LUMO"] = parse_homo_lumo(gauss_log_str)
    except:
        res_dict["HOMO"], res_dict["LUMO"] = np.nan, np.nan

    #res_dict["SMILES"] = sm
    res_dict["est_RI"] = calc_n(alpha, vol)
    # res_dict["alpha"]=alpha
    res_dict["vol"] = vol
    res_dict["est_dens"] = calc_density(mw, vol)

    for k, v in alpha_dict.items():
        res_dict["alpha_"+str(k)] = v

    # exclude nan
    final_dict = {}
    for k, v in res_dict.items():
        if v == v:
            final_dict[k] = v

    return final_dict


def calc_density(mw, volume):
    Na = 6.02
    return mw/volume/Na*10
