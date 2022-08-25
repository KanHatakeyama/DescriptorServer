from ..models import Molecule
from .Autodescriptor import *
import json

rdkit_calculator = RDKitDescriptors()
avfp_calculator = Fingerprint()
mord2d_calculator = MordredDescriptor()
jr_calculator = GroupContMethod()


def process_smiles(smiles_list, option_list):
    smiles_list = smiles_list.split("\r\n")

    res_dict = {}

    for smiles in smiles_list:
        res_dict[smiles] = fetch_descriptor(smiles, option_list)

    return res_dict


def fetch_descriptor(smiles, option_list):

    data_dict = {}

    # find smiles from SQL
    try:
        obj = Molecule.objects.get(SMILES__exact=smiles)
    except:
        # if not found, make new record
        obj = Molecule(SMILES=smiles)
        obj.save()

    data_dict["SMILES"] = smiles

    # descriptors
    # if available, use database values. if not, calculate them

    command_list = [
        ["rdkit", obj.RDKit_desc, rdkit_calculator],
        ["avalonFP", obj.avfp_desc, avfp_calculator],
        ["mordred2d", obj.mord2d_desc, mord2d_calculator],
        ["JR", obj.jr_desc, jr_calculator],
    ]
    for command in command_list:
        name = command[0]
        desc = command[1]
        calculator = command[2]

        if name in option_list:
            if desc is None:
                desc = calculator.calc(smiles)
                obj.save()
            data_dict[name] = desc

    return data_dict
