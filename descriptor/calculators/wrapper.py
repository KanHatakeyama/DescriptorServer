from ..models import Molecule
from .Autodescriptor import *
from .gaussian.django_wrapper import GaussianPM7
import json
import random

rdkit_calculator = RDKitDescriptors()
avfp_calculator = Fingerprint()
mord2d_calculator = MordredDescriptor()
jr_calculator = GroupContMethod()
pm7_calculator = GaussianPM7()


def process_smiles(smiles_list, option_list):
    smiles_list = smiles_list.split("\r\n")

    res_dict = {}

    # register smiles
    for smiles in smiles_list:
        # find smiles from SQL
        try:
            obj = Molecule.objects.get(SMILES__exact=smiles)
            obj.done = False
        except:
            # if not found, make new record
            obj = Molecule(SMILES=smiles)
        obj.save()

    random.shuffle(smiles_list)

    # calculate
    for smiles in smiles_list:
        res_dict[smiles] = fetch_descriptor(smiles, option_list)

    return res_dict


def fetch_descriptor(smiles, option_list):

    data_dict = {}
    data_dict["SMILES"] = smiles

    obj = Molecule.objects.get(SMILES__exact=smiles)

    # get descriptors
    # if available, use database values. if not, calculate them

    command_list = [
        ["rdkit", "RDKit_desc", rdkit_calculator],
        ["avalonFP", "avfp_desc", avfp_calculator],
        ["mordred2d", "mord2d_desc", mord2d_calculator],
        ["JR", "jr_desc", jr_calculator],
        ["PM7", "PM7_desc", pm7_calculator],
    ]
    for command in command_list:
        name = command[0]
        desc = command[1]
        calculator = command[2]

        if name in option_list:
            if getattr(obj, desc) is None:
                dict_desc = calculator.calc(smiles)
                json_desc = json.dumps(dict_desc)
                setattr(obj, desc, json_desc)
                obj.save()
            else:
                dict_desc = json.loads(getattr(obj, desc))

            for k, v in dict_desc.items():
                data_dict[f"{name}_{k}"] = v

    obj.done = True
    obj.save()
    return data_dict
