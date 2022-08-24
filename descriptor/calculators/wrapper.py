from ..models import Molecule
from .Autodescriptor import *

rdkit_calculator = RDKitDescriptors()
avfp_calculator = Fingerprint()
mord2d_calculator = MordredDescriptor()
jr_calculator = GroupContMethod()


def fetch_descriptor(request):

    smiles = (request.GET["SMILES"])
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

    # rdkit
    if "rdkit" in request.GET:
        if obj.RDKit_desc is None:
            obj.RDKit_desc = rdkit_calculator.calc(smiles)
            obj.save()
        data_dict["rdkit"] = obj.RDKit_desc

    # avalon fp
    if "avfp" in request.GET:
        if obj.avfp_desc is None:
            obj.avfp_desc = avfp_calculator.calc(smiles)
            obj.save()
        data_dict["avalonFP"] = obj.avfp_desc

    # mord2d
    if "mord2d" in request.GET:
        if obj.mord2d_desc is None:
            obj.mord2d_desc = mord2d_calculator.calc(smiles)
            obj.save()
        data_dict["Mordred 2D"] = obj.mord2d_desc

    # JR
    if "jr" in request.GET:
        if obj.jr_desc is None:
            obj.jr_desc = jr_calculator.calc(smiles)
            obj.save()
        data_dict["JR"] = obj.jr_desc

    # refresh record
    return data_dict
