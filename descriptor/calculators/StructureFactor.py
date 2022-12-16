from .Autodescriptor import RDKitDescriptors, mol_from_smiles
from sklearn.decomposition import PCA
import numpy as np
from rdkit.Chem import AllChem


class StructureFactor(RDKitDescriptors):
    def __init__(self, dict_mode=True):
        super(RDKitDescriptors, self).__init__()
        self.calculator = calc_factors
        self.desc_list = ["volume", "eps1", "eps2"]
        self.auto_correct = False
        self.dict_mode = dict_mode

    def _desc_calc(self, m):
        return self.calculator(m)


def calc_factors(smiles):
    m = mol_from_smiles(smiles)
    try:
        eps1, eps2 = calc_structure_factors(m)
    except Exception as e:
        print(e)
        eps1, eps2 = 1, 1

    try:
        vol = AllChem.ComputeMolVolume(m)
    except:
        vol = np.nan
    return [vol, eps1, eps2]


def calc_3d_pos(m):
    # AllChem.Compute2DCoords(m)
    AllChem.EmbedMolecule(m, randomSeed=0xf00d)
    pos_array = []
    for c in m.GetConformers():
        pos_array.append(c.GetPositions())

    pos_array = np.array(pos_array)
    pos_array = pos_array.reshape(-1, 3)
    return pos_array


def calc_structure_factors(m):
    pos_array = calc_3d_pos(m)

    # calc PCA pos and eps
    pca = PCA(n_components=3)
    pca_pos = pca.fit_transform(pos_array)

    var_array = []
    for dim in range(3):
        var_array.append(np.var(pca_pos[:, dim]))

    var_array = sorted(var_array, reverse=True)
    eps1 = var_array[1]/var_array[0]
    eps2 = var_array[2]/var_array[0]
    return eps1, eps2
