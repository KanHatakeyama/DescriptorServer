from rdkit import Chem
from rdkit.Chem import AllChem
import re
import joblib
import openbabel as ob
from openbabel import pybel
import os
import hashlib
import random
import string


def check_normal_finish(gauss_log_str):
    flg = gauss_log_str.find("Normal termination of Gaussian 16")

    if flg >= 0:
        return True
    else:
        return False


def load_txt(path):
    with open(path) as f:
        log_string = f.read()
    return log_string


# use openbabel
def calc_geom(sm, project_name):
    g_temp_path = "temp/"+project_name+".inp"
    m = pybel.readstring("smiles", sm)
    m.localopt()
    m.write('gjf', g_temp_path, overwrite=True)
    geom_str = load_txt(g_temp_path)
    geom_str = "\n".join(geom_str.split("\n")[3:])
    # print(geom_str)
    return geom_str


def get_hash(sm):
    return hashlib.md5(sm.encode()).hexdigest()


def auto_gaussian(sm, header, footer,
                  project_name="proj"):
    # calc geometry roughly
    geom_str = calc_geom(sm, project_name)

    # run gaussian
    run_gaussian(project_name, header, geom_str, footer)

    # load result
    gauss_log_str = load_txt("temp/"+project_name+".log")

    fin_flg = check_normal_finish(gauss_log_str)

    res_dict = {}

    if fin_flg:
        alpha_dict = parse_alpha(gauss_log_str)
        res_dict["alpha"] = alpha_dict

    res_dict["SMILES"] = sm
    res_dict["log"] = gauss_log_str
    res_dict["fin"] = fin_flg
    #hs = get_hash(sm)
    #joblib.dump(res_dict, folder_path+"/"+hs+".bin", compress=9)
    return res_dict

    return fin_flg


def run_gaussian(project_name, header, geom_str, footer):
    out_string = header+geom_str+footer
    out_path = "temp/"+project_name+".com"
    with open(out_path, mode='w') as f:
        f.write(out_string)

    # run gaussian
    command = "g16 "+out_path
    os.system(command)

# parse*************


def parse_volumme(gauss_log_str):
    pattern = "Molar volume = (.*?)bohr"
    result = re.search(pattern, gauss_log_str, re.S)
    vol = result.group()
    vol = vol.replace(" ", "")
    for c in ["Molar", "volume", "=", "bohr"]:
        vol = vol.replace(c, "")
    return float(vol)


def parse_alpha(gauss_log_str):
    alpha_dict = {}
    pattern = " Isotropic polarizability (.*?)Bohr"
    for result in re.finditer(pattern, gauss_log_str, re.S):
        al = result.group()
        al = al.replace("Isotropic polarizability for W=    ", "")
        al = al.replace(" Bohr", "")
        temp = al.split(" ")
        alpha_dict[float(temp[1])] = (float(temp[-1]))

    return alpha_dict


def parse_energy(gauss_log_str):
    pattern = "SCF Done:(.*?)A.U."

    result = list(re.finditer(pattern, gauss_log_str, re.S))[-1]
    string = result.group()
    string = string.replace("SCF Done:  E(RB3LYP) =", "")
    string = string.replace("A.U.", "")
    string = string.replace(" ", "")

    return float(string)


def parse_homo_lumo(gauss_log_str):
    pattern = "Alpha  occ. eigenvalues"

    # find pos of occupied orbital energy
    pos = gauss_log_str.rfind(pattern)
    homo_line, lumo_line = gauss_log_str[pos:pos+1000].split("\n")[:2]

    line = homo_line

    def get_numbers(line):
        arrays = line.split(" ")
        num_list = []
        for v in arrays:
            try:
                num_list.append(float(v))
            except:
                pass
        return num_list

    homo_energies = get_numbers(homo_line)
    lumo_energies = get_numbers(lumo_line)

    homo = homo_energies[-1]
    lumo = lumo_energies[0]

    return homo, lumo


def parse_dipole(gauss_log_str):
    pattern = "Dipole moment"

    # find pos of occupied orbital energy
    pos = gauss_log_str.rfind(pattern)
    line = gauss_log_str[pos:pos+1000].split("\n")[1]

    di_X = re.search(r'X=\s+[-\d]+.[-\d]+', line).group()
    di_Y = re.search(r'Y=\s+[-\d]+.[-\d]+', line).group()
    di_Z = re.search(r'Z=\s+[-\d]+.[-\d]+', line).group()
    di_tot = re.search(r'Tot=\s+[-\d]+.[-\d]+', line).group()

    def parse(di):
        di = di[5:].replace(" ", "")
        return float(di)

    return parse(di_X), parse(di_Y), parse(di_Z), parse(di_tot)


def random_name(n):
    randlst = [random.choice(string.ascii_letters + string.digits)
               for i in range(n)]
    return ''.join(randlst)
