from .GaussianWrapper import auto_gaussian, random_name
from .parse_result import parse_result
import os
import glob

# haatree to wavelength
# https://www2.chemistry.msu.edu/faculty/reusch/virttxtjml/cnvcalc.htm


# gaussian header
header = """%CPU=0-7
%Mem=30GB
%chk=test123456.chk
#opt pm7

structure opt pm7"""

# footer
footer = """--Link1--
%CPU=0-7
%Mem=30GB
%chk=test123456.chk
#p pm7 maxdisk=100GB geom=Check Polar CPHF=RdFreq

bunkyoku

0 1

 0.1302 0.0937 0.0773 0.0694 0.0570

"""


class GaussianPM7:
    def __init__(self):
        pass

    def calc(self, sm):
        return calc_pm7_by_gaussian(sm)


def calc_pm7_by_gaussian(sm):
    project_name = random_name(30)
    global header, footer
    header = header.replace("test123456", project_name)
    footer = footer.replace("test123456", project_name)

    sm = sm.replace("*", "H")
    # calc
    try:
        res = auto_gaussian(sm, header, footer, project_name=project_name)
        parsed_res = parse_result(res)
    except:
        parsed_res = {}

    # delete temp files
    for path in [project_name+".chk",
                 f"temp/{project_name}.com",
                 f"temp/{project_name}.inp",
                 f"temp/{project_name}.log", ]:
        if os.path.exists(path):
            os.remove(path)

    # remove other caches
    for path in glob.glob("*.chk"):
        try:
            os.remove(path)
        except:
            pass

    return parsed_res
