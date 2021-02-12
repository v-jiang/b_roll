## code takes a pdb file, native pdb, and list of mutated residues to generate heatmap of mutations
import pandas as pd
%matplotlib inline
import matplotlib.pyplot as plt 
!pip install seaborn
import seaborn as sns

def score_mutant(filename, native_filename, mut_res):
    native_pose = pose_from_pdb(native_filename)
    scorefunction = create_score_function("ref2015_cst")
    pyrosetta.rosetta.core.scoring.constraints.add_fa_constraints_from_cmdline_to_pose(native_pose)
    scorefunction(native_pose)
    res_dic = {}
    pose = pose_from_pdb(filename)
    scorefunction(pose)
    name = []
    energies = []
    native_res = []
    for res_int in mut_res:
        n =str(res_int)+(native_pose.residue(res_int).name())
        native_res.append(n)
        en = (native_pose.energies().residue_total_energy(res_int))
        res_name = (pose.residue(res_int).name())
        name.append(res_name)
        res_en = pose.energies().residue_total_energy(res_int) - en
        energies.append(res_en)
    res_dic["wt_res"] = native_res
    res_dic['mut_res'] = name
    res_dic['ddG'] = energies
    ddG = pd.DataFrame.from_dict(res_dic)
    ddG_table = pd.pivot_table(ddG, index = "wt_res", columns = "mut_res", values="ddG")
    ax = sns.heatmap(ddG_table)
    return ax
