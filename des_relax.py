from pyrosetta import *
init('-in:auto_setup_metals  -constraints:cst_fa_file b-roll.cst')
from pyrosetta.rosetta.core.select import residue_selector as selection
from pyrosetta.rosetta.core.pack.task import operation
from pyrosetta.rosetta.protocols import minimization_packing as pack_min
from pyrosetta.rosetta.core.select.movemap import *
from pyrosetta.rosetta.protocols.relax import FastRelax

def des_relax_pose(pose):
## set up score function and reverse fold tree
    scorefunction = create_score_function("ref2015_cst")
    pyrosetta.rosetta.core.scoring.constraints.add_fa_constraints_from_cmdline_to_pose(pose)
    ft = pose.fold_tree() 
    ft.reorder(pose.total_residue())
    pose.fold_tree(ft)
    scorefunction(pose)
    
    ## added is residues different from the wt b roll
    added =  pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector()
    added.set_index("10,11,21,22,32,33,43,44,54,65,66,77,78,79,80,81,83,84,85")

    # select interaction partners within a 5A radius
    nbr_selector = selection.NeighborhoodResidueSelector()
    nbr_selector.set_focus_selector(added)
    nbr_selector.set_include_focus_in_subset(False)
    nbr_selector2 = selection.NeighborhoodResidueSelector()
    nbr_selector2.set_focus_selector(added)
    nbr_selector2.set_include_focus_in_subset(True)
   
    # tell Rosetta which residues to design and repack
    tf = pyrosetta.rosetta.core.pack.task.TaskFactory()
    tf.clear()
    tf.push_back(operation.InitializeFromCommandline())
    prevent_repacking_rlt = operation.PreventRepackingRLT()
    prevent_subset_repacking = operation.OperateOnResidueSubset(prevent_repacking_rlt, nbr_selector2, True )
    tf.push_back(prevent_subset_repacking)
    restrict_to_repack = operation.RestrictToRepackingRLT()
    prevent_nbr_design = operation.OperateOnResidueSubset(restrict_to_repack, nbr_selector, False )
    tf.push_back(prevent_nbr_design)
    aa_to_design = pyrosetta.rosetta.core.pack.task.operation.RestrictAbsentCanonicalAASRLT()
    aa_to_design.aas_to_keep("ACDEFGHIKLMNPQRSTVWY")
    aa_design = operation.OperateOnResidueSubset(aa_to_design, added, True )
    
    packer_task = tf.create_task_and_apply_taskoperations(pose)
    print(packer_task)
    
    # Create relax move map and set up fast relax
    mmf = MoveMapFactory()
    mmf.add_bb_action(mm_enable, added)
    mmf.add_chi_action(mm_enable, nbr_selector2)

    fr = FastRelax()
    fr.set_scorefxn(scorefunction)
    fr.set_movemap_factory(mmf)
    fr.set_task_factory(tf)
    fr.constrain_relax_to_start_coords(True)
    
    rel_pose = fr.apply(pose)
    return rel_pose
    
def main(pdb_file):
    pose = pose_from_pdb(pdb_file)
    jd = PyJobDistributor("outputs/mut_des_beta_roll", 1, scorefunction)
    jd.native_pose = pose

    rel_pose = Pose()
    while not jd.job_complete:
        rel_pose.assign(pose)
        if not os.getenv("DEBUG"):
            des_relax_pose(rel_pose)
        jd.output_decoy(rel_pose)
