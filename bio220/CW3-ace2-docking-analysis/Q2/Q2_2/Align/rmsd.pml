############################################
# 1. LOAD DATA
############################################
reinitialize everything
load 1R4L_protein_zn.dock4, rs1
load ../../Q2_3/PLIP Input/result.dock4, rs2
load 1R4L_protein_zn_water_cl.dock4, rs3
load 1R4L_protein_zn_water_cl_nag800_nag801.dock4, rs4
load ../../Q2_1/1R4L.pdb, 1R4L
create ref, 1R4L and resn XX5

############################################
# 2. PYTHON LOGIC: FIND LOWEST RMSD
############################################
python
import pymol.cmd as cmd

def find_best_docking_pose(target_objs, ref_obj):
    print("\n" + "="*80)
    print(f"{'Object':<10} | {'Total':<6} | {'Best':<6} | {'RMSD (Heavy)':<12} | {'Status'}")
    print("="*80)

    for obj in target_objs:
        n_states = cmd.count_states(obj)
        best_rmsd = 9999.0
        best_state = -1
        
        for i in range(1, n_states + 1):
            mobile_sel = f"{obj} and state {i} and not element H"
            target_sel = f"{ref_obj}"
            
            try:
                results = cmd.align(mobile_sel, target_sel, cycles=0, transform=0)
                
                curr_rmsd = results[0]
                atom_count = results[1]

                if atom_count < 5: 
                    continue

                if curr_rmsd < best_rmsd:
                    best_rmsd = curr_rmsd
                    best_state = i
            except:
                continue

        if best_state != -1:
            print(f"{obj:<10} | {n_states:<6} | {best_state:<6} | {best_rmsd:.3f}        | OK")
            new_name = f"{obj}_best"
            cmd.create(new_name, obj, best_state, 1)
        else:
            print(f"{obj:<10} | {n_states:<6} | {'--':<6} | {'FAIL':<12} | Mismatch")

docking_results = ["rs1", "rs2", "rs3", "rs4"]
find_best_docking_pose(docking_results, "ref")
python end

############################################
# 3. VISUALIZATION & ORGANIZATION
############################################
# Group the original raw data and hide it
group raw_data, rs1 rs2 rs3 rs4
disable raw_data

# Group the best poses
group best_poses, rs1_best rs2_best rs3_best rs4_best 1R4L ref

# Set up a nice view for comparison
show sticks, ref
show sticks, best_poses

# Color coding for differentiation
color green, ref
color cyan, rs1_best
color magenta, rs2_best
color yellow, rs3_best
color salmon, rs4_best

# Add labels to identify them
set label_size, 20
label ref and name C1, "Crystal Ref"
label rs1_best and name C1, "RS1 Best"
label rs2_best and name C1, "RS2 Best"
label rs3_best and name C1, "RS3 Best"
label rs4_best and name C1, "RS4 Best"

# Center view
zoom ref