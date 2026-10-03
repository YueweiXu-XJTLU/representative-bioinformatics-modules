reinitialize everything
load ../1R4L_protein_zn_water/files/receptor.pdb, receptor
load result.dock4, result
create best_dock, result, 1, 1
select predicted_complex, receptor or best_dock
save predicted_SwissDock.pdb, predicted_complex