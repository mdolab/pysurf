#mv wing_L1_hyp.cgns wing_L1_hyp_symbc.cgns
cgns_utils removebc wing_L1_hyp.cgns
cgns_utils overwritebc wing_L1_hyp.cgns wing_bc.bc
cgns_utils fillOpenBCs wing_L1_hyp.cgns bcoverset overset
