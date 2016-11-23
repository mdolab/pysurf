# Rescale primary meshes
cgns_utils scale primary_meshes/fuse_L1_hyp.cgns 0.0254 fuse_L1_temp.cgns
cgns_utils scale primary_meshes/wing_L1_hyp.cgns 0.0254 wing_L1_temp.cgns

# Coarsen collar grid
#cgns_utils coarsen collar.cgns collar_L1_temp.cgns
cp collar.cgns collar_L1_temp.cgns

# Combine meshes in a single CGNS file
cgns_utils combine fuse_L1_temp.cgns wing_L1_temp.cgns collar_L1_temp.cgns crm_wb.cgns

# Remove temporary meshes
rm fuse_L1_temp.cgns
rm wing_L1_temp.cgns
rm collar_L1_temp.cgns

# Now create the background mesh
# positional arguments:
#  gridFile    Name of input CGNS file
#  dh          Uniform cartesian spacing size
#  hExtra      Extension in "O" dimension
#  nExtra      Number of nodes to use for extension
#  sym         Normal for possible sym plane
#  mgcycle     Minimum MG cycle to enforce
#  outFile     Name of output CGNS file
#cgns_utils simpleOCart crm_wb.cgns 0.6 1500.0 65 y 3 background.cgns

# Now combine everything in a single file
cgns_utils combine crm_wb.cgns background.cgns crm_wb.cgns

# First we use cgns_utils to set all symmetry plane coordinates to a hard zero
cgns_utils symmzero crm_wb.cgns y crm_wb.cgns

# Use cgns_utils to merge contiguous blocks
cgns_utils merge crm_wb.cgns

# Check block-to-block connectivities
cgns_utils connect crm_wb.cgns

# Coarsen collar grid
cgns_utils coarsen crm_wb.cgns crm_wb_L2.cgns
