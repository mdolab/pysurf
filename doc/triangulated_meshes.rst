.. _triangulated_meshes:

Generating triangulated meshes in ICEM
======================================

#. In order to start the triangulated mesh generation, you will need a valid .tin file. You can import an IGES or STEP file in ICEM in order to generate a .tin file.

#. Open the .tin file in ICEM.

#. Check if you have a smooth surface, without any overlapping patches. You might need to adjust/heal the surface if this is case.

#. You might want to assign different parts to split the surface in important segments. For instance, for a wing, you can create separate parts for the upper skin, lower skin, trailing edge, and tip.

#. Then delete all points and curves from the geometry. You should use the appropriate function on the top bar to remove points, and the other one to remove curves. Remember to also delete all dormant points and curves.

#. Then go to Geometry > Repair geometry > Build Diagnostic Topology. Set the tolerance as half the length of the smallest edge of your geometry (e.g. trailing edge). Activate "Single Curve Cleanup" and use the same tolerance. Also activate "Split surface at T-connections", "Join edge curves", "Delete unnatached curves and points", "Keep dormant as dormant", and "Use Local Tolerance". Then select Apply.

#. In the end, the free edges should be yellow, and the interior edges should be red. If you get yellow curves in the interior of the surfaces, then this means that you have a discontinuity in your geometry. You might need to adjust the tolerances or find a way to reduce the gap.

#. Once the topology is correct, you might want to create separate parts for feature curves you find important, such as lines that follows sharp corners. In order to do that, turn off all surfaces and points under the "Geometry" branch of the tree on the left side of the screen. Now you only have the topology curves on your display. Then you can create parts (Right click on the "Parts" branch of the tree and select "Create Part") by selecting the curves you find important. I usually create one part for the upper trailing edge curve, lower trailing edge curve, and leading edge curve. Then ICEM will create bar elements for these curves, which can be used as guide curves during the hyperbolic surface marching procedure.

#. Then go to Mesh > Global Mesh Setup > Global Mesh Size. Select an appropriate size for the Max element (You can adjust this later if the cells are too big). Click Apply.

#. Go to Mesh > Global Mesh Setup > Shell Meshing Parameters. Set "Ignore Size" to the same tolerance you used for the topology, or even smaller. This prevents the collapsing of the trailing edge.

#. Go to Mesh > Compute Mesh > Surface Mesh Only. Select "Overwrite Surface Preset/Default Mesh Type" and choose "All Tri". Select "Overwrite Surface Preset/Default Mesh Method" and choose "Patch Dependent". Then select "Compute". If the resulting mesh is too coarse, you can adjust the Max element size explained two steps ago. But remember that you can also locally refine the mesh, as explained in the next step.

#. If some surfaces become flat (that is, if you notice that no elements are projected into a specific surface), you might need to split it. Go to Geometry > Repair Geometry > Split folded surfaces, and then select the surface that you are having problems with. Then regenerate the topology and repeat all steps from there.

#. If you want to refine a specific region of the surface (e.g. leading edge) after the initial mesh has been generated, you can go to Edit Mesh > Adjust Mesh Density > Refine Selected Mesh. Then you can select the cells you want to refine!

#. After you are happy with your mesh, then you should check if the normals are consistent. The normals of all elements should be pointing to the exterior for the surface. In order to check that you can use a right on click on Mesh > Shells (on the model tree window on the left side of the screen), and then choose "Normals using arrow". Find an element whose arrow is pointing out. Then select Edit Mesh > ReOrient Mesh > Reorient Consistent. Then select the element you found previously. Now all normals should be consistent.

#. Now it is time to export the mesh. Go to Output > Select Solver. Select "CGNS" and "NASTRAN" as options. Activate "Set as Default". Then click Apply.

#. Go to Output > Write Input. Save the project. Select the appropriate .uns file. Then the options should be: "Unstructured", "No", "Face elements", "Yes", "2.4". Then click Done. A new CGNS file should appear in the directory. It can be used as an input to pySurf.
