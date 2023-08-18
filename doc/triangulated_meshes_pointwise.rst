Generating triangulated meshes in Pointwise
===========================================

#. Open the surface mesh in Pointwise. This can ideally be the same Pointwise surface mesh that was used for extrusion by pyHyp.

#. Select the surface mesh (domains) for the surfaces that need to be triangulated.

#. Go to Create > Diagonalize. 

#. Select Diagonalize Method as 'Aligned'. This will align all the triangulated elements in the same way as the structured surface mesh. Keep default check on 'Link' and 'Hide' Visibility/Delete Options.

#. Before exporting, check that the solver attributes are set to File Type: 'adf' and Dimension is set to '2D'.

#. Selected triangulated domains can then be exported by clicking on Export > CAE > cgns file name. Use Data Precision: Double. Other export options can be unchecked.


.. note::

While using this triangulated surface for DVGeoMulti, the terminal will show the number of triangular elements and some details about unspecified boundary conditions on the triangulated surface, which can be ignored.
The structured surface mesh can also be refined before triangulation in specific intersection regions. The resulting triangulated mesh will reflect the same.
