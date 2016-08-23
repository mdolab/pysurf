The main program that should be called in order to generate a mesh is
"computeMesh" defined in tfiMeshTools.py

Even though TFI lib has support for cubic interpolation, the derivatives
estimatimation is not working quite well for now.
So we recomend using linear interpolation only. This is done by setting
derivDict to None when calling computeMesh.
