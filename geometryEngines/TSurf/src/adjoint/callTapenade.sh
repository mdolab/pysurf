# Preprocess some files
cpp -traditional -P -DUSE_TAPENADE ../common/precision.F90 precision.f90
cpp -traditional -P -DUSE_TAPENADE ../common/constants.F90 constants.f90

# INTERSECTION AND CURVE DERIVATIVES

# Generate forward mode
tapenade -d -head "triTriIntersect remesh_main barProjection(x1,x2,x)\(xf)" ../intersections/Intersection.F90 ../utilities/Utilities.F90 ../curveSearch/curveUtils.F90

# Generate backward mode
tapenade -b -head "triTriIntersect remesh_main barProjection(x1,x2,x)\(xf)" ../intersections/Intersection.F90 ../utilities/Utilities.F90 ../curveSearch/curveUtils.F90

# PROJECTION DERIVATIVES

# Generate forward mode
tapenade -d -head "triaProjection(x1,x2,x3,x)\(xf,u,v) quadProjResidual(x1,x2,x3,x4,x,u,v)\(residual) quadProjOutput(x1,x2,x3,x4,u,v)\(xf) triaWeights quadWeights computeNodalNormals" ../ADT/adtProjections.F90 precision.f90 constants.f90

# Generate backward mode
tapenade -b -head "triaProjection(x1,x2,x3,x)\(xf,u,v) quadProjResidual(x1,x2,x3,x4,x,u,v)\(residual) quadProjOutput(x1,x2,x3,x4,u,v)\(xf) triaWeights quadWeights computeNodalNormals" ../ADT/adtProjections.F90 precision.f90 constants.f90

# Adjust AD code to use original modules that should not be differentiated
sed -i -e 's/PRECISION_D/PRECISION/' -e 's/CONSTANTS_D/CONSTANTS/' adtprojections_d.f90
sed -i -e 's/PRECISION_B/PRECISION/' -e 's/CONSTANTS_B/CONSTANTS/' adtprojections_b.f90

# Remove modules that are not useful, but were differentiated anyway
rm precision_d.*
rm precision_b.*
rm constants_d.*
rm constants_b.*

# Remove processed files
rm precision.f90
rm constants.f90
