# INTERSECTION DERIVATIVES

# Generate forward mode
tapenade -d -head "triTriIntersect remesh_main" ../intersections/Intersection.F90 ../utilities/Utilities.F90 #../ADT/adtSearch.F90

# Generate backward mode
tapenade -b -head "triTriIntersect remesh_main" ../intersections/Intersection.F90 ../utilities/Utilities.F90 #../ADT/adtSearch.F90

# PROJECTION DERIVATIVES

# Run preprocessor
cpp -traditional -P ../ADT/adtUtils.F90 adtUtils.f90

# Generate forward mode
tapenade -d -head "triaProjection(x1,x2,x3,x)\(xf,norm)" adtUtils.f90

# Generate backward mode
tapenade -b -head "triaProjection(x1,x2,x3,x)\(xf,norm)" adtUtils.f90

# Remove preprocessed file
rm adtUtils.f90
