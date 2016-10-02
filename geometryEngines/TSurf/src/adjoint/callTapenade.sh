# Generate forward mode
tapenade -d -head "triTriIntersect remesh_main" ../intersections/Intersection.F90 ../utilities/Utilities.F90 ../ADT/adtSearch.F90

# Generate backward mode
tapenade -b -head "triTriIntersect remesh_main" ../intersections/Intersection.F90 ../utilities/Utilities.F90 ../ADT/adtSearch.F90
