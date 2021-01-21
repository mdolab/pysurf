with open("results.txt", "r") as fid:
    lines = fid.readlines()

# Initialize connectivities and coordinate arrays
coor = []
conn = []

# Gather number of lines and compute the number of tets
numPts = len(lines)
numTets = int(numPts / 4)

# Create all connectivities
for ii in range(numTets):
    # Create connectivity list
    conn.append(" ".join(map(str, [4 * ii + 1, 4 * ii + 2, 4 * ii + 3, 4 * ii + 4, "\n"])))

# Define header of the tecplot file
header = """TITLE = "Intersection parents"
VARIABLES = "CoordinateX", "CoordinateY", "CoordinateZ"
ZONE N = %d, E = %d, DATAPACKING = POINT, ZONETYPE = FETETRAHEDRON\n""" % (
    numPts,
    numTets,
)

# Now add all coordinates and connectivities
lines = [header] + lines + ["\n"] + conn

# Now write everything into a file
with open("intTets.dat", "w") as fid:
    fid.writelines(lines)
