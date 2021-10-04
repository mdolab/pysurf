# Use f2py to generate preliminary signature file
f2py cgnsAPI.F90 -m cgnsAPI -h cgnsAPI.pyf --overwrite-signature

# Remove allocatable variables from the signature
sed -i '/allocatable/d' cgnsAPI.pyf

# Move signature file to f2py folder
mv cgnsAPI.pyf ../python/f2py/cgnsAPI.pyf
