# Use f2py to generate preliminary signature file
f2py intersectionAPI.F90 -m intersectionAPI -h intersectionAPI.pyf --overwrite-signature

# Remove allocatable variables from the signature
sed -i '/allocatable/d' intersectionAPI.pyf

# Move signature file to f2py folder
cp intersectionAPI.pyf ../python/f2py/intersectionAPI.pyf
