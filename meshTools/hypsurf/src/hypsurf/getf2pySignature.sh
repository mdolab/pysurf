f2py hypsurfAPI.F90 -m hypsurfAPI -h hypsurfAPI.pyf --overwrite-signature
sed -i '8s/.*/            real(kind=realtype) dimension(3 * numnodes),depend(numnodes),intent(out) :: rnext/' hypsurfAPI.pyf
sed -i '18s/.*/            real(kind=realtype) dimension(3 * numnodes),depend(numnodes),intent(out) :: rnext/' hypsurfAPI.pyf
# sed -i '28s/.*/            real(kind=realtype) dimension(3 * numnodes),depend(numnodes),intent(out) :: rnext/' hypsurfAPI.pyf
sed -i '9s/.*/            real(kind=realtype) dimension(3,numnodes),depend(numnodes),intent(out) :: nnext/' hypsurfAPI.pyf
sed -i '19s/.*/            real(kind=realtype) dimension(3,numnodes),depend(numnodes),intent(out) :: nnext/' hypsurfAPI.pyf
# sed -i '29s/.*/            real(kind=realtype) dimension(3,numnodes),depend(numnodes),intent(out) :: nnext/' hypsurfAPI.pyf
cp hypsurfAPI.pyf ../python/f2py/hypsurfAPI.pyf
