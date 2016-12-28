f2py hypsurfAPI.F90 -m hypsurfAPI -h hypsurfAPI.pyf --overwrite-signature

sed -i '8s/.*/            real(kind=realtype) dimension(3 * numnodes),depend(numnodes),intent(out) :: rnext/' hypsurfAPI.pyf
sed -i '20s/.*/            real(kind=realtype) dimension(3 * numnodes),depend(numnodes),intent(in) :: rnext/' hypsurfAPI.pyf
sed -i '34s/.*/            real(kind=realtype) dimension(3 * numnodes),depend(numnodes),intent(in) :: rnext/' hypsurfAPI.pyf

sed -i '21s/.*/            real(kind=realtype) dimension(3 * numnodes),depend(numnodes),intent(out) :: rnextd/' hypsurfAPI.pyf

sed -i '35s/.*/            real(kind=realtype) dimension(3 * numnodes),depend(numnodes),intent(in) :: rnextb/' hypsurfAPI.pyf

sed -i '9s/.*/            real(kind=realtype) dimension(3,numnodes),depend(numnodes),intent(out) :: nnext/' hypsurfAPI.pyf
sed -i '22s/.*/            real(kind=realtype) dimension(3,numnodes),depend(numnodes),intent(in) :: nnext/' hypsurfAPI.pyf
sed -i '36s/.*/            real(kind=realtype) dimension(3,numnodes),depend(numnodes),intent(in) :: nnext/' hypsurfAPI.pyf

sed -i '23s/.*/            real(kind=realtype) dimension(3,numnodes),depend(numnodes),intent(out) :: nnextd/' hypsurfAPI.pyf

sed -i '37s/.*/            real(kind=realtype) dimension(3,numnodes),depend(numnodes),intent(in) :: nnextb/' hypsurfAPI.pyf

sed -i '10s/.*/            integer(kind=inttype),intent(in) :: storedict/' hypsurfAPI.pyf

sed -i '24s/.*/            integer(kind=inttype),intent(in) :: layerid/' hypsurfAPI.pyf

sed -i '38s/.*/            integer(kind=inttype),intent(in) :: projid/' hypsurfAPI.pyf

sed -i '32s/.*/            real(kind=realtype) dimension(3 * numnodes), intent(in) :: rremeshed/' hypsurfAPI.pyf
sed -i '33s/.*/            real(kind=realtype) dimension(3 * numnodes), intent(out), depend(numnodes) :: rremeshedb/' hypsurfAPI.pyf

cp hypsurfAPI.pyf ../python/f2py/hypsurfAPI.pyf
