import _mod
import f90wrap.runtime
import logging


class Hypsurfapi(f90wrap.runtime.FortranModule):
    """
    Module hypsurfapi


    Defined at hypsurfAPI.F90 lines 18-1028

    """

    @staticmethod
    def computematrices(r0, n0, s0, rm1, sm1, layerindex, theta, sigmasplay, bc1, bc2, numlayers, epse0, f, numnodes):
        """
        computematrices(r0, n0, s0, rm1, sm1, layerindex, theta, sigmasplay, bc1, bc2, \
            numlayers, epse0, f, numnodes)
        
        
        Defined at hypsurfAPI.F90 lines 30-210
        
        Parameters
        ----------
        r0 : float array
        n0 : float array
        s0 : float array
        rm1 : float array
        sm1 : float array
        layerindex : int
        theta : float
        sigmasplay : float
        bc1 : float
        bc2 : float
        numlayers : int
        epse0 : float
        f : float array
        numnodes : int
        
        """
        _mod.f90wrap_computematrices(
            r0=r0,
            n0=n0,
            s0=s0,
            rm1=rm1,
            sm1=sm1,
            layerindex=layerindex,
            theta=theta,
            sigmasplay=sigmasplay,
            bc1=bc1,
            bc2=bc2,
            numlayers=numlayers,
            epse0=epse0,
            f=f,
            numnodes=numnodes,
        )

    @staticmethod
    def matrixbuilder(curr_index, bc1, bc2, r0, rm1, n0, s0, sm1, numlayers, epse0, layerindex, theta, numnodes, k, f):
        """
        matrixbuilder(curr_index, bc1, bc2, r0, rm1, n0, s0, sm1, numlayers, epse0, \
            layerindex, theta, numnodes, k, f)
        
        
        Defined at hypsurfAPI.F90 lines 212-406
        
        Parameters
        ----------
        curr_index : int
        bc1 : float
        bc2 : float
        r0 : float array
        rm1 : float array
        n0 : float array
        s0 : float array
        sm1 : float array
        numlayers : int
        epse0 : float
        layerindex : int
        theta : float
        numnodes : int
        k : float array
        f : float array
        
        """
        _mod.f90wrap_matrixbuilder(
            curr_index=curr_index,
            bc1=bc1,
            bc2=bc2,
            r0=r0,
            rm1=rm1,
            n0=n0,
            s0=s0,
            sm1=sm1,
            numlayers=numlayers,
            epse0=epse0,
            layerindex=layerindex,
            theta=theta,
            numnodes=numnodes,
            k=k,
            f=f,
        )

    @staticmethod
    def dissipationcoefficients(layerindex, r0_xi, r0_eta, dsensor, angle, numlayers, epse0):
        """
        epse, epsi = dissipationcoefficients(layerindex, r0_xi, r0_eta, dsensor, angle, \
            numlayers, epse0)
        
        
        Defined at hypsurfAPI.F90 lines 408-454
        
        Parameters
        ----------
        layerindex : int
        r0_xi : float array
        r0_eta : float array
        dsensor : float
        angle : float
        numlayers : int
        epse0 : float
        
        Returns
        -------
        epse : float
        epsi : float
        
        """
        epse, epsi = _mod.f90wrap_dissipationcoefficients(
            layerindex=layerindex,
            r0_xi=r0_xi,
            r0_eta=r0_eta,
            dsensor=dsensor,
            angle=angle,
            numlayers=numlayers,
            epse0=epse0,
        )
        return epse, epsi

    @staticmethod
    def areafactor(r0, d, nuarea, numareapasses, bc1, bc2, n, s):
        """
        maxstretch = areafactor(r0, d, nuarea, numareapasses, bc1, bc2, n, s)


        Defined at hypsurfAPI.F90 lines 456-522

        Parameters
        ----------
        r0 : float array
        d : float
        nuarea : float
        numareapasses : int
        bc1 : float
        bc2 : float
        n : int
        s : float array

        Returns
        -------
        maxstretch : float

        """
        maxstretch = _mod.f90wrap_areafactor(
            r0=r0, d=d, nuarea=nuarea, numareapasses=numareapasses, bc1=bc1, bc2=bc2, n=n, s=s
        )
        return maxstretch

    @staticmethod
    def smoothing(r, eta, alphap0, numsmoothingpasses, numlayers, n):
        """
        smoothing(r, eta, alphap0, numsmoothingpasses, numlayers, n)


        Defined at hypsurfAPI.F90 lines 524-571

        Parameters
        ----------
        r : float array
        eta : float
        alphap0 : float
        numsmoothingpasses : int
        numlayers : int
        n : int

        """
        _mod.f90wrap_smoothing(
            r=r, eta=eta, alphap0=alphap0, numsmoothingpasses=numsmoothingpasses, numlayers=numlayers, n=n
        )

    @staticmethod
    def qualitycheck(r, numlayers, numnodes, ratios, layerindex=None):
        """
        fail = qualitycheck(r, numlayers, numnodes, ratios[, layerindex])


        Defined at hypsurfAPI.F90 lines 573-699

        Parameters
        ----------
        r : float array
        numlayers : int
        numnodes : int
        ratios : float array
        layerindex : int

        Returns
        -------
        fail : int

        """
        fail = _mod.f90wrap_qualitycheck(
            r=r, numlayers=numlayers, numnodes=numnodes, ratios=ratios, layerindex=layerindex
        )
        return fail

    @staticmethod
    def march(
        py_projection,
        rstart,
        dstart,
        theta,
        sigmasplay,
        bc1,
        bc2,
        plotquality,
        epse0,
        alphap0,
        extension,
        nuarea,
        ratioguess,
        cmax,
        numsmoothingpasses,
        numareapasses,
        numlayers,
        numnodes,
        r,
        ratios,
    ):
        """
        fail = march(py_projection, rstart, dstart, theta, sigmasplay, bc1, bc2, \
            plotquality, epse0, alphap0, extension, nuarea, ratioguess, cmax, \
            numsmoothingpasses, numareapasses, numlayers, numnodes, r, ratios)
        
        
        Defined at hypsurfAPI.F90 lines 701-851
        
        Parameters
        ----------
        py_projection : float
        rstart : float array
        dstart : float
        theta : float
        sigmasplay : float
        bc1 : float
        bc2 : float
        plotquality : int
        epse0 : float
        alphap0 : float
        extension : float
        nuarea : float
        ratioguess : float
        cmax : float
        numsmoothingpasses : int
        numareapasses : int
        numlayers : int
        numnodes : int
        r : float array
        ratios : float array
        
        Returns
        -------
        fail : int
        
        ===========================================================
         Some functions require the area factors of the first-before-last curve
         We will repeat the first curve areas for simplicity.
         rNext, NNext, rm1 for the first iteration are computed at the beginning of the \
             function.
         But we still need to find Sm1
        """
        fail = _mod.f90wrap_march(
            py_projection=py_projection,
            rstart=rstart,
            dstart=dstart,
            theta=theta,
            sigmasplay=sigmasplay,
            bc1=bc1,
            bc2=bc2,
            plotquality=plotquality,
            epse0=epse0,
            alphap0=alphap0,
            extension=extension,
            nuarea=nuarea,
            ratioguess=ratioguess,
            cmax=cmax,
            numsmoothingpasses=numsmoothingpasses,
            numareapasses=numareapasses,
            numlayers=numlayers,
            numnodes=numnodes,
            r=r,
            ratios=ratios,
        )
        return fail

    @staticmethod
    def findradius(r, numnodes):
        """
        radius = findradius(r, numnodes)


        Defined at hypsurfAPI.F90 lines 853-887

        Parameters
        ----------
        r : float array
        numnodes : int

        Returns
        -------
        radius : float

        """
        radius = _mod.f90wrap_findradius(r=r, numnodes=numnodes)
        return radius

    @staticmethod
    def findratio(dmax, d0, numlayers, ratioguess):
        """
        q = findratio(dmax, d0, numlayers, ratioguess)


        Defined at hypsurfAPI.F90 lines 889-927

        Parameters
        ----------
        dmax : float
        d0 : float
        numlayers : int
        ratioguess : float

        Returns
        -------
        q : float

        """
        q = _mod.f90wrap_findratio(dmax=dmax, d0=d0, numlayers=numlayers, ratioguess=ratioguess)
        return q

    @staticmethod
    def matinv3(a, b):
        """
        matinv3(a, b)


        Defined at hypsurfAPI.F90 lines 929-957

        Parameters
        ----------
        a : float array
        b : float array

        """
        _mod.f90wrap_matinv3(a=a, b=b)

    @staticmethod
    def giveangle(r0, r1, r2, n1, angle):
        """
        giveangle(r0, r1, r2, n1, angle)


        Defined at hypsurfAPI.F90 lines 959-997

        Parameters
        ----------
        r0 : float array
        r1 : float array
        r2 : float array
        n1 : float array
        angle : float

        """
        _mod.f90wrap_giveangle(r0=r0, r1=r1, r2=r2, n1=n1, angle=angle)

    @staticmethod
    def cross(a, b, c):
        """
        cross(a, b, c)


        Defined at hypsurfAPI.F90 lines 999-1011

        Parameters
        ----------
        a : float array
        b : float array
        c : float array

        """
        _mod.f90wrap_cross(a=a, b=b, c=c)

    @staticmethod
    def m33det(a):
        """
        det = m33det(a)


        Defined at hypsurfAPI.F90 lines 1013-1028

        Parameters
        ----------
        a : float array

        Returns
        -------
        det : float

        """
        det = _mod.f90wrap_m33det(a=a)
        return det

    _dt_array_initialisers = []


hypsurfapi = Hypsurfapi()
