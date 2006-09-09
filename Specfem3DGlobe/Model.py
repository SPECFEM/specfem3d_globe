#!/usr/bin/env python


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# base class for all models
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


from cig.addyndum.components import Component


class Model(Component):
    
    # parameters for all models
    import pyre.inventory as pyre
    import cig.addyndum.inventory as addyndum
    from TopoBathy import TopoBathy
    ATTENUATION        = pyre.bool("attenuation")
    ELLIPTICITY        = pyre.bool("ellipticity")
    GRAVITY            = pyre.bool("gravity")
    OCEANS             = pyre.bool("oceans")
    ROTATION           = pyre.bool("rotation")
    TOPOGRAPHY         = pyre.bool("topography")
    topoBathy          = addyndum.facility("topo-bathy", factory=TopoBathy)
    
    # configuration
    def _configure(self):
        Component._configure(self)
        if not (self.TOPOGRAPHY or self.OCEANS):
            # We don't need topo-bathy data.
            self.jettison(Model.topoBathy)
        return


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3D isotropic (S20RTS)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


class Model3DIsotropic(Model):

    # name by which the Fortran code knows this model
    MODEL = "s20rts"

    # names by which the user refers to this model
    componentNames = [ "3D_isotropic", "s20rts" ]

    # additional parameters for this model
    import cig.addyndum.inventory as addyndum
    CNtype2            = addyndum.inputFile("CNtype2",           default="DATA/crust2.0/CNtype2.txt")
    CNtype2_key_modif  = addyndum.inputFile("CNtype2_key_modif", default="DATA/crust2.0/CNtype2_key_modif.txt")
    P12                = addyndum.inputFile("P12",               default="DATA/s20rts/P12.dat")
    S20RTS             = addyndum.inputFile("S20RTS",            default="DATA/s20rts/S20RTS.dat")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3D attenuation (Min Chen)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


class Model3DAttenuation(Model):
    
    # name by which the Fortran code knows this model
    MODEL = "Min_Chen"
    
    # names by which the user refers to this model
    componentNames = [ "3D_attenuation", "Min_Chen" ]
    
    # additional parameters for this model
    import cig.addyndum.inventory as addyndum
    Adrem119           = addyndum.inputFile("Adrem119",        default="DATA/Montagner_model/Adrem119")
    glob_prem3sm01     = addyndum.inputFile("glob_prem3sm01",  default="DATA/Montagner_model/glob-prem3sm01")
    globpreman3sm01    = addyndum.inputFile("globpreman3sm01", default="DATA/Montagner_model/globpreman3sm01")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3D anisotropic (Brian Savage)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


class Model3DAnisotropic(Model):

    # name by which the Fortran code knows this model
    MODEL = "Brian_Savage"
    
    # names by which the user refers to this model
    componentNames = [ "3D_anisotropic", "Brian_Savage" ]

    # initialization
    def _init(self):
        Model._init(self)
        # force 'attenuation' to 'True'
        if not self.ATTENUATION:
            self._warning.log("setting 'attenuation' to 'True' for 3D anisotropic model")
            self.ATTENUATION = True
        return


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# other assorted models
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#                                            name by which the Fortran code
#     Python class name                      knows the model                        names by which the user refers to the model
#     -------------------------------------  -------------------------------------  -------------------------------------------------
class ModelIsotropicPrem(Model):             MODEL = "isotropic_prem";              componentNames = [ "isotropic_prem" ]
class ModelTransverslyIsotropicPrem(Model):  MODEL = "transversly_isotropic_prem";  componentNames = [ "transversly_isotropic_prem" ]
class ModelIaspei(Model):                    MODEL = "iasp91";                      componentNames = [ "iaspei", "iasp91" ]
class ModelAK135(Model):                     MODEL = "ak135";                       componentNames = [ "ak135" ]


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# model facility
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


from cig.addyndum.inventory import builtInComponents
builtins = builtInComponents(globals())


def model(*args, **kwds):
    from cig.addyndum.inventory import Facility
    kwds['builtins'] = builtins
    return Facility(*args, **kwds)


# end of file
