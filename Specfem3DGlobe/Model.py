#!/usr/bin/env python


from pyre.components.Component import Component
from pyre.inventory.Facility import Facility
import Specfem3DGlobeCode


class Model(Component):

    class Inventory(Component.Inventory):

        from pyre.inventory import bool, inputFile
        
        ATTENUATION                   = bool("attenuation")
        ELLIPTICITY                   = bool("ellipticity")
        GRAVITY                       = bool("gravity")
        OCEANS                        = bool("oceans")
        ROTATION                      = bool("rotation")
        TOPOGRAPHY                    = bool("topography")

        # The hardwired parameters NX_BATHY, NY_BATHY, and
        # RESOLUTION_TOPO_FILE would also have to be Pyrized for the
        # following item to be truly useful...  however this would
        # force the array 'ibathy_topo' to be allocatable.
        topoBathyFile                 = inputFile("topo-bathy-file", default="DATA/topo_bathy/topo_bathy_etopo4_smoothed_window7.dat")
        
    def __init__(self, name):
        Component.__init__(self, name, "model")
        self.aliases.extend(self.classAliases)
        self.PATHNAME_TOPO_FILE = None

    def _init(self):
        Component._init(self)
        # Access our InputFile inventory items to make sure they're
        # readable.  (They will be reopened by the Fortran code.)
        if self.inventory.TOPOGRAPHY or self.inventory.OCEANS:
            f = self.inventory.topoBathyFile
            self.PATHNAME_TOPO_FILE = f.name
            f.close()


# built-in models

def BuiltInModel(name, *aliases):
    return type(Model)(
        'BuiltInModel(%s)' % name,
        (Model,),
        {'className': name,
        'classAliases': [name] + list(aliases)}
        )

class Model3DIsotropic(Model):
    className = "s20rts"
    classAliases = ["s20rts", "3D_isotropic"]
    class Inventory(Model.Inventory):
        from pyre.inventory import inputFile
        CNtype2            = inputFile("CNtype2",           default="DATA/crust2.0/CNtype2.txt")
        CNtype2_key_modif  = inputFile("CNtype2_key_modif", default="DATA/crust2.0/CNtype2_key_modif.txt")
        P12                = inputFile("P12",               default="DATA/s20rts/P12.dat")
        S20RTS             = inputFile("S20RTS",            default="DATA/s20rts/S20RTS.dat")
    def __init__(self, name):
        Model.__init__(self, name)
        self.CNtype2           = None
        self.CNtype2_key_modif = None
        self.P12               = None
        self.S20RTS            = None
    def _init(self):
        Model._init(self)
        # Access our InputFile inventory items to make sure they're
        # readable.  (They will be reopened by the Fortran code.)
        f = self.inventory.CNtype2;            self.CNtype2           = f.name;  f.close()
        f = self.inventory.CNtype2_key_modif;  self.CNtype2_key_modif = f.name;  f.close()
        f = self.inventory.P12;                self.P12               = f.name;  f.close()
        f = self.inventory.S20RTS;             self.S20RTS            = f.name;  f.close()

class Model3DAttenuation(Model):
    className = "Min_Chen"
    classAliases = ["Min_Chen", "3D_attenuation"]
    class Inventory(Model.Inventory):
        from pyre.inventory import inputFile
        Adrem119           = inputFile("Adrem119",        default="DATA/Montagner_model/Adrem119")
        glob_prem3sm01     = inputFile("glob_prem3sm01",  default="DATA/Montagner_model/glob-prem3sm01")
        globpreman3sm01    = inputFile("globpreman3sm01", default="DATA/Montagner_model/globpreman3sm01")
    def _init(self):
        Model._init(self)
        self.inventory.Adrem119.close()
        self.inventory.glob_prem3sm01.close()
        self.inventory.globpreman3sm01.close()

builtInModelClasses = [
    BuiltInModel("isotropic_prem"),
    BuiltInModel("transversly_isotropic_prem"),
    BuiltInModel("iaspei", "iasp91"),
    BuiltInModel("ak135"),
    Model3DIsotropic,
    BuiltInModel("Brian_Savage", "3D_anisotropic"),
    Model3DAttenuation,
    ]

def retrieveBuiltInModelClass(componentName):
    for cls in builtInModelClasses:
        for name in cls.classAliases:
            if name == componentName:
                return cls
    return None


# model facility

class ModelFacility(Facility):

    def __init__(self, name):
        Facility.__init__(self, name, default="isotropic_prem")
        return
    
    def _retrieveComponent(self, instance, componentName):
        cls = retrieveBuiltInModelClass(componentName)
        if cls is None:
            return Facility._retrieveComponent(self, instance, componentName)
        model = cls(componentName)
        model.aliases.append(self.name)
        import pyre.parsing.locators
        locator = pyre.parsing.locators.simple('built-in')
        return model, locator


# end of file
