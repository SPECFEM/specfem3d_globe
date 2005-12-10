#!/usr/bin/env python


from pyre.components.Component import Component
from pyre.inventory.Facility import Facility
import Specfem3DGlobeCode


class Model(Component):

    class Inventory(Component.Inventory):

        from pyre.inventory import bool
        
        ATTENUATION                   = bool("attenuation")
        ELLIPTICITY                   = bool("ellipticity")
        GRAVITY                       = bool("gravity")
        OCEANS                        = bool("oceans")
        ROTATION                      = bool("rotation")
        TOPOGRAPHY                    = bool("topography")
        
    def __init__(self, name):
        Component.__init__(self, name, "model")
        self.aliases.extend(self.classAliases)


# built-in models

def BuiltInModel(name, *aliases):
    return type(Model)(
        'BuiltInModel(%s)' % name,
        (Model,),
        {'className': name,
        'classAliases': [name] + list(aliases)}
        )

builtInModelClasses = [
    BuiltInModel("isotropic_prem"),
    BuiltInModel("transversly_isotropic_prem"),
    BuiltInModel("iaspei", "iasp91"),
    BuiltInModel("ak135"),
    BuiltInModel("s20rts", "3D_isotropic"),
    BuiltInModel("Brian_Savage", "3D_anisotropic"),
    BuiltInModel("Min_Chen", "3D_attenuation"),
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
