#!/usr/bin/env python


# courtesy Bob Ippolito  bob@redivi.com 

class MetaCategory(type):
     def __new__(cls, name, bases, dct):
         if '__category_of__' in dct:
             return type.__new__(cls, name, bases, dct)
         if not len(bases) == 1 and isinstance(bases[0], cls):
             raise TypeError("Categories may only have a Category(...) as their base")
         cls = bases[0].__category_of__
         for k,v in dct.iteritems():
             if k == '__module__' or k == '__doc__':
                 continue
             setattr(cls, k, v)
         return cls

def Category(cls):
     return MetaCategory(
         'Category(%s)' % (cls.__name__),
         (object,),
         {'__category_of__': cls}
     )


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Property


from pyre.inventory.Property import Property
from pyre.inventory.Trait import Trait
from journal.devices.Renderer import Renderer as PyreRenderer
import sys


class PropertyValueError(ValueError):
     

     class Renderer(PyreRenderer):
          def __init__(self, header=None, format=None, footer=None):
               PyreRenderer.__init__(
                    self,
                    header="[%(filename)s:%(line)s] property '%(name)s': %(error)s\n >> %(src)s",
                    format=format, footer=footer
                    )


     def __init__(self, *args):
          ValueError.__init__(self, *args)
          self.origExc = sys.exc_info()


     def __str__(self):
          return self.origExc[0].__name__ + ": " + self.origExc[1].__str__()


class PropertyPatches(Category(Property)):

    
    """Add location info to ValueError exceptions.
    
    Instead of _cast, call _richCast;
    instead of validator, call _richValidator.
    
    """
    

    def _set(self, instance, value, locator):
        # None is a special value; it means that a property is not set
        if value is not None:
            # convert
            value = self._richCast(value, locator)
            # validate 
            value = self._richValidator(value, locator)

        # record
        return Trait._set(self, instance, value, locator)


    def _getDefaultValue(self, instance):
        """retrieve the default value and return it along with a locator"""

        value = self.default

        # None is a special value and shouldn't go through the _cast
        if value is not None:
            # convert
            value = self._richCast(value)
            # validate
            value = self._richValidator(value)
        
        import pyre.parsing.locators
        locator = pyre.parsing.locators.simple('default')

        return value, locator

    
    def _richCast(self, value, locator=None):
        try:
            value = self._cast(value)
        except ValueError:
             raise PropertyValueError, (self, locator)
        except KeyError: # e.g., Bool
             raise PropertyValueError, (self, locator)
        except NameError: # e.g., Dimensional
             raise PropertyValueError, (self, locator)
        except TypeError: # Dimensional again...
             raise PropertyValueError, (self, locator)
        except: # ...I suppose anything could happen in the 'eval'
             raise PropertyValueError, (self, locator)
        else:
             return value
    

    def _richValidator(self, value, locator=None):
        if self.validator:
            try:
                value = self.validator(value)
            except ValueError:
                raise PropertyValueError, (self, locator)
        return value


# end of file
