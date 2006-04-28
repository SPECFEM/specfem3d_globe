#!/usr/bin/env python


from opal.applications.CGI import CGI
from Specfem import Specfem
from Cheetah.Template import Template
import os


# hocus-pocus so we don't have to write "inventory" so much
class SearchAdapter(object):
    def __init__(self, component):
        self.component = component
    def __getattr__(self, name):
        from pyre.components.Component import Component
        attr = getattr(self.component.inventory, name)
        if isinstance(attr, Component):
            return SearchAdapter(attr)
        return attr


class CGIScript(CGI, Specfem):

    class Inventory(CGI.Inventory, Specfem.Inventory):
        from pyre.inventory import bool, facility, outputFile, str
        output = outputFile("output")
        template = str("template", default="page1.html.tmpl")
        
    def dumpConfiguration(self):
        configuration = self.retrieveUsableConfiguration()
        print >> self.inventory.output, "\n".join(self.weaver.render(configuration))

    def main(self, *args, **kwds):
	# os.environ['REQUEST_METHOD']
        import Specfem3DGlobe
	from os.path import join
        
	scriptDir = Specfem3DGlobe.__path__[0]
        templateName = join(scriptDir, self.inventory.template)
        template = Template(file=templateName, searchList=[SearchAdapter(self)])
        
	print template

    def __init__(self):
        super(CGIScript, self).__init__("Specfem3DGlobe")


if __name__ == '__main__':
    cgiScript = CGIScript()
    cgiScript.run()


# end of file
