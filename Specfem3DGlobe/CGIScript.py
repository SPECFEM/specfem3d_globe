#!/usr/bin/env python


from opal.applications.CGI import CGI
from Specfem import Specfem
from Cheetah.Template import Template
import os


class CGIScript(CGI, Specfem):

    class Inventory(CGI.Inventory, Specfem.Inventory):
        from pyre.inventory import bool, facility, outputFile, str
        output = outputFile("output")
    
    def dumpConfiguration(self):
        configuration = self.retrieveUsableConfiguration()
        print >> self.inventory.output, "\n".join(self.weaver.render(configuration))

    def main(self, *args, **kwds):
        import Specfem3DGlobe
	scriptDir = Specfem3DGlobe.__path__[0]
        simulation_type = self.inventory.solver.inventory.simulation_type
	# os.environ['REQUEST_METHOD']
	from os.path import join
        template = Template(file=join(scriptDir, "page1.html.tmpl"), searchList=[locals()])
	print template

    def __init__(self):
        super(CGIScript, self).__init__("specfem")


if __name__ == '__main__':
    cgiScript = CGIScript()
    cgiScript.run()


# end of file
