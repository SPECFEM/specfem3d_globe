
import os, sys

sys.path.insert(1, "Specfem3DGlobe/tools")

from Cheetah.Template import Template

template = Template(file="Specfem3DGlobe/xspecfem3D.tmpl")
template.interpreter = os.path.join(os.getcwd(), "pyspecfem3D")
script = open("xspecfem3D", "w")
print >> script, template
