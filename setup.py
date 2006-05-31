
import os, sys

sys.path.insert(1, "Specfem3DGlobe/tools")

from Cheetah.Template import Template

wd = os.getcwd()
template = Template(file="Specfem3DGlobe/xspecfem3D.tmpl")
template.interpreter = os.path.join(wd, "pyspecfem3D")
template.srcdir = wd
try:
    script = open("xspecfem3D", "w")
    print >> script, template
except Exception:
    os.remove("xspecfem3D")
    raise
