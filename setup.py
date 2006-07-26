
import os, os.path
from cig_setup import setup

setup(
    name = 'Specfem3DGlobe', 
    version = '3.5',
    url = 'http://www.gps.caltech.edu/~jtromp/research/downloads.html',
    author = 'Dimitri Komatitsch and Jeroen Tromp',
    author_email = 'jtromp AT caltech.edu',
    packages = [ 'Specfem3DGlobe' ],
    
    install_requires = [
    'cig >= 1.0dev-4103, < 2.0, == dev',
    'pythia >= 0.8-1.0dev-r4100, < 0.9, == dev',
    ],
    
    dependency_links = [
    'svn://geodynamics.org/cig/cs/framework/trunk#egg=cig-dev',
    'svn://geodynamics.org/cig/cs/pythia/trunk#egg=pythia-dev',
    ],

    interpreter = os.path.join(os.getcwd(), "pyspecfem3D"),
    entry_points = {
    'console_scripts': [
    'xspecfem3D = Specfem3DGlobe.Specfem:main',
    ],
    },
)
