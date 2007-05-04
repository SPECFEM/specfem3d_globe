
from archimedes import use_merlin
use_merlin()

from merlin import setup

setup(
    name = 'Specfem3DGlobe', 
    version = '3.6',
    url = 'http://www.gps.caltech.edu/~jtromp/research/downloads.html',
    author = 'Dimitri Komatitsch and Jeroen Tromp',
    author_email = 'jtromp AT caltech.edu',
    packages = [ 'Specfem3DGlobe' ],
    
    install_requires = [
    'pythia[mpi] >= 0.8.1.3, < 0.8.2a',
    ],
    
    entry_points = {
    'console_scripts': [
    'xspecfem3D = Specfem3DGlobe.Specfem:main',
    'xcreate_movie_AVS_DX = Specfem3DGlobe.Specfem:create_movie_AVS_DX',
    'xsfdaemon = Specfem3DGlobe.Daemon:main',
    ],
    },
)
