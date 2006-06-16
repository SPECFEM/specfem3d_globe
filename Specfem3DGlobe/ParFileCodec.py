#!/usr/bin/env python


from pyre.inventory.odb.Registry import Registry
from pyre.odb.fs.CodecODB import CodecODB
import pyre.parsing.locators as locators


class ParFileCodec(CodecODB):

    
    def __init__(self, name):
        CodecODB.__init__(self, encoding='specfem')
        self.name = name
        return

    
    def resolve(self, db):
        return db

    
    def _decode(self, shelf):
        root = Registry("root")
        self._parse(shelf.name, root)
        shelf['inventory'] = root
        shelf._frozen = True
        return


    def _parse(self, pathname, root):
        # Technically, the parameters must be in a specific order, but
        # we don't enforce that here.
        from cig.addyndum.util import setPropertyWithPath
        f = open(pathname, "r")
        lineno = 0
        for line in f:
            lineno = lineno + 1
            if line[0] == '#':
                continue
            tokens = line.split()
            if not len(tokens):
                continue
            if tokens[1] != '=':
                raise RuntimeError("%s: line %d: '%s' unexpected" %
                                   (filename, lineno, tokens[1]))
            f90Var = tokens[0]
            path = self.name + '.' + xlate[f90Var]
            value = tokens[2]
            if value == ".true.":
                value = "True"
            elif value == ".false.":
                value = "False"
            elif (f90Var[-11:] == "_IN_DEGREES" or
                  f90Var == 'GAMMA_ROTATION_AZIMUTH'):
                value += "*deg"
            elif f90Var[-11:] == "_IN_MINUTES":
                value += "*minute"
            elif f90Var == 'SIMULATION_TYPE':
                choice = ['0', 'forward', 'adjoint', 'both']
                v = int(value)
                if 0 < v and v < len(choice):
                    value = choice[v]
            value = value.replace(".d", ".", 1)
            locator = locators.file(pathname, lineno)
            setPropertyWithPath(root, path, value, locator)
        f.close()
        return

    
xlate = {
                      'SIMULATION_TYPE': 'solver.simulation-type',
                         'SAVE_FORWARD': 'solver.save-forward',
                              'NCHUNKS': 'mesher.nchunks',
          'ANGULAR_WIDTH_XI_IN_DEGREES': 'mesher.angular-width-xi',
         'ANGULAR_WIDTH_ETA_IN_DEGREES': 'mesher.angular-width-eta',
           'CENTER_LATITUDE_IN_DEGREES': 'mesher.center-latitude',
          'CENTER_LONGITUDE_IN_DEGREES': 'mesher.center-longitude',
               'GAMMA_ROTATION_AZIMUTH': 'mesher.gamma-rotation-azimuth',
                               'NEX_XI': 'mesher.nex-xi',
                              'NEX_ETA': 'mesher.nex-eta',
                             'NPROC_XI': 'mesher.nproc-xi',
                            'NPROC_ETA': 'mesher.nproc-eta',
                                'MODEL': 'model',
                               'OCEANS': 'model.oceans',
                          'ELLIPTICITY': 'model.ellipticity',
                           'TOPOGRAPHY': 'model.topography',
                              'GRAVITY': 'model.gravity',
                             'ROTATION': 'model.rotation',
                          'ATTENUATION': 'model.attenuation',
                 'ABSORBING_CONDITIONS': 'solver.absorbing-conditions',
             'RECORD_LENGTH_IN_MINUTES': 'solver.record-length',
                        'MOVIE_SURFACE': 'solver.movie-surface',
                         'MOVIE_VOLUME': 'solver.movie-volume',
                'NTSTEP_BETWEEN_FRAMES': 'solver.ntstep-between-frames',
                           'HDUR_MOVIE': 'solver.hdur-movie',
                      'SAVE_MESH_FILES': 'mesher.save-files',
                       'NUMBER_OF_RUNS': 'solver.number-of-runs',
                   'NUMBER_OF_THIS_RUN': 'solver.number-of-this-run',
                           'LOCAL_PATH': 'scratch-dir',
           'NTSTEP_BETWEEN_OUTPUT_INFO': 'solver.ntstep-between-output-info',
        'NTSTEP_BETWEEN_OUTPUT_SEISMOS': 'solver.ntstep-between-output-seismos',
              'RECEIVERS_CAN_BE_BURIED': 'solver.receivers-can-be-buried',
           'PRINT_SOURCE_TIME_FUNCTION': 'solver.print-source-time-function',
    }


# end of file
