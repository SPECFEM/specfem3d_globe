#!/usr/bin/env python


from pyre.components.Component import Component
import Specfem3DGlobeCode


class Solver(Component):
    
    class Inventory(Component.Inventory):

        from pyre.inventory import bool, float, int
    
        MOVIE_SURFACE                 = bool("movie-surface")
        MOVIE_VOLUME                  = bool("movie-volume")
        RECEIVERS_CAN_BE_BURIED       = bool("receivers-can-be-buried")
        PRINT_SOURCE_TIME_FUNCTION    = bool("print-source-time-function")
        
        HDUR_MOVIE                    = float("hdur-movie")
        RECORD_LENGTH_IN_MINUTES      = float("record-length")
        
        NTSTEP_BETWEEN_FRAMES         = int("ntstep-between-frames")
        NTSTEP_BETWEEN_OUTPUT_INFO    = int("ntstep-between-output-info")
        NTSTEP_BETWEEN_OUTPUT_SEISMOS = int("ntstep-between-output-seismos")
        NUMBER_OF_RUNS                = int("number-of-runs")
        NUMBER_OF_THIS_RUN            = int("number-of-this-run")
        
    def __init__(self, name):
        Component.__init__(self, name, "solver")


# end of file
