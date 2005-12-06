#!/usr/bin/env python


from pyre.applications.Script import Script as PyreScript


class Script(PyreScript):

    class Inventory(PyreScript.Inventory):

        from pyre.inventory import bool, choice, float, int, str
        
        ANGULAR_WIDTH_ETA_IN_DEGREES  = float("angular-width-eta")
        ANGULAR_WIDTH_XI_IN_DEGREES   = float("angular-width-xi")
        CENTER_LATITUDE_IN_DEGREES    = float("center-latitude")
        CENTER_LONGITUDE_IN_DEGREES   = float("center-longitude")
        GAMMA_ROTATION_AZIMUTH        = float("gamma-rotation-azimuth")
        HDUR_MOVIE                    = float("hdur-movie")
        RECORD_LENGTH_IN_MINUTES      = float("record-length")
        
        NCHUNKS                       = int("nchunks", validator=choice([1,2,3,6]))
        NEX_ETA                       = int("nex-eta")
        NEX_XI                        = int("nex-xi")
        NPROC_ETA                     = int("nproc-eta")
        NPROC_XI                      = int("nproc-xi")
        NTSTEP_BETWEEN_FRAMES         = int("ntstep-between-frames")
        NTSTEP_BETWEEN_OUTPUT_INFO    = int("ntstep-between-output-info")
        NTSTEP_BETWEEN_OUTPUT_SEISMOS = int("ntstep-between-output-seismos")
        NUMBER_OF_RUNS                = int("number-of-runs")
        NUMBER_OF_THIS_RUN            = int("number-of-this-run")

        ABSORBING_CONDITIONS          = bool("absorbing-conditions")
        ATTENUATION                   = bool("attenuation")
        ELLIPTICITY                   = bool("ellipticity")
        GRAVITY                       = bool("gravity")
        MOVIE_SURFACE                 = bool("movie-surface")
        MOVIE_VOLUME                  = bool("movie-volume")
        OCEANS                        = bool("oceans")
        PRINT_SOURCE_TIME_FUNCTION    = bool("print-source-time-function")
        RECEIVERS_CAN_BE_BURIED       = bool("receivers-can-be-buried")
        ROTATION                      = bool("rotation")
        SAVE_MESH_FILES               = bool("save-mesh-files")
        TOPOGRAPHY                    = bool("topography")

        LOCAL_PATH                    = str("local-path")
        MODEL                         = str("model", validator=choice(['isotropic_prem',
                                                                       'transversly_isotropic_prem',
                                                                       'iaspei',
                                                                       's20rts',
                                                                       'Brian_Savage',
                                                                       'Min_Chen',
                                                                       ]))
    
    def readValue(self, name):
        return reduce(getattr, name.split('.'), self)


# end of file
