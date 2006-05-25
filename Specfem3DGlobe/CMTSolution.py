#!/usr/bin/env python


def cmtValidator(f):
    NLINES_PER_CMTSOLUTION_SOURCE = 13 # constants.h
    lineTally = 0
    for line in f:
        lineTally = lineTally + 1
    if lineTally % NLINES_PER_CMTSOLUTION_SOURCE != 0:
        raise ValueError("total number of lines in 'cmt-solution' file '%s' should be a multiple of %d"
                         % (f.name, NLINES_PER_CMTSOLUTION_SOURCE))
    NSOURCES = lineTally / NLINES_PER_CMTSOLUTION_SOURCE
    if NSOURCES < 1:
        raise ValueError("need at least one source in 'cmt-solution' file '%s'" % f.name)
    f.seek(0)
    return f


# end of file
