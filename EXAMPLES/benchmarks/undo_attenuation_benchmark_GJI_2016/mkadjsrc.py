#!/usr/bin/env python3

import argparse
import os

import numpy as np
import matplotlib.pyplot as plt
import obspy
from obspy.core import Stream, Trace, Stats, UTCDateTime
assert hasattr(obspy, '__version__') and obspy.__version__[0] >= '1', \
    'ObsPy is too old to filter correctly. Please install 1.0.0 or newer.'


START_50_100 = 3395
END_50_100 = 3595
START_100_200 = 3390
END_100_200 = 3825
COMPONENT = 'Z'


def read_specfem_seismogram(output_files, network, station, band):
    st = Stream()
    for component in 'ZNE':
        channel = '%sX%s' % (band, component)
        filename = os.path.join(output_files,
                                '%s.%s.%s.sem.ascii' % (network, station,
                                                        channel))
        tmp = np.genfromtxt(filename)

        stats = Stats()
        stats.network = network
        stats.station = station
        stats.channel = channel
        stats.delta = tmp[1, 0] - tmp[0, 0]
        stats.npts = tmp.shape[0]
        stats.starttime = tmp[0, 0]

        tr = Trace(tmp[:, 1], stats)
        st += tr

    return st


def plot_seismo(data, title, markers=(), colours=(), units='m', **kwargs):
    if not args.plot:
        return

    if isinstance(units, str):
        units = [units] * len(data)

    fig = data.plot(type='relative', reftime=UTCDateTime(0), handle=True,
                    **kwargs)

    for ax, unit in zip(fig.get_axes(), units):
        for mark, col in zip(markers, colours):
            ax.axvline(mark, c=col, ls='--', lw=1)
        ax.set_ylabel(unit)

    fig.suptitle(title)

    if output:
        output.savefig(fig)


def bandpass(data, low_period, high_period, markers=(), colours=()):
    filt = data.copy().filter('bandpass',
                              freqmin=1.0 / high_period,
                              freqmax=1.0 / low_period,
                              zerophase=True)
    vel = filt.copy().differentiate()

    plot_seismo(filt, 'Filtered %ds - %ds Displacement' % (low_period,
                                                           high_period),
                markers=markers, colours=colours)
    plot_seismo(vel, 'Filtered %ds - %ds Velocity' % (low_period, high_period),
                units='m/s', markers=markers, colours=colours)

    return filt, vel


def create_adjsrc(all_data, t1, t2, component, EPS=1e-17):
    """
    Create an adjoint source time function.

    Based on the calculation in utils/adjoint_source/traveltime, but without
    any extraneous work. If you want to calculate a traveltime source, then you
    need to input the velocity (derivative) instead of the displacement as
    input.

    Parameters
    ----------

    all_data : obspy.core.stream.Stream
        A Stream of each time series.
    t1, t2 : float
        The start and end times for the window.
    component : str
        Which component to use; an empty string signifies all components.
    EPS : float
        Tolerance below which the norm will be treated as zero.

    Returns
    -------
    dict
        A dictionary of the same input time series components, but with the
        values windowed and normalized within the requested range.

    """
    adj_src = all_data.copy()
    for tr in adj_src:
        if component and tr.stats.channel[-1] != component:
            # A specific component is chosen and this is not it.
            tr.data.fill(0)
            continue

        time = tr.times() + (tr.stats.starttime - UTCDateTime(0))
        data = tr.data

        # Select time range.
        mask = np.where((t1 <= time) & (time <= t2))
        ts = mask[0][0]
        te = mask[0][-1]

        # Create time window (parabola shaped).
        tw = np.zeros_like(time)
        tw[ts:te] = 1 - (2 * (np.arange(ts, te) - ts) / (te - ts) - 1) ** 2

        # Find normalization factor.
        norm = np.sum(tw * data * data)
        print('ID =', tr.id, 'norm =', norm)

        # Create adjoint source.
        if abs(norm) > EPS:
            adj = - data * tw / norm
        else:
            print('norm < EPS for component', tr.stats.channel[-1])
            adj = np.zeros_like(data)

        tr.data = adj

    return adj_src


def generate_adjoint_files(basedir, section, data):
    section_dir = os.path.join(basedir, section)
    try:
        os.mkdir(section_dir)
    except OSError:
        pass

    for name in ['DATA', 'DATABASES_MPI', 'bin']:
        try:
            os.symlink(os.path.join(os.pardir, name),
                       os.path.join(section_dir, name))
        except OSError:
            pass

    try:
        os.mkdir(os.path.join(section_dir, 'OUTPUT_FILES'))
    except OSError:
        pass
    source = os.path.join(os.pardir, os.pardir, 'OUTPUT_FILES',
                          'addressing.txt')
    try:
        os.symlink(source,
                   os.path.join(section_dir, 'OUTPUT_FILES', 'addressing.txt'))
    except OSError:
        pass

    sem_dir = os.path.join(section_dir, 'SEM')
    try:
        os.mkdir(sem_dir)
    except OSError:
        pass

    for tr in data:
        name = os.path.join(sem_dir, '%s.adj' % (tr.id.replace('..', '.'), ))
        output = np.column_stack((
            tr.times() + (tr.stats.starttime - UTCDateTime(0)),
            tr.data))
        np.savetxt(name, output)


# Read any command-line arguments.
parser = argparse.ArgumentParser(description='Generate adjoint source files.')
parser.add_argument('synthetics',
                    help='Top-level directory for synthetics simulation.')
parser.add_argument('adjoint',
                    help='Top-level directory to place adjoint source.')
parser.add_argument('-p', '--plot', action='store_true',
                    help='Plot the seismograms and adjoint sources.')
parser.add_argument('-o', '--output', help='PDF file to place plots.')
args = parser.parse_args()
plot = args.plot
if args.plot and args.output:
    from matplotlib.backends.backend_pdf import PdfPages
    output = PdfPages(args.output)
else:
    output = None

synthetics_dir = os.path.join(args.synthetics, 'OUTPUT_FILES')

# Load output files.
st = read_specfem_seismogram(synthetics_dir, 'SY', 'STA00', 'M')
print('Read in', st)

plot_seismo(st, 'Raw Displacement Seismograms',
            markers=(START_50_100, END_50_100, START_100_200, END_100_200),
            colours='rrgg')

# Detrend (not much change here).
st.detrend('linear')

plot_seismo(st, 'Detrended Displacement Seismograms',
            markers=(START_50_100, END_50_100, START_100_200, END_100_200),
            colours='rrgg')

# Two filtered results:
#   * First band from 50s to 100s
#   * Second band from 100s to 200s

# 50s - 100s
filt_50_100, vel_50_100 = bandpass(st, 50, 100,
                                   markers=(START_50_100, END_50_100),
                                   colours='rr')

# 100s - 200s
filt_100_200, vel_100_200 = bandpass(st, 100, 200,
                                     markers=(START_100_200, END_100_200),
                                     colours='gg')

# Do cutting

cut_50_100 = create_adjsrc(vel_50_100, START_50_100, END_50_100, COMPONENT)
cut_100_200 = create_adjsrc(vel_100_200, START_100_200, END_100_200, COMPONENT)

plot_seismo(cut_50_100, 'Adjoint Source 50s - 100s', units='m/s',
            markers=(START_50_100, END_50_100), colours='rr')
plot_seismo(cut_100_200, 'Adjoint Source 100s - 200s', units='m/s',
            markers=(START_100_200, END_100_200), colours='gg')

# Plot together

st_50_100 = Stream()
for this_st, name in zip([st, filt_50_100, vel_50_100],
                         ['1-raw', '2-filt', '3-vel']):
    for tr in this_st.select(component=COMPONENT):
        tr_new = tr.copy()
        tr_new.stats.location = name
        st_50_100 += tr_new

st_100_200 = Stream()
for this_st, name in zip([st, filt_100_200, vel_100_200],
                         ['1-raw', '2-filt', '3-vel']):
    for tr in this_st.select(component=COMPONENT):
        tr_new = tr.copy()
        tr_new.stats.location = name
        st_100_200 += tr_new

plot_seismo(st_50_100, 'All Stages 50s - 100s', units=('m', 'm', 'm/s'),
            markers=(START_50_100, END_50_100), colours='rr',
            equal_scale=False)
plot_seismo(st_100_200, 'All Stages 100s - 200s', units=('m', 'm', 'm/s'),
            markers=(START_100_200, END_100_200), colours='gg',
            equal_scale=False)

# Write output

generate_adjoint_files(args.adjoint, '50s-100s', cut_50_100)
generate_adjoint_files(args.adjoint, '100s-200s', cut_100_200)

if args.plot:
    if output:
        output.close()
    else:
        plt.show()
