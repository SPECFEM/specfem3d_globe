#!/usr/bin/env python3

import argparse
import os

import numpy as np
import matplotlib.pyplot as plt
import scipy.signal
# ObsPy has a simpler filtering interface
import obspy.signal


START_50_100 = 3325
END_50_100 = 3635
START_100_200 = 3350
END_100_200 = 3900
COMPONENT = 'Z'


def plot_seismo(data, title, markers=(), colours=(), units='m'):
    if not args.plot:
        return

    fig, ax = plt.subplots(3, 1, figsize=(16, 12), sharex=True)

    for i, comp in enumerate('ZEN'):
        ax[i].plot(data[comp][:, 0], data[comp][:, 1])
        ax[i].set_xlim(data[comp][0, 0], data[comp][-1, 0])

        for mark, col in zip(markers, colours):
            ax[i].axvline(mark, c=col, ls='--')

        ax[i].set_ylabel('%s (%s)' % (comp, units))

    ax[2].set_xlabel('Time (s)')
    fig.suptitle(title)

    if output:
        output.savefig(fig)


def bandpass(data, low_period, high_period, dt, markers=(), colours=()):
    filt = {}
    vel = {}

    for c in 'ZEN':
        filt[c] = np.empty_like(data[c])
        filt[c][:, 0] = data[c][:, 0]
        # NOTE: ObsPy has a bandpass filter, but it is broken in current
        # releases. It should be fixed in master and the next (0.11.0) release.
        tmp = obspy.signal.lowpass(data[c][:, 1], 1.0 / low_period,
                                   df=1.0 / dt, zerophase=True)
        filt[c][:, 1] = obspy.signal.highpass(tmp, 1.0 / high_period,
                                              df=1.0 / dt, zerophase=True)

        vel[c] = np.empty_like(detrend[c])
        vel[c][:, 0] = detrend[c][:, 0]
        vel[c][:, 1] = np.gradient(filt[c][:, 1])

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

    all_data : dict
        A dictionary of each time series indexing by component.
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
    adj_src = {}
    for comp in 'ZEN':
        adj_src[comp] = all_data[comp].copy()

        if component and comp != component:
            # A specific component is chosen and this is not it.
            adj_src[comp][:, 1] = 0.0
            continue

        time = all_data[comp][:, 0]
        data = all_data[comp][:, 1]

        # Select time range.
        mask = np.where((t1 <= time) & (time <= t2))
        ts = mask[0][0]
        te = mask[0][-1]

        # Create time window (parabola shaped).
        tw = np.zeros_like(time)
        tw[ts:te] = 1 - (2 * (np.arange(ts, te) - ts) / (te - ts) - 1) ** 2

        # Find normalization factor.
        norm = np.sum(tw * data * data)
        print('component = ', comp, 'norm = ', norm)

        # Create adjoint source.
        if abs(norm) > EPS:
            adj = - data * tw / norm
        else:
            print('norm < EPS for component', comp)
            adj = 0.0

        adj_src[comp][:, 1] = adj

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

    for c in 'ZEN':
        name = os.path.join(sem_dir, 'SY.STA00.MX%s.adj' % (c, ))
        np.savetxt(name, data[c])


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
files = {}

for c in 'ZEN':
    name = os.path.join(synthetics_dir, 'SY.STA00.MX%s.sem.ascii' % (c, ))
    files[c] = np.genfromtxt(name)
dt = files['Z'][1, 0] - files['Z'][0, 0]

plot_seismo(files, 'Raw Displacement Seismograms',
            markers=(START_50_100, END_50_100, START_100_200, END_100_200),
            colours='rrgg')

# Detrend (not much change here).
detrend = {}
for c in 'ZEN':
    detrend[c] = np.empty_like(files[c])
    detrend[c][:, 0] = files[c][:, 0]
    detrend[c][:, 1] = scipy.signal.detrend(files[c][:, 1])

plot_seismo(detrend, 'Detrended Displacement Seismograms',
            markers=(START_50_100, END_50_100, START_100_200, END_100_200),
            colours='rrgg')

# Two filtered results:
#   * First band from 50s to 100s
#   * Second band from 100s to 200s

# 50s - 100s
filt_50_100, vel_50_100 = bandpass(detrend, 50, 100, dt,
                                   markers=(START_50_100, END_50_100),
                                   colours='rr')

# 100s - 200s
filt_100_200, vel_100_200 = bandpass(detrend, 100, 200, dt,
                                     markers=(START_100_200, END_100_200),
                                     colours='gg')

# Do cutting

cut_50_100 = create_adjsrc(vel_50_100, START_50_100, END_50_100, COMPONENT)
cut_100_200 = create_adjsrc(vel_100_200, START_100_200, END_100_200, COMPONENT)

plot_seismo(cut_50_100, 'Adjoint Source 50s - 100s', units='m/s',
            markers=(START_50_100, END_50_100), colours='rr')
plot_seismo(cut_100_200, 'Adjoint Source 100s - 200s', units='m/s',
            markers=(START_100_200, END_100_200), colours='gg')

# Write output

generate_adjoint_files(args.adjoint, '50s-100s', cut_50_100)
generate_adjoint_files(args.adjoint, '100s-200s', cut_100_200)

if args.plot:
    if output:
        output.close()
    else:
        plt.show()
