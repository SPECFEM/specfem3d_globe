#!/usr/bin/env python
#
# reads in an HDF5-format seismogram file
# and outputs infos about content and checks validity
#
from __future__ import print_function

import sys,os

try:
    import h5py
except:
    print("Error importing h5py, check if h5py module is installed")
    print("")
    # python version
    print("python version:")
    print(sys.version)
    print("")
    # import module paths
    print("module paths:")
    for path in sys.path:
        print(path)
    print("")
    sys.exit(1)

try:
    import numpy as np
except:
    print("Error importing numpy, check if numpy module is installed")
    print("")
    # python version
    print("python version:")
    print(sys.version)
    print("")
    # import module paths
    print("module paths:")
    for path in sys.path:
        print(path)
    print("")
    sys.exit(1)

# visualization
import matplotlib.pyplot as plt

######################################################################

# read and plot the seismograms in hdf5 format
def plot_HDF5_data_info(file_in,show_plot,convert_waveforms):
    print("")
    print("reading file: ",file_in)
    print("")

    # initializes data
    time = None
    networks = None
    stations = None
    dists = None
    channel = None
    aux = None
    prov = None
    data = None
    data_label = ""
    keys = []

    # open
    with h5py.File(file_in, "r") as f_in:
        # contents
        # coords                   Dataset {3, 4}
        # displ                    Dataset {4, 5000, 3}
        # network                  Dataset {4}
        # station                  Dataset {4}
        # time                     Dataset {5000}
        keys = list(f_in.keys())

        # file info
        print("file info...")
        print("  keys: ",keys)
        print("")

        # debug layout
        #print("keys:")
        #for key in keys:
        #    print(f"  {key}: shape = {f[key].shape}")
        #print("")

        if "AuxiliaryData" in keys:
            aux = f_in["AuxiliaryData"]
            for item in aux: print("  AuxiliaryData: ",item)

        if "Provenance" in keys:
            prov = f_in["Provenance"]
            for item in prov: print("  Provenance: ",item)


        # load data
        print("  basic data:")

        # time
        if "time" in keys:
            time = f_in["time"][:]
        # networks
        if "networks" in keys:
            networks = f_in["networks"][:]
        # station names
        if "stations" in keys:
            stations = f_in["stations"][:]
        # epicentral distances
        if "dists" in keys:
            dists = f_in["dists"][:]

        # info
        if not stations is None:
            nrec = len(stations)
            print("  number of stations  : ",nrec)

        if not time is None:
            ntime = len(time)
            print("  number of time steps: ",ntime)

        # checks lengths
        if not networks is None:
            if len(networks) != nrec:
                print("Error: invalid data: station and network arrays have different different lengths",nrec,len(networks))
                sys.exit(1)

        print("")
        print("  seismogram data:")

        # channel names
        channels = []     # e.g. ["MXE","MXN","MXZ"]
        for key in keys:
            if key.startswith(("LX","MX","BX","HX","CX","FX")):
                channels.append(key)

        print("  channels: ",channels)
        nchannels = len(channels)
        for cha in channels:
            # channel names
            print(f"    {cha}: shape = {f_in[cha].shape}")

        # check channels
        if nchannels == 0:
            print("Error: no channel data found")
            sys.exit(1)

        # displacement label
        data_label = "displacement (m)"

        # displacement data
        data = np.zeros((nrec,ntime,nchannels))

        # loops over all receivers
        for irec in range(nrec):
            # gets all channel data
            displ = np.zeros((ntime,nchannels))
            for i,cha in enumerate(channels):
                trace = f_in[cha][irec,:]
                displ[:,i] = trace[:]

            # debug data
            #print(f"  displacement: shape = {displ.shape}")

            # stores displacement for receiver
            data[irec,:,:] = displ[:,:]

    # check if data available
    if data is None:
        print("Error: no data found")
        sys.exit(1)

    print("")
    print("file record:")
    print(f"  data: shape = {data.shape}")
    # global min/max
    glob_min = data.min()
    glob_max = data.max()
    print(f"        global min/max = {glob_min:6.2e} / {glob_max:6.2e}")
    print("")

    # infos
    for irec in range(nrec):
        print(f"  station {irec}:")

        # station name
        sta = stations[irec]
        # avoid bytes string issues with strings like b'Hello', converts to text string
        if isinstance(sta, (bytes, bytearray)): sta = sta.decode("utf-8")

        # network name
        net = networks[irec]
        # avoid bytes string issues with strings like b'Hello', converts to text string
        if isinstance(net, (bytes, bytearray)): net = net.decode("utf-8")

        for icomp in range(nchannels):
            # channel
            cha = channels[icomp]
            # avoid bytes string issues with strings like b'Hello', converts to text string
            if isinstance(cha, (bytes, bytearray)): cha = cha.decode("utf-8")

            name = "{}.{}.{}".format(net,sta,cha)

            displ = data[irec,:,icomp]
            val_min = displ.min()
            val_max = displ.max()

            print(f"  {name:<18}  min/max = {val_min:6.2e} / {val_max:6.2e}")
        print("")

    print("")

    # plot
    if show_plot:
        # figure
        # label
        label = data_label
        print("showing: ",label)

        # only plot vertical traces
        preferred_channel = 2

        # global channel min/max
        use_global_range = True
        if use_global_range:
            glob_min = data[:,:,preferred_channel].min()
            glob_max = data[:,:,preferred_channel].max()
            print(f"  channel {channels[preferred_channel]}: global min/max = {glob_min:6.2e} / {glob_max:6.2e}")
            print("")
            # add margin for plotting
            fac_min = 1.05 if glob_min < 0 else 0.95
            fac_max = 1.05 if glob_max > 0 else 0.95
            glob_min *= fac_min
            glob_max *= fac_max

        # plot 5 stations by page
        nstations_per_page = 5

        if nrec % nstations_per_page == 0:
            npages = int(nrec/nstations_per_page)
        else:
            npages = int(nrec/nstations_per_page) + 1

        for ipage in range(npages):
            # single page - plot with nstations_per_page station maximum
            if nrec % nstations_per_page == 0:
                nstas = nstations_per_page
            else:
                if nrec - (ipage+1)*nstations_per_page > 0:
                    nstas = nstations_per_page
                else:
                    nstas = (nrec - ipage*nstations_per_page) % nstations_per_page

            fig, ax = plt.subplots(nstas, figsize=(20, 10))
            for i in range(nstas):
                # station index
                ista = i + ipage*nstations_per_page

                # station name
                sta = stations[ista]
                # avoid bytes string issues with strings like b'Hello', converts to text string
                if isinstance(sta, (bytes, bytearray)): sta = sta.decode("utf-8")

                # network name
                net = networks[ista]
                # avoid bytes string issues with strings like b'Hello', converts to text string
                if isinstance(net, (bytes, bytearray)): net = net.decode("utf-8")

                # channel
                cha = channels[preferred_channel]
                # avoid bytes string issues with strings like b'Hello', converts to text string
                if isinstance(cha, (bytes, bytearray)): cha = cha.decode("utf-8")

                name = "{}.{}.{}".format(net,sta,cha)
                print("  station: ",name)

                # plot z component
                ax[i].plot(time, data[ista, :, preferred_channel], "b")

                # labels
                ax[i].set_title(name)
                ax[i].set_xlabel("time (s)")
                ax[i].set_ylabel(label)
                if use_global_range: ax[i].set_ylim(glob_min,glob_max)

            plt.tight_layout()

            # save as figure
            plt.savefig(f"tmp_figure_{ipage}.jpg")
            print("")
            print(f"  plotted to: tmp_figure_{ipage}.jpg")
            print("")
            
            plt.show()

            # Wait for the user input to terminate the program
            #input("Press any key to terminate the program")

    if convert_waveforms:
        print("converting to ASCII:")

        # takes same output directory as input file
        # for example: ./plot_HDF5_data_info.py OUTPUT_FILES/seismograms.h5
        #              -> output in OUTPUT_FILES/ directory
        out_dir = os.path.dirname(file_in)
        if out_dir == '':
            out_dir = "./"  # -> ./ + DB.STA.XXX.semd
        else:
            out_dir = out_dir + "/"  # OUTPUT_FILES -> OUTPUT_FILES/ + DB.STA.XXX.semd

        # label
        label = data_label

        # loops over all stations
        for irec in range(nrec):
            # station name
            sta = stations[irec]
            # avoid bytes string issues with strings like b'Hello', converts to text string
            if isinstance(sta, (bytes, bytearray)): sta = sta.decode("utf-8")

            # network name
            net = networks[irec]
            # avoid bytes string issues with strings like b'Hello', converts to text string
            if isinstance(net, (bytes, bytearray)): net = net.decode("utf-8")

            # components
            ncomp = 3
            for icomp in range(ncomp):
                # trace data
                # displ w/ shape (nrec,nstep,ndim)
                trace = data[irec, :, icomp]

                length = len(trace)

                # check if time and trace have same length
                if length != len(time):
                    print("Error: different trace / time lengths: {} / {}".format(length,len(time)))
                    sys.exit(1)

                # channel name
                cha = channels[icomp]
                # avoid bytes string issues with strings like b'Hello', converts to text string
                if isinstance(cha, (bytes, bytearray)): cha = cha.decode("utf-8")

                # file ending
                ending = ".sem.ascii.converted"

                name = "{}.{}.{}".format(net,sta,cha)
                print("  station: ",name)

                filename = out_dir + name + ending

                # single file per station and component
                # file output
                f_out = open(filename, "w")

                t0 = time[0]
                # file header
                f_out.write("# STATION %s\n" % (sta))
                f_out.write("# CHANNEL %s\n" % (cha))
                f_out.write("# START_TIME %s\n" % (str(t0)))
                #f_out.write("# SAMP_FREQ %f\n" % (tr.stats.sampling_rate))
                f_out.write("# NDAT %d\n" % (length))

                # data section
                # fills numpy array
                xy = np.empty(shape=(0,2))
                for ii in range(length):
                    t = time[ii]
                    val = trace[ii]
                    xy = np.append(xy,[[t,val]],axis=0)

                # data column
                #print(xy[:,1])
                #print(xy[0,1],xy[1,1],xy[2,1],xy[3,1])

                # saves as ascii
                np.savetxt(f_out, xy, fmt="%f\t%e")
                f_out.close()

                # user output
                print("  written file: ",filename)

    print("")
    print("all done")
    print("")


def usage():
    print("usage: ./plot_HDF5_data_info.py filename[e.g. seismograms.h5] [show] [convert]")
    print(" with")
    print("   filename  - e.g. OUTPUT_FILES/seismograms.h5")
    print("   show      - (optional) plot waveforms")
    print("   convert   - (optional) convert seismograms to ASCII traces")
    sys.exit(1)


if __name__ == "__main__":
    # gets arguments
    do_plot_waveforms = False
    do_convert_waveforms = False

    if len(sys.argv) == 2:
        filename = sys.argv[1]
    elif len(sys.argv) == 3:
        filename = sys.argv[1]
        if sys.argv[2] == "show": do_plot_waveforms = True
    elif len(sys.argv) == 4:
        filename = sys.argv[1]
        if sys.argv[2] == "show": do_plot_waveforms = True
        if sys.argv[3] == "convert": do_convert_waveforms = True
    else:
        usage()
        sys.exit(1)

    plot_HDF5_data_info(filename,do_plot_waveforms,do_convert_waveforms)


