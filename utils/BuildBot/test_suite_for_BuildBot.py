#! /usr/bin/env python3

### Author: Kevin Pouget, Grenoble, France, July 2014

import queue
import sh

SOURCE_DIR = "/home/kevin/cl-specfem3d/SPECFEM3D_GLOBE"
BUILD_DIR = "/home/kevin/specfem-testsuite"
PAR_FILE = "DATA/Par_file"

DEBUG = False

to_import = ("cd", "mkdir", "make", "sed", "grep", "gnuplot", "{}/configure".format(SOURCE_DIR))

def process_output(line):
    print(line)

def bake_echo(name):
    def echo(*kwargs):
        print("| {} {}".format(name, " ".join(kwargs)))
        return "(running {} {})".format(name, " ".join(kwargs))
    return echo

for imp in to_import:
    name = imp.rpartition("/")[-1]
    globals()[name] = sh.__getattr__(imp) if not DEBUG else bake_echo(name)

def lst_to_ctes(*lst):
    for cte in lst:
        globals()[cte.upper()] = cte
    return lst

###############################################################################################

CUDA_CONF = ["CUDA_INC=/usr/local/cuda-5.5.22/include/", "CUDA_LIB=/usr/local/cuda-5.5.22/lib64"]

CONFIGURE_OPT = {
    OPENCL: ["--with-opencl"],
    CUDA: ["--with-cuda=cuda5"] + CUDA_CONF,
    BOTH: ["--with-opencl", "--with-cuda=cuda5"] + CUDA_CONF
    }

FLAGS = ("USE_TEXTURES_FIELDS", "USE_TEXTURES_CONSTANTS", "MANUALLY_UNROLLED_LOOPS")

RUNTIMES = lst_to_ctes("opencl", "cuda", "both")
CONFFILE_DEFAULT_OPT = {
    "GPU_MODE": ".true.",
    "GPU_RUNTIME": "0",
    "GPU_PLATFORM": "NVIDIA",
    "GPU_DEVICE": "Tesla"
    }

GPU_RUNTIME_OPT = {
    OPENCL: "2",
    CUDA: "1",
    BOTH: "0"
    }

def analyze_results(scenarii):
    cd(BUILD_DIR)
    return
    to_plot = []
    seismo_glob = "{}/inout/OUTPUT_FILES/*.ascii".format(VERSIONS_TO_TEST[0])
    to_plot.append("set terminal pdf")
    to_plot.append("set output 'out.pdf'")

    for seismo in sh.glob(seismo_glob):
        plot = "plot "
        for i, version in enumerate(VERSIONS_TO_TEST):
            right_seismo = seismo.replace(VERSIONS_TO_TEST[0], version)
            plot += "'{}' using 1:2 {}".format(right_seismo, "" if i == len(VERSIONS_TO_TEST) - 1 else ", ")
        to_plot.append(plot)
    to_plot.append("quit")
    plotter = gnuplot(_in=[p+"\n" for p in to_plot])

def set_config_options(options=CONFFILE_DEFAULT_OPT):
    for opt, val in options.items():
        sed("-i", "/{}/d".format(opt), PAR_FILE)
        if DEBUG:
            print("(add to {}: {} = {})".format(PAR_FILE, opt, val))
            continue
        with open(PAR_FILE, "a+") as par_file:
            print("{} = {}".format(opt, val), file=par_file)

class Scenario:
    built_runtimes = []
    scenarii = []
    cpt = 0
    def __init__(self, runtime, make_flags=[], config_options={}):
        self.uid = Scenario.cpt
        Scenario.cpt += 1

        self.runtime = runtime
        self.make_flags = make_flags
        self.config_options = config_options

        Scenario.scenarii.append(self)

    def run(self):
        print("------<Test #{}>----------".format(Scenario.cpt))

        scenario.initialize_build_dir()
        scenario.setup_environment()
        scenario.run_execution()
        print("------</Test #{}>----------\n".format(Scenario.cpt))

    def initialize_build_dir(self):
        if self.runtime in Scenario.built_runtimes:
            set_config_options()
            return

        cd(BUILD_DIR)
        mkdir("{}".format(self.runtime), "-p")
        cd("{}".format(self.runtime))
        print("Configure {}".format(self.runtime))
        print(configure(CONFIGURE_OPT[self.runtime]))

        set_config_options()

        Scenario.built_runtimes.append(self.runtime)

    def setup_environment(self):
        print("prepare the example")
        make("prepare-example")
        set_config_options(self.config_options)

        print("building the example")
        cflags = "'{}'".format(" ".join(["-D{}".format(flag) for flag in self.make_flags])) if self.make_flags else ""

        print(make("CFLAGS={}".format(cflags), _out=process_output, _err=process_output))
        print(make("build-example", _out=process_output, _err=process_output))

    def run_execution(self):
        print("run the example")
        print(make("run-example"))

    def save_results(self):
        print("save the results (TODO)")

if __name__ == "__main__":
    mkdir(BUILD_DIR, "-p")

    for rt in RUNTIMES:
        if rt == OPENCL: continue
        if rt == BOTH:
            for each in (OPENCL, CUDA):
                Scenario(rt, config_options={"GPU_RUNTIME" : GPU_RUNTIME_OPT[each]})
            continue

        Scenario(rt)
        for flag in FLAGS:
            Scenario(rt, [flag])

    for scenario in Scenario.scenarii:
        scenario.run()


    analyze_results(Scenario.scenarii)

