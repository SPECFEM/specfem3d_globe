#! /usr/bin/env python3

import subprocess
import re
import sys
import os

SPLITTER = "CUT-CUT-CUT"
GDB_COMMAND = ["gdb", "-ex", "set filename-display absolute", "-ex", 'printf "{}\n"'.format(SPLITTER),  '-quiet']

LIST_KERNELS = ["-ex", 'info functions __device_stub__']
QUIT = ['-ex', 'quit']

def get_cuda_kernel_prototypes(binary):
    os.putenv("LD_PRELOAD", "")
    output = subprocess.Popen(GDB_COMMAND + [binary] + LIST_KERNELS + QUIT,
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()[0]

    functions = " ".join([f for f in str(output).split(SPLITTER)[1].split("\n") \
                              if not f.startswith("All functions") and \
                              not f.startswith("File ") \
                              and f])

    #full_fct_name = functions.split("\\nvoid ")[7].split("(")[0]
    kern_list = [e.split("Pi")[0]\
                     .split("PK")[0]\
                     .split("_implii")[0]\
                     .split("ii")[0] for e in
                 re.findall(r'Z\d+(\S+(?!(?:Pi)))(?:Pf|Pi|Pc)', functions)]

    kern_list = sorted([e if not e.endswith("i") else e[:-1] for e in kern_list])
    kern_list = [e.replace("ILi", "<").replace("EEv", ">") for e in kern_list]

    print_kernels = []
    for kern in kern_list:
        print_kernels += ["-ex", "print {}".format(kern)]

    output = subprocess.Popen(GDB_COMMAND + [binary] + print_kernels + QUIT,
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()[0]
    locations = [e.split(" = ")[1].strip() for e in str(output).split(SPLITTER)[1].replace("\\n", " ").split("$") if " = " in e]
    return locations

def get_cuda_kernel_names(binary):
    for proto in get_cuda_kernel_prototypes(binary):
        name = proto.split(" <")[1].split("(")[0]
        yield name

def get_symbol_location(symbols, binary):
    for symb in symbols:
        print_locations = ["-ex", "python symb = gdb.lookup_global_symbol('{}'); print '8<{{}} {{}}'.format(symb.symtab, str(symb.value().address).split(" ")[0])".format(symb)]

        output, err = subprocess.Popen(GDB_COMMAND + [binary] + print_locations + QUIT,
                                       stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE).communicate()

        location, address = [e.strip("\" '").replace("\\n", "") for e in str(output).split("8<")[1:]][0].split(" ")
        yield symb, location, address

if __name__ == "__main__":
    symbols = get_cuda_kernel_names()
    for symb, loc, address in get_symbol_location(symbols):
        print(symb, loc, address)
