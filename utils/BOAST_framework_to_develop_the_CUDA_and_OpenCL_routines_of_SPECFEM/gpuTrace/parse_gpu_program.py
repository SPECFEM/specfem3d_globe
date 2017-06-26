#! /usr/bin/env python3

import sys
import os
import collections

import gdb_helper

programs = {}
def ocl_prepare_program(uid, lines):
    progr_lines = []
    for line in lines:
        progr_lines += line.split("\n")
    programs[uid] = progr_lines

def ocl_parse(program_uid, kernel_name):
    def get_lines():
        yield from programs[program_uid]

    return parse(get_lines(), kernel_name)

def cuda_parse(source_name, kernel_name):
    def get_lines():
        with open(source_name, "r") as fsource:
            yield from fsource.readlines()

    arglist = parse(get_lines(), kernel_name)

    arguments = list(zip(arglist[0::2], arglist[1::2]))

    return arguments

def parse(lines, kernel_name):
    if "<" in kernel_name:
        kernel_name = kernel_name.split("<")[0]
    else:
        kernel_name = "void {}".format(kernel_name)

    parameters = []
    while True:
        line = lines.send(None)
        if " __attribute__" in line:
            idx = line.index(" __attribute__")
            closeAt = idx
            parents = 0
            inside = False
            while not (inside and parents == 0):
                if line[closeAt] == "(":
                    inside = True
                    parents += 1
                elif line[closeAt] == ")":
                    parents -= 1
                closeAt += 1

            line = line[:idx] + line[closeAt:]

        if kernel_name in line:
            while not "{" in line:
                line += " " + lines.send(None).strip()
            break

    for param in line.split("(")[1].split(")")[0].split(","):
        param = param.strip()
        type_ = " ".join([w for w in param.split()[:-1] if not w.startswith("__")])
        name = param.split()[-1]
        while name.startswith("*"):
            name = name[1:]
            type_ = type_ + " *"
        parameters.append(type_)
        parameters.append(name)

    return parameters

def cuda_get_raw_lookup_table(binary=None):
    lookup_table = collections.OrderedDict()
    for addr, info in cuda_get_lookup_table(binary).items():
        # convert str(addr) to int
        lookup_table[int(addr, 16)] = info

    return lookup_table

def cuda_get_lookup_table(binary=None):
    if binary is None:
        binary = os.readlink("/proc/self/exe")
    lookup_table = collections.OrderedDict()

    symbols = gdb_helper.get_cuda_kernel_names(binary)

    for symb, loc, address in gdb_helper.get_symbol_location(symbols, binary):
        params = cuda_parse(loc, symb)
        lookup_table[address] = symb, params

    return lookup_table


if __name__ == "__main__":
    binary = sys.argv[1] if len(sys.argv) > 1 else "./hello"
    symbols = gdb_helper.get_cuda_kernel_names()

    for symb, loc, address in gdb_helper.get_symbol_location(symbols):
        print(symb)
        print(cuda_parse(loc, symb))
