#!/usr/bin/env -S python3 -Werror
# -*- coding: utf-8 -*-

import os
import argparse
import configparser

class dependency_structure:
    def __init__(self, filename, program, module, name, dependencies, includes):
        self.filename     = filename
        self.program      = program
        self.module       = module
        self.name         = name
        self.dependencies = dependencies
        self.includes     = includes

def get_program_name_from_program(line):
    start_index = line.find("program") + 8
    end_index   = line.find(",")
    
    return line[start_index:end_index]

def get_module_name_from_module(line):
    start_index = line.find("module") + 7
    end_index   = line.find(",")
    
    return line[start_index:end_index]

def get_module_name_from_use(line):
    start_index = line.find("use") + 4
    end_index   = line.find(",")
    
    return line[start_index:end_index]

def get_file_name_from_include(line):
    start_index = line.find('include "') + 9
    end_index   = line.rfind('"')
    
    return line[start_index:end_index]

def canonicalize_path(path):
    path_split = path.split("/")
    return os.path.join(*path_split)

def canonicalize_paths(path_string):
    path_string_split = path_string.split(' ')
    paths = []
    for path_string_split_i in path_string_split:
        paths.append(canonicalize_path(path_string_split_i))
    return paths

parser = argparse.ArgumentParser(description="Generates a Makefile fragment describing the dependencies of a Fortran project.")
parser.add_argument("file", help="input file to read")
args = parser.parse_args()

config = configparser.ConfigParser()
config.read(args.file)

directories          = config['depends']['directories'].split(' ')
no_existence_check   = config['depends']['no_existence_check'].split(' ')
module_name_mismatch = canonicalize_paths(config['depends']['module_name_mismatch'])
skip_indexing        = canonicalize_paths(config['depends']['skip_indexing'])
depends_file         = canonicalize_path(config['depends']['depends_file'])
before_file          = canonicalize_path(config['depends']['before_file'])
allsrc_raw           = config['depends']['allsrc'].split(' ')
not_allsrc_raw       = config['depends']['not_allsrc'].split(' ')
extra_testnml        = config['depends']['extra_testnml'].split(' ')

fail = False

filepaths = []
for directory in sorted(directories):
    # Canonicalize the filepath, so that (for example), this works if `.\` is in front of the path, like PowerShell does.
    directory = os.path.relpath(directory)
    
    if not os.path.isdir(directory):
        print("ERROR: {} in directories is not a directory.".format(directory))
        fail = True
    else:
        for filename in sorted(os.listdir(directory)):
            if filename.endswith(".f90"):
                filepaths.append(os.path.join(directory, filename))

all_source = set()
for filename in sorted(os.listdir("mk")):
    all_source.add(os.path.join("mk", filename))

if not ((len(allsrc_raw) == 1) and (allsrc_raw[0] == "")):
    for filepath in sorted(allsrc_raw):
        # Canonicalize the filepath, so that (for example), this works if `.\` is in front of the path, like PowerShell does.
        filepath = os.path.relpath(filepath)
        
        if not os.path.isfile(filepath):
            print("ERROR: {} in allsrc is not a file.".format(filepath))
            fail = True
        
        all_source.add(filepath)

not_all_source = set()
if not ((len(not_allsrc_raw) == 1) and (not_allsrc_raw[0] == "")):
    for filepath in sorted(not_allsrc_raw):
        # Canonicalize the filepath, so that (for example), this works if `.\` is in front of the path, like PowerShell does.
        filepath = os.path.relpath(filepath)
        
        if not os.path.isfile(filepath):
            print("ERROR: {} in not_allsrc is not a file.".format(filepath))
            fail = True
        
        print("not", filepath)
        not_all_source.add(filepath)

if fail:
    print("Error(s) encountered, stopping.")
    exit(1)

depstructs = []
all_dependencies = set()
for filepath in sorted(filepaths):
    if filepath in skip_indexing:
        continue
    
    if (not filepath in not_all_source):
        all_source.add(filepath)
    
    with open(filepath, "r") as file_handler:
        program      = False
        module       = False
        name         = None
        dependencies = set() # `use` lines
        includes     = set()
        
        for line in file_handler.readlines():
            if line.strip().startswith("program "):
                program = True
                
                name = get_program_name_from_program(line)
            
            if line.strip().startswith("module "):
                if not line.strip().startswith("module procedure "):
                    module = True
                    
                    name = get_module_name_from_module(line)
            
            if line.strip().startswith("use "):
                if get_module_name_from_use(line) == name:
                    print("ERROR: {} can't depend on itself.".format(filepath))
                    fail = True
                
                dependencies.add(get_module_name_from_use(line))
                all_dependencies.add(get_module_name_from_use(line))
                
                if get_module_name_from_use(line).startswith("iso_"):
                    
                    fail = True
            
            if line.strip().startswith("include "):
                # TODO
                # if get_file_name_from_include(line) == filename:
                    # print("ERROR: {} can't include itself.".format(filepath))
                    # fail = True
                
                includes.add(get_file_name_from_include(line))
        
        print("{} program={} module={} name={} dependencies={} includes={}".format(filepath, program, module, name, dependencies, includes))
        
        if (program and module):
            print("ERROR: {} contains both a program and a module. depends.py assumes that a file contains one or the other, not both.".format(filepath))
            fail = True
        
        if not name is None:
            if not filepath in module_name_mismatch:
                if os.path.basename(filepath) != name + ".f90":
                    print("ERROR: {} contains module or program {}, when I require that the two have the same name.".format(filepath, name))
                    fail = True
            
            if "." in name:
                print("ERROR: The module name for {} can't contain a period as this code assumes there are no periods in the module name.".format(filepath))
                fail = True
        
        if module:
            if os.path.split(filepath)[0] != "src":
                print("ERROR: {} contains a module which is not in the src directory. All modules must be in the src directory.".format(filepath))
                fail = True
        
        if (program or module):
            depstructs.append(dependency_structure(filepath, program, module, name, dependencies, includes))
        else:
            assert name is None, "name should not be set for something which is not a program or module"

for dependency in sorted(all_dependencies):
    if dependency in no_existence_check:
        continue
    
    if not os.path.exists(os.path.join("src", dependency+".f90")):
        print("ERROR: Module dependency {} does not exist.".format(dependency))
        fail = True

if fail:
    print("Error(s) encountered, stopping.")
    exit(1)

for depstruct in depstructs:
    # For programs, I need the entire recursive dependency structure (both modules and `include`s).
    if depstruct.program:
        num_dependencies_prev = 0
        num_dependencies      = len(depstruct.dependencies) + len(depstruct.includes)
        
        while (num_dependencies > num_dependencies_prev):
            for sub_depstruct in depstructs:
                if sub_depstruct.name in depstruct.dependencies:
                    depstruct.dependencies = depstruct.dependencies.union(sub_depstruct.dependencies)
                    depstruct.includes     = depstruct.includes.union(sub_depstruct.includes)
            
            num_dependencies_prev = num_dependencies
            num_dependencies      = len(depstruct.dependencies) + len(depstruct.includes)
        
        print("Complete {} program dependencies: {}".format(depstruct.name, depstruct.dependencies))

with open(depends_file, "w") as output_handler:
    output_handler.write("# Automatically generated by depends.py.\n\n")
    
    # Write module dependencies.
    
    output_handler.write("#######################\n")
    output_handler.write("# Module dependencies #\n")
    output_handler.write("#######################\n\n")
    for depstruct in depstructs:
        if depstruct.module:
            directory = os.path.split(depstruct.filename)[0]
            filename  = os.path.split(depstruct.filename)[1]
            basename  = filename.split(".")[0]
            
            print("Writing module dependencies: {}".format(depstruct.filename))
            
            output_handler.write("{}$(DIR_SEP){}.$(OBJEXT):".format(directory, basename))
            
            for dependency in sorted(depstruct.dependencies):
                if dependency != "build":
                    output_handler.write(" src$(DIR_SEP){}.$(OBJEXT)".format(dependency))
                else:
                    output_handler.write(" src$(DIR_SEP)$(BUILD).$(OBJEXT)")
            
            for include in sorted(depstruct.includes):
                output_handler.write(" src$(DIR_SEP){}".format(include))
            
            output_handler.write(" {}$(DIR_SEP){}\n\n".format(directory, filename))
    
    # Write program dependencies.
    
    test_nmls = []
    output_handler.write("########################\n")
    output_handler.write("# Program dependencies #\n")
    output_handler.write("########################\n\n")
    for depstruct in depstructs:
        if depstruct.program:
            directory = os.path.split(depstruct.filename)[0]
            
            print("Writing program dependencies: {}".format(depstruct.name))
            
            output_handler.write("{}$(BINEXT):".format(depstruct.name))
            
            dependency_string = ""
            for dependency in sorted(depstruct.dependencies):
                if dependency != "build":
                    dependency_string = dependency_string + " src$(DIR_SEP){}.$(OBJEXT)".format(dependency)
                else:
                    dependency_string = dependency_string + " src$(DIR_SEP)$(BUILD).$(OBJEXT)"
            
            output_handler.write(dependency_string)
            for include in sorted(depstruct.includes):
                output_handler.write(" src$(DIR_SEP){}".format(include))
            output_handler.write(" {}$(DIR_SEP){}.f90\n".format(directory, depstruct.name))
            output_handler.write("\t$(FC) $(OFLAG) $@ $(FFLAGS) {} {}$(DIR_SEP){}.f90\n\n".format(dependency_string.strip(), directory, depstruct.name))
            
            # For programs starting with `test_`, write code to run test too.
            if depstruct.name.startswith("test_") and ((depstruct.name[5:] in depstruct.dependencies) or (depstruct.name in extra_testnml)):
                output_handler.write("{}.nml: {}$(BINEXT)\n".format(depstruct.name[5:], depstruct.name))
                output_handler.write("\t$(RUN){}$(BINEXT)\n\n".format(depstruct.name))
                test_nmls.append(depstruct.name[5:]+".nml")

with open(before_file, "w") as output_handler:
    output_handler.write("# Automatically generated by depends.py.\n\n")
    
    # Generate `TESTNML`
    print("Writing TESTNML...")
    output_handler.write("##################\n")
    output_handler.write("# Test namelists #\n")
    output_handler.write("##################\n\nTESTNML =")
    for test_nml in sorted(test_nmls):
        output_handler.write(" {}".format(test_nml))
    
    output_handler.write("\n.PRECIOUS: $(TESTNML)\n\n")
    
    print("Writing ALLSRC...")
    output_handler.write("####################\n")
    output_handler.write("# All source files #\n")
    output_handler.write("####################\n\nALLSRC =")
    for source in sorted(all_source):
        directory = os.path.split(source)[0]
        filename  = os.path.split(source)[1]
        
        output_handler.write(" {}$(DIR_SEP){}".format(directory, filename))
