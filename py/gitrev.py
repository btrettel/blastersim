#!/usr/bin/env -S python3 -Werror
# -*- coding: utf-8 -*-

import os
from subprocess import run, PIPE

with open(os.path.join("src", "rev.f90"), "w") as f90_output_handler:
    f90_output_handler.write("! auto-generated\n")
    f90_output_handler.write("module rev\n")
    f90_output_handler.write("character(len=*), public, parameter :: REVISION = \"")
    result = run(["git", "rev-parse", "--short", "HEAD"], stdout=PIPE, stderr=PIPE, universal_newlines=True)
    revision = result.stdout.strip()
    f90_output_handler.write(revision)
    f90_output_handler.write("\"\n")
    
    f90_output_handler.write("character(len=*), public, parameter :: REVISION_DATE = \"")
    result = run(["git", "show", "--no-patch", "--format=%ci", "HEAD"], stdout=PIPE, stderr=PIPE, universal_newlines=True)
    f90_output_handler.write(result.stdout.strip())
    f90_output_handler.write("\"\n")
    
    f90_output_handler.write("logical, public, parameter :: MODIFIED = ")
    result = run(["git", "diff", "--exit-code", "--quiet"], stdout=PIPE, stderr=PIPE, universal_newlines=True)
    if result.returncode == 0:
        f90_output_handler.write(".false.\n")
    else:
        f90_output_handler.write(".true.\n")
    
    # Use Git tags to put version numbers in the program.
    # Example: `git tag -a v0.1.0 -m "version 0.1.0"`
    f90_output_handler.write("character(len=*), public, parameter :: TAG = \"")
    result = run(["git", "describe", "--tags"], stdout=PIPE, stderr=PIPE, universal_newlines=True)
    tag = result.stdout.strip() # This would be `v0.1.0` in the example above.
    if result.returncode != 0:
        tag = revision
    f90_output_handler.write(tag)
    f90_output_handler.write("\"\n")
    f90_output_handler.write("end module rev")
    
    with open(os.path.join("rev.tex"), "w") as tex_output_handler:
        tex_output_handler.write("% auto-generated\n")
        tex_output_handler.write("% chktex-file 8\n")
        tex_output_handler.write(r"\newcommand*{\gittag}{"+tag+"}")
