#!/usr/bin/env -S python3 -Werror
# -*- coding: utf-8 -*-

import xml.etree.ElementTree as ET
import sys
import argparse
from subprocess import run, PIPE

#rule_ids_to_ignore = set()
# COM.PRES.Indent: I don't indent in their expected style.
# COM.PRES.LengthLine: My maximum line length is 132, not 100.
# COM.INST.CodeComment: Commented out code can be useful sometimes.
# COM.DATA.NotUsed: The compilers seem to do a better check on these with no false positives.
# COM.FLOW.Abort: The `stop` procedure is useful.
# COM.DATA.Invariant: Too many false positives.
# COM.DATA.Initialisation: Too many false positives. ("The variable parameter is used before being initialized". No, the variable name is not parameter!)
# F90.REF.Interface: Too many false positives.
# F90.DATA.Float: Too many false positives.
# F90.DATA.Declaration: The compilers seem to do a better check on these with no false positives.
# COM.INST.Brace: Not sure how to satisfy this one.
# COM.FLOW.FileExistence: False positive when this is done with a procedure.
# F90.INST.Only: Why have this?
# F90.DESIGN.IO: Too many false positives.
# COM.PROJECT.Header: You need to add comments precisely how they expect it, which is quite brittle. Too many false positives.
# F90.DESIGN.Interface: It's not clear to me why this warning is appearing. Could be due to using derived types.
# COM.FLOW.CheckCodeReturn: False positives on derived type constructors.
# F90.NAME.KeyWords: Too many false positives.
# F90.TYPE.Integer: Too many false positives.
# COM.FLOW.ExitLoop: There are reasons to have more than one exit per loop.
# COM.FLOW.CheckUser: I have no idea what this is supposed to catch.
# F90.REF.Open, F90.ERR.OpenRead, F90.ERR.Allocate: I decided to not use `iostat`/`stat` in `open`, `allocate`, etc., if I simply want the compiler to terminate the program at that point.
# COM.DESIGN.Alloc: This can't be satisfied for my operators in autodiff.f90.
# F90.INST.Nullify: This leads to false positives for allocated pointer arrays.
# COM.FLOW.Exit: Too many false positives.
# F90.REF.Label: False positives for `end where`, and false positives in general.
# COM.TYPE.Expression: Too many false positives.
# F90.DATA.Constant: I don't care about this one.
# COM.NAME.Homonymy: Too many false positives.
# F90.BLOC.File: Disabled as I often close units in different scopes.
# COM.FLOW.FilePath: I don't care about this one.
# F90.DATA.ConstantFloat: False positives on arrays.
# F90.INST.Intent: False positives for automatic arrays.
# COM.DATA.FloatCompare: False positives...
# F90.DESIGN.Include: I think that `include` is okay.
# COM.FLOW.BooleanExpression: I could see how this is bad, but simple versions should be okay in my view.
# F90.REF.Variable: Seems to have false positives. I'm also not sure the idea of this is a good one.
# F90.TYPE.Real: False positive in genunits_io.f90.
rule_ids_to_ignore = {"COM.PRES.Indent", "COM.PRES.LengthLine", "COM.INST.CodeComment", "COM.DATA.NotUsed", "COM.FLOW.Abort", "COM.DATA.Invariant", "COM.DATA.Initialisation", "F90.REF.Interface", "F90.DATA.Float", "F90.DATA.Declaration", "COM.INST.Brace", "COM.FLOW.FileExistence", "F90.INST.Only", "F90.DESIGN.IO", "COM.PROJECT.Header", "F90.DESIGN.Interface", "COM.FLOW.CheckCodeReturn", "F90.NAME.KeyWords", "F90.TYPE.Integer", "COM.FLOW.ExitLoop", "COM.FLOW.CheckUser", "F90.REF.Open", "F90.ERR.OpenRead", "F90.ERR.Allocate", "COM.DESIGN.Alloc", "F90.INST.Nullify", "COM.FLOW.Exit", "F90.REF.Label", "COM.TYPE.Expression", "F90.DATA.Constant", "COM.NAME.Homonymy", "F90.BLOC.File", "COM.FLOW.FilePath", "F90.DATA.ConstantFloat", "F90.INST.Intent", "COM.DATA.FloatCompare", "F90.DESIGN.Include", "COM.FLOW.BooleanExpression", "F90.REF.Variable", "F90.TYPE.Real"}

# "Variable names should be unique. The variable public is already defined in this file.": False positives when types are declared public.
# "Variable names should be unique. The variable I10 is already defined in this file.": False positives.
#rule_messages_to_ignore = set()
rule_messages_to_ignore = {"Internal i-Code error", "The return code of the function i is not checked", "It misses the declaration SELECTED_REAL_KIND in the initialisation of reals_csv", "It misses the declaration SELECTED_REAL_KIND in the initialisation of headers_csv", "Variable names should be unique. The variable public is already defined in this file.", "Variable names should be unique. The variable I10 is already defined in this file."}
#rule_messages_to_ignore = {'Return code used in arithmetical statement", "The variable selected_real_kind is a keyword in Fortran90 language", "The sequence IMPLICIT NONE must be declared after the method", "The variable public is used before being initialized", "The variable outer must be declared", "The variable inner must be declared", "The return code of the function i is not checked", "The return code of the function lower_index is not checked", "The return code of the function upper_index is not checked", "The variable selected_int_kind is a keyword in Fortran90 language", "The variable check_flag must be declared", "The existences of the file error_unit must be checked'}

parser = argparse.ArgumentParser(description="A wrapper for i-Code CNES (a Fortran linter) to output plain text with meaningful exit codes.")
parser.add_argument('file', nargs='+')
args = parser.parse_args()

command = ['icode']
command.extend(args.file)

result = run(command, stdout=PIPE, stderr=PIPE, universal_newlines=True)

root = ET.fromstring(result.stderr)

warnings_found = 0

for analysisRules in root.findall('analysisRule'):
    analysisRuleId = analysisRules.get('analysisRuleId')
    
    if analysisRuleId in rule_ids_to_ignore:
        continue
    
    for results in analysisRules:
        filename = results.get('fileName')
        line_number = results.get('resultLine')
        display = True
        for resultMessage in results:
            for rule_message_to_ignore in rule_messages_to_ignore:
                if rule_message_to_ignore in resultMessage.text:
                    display = False
            
            if display:
                print(f"{filename}:{line_number} {resultMessage.text} ({analysisRuleId})")
                warnings_found += 1

for f90_file in args.file:
    command = ["flint", "lint", "--flintrc", "f90.yaml", f90_file]
    
    result = run(command, stdout=PIPE, stderr=PIPE, universal_newlines=True)
    
    if len(result.stdout) > 0:
        warnings_found = warnings_found + 1
        print(f90_file + ":")
        print(result.stdout)

if warnings_found > 0:
    sys.exit(1)
