#!/usr/bin/env -S python3 -Werror
# -*- coding: utf-8 -*-

import os
import argparse
import configparser

# Branch statement density is a sort of poor man's cyclomatic complexity. But it's better in one sense, as I've heard that ["the McCabe cyclomatic complexity metric is no better than counting lines of code"](https://shape-of-code.com/2023/09/03/halstead-mccabe-metrics-the-wisdom-of-the-ancients/). This looks at the *density*, so it's not the same information as the total lines of code.

ASSERTION_DENSITY_LOWER_LIMIT = 5.0  # percent
BRANCH_DENSITY_UPPER_LIMIT    = 10.0 # percent (number is arbitrary for the moment)
COMMENT_DENSITY_LOWER_LIMIT   = 5.0 # percent
TEST_RATIO_LOWER_LIMIT        = 0.5

class file_stats_class:
    def __init__(self, filename, assertion_density, branch_density, comment_density, lines, test_subroutine_calls, test_subroutine_definitions):
        self.filename                    = filename
        self.assertion_density           = assertion_density
        self.branch_density              = branch_density
        self.comment_density             = comment_density
        self.lines                       = lines
        self.test_ratio                  = 0.0
        self.test_subroutine_calls       = test_subroutine_calls
        self.test_subroutine_definitions = test_subroutine_definitions

def is_empty_or_comment(line):
    line_clean = line.split("!")[0].strip()
    
    return not len(line_clean) > 0

def canonicalize_path(path):
    path_split = path.split("/")
    return os.path.join(*path_split)

def canonicalize_paths(path_string):
    path_string_split = path_string.split(' ')
    paths = []
    for path_string_split_i in path_string_split:
        paths.append(canonicalize_path(path_string_split_i))
    return paths

parser = argparse.ArgumentParser(description="A linter for Fortran 90 and later intended for Ben Trettel's use.")
parser.add_argument("file", help="input file to read")
args = parser.parse_args()

config = configparser.ConfigParser()
config.read(args.file)

directories              = config['f90lint']['directories'].split(' ')
skip_indexing            = canonicalize_paths(config['f90lint']['skip_indexing'])
ignore_assertion_density = canonicalize_paths(config['f90lint']['ignore_assertion_density'])
ignore_test_sloc_ratio   = canonicalize_paths(config['f90lint']['ignore_test_sloc_ratio'])
ignore_branch_density    = canonicalize_paths(config['f90lint']['ignore_branch_density'])
ignore_comment_density   = canonicalize_paths(config['f90lint']['ignore_comment_density'])

fail = False

filepaths = []
for directory in sorted(directories):
    # Canonicalize the filepath, so that (for example), this works if `.\` is in front of the path, like PowerShell does.
    directory = os.path.relpath(directory)
    
    if not os.path.isdir(directory):
        print("ERROR: {} in directories is not a directory.".format(directory))
        fail = True
    else:
        for filename in os.listdir(directory):
            # Exclude symlinks as those presumably are linted elsewhere.
            if os.path.isfile(os.path.join(directory, filename)) and filename.endswith(".f90"):
                filepaths.append(os.path.join(directory, filename))

assert(len(filepaths) > 0)

if fail:
    print("Error(s) encountered, stopping.")
    exit(1)

exit_code = 0

global_num_lines_code        = 0
global_num_lines_tests       = 0
global_num_assertions        = 0
global_num_branch_statements = 0
global_num_comments          = 0
num_files_analyzed           = 0
file_stats = []
for filename in sorted(filepaths):
    if filename in skip_indexing:
        continue
    
    num_files_analyzed = num_files_analyzed + 1
    
    with open(filename, "r") as file_handler:
        contains_executable_code = False
        for line in file_handler.readlines():
            if line.strip().startswith("program ") or line.strip().startswith("contains"):
                contains_executable_code = True
                break
        
        file_handler.seek(0)
        
        # Ignore files that don't contain executable code.
        if contains_executable_code:
            file_handler.seek(0)
            local_num_lines_code_or_tests = 0
            local_num_lines_tests         = 0
            local_num_assertions          = 0
            local_num_branch_statements   = 0
            local_num_comments            = 0
            line_no                       = 0
            test_subroutine_calls         = set()
            test_subroutine_definitions   = set()
            for line in file_handler.readlines():
                line_no = line_no + 1
                
                # If the line isn't empty or only comments.
                if is_empty_or_comment(line):
                    if "!" in line.strip():
                        global_num_comments = global_num_comments + 1
                        local_num_comments  = local_num_comments  + 1
                else:
                    local_num_lines_code_or_tests = local_num_lines_code_or_tests + 1
                    
                    line_no_comments = line.split("!")[0].strip()
                    
                    # Tests don't have assertions, so they shouldn't count towards the assertion count.
                    if not filename.startswith("test"):
                        global_num_lines_code = global_num_lines_code + 1
                        
                        # `error stop` lines are for error termination, so they are like assertions. They are used where assertions aren't appropriate.
                        if (line.strip().startswith("call assert") or line.strip().startswith("error stop") or ("%check(" in line)):
                            global_num_assertions = global_num_assertions + 1
                            local_num_assertions  = local_num_assertions  + 1
                        
                        if (line.strip().startswith("if ") or line.strip().startswith("select case") or ("merge(" in line)):
                            global_num_branch_statements = global_num_branch_statements + 1
                            local_num_branch_statements  = local_num_branch_statements  + 1
                        
                        # Normal `do` loops return an error unless they have the pragma `! SERIAL` on the same line. The point of this is to encourage me to use `do concurrent` as much as possible to make sure that the code can be parallelized well.
                        if ":" in line_no_comments:
                            line_do = line_no_comments.split(":")[1].strip()
                        else:
                            line_do = line_no_comments
                        
                        if line_do.startswith("do ") or (line_do == "do"):
                            if (not line_do.startswith("do concurrent ")) and (not "! SERIAL" in line):
                                print("{}:{}: do concurrent not used". format(filename, line_no))
                                exit_code = 1
                        
                        # # Require `function`s to be marked as `pure`, `impure`, or `elemental`.
                        # if line_no_comments.startswith("function "):
                            # print("{}:{}: functions must be pure or elemental". format(filename, line_no))
                            # exit_code = 1
                        
                        # # Require `subroutines`s to be marked as `pure`, `impure`, or `elemental`.
                        # if line_no_comments.startswith("subroutine "):
                            # print("{}:{}: subroutines must be marked as pure, elemental, or impure". format(filename, line_no))
                            # exit_code = 1
                        
                        if line_no_comments == "end type":
                            print("{}:{}: end type lacks name". format(filename, line_no))
                            exit_code = 1
                        
                        # Not all `interface` blocks have names.
                        # if line_no_comments == "end interface":
                            # print("{}:{}: end interface lacks name". format(filename, line_no))
                            # exit_code = 1
                        
                        if line_no_comments == "end function":
                            print("{}:{}: end function lacks name". format(filename, line_no))
                            exit_code = 1
                        
                        if line_no_comments == "end subroutine":
                            print("{}:{}: end subroutine lacks name". format(filename, line_no))
                            exit_code = 1
                    else:
                        global_num_lines_tests = global_num_lines_tests + 1
                        
                        if line_no_comments.startswith("call test_") and not "%" in line_no_comments:
                            start_index = line.find(" test_") + 1
                            end_index   = line.find("(")
                            test_subroutine_calls.add(line[start_index:end_index])
                        elif line_no_comments.startswith("subroutine test_"):
                            start_index = line.find(" test_") + 1
                            end_index   = line.find("(")
                            test_subroutine_definitions.add(line[start_index:end_index])
            
            file_stats.append(file_stats_class(
                                filename, 
                                100.0 * local_num_assertions / local_num_lines_code_or_tests,
                                100.0 * local_num_branch_statements / local_num_lines_code_or_tests,
                                100.0 * local_num_comments / (local_num_lines_code_or_tests + local_num_comments),
                                local_num_lines_code_or_tests,
                                test_subroutine_calls,
                                test_subroutine_definitions))
            assert(global_num_assertions < global_num_lines_code)
            assert(global_num_branch_statements < global_num_lines_code)

assert(num_files_analyzed > 0)
assert(global_num_lines_code > 0)

# Calculate `test_ratio`s.
for file_stat in file_stats:
    if not file_stat.filename.startswith("test"):
        for sub_file_stat in file_stats:
            if sub_file_stat.filename.startswith("test"):
                if os.path.basename(sub_file_stat.filename) == "test_" + os.path.basename(file_stat.filename):
                    file_stat.test_ratio = sub_file_stat.lines / file_stat.lines
                    #print(file_stat.filename, sub_file_stat.filename, file_stat.test_ratio)
                    break

for file_stat in file_stats:
    if (file_stat.assertion_density < ASSERTION_DENSITY_LOWER_LIMIT) and (not file_stat.filename.startswith("test")) and (not file_stat.filename in ignore_assertion_density):
        print("{}: Assertion density is {:.2f}% (<{:.2f}%).".format(file_stat.filename, file_stat.assertion_density, ASSERTION_DENSITY_LOWER_LIMIT))
        exit_code = 1
    
    if (file_stat.branch_density > BRANCH_DENSITY_UPPER_LIMIT) and (not file_stat.filename.startswith("test")) and (not file_stat.filename in ignore_branch_density):
        print("{}: Branch density is {:.2f}% (>{:.2f}%).".format(file_stat.filename, file_stat.branch_density, BRANCH_DENSITY_UPPER_LIMIT))
        exit_code = 1
    
    if (file_stat.comment_density < COMMENT_DENSITY_LOWER_LIMIT) and (not file_stat.filename.startswith("test")) and (not file_stat.filename in ignore_comment_density):
        print("{}: Comment density is {:.2f}% (<{:.2f}%).".format(file_stat.filename, file_stat.comment_density, COMMENT_DENSITY_LOWER_LIMIT))
        exit_code = 1
    
    if (file_stat.test_ratio < TEST_RATIO_LOWER_LIMIT) and file_stat.filename.startswith("src") and (not file_stat.filename in ignore_test_sloc_ratio):
        print("{}: Lines-of-tests to lines-of-code ratio is {:.2f} (<{:.2f}).".format(file_stat.filename, file_stat.test_ratio, TEST_RATIO_LOWER_LIMIT))
        exit_code = 1

global_assertion_density = 100.0 * global_num_assertions / global_num_lines_code
if (global_assertion_density < ASSERTION_DENSITY_LOWER_LIMIT):
    print("Global assertion density is {:.2f}% (<{:.2f}%).".format(global_assertion_density, ASSERTION_DENSITY_LOWER_LIMIT))
    exit_code = 1

global_branch_density = 100.0 * global_num_branch_statements / global_num_lines_code
if (global_branch_density > BRANCH_DENSITY_UPPER_LIMIT):
    print("Global branch density is {:.2f}% (>{:.2f}%).".format(global_branch_density, BRANCH_DENSITY_UPPER_LIMIT))
    exit_code = 1

global_comment_density = 100.0 * global_num_comments / (global_num_lines_code + global_num_comments)
if (global_comment_density < COMMENT_DENSITY_LOWER_LIMIT):
    print("Global comment density is {:.2f}% (<{:.2f}%).".format(global_comment_density, COMMENT_DENSITY_LOWER_LIMIT))
    exit_code = 1

global_test_ratio = global_num_lines_tests / global_num_lines_code
if (global_test_ratio < TEST_RATIO_LOWER_LIMIT):
    print("Global lines-of-tests to lines-of-code ratio is {:.2f} (<{:.2f}).".format(global_test_ratio, TEST_RATIO_LOWER_LIMIT))
    exit_code = 1

for file_stat in file_stats:
    if file_stat.filename.startswith("test"):
        if len(file_stat.test_subroutine_definitions - file_stat.test_subroutine_calls) > 0:
            print("{}: Missing test subroutine call(s): {}".format(file_stat.filename, file_stat.test_subroutine_definitions - file_stat.test_subroutine_calls))
            exit_code = 1
        elif len(file_stat.test_subroutine_calls - file_stat.test_subroutine_definitions) > 0:
            print("{}: Missing test subroutine definition(s): {}".format(file_stat.filename, file_stat.test_subroutine_calls - file_stat.test_subroutine_definitions))
            exit_code = 1

print("\n==============")
print("= Statistics =")
print("==============")
print("Number of source files analyzed: {}".format(num_files_analyzed))
print("Global assertion density: {:.2f}%".format(global_assertion_density))
print("Global branch density: {:.2f}%".format(global_branch_density))
print("Global comment density: {:.2f}%".format(global_comment_density))
print("Global lines-of-tests to lines-of-code ratio: {:.2f}".format(global_test_ratio))

exit(exit_code)
