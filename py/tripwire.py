#!/usr/bin/env -S python3 -Werror
# -*- coding: utf-8 -*-

# Run tests with `python3 -m unittest tripwire.py`.

import unittest

TRIPWIRE_DIRECTIVE_PREFIX = "tripwire$"
TRIPWIRE_BEGIN_DIRECTIVE  = TRIPWIRE_DIRECTIVE_PREFIX+" begin"
TRIPWIRE_END_DIRECTIVE    = TRIPWIRE_DIRECTIVE_PREFIX+" end"

def crc32(x):
    from zlib import crc32
    
    checksum_integer = crc32(bytes(x, "utf-8"))
    
    # <https://stackoverflow.com/a/12638477/1124489>
    checksum = f"{checksum_integer:0{8}x}"
    assert len(checksum) == 8
    
    return checksum.upper()

def get_fenced_lines(file, checksum):
    fenced_lines = ""
    
    with open(file, "r") as file_in:
        next_line_in_fence = False
        line_no = 0
        for line in file_in:
            line_no += 1
            
            if next_line_in_fence:
                if TRIPWIRE_END_DIRECTIVE in line:
                    break
                else:
                    fenced_lines = fenced_lines + line
            
            if TRIPWIRE_BEGIN_DIRECTIVE in line:
                line_checksum, message = parse_begin_directive(get_tripwrite_directive(line))
                begin_line_no = line_no
                if checksum == line_checksum:
                    next_line_in_fence = True
    
    return fenced_lines, message, begin_line_no

def get_reported_checksums(file):
    reported_checksums = set()
    duplicate_checksum_line_nos = set()
    
    with open(file, "r") as file_in:
        line_no = 0
        for line in file_in:
            line_no += 1
            
            if TRIPWIRE_BEGIN_DIRECTIVE in line:
                reported_checksum, _ = parse_begin_directive(get_tripwrite_directive(line))
                
                if (reported_checksum in reported_checksums):
                    duplicate_checksum_line_nos.add(line_no)
                
                reported_checksums.add(reported_checksum)
    
    return reported_checksums, duplicate_checksum_line_nos

def get_tripwrite_directive(line):
    assert TRIPWIRE_DIRECTIVE_PREFIX in line
    
    return line.strip()[line.strip().index(TRIPWIRE_DIRECTIVE_PREFIX):]

def parse_begin_directive(directive):
    assert directive.startswith(TRIPWIRE_BEGIN_DIRECTIVE)
    
    directive_split = directive.strip().split(" ")
    assert directive_split[1] == "begin"
    
    checksum = directive_split[2]
    message  = " ".join(directive_split[3:]).strip()
    
    return checksum, message

if __name__ == "__main__":
    import argparse
    import sys
    
    parser = argparse.ArgumentParser(description="Detects when specified lines of code change.")
    parser.add_argument("files", nargs="*", help="input file(s) to read")
    args = parser.parse_args()
    
    exit_code = 0
    for file in args.files:
        reported_checksums, duplicate_checksum_line_nos = get_reported_checksums(file)
        
        for reported_checksum in reported_checksums:
            fenced_lines, message, line_no = get_fenced_lines(file, reported_checksum)
            actual_checksum = crc32(fenced_lines)
            
            if actual_checksum != reported_checksum:
                print("{}:{} (reported {}, actual {}): {}".format(file, line_no, reported_checksum, actual_checksum, message))
                exit_code = 1
        
        for duplicate_checksum_line_no in duplicate_checksum_line_nos:
            print("{}:{}: duplicate checksum, possibly not checked".format(file, line_no))
            exit_code = 1
    
    if (exit_code > 0):
        print("ONE OR MORE TRIPWIRE(S) TRIPPED OR ERRORS ENCOUNTERED.")
    sys.exit(exit_code)

class Tests(unittest.TestCase):
    # <https://docs.python.org/3/library/unittest.html>
    
    def test_crc32(self):
        # <https://asecuritysite.com/hash/crc32>
        self.assertEqual(crc32("The quick brown fox jumps over the lazy dog"), "414FA339")
        self.assertEqual(crc32("Test vector from febooti.com"), "0C877F61")
        self.assertEqual(crc32(""), "00000000")
        
        # <https://lxp32.github.io/docs/a-simple-example-crc32-calculation/>
        self.assertEqual(crc32("123456789"), "CBF43926")
    
    def test_get_tripwrite_directive(self):
        self.assertEqual(get_tripwrite_directive("don't include this tripwire$ begin CHECKSUM Multi-word message"), "tripwire$ begin CHECKSUM Multi-word message")
        self.assertEqual(get_tripwrite_directive("don't include this tripwire$ end"), "tripwire$ end")
    
    def test_parse_begin_directive(self):
        checksum, message = parse_begin_directive("tripwire$ begin CHECKSUM Multi-word message")
        
        self.assertEqual(checksum, "CHECKSUM")
        self.assertEqual(message, "Multi-word message")
