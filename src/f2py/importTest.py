#!/usr/bin/env python
# Standard Python modules
import argparse
import os
import sys

parser = argparse.ArgumentParser()
parser.add_argument(
    "name",
    nargs="+",
    help="One or more library names (example: libpackage.so). Note: This script must be run in the same dir as the library.",
)
args = parser.parse_args()

# Make sure the executing directory is always in the PATH before importing
sys.path.insert(0, os.getcwd())

for libraryName in args.name:
    # Only get the filename without the extension
    moduleName = os.path.splitext(libraryName)[0]
    print(f"Testing if module {moduleName} can be imported...")

    try:
        import_cmd = f"import {moduleName}"
        exec(import_cmd)
    except ImportError as e:
        print(f"Error: {e}")
        print(f"Error: library {libraryName} was not imported correctly")
        sys.exit(1)

    print(f"Module {moduleName} was successfully imported")
