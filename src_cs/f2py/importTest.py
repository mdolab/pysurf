#! /usr/bin/env python
import sys

modules = ["cgnsAPI", "adtAPI", "curveSearchAPI", "intersectionAPI", "utilitiesAPI"]
modules_CS = [name + "_cs" for name in modules]

for name in modules_CS:
    print("")
    print(f"Testing if module {name} can be imported...")
    import_cmd = f"import {name}"
    try:
        exec(import_cmd)
    except ImportError:
        print(f"Error importing {name}")
        sys.exit(1)

    print(f"Module {name} was successfully imported.")
