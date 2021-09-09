#! /usr/bin/env python
import sys

modules = ["cgnsAPI", "adtAPI", "curveSearchAPI", "intersectionAPI", "utilitiesAPI"]

for name in modules:
    print("")
    print(f"Testing if module {name} can be imported...")
    import_cmd = f"import {name}"
    try:
        exec(import_cmd)
    except ImportError:
        print(f"Error importing {name}")
        sys.exit(1)

    print(f"Module {name} was successfully imported.")
