#!/usr/bin/env python3

"""

Program     : SynViser, synteny visualisation tool.
Language    : Python
Authors     : Yehor Tertyshnyk, Mirela Krol
Email       : egor.tertyshnyk@gmail.com
Github      : https://github.com/Nepomorpha/SynViser
License     : MIT License
Citation    : -TBS-

Description : This program is designed to take various inputs from config files and produce synteny plots.

Usage:         python synviser.py config.toml

Instructions: see README.md or dedicated .ipynb notebooks for more details.

"""

import sys

from functions import plot


def main():
    if len(sys.argv) < 2:
        print("usage: python synviser.py config.toml [config2.toml ...]")
        sys.exit(1)
    failed = False 
    for config_path in sys.argv[1:]:
        try:
            plot(config_path)
        except Exception as exc:
            failed = True
            print(f"error in {config_path}: {exc}",file=sys.stderr)

    if failed:
        sys.exit(1)

    print("plotting is done!")



if __name__ == "__main__":
    main()
