#!/usr/bin/env python3

import sys

from functions import plot


def main():
    if len(sys.argv) < 2:
        print("usage: python synviser.py config.toml [config2.toml ...]")
        sys.exit(1)
    for config_path in sys.argv[1:]:
        try:
            plot(config_path)
        except Exception as exc:
            print(f"ERROR in {config_path}: {exc}", file=sys.stderr)
    print("that's it!")



if __name__ == "__main__":
    main()
