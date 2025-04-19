#!/usr/bin/env python3

import json
import sys

def main():
    if len(sys.argv) < 2:
        print("Usage: parse_json.py <path-to-config.json>")
        sys.exit(1)

    json_path = sys.argv[1]
    with open(json_path, "r") as f:
        data = json.load(f)

    compile_time_data = data.get("compileTime", {})

    for key, val in compile_time_data.items():
        # Define macro with no value only for True.
        if val is True:
            # Print just the macro name (no "=...")
            print(key)
        # Skip if the value is False
        elif val is False:
            continue
        elif isinstance(val, str):
            # Strings get quotes
            print(f'{key}="{val}"')
        elif isinstance(val, (int, float)):
            # Numeric values
            print(f"{key}={val}")
        else:
            # If it's null, false, or something else, handle as you prefer
            # e.g., define but no value
            print(key)

if __name__ == "__main__":
    main()
