#!/usr/bin/env bash
set -e

# Shell wrapper to run the build script without installing the module
# You still need to pass a bunch of arguments
PYTHONPATH=src python -m hgvmbuilder.build "${@}"

