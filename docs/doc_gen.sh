#!/usr/bin/env bash

sphinx-apidoc -o . ././../src/
make html
