#!/usr/bin/env bash

sphinx-apidoc -o ./source/ ././../src/
make html
