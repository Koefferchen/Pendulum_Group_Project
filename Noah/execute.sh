#!/bin/bash
cd ./C-Code
make
cd ../Python-Code
python3 create_plots.py
cd ..

