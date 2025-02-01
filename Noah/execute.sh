#!/bin/bash

cd ./C-Code
make

cd ../Python-Code

python3 plotting_simp_pend.py

python3 plotting_doub_poincare.py
python3 plotting_doub_pend.py
python3 animating_doub_pend.py

python3 plotting_trip_pend.py
python3 animating_trip_pend.py
python3 plotting_trip_chaos.py 0
python3 animating_trip_chaos.py 0
python3 plotting_trip_chaos.py 1
python3 animating_trip_chaos.py 1
python3 plotting_trip_compare.py

python3 plotting_num_deviation.py

cd ..

