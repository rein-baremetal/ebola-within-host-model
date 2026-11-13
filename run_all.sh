#!/usr/bin/env bash
# Rebuilds engine, figures, checks and the .docx from this folder (about 8 minutes). Requires g++, python3 + numpy + PIL, rsvg-convert.
set -e
cd "$(dirname "$0")"
g++ -std=c++20 -O2 -pthread -Wall -Wextra code_annotated.cpp -o engine
g++ -std=c++20 -O2 rpt_adapt.cpp -o rpt_adapt
./engine --selftest input.txt
bash check_determinism.sh
bash check_sanitisers.sh
python3 verify_threshold.py input.txt
./engine input.txt sample_trajectory.txt errors.txt
sed 's/^tip_dose *=.*/tip_dose = 0/' input.txt > untreated_input.txt
./engine untreated_input.txt untreated_output.txt errors.txt
./engine --screen screen.txt input.txt errors.txt
python3 verify_figures.py sample_trajectory.svg sample_trajectory.txt untreated_output.txt
python3 verify_figures.py --pareto screen_pareto.svg screen_samples.csv
python3 make_figures.py .
python3 make_figures_analysis.py .
python3 build_docx_figures.py
python3 claims_check.py
