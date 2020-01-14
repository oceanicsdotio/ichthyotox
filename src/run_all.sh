make clean
make

./ichthyotox ichthyotox 100 A
./ichthyotox ichthyotox 101 B
./ichthyotox ichthyotox 102 C
./ichthyotox ichthyotox 103 D

cd Python\ Viz\ Scripts/
python2.7 cyano_ratio_screen.py
python2.7 cyano_growth_screen.py
python2.7 cyano_paths_screen.py
python2.7 cyano_toxin_screen.py
python2.7 dissolved_toxin_screen.py

python2.7 fish_toxin_screen.py
python2.7 fish_growth_screen.py
python2.7 fish_cues_screen.py
python2.7 fish_singlecues_screen.py
python2.7 fish_singlepaths_screen.py
python2.7 fish_msd_screen.py

cd ..