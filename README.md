# clas12DVCS
A python code set to analyze CLAS12 DVCS data.

root2pickle assumes to take root file converted by 
[https://github.com/Sangbaek/convertingHipo](https://github.com/Sangbaek/convertingHipo).

## How to use in ifarm

```
module load python/3.9.1
python3 -m pip install numpy
python3 -m pip install pandas
python3 -m pip install uproot
```
## Usage of scripts
```
python3 root2pickleEpggExp.py -f "/path/to/root/file" -o "output_name.pkl" -s "10000"
```
The optional s flag reads the number of entries to be read. If unused, it converts whole root file. 

For the EpgRec script,

```
python3 root2pickleEpgRec.py -f "/path/to/root/file" -o "output_name.pkl" -s "10000" -gen dvcs
```
we can select the generator between ```dvcs``` and ```pi0```.

## Usage of scripts, MC data

EpggRec: dvpi0 -> (dvpi0 -> 2 gamma)

EpgRec: dvpi0 -> (dvpi0 -> one gamma)

EpgRec: dvcs -> dvcs only

## Usage of scripts, exp data

EpggExp: dvpi0 candidates -> (dvpi0 -> 2 gamma)

EpgExp: dvpi0 candidates -> (dvpi0 -> one gamma)

EpgExp: dvcs candidates -> dvcs + (dvpi0 -> one gamma)