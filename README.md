# genstudies-eft

Clone repository:
```sh
git clone --recurse-submodules https://github.com/fstaeg/genstudies-eft.git
```

Setup (should be sourced at the start of each new session):
```sh
. /cvmfs/sft.cern.ch/lcg/views/LCG_106/x86_64-el8-gcc11-opt/setup.sh
```

## Workflow

To generate the samples with SMEFT weights included, we use the `EFT2Obs` tool, which we include in this repository as a submodule.

### Installing EFT2Obs

```sh
. /cvmfs/sft.cern.ch/lcg/views/LCG_106/x86_64-el8-gcc11-opt/setup.sh

cd EFT2Obs
source env.sh # should be sourced at the start of each new session

# Installing MG5_aMC@NLO and Pythia8
./scripts/setup_mg5.sh

# Installing the SMEFTsim model
./scripts/setup_model_SMEFTsim3.sh
```

### Setting up the event generation

1) First we create a subdirectory `cards/DYjets-SMEFTsim3` with a file `proc_card.dat`:
```
import model SMEFTsim_topU3l_MwScheme_UFO-massless

# Definitions
define ell+ = e+ mu+
define ell- = e- mu-

generate p p > ell+ ell- NP<=1 @0
add process p p > ell+ ell- j NP<=1 @1

output DYjets-SMEFTsim3 -nojpeg
```

- first we import the SMEFTsim3 model
- then we generate pp $\rightarrow$ $\ell^{+} \ell^{-}(+j)$ events with up to one new physics vertices inserted (`NP<=1`)
- the directory specified in the last line has to match the name of the subdirectory in `cards` (`DYjets-SEMFTsim3`)

2) Initialise the process in MadGraph:
```sh
./scripts/setup_process.sh DYjets-SMEFTsim3
```

- this creates a process directory `MG5_aMC-v2_9_16/DYjets-SMEFTsim3`
- take a look at the Feynman diagrams in `MG5_aMC-v2_9_16/DYjets-SMEFTsim3/SubProcesses/*/matrix*.ps`

3) The last step created two new cards in `cards/DYjets-SMEFTsim3`: `pythia8_card.dat` and `run_card.dat`. 

- We make three changes to `run_card.dat`:
```txt
10.0 = xqcut
False = use_syst
none = systematics_program
```

- In `pythia8_card.dat`, we turn off multi-parton interactions to speed things up:
```txt
partonlevel:mpi = off
```

- in `run_card.dat` we can also apply generator level cuts. Later we will probably have to generate events in several bins of mll, to improve the statistics of our samples at high invariant mass. For this we can use the parameters `mmll` and `mmllmax`.

4) Next we create a configuration file that specifies the SMEFT operators we want to study. It is generated with:
```sh
python scripts/make_config.py -p DYjets-SMEFTsim3 \
 -o cards/DYjets-SMEFTsim3/config.json --pars SMEFT:108,109,113,115,117,119,121 \
 --def-val 1.0 --def-sm 0.0 --def-gen 0.0
```

- 108,...,121 are the IDs of the seven operators we are interested in: `clj1, clj3, ceu, ced, cje, clu, cld`

5) Next we use the `config.json` to create two additional cards, `param_card.dat` and `reweight_card.dat`:
```sh
python scripts/make_param_card.py -p DYjets-SMEFTsim3 \
 -c cards/DYjets-SMEFTsim3/config.json -o cards/DYjets-SMEFTsim3/param_card.dat

python scripts/make_reweight_card.py cards/DYjets-SMEFTsim3/config.json \
 cards/DYjets-SMEFTsim3/reweight_card.dat
```

### Generating events

1) Now that all the cards are prepared, we can produce a gridpack:
```sh
./scripts/make_gridpack.sh DYjets-SMEFTsim3
```

2) From the gridpack, we can now generate events:
```sh
python scripts/launch_jobs.py --gridpack gridpack_DYjets-SMEFTsim3.tar.gz \
 -j 10 -e 10000 -s 1 --job-mode interactive --parallel 10 \
 --to-step shower --save-hepmc DYjets-SMEFTsim3
```

- this will generate a total of 100k events split in 10 jobs that run in parallel
- in the last line, we set the option to stop after the showering and creating the HepMC files. By default it would run a Rivet routine.

3) We are done. The generated HepMC files are saved in `DYjets-SMEFTsim3/events_{N}.hepmc.gz`. Move them to the `genstudies-eft` directory:
```sh
cp -r DYjets-SMEFTsim3 ..
cd ..
```

### Processing the samples

0) Install the python packages needed to read HepMC files and run the jet clustering:
```sh
pip install pyhepmc fastjet
```

1) Create a ROOT tree from the HepMC files:
```sh
python hepmc_to_tree.py -d DYjets-SMEFTsim3
# options
# -d {dirname}: directory of input HepMC files and for output ROOT file
# -f {filename.root}: name of the output ROOT file (default: events.root)
# -n {N}: only read the first N input files. Useful to save time when testing code (default: all of them)
```

2) Make plots:
```sh
python make_plots.py -i DYjets-SMEFTsim3/events.root -o plots
```
