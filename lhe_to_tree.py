import ROOT
import os
from array import array
from argparse import ArgumentParser
ROOT.gROOT.ProcessLine('#include "EFT2Obs/LHEF.h"')

parser = ArgumentParser()
parser.add_argument('-d', '--directory', default=None, help='directory of input LHE files and for output ROOT file')
parser.add_argument('-f', '--filename', default='events.root', help='name of output ROOT file')
parser.add_argument('-n', '--nfiles', type=int, default=None, help='only run on the first n input files')
args = parser.parse_args()

# create the root file
outfile = ROOT.TFile(os.path.join(args.directory, args.filename), 'RECREATE') 

# define variables 
# arrays for variables with one entry per event, 
# vectors for variables with multiple entries per event
mll = array('d',[0.0])
yll = array('d',[0.0])
ptll = array('d',[0.0])
sm_weight = array('d',[0.0])
smeft_weight = ROOT.std.vector('double')()
l_pdgid = ROOT.std.vector('int')()
l_px = ROOT.std.vector('double')()
l_py = ROOT.std.vector('double')()
l_pz = ROOT.std.vector('double')()
l_E = ROOT.std.vector('double')()
l_pt = ROOT.std.vector('double')()
l_eta = ROOT.std.vector('double')()
l_phi = ROOT.std.vector('double')()

# create ROOT tree and add branches
# note the different syntax for branches with one or multiple entries per event
ttree = ROOT.TTree('events','tree of generated events')
ttree.Branch('mll', mll, 'mll/D')
ttree.Branch('yll', yll, 'yll/D')
ttree.Branch('ptll', ptll, 'ptll/D')
ttree.Branch('EventWeight_SM', sm_weight, 'EventWeight_SM/D')
ttree.Branch('EventWeight_SMEFT', smeft_weight)
ttree.Branch('Lepton_pdgId',l_pdgid)
ttree.Branch('Lepton_px',   l_px)
ttree.Branch('Lepton_py',   l_py)
ttree.Branch('Lepton_pz',   l_pz)
ttree.Branch('Lepton_E',    l_E)
ttree.Branch('Lepton_pt',   l_pt)
ttree.Branch('Lepton_eta',  l_eta)
ttree.Branch('Lepton_phi',  l_phi)

# get all the events_{n}.lhe.gz files in the input directory
infiles = [f for f in os.listdir(args.directory) if f.startswith('events_') and f.endswith('.lhe.gz')]
infiles = sorted(infiles, key=lambda x: float(x[7:-7])) # sort them
if args.nfiles is not None: 
    infiles = infiles[:args.nfiles]

nevent = 0

# loop over input files
for infile in infiles:
    infile = os.path.join(args.directory, infile)
    os.system('gunzip %s' % infile) # decompress the file
    infile = infile.replace('.gz','')
    reader = ROOT.LHEF.Reader(infile)

    # loop over events
    while reader.readEvent():
        if reader.hepeup.isGroup: subevents = reader.hepeup.subevents
        else: subevents = [reader.hepeup]
        
        for event in subevents:
            nevent += 1
            sm_weight[0] = reader.hepeup.weights[0].first # event weight
            # loop over the event weights from SMEFTsim
            for i in range(1,len(reader.hepeup.weights)):
                smeft_weight.push_back(reader.hepeup.weights[i].first)
            
            leptons = []
            
            # loop over particles in the event
            for i in range(event.NUP):  
                
                # ISTUP: status code 
                # -1 = incoming parton, 1 = final-state parton, 2 = intermediate resonance
                if event.ISTUP[i] != 1: # skip incoming and intermediate particles
                    continue
                
                # PUP: particle momentum
                # 0,1,2,3,4 = px,py,pz,E,m
                particle = ROOT.Math.PxPyPzEVector( # create a 4-vector
                    event.PUP[i][0], 
                    event.PUP[i][1], 
                    event.PUP[i][2], 
                    event.PUP[i][3]
                )

                # IDUP: pdg ID 
                # +/-11 = electron, +/-13 = muon
                # for now we only use electrons and muons
                if abs(event.IDUP[i]) in [11,13]: 
                    leptons.append(particle)

                    l_pdgid.push_back(event.IDUP[i])
                    l_px.push_back(particle.Px())
                    l_py.push_back(particle.Py())
                    l_pz.push_back(particle.Pz())
                    l_E.push_back(particle.E())
                    l_pt.push_back(particle.Pt())
                    l_eta.push_back(particle.Eta())
                    l_phi.push_back(particle.Phi())

            # invariant mass and rapidity of the dilepton system
            mll[0] = (leptons[0]+leptons[1]).M()
            yll[0] = (leptons[0]+leptons[1]).Rapidity()
            ptll[0] = (leptons[0]+leptons[1]).pt()

            # fill the tree and reset variables before moving to the next event 
            ttree.Fill()

            mll[0] = 0.0
            yll[0] = 0.0
            ptll[0] = 0.0
            sm_weight[0] = 0.0
            smeft_weight.clear()
            l_pdgid.clear()
            l_px.clear()
            l_py.clear()
            l_pz.clear()
            l_E.clear()
            l_pt.clear()
            l_eta.clear()
            l_phi.clear()

            # print a status update every 10000 events
            if nevent % 10000 == 0: print('done %s events' % nevent)

    os.system('gzip %s' % infile) # compress the file again

# save the root file
outfile.Write()
outfile.Close()

