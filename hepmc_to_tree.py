import ROOT
import pyhepmc
import fastjet
import os
from array import array
from argparse import ArgumentParser

jet_definition = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.4)

def make_branches(tree, variables):
    for var in variables:
        if isinstance(variables[var], array):
            tree.Branch(var, variables[var], '%s/D' % var)
        else:
            tree.Branch(var, variables[var])

def clear_variables(variables):
    for var in variables:
        if isinstance(variables[var], array):
            variables[var][0] = 0.0
        else:
            variables[var].clear()

parser = ArgumentParser()
parser.add_argument('-d', '--directory', default=None, help='directory of input HepMC files and for output ROOT file')
parser.add_argument('-f', '--filename', default='events.root', help='name of output ROOT file')
parser.add_argument('-n', '--nfiles', type=int, default=None, help='only run on the first n input files')
args = parser.parse_args()

outfile = ROOT.TFile(os.path.join(args.directory, args.filename), 'RECREATE')

infiles = [f for f in os.listdir(args.directory) if f.startswith('events_') and f.endswith('.hepmc.gz')]
infiles = sorted(infiles, key=lambda x: float(x[7:-9])) # sort them
if args.nfiles is not None: 
    infiles = infiles[:args.nfiles]

variables = {
    'mll': array('d',[0.0]),
    'yll': array('d',[0.0]),
    'ptll': array('d',[0.0]),
    'EventWeight_SM': array('d',[0.0]),
    'EventWeight_SMEFT': ROOT.std.vector('double')(),
    'Lepton_pdgId': ROOT.std.vector('double')(),
    'Lepton_pt': ROOT.std.vector('double')(),
    'Lepton_eta': ROOT.std.vector('double')(),
    'Lepton_phi': ROOT.std.vector('double')(),
    'Jet_pt': ROOT.std.vector('double')(),
    'Jet_eta': ROOT.std.vector('double')(),
    'Jet_phi': ROOT.std.vector('double')(),
}

ttree = ROOT.TTree('events','tree of generated events')
make_branches(ttree, variables)

nevents = 0
for infile in infiles:
    infile = os.path.join(args.directory, infile)
    os.system('gunzip %s' % infile) # decompress the file
    infile = infile.replace('.gz','')
    with pyhepmc.open(infile) as f:
        for event in f:
            nevents += 1
            variables['EventWeight_SM'][0] = event.weights[0]
            for j in range(1,len(event.weights)):
                variables['EventWeight_SMEFT'].push_back(event.weights[j])
            particles = list()
            leptons = list()
            for particle in event.particles:
                if particle.status == 1: # final state
                    momentum = particle.momentum
                    particles.append(fastjet.PseudoJet(momentum.px, momentum.py, momentum.pz, momentum.e))
                    if abs(particle.pid) in [11,13]:
                        lepton = ROOT.Math.PxPyPzEVector(
                            momentum.px, momentum.py, momentum.pz, momentum.e
                        )
                        lepton.pid = particle.pid
                        leptons.append(lepton)
            leptons = sorted(leptons, key=lambda x: x.pt(), reverse=True)
            leptons = list(filter(lambda x: x.pt()>10., leptons))
            for lepton in leptons:
                variables['Lepton_pdgId'].push_back(lepton.pid)
                variables['Lepton_pt'].push_back(lepton.pt())
                variables['Lepton_eta'].push_back(lepton.eta())
                variables['Lepton_phi'].push_back(lepton.phi())

            cluster = fastjet.ClusterSequence(particles, jet_definition)
            jets = cluster.inclusive_jets(ptmin=10.)
            jets = sorted(jets, key=lambda x: x.pt(), reverse=True)
            for jet in jets:
                variables['Jet_pt'].push_back(jet.pt())
                variables['Jet_eta'].push_back(jet.eta())
                variables['Jet_phi'].push_back(jet.phi())

            if len(leptons)>1:
                variables['mll'][0] = (leptons[0]+leptons[1]).M()
                variables['yll'][0] = (leptons[0]+leptons[1]).Rapidity()
                variables['ptll'][0] = (leptons[0]+leptons[1]).pt()
                ttree.Fill()

            clear_variables(variables)

            if nevents % 10000 == 0: print('done %s events' % nevents)
    os.system('gzip %s' % infile)

outfile.Write()
outfile.Close()


