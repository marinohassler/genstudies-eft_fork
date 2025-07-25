import ROOT
import os
import hist
import mplhep as hep
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from argparse import ArgumentParser

style = hep.style.CMS
style['font.size'] = 18
plt.style.use(style)

parser = ArgumentParser()
parser.add_argument('-i', '--input', required=True, help='path to ROOT file')
parser.add_argument('-o', '--output', default='plots', help='output directory')
args = parser.parse_args()

def build_lv(pt, eta, phi, mass=0.0):
    px = pt * np.cos(phi)
    py = pt * np.sin(phi)
    pz = pt * np.sinh(eta)
    energy = np.sqrt(px**2 + py**2 + pz**2 + mass**2)
    lv = ROOT.TLorentzVector()
    lv.SetPxPyPzE(px, py, pz, energy)
    return lv

mll_axis = hist.axis.Regular(50, 0, 250, name='mll')
yll_axis = hist.axis.Regular(60, -3, 3, name='yll')
costheta_axis = hist.axis.Regular(50, -1, 1, name='costheta')

def get_mll(tree):
    return tree.mll

def get_yll(tree):
    lv1 = build_lv(tree.Lepton_pt[0], tree.Lepton_eta[0], tree.Lepton_phi[0])
    lv2 = build_lv(tree.Lepton_pt[1], tree.Lepton_eta[1], tree.Lepton_phi[1])
    return (lv1 + lv2).Rapidity()

def get_costheta(tree):
        return tree.cosTheta

hist_mll_yll = hist.Hist(mll_axis, yll_axis, storage=hist.storage.Weight())
hist_mll_costheta = hist.Hist(mll_axis, costheta_axis, storage=hist.storage.Weight())

infile = ROOT.TFile.Open(args.input)
tree = infile.Get('events')

for i in range(tree.GetEntries()):
    tree.GetEntry(i)
    mll_val = get_mll(tree)
    yll_val = get_yll(tree)
    costheta_val = get_costheta(tree)
    weight_sm = tree.EventWeight_SM

    hist_mll_yll.fill(mll=mll_val, yll=yll_val, weight=weight_sm)
    hist_mll_costheta.fill(mll=mll_val, costheta=costheta_val, weight=weight_sm)

os.makedirs(args.output, exist_ok=True)

def plot_2d_hist(histogram, xlabel, ylabel, title, output_filename):
    fig, ax = plt.subplots(figsize=(8, 6), dpi=300)
    hep.cms.label("Private Work", data=True, ax=ax, rlabel='')

    im = ax.imshow(
        histogram.values().T,
        origin='lower',
        extent=[
            histogram.axes[0].edges[0], histogram.axes[0].edges[-1],
            histogram.axes[1].edges[0], histogram.axes[1].edges[-1]
        ],
        aspect='auto',
        interpolation='nearest',
        norm=LogNorm(vmin=1e-1, vmax=histogram.values().max())
    )

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    cbar = fig.colorbar(im, ax=ax, label='Events (weighted)')
    plt.tight_layout()
    fig.savefig(f"{args.output}/{output_filename}", bbox_inches='tight')
    plt.close(fig)

#plot_2d_hist(hist_mll_yll, 'mll [GeV]', 'yll', 'DY SM: mll vs yll', 'mll_vs_yll_SM_2D.png')
plot_2d_hist(hist_mll_costheta, 'mll [GeV]', 'cosTheta', 'DY SM: mll vs cosTheta', 'mll_vs_costheta_SM_2D_1.png')


#2D histogram clearly shows the expected Z boson resonance peak in mll
#and the typical rapidity spread of DY events in SM.

