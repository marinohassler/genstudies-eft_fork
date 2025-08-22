import ROOT
import os
import hist
import mplhep as hep
import numpy as np
from matplotlib import pyplot as plt
from argparse import ArgumentParser

style = hep.style.CMS
style['font.size'] = 18
plt.style.use(style)

parser = ArgumentParser()
parser.add_argument('-i', '--input', default=None, help='path to ROOT file')
parser.add_argument('-o', '--output', default='plots', help='output directory')
args = parser.parse_args()

# define variables and binning
variables = {
    'mll': {
        'function': lambda tree: tree.mll,
        'axis': hist.axis.Regular(25, 50, 300)
    },
    'ptl': {
        'function': lambda tree: tree.Lepton_pt,
        'axis': hist.axis.Regular(30, 0, 150)
    },
    'detall': {
        'function': lambda tree: abs(tree.Lepton_eta[0]-tree.Lepton_eta[1]),
        'axis': hist.axis.Regular(25, 0, 5)
    },
    'ptll': {
        'function': lambda tree: tree.ptll,
        'axis': hist.axis.Regular(30, 0, 150)
    }
}

# create histograms
histos = {}
for k,v in variables.items():
    histos[k] = { 'sm': hist.Hist(v['axis'], storage=hist.storage.Weight()) }
    histos[k]['clj3'] = hist.Hist(v['axis'], storage=hist.storage.Weight())

# loop through the ROOT tree and fill the histograms
infile = ROOT.TFile(args.input, 'READ')
tree = infile.Get('events')

for i in range(tree.GetEntries()):
    tree.GetEntry(i)
    for k,v in variables.items():
        histos[k]['sm'].fill(v['function'](tree), weight=tree.EventWeight_SM)
        histos[k]['clj3'].fill(v['function'](tree), weight=tree.EventWeight_SMEFT[4]) # rw0004
        # TODO: generalize this to look at each Wilson coefficient
        # for now just use weight 4 (clj3=1): compare with your reweight_card.dat, rw0004

# plotting
if not os.path.isdir(args.output):
    os.makedirs(args.output)

for k,v in histos.items():
    fig,ax = plt.subplots(2, 1, sharex=True, gridspec_kw={"height_ratios": [3,1]}, dpi=300)
    fig.tight_layout(pad=-0.5)
    hep.cms.label('Private Work', data=True, ax=ax[0], rlabel='')
    
    # upper panel: histograms
    ax[0].set_yscale('log')
    ax[0].set_ylabel('Events')
    ax[0].set_xlim(v['sm'].axes[0].edges[0], v['sm'].axes[0].edges[-1])

    ax[0].stairs(
        v['sm'].values(), 
        v['sm'].axes[0].edges, 
        label=f"DY SM [{int(round(np.sum(v['sm'].values()), 0))}]", 
        color='cornflowerblue',
        edgecolor='black',
        fill=True,
        linewidth=1.0,
    )

    ax[0].stairs(
        v['clj3'].values(), 
        v['clj3'].axes[0].edges, 
        label=f"DY clj3=1.0 [{int(round(np.sum(v['clj3'].values()), 0))}]", 
        color='maroon',
        fill=False,
        linewidth=1.0,
    )
    
    ax[0].legend(loc="upper right", fontsize=12)

    # lower panel: ratio
    ax[1].set_ylabel('BSM / SM')
    ax[1].set_xlabel(k)

    ax[1].plot(
        v['sm'].axes[0].edges, 
        np.ones_like(v['sm'].axes[0].edges), 
        color='cornflowerblue', 
        linestyle='dashed',
    )

    # we can not divide by 0
    denominator = np.where(v['sm'].values()>0., v['sm'].values(), 1e-6)
    ax[1].stairs(
        v['clj3'].values()/denominator, 
        v['clj3'].axes[0].edges, 
        edgecolor='maroon',
        fill=False,
        linewidth=1.0,
    )

    ymax = np.max(v['clj3'].values()/denominator)
    ymin = np.min(np.where(v['clj3'].values()/denominator>0., v['clj3'].values()/denominator, 1.))
    ylim = max(abs(1-ymax), abs(1-ymin))
    ax[1].set_ylim(1-ylim, 1+ylim)
    
    fig.savefig(f"{args.output}/{k}.png", bbox_inches="tight")


