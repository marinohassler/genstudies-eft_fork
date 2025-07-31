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
parser.add_argument('-i', '--input', required=True, help='path to ROOT file')
parser.add_argument('-o', '--output', default='plots', help='output directory')
args = parser.parse_args()

variables = {
    'mll': {
        'function': lambda t: t.mll,
    },
    'cosTheta': {
        'function': lambda t: t.cosTheta,
        'axis': hist.axis.Regular(50, -1, 1, name='cosTheta')
    }
}

mll_bins = [(0, 200), (200, 400), (400, 600), (600, 800), (800, 1000)]
colors = ['blue', 'orange', 'green', 'red', 'purple']

os.makedirs(args.output, exist_ok=True)
infile = ROOT.TFile.Open(args.input)
tree = infile.Get('events')

cos_axis = variables['cosTheta']['axis']
histos = [hist.Hist(cos_axis, storage=hist.storage.Weight()) for _ in mll_bins]

for i in range(tree.GetEntries()):
    tree.GetEntry(i)
    cos_val = variables['cosTheta']['function'](tree)
    mll_val = variables['mll']['function'](tree)
    weight = tree.EventWeight_SM

    for idx, (xmin, xmax) in enumerate(mll_bins):
        if xmin <= mll_val < xmax:
            histos[idx].fill(cosTheta=cos_val, weight=weight)
            break


fig, ax = plt.subplots(figsize=(8,6), dpi=300)
hep.cms.label("Private Work", data=True, ax=ax, rlabel='')

for idx, (xmin, xmax) in enumerate(mll_bins):
    vals = histos[idx].values()
    ax.stairs(
        vals,
        histos[idx].axes[0].edges,
        label=f'mll in [{xmin}, {xmax}) GeV',
        color=colors[idx],
        linewidth=1.5,
        fill=False,
    )

ax.set_xlabel('cosTheta')
ax.set_ylabel('Events')
ax.set_yscale('log')
ax.legend()
plt.tight_layout()
fig.savefig(f"{args.output}/cosTheta_in_mll_5bins_SM.png", bbox_inches="tight")
plt.close(fig)
