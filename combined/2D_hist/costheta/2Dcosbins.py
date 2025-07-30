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
        'axis': hist.axis.Regular(50, 0, 1000, name='mll')
    },
    'cosTheta': {
        'function': lambda t: t.cosTheta,
        'axis': hist.axis.Regular(50, -1, 1)
    }
}

cos_bins = [(0.0, 0.25), (0.25, 0.5), (0.5, 0.75), (0.75, 1.0)]
colors = ['blue', 'orange', 'green', 'red']

os.makedirs(args.output, exist_ok=True)
infile = ROOT.TFile.Open(args.input)
tree = infile.Get('events')

histos = []
mll_axis = variables['mll']['axis']
for _ in cos_bins:
    histos.append(hist.Hist(mll_axis, storage=hist.storage.Weight()))

for i in range(tree.GetEntries()):
    tree.GetEntry(i)
    cos_val = abs(variables['cosTheta']['function'](tree))
    mll_val = variables['mll']['function'](tree)
    weight = tree.EventWeight_SM
    
    for idx, (ymin, ymax) in enumerate(cos_bins):
        if ymin <= cos_val < ymax:
            histos[idx].fill(mll=mll_val, weight=weight)
            break  

#Plotting
fig, ax = plt.subplots(figsize=(8,6), dpi=300)
hep.cms.label("Private Work", data=True, ax=ax, rlabel='')

for idx, (ymin, ymax) in enumerate(cos_bins):
    vals = histos[idx].values()
    total = np.sum(vals)
    ax.stairs(
        vals,
        histos[idx].axes[0].edges,
        label=f'|cos \u03B8| in [{ymin}, {ymax})',
        color=colors[idx],
        linewidth=1.5,
        fill=False,
    )

ax.set_yscale('log')
ax.set_xlabel('mll [GeV]')
ax.set_ylabel('Events')
ax.legend()
plt.tight_layout()
fig.savefig(f"{args.output}/mll_in_cos_4bins_SM.png", bbox_inches="tight")
plt.close(fig)
