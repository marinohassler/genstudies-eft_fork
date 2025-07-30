
import ROOT
import os
import hist
import mplhep as hep
import numpy as np
from matplotlib import pyplot as plt
from argparse import ArgumentParser
from tqdm import tqdm

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
        'axis': hist.axis.Regular(50, -1, 1, name='cosTheta')
    }
}

cos_bins = [(0.0, 0.5), (0.5, 1.0)]
colors = ['blue', 'red', 'green', 'red']

operators = {
    "clj1": (1, 2),
    "clj3": (3, 4),
    "ceu": (5, 6),
    "ced": (7, 8),
    "cje": (9, 10),
    "clu": (11, 12),
    "cld": (13, 14),
}

c1, c2 = 0.5, 1.0
c = 1.0
C = np.array([[c1, c1**2], [c2, c2**2]])

os.makedirs(args.output, exist_ok=True)
infile = ROOT.TFile.Open(args.input)
tree = infile.Get('events')
n_entries = tree.GetEntries()

histos = {}
for op in operators:
    histos[op] = []
    for _ in cos_bins:
        histos[op].append({
            'sm': hist.Hist(variables['mll']['axis'], storage=hist.storage.Weight()),
            'full': hist.Hist(variables['mll']['axis'], storage=hist.storage.Weight())
        })

print(f"Processing {n_entries} events ...")
for i in tqdm(range(n_entries), desc="Events"):
    tree.GetEntry(i)
    cos_val = abs(variables['cosTheta']['function'](tree))
    mll_val = variables['mll']['function'](tree)
    w0 = tree.EventWeight_SM

    for op, (idx1, idx2) in operators.items():
        w1 = tree.EventWeight_SMEFT[idx1]
        w2 = tree.EventWeight_SMEFT[idx2]
        W_evt = np.array([w1 - w0, w2 - w0])
        A, B = np.linalg.solve(C, W_evt)

        w_reco_full = w0 + A * c + B * c**2

        for idx, (cmin, cmax) in enumerate(cos_bins):
            if cmin <= cos_val < cmax:
                histos[op][idx]['sm'].fill(mll=mll_val, weight=w0)
                histos[op][idx]['full'].fill(mll=mll_val, weight=w_reco_full)
                break

for op in operators:
    fig, ax = plt.subplots(figsize=(8, 6), dpi=300)
    hep.cms.label("Private Work", data=True, ax=ax, rlabel='')

    for idx, (cmin, cmax) in enumerate(cos_bins):
        h_sm = histos[op][idx]['sm']
        h_full = histos[op][idx]['full']
        edges = h_sm.axes[0].edges

        ax.stairs(h_sm.values(), edges, label=f'SM |cos \u03B8| in [{cmin}, {cmax})', color=colors[idx], linewidth=1.5)
        ax.stairs(h_full.values(), edges, label=f'{op} full |cos \u03B8| in [{cmin}, {cmax})', linestyle='--', color=colors[idx])

    ax.set_yscale('log')
    ax.set_xlabel('mll [GeV]')
    ax.set_ylabel('Events')
    ax.legend()
    plt.tight_layout()
    fig.savefig(f"{args.output}/mll_in_costheta_2bins_{op}.png", bbox_inches="tight")
    plt.close(fig)

print(f"All plots saved in '{args.output}'")
