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
            'linear': hist.Hist(variables['mll']['axis'], storage=hist.storage.Weight())
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
        A, _ = np.linalg.solve(C, W_evt)  #keep only linear term
        w_reco_linear = w0 + A * c

        for idx, (cmin, cmax) in enumerate(cos_bins):
            if cmin <= cos_val < cmax:
                histos[op][idx]['sm'].fill(mll=mll_val, weight=w0)
                histos[op][idx]['linear'].fill(mll=mll_val, weight=w_reco_linear)
                break

for op in operators:
    fig, ax = plt.subplots(nrows=2, figsize=(8, 8), dpi=300, sharex=True, gridspec_kw={'height_ratios': [3, 1]})
    hep.cms.label("Private Work", data=True, ax=ax[0], rlabel='')

    all_ratios = []
    for idx, (cmin, cmax) in enumerate(cos_bins):
        h_sm = histos[op][idx]['sm']
        h_lin = histos[op][idx]['linear']
        edges = h_sm.axes[0].edges

        ax[0].stairs(h_sm.values(), edges, label=f'SM |cosθ| in [{cmin}, {cmax})', color=colors[idx], linewidth=1.5)
        ax[0].stairs(h_lin.values(), edges, label=f'{op} linear |cosθ| in [{cmin}, {cmax})', linestyle='--', color=colors[idx])

        denominator = np.where(h_sm.values() > 0., h_sm.values(), 1e-6)
        ratio = h_lin.values() / denominator
        all_ratios.append(ratio)
        ax[1].stairs(ratio, edges, edgecolor=colors[idx], fill=False, linewidth=1.0, label=f'Linear / SM [{cmin},{cmax})')

    ax[0].set_yscale('log')
    ax[0].set_ylabel('Events')
    ax[0].legend()

    ax[1].axhline(1.0, color='gray', linestyle='dashed', linewidth=1)
    ax[1].set_ylabel('BSM / SM')
    ax[1].set_xlabel('mll [GeV]')
    ax[1].set_xlim(0, 1000)
    ax[1].legend(fontsize=10)

    rvalid = np.concatenate([r[np.isfinite(r)] for r in all_ratios])
    if len(rvalid) > 0:
        rmin, rmax = np.min(rvalid), np.max(rvalid)
        pad = 0.1
        ax[1].set_ylim(rmin - abs(rmin)*pad, rmax + abs(rmax)*pad)
    else:
        ax[1].set_ylim(0.5, 1.5)

    plt.tight_layout()
    fig.savefig(f"{args.output}/mll_in_costheta_ratio_linear_{op}.png", bbox_inches="tight")
    plt.close(fig)

print(f"All linear-only plots saved in '{args.output}'")
