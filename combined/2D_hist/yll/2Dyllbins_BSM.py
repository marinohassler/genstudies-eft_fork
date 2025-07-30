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

def build_lv(pt, eta, phi, mass=0.0):
    px = pt * np.cos(phi)
    py = pt * np.sin(phi)
    pz = pt * np.sinh(eta)
    energy = np.sqrt(px**2 + py**2 + pz**2 + mass**2)
    lv = ROOT.TLorentzVector()
    lv.SetPxPyPzE(px, py, pz, energy)
    return lv

def get_yll(tree):
    lv1 = build_lv(tree.Lepton_pt[0], tree.Lepton_eta[0], tree.Lepton_phi[0])
    lv2 = build_lv(tree.Lepton_pt[1], tree.Lepton_eta[1], tree.Lepton_phi[1])
    return abs((lv1 + lv2).Rapidity())

variables = {
    'mll': {
        'function': lambda t: t.mll,
        'axis': hist.axis.Regular(50, 0, 1000, name='mll')
    }
}

operators = {
    "clj1": (1, 2),
    "clj3": (3, 4),
    "ceu":  (5, 6),
    "ced":  (7, 8),
    "cje":  (9,10),
    "clu": (11,12),
    "cld": (13,14),
}

c1, c2 = 0.5, 1.0
c = 1.0
C = np.array([[c1, c1**2], [c2, c2**2]])

yll_bins = [(0.0, 1.0), (1.0, 2.0), (2.0, 2.5)]
colors = ['blue', 'orange', 'green', 'red', 'purple']

os.makedirs(args.output, exist_ok=True)

infile = ROOT.TFile.Open(args.input)
tree = infile.Get('events')
n_entries = tree.GetEntries()

mll_axis = variables['mll']['axis']
histos = {
    op: {
        'sm': [hist.Hist(mll_axis, storage=hist.storage.Weight()) for _ in yll_bins],
        'linear': [hist.Hist(mll_axis, storage=hist.storage.Weight()) for _ in yll_bins],
        'full': [hist.Hist(mll_axis, storage=hist.storage.Weight()) for _ in yll_bins],
    }
    for op in operators
}


print(f"Processing {n_entries} events …")
for i in tqdm(range(n_entries), desc="Events"):
    tree.GetEntry(i)
    yll_val = get_yll(tree)
    mll_val = variables['mll']['function'](tree)
    w0 = tree.EventWeight_SM

    for idx, (ymin, ymax) in enumerate(yll_bins):
        if ymin <= yll_val < ymax:
            for op, (idx1, idx2) in operators.items():
                w1 = tree.EventWeight_SMEFT[idx1]
                w2 = tree.EventWeight_SMEFT[idx2]
                W_evt = np.array([w1 - w0, w2 - w0])
                A, B = np.linalg.solve(C, W_evt)

                w_reco_lin = w0 + A * c
                w_reco_full = w0 + A * c + B * c**2

                histos[op]['sm'][idx].fill(mll=mll_val, weight=w0)
                histos[op]['linear'][idx].fill(mll=mll_val, weight=w_reco_lin)
                histos[op]['full'][idx].fill(mll=mll_val, weight=w_reco_full)
            break


def plot(h, var, op, bin_label):
    fig, ax = plt.subplots(2,1, sharex=True,
        gridspec_kw={"height_ratios": [3,1]}, dpi=300)
    fig.tight_layout(pad=-0.5)

    hep.cms.label('Private Work', data=True, ax=ax[0])

    ax[0].set_yscale('log')
    ax[0].set_ylabel('Events')
    ax[0].set_xlim(h['sm'].axes[0].edges[0], h['sm'].axes[0].edges[-1])

    ax[0].stairs(h['sm'].values(), h['sm'].axes[0].edges,
        label=f"SM [{int(np.sum(h['sm'].values())):,}]",
        color='cornflowerblue', edgecolor='black', fill=True, linewidth=1.0)

    ax[0].stairs(h['linear'].values(), h['linear'].axes[0].edges,
        label=f"{op} linear [{int(np.sum(h['linear'].values())):,}]",
        color='maroon', linestyle='--', fill=False, linewidth=1.2)

    ax[0].stairs(h['full'].values(), h['full'].axes[0].edges,
        label=f"{op} full [{int(np.sum(h['full'].values())):,}]",
        color='green', linestyle='-', fill=False, linewidth=1.0)

    ax[0].legend(loc="upper right", fontsize=10, frameon=False)

    denominator = np.where(h['sm'].values() > 0., h['sm'].values(), 1e-6)
    ratio_lin  = h['linear'].values() / denominator
    ratio_full = h['full'].values() / denominator

    ax[1].set_ylabel('BSM / SM')
    ax[1].set_xlabel(var)
    ax[1].axhline(1.0, color='grey', linestyle='dashed', linewidth=1)

    ax[1].stairs(ratio_lin, h['sm'].axes[0].edges,
        edgecolor='maroon', fill=False, linewidth=1.0, label='Linear/SM')

    ax[1].stairs(ratio_full, h['sm'].axes[0].edges,
        edgecolor='green', fill=False, linewidth=1.0, linestyle='-', label='Full/SM')

    ax[1].legend(loc="upper right", fontsize=10)

    rvalid = np.concatenate([ratio_lin[np.isfinite(ratio_lin)], ratio_full[np.isfinite(ratio_full)]])
    if len(rvalid) > 0:
        rmax = np.max(rvalid)
        rmin = np.min(rvalid)
        pad = 0.1
        ax[1].set_ylim(rmin - abs(rmin)*pad, rmax + abs(rmax)*pad)
    else:
        ax[1].set_ylim(0.5, 1.5)

    fig.savefig(f"{args.output}/{op}_{bin_label}_{var}_BSM.png", bbox_inches="tight")
    plt.close(fig)


for op in operators:
    for idx, (ymin, ymax) in enumerate(yll_bins):
        bin_label = f"yll_{ymin}_{ymax}"
        plot({
            'sm': histos[op]['sm'][idx],
            'linear': histos[op]['linear'][idx],
            'full': histos[op]['full'][idx]
        }, var='mll', op=op, bin_label=bin_label)
