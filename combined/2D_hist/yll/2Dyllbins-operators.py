
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

yll_bins = [(0.0, 1.0), (1.0, 2.0), (2.0, 2.5)]
colors = ['black', 'blue', 'red', 'green', 'orange']

operators = {
    "clj1": (1, 2),
    "clj3": (3, 4),
    "ceu": (5, 6),
    "ced": (7, 8),
    "cje": (9, 10),
    "clu": (11, 12),
    "cld": (13, 14),
}

os.makedirs(args.output, exist_ok=True)
infile = ROOT.TFile.Open(args.input)
tree = infile.Get('events')
n_entries = tree.GetEntries()

c1, c2 = 0.5, 1.0
c = 1.0
C = np.array([[c1, c1**2], [c2, c2**2]])

histos = {}
for op in operators:
    histos[op] = []
    for _ in yll_bins:
        histos[op].append({
            'sm': hist.Hist(variables['mll']['axis'], storage=hist.storage.Weight()),
            'full': hist.Hist(variables['mll']['axis'], storage=hist.storage.Weight())
        })

print(f"Processing {n_entries} events …")
for i in tqdm(range(n_entries), desc="Events"):
    tree.GetEntry(i)
    yll_val = get_yll(tree)
    mll_val = variables['mll']['function'](tree)
    w0 = tree.EventWeight_SM

    for op, (idx1, idx2) in operators.items():
        w1 = tree.EventWeight_SMEFT[idx1]
        w2 = tree.EventWeight_SMEFT[idx2]

        W_evt = np.array([w1 - w0, w2 - w0])
        A, B = np.linalg.solve(C, W_evt)

        w_reco_full = w0 + A * c + B * c**2

        for idx, (ymin, ymax) in enumerate(yll_bins):
            if ymin <= yll_val < ymax:
                histos[op][idx]['sm'].fill(mll=mll_val, weight=w0)
                histos[op][idx]['full'].fill(mll=mll_val, weight=w_reco_full)
                break


for op in operators:
    fig, ax = plt.subplots(figsize=(8, 6), dpi=300)
    hep.cms.label("Private Work", data=True, ax=ax, rlabel='')
    
    for idx, (ymin, ymax) in enumerate(yll_bins):
        h_sm = histos[op][idx]['sm']
        h_full = histos[op][idx]['full']
        edges = h_sm.axes[0].edges

        ax.stairs(h_sm.values(), edges, label=f'SM |yll| ∈ [{ymin}, {ymax})', color=colors[idx], linewidth=1.5)
        ax.stairs(h_full.values(), edges, label=f'{op} full |yll| ∈ [{ymin}, {ymax})', linestyle='--', color=colors[idx])

    ax.set_yscale('log')
    ax.set_xlabel('mll [GeV]')
    ax.set_ylabel('Events')
    ax.legend()
    plt.tight_layout()
    fig.savefig(f"{args.output}/mll_in_yll_3bins_{op}.png", bbox_inches="tight")
    plt.close(fig)

print(f"Plots saved in {args.output}")
