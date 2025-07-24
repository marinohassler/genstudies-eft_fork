import ROOT
import os
import numpy as np
import hist
import mplhep as hep
from matplotlib import pyplot as plt
from argparse import ArgumentParser
from tqdm import tqdm

style = hep.style.CMS
style['font.size'] = 18
plt.style.use(style)

parser = ArgumentParser()
parser.add_argument('-i', '--input', required=True, help='Input ROOT file')
parser.add_argument('-o', '--output', default='plots', help='Output directory')
args = parser.parse_args()

if not os.path.isdir(args.output):
    os.makedirs(args.output)

def build_lv(pt, eta, phi, mass=0.0):
    px = pt * np.cos(phi)
    py = pt * np.sin(phi)
    pz = pt * np.sinh(eta)
    E  = np.sqrt(px**2 + py**2 + pz**2 + mass**2)
    lv = ROOT.TLorentzVector()
    lv.SetPxPyPzE(px, py, pz, E)
    return lv

variables = {
    'mll': {
        'function': lambda t: t.mll,
        'axis': hist.axis.Regular(50, 0, 1000)
    },
    'yll': {
        'function': lambda t: (
            build_lv(t.Lepton_pt[0], t.Lepton_eta[0], t.Lepton_phi[0]) +
            build_lv(t.Lepton_pt[1], t.Lepton_eta[1], t.Lepton_phi[1])
        ).Rapidity(),
        'axis': hist.axis.Regular(60, -2.5, 2.5)
    },
    'cosTheta': {
        'function': lambda t: t.cosTheta,
        'axis': hist.axis.Regular(50, -1, 1)
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

infile = ROOT.TFile(args.input, 'READ')
tree = infile.Get('events')
n_entries = tree.GetEntries()

c1, c2 = 0.5, 1.0
c = 1.0
C = np.array([[c1, c1**2], [c2, c2**2]])

histos = {}
for op in operators:
    histos[op] = {}
    for k, v in variables.items():
        histos[op][k] = {
            'sm': hist.Hist(v['axis'], storage=hist.storage.Weight()),
            'linear': hist.Hist(v['axis'], storage=hist.storage.Weight()),
            'full': hist.Hist(v['axis'], storage=hist.storage.Weight())
        }

print(f"Processing {n_entries} events …")
for i in tqdm(range(n_entries), desc="Events"):
    tree.GetEntry(i)
    w0 = tree.EventWeight_SM

    for op, (idx1, idx2) in operators.items():
        w1 = tree.EventWeight_SMEFT[idx1]
        w2 = tree.EventWeight_SMEFT[idx2]

        W_evt = np.array([w1-w0, w2-w0])
        A,B = np.linalg.solve(C,W_evt)

        w_reco_lin  = w0 + A*c
        w_reco_full = w0 + A*c + B*c**2

        for k, v in variables.items():
            val = v['function'](tree)
            histos[op][k]['sm'].fill(val, weight=w0)
            histos[op][k]['linear'].fill(val, weight=w_reco_lin)
            histos[op][k]['full'].fill(val, weight=w_reco_full)

#Plotting both
def plot(h, var, op):
    v = h
    fig, ax = plt.subplots(2,1, sharex=True,
        gridspec_kw={"height_ratios": [3,1]}, dpi=300)
    fig.tight_layout(pad=-0.5)

    hep.cms.label('Private Work', data=True, ax=ax[0])

    ax[0].set_yscale('log')
    ax[0].set_ylabel('Events')
    ax[0].set_xlim(v['sm'].axes[0].edges[0], v['sm'].axes[0].edges[-1])

    #SM
    ax[0].stairs(
        v['sm'].values(),
        v['sm'].axes[0].edges,
        label=f"SM [{int(np.sum(v['sm'].values())):,}]",
        color='cornflowerblue',
        edgecolor='black',
        fill=True,
        linewidth=1.0,
    )

    #Linear
    ax[0].stairs(
        v['linear'].values(),
        v['linear'].axes[0].edges,
        label=f"{op} rec. BSM (linear) [{int(np.sum(v['linear'].values())):,}]",
        color='maroon',
        linestyle='--',
        fill=False,
        linewidth=1.2,
    )

    #Full
    ax[0].stairs(
        v['full'].values(),
        v['full'].axes[0].edges,
        label=f"{op} rec. BSM (lin+quad) [{int(np.sum(v['full'].values())):,}]",
        color='green',
        linestyle='-',
        fill=False,
        linewidth=1.0,
    )

    ax[0].legend(loc="upper right", fontsize=10, frameon=False)

    #Ratios
    denominator = np.where(v['sm'].values() > 0., v['sm'].values(), 1e-6)
    ratio_lin  = v['linear'].values() / denominator
    ratio_full = v['full'].values() / denominator

    ax[1].set_ylabel('BSM / SM')
    ax[1].set_xlabel(var)

    ax[1].axhline(1.0, color='grey', linestyle='dashed', linewidth=1)

    ax[1].stairs(
        ratio_lin,
        v['sm'].axes[0].edges,
        edgecolor='maroon',
        fill=False,
        linewidth=1.0,
        label='Linear/SM'
    )

    ax[1].stairs(
        ratio_full,
        v['sm'].axes[0].edges,
        edgecolor='green',
        fill=False,
        linewidth=1.0,
        linestyle='-',
        label='Full/SM'
    )

    ax[1].legend(loc="upper right", fontsize=10)

    #Ratio y–range
    rvalid = np.concatenate([ratio_lin[np.isfinite(ratio_lin)], ratio_full[np.isfinite(ratio_full)]])
    if len(rvalid) > 0:
        rmax = np.max(rvalid)
        rmin = np.min(rvalid)
        pad = 0.1
        ylim_low = rmin - abs(rmin)*pad
        ylim_high = rmax + abs(rmax)*pad
        ax[1].set_ylim(ylim_low, ylim_high)
    else:
        ax[1].set_ylim(0.5, 1.5)

    fig.savefig(f"{args.output}/{op}_{var}_both.png", bbox_inches="tight")
    plt.close(fig)

print(f"Saving plots …")
for op in operators:
    for k in variables:
        plot(histos[op][k], k, op)

print(f"All plots saved in `{args.output}`")
