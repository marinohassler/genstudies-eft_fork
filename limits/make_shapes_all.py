import os, math, argparse
import numpy as np
import uproot, awkward as ak
import ROOT
from array import array

ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(False)

# --------------------
#helpers
# --------------------
def to_np(x):
    """awkward -> numpy, flattening jagged if needed."""
    if isinstance(x, ak.highlevel.Array):
        return ak.to_numpy(x)
    return np.asarray(x)

def clean_mask(*arrs):
    """finite mask common to all arrays"""
    mask = np.ones_like(to_np(arrs[0]), dtype=bool)
    for a in arrs:
        v = to_np(a)
        mask &= np.isfinite(v)
    return mask

def whist_with_err(x, w, bins, include_overflow=True):
    """
    Weighted histogram with per-bin errors (sum w^2).
    Returns (sum w, sqrt(sum w^2)).
    Properly folds overflow into the last bin if requested.
    """
    x = to_np(x); w = to_np(w)
    m = clean_mask(x, w)
    x, w = x[m], w[m]

    sw, _  = np.histogram(x, bins=bins, weights=w)
    sw2, _ = np.histogram(x, bins=bins, weights=w*w)

    if include_overflow:
        of_mask = x >= bins[-1]
        if np.any(of_mask):
            sw[-1]  += np.sum(w[of_mask])
            sw2[-1] += np.sum((w[of_mask])**2)

    err = np.sqrt(sw2, dtype=float)
    return sw.astype(float), err.astype(float)

def make_th1(name, contents, errors, edges):
    """TH1D with custom binning + errors."""
    h = ROOT.TH1D(name, name, len(edges)-1, array('d', edges.tolist()))
    h.Sumw2()
    for i, (c,e) in enumerate(zip(contents, errors), start=1):
        h.SetBinContent(i, float(c))
        h.SetBinError(i,   float(e))
    return h

def unit_area_positive(arr, eps=1e-12):
    """Turn (possibly signed) array into positive unit-area shape."""
    a = np.abs(arr).astype(float)
    a[a < eps] = eps
    s = float(np.sum(a))
    return a/s if s > 0 else a

# --------------------
#core per-operator routine
# --------------------
def build_one_operator(fname, treename, mll_branch, w_sm_branch, w_vec_branch,
                       idx_lin, edges, outdir, channel,
                       scale_if_huge=False, scale_threshold=1e6, scale_target=1e4,
                       mll_scale=1.0, debug=False):
    #read tree
    with uproot.open(fname)[treename] as t:
        mll   = t[mll_branch      ].array(library="ak")
        w_sm  = t[w_sm_branch     ].array(library="ak")
        w_vec = t[w_vec_branch    ].array(library="ak")

    #optional unit conversion 
    if mll_scale != 1.0:
        mll = mll * mll_scale

    #debug coverage prints
    if debug:
        try:
            m_min = float(ak.min(mll))
            m_max = float(ak.max(mll))
            n_gt110 = int(ak.sum(mll > 110))
            n_gt200 = int(ak.sum(mll > 200))
            n_gt500 = int(ak.sum(mll > 500))
            print(f"[DEBUG] mll range: {m_min:.3g} – {m_max:.3g}")
            print(f"[DEBUG] counts mll>110,>200,>500: {n_gt110}, {n_gt200}, {n_gt500}")
        except Exception as e:
            print("[DEBUG] mll debug failed:", e)

    #relative -> absolute per-event EFT weights
    w_lin_rel  = w_vec[:, idx_lin]
    w_quad_rel = w_lin_rel**2       

    w_lin_abs  = w_sm * w_lin_rel   
    w_quad_abs = w_sm * w_quad_rel  

    #histograms
    h_sm_np,  h_sm_err  = whist_with_err(mll, w_sm,       edges)
    h_lin_np, _         = whist_with_err(mll, w_lin_abs,  edges)
    h_quad_np, _        = whist_with_err(mll, w_quad_abs, edges)

    I_SM   = float(np.sum(h_sm_np))
    I_LIN  = float(np.sum(h_lin_np))   
    I_QUAD = float(np.sum(h_quad_np))  

    #optional global downscale (keep Combine numerically comfy)
    scale = 1.0
    if scale_if_huge and I_SM > scale_threshold:
        scale = float(scale_target) / I_SM

    if debug:
        print(f"[DEBUG] scale_if_huge={scale_if_huge} threshold={scale_threshold} target={scale_target} -> scale={scale:g}")

    h_sm_np  *= scale
    h_sm_err *= scale
    h_lin_np *= scale
    h_quad_np*= scale
    I_SM   *= scale
    I_LIN  *= scale
    I_QUAD *= scale

    #CED: unit-area positive EFT shapes
    lin_ced  = unit_area_positive(h_lin_np)
    quad_ced = unit_area_positive(h_quad_np)

    #TH1 with errors for SM (others no errors)
    h_sm       = make_th1("sm",       h_sm_np,   h_sm_err, edges)
    h_lin_ced  = make_th1("lin_ced",  lin_ced,   np.zeros_like(lin_ced),  edges)
    h_quad_ced = make_th1("quad_ced", quad_ced,  np.zeros_like(quad_ced), edges)
    h_data     = h_sm.Clone("data_obs")

    #write ROOT
    os.makedirs(outdir, exist_ok=True)
    outroot = os.path.join(outdir, "shapes_tmp.root")
    f = ROOT.TFile(outroot, "RECREATE")
    d = f.mkdir(channel); d.cd()
    for h in (h_sm, h_lin_ced, h_quad_ced, h_data):
        h.Write()
    f.Close()

    return dict(
        outroot=outroot,
        I_SM=I_SM, I_LIN=I_LIN, I_QUAD=I_QUAD, scale=scale
    )

# --------------------
#datacard writer (single-operator model)
# --------------------
def write_datacard(op, channel, outdir, numbers_root, ABS_ILIN, I_QUAD, automc=True):
    card = f"""imax 1
jmax 2
kmax *

shapes * {channel} shapes_{op}.root {channel}/$PROCESS {channel}/$PROCESS_$SYSTEMATIC

bin {channel}
observation -1

# SM absolute; EFT are unit-area shapes (CED)
bin          {channel}   {channel}    {channel}
process      sm          lin_ced      quad_ced
process      1           0            2
rate         -1          1            1

# Use Combine's built-in POI 'r' (scan r over +/-)
# Numbers already include any global scale:
r_lin   rateParam  {channel}  lin_ced   (@0*{ABS_ILIN:.16g})      r
r_quad  rateParam  {channel}  quad_ced  (@0*@0*{I_QUAD:.16g})     r

{"* autoMCStats 0" if not automc else ""}
"""
    #move shapes_tmp.root -> shapes_op.root so the card points to it
    os.replace(numbers_root, os.path.join(outdir, f"shapes_{op}.root"))
    with open(os.path.join(outdir, f"datacard_mll_{op}_scaled.txt"), "w") as dc:
        dc.write(card)

# --------------------
#CLI
# --------------------
if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--file",     required=True, help="path to events.root")
    ap.add_argument("--tree",     default="events")
    ap.add_argument("--mll",      default="mll")
    ap.add_argument("--wsm",      default="EventWeight_SM")
    ap.add_argument("--wvec",     default="EventWeight_SMEFT")
    ap.add_argument("--edges",    default="50,70,90,110,140,200,300,500,800,1200")
    ap.add_argument("--mll-scale", type=float, default=1.0,
                    help="Multiply mll by this (e.g. 0.001 if the branch is in MeV).")
    ap.add_argument("--op",       default="", help="single operator name:index (e.g. ced:8)")
    ap.add_argument("--ops",      default="", help="comma list of name:index (e.g. clj1:2,clj3:4,...)")
    ap.add_argument("--outdir",   default="../combine")
    ap.add_argument("--channel",  default="ch1")
    ap.add_argument("--auto-mcstats", action="store_true",
                    help="leave autoMCStats ON (default is ON)")
    #scaling controls
    ap.add_argument("--scale-if-huge", action="store_true",
                    help="If set, downscale SM to target if I_SM > threshold.")
    ap.add_argument("--scale-threshold", type=float, default=1e6)
    ap.add_argument("--scale-target",    type=float, default=1e4)
    ap.add_argument("--debug", action="store_true")
    args = ap.parse_args()

    edges = np.array([float(x) for x in args.edges.split(",")], dtype=float)

    #build list of (name, index)
    pairs = []
    if args.op:
        name, idx = args.op.split(":")
        pairs.append((name, int(idx)))
    if args.ops:
        for tok in args.ops.split(","):
            name, idx = tok.split(":")
            pairs.append((name, int(idx)))

    if not pairs:
        raise SystemExit("Provide --op NAME:IDX or --ops NAME:IDX,NAME:IDX,...")

    print(f"Reading: {args.file}  tree: {args.tree}\n")

    for name, idx in pairs:
        print(f"=== Operator {name} (index {idx}) ===")
        res = build_one_operator(
            fname=args.file, treename=args.tree,
            mll_branch=args.mll, w_sm_branch=args.wsm, w_vec_branch=args.wvec,
            idx_lin=idx, edges=edges, outdir=args.outdir, channel=args.channel,
            scale_if_huge=args.scale_if_huge,
            scale_threshold=args.scale_threshold, scale_target=args.scale_target,
            mll_scale=args.mll_scale, debug=args.debug
        )

        #numbers (already include any scale)
        ABS_ILIN = abs(res["I_LIN"])
        I_QUAD   = res["I_QUAD"]

        #write shapes_OP.root and card
        write_datacard(
            op=name, channel=args.channel, outdir=args.outdir,
            numbers_root=res["outroot"], ABS_ILIN=ABS_ILIN, I_QUAD=I_QUAD,
            automc=args.auto_mcstats
        )

        print(f" -> Wrote {os.path.join(args.outdir, f'shapes_{name}.root')}")
        expo = (math.log10(res["scale"]) if res["scale"]>0 else 0.0)
        print(f"    Integrals (after global scale 10^{expo:.0f}):")
        print(f"      I_SM   = {res['I_SM']}")
        print(f"      I_LIN  = {res['I_LIN']}   (signed total)")
        print(f"      I_QUAD = {res['I_QUAD']}")
        print(f" -> Wrote {os.path.join(args.outdir, f'datacard_mll_{name}_scaled.txt')}\n")
