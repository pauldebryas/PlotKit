#!/usr/bin/env python
import uproot
import numpy as np
import matplotlib.pyplot as plt
import sys
import os

def plot_systematic_variation(root_path, output_path, branch_name, NPname, DV, channel):
    """
    Reads nominal, Up, and Down histograms from a ROOT file using uproot and plots them.
    Expected paths inside ROOT file:
      signal_region/{branch_name}
      signal_region/{branch_name}_{NPname}Up
      signal_region/{branch_name}_{NPname}Down
    """

    if not os.path.exists(root_path):
        print(f" Error: File not found -> {root_path}")
        sys.exit(1)

    # --- Open ROOT file ---
    with uproot.open(root_path) as file:
        # Try to get histograms
        try:
            h_nom = file[f"signal_region/{branch_name}"]
            h_up = file[f"signal_region/{branch_name}_{NPname}Up"]
            h_down = file[f"signal_region/{branch_name}_{NPname}Down"]
        except KeyError as e:
            print(f" Missing histogram: {e}")
            sys.exit(1)

        # --- Convert to numpy ---
        y_nom, x_nom = h_nom.to_numpy(flow=False)
        y_up, x_nom_up = h_up.to_numpy(flow=False)
        y_down, x_nom_down = h_down.to_numpy(flow=False)

    # --- Compute bin centers ---
    if np.array_equal(x_nom, x_nom_up) and np.array_equal(x_nom, x_nom_down):
        # Sometimes uproot returns bin centers instead of edges
        x_centers = 0.5 * (x_nom[1:] + x_nom[:-1])
    else:
        print(x_nom)
        print(x_nom_up)
        print(x_nom_down)
        raise

    # --- Normalize (optional, comment out if not desired) ---
    # if y_nom.sum() > 0:
    #     y_nom /= y_nom.sum()
    #     y_up /= y_up.sum()
    #     y_down /= y_down.sum()

    #print(y_up)
    #print(y_down)

    # --- Plot ---
    plt.figure(figsize=(8, 6))
    plt.step(x_nom, np.append(y_nom, 0), where="post", color="black", linewidth=0.5, label="Nominal")
    plt.step(x_nom, np.append(y_up, 0) , where="post", color="red", linewidth=0.5, linestyle="--", label=f"{NPname} Up")
    plt.step(x_nom, np.append(y_down, 0), where="post", color="blue", linewidth=0.5, linestyle="--", label=f"{NPname} Down")
    #plt.axhline(y=0.35520994, color="gray", linestyle="--", linewidth=0.5)

    plt.title(f"{branch_name} systematic variation: {NPname}")
    plt.xlabel(f"{DV}")
    plt.ylabel("counts")
    plt.legend()
    plt.grid(True, linestyle=":", alpha=0.7)
    plt.tight_layout()

    # --- Save PDF ---
    out_pdf = output_path + f"_{branch_name}_{NPname}_{channel}.pdf"
    plt.savefig(out_pdf)
    print(f" Saved plot to {out_pdf}")

if __name__ == "__main__":

    channels = ["tee"] #["tee", "tmm", "tte", "ttm", "tem"]
    DV='DNNscore'
    NPname="statFakeFactorsTau2017"

    period = "2017"
    tag= "FinalProd"
    HNLmass= '200'
    branch_name="FakeBackground"

    output_rel_path= "other_scripts/figures/"
    global_path = "/afs/cern.ch/user/p/pdebryas/HNL_analysis/Analysis/PlotKit"
    output_path = f"{global_path}/{output_rel_path}"

    for channel in channels:
        root_path= f"{global_path}/D_RootHist/results/{period}/{tag}/{channel}/FakeRate/SignalRegion/TH1_{DV}_HNLMass{HNLmass}_HNL.root"
        plot_systematic_variation(root_path, output_path, branch_name, NPname, DV, channel)
