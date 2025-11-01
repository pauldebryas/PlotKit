import uproot
import numpy as np
import os
import glob

# Parameters
period = "2016_HIPM"
tag = "FinalProd"

channels = ['tte', 'ttm', 'tee', 'tmm', 'tem']

base_path = f"D_RootHist/results/{period}/{tag}"

print("-" * 76)

for channel in channels:
    folder = os.path.join(base_path, channel, "FakeRate/SignalRegion")
    root_files = glob.glob(os.path.join(folder, "*.root"))

    print(f"channel {channel}")

    for root_file in root_files:
        fname = os.path.basename(root_file)
        found_negative = False
        messages = []

        with uproot.open(root_file) as f:
            # Look for a "signal_region" directory
            if "signal_region" not in f:
                continue
            signal_region = f["signal_region"]

            for hist_name, obj in signal_region.items():
                if not obj.classname.startswith("TH1"):
                    continue

                values, edges = obj.to_numpy()

                for i, val in enumerate(values):
                    if val < 0:
                        low, high = edges[i], edges[i + 1]
                        msg = (f"             Histogram: {hist_name:<50} | "
                               f"Bin {i:<3} | value = {val}")
                        messages.append(msg)

        if messages:
            print(f"       Negative bin found in root file {fname}")
            for m in messages:
                print(m)
            print("")

    print("-" * 76)