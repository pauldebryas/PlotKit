import awkward as ak
import yaml
import os
import uproot
import numpy as np

# ------------------------- Parameters -------------------------
period= '2017'
tag = 'FinalProd'
channel = 'tmm' # could be "tee" or "tmm"
test = True
files_to_split = ['EWKWMinus2Jets_WToLNu_M-50_anatuple.root'] # files to split in test mode Ex: SingleElectron_2016_anatuple, WminusHToTauTau_anatuple, WJetsToLNu_HT-1200To2500_anatuple
# --------------------------------------------------------------
with open(os.path.join(os.getenv("RUN_PATH"), f'common/config/all/config_FakeRate.yaml'), 'r') as f:
    GlobalParms = yaml.safe_load(f)
lepton = GlobalParms['lepton']

def mask_OS_SS(tree, channel):
    l2 = lepton[channel]['lepton2']['name']
    l3 = lepton[channel]['lepton3']['name']
    q2 = tree[f"{l2}_charge"].array(library="ak")
    q3 = tree[f"{l3}_charge"].array(library="ak")
    mask_os = (q2 != q3)
    mask_ss = (q2 == q3)          # <- SS must be equality
    return mask_os, mask_ss

def split_root_files_SS_OS(file, output_file_OS, output_file_SS):
    with uproot.open(file) as src:

        # List trees present
        all_keys = list(src.keys())                  # e.g. ["Events;1","Events_JES_2016_up;1", ... , "bjets;1"]
        event_keys = [k for k in all_keys if k.split(';')[0].startswith("Events")]
        bjets_keys = [k for k in all_keys if k.split(';')[0].startswith("bjets")]

        # Build masks for every Events* tree
        masks = {}
        for k in event_keys:
            masks[k] = mask_OS_SS(src[k], channel)   # tuple (mask_os, mask_ss)

        # Nominal Events mask (for bjets filtering)
        nominal_events_key = next((k for k in event_keys if k.split(';')[0] == "Events"), None)
        nominal_masks = masks.get(nominal_events_key, (None, None))

        # Create the two output ROOT files and write filtered trees
        with uproot.recreate(output_file_OS) as out_os, uproot.recreate(output_file_SS) as out_ss:

            # Write all Events* trees with their own masks
            for k in event_keys:
                tree = src[k]
                mask_os, mask_ss = masks[k]

                arr = tree.arrays(library="ak")      # awkward Record (all branches)
                clean_name = k.split(';')[0]
                out_os[clean_name] = arr[mask_os]
                out_ss[clean_name] = arr[mask_ss]

            # Write bjets using the nominal Events mask (if available), else copy as-is
            for k in bjets_keys:
                tree = src[k]
                arr = tree.arrays(library="ak")  # Awkward Record with fields: n, pt, eta, phi, hadronFlavour, btagDeepFlavB, mass, vetomap

                if nominal_masks[0] is not None:
                    mask_os, mask_ss = nominal_masks

                    # Build dicts of per-branch arrays (this is what uproot expects)
                    os_dict = {name: arr[name][mask_os] for name in arr.fields}
                    ss_dict = {name: arr[name][mask_ss] for name in arr.fields}
                
                    clean_name = k.split(';')[0]
                    if ak.any(mask_os):
                        out_os[clean_name] = os_dict
                    if ak.any(mask_ss):
                        out_ss[clean_name] = ss_dict
                else:
                    # No nominal Events tree found: just copy branches as-is
                    print('Error: No nominal Events tree found')
                    raise

def main():
    eos_path = f'/eos/user/p/pdebryas/HNL/anatuple/{period}/{tag}/'
    folder_tll = os.path.join(eos_path, channel, 'anatuple')
    output_folder_OS = os.path.join(eos_path, f"{channel}_os", 'anatuple')
    output_folder_SS = os.path.join(eos_path, f"{channel}_ss", 'anatuple')
    os.makedirs(output_folder_OS, exist_ok=True)
    os.makedirs(output_folder_SS, exist_ok=True)
    
    files_list = sorted(list(os.listdir(folder_tll)))
    # files_list = ['DYJetsToLL_0J_anatuple.root', 'DYJetsToLL_1J_anatuple.root', 'DYJetsToLL_2J_anatuple.root', 'DYJetsToLL_LHEFilterPtZ-0To50_anatuple.root', 'DYJetsToLL_LHEFilterPtZ-100To250_anatuple.root', 'DYJetsToLL_LHEFilterPtZ-250To400_anatuple.root', 'DYJetsToLL_LHEFilterPtZ-400To650_anatuple.root', 'DYJetsToLL_LHEFilterPtZ-50To100_anatuple.root', 'DYJetsToLL_LHEFilterPtZ-650ToInf_anatuple.root', 'DYJetsToLL_M-50_anatuple.root', 'EWKWMinus2Jets_WToLNu_M-50_anatuple.root', 'EWKWPlus2Jets_WToLNu_M-50_anatuple.root', 'EWKZ2Jets_ZToLL_M-50_anatuple.root', 'GluGluHToTauTau_anatuple.root', 'HNL_VBF_tau_M-1000_anatuple.root', 'HNL_VBF_tau_M-600_anatuple.root', 'HNL_VBF_tau_M-700_anatuple.root', 'HNL_VBF_tau_M-800_anatuple.root', 'HNL_VBF_tau_M-900_anatuple.root', 'HNL_tau_M-1000_anatuple.root', 'HNL_tau_M-100_anatuple.root', 'HNL_tau_M-125_anatuple.root', 'HNL_tau_M-150_anatuple.root', 'HNL_tau_M-200_anatuple.root', 'HNL_tau_M-20_anatuple.root', 'HNL_tau_M-250_anatuple.root', 'HNL_tau_M-300_anatuple.root', 'HNL_tau_M-30_anatuple.root', 'HNL_tau_M-350_anatuple.root', 'HNL_tau_M-400_anatuple.root', 'HNL_tau_M-40_anatuple.root', 'HNL_tau_M-450_anatuple.root', 'HNL_tau_M-500_anatuple.root', 'HNL_tau_M-50_anatuple.root', 'HNL_tau_M-600_anatuple.root', 'HNL_tau_M-60_anatuple.root', 'HNL_tau_M-700_anatuple.root', 'HNL_tau_M-70_anatuple.root', 'HNL_tau_M-75_anatuple.root', 'HNL_tau_M-800_anatuple.root', 'HNL_tau_M-85_anatuple.root', 'HNL_tau_M-900_anatuple.root', 'ST_t-channel_antitop_4f_InclusiveDecays_anatuple.root', 'ST_t-channel_top_4f_InclusiveDecays_anatuple.root', 'ST_tW_antitop_5f_inclusiveDecays_anatuple.root', 'ST_tW_top_5f_inclusiveDecays_anatuple.root', 'SingleElectron_2016_anatuple.root', 'TTTo2L2Nu_anatuple.root', 'TTToHadronic_anatuple.root', 'TTToSemiLeptonic_anatuple.root', 'TTWJetsToLNu_anatuple.root', 'TTWW_anatuple.root', 'TTWZ_anatuple.root', 'TTZToLLNuNu_M-10_anatuple.root', 'TTZZ_anatuple.root', 'W1JetsToLNu_anatuple.root', 'W2JetsToLNu_anatuple.root', 'W3JetsToLNu_anatuple.root', 'W4JetsToLNu_anatuple.root', 'WJetsToLNu_HT-100To200_anatuple.root', 'WJetsToLNu_HT-1200To2500_anatuple.root', 'WJetsToLNu_HT-200To400_anatuple.root', 'WJetsToLNu_HT-2500ToInf_anatuple.root', 'WJetsToLNu_HT-400To600_anatuple.root', 'WJetsToLNu_HT-600To800_anatuple.root', 'WJetsToLNu_HT-70To100_anatuple.root', 'WJetsToLNu_HT-800To1200_anatuple.root', 'WJetsToLNu_anatuple.root', 'WWW_anatuple.root', 'WWZ_anatuple.root', 'WW_anatuple.root', 'WZZ_anatuple.root', 'WZ_anatuple.root', 'WminusHToTauTau_anatuple.root', 'WplusHToTauTau_anatuple.root', 'ZHToTauTau_anatuple.root', 'ZZZ_anatuple.root', 'ZZ_anatuple.root', 'ttHToTauTau_anatuple.root']
    # files_list = ['WJetsToLNu_HT-200To400_anatuple.root', 'WJetsToLNu_HT-2500ToInf_anatuple.root', 'WJetsToLNu_HT-400To600_anatuple.root', 'WJetsToLNu_HT-600To800_anatuple.root', 'WJetsToLNu_HT-70To100_anatuple.root', 'WJetsToLNu_HT-800To1200_anatuple.root', 'WJetsToLNu_anatuple.root', 'WWW_anatuple.root', 'WWZ_anatuple.root', 'WW_anatuple.root', 'WZZ_anatuple.root', 'WZ_anatuple.root', 'WminusHToTauTau_anatuple.root', 'WplusHToTauTau_anatuple.root', 'ZHToTauTau_anatuple.root', 'ZZZ_anatuple.root', 'ZZ_anatuple.root', 'ttHToTauTau_anatuple.root']
    if test == True:
        files_list = files_to_split
    
    for file_name in files_list:
        file_path = os.path.join(folder_tll, file_name)
        output_file_OS = os.path.join(output_folder_OS, file_name)
        output_file_SS = os.path.join(output_folder_SS, file_name)

        print(f"Splitting {file_name}...")
        split_root_files_SS_OS(file_path, output_file_OS, output_file_SS)
        
if __name__ == "__main__":
    main()
