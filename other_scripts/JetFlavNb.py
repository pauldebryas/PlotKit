import os
import yaml
import uproot
import numpy as np
import correctionlib

period= '2018'
tag = 'WithCorrections'
channel = 'tmm'

def get_correction_central(corr, period):
    #load correction from central repo
    POG_path = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/"
    Area_dir = {
        '2016_HIPM': '2016preVFP_UL',
        '2016': '2016postVFP_UL',
        '2017': '2017_UL',
        '2018': '2018_UL'
    }
    f_path = os.path.join(POG_path, 'BTV/' + Area_dir[period] + '/btagging.json.gz')
    ceval = correctionlib.CorrectionSet.from_file(f_path)
    return ceval

# btag wp
cset = get_correction_central('btag', period)
#loose WP value
LooseWP = cset["deepJet_wp_values"].evaluate("L")

#input folder where anatuple are stored
input_folder = os.path.join('/eos/user/p/pdebryas/HNL/anatuple', period, tag, channel, 'anatuple')
#list of MC samples
with open(f'{os.getenv("RUN_PATH")}/config/shared/inputs_MCbackground_2018.yaml', 'r') as f:
    inputs = yaml.safe_load(f)

#check no files is missing
list_MCfiles = {}
list_MCmissingfiles = []
for input in inputs:
    files = input.get('files')
    if files != None:
        for file in files:
            if not os.path.isfile(os.path.join(input_folder, file)):
                print(f'missing file: {file}')
                list_MCmissingfiles.append(file)
            else:
                list_MCfiles[file[:-14]] = os.path.join(input_folder, file)

# Open the ROOT file
file = uproot.open(list_MCfiles['WZ'])

# Get the TTree from the file
tree = file["bjets"] 

# Check if the tree is retrieved successfully
if not tree:
    print(f"Error: Could not retrieve TTree bjets.")
else:
    # Get the branch from the tree
    btagDeepFlavB = np.concatenate([item.flatten() for item in np.array(tree["btagDeepFlavB"]) if len(item) > 0]) 
    hadronFlavour = np.concatenate([item.flatten() for item in np.array(tree["hadronFlavour"]) if len(item) > 0]) 
    jetpt = np.concatenate([item.flatten() for item in np.array(tree["pt"]) if len(item) > 0])
    jeteta =  np.concatenate([item.flatten() for item in np.array(tree["eta"]) if len(item) > 0])
    
    if len(btagDeepFlavB) != len(hadronFlavour):
        print('not same length')

    if len(jetpt) != len(jeteta):
        print('not same length')

    if len(jetpt) != len(btagDeepFlavB):
        print('not same length')

    Jet_Flavors = ['5','4','0']
    eff = {}
    print(f'Total number of jets {len(btagDeepFlavB)}')
    for flavor in Jet_Flavors:
        mask_flav = (hadronFlavour == int(flavor))
        mask_b_tag = (btagDeepFlavB >= LooseWP)
        print(f'   {np.sum(mask_flav)} jet with flavor ' + flavor)
        eff[flavor] = np.sum(mask_flav&mask_b_tag)/np.sum(mask_flav)

print(eff)

# Close the file (optional)
file.close()