import os

#params
channels = ['tmm', 'tee','tem','ttm','tte']
period = '2018'
tag = 'AddJETcorr'
region = 'InvertedBjetsVetoRegion'

# Define the input list
dict = {
    'tmm':[['FPP'], ['PFP', 'PPF', 'PFF'], ['FFP', 'FFF', 'FPF']],
    'tee':[['FPP'], ['PFP', 'PPF', 'PFF'], ['FFP', 'FFF', 'FPF']],
    'tem':[['FPP'], ['PFP'], ['PPF'], ['FFP','PFF','FFF', 'FPF']],
    'ttm':[['FPP', 'PFP', 'FFP'], ['PPF'], ['PFF', 'FFF', 'FPF']],
    'tte':[['FPP', 'PFP', 'FFP'], ['PPF'], ['PFF', 'FFF', 'FPF']],
}
label_dict = {
    'tmm':['fake tau', 'fake muon', 'cross contribution'],
    'tee':['fake tau', 'fake electron', 'cross contribution'],
    'tem':['fake tau', 'fake electron', 'fake muon', 'cross contribution'],
    'ttm':['fake tau', 'fake muon', 'cross contribution'],
    'tte':['fake tau', 'fake electron', 'cross contribution'],
}

print('')
for channel in channels:
    print(f'channel {channel}')
    my_list = dict[channel]
    label = label_dict[channel]

    # Read the file path from arguments
    file_path = os.path.join(os.getenv("RUN_PATH"), f'D_RootHist/results/{period}/{tag}/{channel}/AR_info_{region}.txt') 

    # Parse the text data
    regions = {}

    with open(file_path, 'r') as file:
        data_lines = file.readlines()

    current_region = None
    for line in data_lines:
        line = line.strip()
        if "region:" in line:
            current_region = line.split(" ")[0]
        elif "FFs * FakesProp" in line:
            value = float(line.split("=")[1].strip())
            regions[current_region] = value

    # Compute total sum
    FFs_times_FakesProp_Total = sum(regions.values())

    # Compute proportions
    proportions = []
    for group in my_list:
        group_sum = sum(regions[region] for region in group)
        proportion = (group_sum / FFs_times_FakesProp_Total) * 100
        proportions.append(proportion)

    for i in range(len(my_list)):
        print(f'{label[i]} ({my_list[i]}) : {str(round(proportions[i],2))} %')
    
    print('')