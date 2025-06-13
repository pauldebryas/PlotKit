import ROOT
import os

def change_histograms_names(input_file, output_file, period):
    def process_directory(input_dir, output_dir):
        for key in input_dir.GetListOfKeys():
            obj = key.ReadObj()
            
            # Check if the object is a directory
            if isinstance(obj, ROOT.TDirectory):
                subdir_name = obj.GetName()
                print(f'Processing subdir: {subdir_name}')
                
                # Create output subdirectory and recursively process
                output_subdir = output_dir.mkdir(subdir_name)
                input_dir.cd(subdir_name)
                output_subdir.cd()
                process_directory(obj, output_subdir)
                input_dir.cd()
            
            # Check if the object is a TH1D histogram
            elif isinstance(obj, ROOT.TH1D):
                print(f'Processing histogram: {key.GetName()}')

                # Check if the histogram is in the list of histograms to be renamed
                if key.GetName().startswith('FakeBackground_'):
                    if key.GetName().endswith('Up'):
                        new_name = key.GetName()[:-2] + period + 'Up'
                    elif key.GetName().endswith('Down'):
                        new_name = key.GetName()[:-4] + period + 'Down'
                    else:
                        new_name = key.GetName()  # keep original name if no change
                    print(f'... is being replaced by {new_name}')

                    # Create a new histogram with the new name and fill it with old histogram's data
                    new_hist = ROOT.TH1D(new_name, obj.GetTitle(), obj.GetNbinsX(), obj.GetXaxis().GetXmin(), obj.GetXaxis().GetXmax())
                    new_hist.Add(obj)  # copy contents of obj to new_hist
                    new_hist.SetDirectory(output_dir)
                    new_hist.Write()

                else:
                    # For histograms that don't need renaming, simply clone and write them
                    new_hist = ROOT.TH1D(key.GetName(), obj.GetTitle(), obj.GetNbinsX(), obj.GetXaxis().GetXmin(), obj.GetXaxis().GetXmax())
                    new_hist.Add(obj)  # copy contents of obj to new_hist
                    new_hist.SetDirectory(output_dir)
                    new_hist.Write()

    # Open the input ROOT file
    input_root_file = ROOT.TFile.Open(input_file, "READ")
    if not input_root_file or input_root_file.IsZombie():
        print(f"Error opening file {input_file}")
        return

    # Open the output ROOT file
    output_root_file = ROOT.TFile.Open(output_file, "RECREATE")
    if not output_root_file or output_root_file.IsZombie():
        print(f"Error creating file {output_file}")
        input_root_file.Close()
        return

    # Process the ROOT file starting from the top-level directory
    process_directory(input_root_file, output_root_file)

    # Close the files
    output_root_file.Close()
    input_root_file.Close()

    print(f"Histograms renamed and saved to {output_file}")

# Example usage
period = '2017'
tag = 'AddSyst'
channel = 'ttm'

input_file = f"{os.getenv('RUN_PATH')}/output/{period}/{tag}/{channel}/TH1Hist/TH1_dr_mutau_HNL.root"
output_file = f"{os.getenv('RUN_PATH')}/output/{period}/{tag}/{channel}/TH1Hist/TH1_dr_mutau_HNL_new.root"

change_histograms_names(input_file, output_file, period)
