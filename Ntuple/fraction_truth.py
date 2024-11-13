import ROOT

# Open the ROOT file
file_path = "Truth_eta/Truth_aQGC/WmZ_llqq_ntuple_truth.root"

def calculate_fractions(file_path, variable_names):
    import ROOT

    # Open the ROOT file
    root_file = ROOT.TFile(file_path, "READ")

    # Get the TTree
    tree = root_file.Get("Truth")
    for variable_name in variable_names:
        
        # Check if the TTree exists
        if not tree:
            print(f"Error: No TTree named 'Truth' found in file {file_path}.")
            root_file.Close()
            return

        # Create histograms
        histogram0 = ROOT.TH1D("histogram0", "title", 1, -0.5, 0.5)
        histogram1 = ROOT.TH1D("histogram1", "title", 1, 0.5, 1.5)
        histogram2 = ROOT.TH1D("histogram2", "title", 1, 1.5, 2.5)

        # Fill the histograms with the values of the variable where it is 0, 1, or 2
        tree.Draw(f"{variable_name}>>histogram0", f"{variable_name} == 0")
        tree.Draw(f"{variable_name}>>histogram1", f"{variable_name} == 1")
        tree.Draw(f"{variable_name}>>histogram2", f"{variable_name} == 2")

        # Get the number of entries in the histograms
        num_entries0 = histogram0.GetEntries()
        num_entries1 = histogram1.GetEntries()
        num_entries2 = histogram2.GetEntries()

        # Get the total number of events in the tree
        total_events = num_entries0 + num_entries1 + num_entries2

        # Calculate the fractions of events
        fraction0 = num_entries0 / total_events if total_events != 0 else 0
        fraction1 = num_entries1 / total_events if total_events != 0 else 0
        fraction2 = num_entries2 / total_events if total_events != 0 else 0

        print(f"Fraction of events where {variable_name} is 0: {fraction0}")
        print(f"Fraction of events where {variable_name} is 1: {fraction1}")
        print(f"Fraction of events where {variable_name} is 2: {fraction2}")

    # Close the ROOT file
    root_file.Close()

file_path = "Truth_eta/Truth_eta_40/Truth_SM/WmZ_llqq_ntuple_truth.root"
variable_name = ["FromTruth_W_q_misIDtagjets","misID_tag_jets"]
calculate_fractions(file_path, variable_name)