import glob
import os
from array import array
import ROOT
import itertools
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
import pandas as pd
import numpy as np
import pyAMI.client
import pyAMI.atlas.api as AtlasAPI
client = pyAMI.client.Client('atlas')
AtlasAPI.init()

def extract_EFT_op_proces_dec(conf):
    """
    This function takes conf as a parameter,
    and returns the EFT operator, process, and decay.
    """
    # Split the conf name with the underscore 
    parts = conf.split('_')
    EFT_op,EFT_type, proc, decay = parts[3], parts[4], parts[1], parts[2]
    
    return EFT_op, EFT_type, proc, decay


def cross_section_fb(EFT_op, EFT_type,  proces, dec):
    """
    This function takes an EFT operator, a process, and a decay as parameters,
    and returns the cross section for the corresponding dataset.
    """
    all_ops = ["FM0","FM1","FM2","FM3","FM4","FM5","FM7",
               "FS02","FS1",
               "FT0","FT1","FT2","FT5","FT6","FT7","FT8","FT9"]
    all_ops2 = ["FM","FS","FT"]
    Processes = [proces]
    Decay = [dec]

    organized_info = {}

    if EFT_type == 'QUAD':
        for op in all_ops2:
            for process in Processes:
                for decay in Decay:
                    pattern = f'%MGPy8EG_aQGC{op}%_QUAD_1_{process}_{decay}%'
                    List_ami = AtlasAPI.list_datasets(client, patterns=[pattern], 
                                                    fields=['ldn', 'cross_section', 'dataset_number'], 
                                                    limit=[1, 100], type='EVNT', project='mc16%')

                    if List_ami:
                        for dataset in List_ami:
                            operator = next((op for op in all_ops if op in dataset['ldn']), None)
                            if operator:
                                key = f"{operator}_{EFT_type}_{process}_{decay}"
                                if key not in organized_info:
                                    organized_info[key] = {'name': [], 'cross_section': [], 'dataset_number': []}

                                organized_info[key]['name'].append(dataset['ldn'])
                                organized_info[key]['cross_section'].append(dataset['cross_section']) # Cross section in nb
                                organized_info[key]['dataset_number'].append(dataset['dataset_number'])
                    else:
                        print(f'No datasets found for {op}_{process}_{decay}\n')
                        
    elif (EFT_type == 'SM'):
        for process in Processes:
            for decay in Decay:
                if process == 'WmZ' or process == 'WpZ':
                    process_ = 'WZ'
                pattern = f'%WZjj_llqq_EW6%'
                List_ami = AtlasAPI.list_datasets(client, patterns=[pattern], 
                                                fields=['ldn', 'cross_section', 'dataset_number'], 
                                                limit=[1, 100], type='EVNT', project='mc16%')
                #print(List_ami)
                if List_ami:
                    for dataset in List_ami:
                        key = f"{EFT_op}_{EFT_type}_{process}_{decay}"
                        if key not in organized_info:
                            organized_info[key] = {'name': [], 'cross_section': [], 'dataset_number': []}

                        organized_info[key]['name'].append(dataset['ldn'])
                        organized_info[key]['cross_section'].append(dataset['cross_section'])

    # Construct the key from the EFT operator, process, and decay
    key_ = f"{EFT_op}_{EFT_type}_{proces}_{dec}"

    # Check if the key exists in the organized_info dictionary
    if key_ in organized_info:
        # If the key exists, return the first cross section
        xsection_nb = organized_info[key_]['cross_section'][0] # Return the cross section in nb
        xsection_fb = float(xsection_nb) * 1e6 # Convert the cross section to fb
        return xsection_fb
    else:
        # If the key does not exist, print an error message and return None
        print(organized_info[key_]['cross_section'][0])
        print(f'Error: No datasets found for {key_}')
        return None



def plot_histograms(output_plot,desired_num_bins, file_path, label):
    # Open the first ROOT file to retrieve the list of parameters (branches)
    if not os.path.exists(output_plot):
        os.makedirs(output_plot)
    first_root_file = ROOT.TFile(file_path, "READ")
    keys = [key.GetName() for key in first_root_file.GetListOfKeys()]
    print(keys)
    first_root_file.Close()

    for parameter_to_plot in keys:    
        # Create a new TCanvas
        canvas_name = "canvas_" + parameter_to_plot
        canvas = ROOT.TCanvas(canvas_name, "Stacked Histograms", 800, 600)

        # Create a THStack
        hs = ROOT.THStack("hs", "Distribution of " + parameter_to_plot)
        
        
        # Loop over the files to retrieve and stack the histograms

        root_file = ROOT.TFile(file_path, "READ")
    
        # Retrieve the histogram
        histogram = root_file.Get(parameter_to_plot)
        histogram.SetDirectory(0)  # Detach from the file
        current_num_bins = histogram.GetNbinsX()
        
        rebin_factor = int(current_num_bins / desired_num_bins)
        if rebin_factor > 0:
            histogram.Rebin(rebin_factor)
            
         
        # Add the histogram to the stack
        hs.Add(histogram)
        


        legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)  # x1, y1, x2, y2

        # Set legend style
        legend.SetFillColor(0) 
        legend.SetFillStyle(0)  
        legend.SetBorderSize(0)  
        legend.SetTextSize(0.03)  

        # Add entry to the legend
        legend.AddEntry(histogram, label, "l")
    
        # Close the ROOT file
        root_file.Close()
        
        # Draw the histograms
        #hs.Draw("nostack")
        hs.Draw("hist")
        hs.GetXaxis().SetTitle(parameter_to_plot )
        hs.GetYaxis().SetTitle("Events / bin")
        
        # Set log scale for y-axis
        canvas.SetLogy()
        
        # Draw the legend
        legend.Draw()
        canvas.Draw()
        
        # Update the canvas
        canvas.Update()
        
        # Save the canvas to a file
        canvas.SaveAs(output_plot+parameter_to_plot + "_hist.png")
