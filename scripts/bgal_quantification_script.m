%% Script for quantifying beta-gal + signal validation
dependencies = 'C:\Users\humza\Desktop\github\scripts\dependencies\'; %functions required for segmentation
currFile = 'C:\Users\humza\Desktop\github\raw_data\wella06_channelbfp,gfp,cy5,brightfield_seq0005xy1.tif'; %path of current file
datadirIF = 'C:\Users\humza\Desktop\github\exported\'; %directory to save exported data
mask_dir = 'C:\Users\humza\Desktop\github\mask_dir\'; %directory to save mask files + image validation

% Specify channels in tiff stack 
num_stacks = 6; %number of stacks in tiff image
nuc_channel = 1; %channel for segmenting nucleus
cyto_channel = 2; %channel for segmenting cytoplasm
R = 6; %red channel of RGB image
G = 5; %green channel of RGB image
B = 4; %blue channel of RGB image
cytothresh = 260; %threshold for cytoplasmic binarization
bgal_pct = 5; %percentile of red channel signal distribution to assign beta-gal value

% Execute beta-gal quantification function
IF_nd2_bgal_github(dependencies, currFile, datadirIF, mask_dir, num_stacks, nuc_channel, cyto_channel, R, G, B, cytothresh, bgal_pct)

% Execute beta-gal quantification validation function
betagal_quantification_validation(currFile, datadirIF, mask_dir)