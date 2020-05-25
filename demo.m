clc; clf; clear;
addpath(genpath('utils/'));

%%
mesh_dir =  './data_mesh/';
num_model = 1;

if_output = true;
output_dir = './data_wavelet/';
collection_type = 2;  % 1;

% the number of eigenfunctions and scales for WEDS
numEigs = 300;
weds_scale = 128;


ref = false;  % generating desc for partial matching
if ref   
    ref_name = 'cuts_victoria_shape_0.off'; % for example
    ref_fullname = [mesh_dir, '/', ref_name];
    if ref_name(end-2:end) == 'off'
        shape_ref=read_shape(ref_fullname);
    elseif ref_name(end-2:end) == 'obj'
        shape_ref=read_shape_obj(ref_fullname);
    elseif ref_name(end-2:end) == 'ply'
        shape_ref=read_shape_ply(ref_fullname);
    end
else
    shape_ref = 0;
end

parfor k = 0:num_model-1
    s_name = ['tr_reg_', num2str(k, '%03d'), '.off'];
    WEDS = waveletEnergyDecompositionSignature(mesh_dir, s_name, numEigs, weds_scale, ref, shape_ref, if_output, output_dir, collection_type);
end