%% Series of sections to call my image-based retrospective gating code
%% Load Radial Data

[fid,traj,ImSize,FOV,mypath,twix_obj] = load_radial_data();

%% First, Identify the slice and coil elements

%Pass fids, trajectories, the collected image size (will be subsampled to
%96x96x96 matrix), and the number of projections to use (200 seems to be
%adequate - makes for a fast recon)
%Orig_ImSize = 268;
Desired_Size = 96;
NProj = 200;

Orig_ImSize = ImSize;

%Orig_ImSize = ImSize;
% Desired_Size = 96;
% NProj = 200;

ONE_select_slice_coils(fid,traj,Orig_ImSize,Desired_Size,NProj) 

%% Now, inspect the images to determine the desired slice and coil elements

%You are about to reconstruct >1000 images, so you want to choose only 1 or
%2 coil elements so that it doesn't take days to do. Also, it would be
%infeasible to save all >1000 3D images, so just pick the slice for gating
%(a coronal slice with very clear diaphragm)
slice = 59; %which slice you want
dim = 2; %Which dimension (from imslice)
coils = [14 18]; %Pick 1 or 2 coil elements that really show the diaphragm

%About to reconstruct images with the number of projections specified by
%Sliding_Window. There will be 50% overlap of these image.
Sliding_Window = NProj;

All_Im = TWO_sliding_window_recon(fid,traj,Orig_ImSize,slice,dim,coils,Desired_Size,Sliding_Window);

All_Im = abs(All_Im);
figure('Name','Find indices of line over diaphram')
imagesc(squeeze(abs(All_Im(:,:,1))))
colormap(gray)
axis square;
axis off;

%% Next, visualize Respiratory motion and bin diaphragm posisitons

% Figure out the x and y components of a line going over the diaphragm
X = [42 75];
Y = [61 61];

diaphragm_pos = THREE_diaphragm_motion(All_Im,X,Y);


%% Get indices of gating

%From the diaphragm positions found in the previous step, create an array
%holding projection indices for each diaphram position.
my_index = FOUR_get_gated_indices(diaphragm_pos,Sliding_Window,fid);

numbinned = sum(my_index,2);
my_index(numbinned<12000,:) = [];

my_index(1,:) = [];
my_index(1,:) = [];
my_index(4,:) = [];
my_index(4,:) = [];

numbinned = sum(my_index,2);
[~,final_gate] = max(numbinned);

figure('Name','Check_Gating_Indices')
imagesc(my_index);

save(fullfile(mypath,'Gating_Index.mat'),'my_index');
%% Write out data for BART recon
FIVEC_Write_Data_Coords(fid,traj,Orig_ImSize,my_index,'Exp_data',Orig_ImSize,FOV,final_gate);


%% Cardiac-Gated Imaging?
[chigh,clow] = cardiac_gate(fid(:,:,24),3.5);
[ImHigh,ImLow] = cardiac_gate_step2(fid,traj,my_index(6,:),chigh,clow,400,128);

%% Finally, reconstruct images for every binned diaphragm position - This will give several images 
%Note, take quantitative values with a grain of salt, because each image
%will be created with a different number of projections. 
%FOV = 400;

% If you want lower resolution (but more fully sampled) images, you can
% pass a Subsampled_ImSize. If you don't pass anything, it will just
% default to the original image size.

%I think it makes the most sense to reconstruct the different diaphragm
%bins to a 200 isotropic matrix size (2 mm resolution), and only recon the
%end expiration image at 1 mm resolution

Subsampled_ImSize = Orig_ImSize;
%Subsampled_ImSize = 200;

FIVE_Reconstruct_Images(fid,traj,Orig_ImSize,my_index,mypath,FOV,Subsampled_ImSize)

% Combine a couple of expiratory bins to get a more fully sampled image
combinebins = [4 5];

%FOV = 400;
Subsampled_ImSize = Orig_ImSize;

Patient_Name = 'Xe-052';

FIVEB_Reconstruct_Combined_Image(fid,traj,Orig_ImSize,my_index,mypath,FOV,combinebins,Subsampled_ImSize,twix_obj,Patient_Name)


%% Soft gating Reconstruction - This doesn't work right now. Maybe revisit sometime
approximate_full_exhale = 6;
approximate_full_inhale = 10;

threshold = (approximate_full_inhale-approximate_full_exhale)/4 + approximate_full_exhale;

weights = zeros(size(diaphragm_pos));

weights(diaphragm_pos == approximate_full_exhale) = 1;
weights(diaphragm_pos == approximate_full_exhale+1) = 1;
weights(diaphragm_pos == approximate_full_exhale-1) = 1;

FOV = 400;

Low_res = FOV/Desired_Size;
diaphragm_pos = (diaphragm_pos-approximate_full_exhale)*Low_res;

alpha = 3/((approximate_full_inhale-approximate_full_exhale)*Low_res);

for i = 1:length(diaphragm_pos)
    if weights(i) == 0
        weights(i) = exp(-alpha*(diaphragm_pos(i)));
        if diaphragm_pos(i) < -Low_res
            weights(i) = 0;
        end
    end
end

figure('Name','Evaluate Weighting')
subplot(2,1,1)
plot(diaphragm_pos)

subplot(2,1,2)
plot(weights);

Subsampled_ImSize = Orig_ImSize;

FIVEA_SoftGate_Images(fid,traj,Orig_ImSize,weights,pwd,FOV,Subsampled_ImSize)
