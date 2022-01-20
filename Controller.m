%% Series of sections to call my image-based retrospective gating code

%% First, Identify the slice and coil elements

%Pass fids, trajectories, the collected image size (will be subsampled to
%96x96x96 matrix), and the number of projections to use (200 seems to be
%adequate - makes for a fast recon)
%Orig_ImSize = 268;
%Desired_Size = 96;
%NProj = 200;

Orig_ImSize = 180;
Desired_Size = 90;
NProj = 5000;

ONE_select_slice_coils(fid,traj,Orig_ImSize,Desired_Size,NProj) 

%% Now, inspect the images to determine the desired slice and coil elements

%You are about to reconstruct >1000 images, so you want to choose only 1 or
%2 coil elements so that it doesn't take days to do. Also, it would be
%infeasible to save all >1000 3D images, so just pick the slice for gating
%(a coronal slice with very clear diaphragm)
slice = 64; %which slice you want
dim = 2; %Which dimension (from imslice)
coils = [3]; %Pick 1 or 2 coil elements that really show the diaphragm

%About to reconstruct images with the number of projections specified by
%Sliding_Window. There will be 50% overlap of these image.
Sliding_Window = NProj;

All_Im = TWO_sliding_window_recon(fid,traj,Orig_ImSize,slice,dim,coils,Desired_Size,Sliding_Window);

figure('Name','Find indices of line over diaphram')
imagesc(squeeze(All_Im(:,:,1)))
colormap(gray)
axis square;
axis off;

%% Next, visualize Respiratory motion and bin diaphragm posisitons

% Figure out the x and y components of a line going over the diaphragm
X = [46 59];
Y = [70 70];

diaphragm_pos = THREE_diaphragm_motion(All_Im,X,Y);

%% Get indices of gating

%From the diaphragm positions found in the previous step, create an array
%holding projection indices for each diaphram position.
my_index = FOUR_get_gated_indices(diaphragm_pos,Sliding_Window,fid);

figure('Name','Check_Gating_Indices')
imagesc(my_index);

%% Finally, reconstruct images for every binned diaphragm position - This will give several images 
%Note, take quantitative values with a grain of salt, because each image
%will be created with a different number of projections. 
FOV = 400;

% If you want lower resolution (but more fully sampled) images, you can
% pass a Subsampled_ImSize. If you don't pass anything, it will just
% default to the original image size.

Subsampled_ImSize = Orig_ImSize;
Subsampled_ImSize = 200;

FIVE_Reconstruct_Images(fid,traj,Orig_ImSize,my_index,pwd,FOV,Subsampled_ImSize)

