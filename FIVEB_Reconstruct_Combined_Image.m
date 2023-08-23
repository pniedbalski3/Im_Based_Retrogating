function FIVEB_Reconstruct_Combined_Image(fid,traj,Orig_ImSize,gating,write_path,FOV,comb_bins,Desired_ImSize,twix_obj,Patient_Name)

%% Function to reconstruct images! 
if nargin <8
    Desired_ImSize = Orig_ImSize;
end

if Desired_ImSize ~= Orig_ImSize
    [fid,traj] = ImTools.change_mat_size(fid,traj,Desired_ImSize,Orig_ImSize);
end


%% To avoid potential memory issues, just do one channel at a time and sum of squares combine as I go... Kind of sloppy, but won't cause memory errors
nIter = 10;

which_bins = zeros(1,size(gating,2));

for i = comb_bins
    which_bins(gating(i,:)~=0) = 1;
end

tmptraj = traj;
tmptraj(:,:,which_bins==0) = [];
tmpDCF = Recon.get_DCF(tmptraj,Desired_ImSize,nIter);
%SOS_Im = zeros(Desired_ImSize,Desired_ImSize,Desired_ImSize);
for j = 1:size(fid,3) %Loop over coils
    tmpfid = squeeze(fid(:,:,j));
    tmpfid(:,which_bins==0) = [];
    this_coil_Im = Recon.mem_eff_recon(Desired_ImSize,tmpfid,tmptraj,tmpDCF,j,size(fid,3),[0 0 0],true,1);
    if j == 1
        SOS_Im = zeros(size(this_coil_Im));
    end
    SOS_Im = SOS_Im + (abs(this_coil_Im)).^2;
end
%Wrap it up and write out images to NIFTI
SOS_Im = sqrt(SOS_Im);

imslice(abs(SOS_Im/max(SOS_Im(:))));

niftiwrite(SOS_Im,fullfile(write_path,'End_Expiration_Image'),'Compressed',true);
% [N4_Im,~] = GasExchangeV3.Dissolved_ProtonBiasCorrection(SOS_Im);
% N4_Im = N4_Im/max(N4_Im(:));


%% Write to DICOM first before rotating image for proper Nifti orientation.
if ~isfolder(fullfile(write_path,'Exp_DICOM'))
    mkdir(fullfile(write_path,'Exp_DICOM'))
end
Info = twix_object_to_dicom_metadata(twix_obj);

%Need to fix Info for different resolution and slice thickness than
%specified in the the twix file.
Info.PatientName = Patient_Name;
Info.PatientID = Patient_Name;
Info.SliceThickness = FOV/Desired_ImSize;
Info.SpacingBetweenSlices = FOV/Desired_ImSize;
Info.PixelSpacing = [FOV/Desired_ImSize FOV/Desired_ImSize];

Slice_Spacing = FOV/Desired_ImSize;

%I'm definitely going to need to rotate images...

%I'm not sure that Image Position Patient and Image Orientation Patient are
%correct... But can at least give this a shot.
Slice_Loc = FOV/2;
for i = 1:Desired_ImSize
    Info.SliceLocation = Slice_Loc - (i-1)*Slice_Spacing;
    Info.ImagePositionPatient   = [-FOV/2,FOV/2,Info.SliceLocation];
    %Info.Private_0019_1015      = slicePos(1:3,1);
    %Info.Private_0019_1016      = slicePos(1:3,1);
    
    myslice = squeeze(SOS_Im(:,:,i));
    dicomwrite(myslice,fullfile(write_path,'Exp_DICOM',['dcm_Img_slc_' num2str(i,'%0.3d')]),Info);
end
    

%% Write Image - This should be canonical Orientation.
%Nifti_Info = AllinOne_Tools.nifti_metadata(SOS_Im,FOV/Desired_ImSize,FOV);
% SOS_Im = fliplr(SOS_Im);
% SOS_Im = flip(SOS_Im,3);
% Nifti_Info = AllinOne_Tools.nifti_metadata(SOS_Im,FOV/Desired_ImSize,FOV);
% Nifti_Info.Transform.T(1,1) = 1;
% Nifti_Info.Transform.T(2,2) = 1;
% Nifti_Info.Transform.T(3,3) = 1;
% Nifti_Info.Transform.T(4,4) = 1;
% Nifti_Info.Transform.T(4,1) = -FOV/2;
% Nifti_Info.Transform.T(4,2) = -FOV/2;
% Nifti_Info.Transform.T(4,3) = -FOV/2;
% 
% niftiwrite(SOS_Im,fullfile(write_path,'End_Expiration_Image'),Nifti_Info,'Compressed',true);


    
    
    
    
        
        


