function FIVE_Reconstruct_Images(fid,traj,Orig_ImSize,gating,write_path,FOV,Desired_ImSize)

%% Function to reconstruct images! 
if nargin <7
    Desired_ImSize = Orig_ImSize;
end

if Desired_ImSize ~= Orig_ImSize
    [fid,traj] = ImTools.change_mat_size(fid,traj,Desired_ImSize,Orig_ImSize);
end

%% To avoid potential memory issues, just do one channel at a time and sum of squares combine as I go... Kind of sloppy, but won't cause memory errors
nIter = 10;

for i = 1:size(gating,1) %Loop over diaphragm positions
    %Get trajectory the right size and do DCF only once!
    tmptraj = traj;
    tmptraj(:,:,gating(i,:)==0) = [];
    tmpDCF = get_DCF(tmptraj,Desired_ImSize,nIter);
    
    this_coil_Im = zeros(Desired_ImSize,Desired_ImSize,Desired_ImSize);
    SOS_Im = zeros(Desired_ImSize,Desired_ImSize,Desired_ImSize);
    for j = 1:size(fid,3) %Loop over coils
        tmpfid = squeeze(fid(:,:,j));
        tmpfid(:,gating(i,:)==0) = [];
        this_coil_Im = Recon.mem_eff_recon(Desired_ImSize,tmpfid,tmptraj,tmpDCF,j,size(fid,3),[0 0 0],true);
        SOS_Im = SOS_Im + (abs(this_coil_Im)).^2;
    end
    %Wrap it up and write out images to NIFTI
    SOS_Im = sqrt(SOS_Im);
    Nifti_Info = AllinOne_Tools.nifti_metadata(SOS_Im,FOV/Desired_ImSize,FOV);
    niftiwrite(SOS_Im,fullfile(write_path,['Image_diaphragm_pos_' num2str(i)]),Nifti_Info,'Compressed',true);
end
        
        


