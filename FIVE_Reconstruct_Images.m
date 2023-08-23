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
    disp(['Density Compensation for Image ' num2str(i) ' of ' num2str(size(gating,1))]); 
    tmpDCF = Recon.get_DCF(tmptraj,Desired_ImSize,nIter);
    disp(['Density Compensation for Image ' num2str(i) ' of ' num2str(size(gating,1)) ' Complete']); 
    %this_coil_Im = zeros(Desired_ImSize,Desired_ImSize,Desired_ImSize);
    %SOS_Im = zeros(Desired_ImSize,Desired_ImSize,Desired_ImSize);
    Test_Im = zeros(Desired_ImSize,Desired_ImSize,Desired_ImSize);
    Test_Im2 = zeros(Desired_ImSize,Desired_ImSize,Desired_ImSize);
    Test_Im3 = zeros(Desired_ImSize,Desired_ImSize,Desired_ImSize);
    for j = 1:size(fid,3) %Loop over coils
        tmpfid = squeeze(fid(:,:,j));
        tmpfid(:,gating(i,:)==0) = [];
        [this_coil_Im,k] = Recon.mem_eff_recon(Desired_ImSize,tmpfid,tmptraj,tmpDCF,j,size(fid,3),[0 0 0],true,1);
        % map = bart('caldir 24',k);
        if j == 1
            SOS_Im = zeros(size(this_coil_Im));
        end
        SOS_Im = SOS_Im + (abs(this_coil_Im)).^2;
        % Test_Im = Test_Im + this_coil_Im.*map;
        % Test_Im2 = Test_Im2 + (this_coil_Im.*map).^2;
        % Test_Im3 = Test_Im3 + this_coil_Im.*exp(-sqrt(-1)*map);
    end
    %Wrap it up and write out images to NIFTI
    SOS_Im = sqrt(SOS_Im);
    Nifti_Info = AllinOne_Tools.nifti_metadata(SOS_Im,FOV/Desired_ImSize,FOV);
    SOS_Im = fliplr(SOS_Im);
    SOS_Im = flip(SOS_Im,3);
    Nifti_Info.Transform.T(1,1) = 1;
    Nifti_Info.Transform.T(2,2) = 1;
    Nifti_Info.Transform.T(3,3) = 1;
    Nifti_Info.Transform.T(4,4) = 1;
    Nifti_Info.Transform.T(4,1) = -FOV/2;
    Nifti_Info.Transform.T(4,2) = -FOV/2;
    Nifti_Info.Transform.T(4,3) = -FOV/2;
    niftiwrite(SOS_Im,fullfile(write_path,['Image_diaphragm_pos_' num2str(i)]),Nifti_Info,'Compressed',true);
    disp(['Image Recon and Data Writing Complete for Image ' num2str(i) ' of ' num2str(size(gating,1))]); 
end
        
        


