function FIVEA_SoftGate_Images(fid,traj,Orig_ImSize,weights,write_path,FOV,Desired_ImSize)

%% Function to reconstruct images using soft gating based on the diaphragm position and threshold passed by the user 
if nargin <8
    Desired_ImSize = Orig_ImSize;
end

if Desired_ImSize ~= Orig_ImSize
    [fid,traj] = ImTools.change_mat_size(fid,traj,Desired_ImSize,Orig_ImSize);
end

%% Convert weights (in steps of 200 projections) to projection by projection

NPro = size(fid,2);

final_weight = ones(1,NPro);

for i = 1:length(weights)
    final_weight(((i-1)*100+1):((i-1)*100+200)) = weights(i)* final_weight(((i-1)*100+1):((i-1)*100+200));
end

figure('Name','Check Weights')
plot(final_weight);


%% Reconstruct Images

nIter = 10;
for i = 1:length(final_weight)
    fid(:,i,:) = fid(:,i,:)*final_weight(i);
end

DCF = get_DCF(traj,Desired_ImSize,nIter);

SOS_Im = zeros(Desired_ImSize,Desired_ImSize,Desired_ImSize);
for j = 1:size(fid,3) %Loop over coils
    tmpfid = squeeze(fid(:,:,j));
    this_coil_Im = Recon.mem_eff_recon(Desired_ImSize,tmpfid,traj,DCF,j,size(fid,3),[0 0 0],true);
    SOS_Im = SOS_Im + (abs(this_coil_Im)).^2;
end
SOS_Im = sqrt(SOS_Im);
imslice(SOS_Im);
Nifti_Info = AllinOne_Tools.nifti_metadata(SOS_Im,FOV/Desired_ImSize,FOV);
niftiwrite(SOS_Im,fullfile(write_path,'Soft Gated Image'),Nifti_Info,'Compressed',true);




