function FIVEC_Write_Data_Coords(fid,traj,Orig_ImSize,gating,write_path,Desired_ImSize)

%% Function to reconstruct images! 
if nargin <7
    Desired_ImSize = Orig_ImSize;
end

if Desired_ImSize ~= Orig_ImSize
    [fid,traj] = ImTools.change_mat_size(fid,traj,Desired_ImSize,Orig_ImSize);
end


%% To avoid potential memory issues, just do one channel at a time and sum of squares combine as I go... Kind of sloppy, but won't cause memory errors

for i = 1:size(gating,1) %Loop over diaphragm positions
    %Get trajectory the right size and do DCF only once!
    tmptraj = single(traj);
    tmptraj(:,:,gating(i,:)==0) = [];
    tmpfid = single(fid);
    tmpfid(:,gating(i,:)==0,:) = [];
    
    save_for_cs(tmpfid,tmptraj,fullfile(write_path,['Raw_GatingPosition_' num2str(i) '.mat']));
end
        
        


