function FIVEC_Write_Data_Coords(fid,traj,Orig_ImSize,gating,outfname,Desired_ImSize,FOV,myind)

%% Function to reconstruct images! 
if nargin <7
    Desired_ImSize = Orig_ImSize;
end

if Desired_ImSize ~= Orig_ImSize
    [fid,traj] = ImTools.change_mat_size(fid,traj,Desired_ImSize,Orig_ImSize);
end


%% To avoid potential memory issues, just do one channel at a time and sum of squares combine as I go... Kind of sloppy, but won't cause memory errors

%for i = 1:size(gating,1) %Loop over diaphragm positions
    %Get trajectory the right size and do DCF only once!
tmptraj = traj;
tmptraj(:,:,gating(myind,:)==0) = [];
tmpfid = fid;
tmpfid(:,gating(myind,:)==0,:) = [];
nIter = 10;
myDCF = Recon.get_DCF(tmptraj,Orig_ImSize,nIter);
tmptraj = tmptraj*FOV;
bartdcf = reshape(myDCF,[1 size(myDCF)]);
bartfid = reshape(tmpfid,[1 size(tmpfid)]);

writecfl([outfname,'_traj'],tmptraj);
writecfl([outfname,'_dcf'],bartdcf);
writecfl([outfname,'_data'],bartfid);
writecfl([outfname,'_dcf2'],sqrt(bartdcf));
   % save_for_cs(tmpfid,tmptraj,fullfile(write_path,['Raw_GatingPosition_' num2str(i) '.mat']));
%end
        
        


