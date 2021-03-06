function ONE_select_slice_coils(fid,traj,Orig_ImSize,Desired_Size,N_short)

%In this function, just recon a small portion of the data for every coil
%element so that we can select which slice best shows diaphragm and which
%coil elements are nearest that region.

if nargin<5
    N_short = 200;
end

mytraj = traj(:,:,1:N_short);
myfid = fid(:,1:N_short,:);

%Try gridding onto 64x64x64 matrix. If that doesn't work, we'll go smaller
%Actually, 96 seems to work pretty well, at least for this radial data

[myfid,mytraj] = ImTools.change_mat_size(myfid,mytraj,Desired_Size,Orig_ImSize);

%% Reconstruct
%Make Image to hold everything
Image = zeros(Desired_Size,Desired_Size,Desired_Size,size(myfid,3)+1);

nIter = 10;
myDCF = get_DCF(mytraj,Desired_Size,nIter);

for i = 1:size(myfid,3)
    Image(:,:,:,i+1) = Recon.mem_eff_recon(Desired_Size,myfid(:,:,i),mytraj,myDCF,i,size(myfid,3));
end
Image(:,:,:,1) = Tools.soscoilcombine(Image); 

imslice(abs(Image));
