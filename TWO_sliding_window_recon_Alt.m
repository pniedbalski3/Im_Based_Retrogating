function All_Im = TWO_sliding_window_recon_Alt(fid,traj,Orig_ImSize,slice,dim,coils,Small_size,Sliding_window_size)

%Function to do a sliding window and reconstruct tons! of images. These
%images are 3D, and to save >1000 would be insane. Just save the slice that
%is passed by slice. To save time, only reconstruct the coils closest to
%the diaphragm, which is specified by coils.
if nargin < 8
    Sliding_window_size = 200; %200 seems to work reasonably well.
end

fid = fid(:,:,coils);

if Orig_ImSize ~= Small_size
    [myfid,mytraj] = ImTools.change_mat_size(fid,traj,Small_size,Orig_ImSize);
else
    myfid = fid;
    mytraj = traj;
end

%Overlap by 50% - This may be too slow... We'll see.

%get total number of Images in a lazy way
Ind_Ims_Sw = 1:(Sliding_window_size/2):(size(fid,2)-Sliding_window_size);
Ims = zeros(Small_size,Small_size,length(Ind_Ims_Sw));
tmpIms = zeros(Small_size,Small_size,Small_size,length(coils));
nIter = 10;
counter = 1;
for i = Ind_Ims_Sw
    disp(['Starting Recon ' num2str(counter) ' of ' num2str(length(Ind_Ims_Sw))]);
    tmpfid = myfid(:,i:(i+Sliding_window_size),:);
    tmptraj = mytraj(:,:,i:(i+Sliding_window_size));
    tmptraj = [reshape(tmptraj(1,:,:),1,[])' reshape(tmptraj(2,:,:),1,[])' reshape(tmptraj(3,:,:),1,[])'];
    tmpDCF = Recon.get_DCF_Robertson(tmptraj,[1 1 1]*Small_size,nIter);
    for j = 1:length(coils)
        tmpfid2 = tmpfid(:,:,j);
        tmpfid2 = reshape(tmpfid2,1,[])';
        tmpIms(:,:,:,j) = tmpDCF.reconstruct(tmpfid2,tmptraj);
        %tmpIms(:,:,:,j) = Recon.mem_eff_recon([1 1 1]*Small_size,tmpfid2,tmptraj,tmpDCF,j,size(myfid,3),[0 0 0],false);
    end
    if length(coils)>1
        tmpIm = Tools.soscoilcombine(tmpIms);
    else
        tmpIm = tmpIms;
    end
    if dim==1
        Ims(:,:,counter) = squeeze(tmpIm(slice,:,:));
    elseif dim == 2
        Ims(:,:,counter) = squeeze(tmpIm(:,slice,:));
    elseif dim == 3
        Ims(:,:,counter) = squeeze(tmpIm(:,:,slice));
    else
        error('dim needs to be 1,2, or 3');
    end
    disp(['Done with Recon ' num2str(counter) ' of ' num2str(length(Ind_Ims_Sw))]);
    counter = counter + 1;
end

All_Im = Ims;