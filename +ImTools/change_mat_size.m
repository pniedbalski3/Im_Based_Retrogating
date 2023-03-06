function [myfid,mytraj] = change_mat_size(fid,traj,New_Size,Orig_ImSize)

if length(New_Size)==1
    mytraj = traj*Orig_ImSize/New_Size;
else
    mytraj(1,:,:) = traj(1,:,:)*Orig_ImSize(1)/New_Size(1);
    mytraj(2,:,:) = traj(2,:,:)*Orig_ImSize(2)/New_Size(2);
    mytraj(3,:,:) = traj(3,:,:)*Orig_ImSize(3)/New_Size(3);
end

rad = squeeze(sqrt(mytraj(1,:,:).^2+mytraj(2,:,:).^2+mytraj(3,:,:).^2));

find_rad = rad>0.5;
%I'm sure there's a faster and better way, but this will work.
n_end = size(fid,1);
for i = 1:size(fid,2)
    too_big = find(find_rad,1,'first');
    if too_big<n_end
        n_end = too_big;
    end
end

myfid = fid(1:n_end,:,:);
mytraj = mytraj(:,1:n_end,:);
