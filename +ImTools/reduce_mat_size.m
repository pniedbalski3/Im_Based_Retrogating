function [myfid,mytraj] = change_mat_size(fid,traj,New_Size,Orig_ImSize)

New_Size = 96;
mytraj = traj*Orig_ImSize/New_Size;

rad = squeeze(sqrt(mytraj(1,:,:).^2+mytraj(2,:,:).^2+mytraj(3,:,:).^2));

find_rad = rad>0.5;
%I'm sure there's a faster and better way, but this will work.
n_end = size(myfid,1);
for i = 1:N_short
    too_big = find(find_rad,1,'first');
    if too_big<n_end
        n_end = too_big;
    end
end

myfid = fid(1:n_end,:,:);
mytraj = mytraj(:,1:n_end,:);
