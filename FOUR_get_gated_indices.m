function [my_index] = FOUR_get_gated_indices(diaphragm_pos,sliding_window,fid)

%Get the indices of each diaphragm position

%Just want the number of projections from the fid
NPro = size(fid,2);

first_pos = min(diaphragm_pos);
last_pos = max(diaphragm_pos);

positions = first_pos:last_pos;

my_index = zeros(length(positions),NPro);

index_count = 1;
for i = positions
   locs = find(diaphragm_pos==i);
   for j = 1:length(locs)
       bfr_locs = (locs(j)-1)*sliding_window/2; 
       my_index(index_count,(bfr_locs+1):(bfr_locs+sliding_window))=1;
   end
   index_count = index_count+1;
end