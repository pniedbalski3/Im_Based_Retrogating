function weight = softgate(d_pos,resp_targ,NProj)

%When soft gating, is it better to convert image-based to projections first
%or do after the weighting? Not sure, but it looks a little better
%weighting and then converting the weights to projections

% first_pos = min(d_pos);
% last_pos = max(d_pos);
% 
% positions = first_pos:last_pos;
% 
% index_count = 1;
% for i = positions
%    locs = find(d_pos==i);
%    for j = 1:length(locs)
%        bfr_locs = (locs(j)-1)*NProj/2; 
%        my_index(index_count,(bfr_locs+1):(bfr_locs+NProj))=1;
%    end
%    index_count = index_count+1;
% end    
% 
% for i = 1:size(my_index,1)
%     my_index(i,:) = my_index(i,:)*i;
% end
% resp_sig = mean(my_index,1);
% sresp_sig = smooth(resp_sig,500);

sresp_sig = d_pos;
TR = .0035;
T_resp = 2e3;% ms
N_resp = ceil(T_resp/TR);
resp_max = max_filter(sresp_sig,N_resp);
resp_min = min_filter(sresp_sig,N_resp);
signal_m = mean((resp_max+resp_min)/2);
signal_s = mean(resp_max-resp_min);


ex_mean = mean(sresp_sig);
ex_mean = resp_targ;
ex_std = std(sresp_sig);

iter = 0;
weight = zeros(size(d_pos));
mweight = 0;
beta = 1;
while mweight<.6 || mweight>.8 && iter<50
    iter = iter +1;
    weight = min(1,exp((beta*ex_std-abs(sresp_sig-ex_mean))/(signal_s/4)));
    mweight = mean(weight);
    if mweight <.6
        beta = beta * 1.1;
    elseif mweight >.8
        beta = beta * .9;
    end
end
disp(['Weights calculated using ' num2str(iter) ' iterations']);


newweight(1:NProj/2) = weight(1);
for i = 2:length(weight)
    start_ind = (NProj/2)*(i-1)+1;
    newweight(start_ind:(start_ind+NProj/2-1)) = mean([weight(i),weight(i-1)]);
    endind = start_ind+NProj/2-1;
end
newweight((endind+1):(endind+NProj/2))=weight(end);
newweight = smooth(newweight,500);

weight = newweight;
