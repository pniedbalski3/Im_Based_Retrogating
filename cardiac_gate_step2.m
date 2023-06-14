function [card_high,card_low] = cardiac_gate_step2(fid,traj,resp_gate,card_gate_low,card_gate_high,OrigImSize,FinImSize)

[myfid,mytraj] = ImTools.change_mat_size(fid,traj,FinImSize,OrigImSize);

CardHigh = logical(resp_gate.*card_gate_high);
CardLow = logical(resp_gate.*card_gate_low);

Imfinal = zeros(FinImSize,FinImSize,FinImSize);
Imfinal2 = zeros(FinImSize,FinImSize,FinImSize);
tmpDCF = Recon.get_DCF(mytraj(:,:,CardHigh),FinImSize,10);
tmpDCF2 = Recon.get_DCF(mytraj(:,:,CardLow),FinImSize,10);
for j = 1:size(fid,3)
    tmp = Recon.mem_eff_recon(FinImSize,myfid(:,CardHigh,j),mytraj(:,:,CardHigh),tmpDCF,j,size(myfid,3),[0 0 0],false);
    tmp2 = Recon.mem_eff_recon(FinImSize,myfid(:,CardLow,j),mytraj(:,:,CardLow),tmpDCF2,j,size(myfid,3),[0 0 0],false);
    Imfinal = Imfinal + tmp.^2;
    Imfinal2 = Imfinal2 + tmp2.^2;
end
card_high = sqrt(Imfinal);
card_low = sqrt(Imfinal2);
