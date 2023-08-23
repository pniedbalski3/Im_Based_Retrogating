function [card_high,card_low] = cardiac_gate_step2(mypath,fid,traj,resp_gate,card_gate,OrigImSize,FinImSize,FOV)
%%
[myfid,mytraj] = ImTools.change_mat_size(fid,traj,FinImSize,OrigImSize);

nbins = length(unique(card_gate))-1;

Allims = zeros(FinImSize,FinImSize,FinImSize,nbins);
%%
for i = 1:nbins
    tmpcardgate = card_gate;
    if i > ((nbins+2)/2)
        ind = (-1*(i-(nbins+2)/2+1));
        tmpcardgate(tmpcardgate ~= ind) = 0;
        tmpcardgate(tmpcardgate ~= 0) = 1;
    else
        tmpcardgate(tmpcardgate ~= i) = 0;
        tmpcardgate(tmpcardgate~= 0) = 1;
    end
    tmpcardgate_final = logical(resp_gate.*tmpcardgate);
    tmptraj = mytraj;
    tmptraj(:,:,~tmpcardgate_final) = [];
    tmpfid = myfid;
    tmpfid(:,~tmpcardgate_final,:) = [];
    nIter = 10;
    myDCF = Recon.get_DCF(tmptraj,FinImSize,nIter);
    tmptraj = tmptraj*FinImSize;
    bartdcf = reshape(myDCF,[1 size(myDCF)]);
    bartfid = reshape(tmpfid,[1 size(tmpfid)]);
    if ~isfolder(fullfile('/home/pniedbalski@kumc.edu/Documents/Outputs','Cardiac_Gate_Tmp'))
        mkdir(fullfile('/home/pniedbalski@kumc.edu/Documents/Outputs','Cardiac_Gate_Tmp'));
    end
    outfname = fullfile('/home/pniedbalski@kumc.edu/Documents/Outputs','Cardiac_Gate_Tmp',['Dual_Gated_bin_' num2str(i)]);
    
    writecfl([outfname,'_traj'],tmptraj);
    writecfl([outfname,'_dcf'],bartdcf);
    writecfl([outfname,'_data'],bartfid);
    writecfl([outfname,'_dcf2'],sqrt(bartdcf));

    %copyfile(fullfile(mypath,['Dual_Gated_bin_' num2str(i)]), '/home/pniedbalski@kumc.edu/Documents/Outputs');
    mydir = fullfile('/home/pniedbalski@kumc.edu/Documents/Outputs','Cardiac_Gate_Tmp');
        
    cd '/home/XeAnalysis/bart'
    file_name = fullfile(mydir,['Dual_Gated_bin_' num2str(i)]);

    dims = [num2str(FinImSize) ':' num2str(FinImSize) ':' num2str(FinImSize) ' '];
    % Coil Combine and Reconstruct
    mycmd = system(['./bart fmac ' file_name '_data ' file_name '_dcf ' file_name '_datac']) 
    %Base Recon
    mycmd = system(['./bart nufft -a -d ' dims file_name '_traj ' file_name '_datac ' file_name '_imgL'])
    % FFT
    mycmd = system(['./bart fft 7 ' file_name '_imgL ' file_name '_ksp'])
    % Coil Combination - Use 8 coils to make sure we have adequate memory
    mycmd = system(['./bart cc -M ' file_name '_ksp ' file_name '_ccMatrix'])
    mycmd = system(['./bart ccapply -p 8 ' file_name '_data ' file_name '_ccMatrix ' file_name '_data1'])
    mycmd = system(['./bart ccapply -p 8 ' file_name '_ksp ' file_name '_ccMatrix ' file_name '_ksp1'])
    
    % Estimate coil sensitivies from k-spac>getNextLine (line 38)
    mycmd = system(['./bart caldir 24 ' file_name '_ksp1 ' file_name '_mapsL'])
    
    tic
    mycmd = system(['export OMP_NUM_THREADS=32; ./bart pics -C 50 -i 80 -R T:7:0:0.001 -p ' file_name '_dcf2 -t ' file_name '_traj ' file_name '_data1 ' file_name '_mapsL ' file_name '_picsrecon'])
    toc

    cd(mydir)
    Im = readcfl([file_name '_picsrecon']);
    Im = flip(flip(flip(Im,1),2),3);
    niftiwrite(abs(Im),['Dual_Gated_bin_' num2str(i)],'Compressed',true);
    Allims(:,:,:,i) = Im;
end
%%
testslice = 64;
viewslice = squeeze(Allims(:,testslice,:,:));

figure;imagesc(abs(viewslice(:,:,1)));
roi1 = roipoly;
roi2 = roipoly;
msk = roi1+roi2;

meansig = zeros(1,nbins);
for i = 1:nbins
    tmpim = viewslice(:,:,i);
    meansig(i) = mean(abs(tmpim(msk==1)));
end
figure;
plot(meansig)


% CardHigh = logical(resp_gate.*card_gate_high);
% CardLow = logical(resp_gate.*card_gate_low);
% 
% Imfinal = zeros(FinImSize,FinImSize,FinImSize);
% Imfinal2 = zeros(FinImSize,FinImSize,FinImSize);
% tmpDCF = Recon.get_DCF(mytraj(:,:,CardHigh),FinImSize,10);
% tmpDCF2 = Recon.get_DCF(mytraj(:,:,CardLow),FinImSize,10);
% for j = 1:size(fid,3)
%     tmp = Recon.mem_eff_recon(FinImSize,myfid(:,CardHigh,j),mytraj(:,:,CardHigh),tmpDCF,j,size(myfid,3),[0 0 0],false);
%     tmp2 = Recon.mem_eff_recon(FinImSize,myfid(:,CardLow,j),mytraj(:,:,CardLow),tmpDCF2,j,size(myfid,3),[0 0 0],false);
%     Imfinal = Imfinal + tmp.^2;
%     Imfinal2 = Imfinal2 + tmp2.^2;
% end
% card_high = sqrt(Imfinal);
% card_low = sqrt(Imfinal2);
