function [Avg_Amp,Std_Amp] = bin_wiggles(disfid,gasfid,RBC2Bar,TR,write_path)

%Function to bin wiggles. Pass dissolved and gas fids, alongside
%RBC/Barrier, TR, and a path to write a .mat file to. At the end of the
%function, it will write the high-binned and low-binned indices, and total projections to a .mat file
%for future loading (so that this doesn't have to be repeated every time).
%This also makes it more simple to do manual binning, because you can
%create your binning manually, then save the file, which the subsequent
%analysis will automatically use.

purp = [168 96 168]/255;
[Dis_Shift,~] = Tools.SinglePointDixon_FID(disfid,-RBC2Bar,gasfid);

%We'll need to check that the total number of projections used for binning
%matches what is used.
Tot_Proj = size(disfid,2);

Gask0_1 = abs(gasfid(1,:));

%If we scale the dissolved k0 by the gas k0, that pretty well detrends data
Dis_Shift = Dis_Shift./Gask0_1;

RBC_Fid = real(Dis_Shift);
Bar_Fid = imag(Dis_Shift);

if(mean(RBC_Fid(1,:))<0)
    RBC_Fid = -RBC_Fid;
end
if(mean(Bar_Fid(1,:))<0)
    Bar_Fid = -Bar_Fid;
end

NProj = length(RBC_Fid);
Time_Axis = (0:(NProj-1))*TR/1000;
SampF = 1/(TR/1000);
Freq_Axis = SampF*((-(NProj/2)+1):(NProj/2))/NProj;

%Create a smoothing window
Sm_window = floor(1/(TR/1000)/5);

%Isolate k0
RBC2Gas = RBC_Fid(1,:);
Bar2Gas = Bar_Fid(1,:);

%% I don't think there is any need to fit any data 
%Barrier and Gas both will have similar relaxation behavior to RBC, so just
%scale the data by one of those - For now, we'll stick with gas

%Collier et al. used a bandpass filter with a range of [0.5 2.5] Hz. This
%seems to work pretty well and does minimal violence to the underlying
%wiggles. I think this is better than the combined smoothing/lowpass filter
%I used in the original paper
BP_RBC2Gas = bandpass(RBC2Gas,[0.7 2.5],SampF);
%BandPass filter puts this centered around 0. Let's make this centered
%around the mean RBC/Gas signal
BP_RBC2Gas = BP_RBC2Gas+mean(RBC2Gas);

%% While we're at it, let's get the heart rate of the subject from the wiggles
FFT_RBC = abs(fftshift(fft(BP_RBC2Gas)));
Use_Freq = Freq_Axis(Freq_Axis>0.5 & Freq_Axis<2.5);
[HR_Peak,HR_ind] = max(FFT_RBC(Freq_Axis>0.5 & Freq_Axis<2.5));
HFreq = Use_Freq(HR_ind);
HR = HFreq*60;

%% Get Global Oscillation Amplitude
%Now, similar to the method of Collier et al. Find the minima and maxima.
%We'll assume that the maximum heartrate we'll see is 120 BPM - This means
%maxima should be separated by at least 0.5 s
%We'll also only use the first 6 s of data to avoid possible breath-hold or
%low-SNR issues
Less6ind = find(Time_Axis<7 & Time_Axis>2);
Time6 = Time_Axis(Less6ind);
Data6 = BP_RBC2Gas(Less6ind);

[maxima,maxlocs] = findpeaks(Data6,Time6,'MinPeakDistance',0.5);
[minima,minlocs] = findpeaks(-Data6,Time6,'MinPeakDistance',0.5);
%Perhaps overly simplistic, but let's just get maxima - minima as the mean
%peak separation - make sure vectors are same length
if length(maxima)>length(minima)
    maxima = maxima(1:length(minima));
    maxlocs = maxlocs(1:length(minima));
elseif length(minima)>length(maxima)
    minima = minima(1:length(maxima));
    minlocs = minlocs(1:length(maxima));
end
%Difference and mean - express as a percentage of mean RBC (minima were
%found with a -, so use + to get a difference
%Avg_Amp = mean((maxima+minima)/mean(RBC2Gas)*100);
Avg_Amp = (mean(maxima)+mean(minima))/mean(RBC2Gas)*100;
Std_Amp = std((maxima+minima)/mean(RBC2Gas)*100);
%% Binning
% Before, I stretched data and picked peaks that way. Retrospective gating
% methods may work, but let's try something that may be a little more
% consistent
%Find peaks again - this time for all of the data.
[maximaAll,maxlocsAll] = findpeaks(BP_RBC2Gas,Time_Axis,'MinPeakDistance',0.5);
[minimaAll,minlocsAll] = findpeaks(-BP_RBC2Gas,Time_Axis,'MinPeakDistance',0.5);
%Now, based on TR and HR, decide how many points to use around each peak
%Previously, I used 20% (using a sloppy binning method), so let's use 20%
%again here
HR_Bit = 0.2*(60/HR); %time required (in s) for 20% of one heart beat
NPts = HR_Bit*1000/TR;
%Force NPts to be even to make this easy
if mod(floor(NPts),2) ==0
    NPts = floor(NPts);
else
    NPts = ceil(NPts);
end
%Now - we've found maxima and we know how many points we need. Just need to
%get those points now 
%Probably better ways to do this, but this should work
High_Indx = [];

for i = 1:length(maxlocsAll)
    idx(i) = find(Time_Axis==maxlocsAll(i));
    tmp_indx = (idx(i)-NPts/2):(idx(i)+NPts/2);
    High_Indx = [High_Indx tmp_indx];
end
h_idx = idx;

Low_Indx = [];
for i = 1:length(minlocsAll)
    idx(i) = find(Time_Axis==minlocsAll(i));
    tmp_indx = (idx(i)-NPts/2):(idx(i)+NPts/2);
    Low_Indx = [Low_Indx tmp_indx];
end
l_idx = idx;

    %We might occasionally go beyond the indices allowed... fix that here
    High_Indx(High_Indx<1) = [];
    Low_Indx(Low_Indx<1) = [];
    High_Indx(High_Indx>Tot_Proj) = [];
    Low_Indx(Low_Indx>Tot_Proj) = [];
%Let's force these to have the same number of points - pretty easy to do in
    %this case - pull points from the end since those will be low SNR anyway
    if length(High_Indx)>length(Low_Indx)
        High_Indx = High_Indx(1:length(Low_Indx));
    elseif length(High_Indx)<length(Low_Indx)
        Low_Indx = Low_Indx(1:length(High_Indx));
    end
%We also need the "High" and "Low" RBC/Barrier ratio for separation of
    %images later
    High_RBC2Bar = mean(RBC2Gas(High_Indx))/mean(Bar2Gas(High_Indx));
    Low_RBC2Bar = mean(RBC2Gas(Low_Indx))/mean(Bar2Gas(Low_Indx));
%Not using this right now, so don't want it to trip things up.
try
    
    %Bin some of the intermediate points - Will this give us more/better
    %information?
    Int1_Indx = [];
    Int2_Indx = [];

    Int1a_Indx = [];
    Int2a_Indx = [];
    Int3a_Indx = [];
    Int4a_Indx = [];
    for i = 1:min([length(minlocsAll) length(maxlocsAll)])
        min1 = find(Time_Axis == minlocsAll(i));
        max1 = find(Time_Axis == maxlocsAll(i));
        %Rather than using x axis points, use y axis - specifically points
        %within the middle 20%
        idx = ceil(mean([min1 max1]));
        tmp_indx = (idx-NPts/2):(idx+NPts/2);
        Int1_Indx = [Int1_Indx tmp_indx];
        if length(minlocsAll) > length(maxlocsAll)
            next1 = find(Time_Axis == minlocsAll(i+1));
            next2 = find(Time_Axis == maxlocsAll(i));
        else
            try
            next1 =find(Time_Axis == maxlocsAll(i+1));
            next2 = find(Time_Axis == minlocsAll(i));
            catch
            end
        end
        try
            idx = ceil(mean([next1 next2]));
            tmp_indx = (idx-NPts/2):(idx+NPts/2);
            Int2_Indx = [Int2_Indx tmp_indx];
        catch
        end
        soft = round(linspace(min1,max1,4));
        idx = soft(2);
        tmp_indx = (idx-NPts/2):(idx+NPts/2);
        Int1a_Indx = [Int1a_Indx tmp_indx];
        idx = soft(3);
        tmp_indx = (idx-NPts/2):(idx+NPts/2);
        Int2a_Indx = [Int2a_Indx tmp_indx];
        try
            soft = round(linspace(next1,next2,4));
            idx = soft(2);
            tmp_indx = (idx-NPts/2):(idx+NPts/2);
            Int3a_Indx = [Int3a_Indx tmp_indx];
            idx = soft(3);
            tmp_indx = (idx-NPts/2):(idx+NPts/2);
            Int4a_Indx = [Int4a_Indx tmp_indx];
        catch
        end
    end

    fourbin.bin1 = High_Indx;
    fourbin.bin2 = Int1_Indx;
    fourbin.bin3 = Low_Indx;
    fourbin.bin4 = Int2_Indx;

    %Make sure these get put in the right order
    if maxlocsAll(1)<minlocsAll(1)
        sixbin.bin1 = High_Indx;
        sixbin.bin2 = Int1a_Indx;
        sixbin.bin3 = Int2a_Indx;
        sixbin.bin4 = Low_Indx;
        sixbin.bin5 = Int3a_Indx;
        sixbin.bin6 = Int4a_Indx;
    else
        sixbin.bin1 = Low_Indx;
        sixbin.bin2 = Int1a_Indx;
        sixbin.bin3 = Int2a_Indx;
        sixbin.bin4 = High_Indx;
        sixbin.bin5 = Int3a_Indx;
        sixbin.bin6 = Int4a_Indx;
    end

    fourbin.rbc2bar(1) = mean(RBC2Gas(fourbin.bin1))/mean(Bar2Gas(fourbin.bin1));
    fourbin.rbc2bar(2) = mean(RBC2Gas(fourbin.bin2))/mean(Bar2Gas(fourbin.bin2));
    fourbin.rbc2bar(3) = mean(RBC2Gas(fourbin.bin3))/mean(Bar2Gas(fourbin.bin3));
    fourbin.rbc2bar(4) = mean(RBC2Gas(fourbin.bin4))/mean(Bar2Gas(fourbin.bin4));

    sixbin.rbc2bar(1) = mean(RBC2Gas(sixbin.bin1))/mean(Bar2Gas(sixbin.bin1));
    sixbin.rbc2bar(2) = mean(RBC2Gas(sixbin.bin2))/mean(Bar2Gas(sixbin.bin2));
    sixbin.rbc2bar(3) = mean(RBC2Gas(sixbin.bin3))/mean(Bar2Gas(sixbin.bin3));
    sixbin.rbc2bar(4) = mean(RBC2Gas(sixbin.bin4))/mean(Bar2Gas(sixbin.bin4));
    sixbin.rbc2bar(5) = mean(RBC2Gas(sixbin.bin5))/mean(Bar2Gas(sixbin.bin5));
    sixbin.rbc2bar(6) = mean(RBC2Gas(sixbin.bin6))/mean(Bar2Gas(sixbin.bin6));
catch
    fourbin = 1;
    sixbin = 1;
end

[Alt_Bin,Alt_Min2Min] = Tools.tenbin_wiggles(BP_RBC2Gas,Time_Axis,h_idx,l_idx);
for i = 1:size(Alt_Bin,1)
    Alt_RBC2Bar(i) = mean(RBC2Gas(Alt_Bin(i,:)))/mean(Bar2Gas(Alt_Bin(i,:)));
end
    
%% Summary Figure for Global Oscillations and Binning
WiggleRawFig = figure('Name','Global Oscillation Analysis');
set(WiggleRawFig,'color','white','Units','inches','Position',[1 1 12 8])
tiledlayout(3,3)
%Display Gas
nexttile
plot(Time_Axis,Gask0_1,'.g')
hold on
plot(Time_Axis,abs(disfid(1,:)),'.','Color',purp)
legend('Gas','Dissolved')
title('Gas k0')
%Display Barrier
nexttile
plot(Time_Axis,Bar2Gas,'.b')
title('Barrier k0')
%Display RBC
nexttile
plot(Time_Axis,RBC2Gas,'.r');
title('RBC k0')
%Display RBC/Gas
nexttile
plot(Time_Axis,RBC2Gas,'.r')
hold on
plot(Time_Axis,BP_RBC2Gas,'k','LineWidth',2)
legend('RBC/Gas','Bandpass')
title('RBC/Gas')
%Display Fourier Transform - only show frequences in bandpass range
nexttile
plot(Freq_Axis(Freq_Axis >0.5 & Freq_Axis<2.5),FFT_RBC(Freq_Axis >0.5 & Freq_Axis<2.5),'k')
text(HFreq,HR_Peak,['\leftarrow HR: ' num2str(HR,'%.1f') ' BPM']);
%annotation('textbox',[0.1 0.8 0.8 0.2],'String',['HR: ' num2str(HR,'%.1f') ' BPM'],'FitBoxToText','on');
title('FFT of Bandpass RBC')
%Display Global Oscillation Amplitudes
nexttile
plot(Time6,Data6,'r','LineWidth',2);
hold on
%plot maxima as circles
plot(maxlocs,maxima,'ko','MarkerFaceColor','k')
%plot minima as squares
plot(minlocs,-minima,'ms','MarkerFaceColor','m')
text(1,1.1*max(maxima),['Mean Osc. Ampl: ' num2str(Avg_Amp,'%.1f') '% \pm ' num2str(Std_Amp,'%.1f') '%'  ])
ylim([0.9*min(-minima) 1.15*max(maxima)]);
%plot all maxima and minima on bandpass filtered data
nexttile
plot(Time_Axis,BP_RBC2Gas,'r','LineWidth',2)
hold on
plot(maxlocsAll,maximaAll,'go','MarkerFaceColor','g')
plot(minlocsAll,-minimaAll,'ms','MarkerFaceColor','m')
plot(Time_Axis(High_Indx),BP_RBC2Gas(High_Indx),'.g','MarkerSize',12)
plot(Time_Axis(Low_Indx),BP_RBC2Gas(Low_Indx),'.m','MarkerSize',12)
%plot(Time_Axis(Int1_Indx),BP_RBC2Gas(Int1_Indx),'.c','MarkerSize',12)
%plot(Time_Axis(Int2_Indx),BP_RBC2Gas(Int2_Indx),'.b','MarkerSize',12)
%plot all selected points on raw RBC data
nexttile
plot(Time_Axis,RBC2Gas,'.r','MarkerSize',12)
hold on
plot(Time_Axis(High_Indx),RBC2Gas(High_Indx),'.g','MarkerSize',12)
plot(Time_Axis(Low_Indx),RBC2Gas(Low_Indx),'.m','MarkerSize',12)
%plot(Time_Axis(Int1_Indx),RBC2Gas(Int1_Indx),'.c','MarkerSize',12)
%plot(Time_Axis(Int2_Indx),RBC2Gas(Int2_Indx),'.b','MarkerSize',12)
%we have an extra tile - let's plot all the soft gated sillness as well
% nexttile
% plot(Time_Axis,BP_RBC2Gas,'.r','MarkerSize',12)
% hold on
% plot(Time_Axis(High_Indx),BP_RBC2Gas(High_Indx),'.g','MarkerSize',12)
% plot(Time_Axis(Low_Indx),BP_RBC2Gas(Low_Indx),'.m','MarkerSize',12)
% plot(Time_Axis(Int1a_Indx),BP_RBC2Gas(Int1a_Indx),'.c','MarkerSize',12)
% plot(Time_Axis(Int2a_Indx),BP_RBC2Gas(Int2a_Indx),'.b','MarkerSize',12)
% plot(Time_Axis(Int3a_Indx),BP_RBC2Gas(Int3a_Indx),'.k','MarkerSize',12)
% plot(Time_Axis(Int4a_Indx),BP_RBC2Gas(Int4a_Indx),'.','MarkerSize',12,'Color',purp)

save(fullfile(write_path,'Wiggle_Binning.mat'),'Low_Indx','High_Indx','Tot_Proj','High_RBC2Bar','Low_RBC2Bar','Avg_Amp','Std_Amp','sixbin','fourbin','HR','Alt_Bin','Alt_RBC2Bar','Alt_Min2Min')

saveas(WiggleRawFig,fullfile(write_path,'Global_Wiggle_Summary.png'))
%% This is binning based on retrospective gating
%This works fairly well, but requires much more user input, and there are
%still some issues with it. I think, for now, I like the bandpass filter
%method.

% happy_gate = 0;
% 
% %Set initial values for derivative limits
% hLL_val = -30;
% hUL_val = 30;
% 
% lLL_val = -30;
% lUL_val = 30;
% 
% high_thres = 0.5;
% low_thres = 2;
% 
% while happy_gate == 0
%     h = figure('Name','Wiggle Binning');
%     set(h,'Units','Normalized','Position',[.05 .05 .9 .8],'Color','w')
%     
%     subplot(2,3,1);
%     plot(Time_Axis,RBC2Gas, '.r');
%     title('RBC/Gas')
%     % Display smoothed data 
%     subplot(2,3,2);
%     smRBC = smooth(RBC2Gas,Sm_window);
%     %plot(Time_Axis,BP_RBC2Gas,'.k');
%     plot(Time_Axis,smRBC,'.k')
%     title('Smoothed RBC/Gas');
%     
%     %Get first derivative
%     RBC2Gas_Der1 = diff(smRBC,1);
%     %Leaves off last point - need to fill that here
%     RBC2Gas_Der1(Tot_Proj) = RBC2Gas_Der1(Tot_Proj-1);
%     smRBC2Gas_Der1 = smooth(RBC2Gas_Der1,Sm_window)*1e4;
%     subplot(2,3,3)
%     rectangle('position', [0 lLL_val max(Time_Axis) -lLL_val+lUL_val], 'facecolor',...
%            'b', 'edgecolor', 'w'); hold on;
% 
%     rectangle('position', [0 hLL_val max(Time_Axis) -hLL_val+hUL_val], 'facecolor',...
%            'r', 'edgecolor', 'w'); hold on;
%             
%     plot(Time_Axis,smRBC2Gas_Der1,'k')
%     
%     line([0 max(Time_Axis)], [lLL_val lLL_val], 'color', 'b',...
%         'linewidth', 4, 'linestyle', '--');
% 
%     line([0 max(Time_Axis)], [lUL_val lUL_val], 'color', 'b',...
%         'linewidth', 4, 'linestyle', '--'); hold off; 
%     
%     line([0 max(Time_Axis)], [hLL_val hLL_val], 'color', 'r',...
%         'linewidth', 4, 'linestyle', '--');
% 
%     line([0 max(Time_Axis)], [hUL_val hUL_val], 'color', 'r',...
%         'linewidth', 4, 'linestyle', '--'); hold off; 
%     title('First Derivative')
%     ylim([1.1*min(smRBC2Gas_Der1(100:(end-100))) 1.1*max(smRBC2Gas_Der1(100:(end-100)))]);
%     
%     %Second Derivative
%     RBC2Gas_Der2 = diff(smRBC2Gas_Der1);
%     RBC2Gas_Der2(Tot_Proj) = RBC2Gas_Der2(Tot_Proj-1);
%     smRBC2Gas_Der2 = smooth(RBC2Gas_Der2,Sm_window);
%     
%     subplot(2,3,4);
%     rectangle('position', [0 -100 max(Time_Axis)....
%         100 + high_thres],...
%         'facecolor', [243/255, 151/255, 151/255], 'edgecolor', 'w'); hold on;
%     rectangle('position', [0 low_thres max(Time_Axis) max(smRBC2Gas_Der2)/2],...
%         'facecolor', [171/255, 209/255, 255/255], 'edgecolor', 'w');
%     plot(Time_Axis,smRBC2Gas_Der2, '-k', 'linewidth', 3) 
% 
% 
%     line([0 max(Time_Axis)], [high_thres high_thres], 'color', 'r',...
%         'linewidth', 4, 'linestyle', '--');
% 
%     line([0 max(Time_Axis)], [low_thres low_thres], 'color', 'b',...
%         'linewidth', 4, 'linestyle', '--');
%     title('Second Derivative')
%     ylim([1.1*min(smRBC2Gas_Der2(100:(end-100))) 1.1*max(smRBC2Gas_Der2(100:(end-100)))])
%     
%     StartGateSmoothHigh  = smRBC;
%     StartGateSmoothLow = smRBC;
% 
%     % This helps to make good looking plots
%     StartGateSmoothHigh(smRBC2Gas_Der1 > hLL_val &...
%     smRBC2Gas_Der1 < hUL_val &...
%     smRBC2Gas_Der2 < high_thres) = max(smRBC(:));
%     % This variable contains locations of k0 points at "high" oscillation
%     SelectVectorHigh = smRBC2Gas_Der1 > hLL_val &...
%         smRBC2Gas_Der1 < hUL_val &...
%         smRBC2Gas_Der2 < high_thres;
%     % This helps to make good looking plots
%     StartGateSmoothLow(smRBC2Gas_Der1 > lLL_val &...
%     smRBC2Gas_Der1 < lUL_val &...
%     smRBC2Gas_Der2 > low_thres) =  min(smRBC(:));
% 
% % This variable contains locations of k0 points at "Low" oscillation
%     SelectVectorLow = smRBC2Gas_Der1 > lLL_val &...
%         smRBC2Gas_Der1 < lUL_val &...
%         smRBC2Gas_Der2 > low_thres;
%     disp('Finished Gating at End-Inspiration');
%     
%     % This is strictly for plotting purposes
%     StartGateSmoothHigh(StartGateSmoothHigh   ~= max(StartGateSmoothHigh)) = nan;
%     StartGateSmoothLow(StartGateSmoothLow ~= min(StartGateSmoothLow)) = nan;
% 
%     SelectedProjHigh  = SelectVectorHigh.*smRBC; SelectedProjHigh(SelectedProjHigh == 0) = nan;
%     SelectedProjLow = SelectVectorLow.*smRBC; SelectedProjLow(SelectedProjLow == 0) = nan;
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     subplot(2,3,5);
% 
%     plot(StartGateSmoothHigh, '-k', 'linewidth', 3); hold on;
%     plot(smRBC, '-', 'color', [105/256 105/256 105/256], 'linewidth', 3)
%     plot(SelectedProjHigh, '-r', 'linewidth', 3);
%     plot(StartGateSmoothLow, '-k', 'linewidth', 3);
%     plot(SelectedProjLow, '-b', 'linewidth', 3); hold off;
%     title('Selected Projections')
%     subplot(2,3,6)
%     plot(Time_Axis,RBC2Gas,'.k')
%     hold on
%     plot(Time_Axis(SelectVectorHigh),RBC2Gas(SelectVectorHigh),'.g')
%     plot(Time_Axis(SelectVectorLow),RBC2Gas(SelectVectorLow),'.m')
%     
%     answer = questdlg('Is this Binning okay?','Gating Confirmation','Yes','No','Skip','Yes');       
%     if strcmp(answer,'No') 
%         prompt = {['Low Lower Threshold (Was ' num2str(lLL_val) ')'],['Low Upper Threshold (Was ' num2str(lUL_val) ')'],['High Lower Threshold (Was ' num2str(hLL_val) ')'],['High Upper Threshold (Was ' num2str(hUL_val) ')'],['High Threshold (Was ' num2str(high_thres) ')'],['Low Threshold (Was ' num2str(low_thres) ')']};
%         dlgtitle = 'New Binning Parameters';
%         dims = [1 50];
%         definput = {num2str(lLL_val),num2str(lUL_val),num2str(hLL_val),num2str(hUL_val),num2str(high_thres),num2str(low_thres)};
%         user_input = inputdlg(prompt,dlgtitle,dims,definput);
%         lLL_val = str2num(user_input{1});
%         lUL_val = str2num(user_input{2});
%         hLL_val = str2num(user_input{3});
%         hUL_val = str2num(user_input{4});
%         high_thres = str2num(user_input{5});
%         low_thres = str2num(user_input{6});
%     elseif strcmp(answer,'Yes')
%         happy_gate = 1;
%        % saveas(h,fullfile(write_path,'Wiggle_Binning_Summary.png'))
%     else
%         happy_gate = 1;
%         SelectVectorHigh = ones(size(k0));
%         SelectVectorLow = NaN;
%     end
%     close
% end
% 
% low_indx = SelectVectorLow;
% high_indx = SelectVectorHigh;

%save(fullfile(write_path,'Wiggle_Binning.mat'),'low_indx','high_indx','Tot_Proj')