function finalbin = cardiac_gate(myfid,TR)

%FID and TR in ms
myfid = abs(myfid(1,:));
SampF = 1/(TR/1000);
NPro = size(myfid,2);

Time_Axis = (0:NPro-1)*TR/1000;
Freq_Axis = SampF*((-(NPro/2)+1):(NPro/2))/NPro;

Data_FFt = abs(fftshift(fft(myfid)));

f_display = (Freq_Axis > 0) & (Freq_Axis < 3);

figure('Name','Check_FFT')
tiledlayout(2,2)
nexttile
plot(Freq_Axis(f_display),Data_FFt(f_display));

bp_fid = bandpass(myfid,[1.2 1.4],SampF);

nexttile
plot(bp_fid);

%% Cut off first and last 1000 points
%bp_fid(end-2000) = [];
%bp_fid(1:2000) = [];

% [~,mymax] = findpeaks(bp_fid);
% [~,mymin] = findpeaks(-bp_fid);
%Let's do 5 bins
nbin = 5;

finalbin = zeros(size(bp_fid));
%Let's do 5 s at a time
step = floor(5/(TR/1000));
for i = 2001:step:(length(bp_fid)-2000)
    tmpdata = bp_fid(i:(i+step-1));
    binned = zeros(nbin,length(tmpdata));
    tmpd1 = diff(tmpdata);
    tmpd1(end+1) = tmpd1(end);
    %Let's do 5 bins
    tmpmax = max(tmpdata);
    tmpmin = min(tmpdata);
    mybinstep = (tmpmax-tmpmin)/5;
    for j = 1:nbin
        binned(j,tmpdata>tmpmin + mybinstep*(j-1) & tmpdata<tmpmin + mybinstep*(j)) = 1;
        if j > 1 && j < nbin
            binned(j,:) = binned(j,:).*sign(tmpd1);
        end
        binned(j,:) = binned(j,:)*j;
    end
    tmpbinfinal = sum(binned,1);
    finalbin(i:(i+step-1)) = tmpbinfinal;
    clear tmpd1;
end
%%
nexttile;
plot(Time_Axis, bp_fid);
hold on;
for i = 1:(nbin + nbin-2)
    if i > nbin
        plot(Time_Axis(finalbin==(-1*(i-nbin+1))),bp_fid(finalbin==(-1*(i-nbin+1))),'.');
    else
        plot(Time_Axis(finalbin==i), bp_fid(finalbin==i),'.')
    end
end

%Get 1/10th of cardiac cycle = ~30 pts
% maxind = [];
% minind = [];
% 
% for i = 1:length(mymax)
%     tmp = (mymax(i)-14):(mymax(i)+15);
%     maxind = [maxind tmp];
% end
% for i = 1:length(mymin)
%     tmp = (mymin(i)-14):(mymin(i)+15);
%     minind = [minind tmp];
% end
% 
% maxind(maxind<2000) = [];
% minind(minind<2000) = [];
% 
% 
% maxind(maxind>NPro-2000) = [];
% minind(minind>NPro-2000) = [];

% nexttile
% plot(bp_fid,'k');hold on
% plot(maxind,bp_fid(maxind),'r.');
% plot(minind,bp_fid(minind),'b.');
% 
% tmp = zeros(1,length(myfid));
% 
% highind = tmp;
% highind(maxind) = 1;
% lowind = tmp;
% lowind(minind) = 1;
