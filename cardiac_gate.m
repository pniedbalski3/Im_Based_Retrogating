function [highind,lowind] = cardiac_gate(myfid,TR)

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

[~,mymax] = findpeaks(bp_fid);
[~,mymin] = findpeaks(-bp_fid);

%Get 1/10th of cardiac cycle = ~30 pts
maxind = [];
minind = [];

for i = 1:length(mymax)
    tmp = (mymax(i)-14):(mymax(i)+15);
    maxind = [maxind tmp];
end
for i = 1:length(mymin)
    tmp = (mymin(i)-14):(mymin(i)+15);
    minind = [minind tmp];
end

maxind(maxind<2000) = [];
minind(minind<2000) = [];


maxind(maxind>NPro-2000) = [];
minind(minind>NPro-2000) = [];

nexttile
plot(bp_fid,'k');hold on
plot(maxind,bp_fid(maxind),'r.');
plot(minind,bp_fid(minind),'b.');

tmp = zeros(1,length(myfid));

highind = tmp;
highind(maxind) = 1;
lowind = tmp;
lowind(minind) = 1;
