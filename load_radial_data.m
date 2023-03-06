function [fid,traj,ImSize,FOV,mypath,twix_obj] = load_radial_data()

[myfile,mypath] = uigetfile('*.dat','Select Radial Twix File');


%% Read Twix
twix_obj = mapVBVD_pjn(fullfile(mypath,myfile),'ignoreSeg');

if length(twix_obj)>1
    fid = double(twix_obj{2}.image());
    twix_obj = twix_obj{2};
else
    fid = double(twix_obj.image());
end

fid = permute(fid,[1,3,2]);

%% Hardcoded Parameters for my Radial Setup
ImSize = 400; %Desired Recon size. Specified 268 for 1.5 mm isotropic res, but went way farther out in k-space than that.
FOV = 400;
ADC_Dur = 903; %ADC in us
RUT = 140; %Ramp Up Time in us
Resolution = [FOV/ImSize FOV/ImSize FOV/ImSize]/1000;
%gamma = 11.777;
gamma = 42.6;
MaxGrad = 13.980741193; %Max Gradient in mT/m from simulation ERR file
t_delay = 4; % 6 appears to be the magic number for getting really nice images.
%This number may have changed when I changed the echo time slightly.
TR = 3.5;

%% Create Trajectories
NPro = size(fid,2);
Pts = size(fid,1);
Dur = ADC_Dur-RUT;

Grad = linspace(0,1,RUT);
Grad((RUT):ADC_Dur) = 1;

Grad = Grad * MaxGrad;

ADC_Dur = ADC_Dur*1e-6;
RUT = RUT*1e-6;
Dw = ADC_Dur/Pts;

Arm_untimed = cumtrapz(Grad);
Grad_Time = 0:1e-6:(ADC_Dur-1e-6);

Time = 0:Dw:(ADC_Dur-Dw);
Arm = interp1(Grad_Time,Arm_untimed,Time); %Now we are in mT*s/m

Arm = Arm*gamma/1000;

phi1 = 0.46557123;
phi2 = 0.6823278;
gs = 1;
gr = 1;
gp = 1;

r = zeros(1,NPro);
p = zeros(1,NPro);
s = zeros(1,NPro);
%Rotation code from Jinbang's UTE sequence
for i = 0:(NPro-1)
    kz = (i*phi1-floor(i*phi1))*2-1;
    ts = kz*gs;
    alpha = (i*phi2-floor(i*phi2))*2*pi;
    tr = sqrt(1-kz*kz)*cos(alpha)*gr;
    tp = sqrt(1-kz*kz)*sin(alpha)*gp;
    r(i+1) = tr;
    p(i+1) = tp;
    s(i+1) = -ts;
end

trajx_rot = zeros(length(Arm),NPro);
trajy_rot = trajx_rot;
trajz_rot = trajx_rot;
for i = 1:NPro
    trajx_rot(:,i) = r(i)*Arm';
    trajy_rot(:,i) = p(i)*Arm';
    trajz_rot(:,i) = s(i)*Arm';
end

traj = cat(3,trajx_rot,trajy_rot,trajz_rot);
traj = permute(traj,[3,1,2]);

kFOV_desired = 1./(Resolution);
kMax_desired = kFOV_desired/2;
max_k = max(kMax_desired); %Here, we are in 1/m
traj = traj/max_k/2;

traj = traj_delay_correction(traj,Dw,Dw*t_delay);
