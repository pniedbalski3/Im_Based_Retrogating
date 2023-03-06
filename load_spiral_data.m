function [fid,traj,ImSize,FOV,mypath,twix_obj] = load_spiral_data()

[file,mypath] = uigetfile('.dat','Select FLORET UTE File');

%% Load Data
disp('Reading Data')
twix_obj = mapVBVD(file,'ignoreSeg');

if length(twix_obj)>1
    fid = double(twix_obj{2}.image.unsorted);
    twix_obj = twix_obj{2};
else
    fid = double(twix_obj.image().unsorted);
end
disp('Reading Data Complete');
%% Reshape Data as appropriate
disp('Reshaping');
fid = permute(fid,[1,3,2]);
disp('Reshaping Complete');
%% Get Trajectories
disp('Generating Trajectories');
traj = Tools.seek_spiral_traj(twix_obj);
disp('Generating Trajectories Complete');
%% Fix Trajectory delay
disp('Fix Trajectory Delay');
t_delay = 2.25;
Dw = twix_obj.hdr.MeasYaps.sRXSPEC.alDwellTime{1,1}/1000;
traj = traj_delay_correction(traj,Dw,Dw*t_delay);
disp('Fix Trajectory Delay Complete');
%% Get Image Size
ImSize = twix_obj.hdr.Config.BaseResolution;
FOV = twix_obj.hdr.Config.RoFOV;