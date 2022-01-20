function [CorrectedDissKSpace] = gas_contam_remove(DissKSpace,GasKSpace,dwell_s,GasFreq_DissAcq,GasPhase_DissAcq,GasArea_DissAcq,GasFA,traj,BW_Scale)
%GasPhaseContaminationRemoval - Removes gas phase contamination in dissolved kspace
%  Takes gas phase k-space and modifies it using NMR fits and gas phase k0
%  to produce the expected gas phase contamination k-space data which is
%  then removed from the initial contaminated dissolved phase k-space
%
%  Inputs - DissKSpace - Dissolved kSpace data in [RO, Projections]
%           GasKSpace -  Gas kSpace data in [RO, Projections]
%           dwell_s - imaging dwell time in seconds
%           GasFreq_DissAcq - Gas frequency offset in dissolved spectra acquisition
%           GasPhase_DissAcq - Gas phase in dissolved spectra acquisition
%           GasArea_DissAcq - Gas area in dissolved spectra acquisition
%           GasFA - Gas phase excitation flip angle
%
%  Outpout - CorrectedDissKSpace
%
%   Author: Matthew Willmering
%   Work email: matthew.willmering@cchmc.org
%   Personal email: matt.willmering@gmail.com
%   Website: https://cpir.cchmc.org/
%   Revision history: 
%       December-2019   Finished working version
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Step 0: calculate parameters
time = dwell_s*(0:size(GasKSpace,1)-1); %time each point is acquired

%% Step 1: modulate Contamination (gas) to Dissolved Frequency - first order phase correction
Step1PhaseShift = time' * GasFreq_DissAcq * 2 * pi; %calculate phase accumulation during acquisition
Step1ContamKSpace = GasKSpace .* exp(1i*Step1PhaseShift);%apply varying phase as function of readout

%% Step 2: zero order phase shift of contamination estimation
Step2PhaseShift = GasPhase_DissAcq - rad2deg(mean(angle(GasKSpace(1,:,:)),'all'));
Step2ContamKSpace = Step1ContamKSpace * exp(1i*deg2rad(Step2PhaseShift));%apply same phase to all readout points in projection

%% Step 3: scale contamination 
%get k0 of gas image at last projection
GasK0Mean = movmean(abs(GasKSpace(1,:)),100);  %could use fit to reduce noise but would impart more error if exhaled or other changes
GasArea = GasK0Mean(1,end);
Step3Scale = GasArea_DissAcq/GasArea/sqrt(BW_Scale); %calculate scale of contamination to gas
Step3ContamKSpace = Step2ContamKSpace * Step3Scale / cosd(GasFA); %account for decay from gas pulse

%% Subtract contamination
CorrectedDissKSpace = DissKSpace - Step3ContamKSpace;

%Reconstruct Contamination to debug
% fid = reshape(Step3ContamKSpace,1,[])';
% 
% fid1= reshape(DissKSpace,1,[])';
% 
% Contam = Reconstruction.Dissolved_Phase_LowResRecon(64,fid,traj);
end

