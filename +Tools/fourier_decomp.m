function [Amp,Phase,Freq,R2,BP_ImStack] = fourier_decomp(Im_Stack,Mask,Time_Res,HR)

nKey = size(Im_Stack,4);
t = (0:(nKey-1))*Time_Res;
SampF = 1/(Time_Res);
Freq_Axis = SampF*((-(nKey/2)+1):(nKey/2))/nKey;

Amp = zeros(size(Im_Stack,1),size(Im_Stack,2),size(Im_Stack,3));
Phase = zeros(size(Im_Stack,1),size(Im_Stack,2),size(Im_Stack,3));
Freq = zeros(size(Im_Stack,1),size(Im_Stack,2),size(Im_Stack,3));
R2 = zeros(size(Im_Stack,1),size(Im_Stack,2),size(Im_Stack,3));

BP_ImStack = zeros(size(Im_Stack));

fo = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0,0.5*2*pi,-pi],...
               'Upper',[Inf,1.66*2*pi,pi],...
               'StartPoint',[0.1, HR/60*2*pi,0]);

for i = 1:size(Im_Stack,1)
    for j = 1:size(Im_Stack,2)
        for k = 1:size(Im_Stack,3)
            if Mask(i,j,k)
                Data = squeeze(Im_Stack(i,j,k,:));
                %Found a function online that seems to perform pretty
                %well... And it doesn't require a bandpass filter. And it's
                %really fast.
                Result = Tools.sineFit(t,Data',HR/60,0);
                %Test_Data = Data - mean(Data);
%                 bpData = bandpass(Data,[0.5 1.6],SampF);
%                 BP_ImStack(i,j,k,:) = bpData;
%                 if ~iscolumn(bpData)
%                     bpData = bpData';
%                 end
%                 if ~iscolumn(t)
%                     t = t';
%                 end
%                 [fitObj,gof] = fit(t,bpData,'sin1',fo);
%                 
%                 %fitObj = Spectroscopy.NMR_TimeFit_v(Test_Data,t,1,1,1,1,0,0,length(t)); % first widths lorenzian, 2nd are gauss
%                 %fitObj = fitObj.fitTimeDomainSignal();
%                 Amp(i,j,k) = fitObj.a1;
%                 Phase(i,j,k) = fitObj.c1;
%                 Freq(i,j,k) = fitObj.b1/2/pi;
%                 R2(i,j,k) = gof.rsquare;
                %I think I want to do the amplitude/offset so that we get a
                %percentage:
                Amp(i,j,k) = abs(Result(2)/Result(1));
                Phase(i,j,k) = Result(4);
                Freq(i,j,k) = Result(3);
                R2(i,j,k) = Result(5);
            end
        end
    end
end
%Phase is just arbitrary anyway, so let's make the mean phase of the
%oscillations 0 to better visualize phase differences.
meanPhase = mean(Phase(Mask(:)));
Phase = Phase-meanPhase;
Phase(Phase<-pi) = Phase(Phase<-pi)+2*pi;
