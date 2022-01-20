function [ADCMap,Mean,STD] = ADCFitting_3D(ImageIn,BNMask,bvals)
%function [ADCMapExpFit,ADCMapLinFit,MeanExp,STDExp,MeanLin,STDLin] = ADCFitting(ImageIn,BNMask,bvals)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function to fit ADC images
%Pass a stack of 3D Images and a vector containing the associated b-vals
%
%Author: Peter Niedbalski
%Work Email: pniedbalski@kumc.edu
%Personal Email: pniedbalski3@gmail.com
%Website:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Use Fitting based on equation Im(b) = Im(a)*exp(-b*ADC)

ImageIn = abs(ImageIn);

fo = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0 0],...
               'Upper',[Inf Inf],...
               'StartPoint',[.01 .01]);
fiteq =fittype('exp(-b*x)+c','options',fo);

ADCMap = zeros(size(ImageIn,1),size(ImageIn,2),size(ImageIn,3));
ADCMapLinFit = zeros(size(ImageIn,1),size(ImageIn,2),size(ImageIn,3));

if length(bvals)>2
    bvals = bvals';

    for i = 1:size(ImageIn,1)
        for j = 1:size(ImageIn,2)
            for k = 1:size(ImageIn,3)
                YData = squeeze(ImageIn(i,j,k,:)/ImageIn(i,j,k,1));
                if(BNMask(i,j,k) == 1)
                    ExpFit = fit(bvals,YData,fiteq,fo);
                    ADCMap(i,j,k) = ExpFit.b;
                     linFit = polyfit(-1*bvals,log(YData),1);
                     ADCMapLinFit(i,j,k) = linFit(2);
                end
            end
        end
    end
else
    ADCMap = -1/bvals(2)*log(ImageIn(:,:,:,2)./ImageIn(:,:,:,1)).*BNMask;
end

ADCMap(ADCMap > 0.2) = 0;
ADCMap(ADCMap < 0) = 0;
AveExp = ADCMap;
AveExp(BNMask == 0) = [];
AveExp(AveExp<0) = []; %ADC can't be zero, so kill any points that are negative
Mean = mean(AveExp);
STD = std(AveExp);

