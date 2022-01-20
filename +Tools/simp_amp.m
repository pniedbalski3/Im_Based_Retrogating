function [Amp,Amp1] = simp_amp(Im_Stack,Mask)

%% Simplistic function to get the oscilation amplitude
%Find maxima and take average
%Find minima and take average
%average(max)-average(min)

Amp = zeros(size(Im_Stack,1),size(Im_Stack,2),size(Im_Stack,3));
Amp1 = zeros(size(Im_Stack,1),size(Im_Stack,2),size(Im_Stack,3));
for i = 1:size(Im_Stack,1)
    for j = 1:size(Im_Stack,2)
        for k = 1:size(Im_Stack,3)
            if Mask(i,j,k)
                Data = squeeze(Im_Stack(i,j,k,:));
                [maxima,maxlocs] =findpeaks(Data);
                [minima,minlocs] =findpeaks(-Data);
                %Report as % difference
                Amp(i,j,k) = ((mean(maxima) - mean(-minima))/mean(Data))*100; 
                Amp1(i,j,k) = ((mean(maxima) - mean(-minima))); 
            end
        end
    end
end

                