function output_txt = wiggledisplayfunc(empt,event_obj,fighandle1,fighandle2,xdata,ydata,HR)

% Display an observation's Y-data and label for a data tip
% empt          Not sure what this does... it's not used, but it's in the
%              original help, so I'm going to leave it in.
% fighandle1   figure showing T2*
% fighandle2   figure to which I want to plot data
% event_obj    Handle to event object
% xdata        x data
% ydata        4-D matrix of log of image data
% output_txt   Datatip text (string or string cell array)
% This datatip function spits out T2* as well as fitting information and
% shows the data from which T2* is calculated from the selected datapoint.


%Get location of datatip
pos = get(event_obj,'Position');

%Get the dimension and slice from imslice figure
h= findobj(fighandle1, 'type', 'uicontrol');
%Get the dimension that is selected
Dim = round(h(2).Value);
%Get the slice that is selected
Slice = round(h(3).Value);

%Correctly set x, y, and z based on the dimension specified in imslice
if Dim == 1
    x = Slice;
    y = pos(2);
    z = pos(1);
elseif Dim == 2
    x = pos(1);
    y = Slice;
    z = pos(2);
elseif Dim == 3
    x = pos(2);
    y = pos(1);
    z = Slice;
end

coords = [x y z];
x1 = xdata;
if isrow(x1)
    x1 = x1';
end
y1 = squeeze(ydata(coords(1),coords(2),coords(3),:));
if isrow(y1)
    y1 = y1';
end

xfit = linspace(0,max(xdata));

%Edit start point of FitOpts
%FitOpts.StartPoint(1) = y1(1);

FitAns = Tools.sineFit(x1,y1,HR/60,0);

%FitAns = fit(x1,y1,FitFunc,fitopts);
yfit = FitAns(1)+FitAns(2)*sin(2*pi*FitAns(3)*xfit+FitAns(4));%FitAns(xfit);

%Plot data and fit to figure
figure(fighandle2)
plot(x1,y1,'-ok')
hold on
plot(xfit,yfit,'LineWidth',2);
hold off
ylabel('Image Intensity');
xlabel('time (ms)');

% FitFields = fieldnames(FitAns);
% 
% FitResTxt = [];
% for i = 1:length(FitFields)
%     FitResults = FitAns.(FitFields{i});
%     FitResTxt = [FitResTxt FitFields{i} '=' num2str(FitResults) ' '];
% end

FitResTxt = ['Off = ' num2str(FitAns(1)), ', Amp = ' num2str(abs(FitAns(2)/FitAns(1))) ', Freq = ' num2str(FitAns(3)) ', Phase = ' num2str(FitAns(4))];


output_txt = {num2str(pos),['(x,y,z) = (',num2str(x),',',num2str(y),',',num2str(z),')'],FitResTxt};
