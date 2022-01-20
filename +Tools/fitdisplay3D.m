function fitdisplay3D(RawDataX,RawDataY,FitVal,FitFunc,fitopts)

if nargin < 5
    fitopts = fitoptions('Method','NonlinearLeastSquares');
end

myFig = imslice(FitVal,'Fitting Result');

NCmap = parula;
NCmap(1,:) = [0 0 0];
colormap(myFig,NCmap)

rawfig = figure('Name','Data from which Fit is Derived');

dcm_obj = datacursormode(myFig);
set(dcm_obj,'UpdateFcn',{@Tools.fittingdisplayfunc,myFig,rawfig,RawDataX,RawDataY,FitFunc,fitopts})


