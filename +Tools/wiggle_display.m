function wiggle_display(RawDataX,RawDataY,FitVal,HR)

myFig = imslice(FitVal,'Fitting Result');

NCmap = parula;
NCmap(1,:) = [0 0 0];
colormap(myFig,NCmap)

rawfig = figure('Name','Data from which Fit is Derived');

dcm_obj = datacursormode(myFig);
set(dcm_obj,'UpdateFcn',{@Tools.wiggledisplayfunc,myFig,rawfig,RawDataX,RawDataY,HR})


