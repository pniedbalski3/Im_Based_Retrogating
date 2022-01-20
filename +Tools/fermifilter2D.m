function fermi2 = fermifilter2D(GridSize, Wf)
% pull the gridded k-space dimensions to create the fermi filter
kx = GridSize(1); ky = GridSize(2); 
%GE default for Wf is 10
% hard code 1D BA - 07Aug2019 -- add
% 1D x    
Fermi_xx = Tools.fermifilter1D(kx,Wf);% 1D filter in x dim
Fermi_xx = repmat(Fermi_xx, ky,1);
Fermi_xx = rot90(reshape(Fermi_xx, [kx, ky]));
%figure; mesh(Fermi_xx)
% 1D y
Fermi_yy = Tools.fermifilter1D(ky,Wf);% 1D filter in x dim
Fermi_yy = repmat(Fermi_yy, kx,1);
Fermi_yy =(reshape(Fermi_yy, [ky, kx]));
%figure; mesh(Fermi_yy)

fermi2=Fermi_xx.*Fermi_yy; 

end