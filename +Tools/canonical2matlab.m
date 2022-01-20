function Im = canonical2matlab(In_Im)

%Image = permute(In_Image,[2,3,1])
Im = permute(In_Im,[3,1,2]);
Im = rot90(Im,2);