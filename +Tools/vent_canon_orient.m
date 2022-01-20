function Image = vent_canon_orient(In_Image)

Image = permute(In_Image,[2,3,1]);
Image = flip(Image,3);
Image = flipud(Image);
Image = fliplr(Image);