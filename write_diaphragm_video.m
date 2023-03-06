function write_diaphragm_video(All_Im)


G(size(All_Im,3)) = struct('cdata',[],'colormap',[]);

writerObj = VideoWriter('Diaphragm_Motion.avi','Motion JPEG AVI');
writerObj.Quality = 95;
writerObj.FrameRate = 5; % How many frames per second.
open(writerObj); 

My_Fig = figure('Name','Diaphragm_Motion_Figure');
set(My_Fig,'color','w','Units','inches','Position',[0.25 0.25 7 7])
for k = 1:size(All_Im,3)
    imagesc(squeeze(All_Im(:,:,k)));
    colormap(gray)
    axis square
    axis off
    G(k) = getframe;
    writeVideo(writerObj,G(k));
end

close(writerObj);