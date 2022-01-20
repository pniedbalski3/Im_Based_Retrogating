function diaphragm_pos = THREE_diaphragm_motion(All_Ims,X,Y)

%Function to evaluate all 2D images to get a measure of diaphragm position.
%Ultimately, it would be cool to include some of the neat 2D correlation
%matrix stuff, but for now just do a simple line (specified by X and Y)
%over the diaphragm.

%Start by making sure the line is in the right place
figure('Name','Check Line')
imagesc(All_Ims(:,:,1));
colormap(gray)
axis square
axis off
hold on
plot(X,Y,'g','LineWidth',2);
hold off


Lines = zeros(max([X(2)-X(1),Y(2)-Y(1)]),size(All_Ims,3));
dLines = zeros((max([X(2)-X(1),Y(2)-Y(1)]))-1,size(All_Ims,3));
diaphragm_pos = zeros(1,size(All_Ims,3));
for i = 1:size(All_Ims,3)
    pts = Lines(:,i);
    if Y(1) == Y(2)
        Lines(:,i) = All_Ims(Y(1):Y(2),(X(1)+1):X(2),i);
    else
        Lines(:,i) = All_Ims((Y(1)+1):Y(2),X(1):X(2),i);
    end
    diffline = abs(diff(Lines(:,i)));
%     figure('Name','Check Line and Diff')
%     plot(Lines(:,i));
%     hold on
%     plot(diffline);
%     hold off;
%     legend('Line','Difference')
    dLines(1:length(diffline),i) = abs(diff(Lines(:,i)));
    [~,diaphragm_pos(i)] = max(dLines(:,i));
end

figure('Name','Diaphragm Motion')
imagesc(Lines)
colormap(gray)
axis off
hold on
plot(1:size(All_Ims,3),diaphragm_pos,'g','LineWidth',2);
hold off








