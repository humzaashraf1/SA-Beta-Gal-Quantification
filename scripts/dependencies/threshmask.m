function mask=threshmask(image,blurradius)
blur=imfilter(image,fspecial('disk',blurradius),'symmetric'); %10x:3 20x:6
normlog=mat2gray(log(blur));
thresh=graythresh(normlog); %thresh for this image is 0.333.
%first attempt: increase thresh to 0.5
%thresh = 0.4;
%thresh=graythresh(normlog).*1.5;
%20161205 edit: for the Fab movies, the thresh is being set too high.
%Instead of using graythresh, find the mean and use that
mask=im2bw(normlog,thresh);
mask=imfill(mask,'holes');
end