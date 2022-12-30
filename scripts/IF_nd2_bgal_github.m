function IF_nd2_bgal_github(dependencies, currFile, datadirIF, mask_dir, num_stacks, nuc_channel, cyto_channel, R, G, B, cytothresh, bgal_pct)
addpath(dependencies)
%% Resolve file name
[experimentpath,fn] = fileparts(currFile);

%For tiffs exported from ND2 files via elements:
startIdx = regexp(fn,'[0-9][0-9]');
rows = {'a','b','c','d','e','f','g','h'}; row_idx = startIdx(1) - 1;
row_name = (fn(row_idx)); WellRow=find(ismember(rows,row_name));

col_idx = startIdx(1); col_idx2 = col_idx+1; 
col_str = fn(col_idx); col_str2 = fn(col_idx2);
col_cat = strcat(col_str,col_str2); WellCol = str2num(col_cat);

fn_size = size(fn,2); WellSite = str2num(fn(fn_size));
%% load image data 
timetotal=tic;

raw_IF = [];
for L = 1:num_stacks
    raw_IF(:,:,L) = imread(currFile,L);
end
%% segment nuclei
raw_IF = double(raw_IF);
nucr=12;
debrisarea=100;
boulderarea=1500; 
blurradius=3;

nuc_mask=threshmask(raw_IF(:,:,nuc_channel),blurradius);
nuc_mask=markershed(nuc_mask,round(nucr*2/3));
nuc_mask=secondthresh(raw_IF(:,:,nuc_channel),blurradius,nuc_mask,boulderarea*2);
nuc_mask=bwareaopen(nuc_mask,debrisarea);
nuc_mask=segmentdeflections_bwboundaries(nuc_mask,nucr,debrisarea);
nuc_mask=excludelargeandwarped_3(nuc_mask,boulderarea,0.90);

nuc_perim = bwperim(nuc_mask);
raw_nuc = raw_IF(:,:,nuc_channel);
figure()
himage = imshowpair(nuc_perim,raw_nuc);
saveas(himage,[mask_dir, num2str(WellRow), '_', num2str(WellCol), '_', num2str(WellSite),'_nuc.tif']);
cla reset
%% segment cytoplasm
real_IF = raw_IF;
cyto = real_IF(:,:,cyto_channel);

cyto_thresh = cyto > cytothresh;
cyto = mat2gray(cyto);
% clear small objects
cyto_nuc_cleared = bwareaopen(cyto_thresh, 600);
cyto_mask = cyto_nuc_cleared;

% Watershed
seeds = nuc_mask;
%cyto_smooth = imgaussfilt(cyto,1); % don't smooth too much for watersheding
cyto_min = imimposemin(max(cyto(:))-cyto,seeds); % set locations of seeds to be -Inf (cause matlab watershed)
cyto_ws = watershed(cyto_min);
cyto_ws(cyto_mask==0)=0; % remove areas that aren't in our cyto mask

% Fill holes
labelled_cyto = imfill(cyto_ws,'holes');
cyto_info = regionprops(labelled_cyto,'PixelIdxList','Centroid');  
a = imshow(imoverlay(uint8(cyto*255),bwperim(labelled_cyto),'r'),[]);
hold on
% Display color overlay

labelled_cyto_rgb = label2rgb(uint32(labelled_cyto), 'jet', [1 1 1], 'shuffle');
himage = imshow(labelled_cyto_rgb,[]);
himage.AlphaData = cyto_mask*.3;
saveas(himage,[mask_dir, num2str(WellRow), '_', num2str(WellCol), '_', num2str(WellSite),'_cyto.tif']);
close all
%% feature extraction
nuc_label=bwlabel(nuc_mask);
nuc_info=struct2cell(regionprops(nuc_mask,real_IF(:,:,1),'Area','Centroid','MeanIntensity')');
nuc_area=squeeze(cell2mat(nuc_info(1,1,:)));
nuc_info=regionprops(nuc_label,'PixelIdxList','Centroid','PixelList');
cyto_info = regionprops(labelled_cyto,'PixelIdxList','PixelList');
numcells=numel(nuc_area);

real4 = real_IF(:,:,R);
real5 = real_IF(:,:,G);
real6 = real_IF(:,:,B);

i=1;
for L = 1:size(cyto_info)
    red{1,i} = real6(cyto_info(i).PixelIdxList);
    green{1,i} = real5(cyto_info(i).PixelIdxList);
    blue{1,i} = real4(cyto_info(i).PixelIdxList);
    i=i+1;
end
red = transpose(red);
blue = transpose(blue);
green = transpose(green);

for cc=1:numcells
    cells(cc).bgal_red = red(cc);
    cells(cc).bgal_green = green(cc);
    cells(cc).bgal_blue = blue(cc);
end

channel = red;
index = 1;
red_bgal = [];
for n = 1:size(red,1)
z = 1;
Y = [];
    for n = 1:size(size(channel{index,1},1))
        Y{z,1} = 1/(prctile(channel{index,1},bgal_pct));
        z=z+1;
        X = cell2mat(Y);
        red_bgal{index,1} = X;
    end
index = index +1;
end

IFdata_info_struc = nuc_info;
PixelList_cyto = [];
z=1;
for L = 1:size(cyto_info,1)
    PixelList_cyto{z,1} = cyto_info(z,1).PixelList;
    PixelIDXList_cyto{z,1} = cyto_info(z,1).PixelIdxList;
    z=z+1;
end
[IFdata_info_struc(:).PixelList_cyto] = PixelList_cyto{:}; 
[IFdata_info_struc(:).PixelIDXList_cyto] = PixelIDXList_cyto{:}; 
[IFdata_info_struc(:).bgal] = red_bgal{:};   
%% store data 
if ~isempty(nuc_area)
    save([datadirIF,num2str(WellRow), '_', num2str(WellCol), '_', num2str(WellSite),'_IF.mat'],'cells','IFdata_info_struc');
end

toc(timetotal);
