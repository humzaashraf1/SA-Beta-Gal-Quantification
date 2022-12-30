function betagal_quantification_validation(currFile, datadirIF, mask_dir)
timetotal=tic;

[experimentpath,fn] = fileparts(currFile);
startIdx = regexp(fn,'[0-9][0-9]');
rows = {'a','b','c','d','e','f','g','h'}; row_idx = startIdx(1) - 1;
row_name = (fn(row_idx)); WellRow=find(ismember(rows,row_name));

col_idx = startIdx(1); col_idx2 = col_idx+1; 
col_str = fn(col_idx); col_str2 = fn(col_idx2);
col_cat = strcat(col_str,col_str2); WellCol = str2num(col_cat);

fn_size = size(fn,2); WellSite = str2num(fn(fn_size));

shot = [num2str(WellRow),'_',num2str(WellCol),'_',num2str(WellSite)];
z=1;
for L = 1:6
    raw_IF(:,:,z) = imread(currFile,z);
    z=z+1;
end
raw_IF = double(raw_IF);
load([datadirIF,shot,'_IF','.mat'],'cells','IFdata_info_struc');

nuc_mask_test = zeros([2044 2048]);
count = 1;
for L = 1:size(IFdata_info_struc,1)
    z = 1;
    for N = 1:size(IFdata_info_struc(count).PixelList_cyto,1)
        nuc_mask_test(((IFdata_info_struc(count).PixelList_cyto(z,2))),(IFdata_info_struc(count).PixelList_cyto(z,1))) = 1;
        z=z+1;
    end
    count=count+1;
end
% Reconstruct cyto mask and apply to RGB image
nuc_mask_test = uint16(nuc_mask_test);

in = 25000;
out = 50000;
r = uint16(raw_IF(:,:,6));
g = uint16(raw_IF(:,:,5));
b = uint16(raw_IF(:,:,4));
in_range = in/65536;
out_range = out/65536;
r_ad = imadjust(r,[in_range out_range]);
g_ad = imadjust(g,[in_range out_range]);
b_ad = imadjust(b,[in_range out_range]);

rgbImage = cat(3, r, g, b);
rgbImage2 = cat(3, r_ad, g_ad, b_ad);

labelled_cyto = nuc_mask_test.*rgbImage2;
figure()
himage = imshow((labelled_cyto));
saveas(himage,[mask_dir, num2str(WellRow), '_', num2str(WellCol), '_', num2str(WellSite),'_bgal.tif']);
cla reset

bgal = [IFdata_info_struc(:).bgal]';
bgal = rescale(bgal);
PixelidxCyto = {IFdata_info_struc(:).PixelIDXList_cyto}'; 

nuc_mask_test = zeros([2044 2048]);
for L = 1:size(PixelidxCyto,1)
    x = bgal(L);
    nuc_mask_test(PixelidxCyto{L}) = x;
end

himage = imshow(nuc_mask_test);
colormap turbo
caxis([0 0.3])
saveas(himage,[mask_dir, num2str(WellRow), '_', num2str(WellCol), '_', num2str(WellSite),'_bgalsignal.tif']);
close all

toc(timetotal);
end
