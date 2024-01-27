%%Read and preprocess the image
im_orig = imread(full_path); % read an image
im_d = 255*im2double(im_orig); % convert to double in range 0, 255
medi1 = median(im_d(:)); % calculate the median
im = medi1-im_orig; % subtract the median from the original image
im = im2double(im); % convert to double
im_dens = 1-im2double(im_orig); % invert image
im_smoothed = imgaussfilt(im,4); % apply a gaussian filter with std = 4
BW = im2bw(im_smoothed,0.08); % convert to binary with threshold 0.08

%%Remove objects with area < 4000 pixels
thres = 4000; % set the threshold
CC = bwconncomp(BW); % compute the connected components
if CC.NumObjects > 1
    for ii = 1:CC.NumObjects
        numPixels = CC.PixelIdxList{ii};
        if length(numPixels) < thres
            BW(CC.PixelIdxList{ii}) = 0;
        end
    end
end

%%Calculate region properties (eccentricity). Keep only the object with
%%lowest eccentricity
stats = regionprops(BW,'Eccentricity','PixelIdxList'); 
if length(stats) > 1
    [val,ind]= min([stats.Eccentricity]);
    BW_filt = BW;
    for j = 1:length(stats)
        if j ~= ind
            BW_filt(stats(j).PixelIdxList) = 0;
        end
    end
    BW = BW_filt;
end

%%Calculate the center of mass of the remaining binary image
[rows,cols,vals] = find(BW == 1);
center = [round(mean(rows)), round(mean(cols))];

%%Calculate the density and distance from tumor center for all pixels (our images are 1152 by 1536 pixels)
pair = zeros(1152*1536,2); % create an array which stores [pixelID, image density]
for k = 1:1152
    for j = 1:1536
        distance = sqrt((k-center(1))^2+(j-center(2))^2); % distance between tumor center and current pixel
        pair(j + (k-1)*1536,:) = [distance, im_dens(k,j)];% store the [distance, density] pair
    end
end

%%Discretize domain into bins of 2 pixels length each.
density = [];
delta_dist = 2; % bin size (2 pixels). 
max_dist = 1200; % maximum allowed distance between tumor center and a given pixel
for dist = 1:delta_dist:max_dist
    [ind1,val] = find(pair(:,1) > dist);
    [ind2,val2] = find(pair(:,1) <= dist+delta_dist);
    ind = intersect(ind1,ind2);
    count = sum(pair(ind,2));
    density = [density; count/length(ind)];
end
% the resulting spatial resolution is that of two pixels (delta_dist = 2)