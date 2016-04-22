function lung3d = fLungSegmKohlm(dImg,nGate,bChan)
% Function to segment Lungs according to paper of Kohlmann et al. 2014
% Determination of seed points by two-step histogram-threshold segmentation
% and subsequent region growing from the seed.
%
% Region growing mask can be used as input for Chan-Vese segmentation
%
% input:    
% dImg          Image data (4-D double)
% nGate         gate number (int)
% bChan         flag if Chan-Vese algorithm shall be used (binary)
% 
% output:   
% lung3D        3-D logical containing the voxel, that are associated with lung tissue
%
% -------------------------------------------------------------------------
% (c) 2015: Thomas Kuestner, Verena Neumann
% -------------------------------------------------------------------------

% flag to show results in figure
showfig = 0;

% parameter for region growing algorithm:
% dMaxDiff = 0.0033;
dMaxDiff = 0.0025;
% Chan-Vese parameter
cParams = {1, 1, 0.2, 10, 2};   % lambda1: inside, lambda2: outside
                                % dNu:     , dDeltT:     , dNIter:
% take middle slice
middleSlice = floor(size(dImg,3)/2);
dImg_sl = squeeze(dImg(:,:,middleSlice,nGate));

% compute histogram
[y1,x1] = hist(dImg_sl(:),150);

% compute numerical derivative
deriHist = diff(y1);

% get first index for derivative passing zero
grZero = find(deriHist>0);
% minimum value for the thorax mask
minVal_th = x1(grZero(1));     

% generate the thorax mask
% a) get raw mask
thoraxMask_wh = (dImg_sl>minVal_th);

% b) fill holes
thoraxMask_f = imfill(thoraxMask_wh,'holes');

% c) dilate Mask to fill gaps
SE = strel('rectangle',[3 3]);
thoraxMask_d = imdilate(thoraxMask_f,SE);

% fill again
thoraxMask_f = imfill(thoraxMask_d,'holes');

% d) take only largest connected component
thoraxMask = getLargestCc(thoraxMask_f);

% mask image with thoraxMask
dImg_m = dImg_sl.*thoraxMask;

% according to Kohlmann smooth image for lugn mask:
dImg_m = medfilt2(dImg_m);

% get histogram of masked image and numerical derivative
[y1,x1] = hist(dImg_m(:),150);

% % smooth the histogram (Kohlmann)
% y12 = filter(ones(1,3)/3,1,y1);
% y1 = circshift(y12',-2);

deriHist = diff(y1);
% get first index for derivative passing zero, skipping the first k=5 bins
k = 5;
deriHist(1:k) = 0;
grZero = find(deriHist>0);
minVal_lu = x1(grZero(1));

% generate lung mask
% a) set outer area to minVal_lu
dImg_m(dImg_m == 0)= minVal_lu;

% b) take only values smaller than minVal_lu
lungMask_wh = (dImg_m < minVal_lu);

% c) Opening: erode and dilate image (delete connection between left and
% right part of the lung)
SE = strel('rectangle',[4 4]);
lungMask_e = imerode(lungMask_wh,SE);
lungMask_e = imdilate(lungMask_e,SE);

% d) take 2 largest segments
lungMask_right = getLargestCc(lungMask_e,[],1);
lungMask_left = getLargestCc(lungMask_e,[],2)-lungMask_right;

% e) fill the segments
lungMask_left = imfill(lungMask_left,'holes');
lungMask_right = imfill(lungMask_right,'holes');

% get Euclidean distance to border for every point of the mask
dist_lungMask_left = bwdist(imcomplement(lungMask_left));
dist_lungMask_right = bwdist(imcomplement(lungMask_right));

% get Maximum indices
[vl,ind] = max(dist_lungMask_left(:));
[yl,xl] = ind2sub(size(dist_lungMask_left),ind);
[vr,ind] = max(dist_lungMask_right(:));
[yr,xr] = ind2sub(size(dist_lungMask_right),ind);

% get darkest point in ROI around maximum distance points
% a) circular ROIs
ROIrad = min([vl,vr])/2;
[yy,xx] = ndgrid((1:size(dImg,2))-yr,(1:size(dImg,1))-xr);
ROIcircle_r = (xx.^2 + yy.^2)< ROIrad^2;
[yy,xx] = ndgrid((1:size(dImg,2))-yl,(1:size(dImg,1))-xl);
ROIcircle_l = (xx.^2 + yy.^2)< ROIrad^2;

% b) darkest point inside
leftROI = ROIcircle_l.*dImg_sl;
leftROI(leftROI == 0)= 1;
[~,ind] = min(leftROI(:));
[yld,xld] = ind2sub(size(leftROI),ind);
rightROI = ROIcircle_r.*dImg_sl;
rightROI(rightROI == 0)= 1;
[~,ind]= min(rightROI(:));
[yrd,xrd] = ind2sub(size(rightROI),ind);

lSeed = [yld,xld,middleSlice];
rSeed = [yrd,xrd,middleSlice];

% region growing from the seed points
% normalize image data
dImg_vol = dImg(:,:,:,nGate);
maxImg = max(dImg_vol(:));
normImg = (dImg_vol./maxImg);


% right lung
lung3dr = RegionGrowing(normImg, dMaxDiff, rSeed);
SE = strel('rectangle',[5 5]);
lung3dr = imdilate(lung3dr,SE);
lung3dr = imerode(lung3dr,SE);
% left lung
lung3dl = RegionGrowing(normImg, dMaxDiff, lSeed);
lung3dl = imdilate(lung3dl,SE);
lung3dl = imerode(lung3dl,SE);

% both lungs in one array, fill inner holes in plane
lung3d = logical(lung3dr+lung3dl);
for iI = 1:size(lung3dr,3)
    lung3d(:,:,iI) = imfill(lung3d(:,:,iI),'holes');
end

if bChan
lung3d = ChanVese('exe', dImg_vol, lung3d, 0, cParams);
end

if showfig
    dImgGate = dImg(:,:,:,nGate);
    normImg_m = dImgGate./max(dImgGate(:))*2;
    rgbImg = (normImg_m(:,:,:,[1 1 1])); % make the first a grey-scale image with three channels so it will not be affected by the colormap later on
    rgbImg_r = (lung3d(:,:,:)).*(rgbImg(:,:,:,1)+ones(size(lung3d,1),size(lung3d,2),size(lung3d,3))/2);
    rgbImg_pm = rgbImg;
    rgbImg_pm(:,:,:,2) = rgbImg(:,:,:,2) + rgbImg_r;
    h(nGate).fig = figure;
    h(nGate).num2str(nGate) = imshow(reshape(rgbImg_pm(:,:,middleSlice,:),256,256,3));
    hold on;
    plot(xld,yld,'y+', 'Linewidth', 1.5,'MarkerSize', 20);
    plot(xrd,yrd,'y+','Linewidth', 1.5, 'MarkerSize', 20);
    xX = 1:256; yY = 100*ones(256,1);
    h(nGate).plot = plot(xX,yY,'r');
    
    hsl = uicontrol('Style','slider','Min',1,'Max',size(dImg,3),...
                'SliderStep',[1 1]./72,'Value',36,...
                'Position',[20 20 200 20]);
    set(hsl,'Callback',@(hObject,eventdata) ...
    set(h(nGate).num2str(nGate),'CData',imshow(reshape(rgbImg_pm(:,:,round(get(hObject,'Value')),:),256,256,3))))
    
    hsl2 = uicontrol('Style','slider','Min',1,'Max',size(dImg,2),...
                'SliderStep',[1 1]./256,'Value',156,...
                'Position',[20 120 20 200]);
    set(hsl2,'Callback',@(hObject,eventdata) set(h(nGate).plot,'YData',256-round(get(hObject,'Value')*ones(256,1))))
        
end


end