function euclid = fComparePoints(SRegiResult,Labels)
% compare points in LandmarkGUI of original images
%
% input:
% SRegiResult       struct with registration info
% Labels            assigned label positions
%
% output:
% euclid            euclidean distance between points in reference and moving image
%
% -------------------------------------------------------------------------
% (c) 2015: Thomas Kuestner, Verena Neumann
% -------------------------------------------------------------------------

results = Labels;
DefField = SRegiResult.SDeform;

dVoxSz = SRegiResult.SGeo.dVoxelsize;

% euclid = cell(3,1);

% get the coordinates of the fixed points that shall be compared
G01 = results.fixedLMP.G01; x01 = G01(:,2); y01 = G01(:,3); z01 = G01(:,4);
G02 = results.fixedLMP.G02; x02 = G02(:,2); y02 = G02(:,3); z02 = G02(:,4);
G03 = results.fixedLMP.G03; x03 = G03(:,2); y03 = G03(:,3); z03 = G03(:,4);
G04 = results.fixedLMP.G04; x04 = G04(:,2); y04 = G04(:,3); z04 = G04(:,4);


% transform the coordinates of the moving images according to
% registration result
for iI = 2:4
    dBx(:,:,:,iI) = DefField(iI).dBx;
    dBy(:,:,:,iI) = DefField(iI).dBy;
    dBz(:,:,:,iI) = DefField(iI).dBz;
end
for i = 1:length(x02)
    x02n(i) = x02(i) + dBx(x02(i),y02(i),z02(i),2);
    y02n(i) = y02(i) + dBy(x02(i),y02(i),z02(i),2);
    z02n(i) = z02(i) + dBz(x02(i),y02(i),z02(i),2);
    
    x03n(i) = x03(i) + dBx(x03(i),y03(i),z03(i),3);
    y03n(i) = y03(i) + dBy(x03(i),y03(i),z03(i),3);
    z03n(i) = z03(i) + dBz(x03(i),y03(i),z03(i),3);
    
    x04n(i) = x04(i) + dBx(x04(i),y04(i),z04(i),4);
    y04n(i) = y04(i) + dBy(x04(i),y04(i),z04(i),4);
    z04n(i) = z04(i) + dBz(x04(i),y04(i),z04(i),4);
end

% compute the euclidean distance
euclid(1,:) = sqrt(((x01-x02n')*dVoxSz(1)).^2+((y01-y02n')*dVoxSz(2)).^2+((z01-z02n')*dVoxSz(3)).^2);
euclid(2,:) = sqrt(((x01-x03n')*dVoxSz(1)).^2+((y01-y03n')*dVoxSz(2)).^2+((z01-z03n')*dVoxSz(3)).^2);
euclid(3,:) = sqrt(((x01-x04n')*dVoxSz(1)).^2+((y01-y04n')*dVoxSz(2)).^2+((z01-z04n')*dVoxSz(3)).^2);

