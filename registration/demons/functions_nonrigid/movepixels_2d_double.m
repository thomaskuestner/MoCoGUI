function dImgOut = movepixels_2d_double(dImgIn, dTx, dTy, iMode)
% This function movepixels, will translate the pixels of an image
%  according to x and y translation images (bilinear interpolated).
%
%  Iout = movepixels_2d_double(I,Tx,Ty,mode);
%
% Inputs;
%   Tx, Ty: The transformation images, describing the
%             (backwards) translation of every pixel in x and y direction.
%   mode: If 0: linear interpolation and outside pixels set to nearest pixel
%            1: linear interpolation and outside pixels set to zero
%            (cubic interpolation only supported by compiled mex file)
%            2: cubic interpolation and outsite pixels set to nearest pixel
%            3: cubic interpolation and outside pixels set to zero
%
% Outputs,
%   Iout : The transformed image
%
% Function is written by D.Kroon University of Twente (February 2009)

if (any(size(dImgIn) ~= size (dTx))) || (any(size(dImgIn) ~= size (dTy)))
    error('Movepixels->Inputs 1 to 3 must have same dimentions!');
end

% Make all x,y indices
[dY, dX] = ndgrid(0:size(dImgIn, 1) - 1, 0:size(dImgIn, 2) - 1);
dImgOut = zeros(size(dImgIn));

for iI = 1:size(dImgIn, 3);

    % Calculate the Transformed coordinates
    dXInd = dX + dTx(:,:,iI);
    dYInd = dY + dTy(:,:,iI);
    
    % All the neighbour pixels involved in linear interpolation.
    iXFlor = floor(dXInd);
    iYFlor = floor(dYInd);
    iXCeil = iXFlor + 1;
    iYCeil = iYFlor + 1;
    
    % Linear interpolation constants (percentages)
    dXRem = dXInd - double(iXFlor);
    dYRem = dYInd - double(iYFlor);
    
    dPerc0 = (1 - dXRem) .* (1 - dYRem);
    dPerc1 = (1 - dXRem) .*      dYRem ;
    dPerc2 =      dXRem  .* (1 - dYRem);
    dPerc3 =      dXRem  .*      dYRem ;
    
    % limit indexes to boundaries
    lMaskXBas0 = (iXFlor < 0)  |  (iXFlor > (size(dImgIn,1) - 1));
    lMaskYBas0 = (iYFlor < 0)  |  (iYFlor > (size(dImgIn,2) - 1));
    lMaskXBas1 = (iXCeil < 0)  |  (iXCeil > (size(dImgIn,1) - 1));
    lMaskYBas1 = (iYCeil < 0)  |  (iYCeil > (size(dImgIn,2) - 1));
    
    iXCeil(lMaskXBas1) = 0;
    iYCeil(lMaskYBas1) = 0;
    iXFlor(lMaskXBas0) = 0;
    iYFlor(lMaskYBas0) = 0;
        
    % Get the intensities
    dThisImageIn = dImgIn(:,:,iI);
    dSum0 = dThisImageIn(1 + iYFlor + iXFlor*size(dImgIn, 1));
    dSum1 = dThisImageIn(1 + iYFlor + iXCeil*size(dImgIn, 1));
    dSum2 = dThisImageIn(1 + iYCeil + iXFlor*size(dImgIn, 1));
    dSum3 = dThisImageIn(1 + iYCeil + iXCeil*size(dImgIn, 1));
    
    % Make pixels before outside Ibuffer mode
    if (iMode == 1 || iMode == 3)
        dSum0(lMaskXBas0 | lMaskYBas0) = 0;
        dSum1(lMaskXBas0 | lMaskYBas1) = 0;
        dSum2(lMaskXBas1 | lMaskYBas0) = 0;
        dSum3(lMaskXBas1 | lMaskYBas1) = 0;
    end
    dThisImgOut = dSum0.*dPerc0 + dSum1.*dPerc1 + dSum2.*dPerc2 + dSum3.*dPerc3;
    dImgOut(:,:,iI) = reshape(dThisImgOut, [size(dImgIn,1) size(dImgIn,2)]);
end