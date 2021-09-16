function [dFix, dMove, SDeform, dImgReg, cVoxelInterp] = fRegDemons(dFix, dMove, sParafile, iDim, SGeo)
% function to run demons registration method
%
% input:
% dFix          reference/fixed image (2D/3D)
% dMove         moving images (2D/3D): x-y-(z)-t
% sParafile     (string) path to parameter file
% iDim          registration dimensionality (1=2D, 2=3D)
% SGeo          (struct) geometric image information (voxelsize, orientation, ...)
%
% output:
% dFix          reference/fixed image (2D/3D)
% dMove         moving images (2D/3D): x-y-(z)-t
% SDeform       (struct) deformation field (time x 1) with forward (F) and backward (B) fields
% dImgReg       transformed image from dMove towards dFix according to SDeform
% cVoxelInterp  image resolution of (interpolated) output images
%
% -------------------------------------------------------------------------
% (c) 2015: Thomas Kuestner, Verena Neumann
% -------------------------------------------------------------------------

%%
h = fwaitbar(0,'Running Registration. Please wait!'); st=0;

%% prepare images
[SOptions, dDeformRes] = fReadDemonsParam(sParafile);
nDimImg = ndims(dFix);

if dDeformRes ~= 0
    % interpolate to isotropic resolution
    % define MR Coordinate System
    SGeo.dMRCoord1 = ((1:size(dFix, 1))' - size(dFix, 1)./2 - 0.5).*SGeo.cVoxelsize{1}(1);
    SGeo.dMRCoord2 = ((1:size(dFix, 2))' - size(dFix, 2)./2 - 0.5).*SGeo.cVoxelsize{1}(2);
    SGeo.dMRCoord3 = ((1:size(dFix, 3))' - size(dFix, 3)./2 - 0.5).*SGeo.cVoxelsize{1}(3);
    [SGeo.dMRCoord1, SGeo.dMRCoord2, SGeo.dMRCoord3] = ndgrid(SGeo.dMRCoord1, SGeo.dMRCoord2, SGeo.dMRCoord3);
    % The interpolation coordinate system
    dXi = SGeo.dMRCoord1(1,1,1):sign(SGeo.dMRCoord1(2,1,1)-SGeo.dMRCoord1(1,1,1))*dDeformRes:SGeo.dMRCoord1(end,1,1);
    dYi = SGeo.dMRCoord2(1,1,1):sign(SGeo.dMRCoord2(1,2,1)-SGeo.dMRCoord2(1,1,1))*dDeformRes:SGeo.dMRCoord2(1,end,1);
    if(nDimImg == 2)
%         dZi = 0;
        [dXIso, dYIso] = ndgrid(dXi, dYi);
        dImgIsoTmp = zeros(size(dXIso, 1), size(dXIso, 2), 1, size(dMove, 4));
        % The Interpolation
        fprintf('\nInterpolating MR Images to deformation %1.2f mm isotropic resolution.', dDeformRes);
        dFix = interpn(SGeo.dMRCoord1, SGeo.dMRCoord2, dFix, dXIso, dYIso, '*linear');    
        for iI = 1:size(dMove, 4)
            dImgIsoTmp(:,:,:,iI) = interpn(SGeo.dMRCoord1, SGeo.dMRCoord2, dMove(:,:,:,iI), dXIso, dYIso, '*linear');
            fprintf(1, '.');
        end
        fprintf('\n');
        dMove = dImgIsoTmp;
    else %3D
        dZi = SGeo.dMRCoord3(1,1,1):sign(SGeo.dMRCoord3(1,1,2)-SGeo.dMRCoord3(1,1,1))*dDeformRes:SGeo.dMRCoord3(1,1,end);    
        [dXIso, dYIso, dZIso] = ndgrid(dXi, dYi, dZi);
        dImgIsoTmp = zeros(size(dXIso, 1), size(dXIso, 2), size(dXIso, 3), size(dMove, 4));
        % The Interpolation
        fprintf('\nInterpolating MR Images to deformation %1.2f mm isotropic resolution.', dDeformRes);
        dFix = interpn(SGeo.dMRCoord1, SGeo.dMRCoord2, SGeo.dMRCoord3, dFix, dXIso, dYIso, dZIso, '*linear');    
        for iI = 1:size(dMove, 4)
            dImgIsoTmp(:,:,:,iI) = interpn(SGeo.dMRCoord1, SGeo.dMRCoord2, SGeo.dMRCoord3, dMove(:,:,:,iI), dXIso, dYIso, dZIso, '*linear');
            fprintf(1, '.');
        end
        fprintf('\n');
        dMove = dImgIsoTmp;
    end
    clear 'dImgIsoTmp';
    cVoxelInterp = mat2cell(repmat(dDeformRes,size(dMove,4)+1,3),ones(1,size(dMove,4)+1),3);
else
    cVoxelInterp = SGeo.cVoxelsize;
end

if(nDimImg == 2)
    iNGates = size(dMove, 3) + 1;
    iN3D = 1;
    dImgReg = zeros(size(dFix,1), size(dFix,2), 1, iNGates); % always 4D array: x-y-z-t
elseif(nDimImg == 3)
    iNGates = size(dMove, 4) + 1;
    if(iDim == 1) % 2D reg
        iN3D = size(dFix,3);
    else % 3D reg
        iN3D = 1;
    end
    dImgReg = zeros(size(dFix,1), size(dFix,2), size(dFix,3), iNGates);
end
steps = iN3D * (iNGates-1); % for waitbar

dFix = (dFix - min(dFix(:)))./(max(dFix(:)) - min(dFix(:)));
dMove = (dMove - min(dMove(:)))./(max(dMove(:)) - min(dMove(:)));

%% start registration
for iJ = 1:iN3D
    for iI = 2:iNGates
        fprintf(1, 'Registering gate %u...', iI);
        
        if(~exist('SDeform','var') || length(SDeform) < iI || ~isfield(SDeform(iI),'dFy'))
            SDeform(iI).dFy = zeros(size(dImgReg,1),size(dImgReg,2),size(dImgReg,3));
            SDeform(iI).dFx = zeros(size(dImgReg,1),size(dImgReg,2),size(dImgReg,3));
            SDeform(iI).dFz = zeros(size(dImgReg,1),size(dImgReg,2),size(dImgReg,3));

            SDeform(iI).dBy = zeros(size(dImgReg,1),size(dImgReg,2),size(dImgReg,3));
            SDeform(iI).dBx = zeros(size(dImgReg,1),size(dImgReg,2),size(dImgReg,3));
            SDeform(iI).dBz = zeros(size(dImgReg,1),size(dImgReg,2),size(dImgReg,3));                
        end
        if(iDim == 1) % 2D reg
            if(nDimImg == 2) % 2D image
                dImgReg(:,:,iJ,1) = dFix;
                [dImgReg(:,:,iJ,iI), dBx, dBy, dFx, dFy] = register_images(dMove(:,:,iI-1), dFix, SOptions);                
   
            elseif(nDimImg == 3) % 3D image
                dImgReg(:,:,iJ,1) = dFix(:,:,iJ);
                [dImgReg(:,:,iJ,iI), dBx, dBy, dFx, dFy] = register_images(dMove(:,:,iJ,iI-1), dFix(:,:,iJ), SOptions);   
            end
        
            if dDeformRes ~= 0
                SDeform(iI).dFy(:,:,iJ) = dFx.*dDeformRes;
                SDeform(iI).dFx(:,:,iJ) = dFy.*dDeformRes;
                
                SDeform(iI).dBy(:,:,iJ) = dBx.*dDeformRes;
                SDeform(iI).dBx(:,:,iJ) = dBy.*dDeformRes;
            else
                SDeform(iI).dFy(:,:,iJ) = dFx;
                SDeform(iI).dFx(:,:,iJ) = dFy;
                
                SDeform(iI).dBy(:,:,iJ) = dBx;
                SDeform(iI).dBx(:,:,iJ) = dBy;
            end
        
        elseif(iDim == 2) % 3D reg => just 3D image
            dImgReg(:,:,:,1) = dFix;
            [dImgReg(:,:,:,iI), dBx, dBy, dBz, dFx, dFy, dFz] = register_volumes(dMove(:,:,:,iI-1), dFix, SOptions);
                        
            if dDeformRes ~= 0
                SDeform(iI).dFy = dFx.*dDeformRes;
                SDeform(iI).dFx = dFy.*dDeformRes;
                SDeform(iI).dFz = dFz.*dDeformRes;
                
                SDeform(iI).dBy = dBx.*dDeformRes;
                SDeform(iI).dBx = dBy.*dDeformRes;
                SDeform(iI).dBz = dBz.*dDeformRes;
            else
                SDeform(iI).dFy = dFx;
                SDeform(iI).dFx = dFy;
                SDeform(iI).dFz = dFz;
                
                SDeform(iI).dBy = dBx;
                SDeform(iI).dBx = dBy;
                SDeform(iI).dBz = dBz;
            end      
        end
        
        st= st+1;fwaitbar(st/steps,h);
    end
end
try close(h); catch; end;
