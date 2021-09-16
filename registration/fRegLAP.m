function [dFix, dMove, SDeform, dImgReg, cVoxelInterp] = fRegLAP(dFix, dMove, sParafile, iDim, SGeo)
% function to run LAP registration method
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
if(ischar(sParafile))
    [ SOptions, dDeformRes ] = fReadLAPParam(sParafile);
else
    dDeformRes = 0;
    SOptions = sParafile;
end

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
    clear 'dImgIsoTmp';
    cVoxelInterp = mat2cell(repmat(dDeformRes,size(dMove,4)+1,3),ones(1,size(dMove,4)+1),3);
else
    cVoxelInterp = SGeo.cVoxelsize;
end

nDimImg = ndims(dFix);
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
currpath = fileparts(mfilename('fullpath'));
if(iDim == 1) % 2D reg
    addpath(genpath([currpath,filesep,'LAP',filesep,'2D']));
    rmpath(genpath([currpath,filesep,'LAP',filesep,'3D']));
else % 3D reg
    addpath(genpath([currpath,filesep,'LAP',filesep,'3D']));
    rmpath(genpath([currpath,filesep,'LAP',filesep,'2D']));
end
% dImgReg(:,:,:,1) = dFix;
steps = iN3D * (iNGates-1); % for waitbar

dFix = (dFix - min(dFix(:)))./(max(dFix(:)) - min(dFix(:)));
dMove = (dMove - min(dMove(:)))./(max(dMove(:)) - min(dMove(:)));
dImgReg(:,:,:,1) = dFix;

% Parameters for optical flow estimation:
FilterSizes =  2.^(SOptions.PyrMax:-1:SOptions.PyrMin);   % scale levels (assumes maximum flow = 2^4 pixels)

%% start registration
for iJ = 1:iN3D
    for iI = 2:iNGates
        fprintf('\nRegistration of image %g/%g to reference image.\n',iI-1,size(dImgReg,4)-1);

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
                uLAP = CG_MultiScale_LAP(dFix, dMove(:,:,iI-1), SOptions.PyrMax, SOptions.Nf);
                dImgReg(:,:,iJ,iI) = imshift(dMove(:,:,iI-1),{-uLAP{1},-uLAP{2}});    
            elseif(nDimImg == 3) % 3D image
                dImgReg(:,:,iJ,1) = dFix(:,:,iJ);
                uLAP = CG_MultiScale_LAP(dFix, dMove(:,:,iJ,iI-1), SOptions.PyrMax, SOptions.Nf);
                dImgReg(:,:,iJ,iI) = imshift(dMove(:,:,iJ,iI-1),{-uLAP{1},-uLAP{2}});    
            end
            
            if dDeformRes ~= 0
                SDeform(iI).dFy(:,:,iJ) = real(uLAP).*dDeformRes;
                SDeform(iI).dFx(:,:,iJ) = imag(uLAP).*dDeformRes;
                
                SDeform(iI).dBy(:,:,iJ) = -real(uLAP).*dDeformRes;
                SDeform(iI).dBx(:,:,iJ) = -imag(uLAP).*dDeformRes;
            else
                SDeform(iI).dFy(:,:,iJ) = real(uLAP);
                SDeform(iI).dFx(:,:,iJ) = imag(uLAP);
                
                SDeform(iI).dBy(:,:,iJ) = -real(uLAP);
                SDeform(iI).dBx(:,:,iJ) = -imag(uLAP);
            end
            
        else % 3D reg => just 3D image
            dImgReg(:,:,:,1) = dFix;
            if strcmp(SOptions.Algorithm, 'fast')
                uLAP = CG_MultiScale_LAP3D_Fast(dFix, dMove(:,:,:,iI-1), FilterSizes, SOptions.PreFilt, SOptions.MedFilt);
                dImgReg(:,:,:,iI) = imshift_3D(dMove(:,:,:,iI-1),{-uLAP{1},-uLAP{2},-uLAP{3}},SOptions.Interpolation);    
            else
                % slower
                uLAP = CG_MultiScale_LAP3D(dFix, dMove(:,:,:,iI-1), FilterSizes, SOptions.PreFilt, SOptions.MedFilt);
                dImgReg(:,:,:,iI) = imshift_3D(dMove(:,:,:,iI-1),{-uLAP{1},-uLAP{2},-uLAP{3}},SOptions.Interpolation);
                % forward
%                 dImgReg(:,:,:,iI) = imshift_3D(dFix,uLAP,SOptions.Interpolation);
            end
            
            if dDeformRes ~= 0
                SDeform(iI).dFy = uLAP{1}.*dDeformRes;
                SDeform(iI).dFx = uLAP{2}.*dDeformRes;
                SDeform(iI).dFz = uLAP{3}.*dDeformRes;
                
                SDeform(iI).dBy = -uLAP{1}.*dDeformRes;
                SDeform(iI).dBx = -uLAP{2}.*dDeformRes;
                SDeform(iI).dBz = -uLAP{3}.*dDeformRes;
            else
                SDeform(iI).dFy = uLAP{1};
                SDeform(iI).dFx = uLAP{2};
                SDeform(iI).dFz = uLAP{3};
                
                SDeform(iI).dBy = -uLAP{1};
                SDeform(iI).dBx = -uLAP{2};
                SDeform(iI).dBz = -uLAP{3};
            end          
        end
        
        
%         % backwards deformation field from registration
    %     if strcmp(SOptions.Algorithm, 'fast')
    %         uB = CG_MultiScale_LAP3D_Fast(Imoving, Ifixed, MaxLevel, MinLevel);
    %     else
    %         % slower
    %         uB = CG_MultiScale_LAP3D(Imoving, Ifixed, MaxLevel, MinLevel);
    %     end
    %     
    %     if dDeformRes ~= 0
    %         SDeform(iI).dBy = uB{1}.*dDeformRes;
    %         SDeform(iI).dBx = uB{2}.*dDeformRes;
    %         SDeform(iI).dBz = uB{3}.*dDeformRes;
    %     else
    %         SDeform(iI).dBy = uB{1};
    %         SDeform(iI).dBx = uB{2};
    %         SDeform(iI).dBz = uB{3};
    %     end    

%         % backwards def field with b2f function
%         [SDeform(iI).dBx, SDeform(iI).dBy, SDeform(iI).dBz] = backwards2forwards(SDeform(iI).dFx, SDeform(iI).dFy, SDeform(iI).dFz);

        st= st+1;fwaitbar(st/steps,h);
    end
end
try close(h); catch; end;
