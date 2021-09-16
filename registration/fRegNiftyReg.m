function [SDeform, dImgReg] = fRegNiftyReg(dFix, dMove, sPara, iDim, SGeo, sResultPath, sRegPath)
% function to run elastix registration method
%
% input:
% dFix          reference/fixed image (2D/3D)
% dMove         moving images (2D/3D): x-y-(z)-t
% sPara         (string) parameters to be used
% iDim          registration dimensionality (1=2D, 2=3D)
% SGeo          (struct) geometric image information (voxelsize, orientation, ...)
% sResultPath   (string) path to the folder where results will be saved
% sRegPath      (string) path to registration executables
%
% output:
% SDeform       (struct) deformation field (time x 1) with forward (F) and backward (B) fields
% dImgReg       transformed image from dMove towards dFix according to SDeform
%
% -------------------------------------------------------------------------
% (c) 2015: Thomas Kuestner, Verena Neumann
% -------------------------------------------------------------------------

%%

if(isempty(sPara))
%     sPara = ' --nmi -be 0.0001 -jl 0.8 -ln 3  -sx 8'; % CMRA Teresa
%     sPara = ' --nmi -be 0.0005 -sx 14'; % CMRA Aurelien
    sPara = ' --lncc 3.0 -be 0.0005 -ln 5 -sx 14'; % Abdomen Gastao
end

if(~exist(sResultPath,'dir'))
    mkdir(sResultPath);
end
sCurrDir = pwd;
cd(sResultPath);

if(~ispc)
%    setenv('LD_LIBRARY_PATH', '/usr/local/MATLAB/R2017a/sys/os/glnxa64');
    setenv('LD_LIBRARY_PATH', '');
end

h = fwaitbar(0,'Running Registration. Please wait!'); st=0; 

%% prepare images
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
dImgReg(:,:,:,1) = dFix;
steps = iN3D * (iNGates-1); % for waitbar

% scale to [0,1]
dFix = (dFix - min(dFix(:)))./(max(dFix(:)) - min(dFix(:)));
dMove = (dMove - min(dMove(:)))./(max(dMove(:)) - min(dMove(:)));


%% start registration
sRegPath = [sRegPath,filesep,'registration',filesep,'NIFTYReg',filesep,'bin',filesep];

% -------------------------------------------------------------------------
% Do the registration and read back the data into SDeform (forward and
% backwards)
for iJ = 1:iN3D
    for iI = 2:iNGates
        % register images
        
        if(iDim == 1) % 2D reg
            [dImgReg(:,:,iJ,iI), displacement_field] = nifty_reg(dMove(:,:,iJ,iI-1),dFix(:,:,iJ),sPara,'./',true,sRegPath);

            displacement_field = squeeze(displacement_field);
            if(~exist('SDeform','var') || length(SDeform) < iI || ~isfield(SDeform(iI),'dFy'))
                SDeform(iI).dFy = zeros(size(dImgReg,1),size(dImgReg,2),size(dImgReg,3));
                SDeform(iI).dFx = zeros(size(dImgReg,1),size(dImgReg,2),size(dImgReg,3));
                SDeform(iI).dFz = zeros(size(dImgReg,1),size(dImgReg,2),size(dImgReg,3));
            end
            SDeform(iI).dFy(:,:,iJ) = -displacement_field(:,:,1); % in px | without ./dVoxelsize: in mm
            SDeform(iI).dFx(:,:,iJ) = -displacement_field(:,:,2); % in px
            
            if(~exist('SDeform','var') || length(SDeform) < iI || ~isfield(SDeform(iI),'dBy'))
                SDeform(iI).dBy = zeros(size(dImgReg,1),size(dImgReg,2),size(dImgReg,3));
                SDeform(iI).dBx = zeros(size(dImgReg,1),size(dImgReg,2),size(dImgReg,3));
                SDeform(iI).dBz = zeros(size(dImgReg,1),size(dImgReg,2),size(dImgReg,3));
            end
            SDeform(iI).dBy(:,:,iJ) = displacement_field(:,:,1); % in px | without ./dVoxelsize: in mm
            SDeform(iI).dBx(:,:,iJ) = displacement_field(:,:,2);
        elseif(iDim == 2) % 3D reg
            [dImgReg(:,:,:,iI), displacement_field] = nifty_reg(dMove(:,:,:,iI-1),dFix,sPara,'./',true,sRegPath);

            displacement_field = squeeze(displacement_field);
            
            SDeform(iI).dFy = -displacement_field(:,:,:,1); % in px | without ./dVoxelsize: in mm
            SDeform(iI).dFx = -displacement_field(:,:,:,2); % in px
            SDeform(iI).dFz = -displacement_field(:,:,:,3); % in px
            
            SDeform(iI).dBy = displacement_field(:,:,:,1); % in px
            SDeform(iI).dBx = displacement_field(:,:,:,2); % in px
            SDeform(iI).dBz = displacement_field(:,:,:,3); % in px
        end

        try st= st+1;fwaitbar(st/steps,h); catch; end;
    end
end
try close(h); catch; end;

cd(sCurrDir);

