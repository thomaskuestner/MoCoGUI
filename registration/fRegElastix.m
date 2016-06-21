function [SDeform, dImgReg] = fRegElastix(dFix, dMove, sParafile, iDim, SGeo, sResultPath, sRegPath)
% function to run elastix registration method
%
% input:
% dFix          reference/fixed image (2D/3D)
% dMove         moving images (2D/3D): x-y-(z)-t
% sParafile     (string) path to parameter file
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

if ispc
    elastix = 'elastix.exe';
    transformix = 'transformix.exe';
else
    elastix = 'elastix';
    transformix = 'transformix';
end

if(~exist(sResultPath,'dir'))
    mkdir(sResultPath);
end
sCurrDir = pwd;
cd(sResultPath);

h = fwaitbar(0,'Running Registration. Please wait!'); st=0; 

%% prepare images
nDimImg = ndims(dFix);
if(nDimImg == 2)
    iNGates = size(dMove, 3) + 1;
    iN3D = 1;
    dImgReg = zeros(size(dFix,1), size(dFix,2), 1, iNGates); % always 4D array: x-y-z-t
    SImg.size = [size(dFix, 1), size(dFix, 2)];
    SImg.orientation = eye(2);
    SImg.origin = [0; 0];
elseif(nDimImg == 3)
    iNGates = size(dMove, 4) + 1;
    if(iDim == 1) % 2D reg
        iN3D = size(dFix,3);
        SImg.size = [size(dFix, 1), size(dFix, 2)];
        SImg.orientation = eye(2);
        SImg.origin = [0; 0];
    else % 3D reg
        iN3D = 1;
        SImg.size = [size(dFix, 1), size(dFix, 2), size(dFix, 3)];
        SImg.orientation = eye(3);
        SImg.origin = [0; 0; 0];
    end
    dImgReg = zeros(size(dFix,1), size(dFix,2), size(dFix,3), iNGates);    
end
dImgReg(:,:,:,1) = dFix;
steps = iN3D * (iNGates-1); % for waitbar

sPrecision = fReadText(sParafile,'FixedInternalImagePixelType');
sPrecision = sPrecision(31:end-2);
SImg.spacing = SGeo.cVoxelsize{1};

cImgFix = cell(iN3D,1);
cImgMove = cell(iN3D,(iNGates-1));
if(strcmp(sPrecision,'short'))
    dImgFix = uint16(dFix./max(dFix(:)).*(2.^16 - 1));
    dImgMove = uint16(dMove./max(dMove(:)).*(2.^16 - 1));
    sPrecWrite = 'int16';
elseif(strcmp(sPrecision,'float'))
    dImgFix = dFix;
    dImgMove = dMove;
    sPrecWrite = 'single';
end
    
for iJ=1:iN3D
    if(iDim == 1) % 2D
        SImg.data = dImgFix(:,:,iJ,1);
        cImgFix{iJ} = sprintf('imgFix%02u.mhd',iJ);
        write_mhd(cImgFix{iJ}, SImg, 'ElementType', sPrecWrite);
    else % 3D
        SImg.data = dImgFix;
        cImgFix{iJ} = 'imgFix.mhd';
        write_mhd(cImgFix{iJ}, SImg, 'ElementType', sPrecWrite);
    end

    for iI=2:iNGates % moving images
        if(nDimImg == 2) % iDim = 1 forced -> just 2D reg
            SImg.data = dImgMove(:,:,1,iI-1);
            cImgMove{iJ,iI} = sprintf('imgMove_G%02u.mhd', iI);
        elseif(nDimImg == 3)
            if(iDim == 1) % 2D
                SImg.data = dImgMove(:,:,iJ,iI-1);
                cImgMove{iJ,iI} = sprintf('imgMove%02u_G%02u.mhd', iJ, iI);
            else % 3D
                SImg.data = dImgMove(:,:,:,iI-1);
                cImgMove{iJ,iI} = sprintf('imgMove_G%02u.mhd', iI);
            end
        end
        SImg.spacing = SGeo.cVoxelsize{iI};
        write_mhd(cImgMove{iJ,iI}, SImg, 'ElementType', sPrecWrite);
   end
end


%% start registration
sRegPath = [sRegPath,filesep,'registration',filesep,'elastix'];

% write inverse parameter file
writeInvParamFile(sParafile);
[pathstr,name,ext] = fileparts(sParafile);
sInvParafile = [pathstr,filesep,name,'_inv',ext];


% -------------------------------------------------------------------------
% Do the registration and read back the data into SDeform (forward and
% backwards)
for iJ = 1:iN3D
    for iI = 2:iNGates
        % register images
        system([sRegPath, filesep, elastix,' -f ',cImgFix{iJ},' -m ', cImgMove{iJ,iI}, ' -out ',sResultPath,' -p ', sParafile]);

        % get deformation field from transformation
        %     system([sRegPath, filesep, 'transformix -in ', sprintf('Gate%02u.mhd', iI), ' -out ', sResultPath, ' -tp ',sElastixPath, filesep, 'TransformParameters.0.txt'])
        system([sRegPath, filesep, transformix,' -def all -out ', sResultPath, ' -tp ',sResultPath, filesep, 'TransformParameters.0.txt']);

        % read and save registered image and deformation field
        SRegResult = read_mhd([sResultPath, filesep, 'result.0.mhd']);
        if(nDimImg == 2) % 2D image => 2D reg
            dImgReg(:,:,1,iI) = permute(SRegResult.data,[1 2 4 3]);
        elseif(nDimImg == 3)
            if(iDim == 1) % 2D reg
                dImgReg(:,:,iJ,iI) = SRegResult.data;
            else % 3D reg
                dImgReg(:,:,:,iI) = SRegResult.data;
            end
        end

        SData = read_mhd([sResultPath, filesep, 'deformationField.mhd']);
        if(iDim == 1) % 2D reg
            if(~exist('SDeform','var') || length(SDeform) < iI || ~isfield(SDeform(iI),'dFy'))
                SDeform(iI).dFy = zeros(size(dImgReg,1),size(dImgReg,2),size(dImgReg,3));
                SDeform(iI).dFx = zeros(size(dImgReg,1),size(dImgReg,2),size(dImgReg,3));
                SDeform(iI).dFz = zeros(size(dImgReg,1),size(dImgReg,2),size(dImgReg,3));
            end
            SDeform(iI).dFy(:,:,iJ) = SData.datax./SGeo.cVoxelsize{1}(1); % in px | without ./dVoxelsize: in mm
            SDeform(iI).dFx(:,:,iJ) = SData.datay./SGeo.cVoxelsize{1}(2); % in px
        elseif(iDim == 2) % 3D reg
            SDeform(iI).dFy = SData.datax./SGeo.cVoxelsize{1}(1); % in px | without ./dVoxelsize: in mm
            SDeform(iI).dFx = SData.datay./SGeo.cVoxelsize{1}(2); % in px
            SDeform(iI).dFz = SData.dataz./SGeo.cVoxelsize{1}(3); % in px
        end

    %     % for debugging det-jacobian-field in EvalGUI compare to elastix' field
    %     system([sRegPath, filesep, transformix, ' -jac all -out ', sElastixPath, ' -tp ',sElastixPath, filesep, 'TransformParameters.0.txt']);
    %     SJac = read_mhd([sElastixPath, filesep, 'spatialJacobian.mhd']);
    %     spatialJac = SJac.data;
    %     imagine(spatialJac);


        % compute backwards deformation field

    %     % a) simple but time consuming b2f function
    %     [SDeform(iI).dBx, SDeform(iI).dBy, SDeform(iI).dBz] = backwards2forwards(SDeform(iI).dFx, SDeform(iI).dFy, SDeform(iI).dFz);

        % b) special elastix call
        % gate01 to gate01 with inverse parameter file and 'no initial
        % transform' for transformix call (see elastix manual chapter 6.1.6)
        system([sRegPath, filesep, elastix,' -f ',cImgFix{iJ},' -m ', cImgFix{iJ},' -out ', sResultPath,' -p ', sInvParafile,' -t0 ', sResultPath, filesep, 'TransformParameters.0.txt']);
        fReplaceText([sResultPath, filesep, 'TransformParameters.0.txt'], 'NoInitialTransform');
        system([sRegPath,filesep,transformix,' -def all -out ',sResultPath,' -tp ',sResultPath,filesep,'TransformParameters.0.txt']);

        SData = read_mhd([sResultPath, filesep,'deformationField.mhd']);
        if(iDim == 1) % 2D reg
            if(~exist('SDeform','var') || length(SDeform) < iI || ~isfield(SDeform(iI),'dBy'))
                SDeform(iI).dBy = zeros(size(dImgReg,1),size(dImgReg,2),size(dImgReg,3));
                SDeform(iI).dBx = zeros(size(dImgReg,1),size(dImgReg,2),size(dImgReg,3));
                SDeform(iI).dBz = zeros(size(dImgReg,1),size(dImgReg,2),size(dImgReg,3));
            end
            SDeform(iI).dBy(:,:,iJ) = SData.datax./SGeo.cVoxelsize{1}(1); % in px | without ./dVoxelsize: in mm
            SDeform(iI).dBx(:,:,iJ) = SData.datay./SGeo.cVoxelsize{1}(2); % in px
        elseif(iDim == 2) % 3D reg
            SDeform(iI).dBy = SData.datax./SGeo.cVoxelsize{1}(1); % in px
            SDeform(iI).dBx = SData.datay./SGeo.cVoxelsize{1}(2); % in px
            SDeform(iI).dBz = SData.dataz./SGeo.cVoxelsize{1}(3); % in px
        end

        try st= st+1;fwaitbar(st/steps,h); catch; end;
    end
end
try close(h); catch; end;

cd(sCurrDir);

