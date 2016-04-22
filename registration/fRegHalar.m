function [SDeform, dImgReg] = fRegHalar(dFix, dMove, sParafile, iDim, SGeo, sResultPath, sRegPath)
% function to run halar registration method
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

h = fwaitbar(0,'Running Registration. Please wait!'); st=0;

%% prepare images
nDimImg = ndims(dFix);
iNGates = size(dMove, ndims(dMove)) + 1;
if(nDimImg == 2)
    iN3D = 1;
    dImgReg = zeros(size(dFix,1), size(dFix,2), 1, iNGates); % always 4D array: x-y-z-t
    SImg.size = [size(dFix, 1), size(dFix, 2)];
    SImg.orientation = eye(2);
    SImg.origin = [0; 0];
elseif(nDimImg == 3)
    if(iDim == 1) % 2D reg
        iN3D = size(dFix,3);
    else % 3D reg
        iN3D = 1;
    end
    dImgReg = zeros(size(dFix,1), size(dFix,2), size(dFix,3), iNGates);
    SImg.size = [size(dFix, 1), size(dFix, 2), size(dFix, 3)];
    SImg.orientation = eye(3);
    SImg.origin = [0; 0; 0];
end
dImgReg(:,:,:,1) = dFix;
RORes          = SGeo.cVoxelsize{1}(1);
PERes          = SGeo.cVoxelsize{1}(2);
SLRes          = SGeo.cVoxelsize{1}(3);
steps = iN3D * (iNGates-1); % for waitbar

cImgFix = cell(iN3D,1);
cImgMove = cell(iN3D,(iNGates-1));

dImgFix = uint16(dFix./max(dFix(:)).*(2.^16 - 1));
dImgMove = uint16(dMove./max(dMove(:)).*(2.^16 - 1));
    
    
for iJ=1:iN3D
    if(iDim == 1) % 2D reg
        iImg = dImgFix(:,:,iJ);
        cImgFix{iJ} = sprintf('imgFix%02u.gipl',iJ);
        giplwrite([sResultPath,filesep,cImgFix{iJ}],iImg,'short_noscale');
    else % 3D
        iImg = dImgFix;
        cImgFix{iJ} = 'imgFix.gipl';
        giplwrite([sResultPath,filesep,cImgFix{iJ}],iImg,'short_noscale');
    end

    for iI=2:iNGates % moving images
        if(nDimImg == 2) % iDim = 1 forced -> just 2D reg
            iImg = dImgMove(:,:,iI-1);
            cImgMove{iJ,iI} = sprintf('imgMove_G%02u.gipl', iI);
        elseif(nDimImg == 3)
            if(iDim == 1) % 2D
                iImg = dImgMove(:,:,iJ,iI-1);
                cImgMove{iJ,iI} = sprintf('imgMove%02u_G%02u.gipl', iJ, iI);
            else % 3D
                iImg = dImgMove(:,:,:,iI-1);
                cImgMove{iJ,iI} = sprintf('imgMove_G%02u.gipl', iI);
            end
        end
        giplwrite([sResultPath,filesep,cImgMove{iJ,iI}],iImg,'short_noscale');
   end
end

%% start registration
sRegPath = [sRegPath,filesep,'registration',filesep,'halar'];

if(ispc), slreg = 'lreg.exe'; else slreg = 'lreg'; end;

% register 3D frames and warp 3D frame to end-exhale position
for iJ = 1:iN3D
    for iI = 2:iNGates    
        tic
        % register 3D frames (backward df)
        system([sRegPath,filesep,slreg,' -R ',sResultPath,filesep,cImgMove{iJ,iI},' ',sResultPath,filesep,cImgFix{iJ},' -parin ',sParafile,' -dofout ',sResultPath,filesep,'dof_01_to_',sprintf('%02u.dof',iI)]);

        % transform 3D frame and save deformation field
        system([sRegPath,filesep,slreg,' -T ',sResultPath,filesep,cImgMove{iJ,iI},' ',sResultPath,filesep,'Gate',sprintf('%02u_deformed.gipl ',iI),' -dofin ',sResultPath,filesep,'dof_01_to_',sprintf('%02u.dof',iI),' -mx ',sResultPath,filesep,sprintf('df_01_to_%02u_b_x.gipl',iI),' -my ',sResultPath,filesep,sprintf('df_01_to_%02u_b_y.gipl',iI),' -mz ',sResultPath,filesep,sprintf('df_01_to_%02u_b_z.gipl',iI)]);
        disp(['time for inversion ', num2str(toc)]) 

        % compute forward deformation field by inversion
        tic
        system([sRegPath,filesep,slreg,' -T ',sResultPath,filesep,cImgMove{iJ,iI},' ',sResultPath,filesep,'Gate_dummy.gipl -dofin ',sResultPath,filesep,'dof_01_to_',sprintf('%02u.dof',iI),' -mx ',sResultPath,filesep,sprintf('df_01_to_%02u_f_x.gipl',iI),' -my ',sResultPath,filesep,sprintf('df_01_to_%02u_f_y.gipl',iI),' -mz ',sResultPath,filesep,sprintf('df_01_to_%02u_f_z.gipl',iI)]);
        disp(['time for inversion ', num2str(toc)])

        % save the registered image
        dReg = giplread([sResultPath,filesep,sprintf('Gate%02u_deformed.gipl',iI)]);
        if(nDimImg == 2) % 2D image => 2D reg
            dImgReg(:,:,1,iI) = permute(dReg,[1 2 4 3]);
        elseif(nDimImg == 3)
            if(iDim == 1) % 2D reg
                dImgReg(:,:,iJ,iI) = dReg;
            else % 3D reg
                dImgReg(:,:,:,iI) = dReg;
            end
        end

        % save deformation fields
        % forward df: remaining (target) to end-exhale (source) frames
        % backward df: end-exhale (source) to remaining (target) frames
        % deformation fields from source to target frames 
        % [mm] (with .* SGeo.X) | [px] (without SGeo.X)
        if(iDim == 1) % 2D reg
            if(~exist('SDeform','var') || ~isfield(SDeform(iI),'dFy'))
                SDeform(iI).dFy = zeros(size(dImgReg,1),size(dImgReg,2),size(dImgReg,3));
                SDeform(iI).dFx = zeros(size(dImgReg,1),size(dImgReg,2),size(dImgReg,3));
                SDeform(iI).dFz = zeros(size(dImgReg,1),size(dImgReg,2),size(dImgReg,3));
                
                SDeform(iI).dBy = zeros(size(dImgReg,1),size(dImgReg,2),size(dImgReg,3));
                SDeform(iI).dBx = zeros(size(dImgReg,1),size(dImgReg,2),size(dImgReg,3));
                SDeform(iI).dBz = zeros(size(dImgReg,1),size(dImgReg,2),size(dImgReg,3));                
            end
            SDeform(iI).dFy(:,:,iJ) = giplread([sResultPath,filesep,sprintf('df_01_to_%02u_f_y.gipl',iI)]); % .* RORes;
            SDeform(iI).dFx(:,:,iJ) = giplread([sResultPath,filesep,sprintf('df_01_to_%02u_f_x.gipl',iI)]); % .* PERes;
            
            SDeform(iI).dBx = giplread([sResultPath,filesep,sprintf('df_01_to_%02u_b_x.gipl',iI)]); % .* PERes;
            SDeform(iI).dBy = giplread([sResultPath,filesep,sprintf('df_01_to_%02u_b_y.gipl',iI)]); % .* RORes;
        elseif(iDim == 2) % 3D reg
            SDeform(iI).dFy = giplread([sResultPath,filesep,sprintf('df_01_to_%02u_f_y.gipl',iI)]); % .* RORes;
            SDeform(iI).dFx = giplread([sResultPath,filesep,sprintf('df_01_to_%02u_f_x.gipl',iI)]); % .* PERes;
            SDeform(iI).dFz = giplread([sResultPath,filesep,sprintf('df_01_to_%02u_f_z.gipl',iI)]); % .* SLRes;
            
            SDeform(iI).dBx = giplread([sResultPath,filesep,sprintf('df_01_to_%02u_b_x.gipl',iI)]); % .* PERes;
            SDeform(iI).dBy = giplread([sResultPath,filesep,sprintf('df_01_to_%02u_b_y.gipl',iI)]); % .* RORes;
            SDeform(iI).dBz = giplread([sResultPath,filesep,sprintf('df_01_to_%02u_b_z.gipl',iI)]); % .* SLRes;
        end
        
        st= st+1;fwaitbar(st/steps,h);
    end
end
try close(h); catch; end;
