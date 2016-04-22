function bwe_def = fDeformLine(bwe_open, DefField, gateNo, sliceNo, SCT)
% deform line segment (LandmarkGUI) of original images
%
% input:
% bwe_open          edge image
% DefField          deformation field
% gateNo            current moving state
% sliceNo           current slice (z) position
% SCT               sag/cor/tra image orientation
%
% output:
% bwe_def           deformed/transformed edge image
%
% -------------------------------------------------------------------------
% (c) 2015: Thomas Kuestner, Verena Neumann
% -------------------------------------------------------------------------

warning('off','MATLAB:maxNumCompThreads:Deprecated');
% threshold for pixel
thrs = 0.5;

if strcmp(SCT,'cor')
    % get deformation field that corresponds to gate and slice number
    u(:,:,1) = DefField(gateNo).dFy(:,:,sliceNo);
    u(:,:,2) = DefField(gateNo).dFx(:,:,sliceNo);
    u(:,:,3) = DefField(gateNo).dFz(:,:,sliceNo);
       
    % apply deformation field to the line
%     bwe_def = movepixels(bwe_open,u,1);
    bwe_def = movepixels(bwe_open,u(:,:,1),u(:,:,2),u(:,:,3),1);
    bwe_def(bwe_def>thrs) = 1;
    bwe_def(bwe_def<=thrs) = 0;
    
elseif strcmp(SCT,'sag')
    
    % get deformation field that corresponds to gate and slice number
    dF1 = permute(DefField(gateNo).dFy(:,:,:),[1,3,2]);
    dF2 = permute(DefField(gateNo).dFz(:,:,:),[1,3,2]);
    dF3 = permute(DefField(gateNo).dFx(:,:,:),[1,3,2]);    
    u(:,:,1) = dF1(:,:,sliceNo);
    u(:,:,2) = dF2(:,:,sliceNo);
    u(:,:,3) = dF3(:,:,sliceNo);
    
    % apply deformation field to the line
%     bwe_def = movepixels(bwe_open,u,1);
    bwe_def = movepixels(bwe_open,u(:,:,1),u(:,:,2),u(:,:,3),1);
    bwe_def(bwe_def>thrs) = 1;
    bwe_def(bwe_def<=thrs) = 0;
    
elseif strcmp(SCT,'tra')
    % get deformation field that corresponds to gate and slice number
    dF1 = permute(DefField(gateNo).dFy(:,:,:),[3,2,1]);
    dF2 = permute(DefField(gateNo).dFz(:,:,:),[3,2,1]);
    dF3 = permute(DefField(gateNo).dFx(:,:,:),[3,2,1]);    
    u(:,:,1) = dF1(:,:,sliceNo);
    u(:,:,2) = dF2(:,:,sliceNo);
    u(:,:,3) = dF3(:,:,sliceNo);
    
    % apply deformation field to the line
%     bwe_def = movepixels(bwe_open,u,1);
    bwe_def = movepixels(bwe_open,u(:,:,1),u(:,:,2),u(:,:,3),1);
    bwe_def(bwe_def>thrs) = 1;
    bwe_def(bwe_def<=thrs) = 0;
else
    error('The orientation you have given is not implemented')
end
bwe_def(isnan(bwe_def)) = 0;