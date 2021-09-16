function u = fDeformToU( SDeform, iDir )
% SDeform struct to flow matrix u
% Fixed --- forward field --> Moving
% Moving --- backward field --> Fixed

if(nargin < 2)
    % iDir == 1 (forward), iDir == 2 (backward)
    iDir = 1;
end

if(~isempty(SDeform(1).dFx))
    iSize = size(SDeform(1).dFx);
elseif(~isempty(SDeform(1).dBx))
    iSize = size(SDeform(1).dBx);
elseif(~isempty(SDeform(2).dFx))
    iSize = size(SDeform(2).dFx);
else
    iSize = size(SDeform(2).dBx);
end
   
% spatial size x time x flowDir (3D)
u = zeros([iSize, size(SDeform,2), 3]);

for iTime=1:size(SDeform,2)
    if iDir == 1
        if(isempty(SDeform(iTime).dFx))
            continue;
        end
        u(:,:,:,iTime,1) = SDeform(iTime).dFx;
        u(:,:,:,iTime,2) = SDeform(iTime).dFy;
        u(:,:,:,iTime,3) = SDeform(iTime).dFz;
    else
        if(isempty(SDeform(iTime).dBx))
            continue;
        end
        u(:,:,:,iTime,1) = SDeform(iTime).dBx;
        u(:,:,:,iTime,2) = SDeform(iTime).dBy;
        u(:,:,:,iTime,3) = SDeform(iTime).dBz;
    end    
end

