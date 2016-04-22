function A = fOwnDiff(dFx2,dim,fac)
% compute difference along dimension

A = 0*dFx2;
if dim == 1
    A(1:end-1,:,:) = diff(dFx2,1,dim)/fac;
    A(end,:,:) = A(end-1,:,:);
elseif dim == 2 
    A(:,1:end-1,:) = diff(dFx2,1,dim)/fac;
    A(end,:,:) = A(end-1,:,:);
elseif dim == 3
    A(:,:,1:end-1) = diff(dFx2,1,dim)/fac;
    A(end,:,:) = A(end-1,:,:);  
else
    error('fOwnDiff is not implemented for this dimension')
end