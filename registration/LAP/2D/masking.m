function mask=masking(u,crit, p)
if nargin == 2,
    p=0.5;              % percentage of expected 'good' values
end
good=round(p*length(u(:)));

% Identification of 1-p percent of the image to be corrupted
switch crit
    case 'Laplacian'
        v=imfilter(u,[0 1 0;1 -4 1;0 1 0],'symmetric');
        [~,j]=sort(abs(v(:)));
        mask=zeros(size(u));
        mask(j(1:good))=1;
        mask=mask.*(~isnan(u));
    case 'genLaplacian'
%         K=4;
%         v=u-imfilter(u,ones(2*K+1,2*K+1)/(2*K+1)^2,'symmetric');
        
        f = [0.5, 1, 0.5; 1, -6, 1; 0.5, 1, 0.5];
        v=u-imfilter(u,f,'symmetric');
        [~,j]=sort(abs(v(:)));
        mask=zeros(size(u));
        mask(j(1:good))=1;
        mask=mask.*(~isnan(u));
    case 'Gradient'
        v1=imfilter(u,[1 2 1;0 0 0;-1 -2 -1],'symmetric');
        v2=imfilter(u,[1 2 1;0 0 0;-1 -2 -1],'symmetric');
        [~,j]=sort(abs(v1(:)).^2+abs(v2(:)).^2);
        mask=zeros(size(u));
        mask(j(1:good))=1;
        mask=mask.*(~isnan(u));
    case 'Curl'
        v1=imfilter(imag(u),[1 2 1;0 0 0;-1 -2 -1],'symmetric');
        v2=imfilter(real(u),[1 2 1;0 0 0;-1 -2 -1],'symmetric');
        [~,j]=sort(abs(v1(:)-v2(:)).^2);
        mask=zeros(size(u));
        mask(j(1:good))=1;
        mask=mask.*(~isnan(u));
        
    case 'Median'
        B1 = 11;
        
        % fine scale:
        v = medfilt2(real(u),[3,3], 'symmetric') + 1i.*medfilt2(imag(u),[3,3], 'symmetric');
        % coarse scale:
        v = medfilt2(real(v),[B1, B1], 'symmetric') + 1i.*medfilt2(imag(v),[B1, B1], 'symmetric');
        
        % sort difference:
        [~,j]=sort(abs(v(:) - u(:)));
        mask=zeros(size(u));
        mask(j(1:good))=1;
        mask=mask.*(~isnan(u));
        
    case 'Max'
        [~,j]=sort(abs(u(:)));
        mask=zeros(size(u));
        mask(j(1:good))=1;
        mask=mask.*(~isnan(u));
        
end

