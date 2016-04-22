function [u0,mask]=cleanOF(u,mask)
if nargin<2
    mask=masking(u,'Laplacian',0.99);
end
mask=double(mask);

badidx=find(mask==0);
u(badidx)=0;

% Diffusion
window=3;
h=[0 1 0;1 0 1;0 1 0];
h=conv2(h,h);
u0=u;
encore=1;
mask0=mask;
while encore
    un=imfilter(mask0,h,'symmetric');
    f=imfilter(mask0.*u0,h,'symmetric');
    n=find(mask==0&un~=0);
    u0(n)=f(n)./un(n);
    mask0=double(un~=0);
    encore=double(min(un(:))==0);
end

% % Computation of averages
% window=2;
% u0=u;
% encore=1;
% todo=ones(1,length(badidx));
% while encore
%     window=2*window-1;
%     h=ones(window,window);
%     un=imfilter(mask,h,'symmetric');
%     f=imfilter(mask.*u,h,'symmetric');
% 
%     todonow=todo.*(un(badidx)>=10)';
%     k0=find(todonow==1);
%     todo=todo.*(1-todonow);
%     encore=(sum(todo)>=1);
%     
%     if ~isempty(k0)
%         u0(badidx(k0))=f(badidx(k0))./un(badidx(k0));
%     end
% end 

% [M,N]=size(u);
% x=(1:M)'*ones(1,N);
% y=ones(M,1)*(1:N);
% good=find(mask==1);
% fx=TriScatteredInterp(x(good),y(good),real(u(good)));
% fy=TriScatteredInterp(x(good),y(good),imag(u(good)));
% u0(badidx)=fx(x(badidx),y(badidx))+i*fy(x(badidx),y(badidx));

%u0=zeroLap(u,mask);