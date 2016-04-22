function [J,offset]=interp(x,y,I,inttype, exttype)
% [J,offset]=interp(x,y,I,[inttype])
% Interpolate the image I at the (in general, noninteger) positions x and y
% (usual Matlab convention for images). x and y are 2D matrices with same
% size. The optional argument inttype can take the values
%       * 'nearestneighbor'
%       * 'bilinear' (default)
%       * 'keys'
%       * 'cubicspline'
%       * 'shiftedlinear'
% Assumes images in double precision

% 1. standard image grid
% [s1,s2]=size(I);x=1:s1;x=x'*ones(1,s2);y=1:s2;y=ones(s1,1)*y;
%
% 2. normalized [-0.5,0.5]x[-0.5,0.5] Cartesian grid
% x=-0.5:0.1:0.5;sx=length(x);y=-0.5:0.01:0.5;sy=length(y);x=ones(sy,1)*x;y=y'*ones(1,sx);
%
% Example of use
% [s1,s2]=size(I);x=((1:s2)-(s2+1)/2)/(s2-1);x=ones(s1,1)*x;y=((1:s1)-(s1+1)/2)/(s1-1);y=y'*ones(1,s2);
% z=x+i*y;z=1./(1-z);x=real(z);y=imag(z);[xx,yy]=xy2ij(x,y,size(I));[J,offset]=interp(xx,yy,I,'bilinear');

if nargin<4
    inttype='bilinear';
end

if nargin < 5,
    % Image extension to fill the unknown pixels
    exttype='symh';                     
    
    % Possible options are: 
    % 'zpd' (zero-padding), 'symh' (half-point symmetry), 
    %       'symw' (whole-point symmetry), 'ppd' (periodization)
end

[a,b]=size(I);

% Determination of the interpolation function
switch inttype
    
    case 'nearestneighbor'
        L1=-0.5;    % Support of the interpolation kernel
        L2=+0.5;
        phi=@nearest;
        
    case 'bilinear'
        L1=-1;      
        L2=+1;
        phi=@linspline;
        
    case 'keys'
        L1=-2;      
        L2=+2;
        phi=@keys;
        
    case 'cubicspline'
        L1=-2;      
        L2=+2;
        phi=@cubicspline;
        
    case 'shiftedlinear'
        tau=1/2*(1-sqrt(3)/3);
        L1=floor(-1+tau);      
        L2=ceil(1+tau);
        phi=@(x)linspline(x-tau);        
end

% Minimum and maximum row index needed in the interpolation formula
k0=floor(min(x(:))-L2+1);
k1=floor(max(x(:))-L1);
l0=floor(min(y(:))-L2+1);
l1=floor(max(y(:))-L1);

offset=[1-k0 1-l0];

% Smallest box enclosing the image and the (x,y) positions
kk0=min(k0,1);
kk1=max(k1,a);
ll0=min(l0,1);
ll1=max(l1,b);

% Indices used in the interpolation formula
k=floor(x-L2+1);
l=floor(y-L2+1);

I0=ext(I,exttype,[1-kk0 kk1-a 1-ll0 ll1-b]);
I0=I0(1-kk0+(k0:k1),1-ll0+(l0:l1));
[a0,b0]=size(I0);

% Prefiltering when needed
switch inttype
    
    case 'cubicspline'
        z0=-2+sqrt(3);
        A=(1-z0)/(1+z0);
        % along columns first
        I0=A*(filtering(1,[1 -z0],I0,'causal')+filtering(1,[1 -z0],I0,'anticausal')-I0);
        % then lines
        I0=I0.';
        I0=A*(filtering(1,[1 -z0],I0,'causal')+filtering(1,[1 -z0],I0,'anticausal')-I0);
        if isreal(I)
            I0=real(I0.');
        else
            I0=I0.';
        end
        
    case 'shiftedlinear'
        z0=tau/(1-tau);
        % along columns first
        I0=1/(1-tau)*filtering(1,[1 z0],I0,'causal');
        % then lines
        I0=I0.';
        I0=1/(1-tau)*filtering(1,[1 z0],I0,'causal');
        I0=I0.';
end

% Kernel-based interpolation formula

J=zeros(size(x));
for dk=0:(L2-L1-1)
    for dl=0:(L2-L1-1)
        ind=k+dk+offset(1)+a0*(l+dl+offset(2)-1);       % matrices are ordered along columns in Matlab
        J=J+phi(x-k-dk).*phi(y-l-dl).*I0(ind);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Embedded functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function J=ext(I,exttype,extsize)
[a,b]=size(I);
newa=a+extsize(1)+extsize(2);
newb=b+extsize(3)+extsize(4);

if extsize(1)>extsize(2)
    J=wextend(2,exttype,I,extsize(1),'bn');
    J=J(1:newa,:);
else
    J=wextend(2,exttype,I,extsize(2),'bn');
    J=J(end+(1:newa)-newa,:);
end

if extsize(3)>extsize(4)
    J=wextend(2,exttype,J,extsize(3),'nb');
    J=J(:,1:newb);
else
    J=wextend(2,exttype,J,extsize(4),'nb');
    J=J(:,end+(1:newb)-newb);
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u=nearest(x)
u=(x<0.5&x>=-0.5);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u=linspline(x)
u=1-abs(x);
u=u.*(u>=0);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u=keys(x)
a=-1/2;
x=abs(x);
x2=x.*x;
x3=x.*x2;
u=((a+2)*x3-(a+3)*x2+1).*(x<=1)+...
    (a*x3-5*a*x2+8*a*x-4*a).*(x>1&x<=2);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u=cubicspline(x)
xi=2-abs(x);
xi2=xi.*xi;
xi3=xi.*xi2;
u=(2/3-2*xi+2*xi2-1/2*xi3).*(xi>1&xi<=2)+...
    (1/6*xi3).*(xi>0&xi<=1);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function J=filtering(b,a,I,type)
switch type
    case 'causal'
        J=bsxfun(@plus,filter(b,a,bsxfun(@minus,I,I(1,:))),I(1,:)*sum(b)/sum(a));
    case 'anticausal'
        I=I(end:-1:1,:);
        J=bsxfun(@plus,filter(b,a,bsxfun(@minus,I,I(1,:))),I(1,:)*sum(b)/sum(a));
        J=J(end:-1:1,:);
end
return

