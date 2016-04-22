function [uest,coeffs,err]=optiflowFilter2(I1,I2,basis, K1,I10,I20)

[s1,s2]=size(I1);

if length(size(basis))<=2
    [M,N]=size(basis);
    K0=round(sqrt(M));
    if M~=K0*K0
        error('Could not identify the size of the basis filters.')
    end
    basis=reshape(basis,[K0,K0,N]);
end
[K0,L0,N]=size(basis);
K=(K0-1)/2;
L=(L0-1)/2;
if round(K)~=K|round(L)~=L
    error('Basis filters are not centered.')
end

calculation=2;  	% formula used to calculate velocity from filter
if nargin == 3,
    K1=2*1*K+1;         % block size first dimension
end
K2=K1;              % block size second dimension

% display(['Velocity estimation with ' num2str(N) ' different ' num2str(2*K+1) 'x' num2str(2*K+1) ' filters'])
% display(['Block size = ' num2str(K1) 'x' num2str(K2) ' pixels'])
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% filtering with the basis %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
II=[];
            
for n=1:N
    B1=basis(:,:,n);
    B2=B1(end:-1:1,end:-1:1);
    if nargin>4&n==1
        J=imfilter(I20,B2,'symmetric')-imfilter(I10,B1,'symmetric');
    else
        J=imfilter(I2,B2,'symmetric')-imfilter(I1,B1,'symmetric');
    end
    II=[II J(:)];
end
J=II(:,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% matrices needed %%%%%%%%%
%%%%%% in the linear system %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
A=zeros(N-1,N-1,s1*s2);
b=zeros(N-1,s1*s2);
for k=1:N-1
    for l=k:N-1
        A(k,l,:)=average(II(:,k+1).*II(:,l+1));
        A(l,k,:)=A(k,l,:);
    end
    b(k,:)=-average(II(:,k+1).*J);
end

coeffs=zeros(s1*s2,N-1);
% for k=1:s1*s2
%     coeffs(k,:)=(A(:,:,k)\b(:,k))'; 
% end


for k=1:(N-1)
    for l=(k+1):(N-1)
        c=A(l,k,:)./A(k,k,:);
        for m=(k+1):(N-1)
            A(l,m,:)=A(l,m,:)-c.*A(k,m,:);
        end
        A(l,k,:)=0;
        b(l,:)=b(l,:)-shiftdim(c,1).*b(k,:);
    end
end
for k=(N-1):-1:1
    coeffs(:,k)=shiftdim(b(k,:));
    for m=(k+1):(N-1)
        coeffs(:,k)=coeffs(:,k)-shiftdim(A(k,m,:)).*coeffs(:,m);
    end
    coeffs(:,k)=coeffs(:,k)./shiftdim(A(k,k,:));
end


coeffs=[ones(s1*s2,1) coeffs];

% Exclusion of the boundaries
p=bndindex(I1,[K+(K1-1)/2,K+(K2-1)/2]);
coeffs(p,:)=NaN;

if nargout>2
    err=zeros(s1,s2);
    err(:)=sum(coeffs.*II,2);
end
            
k0=(-K:K)'*ones(1,2*K+1);l0=k0';
basis=reshape(basis,[(2*K+1)^2 N]);
switch calculation
    case 1
        z1=zeros(s1,s2);
        z2=zeros(s1,s2);
        H1(K+1+(2*K+1)*K)=1;
        for n=1:N
            z1(:)=z1(:)+(basis(:,n)'*exp(-2*i*pi*k0(:)/(2*K+1)))*coeffs(:,n);
            z2(:)=z2(:)+(basis(:,n)'*exp(-2*i*pi*l0(:)/(2*K+1)))*coeffs(:,n);
        end
        uest=(2*K+1)/2/pi*(angle(z1./conj(z1))+i*angle(z2./conj(z2)));
        
    case 2
        u1=zeros(s1,s2);
        u11=zeros(s1,s2);
        u2=zeros(s1,s2);
        u22=zeros(s1,s2);
        for n=1:N
            u1(:)=u1(:)-(basis(:,n)'*k0(:))*coeffs(:,n);
            u11(:)=u11(:)+sum(basis(:,n))*coeffs(:,n);
            u2(:)=u2(:)-(basis(:,n)'*l0(:))*coeffs(:,n);
            u22(:)=u22(:)+sum(basis(:,n))*coeffs(:,n);
        end
        uest=2*(u1./(u11)+i*u2./(u22));
        uest(find(abs(uest)>K))=NaN;

    case 3
        theta=2*pi/10;
        z1=zeros(s1,s2);
        z2=zeros(s1,s2);
        H1(K+1+(2*K+1)*K)=1;
        for n=1:N
            z1(:)=z1(:)+(basis(:,n)'*exp(-i*theta*k0(:)))*coeffs(:,n);
            z2(:)=z2(:)+(basis(:,n)'*exp(-i*theta*l0(:)))*coeffs(:,n);
        end
        uest=1/theta*(angle(z1./conj(z1))+i*angle(z2./conj(z2)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Embedded functions %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function J=average(I)
K1=evalin('caller','K1');
K2=evalin('caller','K2');
s1=evalin('caller','s1');
s2=evalin('caller','s2');
J=imfilter(reshape(I,s1,s2),ones(K1,K2),'symmetric')/(K1*K2);
J=J(:);
return

function p=bndindex(I,K)
if length(K)==1
    K=[K,K];
end
[M,N]=size(I);
m=(1:M)'*ones(1,N);
n=ones(M,1)*(1:N);
p=find(~(m>=K(1)+1&m<=M-K(1)&n>=K(2)+1&n<=N-K(2)));
return

function [a,b]=affine(I1,I2)
K=evalin('caller','K');
K1=evalin('caller','K1');
K2=evalin('caller','K2');
s1=evalin('caller','s1');
s2=evalin('caller','s2');
u1=average(I1);
u2=average(I2);
u11=average(I1.*I1);
u12=average(I1.*I2);
u22=average(I2.*I2);
d=u22-u2.*u2;

a=(u12-u1.*u2)./d;
b=(-u12.*u2+u22.*u1)./d;
return