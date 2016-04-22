function basis=loadbasis(btype,K)
if nargin==0
    btype=1;
end

switch btype
    case 1
        basis=zeros(5,5,3);
        basis(:,:,1)=[0 0 1 0 0;0 15 30 15 0;1 30 0 30 1;0 15 30 15 0;0 0 1 0 0];
        basis(:,:,2)=[0 0 2 0 0;0 15 30 11 0;0 4 0 -4 0;0 -11 -30 -15 0;0 0 -2 0 0];
        basis(:,:,3)=basis(:,:,2)';
        basis=reshape(basis,25,3);
    case 2
        s=1+(K-2)/4;
        g=@(x)exp(-x.*x/2/s^2);
        
        K0=ceil(K);
        k=-K0:K0;k=k'*ones(1,2*K0+1);l=k';
        
        basis=zeros(2*K0+1,2*K0+1,6);
        basis(:,:,1)=g(k).*g(l);basis(:,:,1)=basis(:,:,1)/sum(sum(basis(:,:,1)));
        basis(:,:,2)=g(k).*l.*g(l);basis(:,:,2)=basis(:,:,2)/sum(sum(l.*basis(:,:,2)));
        basis(:,:,3)=basis(:,:,2)';
        basis(:,:,4)=g(k).*(k.*k+l.*l).*g(l)/s^2;
        basis(:,:,4)=basis(:,:,4)-basis(:,:,1)*sum(sum(basis(:,:,4)))/sum(sum(basis(:,:,1)));
        basis(:,:,5)=g(k).*(k.*k-l.*l).*g(l)/s^2;
        basis(:,:,6)=g(k).*(k.*l).*g(l)/s^2;
        basis=reshape(basis,(2*K0+1)^2,6);
        
    case 3
        if nargin==1
            K=2;
        end
        s=pi/2+(K-2)/4;
        g=@(x) (s/2)./(pi.*(x.^2 + (s/2)^2));
        gd=@(x) -(16.*x.*s)./(pi.*(4.*x.^2 + (s)^2).^2);
        
        K0=ceil(K);
        k=-K0:K0;k=k'*ones(1,2*K0+1);l=k';
        
        basis=zeros(2*K0+1,2*K0+1,3);
        basis(:,:,1)=g(k).*g(l);
        basis(:,:,2)=g(k).*gd(l);
        basis(:,:,3)=basis(:,:,2)';
%         basis(:,:,4)=finiteDiff(2*K0+1, 2) + finiteDiff(2*K0+1, 2).';
        basis=reshape(basis,(2*K0+1)^2,3);
end
end


% function D = finiteDiff(pts, n)
% % calculates finite difference coefficients for pts points of order n
% 
% N=2*pts-1;  p1=pts-1;
% 
% A=repmat((0:p1)',1,N);      B=repmat((-p1:p1),pts,1);
% 
% M=(B.^A)./gamma(A+1); 
% 
% rhs=zeros(pts,1);   rhs(n+1)=1;
% 
% W=zeros(pts,pts);
% 
% for k=1:pts
%     W(:,k)=M(:,(0:p1)+k)\rhs;
% end
% 
% W=W';   W(1:pts,:)=W(pts:-1:1,:);
% D= W;
% % D = ones(pts,1)*W(ceil(pts/2), :)/pts;
% 
% end