function J=imshift(I,u, exttype)

if nargin == 2,
    % Image extension to fill the unknown pixels
    exttype='symh';                     
    
    % Possible options are: 
    % 'zpd' (zero-padding), 'symh' (half-point symmetry), 
    %       'symw' (whole-point symmetry), 'ppd' (periodization)
end

[M,N]=size(I);
integer_shift=(abs(round(u)-u)<=1e-6);

% if integer_shift
%     u=round(u);
%     xshift=abs(real(u));
%     yshift=abs(imag(u));tic
%     switch sign(real(u))
%         case -1
%             I1=[I;I(end-(1:xshift),:)];
%         case 0
%             I1=I;
%         case 1
%             I1=[I(1+(xshift:-1:1),:);I];
%     end
%     switch sign(imag(u))
%         case -1
%             I1=[I1 I1(:,end-(1:yshift))];
%         case 0
%         case 1
%             I1=[I1(:,1+(yshift:-1:1)) I1];
%     end
% 
%     J=I1((1:M)+max(0,-real(u)),(1:N)+max(0,-imag(u)));
% else
    x=(1:M)'*ones(1,N);
    y=ones(M,1)*(1:N);
    J=interp(x-real(u),y-imag(u),I,'cubicspline', exttype);
%     J=interp(x-real(u),y-imag(u),I,'shiftedlinear', exttype);
% end