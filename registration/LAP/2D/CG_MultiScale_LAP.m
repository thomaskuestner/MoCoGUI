function u_est = CG_MultiScale_LAP(I1, I2, Level_Num, Nf)
% The function implements a multi-scale framework for the LAP optical flow 
% algorithm. Instead of downsampling the images, the framework changes the 
% size of the all-pass filters used in the LAP algorithm. The filter basis
% used in the LAP algorithm spans the derivatives of a Gaussian filter.
% Note that this implementation is for greyscale images only.
% 
% Input parameters:
%       I1 and I2   -> Input images of size M by N (Greyscale)
%       Level_Num   -> Controls the maximum size of the filters,
%        i.e. filter size is [2^(Level_Num+1) + 1] by [2^(Level_Num+1) + 1]
%       Nf          -> Number of filters used in the basis (either 3 or 6)
%        i.e. Nf = 3 -> 1st derivatives of Gaussian, Nf = 6 -> 1st + 2nd
%
% Outputs:
%       u_est       -> Estimate of the optical flow (size M by N) 
%        
% Note that we use a complex representation for the optical flow. Therefore
% we convert a vector field u = [ux, uy]' into u = ux + j*uy (j = sqrt(-1))
%
% Date: 06/05/2015, Author: Chris Gilliam
% Code relates to: "Local All-Pass Filters for Optical Flow Estimation" by
% C. Gilliam and T. Blu, in IEEE ICASSP 2015.
%

% Obtain the dimensions of the images:
[L1, L2] = size(I1);

% Initial optical flow estimate
u_holder = zeros(L1,L2);

% define half support of filters (i.e. R)
amp_array = 2.^([Level_Num:-1:1,1]);

% Initialise local counter:
num_level = 0;

% Start estimating the optical flow
for l = 1:length(amp_array),
    num_level = num_level + 1;
    disp(['Level ', int2str(num_level), '/', int2str(Level_Num +1)]);
    
    % Define filter parameter R at each iteration
    amp_size = amp_array(l);

    % Load Filter Basis:
    Basis_Set = loadbasis(2,amp_size);
    Basis_Set = Basis_Set(:,1:Nf);
    
    if l == 1,
        % No warping for first iteration:
        I1_shift = I1;
    else
        % Warp I1 closer to I2 using current optical flow estimate
        I1_shift = imshift(I1,u_holder);
    end
    
    % Local I2 variable:
    im2 = I2;

    % Using basis functions estimate optical flow on a local scale:
    [uest_Orig, ~,~] = optiflowFilter2(I1_shift, im2, Basis_Set);

    % Clean optical flow
    % Stage 1: Remove nan's in the optical flow using inpainting:
    [uest_Clean, ~] = cleanOF(uest_Orig, not(isnan(uest_Orig)));

    % Stage 2: Remove flow elements that corresponding to large warping
    % errors. 
    scale1 = round(2*amp_size);
    scale = round(4*amp_size);
    uest_Clean = imfilter(uest_Clean, ones(scale1,scale1)./(scale1^2),'symmetric');
    uest_Clean = imfilter(uest_Clean, ones(scale,scale)./(scale^2),'symmetric');

    % Rescale optical flow and add to estimate
    u_holder = u_holder + uest_Clean;

    % Refinement of OF at highest level based on errors in shifted 
    %           images 
    if amp_size <= 2,
        u_holder = Cleaning_Procedure(u_holder);
    end

end 
        
u_est = u_holder;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Embedded functions %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function u_out = Cleaning_Procedure(u_in)
% function cleans the estimate of the optical flow. The cleaning process
% comprises a robust smoothing stage using two Median filters (of different
% sizes)
%

% Define size of median filters
B1 = 11;
B2 = 3;
        
% Two part median filtering:
% fine scale
u_out = medfilt2(real(u_in),[B2,B2], 'symmetric') + 1i.*medfilt2(imag(u_in),[B2,B2], 'symmetric');
% coarse scale
u_out = medfilt2(real(u_out),[B1, B1], 'symmetric') + 1i.*medfilt2(imag(u_out),[B1, B1], 'symmetric');
end

