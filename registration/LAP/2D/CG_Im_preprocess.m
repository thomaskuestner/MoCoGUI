function [I1_out, I2_out] = CG_Im_preprocess(I1_in, I2_in, process_type)
% function performs preprocessing on the input images. The type of
% preprocessing is determined by process_type. The current options are:
%               1: texture-structure decomposition used in Sun et al
%               2: high-pass filtering by removing local averages

switch process_type
    case 1
        % Texture-structure decomposition used in Sun et al
        
        % parameters
        theta   = 1/8; 
        nIters  = 100;           % iteration number
        alp     = 0.95;  
        
        [I1_out, I2_out, ~, ~] = image_decomposition(I1_in, I2_in, theta, nIters, alp);
        
    case 2
        % high pass filtering:
               
        % Find local average of the images:
        h = ones(5);
%         f = daubcqf(4,'min');
        I1_ave = imfilter(I1_in,h./sum(h(:)),'symmetric');
        I2_ave = imfilter(I2_in,h./sum(h(:)),'symmetric');

        % Remove local average:
        I1_out = I1_in - I1_ave;
        I2_out = I2_in - I2_ave;
        
        
%         [y_2,~] = mdwt(I2_in- mean(I2_in(:)),f,2);
%         y_2(1:97,1:146) = y_2(1:97,1:146) - imfilter(y_2(1:97,1:146), h./sum(h(:)),'symmetric');
%         [y_1,~] = mdwt(I1_in- mean(I1_in(:)),f,2);
%         y_1(1:97,1:146) = y_1(1:97,1:146) - imfilter(y_1(1:97,1:146), h./sum(h(:)),'symmetric');
%         [I2_out, ~] = midwt(y_2,f,2);
%         [I1_out, ~] = midwt(y_1,f,2);

%         % combine images into tensor:
%         im(:,:,1) = I1_in;
%         im(:,:,2) = I2_in;
% 
%         % Rescale the input image to [-1 1]
%         im = scale_image(im, 0,255);
%         
%         % Find local average of the images:
%         I1_ave = imfilter(im(:,:,1),h./sum(h(:)),'symmetric');
%         I2_ave = imfilter(im(:,:,2),h./sum(h(:)),'symmetric');
% 
%         % Remove local average:
%         I1_out = I1_in - I1_ave;
%         I2_out = I2_in - I2_ave;

    case 3
            % Gaussian filtering
            
            % paramaters:
%             sd = 0.3;
%             window_size = 5;
%             
%             % generate filter
%             f = fspecial('gaussian', [window_size window_size], sd); % std = 1 is better
            
            f = ([1, 4, 6, 4, 1]./16).'*([1, 4, 6, 4, 1]./16);
            
            % Tensor images:
            im(:,:,1) = I1_in;
            im(:,:,2) = I2_in;
            
            % Subtract filter versions from the original images:
            images  = imfilter(im, f, 'symmetric');
%             images  = scale_image(images, 0, 255);
            
            I1_out = images(:,:,1);
            I2_out = images(:,:,2);
    
    case 4
            % Genealised Laplacian            
%             % paramaters:
%             alp = 0;
%             
%             % generate filter
%             f = fspecial('laplacian', alp);
            f = [0.5, 1, 0.5; 1, -6, 1; 0.5, 1, 0.5];

            % Tensor images:
            im(:,:,1) = I1_in;
            im(:,:,2) = I2_in;
            
            % Subtract filter versions from the original images:
            images  = imfilter(im,f,'symmetric');            
            images  = scale_image(images, 0, 255);
            
            I1_out = images(:,:,1);
            I2_out = images(:,:,2);
            
    case 5
        % Genealised Laplacian            
            % paramaters:
            alp = 0.95;
            
%             % generate filter
%             f = fspecial('laplacian', alp);
            f = [0.5, 1, 0.5; 1, -6, 1; 0.5, 1, 0.5];

%             Tensor images:
            im(:,:,1) = I1_in;
            im(:,:,2) = I2_in;
            
            % Subtract filter versions from the original images:
            Ims  = imfilter(im,f,'symmetric');   
%             Ims  = scale_image(Ims, 0, 255);

            I1_out = Ims(:,:,1) + (1-alp).*(I1_in -  Ims(:,:,1));
            I2_out = Ims(:,:,2) + (1-alp).*(I2_in -  Ims(:,:,2));
            
%             % mean:
%             imean = imfilter(images,ones(K,K),'symmetric')./(K*K);
%             
%             % sd
%             isd = imfilter(abs(images - imean),ones(K,K),'symmetric')./(K*K);
%             
%             
%             images(:,:,1) = (images(:,:,1))./max(images,[],3);
%             images(:,:,2) = (images(:,:,2))./max(images,[],3);
%             images  = scale_image(images, 0, 255);
            
%             images  = scale_image(Ims, 0, 255, -std(Ims(:))*1.5, std(Ims(:))*1.5);
%             
%             I1_out = images(:,:,1);
%             I2_out = images(:,:,2);
    case 6
        % more high pass filtering:
        
        h = ones(5);
        I1_hat = exp(log(I1_in) - imfilter(log(I1_in),h./sum(h(:)),'symmetric'));
%         I1_hat = I1_hat - imfilter(I1_hat,h./sum(h(:)),'symmetric');

        I2_hat = exp(log(I2_in) - imfilter(log(I2_in),h./sum(h(:)),'symmetric'));
%         I2_hat = I2_hat - imfilter(I2_hat,h./sum(h(:)),'symmetric');

%         f = [0.5, 1, 0.5; 1, -6, 1; 0.5, 1, 0.5];
% 
%         I1_hat = exp(log(I1_in) - imfilter(log(I1_in),h./sum(h(:)),'symmetric'));
% %         I1_hat = exp(imfilter(log(I1_in),f,'symmetric'));
%         I1_hat = imfilter(I1_hat,f,'symmetric');
%         I2_hat = exp(log(I2_in) - imfilter(log(I2_in),h./sum(h(:)),'symmetric'));
% %         I2_hat = exp(imfilter(log(I2_in),f,'symmetric'));
%         I2_hat = imfilter(I2_hat,f,'symmetric');
        
        images(:,:,1) = I1_hat;
        images(:,:,2) = I2_hat;
        images  = scale_image(images, 0, 255);
            
        I1_out = images(:,:,1);
        I2_out = images(:,:,2);
end
end
        
