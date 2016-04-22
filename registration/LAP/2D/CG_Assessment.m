function Error_Data = CG_Assessment(u_GT, OF_index, u_est)

[M,N] = size(u_GT);

% ******* Absolute Error Analysis *******
% Calculate the error:
Error_Data.Abs_Error = abs(u_GT(OF_index) - u_est(OF_index));

% Median Error:
Error_Data.amp = median(Error_Data.Abs_Error(:));

% Mean Error:
Error_Data.amp_mean = mean(Error_Data.Abs_Error(:));

% Percentage of error greater than 0.1:
Error_Data.AE_10 = 100*sum(logical(Error_Data.Abs_Error(:) >=0.1))/(N*M);

% Mean of the top 10% highest errors:
Error_Data.rank_AE = sort(Error_Data.Abs_Error(:), 'descend');
if floor(N*M*0.1) > length(OF_index),
    Error_Data.MAE_10P = mean(Error_Data.rank_AE(1:floor(length(OF_index)*0.1)));
else
    Error_Data.MAE_10P = mean(Error_Data.rank_AE(1:floor(N*M*0.1)));
end
    
% ***************************************


% ******* Angular Error Analysis *******
% Calculate the error:
Error_Data.Angle_Error = real(acos((1 + real(conj(u_est(OF_index)).*u_GT(OF_index)))./(sqrt(1 + abs(u_GT(OF_index)).^2).*sqrt(1 + abs(u_est(OF_index)).^2))));

% Median Error:
Error_Data.angle = median(Error_Data.Angle_Error(:));

% Mean Error:
Error_Data.angle_mean = mean(Error_Data.Angle_Error(:));

% Percentage of error greater than 0.1:
Error_Data.rank_AnE = sort(Error_Data.Angle_Error(:), 'descend');

% Mean of the top 10% highest errors:
if floor(N*M*0.1) > length(OF_index),
    Error_Data.MAnE_10P = mean(Error_Data.rank_AnE(1:floor(length(OF_index)*0.1)));
else
    Error_Data.MAnE_10P = mean(Error_Data.rank_AnE(1:floor(N*M*0.1)));
end
% ***************************************

% ******* PSNR of the estimate *******
Error_Data.PSNR = CG_PSNR(u_GT(~isnan(u_est)), u_est(~isnan(u_est)));
% ************************************

end