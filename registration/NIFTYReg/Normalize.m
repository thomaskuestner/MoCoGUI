function [ NormalizedMatrix ] = Normalize(Matrix, minValue, maxValue)
% Normalize maps the matrix elements to the user defined range of
% [minValue,maxValue], thus normalizing the matrix.

% Precision !!!
Matrix = double(Matrix);

Matrixmin = min(Matrix(:));
Matrixmax = max(Matrix(:));
% For the CSM
%Matrixmin = 0;
%Matrixmax = abs(sum(Matrix(:)));

m = double((maxValue-minValue))/double((Matrixmax-Matrixmin)); % Slope.

b1 = minValue-(m*Matrixmin);
b2 = maxValue-(m*Matrixmax);
b = (b1+b2)/2; % Y-intercept.

NormalizedMatrix = (Matrix*m)+b;

end

