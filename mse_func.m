function [R] = mse_func(a,b)
% [D] = mse_func(a,b)
%
%Computes the means squared error of the image
%D is the MSE
%R is the RMSE
%P is the Relative RMSE
%a is the original vector
%b is the estimate of the original vector
%Phillip K Poon 2012

D = mean( (a(:)-b(:)).^2 );
R = sqrt(D);
temp = [a b];
P = 100 * (R/max(temp(:)));

end

