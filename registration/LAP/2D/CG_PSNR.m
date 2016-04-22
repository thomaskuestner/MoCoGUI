function PSNR = CG_PSNR(Base_image, Recon_image)

[m,n] = size(Recon_image);

Max_I = max(max(abs(Base_image)));

MSE = sum(sum( (abs(Base_image-Recon_image).^2) ) )/(n*m);

PSNR = 10*log10( (Max_I).^2/MSE );
return