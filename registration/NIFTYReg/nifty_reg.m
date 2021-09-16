function [registered_image, displacement_field] = nifty_reg(reference, moving, options, tmp_path, clean_temp, niftyreg)
% Calls nifty reg for non-rigid registration and calculates displacement
% field.
%
% Input: 
% reference, moving: images to be registered
% options: string containing any options for reg_3fd (default empty)
% tmp_path: path where temporary nifti files will be written (default
% current directory)
% clean_temp: flag, 1 to delete temp files afterwards, 0 to keep them

% niftyreg = './NIFTI_code/niftyreg/reg-apps/';

if nargin < 3
    options = [];
end

if nargin < 4
    tmp_path = './Temps/';
end

if nargin < 5
    clean_temp = true;
end

if nargin < 6
    niftyreg = './bin/';
end

% Normalize inputs
reference = Normalize(abs(reference),0,1);
moving = Normalize(abs(moving),0,1);

% Save inputs as nifti
nii1 = make_nii(reference);
save_nii(nii1, [tmp_path, 'im1']);

nii2 = make_nii(moving);
save_nii(nii2, [tmp_path, 'im2']);

% Do registration
flag = system([niftyreg 'reg_f3d' ' -ref ' tmp_path, 'im1 '...
    ' -flo ' tmp_path, 'im2 ' ...
    ' -res ' tmp_path, 'registered'...
    ' -cpp ' tmp_path, 'cpp',...
    options]  );

if flag
    error('nifty reg failed')
end

registered = load_nii([tmp_path, 'registered']);

registered_image = registered.img;

% Calculate displacement field

system([niftyreg 'reg_transform', ' -ref ' tmp_path, 'im1 '...
    ' -disp ' [tmp_path, 'cpp.nii '],...
    [tmp_path, 'disp.nii ']
    ]);

disp = load_nii([tmp_path, 'disp.nii']);
displacement_field = disp.img;

% Cleanup temp files

if clean_temp
    
    system(['rm ' tmp_path, 'im1.hdr']);
    system(['rm ' tmp_path, 'im1.img']);
    system(['rm ' tmp_path, 'im1.mat']);   
    
    system(['rm ' tmp_path, 'im2.hdr']);
    system(['rm ' tmp_path, 'im2.img']);
    system(['rm ' tmp_path, 'im2.mat']);
    
    system(['rm ' tmp_path, 'registered.hdr']);
    system(['rm ' tmp_path, 'registered.img']);

    system(['rm ' tmp_path, 'cpp.nii']);

    system(['rm ' tmp_path, 'disp.nii'])

end

