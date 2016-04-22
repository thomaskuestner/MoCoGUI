function  giplwrite(filename, dat, type, vdims)
%	Write a gipl file to file.
%
%	vdims  voxel dimensions [y x z t] (unspecified dimensions default to 1) 
%        
%	type can be:
%	'complex_float'
%	'float'
%	'short'           (Data is scaled and rounded).
%	'short_noscale' % this is what we use for lreg registration
%	'u_char', 'uint8' (Data is scaled and rounded unless it is
%				of type uint8 upon input).
%	'uint16', 'u_short'  No scaling or rounding
%
%   Writes the file big endian (UNIX style), whether on a PC or Unix.
%

numDims = ndims(dat) ;

if ~isreal(dat) & (strcmp(type,'complex_float')==0)
  disp([ 'Input data is complex but you are not writing to complex float',...
         ' using modulus values.' ])
  dat = abs(dat) ;
end
 
for idim = 1:4
  siz(idim) = size(dat,idim) ;
end

[fid, message] = fopen(filename,'w','ieee-be') ;
if fid == -1
  a=['unable to open file ',filename,' for writing.'];
  error(a);
end

if nargin < 4
  vdims = [ 1 1 1 1] ;
end

vox_dims = [ 1 1 1 1] ;
vox_dims(1:length(vdims)) = vdims(:) ;

Y = 1 ; X = 2 ; Z = 3 ; T = 4 ;

dmin = min(min(min(min(dat)))) ;
dmax = max(max(max(max(dat)))) ;

switch type
  case 'complex_float'
    image_type = 192 ;
  case 'float'
    image_type = 64 ;
  case 'short'
    image_type = 15 ;
    warning([ 'Doubles will be written scaled to write as shorts.'])
    
    shrt = 16383 ;
    if dmin >= 0   % input data is positive
      dmin_dat = dmin ;
      dat  = dat - dmin_dat ;
      dmin = 0 ;
      dmax = dmax - dmin_dat ;
      
      dat = dat ./ dmax .* shrt ;
      dmin = dmin/ dmax  * shrt ;
      dmax = dmax / dmax * shrt ;
      
    else  % data contains negative numbers
          % keep 0 as mid gray
      extrem = max([dmax abs(dmin)]) ;
      dat = dat .* (shrt / extrem) ;
      dmin = dmin * (shrt / extrem) ;
      dmax = dmax * (shrt / extrem) ;
    end
    dat = round(dat) ;
    dmin = round(dmin) ;
    dmax = round(dmax) ;
  case 'short_noscale'
    image_type = 15 ;
  case {'u_char','uint8'}
    image_type = 8 ;
    if ~isa(dat,'uint8')
      dat = double(dat) ;
      dmin = min(min(min(min(dat)))) ;
      dmax = max(max(max(max(dat)))) ;
      rang = dmax - dmin ;
      ul = -0.499 ; uu= 255.499 ;
      dat = (dat - dmin)/rang *(uu-ul) + ul ;
      dat = round(dat) ;
      dmin = min(min(min(min(dat)))) ;
      dmax = max(max(max(max(dat)))) ;
    end
 case {'u_short','uint16'}
  image_type = 16 ;
  if dmin < 0 | dmax > 65535
    warning([' Data out of range of uint16 !!'])
  end
  
  otherwise
    error([ 'Unknown output image type ',type])
end

%write some of the giplheader

fwrite(fid,siz(X),'short') ;  % nx

fwrite(fid,siz(Y),'short') ;  % ny

fwrite(fid,siz(Z),'short') ;  % nz

fwrite(fid,siz(T),'short') ;  % nt

fwrite(fid,image_type,'short') ;

fwrite(fid,vox_dims(X),'float') ; % pixel dimensions
fwrite(fid,vox_dims(Y),'float') ; % pixel dimensions
fwrite(fid,vox_dims(Z),'float') ; % pixel dimensions
fwrite(fid,vox_dims(T),'float') ; % pixel dimensions

space = ' ' ;
for ibyte = 1:80 
 fwrite(fid,space,'char') ; % patDesc
end

track=eye(4)
fwrite(fid,track(1,1),'float') ;        % matrix elements
fwrite(fid,track(1,2),'float') ;        % matrix elements
fwrite(fid,track(1,3),'float') ;        % matrix elements
fwrite(fid,track(1,4),'float') ; % matrix elements
fwrite(fid,track(2,1),'float') ;        % matrix elements
fwrite(fid,track(2,2),'float') ;        % matrix elements
fwrite(fid,track(2,3),'float') ;        % matrix elements
fwrite(fid,track(2,4),'float') ; % matrix elements
fwrite(fid,track(3,1),'float') ;        % matrix elements
fwrite(fid,track(3,2),'float') ;        % matrix elements
fwrite(fid,track(3,3),'float') ;        % matrix elements
fwrite(fid,track(3,4),'float') ; % matrix elements

fwrite(fid,0,'float') ; % identifier

for ibyte = 1:28 
 fwrite(fid,space,'char') ; % spare  
end

fwrite(fid,00,'char') ; % orientation
fwrite(fid,2,'char') ; % flag1



fwrite(fid,dmin,'double') ; % data minimum, abs if complex
fwrite(fid,dmax,'double') ; % data maximum, abs if complex

fwrite(fid,0,'double') ; % origins
fwrite(fid,0,'double') ; % origins
fwrite(fid,0,'double') ; % origins
fwrite(fid,0,'double') ; % origins

fwrite(fid,0,'float') ; % pixval offset
fwrite(fid,1,'float') ; % pixval cal   
fwrite(fid,0,'float') ; % user def 1   
fwrite(fid,0,'float') ; % user def 2   
fwrite(fid,719555000,'long') ; % magic number DONT FIDDLE WITH THIS!!!!

switch type
  case 'complex_float'
    %pre-allocate matrices to speed up code
    ridat = zeros([siz(Y) 2*siz(X) siz(Z) siz(T)]) ;

    % compile ridat from dat, then write ridat.
    % Matlab matrices are row-column order and writen/read
    % in column order.
	
    ridat(:, 1:2:2*siz(X) , :, :) = real(dat) ;
    ridat(:, 2:2:2*siz(X) , :, :) = imag(dat) ;
    
    ridat = permute(ridat,[2 1 3 4]) ;
    count = fwrite(fid,ridat,'float') ;
  case 'float'
    wdat = permute(dat,[2 1 3 4]) ;
    count = fwrite(fid,wdat,'float') ;
  case 'short'
    wdat = permute(dat,[2 1 3 4]) ;
    count = fwrite(fid,wdat,'short') ;
  case 'short_noscale'
    wdat = permute(dat,[2 1 3 4]) ;
    count = fwrite(fid,wdat,'short') ;
  case {'u_char','uint8'}
    wdat = permute(dat,[2 1 3 4]) ;
    count = fwrite(fid,wdat,'uchar') ;
  case {'u_short','uint16'}
    wdat = permute(dat,[2 1 3 4]) ;
    count = fwrite(fid,wdat,'uint16') ;
end

fclose(fid) ;



