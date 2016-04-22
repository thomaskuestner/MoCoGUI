function [bwe_open] = fCoords2Line(I,r,c)
% convert coordinate points to line segment (LandmarkGUI)
%
% input:
% I                 corresponding image
% r,c               row and column pixel positions of ROI
%
% output:
% bwe_open          edge image after morphological opening
%
% -------------------------------------------------------------------------
% (c) 2015: Thomas Kuestner, Verena Neumann
% -------------------------------------------------------------------------

% get mask
bw = roipoly(I,r,c);
% get edge of the mask
bwe = edge(bw);   

% second polygon to eliminate 'closure' line
midpt = (c(round(length(r)/2))+c(1)); % some value not on the actual line
bwe_c = edge(roipoly(I,r([1,round(length(r)/2),end]),[c(1),midpt,c(end)]));

bwe_open = bwe.*(1-bwe_c);