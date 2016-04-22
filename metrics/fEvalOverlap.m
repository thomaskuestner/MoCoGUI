function [jaccard, dice, nOverlap, nMeanPixel] = fEvalOverlap(fixedROI, movingROI)
% evaluate overlap measures (LandmarkGUI,EvalGUI)
%
% input:
% fixedROI          ROI of reference image
% movingROI         ROIs of moving images
% gateNo            current moving state
% sliceNo           current slice (z) position
% SCT               sag/cor/tra image orientation
%
% output:
% jaccard           Jaccard coefficient
% dice              DICE coefficient
% nOverlap          amount of overlapping pixels
% nMeanPixel        mean amount of pixels in image ROIs
%
% -------------------------------------------------------------------------
% (c) 2015: Thomas Kuestner, Markus Bednarzyk, Verena Neumann
% -------------------------------------------------------------------------

% Jaccard coefficient
overlap_jac = fixedROI & movingROI;
num_overlap_jac = numel(nonzeros(overlap_jac));
jaccard = num_overlap_jac / (numel(nonzeros(fixedROI | movingROI)));

% Dice coefficient
overlap_dice = fixedROI & movingROI;
num_overlap_dice = numel(nonzeros(overlap_dice));
dice = (2*num_overlap_dice)/(numel(nonzeros(fixedROI))+numel(nonzeros(movingROI)));

% number of overlapping pixel
nOverlap = numel(nonzeros(fixedROI & movingROI));
nMeanPixel = (numel(nonzeros(fixedROI))+numel(nonzeros(movingROI)))/2;
end
