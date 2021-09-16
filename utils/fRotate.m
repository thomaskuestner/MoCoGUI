function [dImg, dFlow] = fRotate(dImg, dFlow, sInputOri, sOutputOri)
% dImg          2D to 5D image
% dFlow         same number of dimensions as dImg with last dimension storing flow in x/y/z direction, i.e. if ndims(dImg)==3 -> size(dFlow) = [size(dImg,1),size(dImg,2),size(dImg,3), 2 or 3], last dim holds cat(4, ux, uy, uz)
% sInputOri     anatomical input orientation (e.g. tra) 
% sOutputOri    desired output orientation (e.g. sag)


if(isnumeric(sInputOri))  % directly supplied iDir 
    iDirs = sInputOri;
    sInputOri = 'cor'; % assumed
else
    iDirs = [];
end

if(length(sInputOri) == 6)
    sOri = sInputOri;
    if(strcmp(sOri(1:2),'LR') || strcmp(sOri(1:2),'RL'))
        if(strcmp(sOri(3:4),'HF') || strcmp(sOri(3:4),'FH'))
            sInputOri = 'cor';
            sExpectedOri = 'RLHFAP';
        else
            sInputOri = 'tra';
            sExpectedOri = 'RLAPFH';
        end
    else
        sInputOri = 'sag';
        sExpectedOri = 'APHFLR';
    end
else
    if(strcmp(sInputOri, 'cor'))
        sOri = 'RLHFAP';
    elseif(strcmp(sInputOri,'tra'))
        sOri = 'RLAPFH';
    elseif(strcmp(sInputOri,'sag'))
        sOri = 'APHFLR';
    end
%     iMult = ones(1,3);
end

if(isempty(iDirs))
    if(strcmp(sInputOri, 'cor'))
        switch sOutputOri
            case 'cor'
                iDirs = 0;
            case 'tra'
                iDirs = 4;
            case 'sag'
                iDirs = 2;
        end
    elseif(strcmp(sInputOri, 'tra'))
        switch sOutputOri
            case 'cor'
                iDirs = 3;
            case 'tra'
                iDirs = 0;
            case 'sag'
                iDirs = [2, 6];
        end
    elseif(strcmp(sInputOri, 'sag'))
        switch sOutputOri
            case 'cor'
                iDirs = 1;
            case 'tra'
                iDirs = [3, 5];
            case 'sag'
                iDirs = 0;
        end
    end
end
    

for iDir = iDirs
    
    if iDir == 0  % do nothing
       continue;
    end

    % - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Identify orientation change
    % Normal, left mouse button -> volume rotation operation
    if iDir == 1 % Moved mouse to right
    %     SData(i).iActiveImage = uint16(SMouse.dAxesStartPos(1, 1));
        iPermutation = [1 3 2]; iFlipdim = 2; iFlipdimFlow = 2; 
        u = [0 1 0]; % around y axis
        theta = -pi/2;
        sOri = sOri([6 5 3 4 1 2]);
    end
    if iDir == 2 % Moved mouse to left
%         SData(i).iActiveImage = uint16(size(SData(i).dImg, 2) - SMouse.dAxesStartPos(1, 1) + 1);
        iPermutation = [1 3 2]; iFlipdim = 3; iFlipdimFlow = 3; 
        u = [0 1 0]; % around y axis
        theta = pi/2;
        sOri = sOri([5 6 3 4 2 1]); 
    end
    if iDir == 3 % Moved mouse up
%         SData(i).iActiveImage = uint16(size(SData(i).dImg, 1) - SMouse.dAxesStartPos(1, 2) + 1);
        iPermutation = [3 2 1]; iFlipdim = 3; iFlipdimFlow = 3; 
        u = [1 0 0]; % around x axis
        theta = -pi/2;
        sOri = sOri([1 2 5 6 4 3]);
    end
    if iDir == 4 % Moved mouse down
%         SData(i).iActiveImage = uint16(SMouse.dAxesStartPos(1, 2));
        iPermutation = [3 2 1]; iFlipdim = 1; iFlipdimFlow = 1; 
        u = [1 0 0]; % around x axis
        theta = pi/2;
        sOri = sOri([1 2 6 5 3 4]);
    end

    % - - - - - - - - - - - - - - - - - - - - - - - - -
    % Shift key or right mouse button -> rotate in-plane

    if iDir == 5 % 90° clock-wise
        iPermutation = [2 1 3]; iFlipdim = 2; iFlipdimFlow = 2; lInPlane = true; 
        u = [0 0 1]; % around z axis
        theta = -pi/2;
%                                 SData(i).iRot(1:2) = [SData(i).iRot(2), mod(SData(i).iRot(2)-90,360)]; % [from, to]
        sOri = sOri([4 3 1 2 5 6]); 
    end
    if iDir == 6 % 90° counter clock-wise
        iPermutation = [2 1 3]; iFlipdim = 1; iFlipdimFlow = 1; lInPlane = true; 
        u = [0 0 1]; % around z axis
        theta = pi/2;
%                                 SData(i).iRot(1:2) = [SData(i).iRot(2), mod(SData(i).iRot(2)+90,360)];
        sOri = sOri([3 4 2 1 5 6]); 
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - -


    % - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Apply the transformation
    switch ndims(dImg) 
        case 2
            dImg =  flipdim(permute(dImg,  iPermutation), iFlipdim);
            iPermutationFlow = [iPermutation(1:2), 3];
        case 3
            dImg =  flipdim(permute(dImg,  iPermutation), iFlipdim);
            iPermutationFlow = [iPermutation, 4];
        case 4
            dImg = flipdim(permute(dImg, [iPermutation, 4, 5]), iFlipdim);
            iPermutationFlow = [iPermutation, 4, 5];           
        case 5
            dImg = flipdim(permute(dImg, [iPermutation, 4, 5]), iFlipdim);
            iPermutationFlow = [iPermutation, 4, 5, 6];
    end
    
%     lMask = flipdim(permute(lMask, iPermutation), iFlipdim);

    % check which dimension was swapped
%     iMultPrev = iMult;
%     iMult = ones(1,3);
%     if(strcmp(sOri(1:2), 'LR') || strcmp(sOri(1:2), 'FH') || strcmp(sOri(1:2), 'PA'))
%         iMult(1) = -1;
%     end
%     if(strcmp(sOri(3:4), 'LR') || strcmp(sOri(3:4), 'FH') || strcmp(sOri(3:4), 'PA'))
%         iMult(2) = -1;
%     end
%     if(strcmp(sOri(5:6), 'LR') || strcmp(sOri(5:6), 'FH') || strcmp(sOri(5:6), 'PA'))
%         iMult(3) = -1;
%     end    
%     iMultTurn = iMult;

%     fprintf('%s: prev = [%d, %d, %d] * curr = [%d, %d, %d] --> [%d, %d, %d]\n', sOri, iMultPrev(1), iMultPrev(2), iMultPrev(3), iMult(1), iMult(2), iMult(3), iMultTurn(1), iMultTurn(2), iMultTurn(3));
    
    if(~isempty(dFlow))
        switch ndims(dFlow) 
            case 3
                iPermutationFlow = [iPermutation(1:2), 3];
            case 4
                iPermutationFlow = [iPermutation, 4];
            case 5
                iPermutationFlow = [iPermutation, 4, 5];           
            case 6
                iPermutationFlow = [iPermutation, 4, 5, 6];
        end
        
        dFlow = flipdim(permute(dFlow, iPermutationFlow), iFlipdimFlow);
        
        switch find(u)
            case 1 % X axis rotation matrix ; u = i = [1 0 0]'
                Rm = [1          0           0;
                      0 cos(theta) -sin(theta);
                      0 sin(theta)  cos(theta)];
            case 2 % Y axis rotation matrix ; u = j = [0 1 0]'
                Rm = [cos(theta)   0  -sin(theta);
                      0            1           0;
                      sin(theta)  0  cos(theta)];
            case 3 % Z axis rotation matrix ; u = k = [0 0 1]'
                Rm = [cos(theta) -sin(theta) 0;
                      sin(theta)  cos(theta) 0;
                      0           0          1];
            otherwise % Any u axis rotation matrix
                
                u = u/norm(u);
                
                Rm = [u(1,1)^2+cos(theta)*(1-u(1,1)^2) (1-cos(theta))*u(1,1)*u(2,1)-u(3,1)*sin(theta) (1-cos(theta))*u(1,1)*u(3,1)+u(2,1)*sin(theta);
                      (1-cos(theta))*u(1,1)*u(2,1)+u(3,1)*sin(theta) u(2,1)^2+cos(theta)*(1-u(2,1)^2) (1-cos(theta))*u(2,1)*u(3,1)-u(1,1)*sin(theta);
                      (1-cos(theta))*u(1,1)*u(3,1)-u(2,1)*sin(theta) (1-cos(theta))*u(2,1)*u(3,1)+u(1,1)*sin(theta) u(3,1)^2+cos(theta)*(1-u(3,1)^2)];
        end
        
        dFlow = shiftdim(dFlow, ndims(dFlow)-1);
        dFlow = mtimesx(Rm, dFlow);
        dFlow = shiftdim(dFlow, 1);
        
        
%         iMinusDim = find(iMultTurn == -1);
%         iMinusDimPrev = find(iMultPrev == -1);
        
%         iMinusDimSwp = find(iSwapComp < 0);

%         if(ndims(dFlow) == 6)
% %             if(~isempty(iMinusDimPrev)), dFlow(:,:,:,:,:,iMinusDimPrev) = -dFlow(:,:,:,:,:,iMinusDimPrev); end;
%             dFlow = dFlow(:,:,:,:,:,abs(iSwapComp));
% %             if(~isempty(iMinusDim)), dFlow(:,:,:,:,:,iMinusDim) = -dFlow(:,:,:,:,:,iMinusDim); end;
%         elseif(ndims(dFlow) == 5)
% %             if(~isempty(iMinusDimPrev)), dFlow(:,:,:,:,iMinusDimPrev) = -dFlow(:,:,:,:,iMinusDimPrev); end;
%             dFlow = dFlow(:,:,:,:,abs(iSwapComp));
% %             if(~isempty(iMinusDim)), dFlow(:,:,:,:,iMinusDim) = -dFlow(:,:,:,:,iMinusDim); end;
%         elseif(ndims(dFlow) == 4)
% %             if(~isempty(iMinusDimPrev)), dFlow(:,:,:,iMinusDimPrev) = -dFlow(:,:,:,iMinusDimPrev); end;
%             dFlow = dFlow(:,:,:,abs(iSwapComp));
%             if(~isempty(iMinusDimSwp)), dFlow(:,:,:,iMinusDimSwp) = -dFlow(:,:,:,iMinusDimSwp); end;
% %             if(~isempty(iMinusDim)), dFlow(:,:,:,iMinusDim) = -dFlow(:,:,:,iMinusDim); end;
%         else % 3
% %             if(~isempty(iMinusDimPrev)), dFlow(:,:,iMinusDimPrev) = -dFlow(:,:,iMinusDimPrev); end;
%             dFlow = dFlow(:,:,abs(iSwapComp(1:2)));
% %             if(~isempty(iMinusDim)), dFlow(:,:,iMinusDim) = -dFlow(:,:,iMinusDim); end;
%         end
    end
%     SData(i).iMult = iMult;
end



