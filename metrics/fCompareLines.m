function [dice,dice0,dice1] = fCompareLines(SRegiResult,Labels,hImg, names, iNGates)
% compare line segments in LandmarkGUI of original images
%
% input:
% SRegiResult       struct with registration info
% Labels            assigned label positions
% hImg              handle to insert output evaluation images
%
% output:
% dice              combined DICE coefficient of dice0+dice1
% dice0             standard DICE cofficient
% dice1             DICE coefficient of by 1 pixel shifted line segments
%
% -------------------------------------------------------------------------
% (c) 2015: Thomas Kuestner, Verena Neumann
% -------------------------------------------------------------------------

% delete demons folder from path as there is a deprecated version of
% moveimage from D.Kroon there.
warning('off','MATLAB:rmpath:DirNotFound')
rmpath(genpath(['.', filesep,'registration',filesep,'demons']));

results = Labels;
DefField = SRegiResult.SDeform;

% names={'G01','G02','G03','G04'};
[xline, yline] = find(~cellfun(@isempty,results.cLines.G01));

for lineNo = 1:size(xline,1)
    if isempty(results.cLines.G01{xline(lineNo),yline(lineNo)}) || isempty(results.cLines.G02{xline(lineNo),yline(lineNo)})...
            ||isempty(results.cLines.G03{xline(lineNo),yline(lineNo)}) || isempty(results.cLines.G04{xline(lineNo),yline(lineNo)})
        disp('Error: One of the Lines is missing')
    else
        if range(results.cLines.G01{xline(lineNo),yline(lineNo)}(:,1)) == 0
            SCT = 'sag';
            rw=3; cl=2;
            I=zeros(size(DefField(2).dFy,1),size(DefField(2).dFy,3));
        elseif range(results.cLines.G01{xline(lineNo),yline(lineNo)}(:,3)) == 0
            SCT = 'cor';
            rw=1; cl=2;
            I=zeros(size(DefField(2).dFy,1),size(DefField(2).dFy,2));
        elseif range(results.cLines.G01{xline(lineNo),yline(lineNo)}(:,2)) == 0
            SCT = 'tra';
            rw=2; cl=2;
            I=zeros(size(DefField(2).dFy,3),size(DefField(2).dFy,2));
        else
            error('Data not valid for this application')
        end
        
        rRef = results.cLines.(names{1}){xline(lineNo),yline(lineNo)}(:,rw);
        cRef = results.cLines.(names{1}){xline(lineNo),yline(lineNo)}(:,cl);
        RefLine = fCoords2Line(I,rRef,cRef);
               
        for gateNo = 2:iNGates
            
            if strcmp(SCT,'sag')
                sliceNo = results.cLines.(names{gateNo}){xline(lineNo),yline(lineNo)}(1,1);
                MRImage = squeeze(SRegiResult.dImgReg(:,sliceNo,:,1));
            elseif strcmp(SCT,'cor')
                sliceNo = results.cLines.(names{gateNo}){xline(lineNo),yline(lineNo)}(1,3);
                MRImage = SRegiResult.dImgReg(:,:,sliceNo,1);
            else % tra
                sliceNo = results.cLines.(names{gateNo}){xline(lineNo),yline(lineNo)}(1,2);
                MRImage = squeeze(SRegiResult.dImgReg(sliceNo,:,:,1));
            end
            
            rMov = results.cLines.(names{gateNo}){xline(lineNo),yline(lineNo)}(:,rw);
            cMov = results.cLines.(names{gateNo}){xline(lineNo),yline(lineNo)}(:,cl);
            MovLine = fCoords2Line(I,rMov,cMov);
            
            RegLine = fDeformLine(MovLine, DefField, gateNo, sliceNo, SCT);
            LineImage=zeros([size(MovLine),3]);
            MRImage = scaleImg(MRImage);
            % delete gray values where colored lines are plotted
            MRImage(RefLine == 1) = 0;
            MRImage(MovLine == 1) = 0;
            MRImage(RegLine == 1) = 0;
            % transform to rgb image
            MRImage = MRImage(:,:,[1 1 1]);
            % put MR image and Lines together
            LineImage(:,:,1) = RefLine;
            LineImage(:,:,2) = MovLine;
            LineImage(:,:,3) = RegLine;
            LineImage = 6*LineImage + 3*MRImage;
            
            % Get all rows and columns where the image is nonzero
            [nonZeroRows nonZeroColumns] = find(RefLine);
            % Get the cropping parameters
            topRow = min(nonZeroRows(:))-20; bottomRow = max(nonZeroRows(:))+10;
            leftColumn = min(nonZeroColumns(:))-10; rightColumn = max(nonZeroColumns(:))+10;
            
            % Extract a cropped image from the original.
            LineImage_c = LineImage(topRow:(min(bottomRow,size(LineImage,1))), leftColumn:(min(rightColumn,size(LineImage,2))),:);
            
            
            %% Evaluation
            [~,dice0(lineNo,gateNo-1),nOverl0(lineNo,gateNo-1),nMeanPx(lineNo,gateNo-1)] = ...
                fEvalOverlap(RefLine,RegLine);
            
            % clean registered line from pixel that have overlapped
            RegLine_clean = RegLine - (RegLine&RefLine);
            
            % shift registered line up and down
            RegLine_up = [RegLine_clean(2:end,:);zeros(1,size(MovLine,2))];
            RegLine_down = [zeros(1,size(MovLine,2));RegLine_clean(1:end-1,:)];
            
            % count the pixel that overlap then
            [~,~,nOverlUp(lineNo,gateNo-1)] = fEvalOverlap(RefLine,RegLine_up);
            [~,~,nOverlDown(lineNo,gateNo-1)] = fEvalOverlap(RefLine,RegLine_down);
            
            % delete pixel in shifted registered line that overlapped with the reference line
            RegLine_up_clean = RegLine_up - (RegLine_up & RefLine);
            RegLine_down_clean = RegLine_down - (RegLine_down & RefLine);
            % move line back to origin
            RegLine_clean_ud = [zeros(1,size(MovLine,2));RegLine_up_clean(1:end-1,:)];
            RegLine_clean_du = [RegLine_down_clean(2:end,:);zeros(1,size(MovLine,2))];
            % combine both lines
            RegLine_clean2 = RegLine_clean_ud .* RegLine_clean_du;
            
            % shift pixel left and right
            RegLine_left = [RegLine_clean2(:,2:end),zeros(size(MovLine,1),1)];
            RegLine_right = [zeros(size(MovLine,1),1),RegLine_clean2(:,1:end-1)];
            % count the pixel that overlap then
            [~,~,nOverlLeft(lineNo,gateNo-1)] = fEvalOverlap(RefLine,RegLine_left);
            [~,~,nOverlRight(lineNo,gateNo-1)] = fEvalOverlap(RefLine,RegLine_right);
            
            nOverl1 = nOverlLeft+nOverlDown+nOverlRight+nOverlRight+nOverlUp;
            dice1 = nOverl1 ./ nMeanPx;
            
            if ~isempty(hImg)
                fHandle = hImg{1};
                fHandle(LineImage_c, [gateNo,lineNo], SCT, [dice0(lineNo,gateNo-1),dice1(lineNo,gateNo-1)]);
%                 if gateNo == 2
%                     hplot = figure;
%                     subtightplot(2,2,1)
%                     imshow(LineImage_c,'InitialMagnification', 400)
%                     text(size(LineImage_c,2)-2,5,'Gate 02',...
%                        'color','w','HorizontalAlignment','right')
%                 else
% %                     figure
%                     subtightplot(2,2,gateNo-1)
%                     imshow(LineImage_c,'InitialMagnification', 400)
%                     text(size(LineImage_c,2)-2,5,['Gate 0',num2str(gateNo)],...
%                        'color','w','HorizontalAlignment','right')
%                 end
%                 if strcmp(SCT,'sag')
%                     daspect([1.95 5 1.95])
%                 end
%                 text(5,5,['Dice coef 0: ',num2str(dice0(lineNo,gateNo-1))],'color','w')
%                 text(5,10,['Dice coef 1: ',num2str(dice1(lineNo,gateNo-1))],'color','w')
            end
            
        end
        
        nOverl1 = nOverlLeft+nOverlDown+nOverlRight+nOverlRight+nOverlUp;
        dice1 = nOverl1 ./ nMeanPx;
        
        dice = dice1+dice0;
        
    end
    if ~isempty(hImg)
        fHandle = hImg{2};
        fHandle(0);
        
%         rgbimg = zeros(72,150,3);
%         rgbimg(11:13,11:30,1) = 1;
%         rgbimg(17:19,11:30,3) = 1;
%         rgbimg(23:25,11:30,2) = 1;
%         rgbimg(35:37,11:30,[1,3]) = 1;
%         rgbimg(41:43,11:30,2:3) = 1;
%         rgbimg(47:49,11:30,1:2) = 1;
%         rgbimg(59:61,11:30,1:3) = 1;
%         
%         subtightplot(2,2,4)
%         imshow(rgbimg,'InitialMagnification', 400)
%         text(35,11,'reference mask','color','w')
%         text(35,17,'registered mask','color','w')
%         text(35,23,'moving mask','color','w')
%         text(35,35,'reference and registered mask','color','w')
%         text(35,41,'moving and registered mask','color','w')
%         text(35,47,'moving and reference mask','color','w')
%         text(35,59,'all masks','color','w')
%         set(hplot,'units','normalized','outerposition',[0 0.5 0.9 0.5])
    end
end
end
