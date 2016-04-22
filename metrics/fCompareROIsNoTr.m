function dice = fCompareROIsNoTr(SRegiResult,Labels,hImg,names,iNGates)
% compare ROIs in LandmarkGUI of original images
%
% input:
% SRegiResult       struct with registration info
% Labels            assigned label positions
% hImg              handle to insert output evaluation images
%
% output:
% dice              DICE coefficient for measuring regional overlap
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
[xline, yline] = find(~cellfun(@isempty,results.cROIs.G01));

for lineNo = 1:size(xline,1)
    if isempty(results.cROIs.G01{xline(lineNo),yline(lineNo)}) || isempty(results.cROIs.G02{xline(lineNo),yline(lineNo)})...
            ||isempty(results.cROIs.G03{xline(lineNo),yline(lineNo)}) || isempty(results.cROIs.G04{xline(lineNo),yline(lineNo)})
        disp('Error: One of the Lines is missing')
    else
        if range(results.cROIs.G01{xline(lineNo),yline(lineNo)}(:,1)) == 0
            SCT = 'sag';
            rw=3; cl=2;
            I=zeros(size(DefField(2).dFy,1),size(DefField(2).dFy,3));
        elseif range(results.cROIs.G01{xline(lineNo),yline(lineNo)}(:,3)) == 0
            SCT = 'cor';
            rw=1; cl=2;
            I=zeros(size(DefField(2).dFy,1),size(DefField(2).dFy,2));
        elseif range(results.cROIs.G01{xline(lineNo),yline(lineNo)}(:,2)) == 0
            SCT = 'tra';
            rw=2; cl=2;
            I=zeros(size(DefField(2).dFy,3),size(DefField(2).dFy,2));
        else
            error('Data not valid for this application')
        end
        
        rRef = results.cROIs.(names{1}){xline(lineNo),yline(lineNo)}(:,rw);
        cRef = results.cROIs.(names{1}){xline(lineNo),yline(lineNo)}(:,cl);
        
        RefROI = roipoly(I,rRef,cRef);
        
        for gateNo = 2:iNGates
            
            if strcmp(SCT,'sag')
                sliceNo = results.cROIs.(names{gateNo}){xline(lineNo),yline(lineNo)}(1,1);
            elseif strcmp(SCT,'cor')
                sliceNo = results.cROIs.(names{gateNo}){xline(lineNo),yline(lineNo)}(1,3);
            else % tra
                sliceNo = results.cROIs.(names{gateNo}){xline(lineNo),yline(lineNo)}(1,2);
            end
            
            rMov = results.cROIs.(names{gateNo}){xline(lineNo),yline(lineNo)}(:,rw);
            cMov = results.cROIs.(names{gateNo}){xline(lineNo),yline(lineNo)}(:,cl);
            MovROI = roipoly(I,rMov,cMov);
            
            RegROI = MovROI;
            %%
            MRImage = SRegiResult.dImg(:,:,sliceNo,1);
            MRImage = scaleImg(MRImage);
            MRImage = MRImage(:,:,[1 1 1]);
            ROIImage(:,:,1) = RefROI;
            
            ROIImage(:,:,3) = RegROI;
            % Get all rows and columns where the image is nonzero
            [nonZeroRows nonZeroColumns] = find(RefROI);
            % Get the cropping parameters
            topRow = min(nonZeroRows(:))-10; bottomRow = max(nonZeroRows(:))+10;
            leftColumn = min(nonZeroColumns(:))-10; rightColumn = max(nonZeroColumns(:))+10;
            % Extract a cropped image from the original.
            LineImage = 5*ROIImage.*MRImage + MRImage;
            LineImage_c = LineImage(topRow:bottomRow, leftColumn:rightColumn,:);
            iptsetpref('ImshowBorder','tight');
%             if images
%                 if gateNo == 2
%                     hplot = figure;
%                     subtightplot(1,4,1)
%                     imshow(LineImage_c,'InitialMagnification', 400)
%                     text(size(LineImage_c,2)-5,5,'Gate 02',...
%                        'color','w','HorizontalAlignment','right')
%                 else
%                     subtightplot(1,4,gateNo-1)
%                     imshow(LineImage_c,'InitialMagnification', 400)
%                     text(size(LineImage_c,2)-5,5,['Gate 0',num2str(gateNo)],...
%                        'color','w','HorizontalAlignment','right')
%                 end
%                 if strcmp(SCT,'sag')
%                     daspect([1.95 5 1.95])
%                 end
%             end
            %%
            [~,dice(lineNo,gateNo-1)] = fEvalOverlap(RefROI,RegROI);
%             if images
%                 text(5,5,['Dice coef: ',num2str(dice(lineNo,gateNo-1))],'color','w')
%             end
            if ~isempty(hImg)
                fHandle = hImg{1};
                fHandle(LineImage_c, [gateNo,lineNo], SCT, dice(lineNo,gateNo-1));
            end
        end
        %%
        if ~isempty(hImg)
            fHandle = hImg{2};
            fHandle(1);
        end
%         if images
%         rgbimg = zeros(72,150,3);
%         rgbimg(23:25,11:30,1) = 1;
%         rgbimg(29:31,11:30,3) = 1;
%         
%         rgbimg(41:43,11:30,[1,3]) = 1;
%         
%         subtightplot(1,4,4)
%         imshow(rgbimg,'InitialMagnification', 400)
%         text(35,23,'reference mask','color','w')
%         text(35,29,'registered mask','color','w')
%         
%         text(35,41,'reference and registered mask','color','w')
%         
%         set(hplot,'units','normalized','outerposition',[0 0 0.9 0.5])
%         end
  %%
    end
end


