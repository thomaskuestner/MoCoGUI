function metrics = fMetricROIs(SRegiResult,Labels,hImg,names,iNGates,cString,lDeform)
% calculate similarity metrics for ROIs in LandmarkGUI of original images
%
% input:
% SRegiResult       struct with registration info
% Labels            assigned label positions
% hImg              handle to insert output evaluation images
%
% output:
% metrics           similarity metrics
%
% -------------------------------------------------------------------------
% (c) 2016: Thomas Kuestner
% -------------------------------------------------------------------------

% delete demons folder from path as there is a deprecated version of
% moveimage from D.Kroon there.
warning('off','MATLAB:rmpath:DirNotFound')
rmpath(genpath(['.', filesep,'registration',filesep,'demons']));

results = Labels;
DefField = SRegiResult.SDeform;

% names={'G01','G02','G03','G04'};
[xline, yline] = find(~cellfun(@isempty,results.cROIs.G01));

metrics = cell(size(xline,1),iNGates-1);

lTmp = num2cell(repmat('''',length(cString),1));
lComma = num2cell(repmat(',',length(cString),1));
sString = cell2mat(strcat(lTmp,cString,lTmp,lComma).');
sString = sString(1:end-1);

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
            
            if(lDeform)
                RegROI = fDeformLine(MovROI, DefField, gateNo, sliceNo, SCT);
            else
                RegROI = MovROI;
            end
            % Get all rows and columns where the image is nonzero
            [nonZeroRows, nonZeroColumns] = find(RefROI | RegROI);
            % Get the cropping parameters
            topRow = min(nonZeroRows(:))-10; bottomRow = max(nonZeroRows(:))+10;
            leftColumn = min(nonZeroColumns(:))-10; rightColumn = max(nonZeroColumns(:))+10;
            % Extract a cropped image from the original.
            dImgRef = SRegiResult.dImg(:,:,sliceNo,1) .* RefROI;
            dImgRef = dImgRef(topRow:bottomRow, leftColumn:rightColumn);
            dImgReg = SRegiResult.dImg(:,:,sliceNo,gateNo) .* RegROI;
            dImgReg = dImgReg(topRow:bottomRow, leftColumn:rightColumn);
            
            eval(sprintf('metrics{lineNo,gateNo-1} = similarity_measure(dImgRef, dImgReg, %s);',sString));

            
        end
    end
end

for gateNo = 2:iNGates
    for lineNo = 1:size(xline)
        if ~isempty(hImg)
%             if(lineNo == 1) 
%                 lCreate = true;
%                 hlistbox = [];
%             else
%                 lCreate = false;
%             end
            hImg(metrics{lineNo,gateNo-1}, [gateNo,lineNo], cString);
        end
    end
end

