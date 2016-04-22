function [hfInsertPoint, hfInsertLineROI, hfPlotLabel, hfPlotMetrics] = fLandmarkCompare(iNPlots, dVoxelsize)
% comparison figure in LandmarkGUI
%
% input:
% iNPlots           amount of to be plotted figures (Points=#Gates, Lines, ROIs)
% dVoxelsize        voxelsize
%
% output:
% hfInsertPoint     handle to insert points
% hfInsertLineROI   handle to insert lines and ROIs
% hfPlotLabel       handle to plot legend image
%
% -------------------------------------------------------------------------
% (c) 2015: Thomas Kuestner
% -------------------------------------------------------------------------


hfig = figure;
fullscreen = get(0,'ScreenSize');
set(hfig,'position',[fullscreen(1)+10 fullscreen(2)+70 fullscreen(3)-20 fullscreen(4)-150]);
set(hfig, 'Color', 'w', 'Name', 'Feature-based landmark evaluation', 'NumberTitle','off');

try
    set(hfig, 'WindowScrollWheelFcn' , @fLandmarkCompare_WindowScrollWheelFcn);
catch
    warning('No scroll wheel functionality!');
end

panel1 = uipanel('Parent',hfig);
panel2 = uipanel('Parent',panel1);

set(panel1,'Position',[0 0 0.98 1],'BackgroundColor','k');
set(panel2,'Position',[0 -1 1 2],'BackgroundColor','k');
% set(gca,'Parent',panel2);

hSlider = uicontrol('Style','Slider','Parent',1,...
      'Units','normalized','Position',[0.98 0 0.02 1],...
      'Value',1,'Callback',{@slider_callback1,panel2});
  
% set some constants
iMaxCol = 4;
% iGapVals = [1,1];
iInitialMag = 400;
iNGatesPlot = iNPlots(1);

iNCols = min([iNGatesPlot+1, iMaxCol]); % + one for the legend
if(iNGatesPlot >= iMaxCol)
    iFac = ceil(iNGatesPlot/(iNCols-1));
else
    iFac = 1;
end
iNLines = iFac*(1+iNPlots(2)+iNPlots(3)+iNPlots(4));% ceil(max([iNCols, iNGatesPlot+1])/iMaxCol);

if(iNPlots(4) > 0)
    hlistbox = cell(1,iNGatesPlot);
    iStartIdx = (iNLines-1)*iNCols + 1;
    for iI=1:iNGatesPlot
        dPos = fCalcPos(iNLines, iNCols, iStartIdx);
        hlistbox{iI} = uicontrol('Style','listbox','Parent',panel2, 'Units','normalized',...
                            'Position',dPos, 'String','', 'BackgroundColor', 'k', 'ForegroundColor', 'w');
        uistack(hlistbox{iI},'top');
        iStartIdx = iStartIdx + 1;          
    end
end

iCurrPlotIdx = 1;
hfInsertPoint = @fInsertPoint; 
hfInsertLineROI = @fInsertLineROI;
hfPlotLabel = @fPlotLabel;
hfPlotMetrics = @fPlotMetrics;

    % slider callback
	function slider_callback1(src,eventdata,arg1)
        val = get(src,'Value');

        set(arg1,'Position',[0 -val 1 2])
    end

    function fLandmarkCompare_WindowScrollWheelFcn(hObject, eventdata, handles)
        dSliderstep = get(hSlider, 'SliderStep');
        if eventdata.VerticalScrollCount > 0
            newVal = max([0, get(hSlider,'Value') - dSliderstep(2)]);
        elseif(eventdata.VerticalScrollCount < 0)
            newVal = min([get(hSlider,'Max'), get(hSlider,'Value') + dSliderstep(2)]); 
        end
        if(exist('newVal','var'))
            set(hSlider, 'Value', newVal);
            set(panel2,'Position',[0 -newVal 1 2]);
        end
    end

    function fInsertPoint(dData, iGate)
        hAx = subtightplot(iNLines, iNCols, iCurrPlotIdx); %, iGapVals);    
        set(hAx,'Parent',panel2, 'Color', [0 0 0]);
        text(1,1,sprintf('Points (all) - image %02d\n %.4f +/- %.4f', iGate, dData(1), dData(2)), 'FontSize', 16, 'Color', [1 1 1]);
        axis(hAx,[0.9 1.5 0 1.1]);
%         axis(hAx,'tight');
        set(hAx,'XTick',[]); set(hAx,'XTickLabel',[]);
        set(hAx,'YTick',[]); set(hAx,'YTickLabel',[]);
        
        iCurrPlotIdx = iCurrPlotIdx + 1;
    end

    function fInsertLineROI(dData, iIdx, SCT, dMetric)
        iGate = iIdx(1);
        if(iGate == 2) % start a new line
            iCurrPlotIdx = ceil(iCurrPlotIdx/iNCols)*iNCols+1;
        end
        hAx = subtightplot(iNLines, iNCols, iCurrPlotIdx); %, iGapVals);
        set(hAx,'Parent',panel2);
        imshow(dData,'InitialMagnification', iInitialMag);
        
        if(strcmp(SCT,'sag')) % sag
            daspect([dVoxelsize(2) dVoxelsize(3) 1]);
        elseif strcmp(SCT,'tra') % tra
            daspect([dVoxelsize(3) dVoxelsize(1) 1]);
        else % cor
            daspect([dVoxelsize(1) dVoxelsize(2) 1]);
        end        
        
        if(length(dMetric) == 2) % lines
            text(size(dData,2)-1,3,sprintf('L%02d-img%02d',iIdx(2),iGate),'Color','w','HorizontalAlignment','right');
            text(1,3,['Dice: ',num2str(dMetric(1))],'color','w')
            text(1,8,['Dice_{mod}: ',num2str(dMetric(2))],'color','w')
            text(1+0.28*size(dData,2),5,['Dice_t: ',num2str(sum(dMetric))],'color','w')
        else % ROIs
            text(size(dData,2)-1,3,sprintf('R%02d-img%02d',iIdx(2),iGate),'Color','w','HorizontalAlignment','right');
            text(1,3,['Dice: ',num2str(dMetric(1))],'color','w')
        end
        iCurrPlotIdx = iCurrPlotIdx + 1;
    end

    function fPlotLabel(type)
        if(type == 0) % from original
            rgbimg = zeros(72,150,3);
            rgbimg(11:13,11:30,1) = 1;
            rgbimg(17:19,11:30,3) = 1;
            rgbimg(23:25,11:30,2) = 1;
            rgbimg(35:37,11:30,[1,3]) = 1;
            rgbimg(41:43,11:30,2:3) = 1;
            rgbimg(47:49,11:30,1:2) = 1;
            rgbimg(59:61,11:30,1:3) = 1;

            hAx = subtightplot(iNLines, iNCols, iCurrPlotIdx); %, iGapVals);
            set(hAx,'Parent',panel2);
            imshow(rgbimg,'InitialMagnification', iInitialMag)
            text(35,11,'reference mask','color','w')
            text(35,17,'transformed mask','color','w')
            text(35,23,'moving mask','color','w')
            text(35,35,'reference and transformed mask','color','w')
            text(35,41,'moving and transformed mask','color','w')
            text(35,47,'moving and reference mask','color','w')
            text(35,59,'all masks','color','w')
        elseif(type == 1) % from transformed
            rgbimg = zeros(72,150,3);
            rgbimg(23:25,11:30,1) = 1;
            rgbimg(29:31,11:30,3) = 1;

            rgbimg(41:43,11:30,[1,3]) = 1;

            hAx = subtightplot(iNLines, iNCols, iCurrPlotIdx); %, iGapVals);
            set(hAx,'Parent',panel2);
            imshow(rgbimg,'InitialMagnification', iInitialMag)
            text(35,23,'reference mask','color','w')
            text(35,29,'transformed mask','color','w')

            text(35,41,'reference and transformed mask','color','w')
        end
    end

    function fPlotMetrics(dData, iIdx, cString)
        % create print string
        lTmp = mat2cell(repmat(': ',length(cString),1),ones(length(cString),1),2);
        cString = strcat(cString,lTmp,cellfun(@(x) sprintf('%.4f',x), dData, 'UniformOutput', false).');

        cPrint = cell(length(cString)+2,1);
        cPrint{1} = sprintf('%%%% R%02d - img%02d %%%%', iIdx(2), iIdx(1));
        cPrint(2:end-1) = cString;
        cPrint{end} = '';
        cOldPrint = get(hlistbox{iIdx(1)-1},'String');
        cNewPrint = [cOldPrint;cPrint];
        set(hlistbox{iIdx(1)-1},'String',cNewPrint);
    end

    function pos_vec = fCalcPos(m,n,p,gap,marg_h,marg_w)        
        if (nargin<4) || isempty(gap),    gap=0.01;  end
        if (nargin<5) || isempty(marg_h),  marg_h=0.05;  end
        if (nargin<6) || isempty(marg_w),  marg_w=marg_h;  end
        if isscalar(gap),   gap(2)=gap;  end
        if isscalar(marg_h),  marg_h(2)=marg_h;  end
        if isscalar(marg_w),  marg_w(2)=marg_w;  end
        gap_vert   = gap(1);
        gap_horz   = gap(2);
        marg_lower = marg_h(1);
        marg_upper = marg_h(2);
        marg_left  = marg_w(1);
        marg_right = marg_w(2);

        %note n and m are switched as Matlab indexing is column-wise, while subplot indexing is row-wise :(
        [subplot_col,subplot_row]=ind2sub([n,m],p);  

        % note subplot suppors vector p inputs- so a merged subplot of higher dimentions will be created
        subplot_cols=1+max(subplot_col)-min(subplot_col); % number of column elements in merged subplot 
        subplot_rows=1+max(subplot_row)-min(subplot_row); % number of row elements in merged subplot   

        % single subplot dimensions:
        %height=(1-(m+1)*gap_vert)/m;
        %axh = (1-sum(marg_h)-(Nh-1)*gap(1))/Nh; 
        height=(1-(marg_lower+marg_upper)-(m-1)*gap_vert)/m;
        %width =(1-(n+1)*gap_horz)/n;
        %axw = (1-sum(marg_w)-(Nw-1)*gap(2))/Nw;
        width =(1-(marg_left+marg_right)-(n-1)*gap_horz)/n;

        % merged subplot dimensions:
        merged_height=subplot_rows*( height+gap_vert )- gap_vert;
        merged_width= subplot_cols*( width +gap_horz )- gap_horz;

        % merged subplot position:
        merged_bottom=(m-max(subplot_row))*(height+gap_vert) +marg_lower;
        merged_left=(min(subplot_col)-1)*(width+gap_horz) +marg_left;
        pos_vec=[merged_left merged_bottom merged_width merged_height];
    end
end

