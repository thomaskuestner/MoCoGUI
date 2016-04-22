function CG_PlotGraphs(u_GT, u_est, Error_Data, type)

switch type.data
    case 1
        % Plot visual comparison between u_GT and u_est:
        
        x = type.x;
        y = type.y;
        max_amp = type.max_amp;
        amp = type.amp;
        
        f1 = figure; set(f1,'position',get(0,'screensize')); subplot(1,2,1);  plotflow(x,y, u_GT, max_amp);
        set(gca, 'fontsize', 12);
        title(['Ground Truth Optical Flow,    Max Amplitude = ', num2str(amp)]);
        xlabel('x component of flow');
        ylabel('y component of flow');

        subplot(1,2,2); plotflow(x,y, u_est, max_amp);
        set(gca, 'fontsize', 12);
        title(['Estimate of Optical Flow using ', type.title]);
        xlabel('x component of flow');
        ylabel('y component of flow');

        % Create textbox
        annotation(f1,'textbox',...
            [0.17 0.133116496065155 0.47510980966325 0.0521472392638037],...
            'String',[{['Median Absolute Error = ', num2str(Error_Data.amp,4),',    Median Angular Error = ', num2str(Error_Data.angle,4), ',    PSNR = ', num2str(Error_Data.PSNR,4),'dB']};...
            {['Percent of Pixels with an Absolute Error > 0.1 = ', num2str(Error_Data.AE_10,4), '%,    Mean of the highest 10% of Absolue Errors = ', num2str(Error_Data.MAE_10P,4)]}],...
            'FontSize',14,...
            'FontName','Arial',...
            'FitBoxToText','on',...
            'HorizontalAlignment','center',...
            'BackgroundColor',[1 1 1]);
    
    case 2
        
        c_max = type.c_max;
        map = zeros(size(u_GT));
        map(type.OF_index) = 1;
        
        f1 = figure; set(f1,'position',get(0,'screensize')); imagesc(abs(u_GT - u_est).*map);
        colorbar('ylim', [0,c_max]); 
        caxis([0, c_max]);
        
        set(gca, 'fontsize', 12);
        title(['Absolute Error in flow using ', type.title]);
        xlabel('x axis');
        ylabel('y axis');
    
    case 3
        color_data = [0, 0, 1; 1, 0, 0; 0, 1, 0; 1, 1, 0; 0, 1, 1; 1, 0, 1];
        max_amp = type.max_amp;
        [M,N] = size(u_GT);
        
        % generate figure for absolute error
        f1 = figure; set(f1,'position',get(0,'screensize'));
        axes1 = axes('Parent', f1,'YScale','log','YMinorTick','on',...
            'Position',[0.0607613469985359 0.11 0.42606149341142 0.815]);
        axis([0, N*M*1.05, 0, ceil(max_amp)]);
        box(axes1,'on');
        hold(axes1,'all');
        
        Data_holder = [];
        for l = 1:length(Error_Data),
            Data_holder = [Data_holder, sort(Error_Data(l).Abs_Error(:))];
        end

        % Y-axis in log scale
        semilogy1 = semilogy(Data_holder,'Parent',axes1);
        set(gca, 'fontsize', 12);
        for l = 1:length(Error_Data),
            color = color_data(l,:);
            tag = ['set(semilogy2(', int2str(l),'), ''Color'', color,''DisplayName'', type.name',int2str(l),');'];
            eval(tag);
        end
%         set(semilogy2(2),'Color',[1 0 0],'DisplayName',type.name2);

        % Create xlabel
        xlabel('Image Pixels');

        % Create ylabel
        ylabel('Log Absolute Error, log_{10}|u - u_{est}|');

        % Create title
        title('Absolue error between the estimate and actual Flow, in ascending order');

        % Create legend
        legend1 = legend(axes1,'show');
        set(legend1,...
            'Position',[0.334797462176672 0.134834167783184 0.142752562225476 0.0828220858895706]);

        % generate figure for angular error
        axes2 = axes('Parent',f1,'YScale','log','YMinorTick','on',...
            'Position',[0.550512445095168 0.11 0.420204978038067 0.815]);
        axis([0, N*M*1.05, 0, 100]);
        box(axes2,'on');
        hold(axes2,'all');
        
        Data_holder = [];
        for l = 1:length(Error_Data),
            Data_holder = [Data_holder, sort(Error_Data(l).Angle_Error(:))];
        end

        % Create multiple lines using matrix input to semilogy
        semilogy2 = semilogy(Data_holder,'Parent',axes2);
        set(gca, 'fontsize', 12);
        for l = 1:length(Error_Data),
            color = color_data(l,:);
            tag = ['set(semilogy1(', int2str(l),'), ''Color'', color,''DisplayName'', type.name',int2str(l),');'];
            eval(tag);
        end
%         set(semilogy1(2),'Color',[1 0 0],'DisplayName',type.name2);

        % Create xlabel
        xlabel('Image Pixels');

        % Create ylabel
        ylabel('Log Angular Error, (log degrees)');

        % Create title
        title('Angular error between the estimate and actual Flow, in ascending order');

        % Create legend
        legend2 = legend(axes2,'show');
        set(legend2,...
            'Position',[0.818692044899952 0.136476203907978 0.142752562225476 0.0828220858895706]);
           
 case 4
        % Plot visual comparison between u_GT and u_est:
        
        x = type.x;
        y = type.y;
        max_amp = type.max_amp;
        amp = type.amp;
        
        f1 = figure; set(f1,'position',get(0,'screensize')); subplot(1,2,1);  plotflow(x,y, u_GT, max_amp);
        set(gca, 'fontsize', 12);
        title(['Ground Truth Optical Flow,    Max Amplitude = ', num2str(amp)]);
        xlabel('x component of flow');
        ylabel('y component of flow');
        axis xy;

        subplot(1,2,2); plotflow(x,y, u_est, max_amp);
        set(gca, 'fontsize', 12);
        title(['Estimate of Optical Flow using ', type.title]);
        xlabel('x component of flow');
        ylabel('y component of flow');
        axis xy;

        % Create textbox
        annotation(f1,'textbox',...
            [0.17 0.133116496065155 0.47510980966325 0.0521472392638037],...
            'String',[{['Median Absolute Error = ', num2str(Error_Data.amp,4),',    Median Angular Error = ', num2str(Error_Data.angle,4), ',    PSNR = ', num2str(Error_Data.PSNR,4),'dB']};...
            {['Percent of Pixels with an Absolute Error > 0.1 = ', num2str(Error_Data.AE_10,4), '%,    Mean of the highest 10% of Absolue Errors = ', num2str(Error_Data.MAE_10P,4)]}],...
            'FontSize',14,...
            'FontName','Arial',...
            'FitBoxToText','on',...
            'HorizontalAlignment','center',...
            'BackgroundColor',[1 1 1]);
        
     case 5
    
        max_amp = type.max_amp;
        pixel_num = max(length(Error_Data(1).Abs_Error(:)), length(Error_Data(2).Abs_Error(:)));
        
        % generate figure for absolute error
        f1 = figure; set(f1,'position',get(0,'screensize'));
        axes1 = axes('Parent', f1,'YScale','log','YMinorTick','on',...
            'Position',[0.0607613469985359 0.11 0.42606149341142 0.815]);
        axis([0, pixel_num*1.05, 0, ceil(max_amp)]);
        box(axes1,'on');
        hold(axes1,'all');

        % Y-axis in log scale
        semilogy1 = semilogy([sort(Error_Data(1).Abs_Error(:)), sort(Error_Data(2).Abs_Error(:))],'Parent',axes1);
        set(gca, 'fontsize', 12);
        set(semilogy1(1),'DisplayName',type.name1);
        set(semilogy1(2),'Color',[1 0 0],'DisplayName',type.name2);

        % Create xlabel
        xlabel('Image Pixels');

        % Create ylabel
        ylabel('Log Absolute Error, log_{10}|u - u_{est}|');

        % Create title
        title('Absolue error between the estimate and actual Flow, in ascending order');

        % Create legend
        legend1 = legend(axes1,'show');
        set(legend1,...
            'Position',[0.334797462176672 0.134834167783184 0.142752562225476 0.0828220858895706]);

        % generate figure for angular error
        axes2 = axes('Parent',f1,'YScale','log','YMinorTick','on',...
            'Position',[0.550512445095168 0.11 0.420204978038067 0.815]);
        axis([0, pixel_num*1.05, 0, 100]);
        box(axes2,'on');
        hold(axes2,'all');

        % Create multiple lines using matrix input to semilogy
        semilogy2 = semilogy([sort(Error_Data(1).Angle_Error(:)), sort(Error_Data(2).Angle_Error(:))],'Parent',axes2);
        set(gca, 'fontsize', 12);
        set(semilogy2(1),'DisplayName', type.name1);
        set(semilogy2(2),'Color',[1 0 0],'DisplayName',type.name2);

        % Create xlabel
        xlabel('Image Pixels');

        % Create ylabel
        ylabel('Log Angular Error, (log degrees)');

        % Create title
        title('Angular error between the estimate and actual Flow, in ascending order');

        % Create legend
        legend2 = legend(axes2,'show');
        set(legend2,...
            'Position',[0.818692044899952 0.136476203907978 0.142752562225476 0.0828220858895706]);
end
end