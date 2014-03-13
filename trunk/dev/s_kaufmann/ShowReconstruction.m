function ShowReconstruction(IMG, ZPosition, ElectrodeOffset)
    IMG.calc_colours.npoints = 320;
    IMG.calc_colours.cb_shrink_move = [1, 1, -0.1500]; % Leitfähigkeitsscala Balken (Width, Height, Position)

    Title = regexprep(IMG.name, '_', '\\_');

    figure(); set(gcf, 'Name', 'Reconstruction Results');

    posn= [ inf,inf,ZPosition-3*ElectrodeOffset,1; ...
            inf,inf,ZPosition-ElectrodeOffset,1; ...
            inf,inf,ZPosition,1 ; ...
            inf,inf,ZPosition+ElectrodeOffset,1; ...
            inf,inf,ZPosition+3*ElectrodeOffset,1];

    subplot(1, 2, 2); show_slices(IMG, posn); title('Cut planes');
    subplot(1, 2, 1); show_fem(IMG, [1, 1]); title(Title);

    % Create textbox
    annotation(gcf,'textbox',...
        [0.8 0.169651162790697 0.4 0.0494186046511628],...
        'String',{'Electrode',['Layer + ' num2str(3*ElectrodeOffset) 'mm']},...
        'FontWeight','bold',...
        'FontSize',10,...
        'FitBoxToText','off',...
        'LineStyle','none');

    % Create textbox
    annotation(gcf,'textbox',...
        [0.8 0.821020930232558 0.4 0.0494186046511628],...
        'String',{'Electrode',['Layer - ' num2str(3*ElectrodeOffset) 'mm']},...
        'FontWeight','bold',...
        'FontSize',10,...
        'FitBoxToText','off',...
        'LineStyle','none');

    % Create textbox
    annotation(gcf,'textbox',...
        [0.8 0.492732558139535 0.4 0.0494186046511628],...
        'String',{'Electrode Layer'},...
        'FontWeight','bold',...
        'FontSize',10,...
        'FitBoxToText','off',...
        'LineStyle','none');

    % Create textbox
    annotation(gcf,'textbox',...
        [0.8 0.340813953488372 0.4 0.0494186046511627],...
        'String',{'Electrode',['Layer + ' num2str(ElectrodeOffset) 'mm']},...
        'FontWeight','bold',...
        'FontSize',10,...
        'FitBoxToText','off',...
        'LineStyle','none');

    % Create textbox
    annotation(gcf,'textbox',...
        [0.8 0.660232558139535 0.4 0.0494186046511628],...
        'String',{'Electrode',['Layer - ' num2str(ElectrodeOffset) 'mm']},...
        'FontWeight','bold',...
        'FontSize',10,...
        'FitBoxToText','off',...
        'LineStyle','none');
    
    figure(); set(gcf, 'Name', 'Reconstruction on the electrode level');
    posn= [ inf, inf, ZPosition, 1 ];
    show_slices(IMG, posn); title('Cut planes');    
    
end