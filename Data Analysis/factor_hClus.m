function hClusters = factor_hClus(facCorr,numFacs,plotFig)
% Function to hierarchically cluster NNMF factors and visualize that
% clustering.

%% Parameters
if nargin == 1
    numFacs = 2;
    plotFig = true;
elseif nargin == 2
    plotFig = true;
end

%%
distance = 1-facCorr; % Compute the distance between any two factors as one minus the correlation between them.
distance = squareform(distance); % squareform needs to be used to convert the dissimilarity matrix ("distance") into the same format as the output of Matlab function pdist. This format is a vector containing all elements under the diagonal of the dissimilarity matrix.

hClusters = linkage(distance,'average');

clusThr = 0.5;

if plotFig
    
    [dendro,~,outperm] = dendrogram(hClusters,size(facCorr,1),'orientation','right');
    
    % Getting axes and figure. The dendrogram function likes to create its own
    % figure and axes, unfortunately.
    ax = gca; hold on;
    f = gcf;
    f.Position(3) = f.Position(3)*2;
    
    ax.Position = [0.3125    0.1160    0.3062    0.8150];
    ax.TickDir = 'out';
    ax.YTick = [];
    ax.XLim(1) = 0;
    ax.YLim = [ax.YLim(1)+0.5 ax.YLim(2)-0.5];
    ax.XLabel.String = 'Dissimilarity (1 - Correlation)';
    ax.FontSize = 12;
    
    ln = plot(repelem(clusThr,2),ax.YLim,'--k');
    ln.LineWidth = 1;
    
    ax2 = axes;
    imagesc(facCorr(outperm,outperm));%,'AlphaData',logical(tril(facCorr))); 
     
    axis image;
    ax2.YTick = [1:numel(outperm)];
    ax2.XTick = [1:numel(outperm)];
    ax2.YTickLabel = [];
    ax2.XTickLabel = [];
    ax2.TickDir = 'out';
    ax2.Box = 'off';
    ax2.XLabel.String = 'Factors across animals';
    ax2.YLabel.String = 'Factors across animals';
    ax2.Title.String = 'Correlation Among All Factors';
    ax2.Position = [0.0006 0.1160 0.3681 0.8150];
    ax2.YRuler.Color = 'none';
    ax2.FontSize = 12;
    ax2.XDir = 'normal';    
    
    ax.Title.String = ['Hierarchical Clustering (' num2str(numFacs) ' factors per experiment)'];
end




%% Nested functions



    function editThr_Callback(hObject, eventdata, handles)
        % hObject    handle to edit1 (see GCBO)
        % eventdata  reserved - to be defined in a future version of MATLAB
        % handles    structure with handles and user data (see GUIDATA)
        
        % Hints: get(hObject,'String') returns contents of edit1 as text
        % str2double(get(hObject,'String')) returns contents as a double
        input = str2double(get(hObject,'String'));
        if isnan(input)
            errordlg('You must enter a numeric value','Invalid Input','modal')
            uicontrol(hObject)
            return
        end
        
        clusThr = input; % Update threshold.
        ln.XData = repelem(clusThr,2); % Update line.
        
        [hClus,clusTable] = assessClus(); % Update clusters.
        tb1.Data = clusTable; % Update table.     
        
        
    end

end