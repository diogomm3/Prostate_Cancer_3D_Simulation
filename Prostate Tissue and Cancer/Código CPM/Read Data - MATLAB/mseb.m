function H = mseb(x,y,errBar,lineProps,transparent)
% Inputs
% x - vector of x values [optional, can be left empty]
% y - vector of y values or a matrix of C x N, where C is the number of
%     lines to be plotted and N is the number of samples in each line and
%     should be equal to length(X)
% errBar - if a vector we draw symmetric errorbars. If it has a
%          size of [C,length(x),2] then asymmetric error bars are drawn,
%          with row 1 being the upper bar and row 2 being the lower
%          bar. In the present version errBar does not support two function
%          handles.
% lineProps  - [optional. Can also be set as "[]" for default settings].
%              Struct containing fields that define lineproperties for the 
%              plot function. It is possible to only define some of the
%              fields.
%       .col - cell array, where each element defines the colour of each
%              line. This can be done using either strings or three-element
%              RGB vectors e.g., either 'b' or [0 0 1] for blue.
%     .style - linestyle of the lines from y-data. Default is '-'.
%     .width - linewidth of the lines from y-data. Default is 2.
% .edgestyle - linestyle of edges that are overlapped by errorbars from
%              other lines. 
%             
% transparent - [optional, 0 by default] if ==1 the shaded error
%               bar is made transparent, which forces the renderer
%               to be openGl. However, if this is saved as .eps the
%               resulting file will contain a raster not a vector
%               image. openGL does not support having logarithmic axes.
%
% Outputs
% H - structure with an element for each line entry containing handles to 
%     the generated plot objects (e.g. H(c) contains the handles to the
%     c'th line entry. 
%
%
% mseb([],y_mean,y_std,lineProps);ylim([-50 150])
error(nargchk(3,5,nargin))
[C, N] = size(y);
if N==1
    C = 1;
    N = length(y);
    y = y';
end
% Cheking the x data
if isempty(x)
    x=repmat(1:N,[C,1]);
elseif length(x(:))==N && C>1
    x=repmat(x(:)',[C,1]);
end
if (size(x,1) ~= C) || (size(x,2)~=N)
    error('inputs x and y do not have same dimensions')
end
if (size(errBar,1) ~= C) || (size(errBar,2)~=N)
    if size(errBar,1)==N && size(errBar,2)==2 && size(errBar,3)~=1 && C==1
        errBar = permute(errBar,[3,1,2]);
    elseif size(errBar,1)==N && size(errBar,2)==1 && size(errBar,3)~=1 && C==1
        errBar=repmat(errBar,[1,2]);
        errBar = permute(errBar,[3,1,2]);
    else
        error('inputs errBar and y do not have same dimensions')
    end
end
if size(errBar,3)==1
    errBar = repmat(errBar,[1,1,2]);
end
if nargin<4
    lineProps = [];
end
if nargin<5 || ~isnumeric(transparent)
    transparent=0;
end
if ~isfield(lineProps,'col')
    for c = 1:C
        colours = 'brkgmcy';
        col_ind = rem(c-1,length(colours))+1;
        lineProps.col{c} = colours(col_ind);
    end
end
if ~isfield(lineProps,'style'),	lineProps.style = '-'; end
if ~isfield(lineProps,'width'),	lineProps.width = 2; end
if ~isfield(lineProps,'edgestyle'),	lineProps.edgestyle = '--'; end
holdStatus=ishold;
if ~holdStatus, hold on,  end
for c = C:-1:1
    H(c).mainLine=plot(x(c,:),y(c,:),'color',lineProps.col{c});
    
    
    
    col=get(H(c).mainLine,'color');
    edgeColor=col+(1-col)*0.55;
    patchSaturation=0.15; 
    if transparent
        faceAlpha=patchSaturation;
        patchColor=col;
        set(gcf,'renderer','openGL')
    else
        faceAlpha=1;
        patchColor=col+(1-col)*(1-patchSaturation);
        set(gcf,'renderer','painters')
    end
    
    uE=y(c,:)+errBar(c,:,1);
    lE=y(c,:)-errBar(c,:,2);
    
    yP=[lE,fliplr(uE)];
    xP=[x(c,:),fliplr(x(c,:))];
    
    xP(isnan(yP))=[];
    yP(isnan(yP))=[];
    
    
    H(c).patch=patch(xP,yP,1,'facecolor',patchColor,...
        'edgecolor','none',...
        'facealpha',faceAlpha);
    
    
    H(c).edge(1)=plot(x(c,:),lE,'-','color',edgeColor);
    H(c).edge(2)=plot(x(c,:),uE,'-','color',edgeColor);
end
for c = C:-1:1
    col=get(H(c).mainLine,'color');
    edgeColor=col+(1-col)*0.55;
    lE = get(H(c).edge(1), 'ydata');
    uE = get(H(c).edge(2), 'ydata');
    
    H(c).edgeoverlap(1)=plot(x(c,:),lE,lineProps.edgestyle,'color',edgeColor);
    H(c).edgeoverlap(2)=plot(x(c,:),uE,lineProps.edgestyle,'color',edgeColor);
    
end
for c = 1:C
    delete(H(c).mainLine)
    H(c).mainLine=plot(x(c,:),y(c,:),lineProps.style,'color',lineProps.col{c},...
        'linewidth',lineProps.width);
    
    
    set(get(get(H(c).patch,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
    set(get(get(H(c).edge(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
    set(get(get(H(c).edge(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
    set(get(get(H(c).edgeoverlap(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
    set(get(get(H(c).edgeoverlap(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
    set(get(get(H(c).mainLine,'Annotation'),'LegendInformation'),'IconDisplayStyle','on')
    
end
if ~holdStatus, hold off, end
end