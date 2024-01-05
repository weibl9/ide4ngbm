function gboxplot(MatCel,xlab,colrs,lxt,lcn)
% MatCel: a cell array including all matrices. All matrices should have
%         equal number of columns. 
% xlab:   the x-label ticks in each groups.
% colrs:  matrix specifies the color of boxes.
% lxt:    text specifies the contents of legend.
% lcn:    text specifies the location of legend. 
% 
% 
% Example: 
% rng default
% ape_rm= randn(500,3);
% ape_gm = 5*randn(500,3); 
% figure;
% gbp = {ape_rm, ape_gm};
% xlabel = 1:size(gbp{1},2);
% colr = [1 0 0 .5; 0 1 0 .5 ]';
% lgd = {'reduced model', 'grey model'};
% gboxplot(gbp, xlabel, colr, lgd,'northeast')
% gboxplot(gbp, xlabel, colr)

%% Check inputs
if ~iscell(MatCel)
    error('Input data is not even a cell array!');
elseif min(size(MatCel)) ~= 1
    error('Input should be a row or column array!');
end

ngs = size(MatCel{1},2);        % number of group boxes
nbs = length(MatCel);           % number of boxes in each group

if size(colrs,2)~=nbs
    error('Wrong amount of colors!');
end

if length(xlab)~=ngs
    error('Wrong amount of X labels given');
end

%% Check legends 
if nargin == 5 
    lgds = 1; % 'lxt' and 'lcn' included
else 
    lgds = 0;
end    

%% Global configuratuion
box_width = 0.15;                % widths of the individual box
outlier_size = 5;               % size of outliers in each plot

%%

% Prepare data 
gdata = cell(ngs,nbs); 
for i = 1:ngs
    for j = 1:nbs
        gdata{i,j} = MatCel{j}(:,i);
    end
end 

% Calculate positions of boxes
pos = 1:0.25:nbs*ngs*0.25+1+0.25*ngs;
pos(1:nbs+1:end) = [];

% Extract data and label it in the group
bxv = [];
bxg = [];
for ii=1:ngs
    for jj=1:nbs
        aux = gdata{ii,jj};
        bxv = vertcat(bxv,aux(:)); %#ok<*AGROW>
        bxg = vertcat(bxg,ones(size(aux(:)))*jj+(ii-1)*nbs); %#ok<*AGROW>
    end
end

% Plot group boxs
boxplot(bxv, bxg, 'positions',pos, ... 
    "Colors","w", "Notch","on", "Widths",box_width, ...  % configuration of box body
    "Symbol","w+", "OutlierSize",outlier_size)  % configuration of outliers 

% Fill boxes with colors
color = repmat(colrs, 1, ngs);

h = findobj(gca,'Tag','Box');
for jj=1:length(h)
   patch(get(h(jj),'XData'),get(h(jj),'YData'),color(1:3,jj)',...
       'FaceAlpha',color(4,jj), 'EdgeColor','none');
end


% locate tehe legends
if lgds
    legend(fliplr(lxt),'location',lcn,'fontsize',12);
end

%
hold on
boxplot(bxv, bxg, 'positions',pos, ... 
    "Colors","k", "Notch","on", "Widths",box_width, ...  % configuration of box body
    "Symbol","r+", "OutlierSize",outlier_size)  % configuration of outliers 
hold off

% Set x-labels
aux = reshape(pos,nbs,[]);
xlabpos = sum(aux,1)./nbs;
set(gca, 'xtick',xlabpos, 'xticklabel',xlab);


end

