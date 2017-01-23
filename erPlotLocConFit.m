function erPlotLocConFit(view,overlayNum,scan,x,y,z,roi)
% erPlotLocCon.m
%
%       $Id$	
%      usage: eventRelatedPlot()
%         by: minyoung lee
%       date: 12/14/16
%    purpose: 
%

% check arguments
if ~any(nargin == [1:7])
  help erPlotLocCon
  return
end

% select the window to plot into
fignum = selectGraphWin;

% turn off menu/title etc.
set(fignum,'NumberTitle','off');
set(fignum,'Name','eventRelatedPlot LocCon');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the hemodynamic response for voxel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ER analysis
groupNum = viewGet(view,'groupnum',viewGet(view,'curGroup'));
concatInfo = viewGet(view,'concatInfo',scan,groupNum);
tr = viewGet(view,'tr', scan, groupNum);
framePeriod = viewGet(view,'framePeriod',scan,groupNum);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
contrastNames = {'25% contrast','50% contrast','100% contrast'};
locationNames = {'LVF -14','LVF -8','RVF 8','RVF 14'};

hdrlen = 25;

loc = [-14 -8 8 14];
con = [.25 .5 1];

% locstr = {'location=[-14]','location=[-8]','location=[8]','location=[14]'};
% constr ={'contrast=[0.25]','contrast=[0.5]','contrast=[1]'};

[stimvolLoc, stimNamesLoc, var] = getStimvol(view, 'location', 'taskNum=2','phaseNum=1','segmentNum=1');
[stimvolCon, stimNamesCon, var] = getStimvol(view, 'contrast', 'taskNum=2','phaseNum=1','segmentNum=1');
d = {};
for locnum = 1:length(loc)
    for connum = 1:length(con)
        stimvol{locnum}{connum} = intersect(stimvolLoc{locnum},stimvolCon{connum});
    end

    % get the time series
    tSeries = squeeze(loadTSeries(view,scan,z,[],x,y));
    d{locnum} = fitTimecourse(tSeries,stimvol{locnum},framePeriod,'concatInfo', concatInfo, ...
            'fitType=deconv','amplitudeType=fit2','displayFit=0', 'hdrlen', hdrlen);
        
    subplot(2,2,locnum)
    plotEhdr(d{locnum}.deconv.time, d{locnum}.deconv.ehdr, d{locnum}.deconv.ehdrste);
    title(sprintf('Voxel[%i %i %i] %s', x,y,z,locationNames{locnum}));
    lhandle = legend(contrastNames);
    set(lhandle,'Interpreter','none');
%    set(lhandle, 'FontSize', 8);
    box off;
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%% function to plot ehdr
%%%%%%%%%%%%%%%%%%%%%%%%%
function plotEhdr(time,ehdr,ehdrste,lineSymbol)

% whether to plot the line inbetween points or not
if ~exist('lineSymbol','var'),lineSymbol = '-';end

brewer = brewermap(6,'Blues');
colors = [brewer(4,:);brewer(5,:);brewer(6,:)];
symbols = [1,1,1];

% display ehdr
for i = 1:size(ehdr,1)
  if ieNotDefined('ehdrste')
    h=plot(time,ehdr(i,:),getsymbol(symbols(i),lineSymbol),'Color', colors(i,:),'MarkerSize',6, 'LineWidth',1.5);  %getcolor(i,getsymbol(i,lineSymbol)),'MarkerSize',8);
  else
    h=errorbar(time,ehdr(i,:),ehdrste(i,:),ehdrste(i,:),getsymbol(symbols(i),lineSymbol),'Color', colors(i,:),'MarkerSize',8,'MarkerEdgeColor','w');%getcolor(i,getsymbol(i,lineSymbol)),'MarkerSize',8);
  end
  set(h,'MarkerFaceColor',colors(i,:));
  hold on
end 

xlabel('Time (sec)');
ylabel('% Signal change');
axis tight;
yaxis(-3,3.5);
% getsymbol(symbols(i),lineSymbol)
% ['o',lineSymbol]
