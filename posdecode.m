% posdecode.m
%
%        $Id:$ 
%      usage: posdecode(v,roiName)
%         by: justin gardner
%       date: 01/04/17
%    purpose: test position decode
%
function retval = posdecode(v,roiName)

% check arguments
if ~any(nargin == [2])
  help posdecode
  return
end

preload = 1;
preloadName = sprintf('posdecode_%s',roiName);
if preload && isfile(setext(preloadName,'mat'))
  disppercent(-inf,'(posdecode) preloading structure');
  load(preloadName);
  disppercent(inf);
else
% load roi and attach pRF coords
roi = loadROITSeries(v,roiName,1,'Concatenation');
roi.subjectID = 's0389';
roi = mlrAnatDBGetPRFROI(v,roi,'noPull=1');

% get event related structures
v = viewSet(v,'curGroup',3);
v = viewSet(v,'curScan',1);
[stimvol stimNames var] = getStimvol(v,'location_x_contrast','taskNum=2');
concatInfo = viewGet(v,'concatInfo');
framePeriod = viewGet(v,'framePeriod');
roi.fit = fitTimecourse(roi.tSeries,stimvol,framePeriod,'concatInfo',concatInfo,'fitType=deconv');
save(preloadName);
end

% pull out r2 values
for iVox = 1:roi.n
  roi.r2(iVox) = roi.fit(iVox).deconv.r2;
end

mlrSmartfig('posdecode','reuse');clf;
for iCond = 1:8
  % get stimulus location from desription
  locStart = findstr('location=',stimNames{iCond});
  stimX = str2num(strtok(stimNames{iCond}(locStart+length('location='):end)));
  stimY = 0;
  stimWidth = 6;
  
  if isodd(iCond)
    color = 'k';
    rowNum = 4;colNum = ceil(iCond/2);
  else
    rowNum = 3; colNum = iCond/2;
    color = 'r';
  end
  subplot(4,4,(rowNum-1)*4+colNum);hold on
  title(sprintf('%i: %s',iCond,stimNames{iCond}));
  for iVox = 1:roi.n
    amplitude(iVox) = roi.fit(iVox).amplitude(iCond);
    x = roi.pRF.params(1,iVox);  
    y = roi.pRF.params(2,iVox);  
    sigma = roi.pRF.params(3,iVox);
    plot(x,y,'ko','MarkerSize',amplitude(iVox)*5,'MarkerFaceColor',color,'MarkerEdgeColor','none');
  end
  hline(0);
  vline(0);
  vline(-14);
  vline(-8);
  vline(14);
  vline(8);
  xaxis(-20,20);
  yaxis(-20,20);
  if isodd(iCond)
    subplot(4,4,4+colNum);hold on
    for iVox = 1:roi.n
      % calculate amount of overlap
      x = roi.pRF.params(1,iVox);  
      y = roi.pRF.params(2,iVox);  
      sigma = roi.pRF.params(3,iVox);
%      if x>0,keyboard,end
      p(iVox) = percentOverlap(x,y,sigma,stimX,stimY,stimWidth);
      p(iVox) = max(p(iVox),0.01);
      color = 'k';
      plot(x,y,'ko','MarkerSize',p(iVox)*10,'MarkerFaceColor',color,'MarkerEdgeColor','none');
    end
    xaxis(-20,20);
    yaxis(-20,20);
    title(sprintf('pRF overlap: %s',stimNames{iCond}));
  end
  subplot(4,4,colNum);hold on
  plot(p,amplitude,'k.');
  xlabel('p');
  ylabel('amplitude');
  title(sprintf('pRF/amp correlation: %s',stimNames{iCond}));
  drawnow
end

% test condition
corrVal = [];
amplitude = {};
p = {};
x = [];
y = [];

for iCond = [3 4 5 6];
  stimXvals = -23:3:23;
  stimXvals = [-8 8];
  nVals = length(stimXvals);
  nRows = floor(sqrt(nVals));
  nCols = ceil(nVals/nRows);
  mlrSmartfig(sprintf('posdecode2_%i',iCond),'reuse');clf;
  % get voxel select
  r2cutoff = 0.1;
  pRFcutoff = 0.4;
  voxelSelect = find((roi.r2>r2cutoff)  & (roi.pRF.r>pRFcutoff)');
  for iX = 1:length(stimXvals);
    subplot(nRows,nCols,iX);
    hold on
    stimX = stimXvals(iX);
    for iVoxSelect = 1:length(voxelSelect)
      iVox = voxelSelect(iVoxSelect);
      % calculate amount of overlap
      x(iVoxSelect) = roi.pRF.params(1,iVox);  
      y(iVoxSelect) = roi.pRF.params(2,iVox);  
      sigma = roi.pRF.params(3,iVox);
      p{iX}(iVoxSelect) = percentOverlap(x(iVoxSelect),y(iVoxSelect),sigma,stimX,stimY,stimWidth);
      amplitude{iCond}(iVoxSelect) = roi.fit(iVox).amplitude(iCond);
    end
    % get rid of nan p
    amplitude{iCond} = amplitude{iCond}(~isnan(p{iX}));
    x = x(~isnan(p{iX}));
    y = y(~isnan(p{iX}));
    p{iX} = p{iX}(~isnan(p{iX}));
    % correlation
    corrVal(iCond,iX) = corr(p{iX}',amplitude{iCond}');
    plot(p{iX},amplitude{iCond},'k.');
    if iX == 1
      title(sprintf('%s x=%0.1f r=%0.4f',stimNames{iCond},stimXvals(iX),corrVal(iCond,iX)));
    else
      title(sprintf('x=%0.1f r=%0.4f',stimXvals(iX),corrVal(iX)));
    end
    drawnow
  end
end

highContrastColor = [0.9 0.03 0.05];
lowContrastColor = [0.5 0.03 0.05];
mlrSmartfig('posdecode_summary');
meanCorr(1,1) = mean([corrVal(4,1) corrVal(6,2)]);
meanCorr(2,1) = mean([corrVal(4,2) corrVal(6,1)]);
meanCorr(1,2) = mean([corrVal(3,1) corrVal(5,2)]);
meanCorr(2,2) = mean([corrVal(3,2) corrVal(5,1)]);
mybar(meanCorr,'groupLabels',{'Same location','Opposite location'},'withinGroupLabels',{'High contrast','Low Contrast'},'yLabelText','Response correlation to pRF model','withinGroupColors',{highContrastColor lowContrastColor});

mlrSmartfig('posdecode_summary2','reuse');clf

subplot(3,2,1);
maxScaleFactor = 12;
iX = 1;
hold on
for iVox = 1:length(x)
  plot(x(iVox),y(iVox),'ko','MarkerFaceColor','k','MarkerEdgeColor','none','MarkerSize',maxScaleFactor*p{iX}(iVox));
end
xaxis(-12,12);yaxis(-12,12);
xaxis(-15,15);yaxis(-15,15);
hline(0);vline(0);
title('pRF Model Prediction');
drawPublishAxis;

subplot(3,2,2);
maxScaleFactor = 12;
iX = 2;
hold on
for iVox = 1:length(x)
  plot(x(iVox),y(iVox),'ko','MarkerFaceColor','k','MarkerEdgeColor','none','MarkerSize',maxScaleFactor*p{iX}(iVox));
end
xaxis(-12,12);yaxis(-12,12);
xaxis(-15,15);yaxis(-15,15);
hline(0);vline(0);
title('pRF Model Prediction');
drawPublishAxis;

subplot(3,2,3);
scaleFactor = maxScaleFactor/max(amplitude{4});
hold on
iCond = 4;
for iVox = 1:length(x)
  plot(x(iVox),y(iVox),'ko','MarkerFaceColor',highContrastColor,'MarkerEdgeColor','none','MarkerSize',scaleFactor*amplitude{iCond}(iVox));
end
xaxis(-12,12);yaxis(-12,12);
xaxis(-15,15);yaxis(-15,15);
hline(0);vline(0);
title('High contrast response');
drawPublishAxis;

subplot(3,2,5);
hold on
iCond = 3;
for iVox = 1:length(x)
  plot(x(iVox),y(iVox),'ko','MarkerFaceColor',lowContrastColor,'MarkerEdgeColor','none','MarkerSize',scaleFactor*amplitude{iCond}(iVox));
end
xaxis(-15,15);yaxis(-15,15);
hline(0);vline(0);
title('Low contrast response');
drawPublishAxis;

subplot(3,2,4);
scaleFactor = maxScaleFactor/max(amplitude{6});
hold on
iCond = 6;
for iVox = 1:length(x)
  plot(x(iVox),y(iVox),'ko','MarkerFaceColor',highContrastColor,'MarkerEdgeColor','none','MarkerSize',scaleFactor*amplitude{iCond}(iVox));
end
xaxis(-12,12);yaxis(-12,12);
xaxis(-15,15);yaxis(-15,15);
hline(0);vline(0);
title('High contrast response');
drawPublishAxis;

subplot(3,2,6);
hold on
iCond = 5;
for iVox = 1:length(x)
  plot(x(iVox),y(iVox),'ko','MarkerFaceColor',lowContrastColor,'MarkerEdgeColor','none','MarkerSize',scaleFactor*amplitude{iCond}(iVox));
end
xaxis(-15,15);yaxis(-15,15);
hline(0);vline(0);
title('Low contrast response');
drawPublishAxis;



keyboard

%%%%%%%%%%%%%%%%%%%%%%%%
%    percentOverlap    %
%%%%%%%%%%%%%%%%%%%%%%%%
function p = percentOverlap(x,y,sigma,stimX,stimY,stimSigma)

% parameters that control resolution
screenWidth = 48;
screenHeight = 48;
pixelsPerDegree = 10;

% make rf
rf = mglMakeGaussian(screenWidth,screenHeight,sigma,sigma,x,y,pixelsPerDegree,pixelsPerDegree);
%cutoff = 0.8;
%rf(rf > cutoff) = 1;
%rf(rf <= cutoff) = 0;

% make stim
stim = mglMakeGaussian(screenWidth,screenHeight,stimSigma,stimSigma,stimX,stimY,pixelsPerDegree,pixelsPerDegree);
%stim(stim > cutoff) = 1;
%stim(stim <= cutoff) = 0;

% get overlap
overlap = rf .* stim;

% calculate percent overlap
p = sum(overlap(:))/sum(rf(:));

