% load all rois into the view structure. (don't need to....)
% retval = posdecodeCrossVal(v, <roiName>)
% if <roiName> not defined, run for all ROIs (X)

% combine rois across hemispheres

function [roi,p] = posdecodeCrossVal(v, roiName, varargin)

% 1. load roi tseries *combined roi
% 2. get sid. import prf info from AnatDB
% 3. 

preload = []; n=[]; nFold=[]; sumplot=[];
% get type of instance getting
[argNames argValues args] = getArgs(varargin,{'preload=0', 'n=100', 'nFold=5', 'sumplot=0'});

v = viewSet(v,'curGroup',3);
v = viewSet(v,'curScan',1);
roi = posdecode(v,roiName,'preload', preload,'sumplot', sumplot);

% get event related structures
[stimvol stimNames var] = getStimvol(v,'location_x_contrast','taskNum=2');

pRFcutoff = .4;
% pull out r2 values
for iVox = 1:roi.n
  roi.r2(iVox) = roi.fit(iVox).deconv.r2;
end
[r2sorted r2sortindex] = sort(roi.r2(:),1,'descend');
prfr_sorted = roi.pRF.r(r2sortindex);
% only keep voxels whose r from pRF is >= pRFcutoff
r2sorted = r2sorted(prfr_sorted>=pRFcutoff);
r2sortindex = r2sortindex(prfr_sorted>=pRFcutoff);

roi.sortindex = r2sortindex(1:n);
% deconvolution HERE (selected voxels only)-> X do amplitude fit
% single trial fit????????????????????????????????????????????

% allStimvol = cell2mat(stimvol);
% allStimvol = num2cell(allStimvol);
concatInfo = viewGet(v,'concatInfo');
framePeriod = viewGet(v,'framePeriod');
% roi.deconv = fitTimecourse(roi.tSeries(roi.sortindex,:),stimvol,framePeriod,'concatInfo',concatInfo,'fitType=deconv','amplitudeType=none');
% -> all trials just averaged/// sorted voxels.....

roi = getInstances(v,roi,stimvol);%'concatInfo',concatInfo);% 'sortindexType', 'user');
crossVal = getCrossValSets(roi{1}.classify.instances, 'nFold',nFold);


% crossVal.trainInstances 
% crossVal.testInstances
% -> remove testInstances from stimvol (separate stimvol for train and test)
for iFold = 1:nFold
    for cond = 1:8
        stimvolTrain{cond} = stimvol{cond}(crossVal.trainInstances{cond,iFold});
        stimvolTest{cond} = stimvol{cond}(crossVal.testInstances{cond,iFold});
    end
    % do deconvolutions.....
    roi{1}.deconvTrain{iFold} = fitTimecourse(roi{1}.tSeries(roi{1}.sortindex,:),stimvolTrain,framePeriod,'concatInfo',concatInfo,'fitType=deconv','amplitudeType=none');
    roi{1}.deconvTest{iFold} = fitTimecourse(roi{1}.tSeries(roi{1}.sortindex,:),stimvolTest,framePeriod,'concatInfo',concatInfo,'fitType=deconv','amplitudeType=none');

    % fit pRF to Train set (roi.deconvTrain{iFold})
    roi{1} = getModelResponse(v,roi{1},iFold);

    resid = [];
    for cond = 1:8
%         if ~isodd(cond)
            resid = [resid roi{1}.pRFModel.residual{cond}];
%         end
    end
    covMat = resid * resid';
    
    time = 0:.5:40;
    for iVox = 1:n
    for cond = 1:8
        thismodelResponse = convolveModelResponseWithHRF(roi{1}.pRFModel.modelResponse{cond}(:,iVox),roi{1}.pRFModel.hrf);
        thismodelResponse = thismodelResponse(1:41);
        thismodelResponse = thismodelResponse - mean(thismodelResponse);
        thismodelResponse = scaleAndOffset(thismodelResponse', roi{1}.pRFModel.beta{cond}(:,iVox));
        modelResponse{cond}(iVox,:) = thismodelResponse;
        
        thisehdr = squeeze(roi{1}.deconvTest{iFold}.deconv.ehdr(iVox,cond,:));
        
        mu{cond}(iVox,1) = mean(thismodelResponse(9:13));
        xval{cond}(iVox,1) = mean(thisehdr(9:13));
    end
    end
    
%     m.mu = mu; m.xval = xval;
% m.nVox = m.roi{1}.newN + m.roi{2}.newN;
% m.nVox_org = m.roi{1}.n + m.roi{2}.n;

    numloc=4; numcon=2;
    ps = repmat(1/(numloc*numcon), 1, numloc*numcon);
    for cond = 1:8
        pbsCorrect{iFold}(cond) = mvnpdf(xval{cond},mu{cond}, covMat/sqrt(length(stimvolTrain{cond})));
    
        for cc = 1:8
            pbs{iFold}(cond,cc) = mvnpdf(xval{cond}, mu{cc}, covMat/sqrt(length(stimvolTrain{cond})));
        end

    
        pb{iFold}(cond) = sum(pbs{iFold}(cond,:));

    end
% for cond = 1:8; subplot(4,2,cond); yaxis([0 max(ymax)]); end

% figure('Name','p stimulus given voxel activation pattern')
% figure(2);
    for cond = 1:8
        for cc = 1:8
            psb{iFold}(cond,cc) = exp((log(pbs{iFold}(cond,cc)) + log(ps(cc))) - log(pb{iFold}(cond)));
            psb2{iFold}(cond,cc) = exp(log(pbs{iFold}(cond,cc)) + log(ps(cc))) / pb{iFold}(cond);
        end

    end
end
    p.ps = ps;
    p.psb = zeros(8); p.pb = zeros(1,8); p.pbs = zeros(8);
    p.tmp.pb = []; p.tmp.pbs = []; p.tmp.psb = [];
for iFold = 1:nFold
    p.pbs = p.pbs + pbs{iFold};
    p.pb = p.pb + pb{iFold};
    p.psb = p.psb + psb{iFold};
    
    p.tmp.pb(:,:,iFold) = pb{iFold};
    p.tmp.pbs(:,:,iFold) = pbs{iFold};
    p.tmp.psb(:,:,iFold) = psb{iFold};
end
p.pbs = p.pbs/5;
p.pb = p.pb/5;
p.psb = p.psb/5;

p.std.pbs = std(p.tmp.pbs,0,3);
p.std.pb = std(p.tmp.pb,0,3);
p.std.psb = std(p.tmp.psb,0,3);

condNames = {'Loc1 Low contrast','Loc1 High contrast', 'Loc2 Low contrast', 'Loc2 High contrast',...
    'Loc3 Low contrast', 'Loc3 High contrast', 'Loc4 Low contrast', 'Loc4 High contrast'};
% figure('Name','p voxel activation given model')
f1 = mlrSmartfig('p voxel activation given model','reuse');clf;
f2 = mlrSmartfig('p stimulus given voxel activation pattern','reuse');clf;
for cond = 1:8
    figure(f1);
    subplot(4,2,cond)
    bar(p.pbs(cond,:))
    myerrorbar(1:8, p.pbs(cond,:),'yError',p.std.pbs(cond,:),'Symbol','o','MarkerSize',0.1);
    title(sprintf('%s: %s (N=%i/%i)', roiName, condNames{cond}, n, roi{1}.n))
    ylabel('p(b|s)')
    set(gca, 'xTickLabel', {'1L','1H','2L','2H','3L','3H','4L','4H'});
    box off
        ymax(cond) = max(get(gca,'yLim'));
end
for cond = 1:8; subplot(4,2,cond); yaxis([0 max(ymax)]); end
ymax = [];
for cond = 1:8
    
    figure(f2);
    subplot(4,2,cond)
    bar(p.psb(cond,:))
    myerrorbar(1:8, p.psb(cond,:),'yError',p.std.psb(cond,:),'Symbol','o','MarkerSize',0.1);
    title(sprintf('%s: %s (N=%i/%i)', roiName, condNames{cond}, n, roi{1}.n))
    ylabel('p(s|b)')
    set(gca, 'xTickLabel', {'1L','1H','2L','2H','3L','3H','4L','4H'});
    box off
        ymax(cond) = max(get(gca,'yLim'));
end
for cond = 1:8; subplot(4,2,cond); yaxis([0 max(ymax)]); end
ymax = [];


f3 = mlrSmartfig('p by contrast','reuse');clf;
psblow = p.psb(1:2:7, 1:2:7);
psbhigh = p.psb(2:2:8,2:2:8);
psblow_std = p.std.psb(1:2:7, 1:2:7);
psbhigh_std = p.std.psb(2:2:8,2:2:8);
brewer = brewermap(5,'*PRGn');
for cond = 1:4 % locations
    correct = cond;
    if mod(cond,2) == 1
        neighbor = cond+1;
    else
        neighbor = cond-1; 
    end
    if cond <= 2
        contra = [3 4];
    else
        contra = [2 1];
    end
    subplot(4,1,cond)
    
%     hArray = bar([psblow(cond,[correct neighbor contra]); psbhigh(cond,[correct neighbor contra])], ...
%          'grouped');
%      set(hArray(1), 'FaceColor', brewer(1,:));
%      set(hArray(2), 'FaceColor', brewer(2,:));
%      set(hArray(3:4), 'FaceColor',brewer(4,:));
%      ylabel('p(s|b)')
%      set(gca,'xTickLabel', {'Low Contrast','High Contrast'});
%      lh=legend('Correct','Neighbor','Contra-Inner','Contra-Outer');
%      set(lh,'FontSize',10,'Color','none', 'Location', 'BestOutside');
%      legend boxoff
    mybar([psblow(cond,[correct neighbor contra]); psbhigh(cond,[correct neighbor contra])],...
        'groupLabels',{'Low Contrast','High Contrast'}, ...
        'withinGroupLabels',{'Correct','Neighbor','Contra-Inner','Contra-Outer'},'yAxisMin=0','dispValues=0',...
        'withinGroupColors',{brewer(1,:) brewer(2,:) brewer(4,:) brewer(4,:)},...
        'yError', [psblow_std(cond,[correct neighbor contra]); psbhigh_std(cond,[correct neighbor contra])]);
    drawPublishAxis
    title(sprintf('%s: Loc%i (N=%i/%i)', roiName, cond, n, roi{1}.n))
    ylabel('p(s|b)');
    ymax(cond) = max(get(gca,'yLim'));
    box off
 lh=legend('Correct','Neighbor','Contra-Inner','Contra-Outer');
     set(lh,'FontSize',10,'Color','none', 'Location', 'BestOutside');
     legend boxoff

end
for cond = 1:4; subplot(4,1,cond); yaxis([0 max(ymax)]); end
    
    
    %         subplot(4,2,cond)
%         bar(psb(cond,:))
%         title(sprintf('%s: %s (N=%i/%i)',roiName, condNames{cond}, m.nVox, m.nVox_org))
%         ylabel('p(s|b)')
%         %xlabel('location X contrast')
%         set(gca, 'xTickLabel', {'1L','1H','2L','2H','3L','3H','4L','4H'});
%         box off
%         ymax(cond) = max(get(gca,'yLim'));
    
%         subplot(4,2,cond)
%         bar(pbs(cond,:))
%         title(sprintf('%s: %s (N=%i/%i)',roiName,condNames{cond}, m.nVox, m.nVox_org))
%         ylabel('p')
%         %xlabel('location X contrast')
%         set(gca, 'xTickLabel', {'1L','1H','2L','2H','3L','3H','4L','4H'});
%         box off
%     ymax(cond) = max(get(gca,'yLim'));    

% keyboard
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   convolveModelResponseWithHRF   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function modelTimecourse = convolveModelResponseWithHRF(modelTimecourse,hrf)

% n = length(modelTimecourse);
modelTimecourse = conv(modelTimecourse,hrf.hrf);
modelTimecourse = modelTimecourse(1:length(hrf.time));
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%    scaleAndOffset    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%
function modelResponse = scaleAndOffset(modelResponse,beta)

designMatrix = modelResponse;
designMatrix(:,2) = 1;

modelResponse = designMatrix*beta;
%%

% start for loop here and do pRF fit each iteration

% roi = getSortPRFFit(v, roi);

% [stimvol stimNames var] = getStimvol(v,'location_x_contrast','taskNum=2');
