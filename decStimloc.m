% remove NaNs -> new roi Coords
% Low/High separately
function [p,m] = decStimloc(v, roiName, m)
hemi = {'l','r'};
if ieNotDefined('v')
    v = newView;
end
v = viewSet(v,'currentGroup', 'Concatenation');
roiList = viewGet(v, 'roiNames');
roisToLoad = cell(1,2);
for r = 1:length(hemi)
    if ~any(strcmp([hemi{r},roiName],roiList))
        roisToLoad{r} = [hemi{r},roiName];
    end
end
if ~any(cellfun('isempty',roisToLoad))
    roisToLoad = roisToLoad(~cellfun('isempty',roisToLoad));
    % get subject ID
    sid = viewGet(v, 'subject');
    subjectID = mlrAnatDBSubjectID(sid);
    if isempty(subjectID), return, end
    % get the local repos
    localRepoSubject = mlrAnatDBGetRepo(subjectID);
    if isempty(localRepoSubject), return, end
    % loadROI
    v = loadROI(v, roisToLoad, [], fullfile(localRepoSubject,'mlrROIs'));
end
if ieNotDefined('m') || ~strcmp(m.roi{1}.name, [hemi{1},roiName])
    m = estimatePRFBetaFromER(v, roiName);
end
% 
% % find voxels without pRF info & remove
% lhRemove = find(isnan(m.beta{1}{1}{1}(1,:)));
% rhRemove = find(isnan(m.beta{2}{1}{1}(1,:)));
% keep{1} = find(~isnan(m.beta{1}{1}{1}(1,:)));
% keep{2} = find(~isnan(m.beta{2}{1}{1}(1,:)));
% 
% lhNewCoords = m.roi{1}.scanCoords(:,keep{1});
% rhNewCoords = m.roi{2}.scanCoords(:,keep{2});

scanCoordsAll = [m.roi{1}.newScanCoords m.roi{2}.newScanCoords];

residTemp = cell(1,2);
for h = 1:length(hemi)
    for loc = 1:length(m.stim.xpos)
        for con = 1:length(m.stim.contrast)
            if mod(con,2) ==0
            residTemp{h} = [residTemp{h} m.residual{h}{loc}{con}];
            end
        end
    end
end
numcon = length(m.stim.contrast); numloc = length(m.stim.xpos);
for loc = 1:numloc
    for con = 1:numcon
        modelRespPrepAll{numcon*(loc-1)+con} = [m.modelResponse{1}{loc}{con} m.modelResponse{2}{loc}{con}];
        betaAll{numcon*(loc-1)+con} = [m.beta{1}{loc}{con} m.beta{2}{loc}{con}];
    end
end
% create covariance matrix
aggregateResidual = [residTemp{1};residTemp{2}];
% for s = 1:size(aggregateResidual,1)
%     aggregateResidual(s,:) = aggregateResidual(s,:) - mean(aggregateResidual(s,:));
% end
covMat = aggregateResidual * aggregateResidual';
% low con =1 , high con = 2
% for each actual stimulus location/contrast, compute 8 (4 loc * 2con)
% probability values
for voxnum = 1:size(scanCoordsAll,2)
    for cond = 1:numloc*numcon
        thismodelResponse = convolveModelResponseWithHRF(modelRespPrepAll{cond}(:,voxnum),m.hrf);
        thismodelResponse = thismodelResponse(1:50);
        thismodelResponse = thismodelResponse - mean(thismodelResponse);
        thismodelResponse = scaleAndOffset(thismodelResponse', betaAll{cond}(:,voxnum));
        modelResponse{cond}(voxnum,:) = thismodelResponse;
                
        x = scanCoordsAll(1,voxnum);
        y = scanCoordsAll(2,voxnum);
        s = scanCoordsAll(3,voxnum);
        [ehdr, time] = gethdr(m.erd,x,y,s);
        
        mu{cond}(voxnum,1) = mean(thismodelResponse(9:13));
        xval{cond}(voxnum,1) = mean(ehdr(cond,9:13));
    end
end
m.mu = mu; m.xval = xval;
m.nVox = m.roi{1}.newN + m.roi{2}.newN;
m.nVox_org = m.roi{1}.n + m.roi{2}.n;
condNames = {'Loc1 Low contrast','Loc1 High contrast', 'Loc2 Low contrast', 'Loc2 High contrast',...
    'Loc3 Low contrast', 'Loc3 High contrast', 'Loc4 Low contrast', 'Loc4 High contrast'};
% figure('Name','p voxel activation given model')
figure(1)
ps = repmat(1/(numloc*numcon), 1, numloc*numcon);
for cond = 1:8
    pbsCorrect(cond) = mvnpdf(xval{cond},mu{cond}, covMat/sqrt(length(m.erd.stimvol{cond})));
    
    for cc = 1:8
        pbs(cond,cc) = mvnpdf(xval{cond}, mu{cc}, covMat/sqrt(length(m.erd.stimvol{cond})));
    end
    subplot(4,2,cond)
    bar(pbs(cond,:))
    title(sprintf('%s: %s (N=%i/%i)',roiName,condNames{cond}, m.nVox, m.nVox_org))
    ylabel('p')
    %xlabel('location X contrast')
    set(gca, 'xTickLabel', {'1L','1H','2L','2H','3L','3H','4L','4H'});
    box off
    
    pb(cond) = sum(pbs(cond,:));
    ymax(cond) = max(get(gca,'yLim'));
end
for cond = 1:8; subplot(4,2,cond); yaxis([0 max(ymax)]); end

% figure('Name','p stimulus given voxel activation pattern')
figure(2);
for cond = 1:8
    for cc = 1:8
        psb(cond,cc) = exp((log(pbs(cond,cc)) + log(ps(cc))) - log(pb(cond)));
        psb2(cond,cc) = exp(log(pbs(cond,cc)) + log(ps(cc))) / pb(cond);
    end
    subplot(4,2,cond)
    bar(psb(cond,:))
    title(sprintf('%s: %s (N=%i/%i)',roiName, condNames{cond}, m.nVox, m.nVox_org))
    ylabel('p(s|b)')
    %xlabel('location X contrast')
    set(gca, 'xTickLabel', {'1L','1H','2L','2H','3L','3H','4L','4H'});
    box off
    ymax(cond) = max(get(gca,'yLim'));
end
for cond = 1:8; subplot(4,2,cond); yaxis([0 max(ymax)]); end
ymax = [];

% figure('Name','p by contrast')
figure(3);
psblow = psb(1:2:7, 1:2:7);
psbhigh = psb(2:2:8,2:2:8);
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
    
    hArray = bar([psblow(cond,[correct neighbor contra]); psbhigh(cond,[correct neighbor contra])], ...
         'grouped');
     set(hArray(1), 'FaceColor', brewer(1,:));
     set(hArray(2), 'FaceColor', brewer(2,:));
     set(hArray(3:4), 'FaceColor',brewer(4,:));
     ylabel('p(s|b)')
     set(gca,'xTickLabel', {'Low Contrast','High Contrast'});
     lh=legend('Correct','Neighbor','Contra-Inner','Contra-Outer');
     set(lh,'FontSize',10,'Color','none', 'Location', 'BestOutside');
     legend boxoff
%     mybar([psblow(cond,[correct neighbor contra]); psbhigh(cond,[correct neighbor contra])],...
%         'groupLabels',{'Low Contrast','High Contrast'}, ...
%         'withinGroupLabels',{'Correct','Neighbor','Contra-Inner','Contra-Outer'},'yAxisMin=0','dispValues=0',...
%         'withinGroupColors',{brewer(1,:) brewer(2,:) brewer(4,:) brewer(4,:)});
%     drawPublishAxis
    title(sprintf('%s: Loc%i (N=%i/%i)', roiName,cond, m.nVox, m.nVox_org))
    ylabel('p(s|b)');
    ymax(cond) = max(get(gca,'yLim'));
    box off
%     % low contrast
%     subplot(4,2,2*cond-1)
%     bar(psblow(cond,[correct neighbor contra]))
%     title(sprintf('Loc%i Low contrast (N=%i/%i)', cond, m.nVox, m.nVox_org))
%     set(gca, 'xTickLabel', {'correct', 'neighbor', 'contra_in', 'contra_out'});
%     box off
%     % high contrast
%     subplot(4,2,2*cond)
%     bar(psbhigh(cond,[correct neighbor contra]))
%     title(sprintf('Loc%i High contrast (N=%i/%i)', cond, m.nVox, m.nVox_org))
%     set(gca, 'xTickLabel', {'correct', 'neighbor', 'contra_in', 'contra_out'});
%     box off
end
for cond = 1:4; subplot(4,1,cond); yaxis([0 max(ymax)]); end

    p.ps = ps; p.pbsCorrect = pbsCorrect; p.pbs = pbs; p.pb = pb;
    p.psb = psb; p.psb2 = psb2; p.psblow = psblow; p.psbhigh = psbhigh;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   convolveModelResponseWithHRF   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function modelTimecourse = convolveModelResponseWithHRF(modelTimecourse,hrf)

% n = length(modelTimecourse);
modelTimecourse = conv(modelTimecourse,hrf.hrf);
modelTimecourse = modelTimecourse(1:length(hrf.time));

function modelResponse = scaleAndOffset(modelResponse,beta)

designMatrix = modelResponse;
designMatrix(:,2) = 1;

modelResponse = designMatrix*beta;
