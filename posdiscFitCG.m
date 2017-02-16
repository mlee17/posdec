function fit = posdiscFitCG(stimfiles)

for s = 1:length(stimfiles)
    load(stimfiles{s});
    [posDiff, contrast{s}, nTrial{s}, nChoice{s}] = getPercentChoice(task{1}{1});
    if s == 1
        stimType = stimulus.stimType;
        ecc = stimulus.eccentricity;
    else
        if ~strcmpi(stimulus.stimType,stimType) || (stimulus.eccentricity~=ecc)
            disp('stim type not matching')
            return
        end
    end
end
%Check N of contrast conditions
ntrials = cellfun('size',nTrial,1);
if max(ntrials) ==5
    ncol=3;
    if ~all(ntrials == 5)
        add = find(ntrials~=5);
        for i = 1:length(add)
            nTrial{add(i)} = [zeros(1,length(posDiff));nTrial{add(i)}];
            nChoice{add(i)} = [zeros(1,length(posDiff));nChoice{add(i)}];
        end
        contrast = contrast{find(ntrials==5,1)};
    else
        contrast = contrast{1};
    end
else % 4
    ncol=2;
    contrast = contrast{1};
end

nTrial = sum(cat(3,nTrial{:}),3);
nChoice = sum(cat(3,nChoice{:}),3);
pInt1 = nChoice./nTrial;

stimSize = stimulus.width;
nCon = length(contrast);
nrow=2;ncol=3;
xs = -6:.01:6;
figure;
for iCon = 1:nCon
    fit{iCon} = fitCG(posDiff, nTrial(iCon,:), nChoice(iCon,:));
    mu = fit{iCon}.fitparams(1);
    sigma = fit{iCon}.fitparams(2);
    lambda = fit{iCon}.fitparams(3);
    B = fit{iCon}.fitparams(4);
    cgfunc = fit{iCon}.cgfunc;
    fitVal = arrayfun(@(x) cgfunc(mu,sigma,lambda,B,x), xs);
    pfitVal = arrayfun(@(x) cgfunc(mu,sigma,lambda,B,x), posDiff);
    error = sqrt(pfitVal.*(1-pfitVal)./nTrial(iCon,:));
    
    % plot
    subplot(nrow,ncol,iCon)
    myerrorbar(posDiff, pInt1(iCon,:), 'yError', error, 'Symbol','o');
    hold on;
    plot(xs, fitVal, 'k')
    title(sprintf('%s Contrast=%0.4f (ecc=%0.1f,size=%0.0f)\nmu=%0.3f, sigma=%0.3f, lambda=%0.3f, B=%0.3f\n%s', ...
        stimType,contrast(iCon),ecc,stimSize, mu,sigma,lambda,B,mat2str(nTrial(iCon,:))));
    ylabel('Percent choices Interval 1 (%)');
    xlabel('Postition difference of targets between interval 1 and 2 (deg)');
    axis([-6 6 0 1]); box off; drawnow;
end

function [posDiff, contrast, nTrial, nChoice] = getPercentChoice(task)
posDiff = task.parameter.posDiff;
nDiff = length(posDiff);
% percent Interval 1
% n Contrast levels
contrast = task.parameter.contrast;
nCon = length(contrast);
for iCon = 1:nCon
    for iDiff = 1:nDiff
        resp{iCon}{iDiff} = task.randVars.resp(task.randVars.con == contrast(iCon) & task.randVars.diff == posDiff(iDiff));
        nTrial(iCon,iDiff) = sum(ismember(resp{iCon}{iDiff},[1 2]));
        nChoice(iCon,iDiff) = sum(resp{iCon}{iDiff} == 1);
    end
end