function fit = posdiscFitCG(stimfiles)

for s = 1:length(stimfiles)
    load(stimfiles{s});
    [posDiff, contrast, nTrial{s}, nChoice{s}] = getPercentChoice(task{1}{1});
end
nTrial = sum(cat(3,nTrial{:}),3);
nChoice = sum(cat(3,nChoice{:}),3);
pInt1 = nChoice./nTrial;

stimType = stimulus.stimType;
stimSize = stimulus.width;
ecc = stimulus.eccentricity;
nCon = length(contrast);

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
    subplot(2,3,iCon)
    myerrorbar(posDiff, pInt1(iCon,:), 'yError', error, 'Symbol','o');
    hold on;
    plot(xs, fitVal, 'k')
    title(sprintf('%s Contrast=%0.4f (ecc=%0.1f,size=%0.0f)\nmu=%0.3f, sigma=%0.3f, lambda=%0.3f, B=%0.3f\n', ...
        stimType,contrast(iCon),ecc,stimSize, mu,sigma,lambda,B));
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