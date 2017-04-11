% posdecode.m
%
%        $Id:$ 
%      usage: posdecode(v,roi,crossValSet))
%         by: minyoung lee
%       date: 01/25/17
%    purpose: 
%
function [roi] = getModelResponse(v, roi)

% pRF x y sigma
% 2 contrasts * 4 locations

% create stim images
stim.xpos = [-14 -8 8 14];
stim.ypos = 0; stimY=stim.ypos;
stim.contrast = [.25 1];
stim.width = 12;
stim.sigma = stim.width/2; stimSigma = stim.sigma;
stim.length = 2; %2 seconds
% for p = 1:length(stim.xpos)
%     stim.im{p} = exp(-(((pRFd.stimX-stim.xpos(p)).^2)/(2*(stim.std^2))+((pRFd.stimY-0).^2)/(2*(stim.std^2))));
% end

hdrlen = 20;
framePeriod = viewGet(v,'framePeriod', 1);
% cutoffr2_ER = 0.1; %cutoffr2 = viewGet(view,'overlayMin');

p.timelag = 1; p.tau=0.6; p.exponent=6; p.diffOfGamma=0;
hrf = getCanonicalHRF(p, framePeriod, hdrlen);

nVox = length(roi.sortindex);
% roin=0;
for iVox = 1:nVox
    thisVox = roi.sortindex(iVox);
    % pRF params
    x = roi.pRF.params(1,thisVox);
    y = roi.pRF.params(2,thisVox);
    sigma = roi.pRF.params(3,thisVox);

    % only take voxels with intact prf params (already done in posdecodeCrossVal)
%     if ~any(isnan(roi.pRF.params(:,voxnum))) && roi.pRF.r(thisVox) >= .4 %&& (erd.r2(x,y,s) >=cutoffr2_ER) && (roi.pRF.r(voxnum)^2>=0.1)
%         roin = roin+1;
%         roi.newScanCoords(:,roin) = [x;y;s];
        
%     rfModel = exp(-(((pRFd.stimX-rfX).^2)/(2*(rfHalfWidth^2))+((pRFd.stimY-rfY).^2)/(2*(rfHalfWidth^2))));
    % init model response
    thisModelResponse = [];
    for locnum = 1:length(stim.xpos)
        for connum = 1:length(stim.contrast)
            stimX = stim.xpos(locnum);
%             if ~any(isnan(roi.pRF.params(:,voxnum)))
            thisModelResponse = convolveModelWithStimulus(x,y,sigma,stimX,stimY,stimSigma);
            modelResponse{locnum}{connum}(:,iVox) = repmat(thisModelResponse,stim.length/framePeriod,1);
            
            modelRespConvolved = convolveModelResponseWithHRF(modelResponse{locnum}{connum}(:,iVox), hrf);
            modelRespConvolved = modelRespConvolved(1:(hdrlen/framePeriod)+1);
            %removing the mean
            modelRespConvolved = modelRespConvolved - mean(modelRespConvolved);
%             % apply concat filtering
%             if isfield(fitParams,'applyFiltering') && fitParams.applyFiltering
%                 thisModelResponse = applyConcatFiltering(thisModelResponse,fitParams.concatInfo,i);
%             else
%             % with no filtering, just remove mean
%                 thisModelResponse = thisModelResponse - mean(thisModelResponse);
%             end
             thisehdr = roi.fit(thisVox).deconv.ehdr(length(stim.contrast)*(locnum-1)+connum,:);%ehdr(length(stim.contrast)*(locnum-1)+connum,:);
            [modelRespConvolved, residual{locnum}{connum}(iVox,:), beta{locnum}{connum}(:,iVox)] ...
                = scaleAndOffset(modelRespConvolved', thisehdr');
%             else
%                 beta{h}{locnum}{connum}(:,roin) = NaN;
%                 residual{h}{locnum}{connum}(roin,:) = NaN;
%             end
            
        end
    end
%     end
 
end
% roi.newN = roin;

%saving vars into a structure
m.beta = beta;
m.modelResponse = modelResponse;
m.residual = residual;
m.hrf = hrf;
m.framePeriod = framePeriod;
m.stim = stim;
roi.pRFbetaSort=m;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   convolveModelWithStimulus   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function modelResponse = convolveModelWithStimulus(x,y,sigma,stimX,stimY,stimSigma)
% parameters that control resolution
screenWidth = 48;
screenHeight = 48;
pixelsPerDegree = 10;
rfModel = mglMakeGaussian(screenWidth,screenHeight,sigma,sigma,x,y,pixelsPerDegree,pixelsPerDegree);
stim = mglMakeGaussian(screenWidth,screenHeight,stimSigma,stimSigma,stimX,stimY,pixelsPerDegree,pixelsPerDegree);

  % multipy the stimulus frame by frame with the rfModel
  % and take the sum
    modelResponse = sum(sum(rfModel.*stim));

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%    scaleAndOffset    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [modelResponse residual beta] = scaleAndOffset(modelResponse,tSeries)

designMatrix = modelResponse;
designMatrix(:,2) = 1;

% get beta weight for the modelResponse
if ~any(isnan(modelResponse))
  beta = pinv(designMatrix)*tSeries;
  beta(1) = max(beta(1),0);
  modelResponse = designMatrix*beta;
  residual = tSeries-modelResponse;
else
  residual = tSeries;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   convolveModelResponseWithHRF   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function modelTimecourse = convolveModelResponseWithHRF(modelTimecourse,hrf)

% n = length(modelTimecourse);
modelTimecourse = conv(modelTimecourse,hrf.hrf);
modelTimecourse = modelTimecourse(1:length(hrf.time));

%%%%%%%%%%%%%%%%%%%%%%%%%
%%   getCanonicalHRF   %%
%%%%%%%%%%%%%%%%%%%%%%%%%
function hrf = getCanonicalHRF(params,sampleRate, hdrlen)

hrf.time = 0:sampleRate:hdrlen;
hrf.hrf = getGammaHRF(hrf.time,params);

% normalize to amplitude of 1
hrf.hrf = hrf.hrf / max(hrf.hrf);
%%%%%%%%%%%%%%%%%%%%%
%%   getGammaHRF   %%
%%%%%%%%%%%%%%%%%%%%%
function fun = getGammaHRF(time,p)
if ~isfield(p,'offset'); p.offset = 0;end
if ~isfield(p,'offset2'); p.offset2 = 0;end
fun = thisGamma(time,1,p.timelag,p.offset,p.tau,p.exponent)/100;
% add second gamma if this is a difference of gammas fit
if p.diffOfGamma
  fun = fun - thisGamma(time,p.amplitudeRatio,p.timelag2,p.offset2,p.tau2,p.exponent2)/100;
end

%%%%%%%%%%%%%%%%%%%
%%   thisGamma   %%
%%%%%%%%%%%%%%%%%%%
function gammafun = thisGamma(time,amplitude,timelag,offset,tau,exponent)

exponent = round(exponent);
% gamma function
gammafun = (((time-timelag)/tau).^(exponent-1).*exp(-(time-timelag)/tau))./(tau*factorial(exponent-1));
% gammafun = (((hrf.time-2)/1.2).^(6-1).*exp(-(hrf.time-2)/1.2))./(1.2*factorial(6-1));

% negative values of time are set to zero,
% so that the function always starts at zero
gammafun(find((time-timelag) < 0)) = 0;

% normalize the amplitude
if (max(gammafun)-min(gammafun))~=0
  gammafun = (gammafun-min(gammafun)) ./ (max(gammafun)-min(gammafun));
end
gammafun = (amplitude*gammafun+offset);


%%%%%%%%%%%%%%%%%%%%%%
%%   getThisGamma   %%
%%%%%%%%%%%%%%%%%%%%%%
function fun = getThisGamma(time,p,i)

fun = thisGamma(time,p.amplitude(i),p.timelag(i),p.offset(i),p.tau(i),p.exponent(i))/100;
% add second gamma if this is a difference of gammas fit
if p.diffOfGammaFit
  fun = fun + thisGamma(time,p.amplitude2(i),p.timelag2(i),p.offset2(i),p.tau2(i),p.exponent2(i))/100;
end
% %%%%%%%%%%%%%%%%%%%
% %%   thisGamma   %%
% %%%%%%%%%%%%%%%%%%%
% function gammafun = thisGamma(time,amplitude,timelag,offset,tau,exponent)
% 
% exponent = round(exponent);
% % gamma function
% gammafun = (((time-timelag)/tau).^(exponent-1).*exp(-(time-timelag)/tau))./(tau*factorial(exponent-1));
% 
% % negative values of time are set to zero,
% % so that the function always starts at zero
% gammafun(find((time-timelag) < 0)) = 0;
% 
% % normalize the amplitude
% if (max(gammafun)-min(gammafun))~=0
%   gammafun = (gammafun-min(gammafun)) ./ (max(gammafun)-min(gammafun));
% end
% gammafun = (amplitude*gammafun+offset);
