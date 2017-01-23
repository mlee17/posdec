function m = estimatePRFBetaFromER(v, roiName)
scan = viewGet(v,'curScan');
hemi = {'l','r'};
for r = 1:2
    roiBoth{r} = loadROITSeries(v,[hemi{r},roiName],1,'Concatenation', 'loadType=none');
    roiBoth{r} = mlrAnatDBGetPRFROI(v,roiBoth{r});
end
% pRF x y sigma
% 2 contrasts * 4 locations

sid = viewGet(v, 'subject');
% get pRF fit params
pRF = mlrAnatDBGetPRF(sid);
fitParams = pRF.params.pRFFit;
pRFd = pRF.d{1};

% create stim images
stim.xpos = [-14 -8 8 14];
stim.contrast = [.25 1];
stim.width = 12;
stim.std = stim.width/2;
stim.length = 2; %2 seconds
for p = 1:length(stim.xpos)
    stim.im{p} = exp(-(((pRFd.stimX-stim.xpos(p)).^2)/(2*(stim.std^2))+((pRFd.stimY-0).^2)/(2*(stim.std^2))));
end

analysis = viewGet(v, 'analysis');
if isempty(analysis)
    disp('(estimatePRFBetaFromER) Cannot find an analysis in the view. Select a file to load');
    v = loadAnalysis(v);
    analysis = viewGet(v, 'analysis');
end
if ~isfield(analysis,'d') || (length(analysis.d) < scan) || isempty(analysis.d)
  disp(sprintf('(estimatePRFBetaFromER) Event related not for scan %i',scan));
  return
end
erd = analysis.d{scan};
if isempty(erd)
  mrWarnDlg(sprintf('(estimatePRFBetaFromER) Could not find d structure for scan %i. Has eventRelated been run for this scan?',scan));
  return
end
erd.r2 = analysis.overlays(1).data{scan};

hdrlen = 25;
framePeriod = viewGet(v,'framePeriod', 1);
cutoffr2_ER = 0.1; %cutoffr2 = viewGet(view,'overlayMin');

hrf = getCanonicalHRF(fitParams, framePeriod, hdrlen);
for h = 1:length(roiBoth)
    roin=0;
    roi = roiBoth{h};
for voxnum = 1:roi.n
    x = roi.scanCoords(1,voxnum);
    y = roi.scanCoords(2,voxnum);
    s = roi.scanCoords(3,voxnum);
    % only take voxels with intact prf params 
    if ~any(isnan(roi.pRF.params(:,voxnum))) && (erd.r2(x,y,s) >=cutoffr2_ER) && (roi.pRF.r(voxnum)^2>=0.1)
        roin = roin+1;
        roi.newScanCoords(:,roin) = [x;y;s];
        
        [ehdr, time] = gethdr(erd,x,y,s);
  
    rfX = roi.pRF.params(1,voxnum);
    rfY = roi.pRF.params(2,voxnum);
    rfHalfWidth = roi.pRF.params(3,voxnum);
    rfModel = exp(-(((pRFd.stimX-rfX).^2)/(2*(rfHalfWidth^2))+((pRFd.stimY-rfY).^2)/(2*(rfHalfWidth^2))));
    % init model response
    thisModelResponse = [];
    for locnum = 1:length(stim.xpos)
        for connum = 1:length(stim.contrast)
            if ~any(isnan(roi.pRF.params(:,voxnum)))
            thisModelResponse = convolveModelWithStimulus(rfModel, stim.im{locnum} * stim.contrast(connum));
            modelResponse{h}{locnum}{connum}(:,roin) = repmat(thisModelResponse,stim.length/framePeriod,1);
            
            modelRespConvolved = convolveModelResponseWithHRF(modelResponse{h}{locnum}{connum}(:,roin), hrf);
            modelRespConvolved = modelRespConvolved(1:hdrlen/framePeriod);
            %removing the mean
            modelRespConvolved = modelRespConvolved - mean(modelRespConvolved);
%             % apply concat filtering
%             if isfield(fitParams,'applyFiltering') && fitParams.applyFiltering
%                 thisModelResponse = applyConcatFiltering(thisModelResponse,fitParams.concatInfo,i);
%             else
%             % with no filtering, just remove mean
%                 thisModelResponse = thisModelResponse - mean(thisModelResponse);
%             end
             thisehdr = ehdr(length(stim.contrast)*(locnum-1)+connum,:);
            [modelRespConvolved, residual{h}{locnum}{connum}(roin,:), beta{h}{locnum}{connum}(:,roin)] ...
                = scaleAndOffset(modelRespConvolved', thisehdr');
            else
                beta{h}{locnum}{connum}(:,roin) = NaN;
                residual{h}{locnum}{connum}(roin,:) = NaN;
            end
            
        end
    end
    end
 
end
roi.newN = roin;
roiBoth{h} = roi;
end

%saving vars into a structure
m.beta = beta;
m.modelResponse = modelResponse;
m.residual = residual;
m.hrf = hrf;
m.framePeriod = framePeriod;
m.stim = stim;
m.roi = roiBoth;
m.erd = erd;

%%%%%%%%%%%%%%%%%%%%%%%%
%    scaleAndOffset    %
%%%%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   convolveModelWithStimulus   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function modelResponse = convolveModelWithStimulus(rfModel,stim)

  % multipy the stimulus frame by frame with the rfModel
  % and take the sum
    modelResponse = sum(sum(rfModel.*stim));
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