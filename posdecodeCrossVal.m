% load all rois into the view structure. 
% retval = posdecodeCrossVal(v, <roiName>)
% if <roiName> not defined, run for all ROIs

% combine rois across hemispheres

function retval = posdecodeCrossVal(v, roiName)

% 1. load roi tseries
% 2. get sid. import prf info from AnatDB
% 3. 

% concatInfo = viewGet(v, 'concatInfo');
framePeriod = viewGet(v, 'framePeriod');

rois 
% [stimvol stimNames var] = getStimvol(v,'location_x_contrast','taskNum=2');
rois = getInstances(v,rois,stimvol,'type','deconv')

function xxx(v, roiName, scan)

%