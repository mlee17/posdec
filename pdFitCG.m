function [fit, fieldNames] = pdFitCG(sbjID)
	sbjDir = sprintf('~/data/posdisc/%s', sbjID);
	files = dir(sbjDir);
	name = {files.name}; 
	name = cellfun(@(x) x(1), name);
	hid = strfind(name,'.');
	isDir = find([files.isdir]);
	rm = union(hid, isDir);
	files(rm)=[];
	nStim = length(files);
	stimlist = struct('taskType', [],'eccentricity', [], 'contrast', [],'stimType', {}, 'task',{}, 'stimulus',{}, 'stimname',{});
	taskName = {'stairCon','fixCon', 'fixEcc', 'constant'};

	for i = 1:nStim
		clearvars task stimulus myscreen
		stimname = files(i).name;
		load([sbjDir,'/',stimname]);
		stimlist(i).taskType = stimulus.taskType;
		if stimulus.fixOrigin(1) == -5
			stimlist(i).eccentricity = stimulus.eccentricity-stimulus.fixOrigin(1);
		else 
			stimlist(i).eccentricity = stimulus.eccentricity;
		end
		stimlist(i).contrast = stimulus.contrast;
		stimlist(i).stimType = stimulus.stimType;
		stimlist(i).task = task{1}{1};
		stimlist(i).stimulus = stimulus;
		stimlist(i).stimname = stimname;
	end 

	[stimlist(find([stimlist.taskType] == 1)).contrast] = deal(NaN);

	stim = struct([]);
	for i = 1:nStim
		switch stimlist(i).taskType
		case 1
			lab = sprintf('%s_Ecc%02d', taskName{1}, stimlist(i).eccentricity);
			if isfield(stim, lab)
					stim(size(stim,1)+1).(lab) = stimlist(i);
			else 
				stim(1).(lab) = stimlist(i);
			end
		case {2,3}
			lab = sprintf('Ecc%d_Con%d', stimlist(i).eccentricity, stimlist(i).contrast*1000);
			if isfield(stim, lab)
					stim(size(stim,1)+1).(lab) = stimlist(i);
			else 
				stim(1).(lab) = stimlist(i);
			end
		end
	end

	fieldNames = fieldnames(stim);
	fieldNames = sort(fieldNames);
	nrow = 4;
	ncol = ceil(length(fieldNames)/nrow);
	stairNames = {'con','dist','dist'};
	figure(2);figure(1);
	for i = 1:length(fieldNames)
		if i==7; figure(2); end
		numSession = length([stim.(fieldNames{i})]);
		stimStr = [stim.(fieldNames{i})];
		Task = cell(1,numSession);
		taskType = stimStr(1).taskType;
		respCombined = [];
		correctCombined =[];
		whichCombined = [];
		con=[];
		posdiff =[];
		diff = [];
		s = cell(1,numSession);
		for iSession = 1:numSession
			s{iSession} = stimStr(iSession).stimulus.stair.(stairNames{taskType});
			Task{iSession} = stimStr(iSession).task;
			nTrial = s{iSession}.trialNum;
			respCombined = [respCombined, Task{iSession}.randVars.resp(1:nTrial)];
			correctCombined = [correctCombined, Task{iSession}.randVars.correct(1:nTrial)];
			whichCombined = [whichCombined, Task{iSession}.randVars.whichint(1:nTrial)];
			con = [con, Task{iSession}.randVars.con(1:nTrial)];
			posdiff = [posdiff, Task{iSession}.randVars.diff(1:nTrial)];
		end
		sCombined = doStaircase('combine', s);
		nTrialCombined = sCombined.trialNum;
		switch taskType
		case 1
			testVal = unique(sCombined.testValues);
			nTotal = arrayfun(@(x) sum(con==x), testVal);
			whichTrial = arrayfun(@(x) find(con==x), testVal, 'UniformOutput', false);
			nC = cellfun(@(x) sum(correctCombined(x)), whichTrial);
			xs = 0:.001:0.6;
		case {2,3}
			testVal = sCombined.testValues;
			testVal(whichCombined==2) = -testVal(whichCombined==2);
			testVal = unique(testVal);
			posdiff(whichCombined==2) = -posdiff(whichCombined==2);
			nTotal = arrayfun(@(x) sum(posdiff==x), testVal);
			whichTrial = arrayfun(@(x) find(posdiff==x), testVal, 'UniformOutput', false);
			nC= cellfun(@(x) sum(respCombined(x)==1), whichTrial);
			xs = -6:.01:6;
		end
		pC = nC./nTotal;
		fit{i} = fitCG(testVal, nTotal, nC);
		mu = fit{i}.fitparams(1);
		sigma = fit{i}.fitparams(2);
		lambda = fit{i}.fitparams(3);
		cgfunc = fit{i}.cgfunc;
		fitVal = arrayfun(@(x) cgfunc(mu,sigma,lambda,x), xs);
		pfitVal = arrayfun(@(x) cgfunc(mu,sigma,lambda,x), testVal);
		error = sqrt(pfitVal.*(1-pfitVal)./nTotal);
		if i < 7
			subplot(nrow/2,ncol,i)
		else
			subplot(nrow/2,ncol,i-6)
		end
		myerrorbar(testVal, pC, 'Symbol','o');
		hold on;
    	plot(xs, fitVal, 'k')
    	if taskType==1
    		title(sprintf('Stair Contrast: Ecc=%d   (N=%d) \n \\mu=%0.2f \\sigma=%0.2f \\lambda=%0.2f \n', ...
    			stimStr(1).eccentricity, nTrialCombined, mu,sigma,lambda), 'Interpreter','tex');
    		ylabel('Percent Correct');
    		xlabel('Contrast');
    		axis([0 0.6 0 1]); box off; drawnow;
    	else
    		title(sprintf('Ecc=%d   Con=%0.3f   (N=%d) \n \\mu=%0.2f \\sigma=%0.2f \\lambda=%0.2f \n', ...
    			stimStr(1).eccentricity, stimStr(1).stimulus.contrast, nTrialCombined, mu,sigma,lambda), 'Interpreter','tex');
    		ylabel('Percent Choice Interval 1');
    		xlabel('Postition difference of targets between interval 1 and 2 (deg)');
    		axis([-6 6 0 1]); box off; drawnow;
    	end

	end



