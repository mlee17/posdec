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
		if stimulus.fixOrigin(1) ~= 0
			stimlist(i).eccentricity = stimulus.eccentricity-stimulus.fixOrigin(1);
		else 
			stimlist(i).eccentricity = stimulus.eccentricity;
		end
		stimlist(i).contrast = stimulus.contrast;
        stimlist(i).ISI = stimulus.ISI;
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
			lab = sprintf('%s_Ecc%02d_ISI%d', taskName{1}, stimlist(i).eccentricity, stimlist(i).ISI*1000);
			if isfield(stim, lab)
					stim(size(stim,1)+1).(lab) = stimlist(i);
			else 
				stim(1).(lab) = stimlist(i);
			end
		case {2,3}
			lab = sprintf('Ecc%d_Con%d_ISI%d', stimlist(i).eccentricity, stimlist(i).contrast*1000, stimlist(i).ISI*1000);
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
    fit = cell(length(fieldNames),1);
	figure(2);figure(1);
	for i = 1:length(fieldNames)
        clearvars nC testVal nTotal whichTrial
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
			nTrial = length(Task{iSession}.randVars.whichint(~isnan(Task{iSession}.randVars.whichint)));
			respCombined = [respCombined, Task{iSession}.randVars.resp(1:nTrial)];
			correctCombined = [correctCombined, Task{iSession}.randVars.correct(1:nTrial)];
			whichCombined = [whichCombined, Task{iSession}.randVars.whichint(1:nTrial)];
			con = [con, Task{iSession}.randVars.con(1:nTrial)];
			posdiff = [posdiff, Task{iSession}.randVars.diff(1:nTrial)];
        end
        
		sCombined = doStaircase('combine', s);
        h = doStaircase('getHistory', sCombined);
       
		nTrialCombined = sCombined.trialNum;
        if nTrialCombined < nTrial
            respCombined = respCombined(1:nTrialCombined);
            correctCombined = correctCombined(1:nTrialCombined);
            whichCombined = whichCombined(1:nTrialCombined);
            con = con(1:nTrialCombined);
            posdiff = posdiff(1:nTrialCombined);
            nTrial = nTrialCombined;
        elseif nTrialCombined > nTrial
            nTrialCombined = nTrial;
        end
        nanInd = isnan(respCombined) | isnan(correctCombined);
        whichCombined(nanInd) = [];
        respCombined(isnan(respCombined)) = [];
        correctCombined(isnan(correctCombined)) = [];
        testValCombined = sCombined.testValues;
        testValCombined(nanInd(1:nTrialCombined)) = [];
        con(nanInd) = [];
        posdiff(nanInd) = [];
		switch taskType
		case 1
			testVal = unique(testValCombined);
			nTotal = arrayfun(@(x) sum(con==x), testVal);
            testVal(nTotal==0) = [];
            nTotal(nTotal==0) = [];
			whichTrial = arrayfun(@(x) find(con==x), testVal, 'UniformOutput', false);
            whichTrial(cellfun(@(x) isempty(x), whichTrial)) = [];
			nC = cellfun(@(x) sum(correctCombined(x)), whichTrial);
%             discard = (nTotal==1);
%             testVal(discard) = []; nTotal(discard)=[]; nC(discard)=[];
			xs = 0:.001:0.6;
            fit{i} = fitCGCon(testVal, nTotal, nC);
            		pC = nC./nTotal;
		mu = fit{i}.fitparams(1);
		sigma = fit{i}.fitparams(2);
		lambda = fit{i}.fitparams(3);
%         gamma = fit{i}.fitparams(4);
		cgfunc = fit{i}.cgfunc;
		fitVal = arrayfun(@(x) cgfunc(mu,sigma,lambda,x), xs);
		pfitVal = arrayfun(@(x) cgfunc(mu,sigma,lambda,x), testVal);
		case {2,3}
			testVal = testValCombined;
			testVal(whichCombined==2) = -testVal(whichCombined==2);
			testVal = unique(testVal);
			posdiff(whichCombined==2) = -posdiff(whichCombined==2);
			nTotal = arrayfun(@(x) sum(posdiff==x), testVal);
            testVal(nTotal==0) = [];
            nTotal(nTotal==0) = [];
			whichTrial = arrayfun(@(x) find(posdiff==x), testVal, 'UniformOutput', false);
            whichTrial(cellfun(@(x) isempty(x), whichTrial)) = [];
			nC= cellfun(@(x) sum(respCombined(x)==1), whichTrial);
%             nC= cellfun(@(x) sum(correctCombined(x)), whichTrial);
%             discard = (nTotal==1);
%             testVal(discard) = []; nTotal(discard)=[]; nC(discard)=[];
			xs = -6:.01:6;
            fit{i} = fitCG(testVal, nTotal, nC);
            pC = nC./nTotal;
		mu = fit{i}.fitparams(1);
		sigma = fit{i}.fitparams(2);
		lambda = fit{i}.fitparams(3);
		cgfunc = fit{i}.cgfunc;
		fitVal = arrayfun(@(x) cgfunc(mu,sigma,lambda,x), xs);
		pfitVal = arrayfun(@(x) cgfunc(mu,sigma,lambda,x), testVal);
		end

		error = sqrt(pfitVal.*(1-pfitVal)./nTotal);
        fit{i}.xs = xs;
        fit{i}.fitVal = fitVal;
%         C = my_norminv(1-0.1,0,1) - my_norminv(0.1,0,1);
		if i < 7
			subplot(nrow/2,ncol,i)
		else
			subplot(nrow/2,ncol,i-6)
		end
		myerrorbar(testVal, pC, 'Symbol','o', 'MarkerSize',6);
		hold on;
    	plot(xs, fitVal, 'k')
    	if taskType==1
    		title(sprintf('Stair Contrast: Ecc=%d   ISI=%d  (N=%d) \n \\mu=%0.2f \\sigma=%0.2f \\lambda=%0.2f \n', ...
    			stimStr(1).eccentricity, stimStr(1).ISI*1000, nTrialCombined, mu,sigma,lambda), 'Interpreter','tex');
    		ylabel('Percent Correct');
    		xlabel('Contrast');
    		axis([0 0.6 0 1]); box off;
            p75 = .75*(1-lambda) + lambda/2;
            [c, ind75] = min(abs(fitVal-p75));
            line([0 xs(ind75)], [p75 p75], 'Color','b', 'LineStyle', ':');
            if mu+sigma<=0.6 && mu+sigma >=0
            line([xs(ind75) xs(ind75)], [p75 0], 'Color','b', 'LineStyle', ':');
            text(xs(ind75)+0.05, 0.05, sprintf('%.2f', xs(ind75)));
            end
            line([0 mu], [.5 .5], 'Color','r', 'LineStyle', ':');
            if mu <=0.6 && mu>=0
            line([mu mu], [.5 0], 'Color','r', 'LineStyle', ':');
            
            end
            
            drawnow;
    	else
    		title(sprintf('Ecc=%d   Con=%0.3f   (N=%d) \n \\mu=%0.2f \\sigma=%0.2f \\lambda=%0.2f \n', ...
    			stimStr(1).eccentricity, stimStr(1).stimulus.contrast, nTrialCombined, mu,sigma,lambda), 'Interpreter','tex');
    		ylabel('Percent Choice Interval 1');
    		xlabel('Postition difference of targets between interval 1 and 2 (deg)');
    		axis([-6 6 0 1]); box off; 
            p75 = .75*(1-lambda) + lambda/2;
            p25 = .25*(1-lambda) + lambda/2;
            [c, ind75] = min(abs(fitVal-p75));
            [c, ind25] = min(abs(fitVal-p25));
            line([-6 xs(ind75)], [p75 p75], 'Color','b', 'LineStyle', ':');
            line([-6 xs(ind25)], [p25 p25], 'Color','b', 'LineStyle', ':');
            if mu+sigma<=6 && mu+sigma>=-6
            line([xs(ind75) xs(ind75)], [p75 0], 'Color','b', 'LineStyle', ':');
            text(xs(ind75)+0.5, 0.05, sprintf('%.2f', xs(ind75)));
            line([xs(ind25) xs(ind25)], [p25 0], 'Color','b', 'LineStyle', ':');
            text(xs(ind25)-2, 0.05, sprintf('%.2f', xs(ind25)));
            end
            line([-6 mu], [.5 .5], 'Color','r', 'LineStyle', ':');
            if mu<=6 && mu>=-6
            line([mu,mu], [.5 0], 'Color','r', 'LineStyle', ':');
            
            end
            
            drawnow;
            
        end
 
    end

