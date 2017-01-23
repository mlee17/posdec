function testlocStimCal(session)
% creating stimvol files

session = 's035020161212';
datadir = '~/data/testloc/';
stimdir = sprintf(['%s', session, '/Etc/'],datadir);
stimfiles = dir(sprintf('%s1*',stimdir));

for i = 1:length(stimfiles)

    display(stimfiles(i).name)
    
    s=load(sprintf('%s%s', stimdir, stimfiles(i).name));
    e = getTaskParameters(s.myscreen, s.task{2});

    % discarding last trial if 
    if any(isnan(e.trials(end).volnum))
        clear tmp
        tmp=find(s.myscreen.events.volnum==e.trials(end).volnum(1)-1);
        s.myscreen.events.n=tmp(end);
    end
    e = getTaskParameters(s.myscreen, s.task{2});
    
 loccon = zeros(1,e.nTrials)-1;
 
loccon(find(e.parameter.location==-14 & e.parameter.contrast==.25)) = 11;
loccon(find(e.parameter.location==-14 & e.parameter.contrast==.5)) = 12;
loccon(find(e.parameter.location==-14 & e.parameter.contrast==1)) = 13;
loccon(find(e.parameter.location==-8 & e.parameter.contrast==.25)) = 21;
loccon(find(e.parameter.location==-8 & e.parameter.contrast==.5)) = 22;
loccon(find(e.parameter.location==-8 & e.parameter.contrast==1)) = 23;
loccon(find(e.parameter.location==8 & e.parameter.contrast==.25)) = 31;
loccon(find(e.parameter.location==8 & e.parameter.contrast==.5)) = 32;
loccon(find(e.parameter.location==8 & e.parameter.contrast==1)) = 33;
loccon(find(e.parameter.location==14 & e.parameter.contrast==.25)) = 41;
loccon(find(e.parameter.location==14 & e.parameter.contrast==.5)) = 42;
loccon(find(e.parameter.location==14 & e.parameter.contrast==1)) = 43;

 eval(sprintf('save %s -struct s',sprintf('%s%s', stimdir, stimfiles(i).name)));
    addCalculatedVar('loccon',loccon, sprintf('%s%s', stimdir, stimfiles(i).name),'taskNum=2','phaseNum=1')
    
    
end