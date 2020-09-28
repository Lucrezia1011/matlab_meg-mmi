clc
ff = dir;


for ii = 1:length(ff)
    if strncmp(ff(ii).name,'s2',2)
        sdan = ff(ii).name(2:end);
        if ~exist(['sub-',sdan],'dir')
            mkdir(['sub-',sdan])
        end
        eval(sprintf('!mv s%s/*/* sub-%s/',sdan,sdan ))
    end
end