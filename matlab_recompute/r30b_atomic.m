function r30b_atomic(file,name,value)
% All checkpoints/results are new-engine files; install only complete writes.
folder=fileparts(file);if isempty(folder),folder=pwd;end
temp=[tempname(folder) '.mat'];D=struct();D.(name)=value;
save(temp,'-struct','D','-v7');[ok,msg]=movefile(temp,file,'f');assert(ok,'R30B:Save','%s',msg);
end
