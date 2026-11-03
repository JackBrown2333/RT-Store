function X=r30full_read(outDir,indices,tooth)
% Read exact raw intervals; absent sparse rows mean zero removal/load.
D=load(fullfile(outDir,'full_run_plan.mat'),'meta','runToken');m=D.meta;token=D.runToken;
indices=double(indices(:));assert(all(indices>=1&indices<=m.intervalCount&mod(indices,1)==0));
if nargin<3,folder=fullfile(outDir,'global_chunks');X=zeros(numel(indices),16);
else,assert(ismember(tooth,m.options.toothIds));folder=fullfile(outDir,sprintf('tooth_%02d',tooth),'chunks');X=zeros(numel(indices),13,2);end
groups=floor((indices-1)/m.options.chunkIntervals)+1;
for j=unique(groups)'
 D=load(fullfile(folder,sprintf('chunk_%06d.mat',j)),'chunk');c=D.chunk;
 assert(strcmp(c.runToken,token)&&c.chunkIndex==j);
 if nargin>=3,assert(c.tooth==tooth);end
 q=find(groups==j);local=indices(q)-c.firstIndex+1;
 [found,where]=ismember(local,double(c.rowIndex));
 if nargin<3,X(q(found),:)=c.values(where(found),:);
 else,X(q(found),:,:)=c.values(where(found),:,:);end
end
end
