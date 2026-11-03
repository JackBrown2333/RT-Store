function R=run_r30_full_force(sourceMat,coreDir,outDir,planFile,user)
% All selected teeth, from nominal stock, on one common full-stroke time grid.
% Two independent flank states in memory. Sparse load chunks on disk.
% Core geometry/force law are called unchanged; no periodic load copying.
if nargin<5,user=struct();end
addpath(coreDir);D=load(sourceMat,'M');M=D.M;clear D;
D=load(planFile,'plan');plan=D.plan;clear D;
assert(isequaln(plan.parameters,M.parameters)&&abs(plan.modelTimeRange_s(1))<1e-12,...
 'R30FULL:Plan','需要同一工况、从模型t=0开始的预检计划。');
o=struct('toothIds',1:30,'ns',M.options.ns,'nz',M.options.nz,'order',M.options.order,...
 'dt_s',M.forceDt_s,'substeps',M.options.substeps,'forceStripWidth_mm',M.options.forceStripWidth_mm/2,...
 'endTime_s',plan.modelTimeRange_s(2),'chunkIntervals',2048,'checkpointEvery_s',30,...
 'progressEvery_s',5,'maxChunksPerCall',Inf,'requireClearStart',true,'plot',true,...
 'targetTolerance_mm',.001,'edgeBand_mm',1,'overviewBin_s',.01,'loadThreshold_N',.01);
names=fieldnames(user);for j=1:numel(names),assert(isfield(o,names{j}));o.(names{j})=user.(names{j});end
o.toothIds=double(o.toothIds(:)');
assert(all(ismember(o.toothIds,1:30))&&numel(unique(o.toothIds))==numel(o.toothIds)&&~isempty(o.toothIds));
assert(o.endTime_s>0&&o.chunkIntervals>=1&&mod(o.chunkIntervals,1)==0&&o.maxChunksPerCall>0);
assert(o.targetTolerance_mm>0&&o.edgeBand_mm>0&&o.overviewBin_s>=o.dt_s);
P=M.parameters;coefficients=M.coefficients;clear M;
assert(P.p.gear.teeth==30&&P.p.wheel.starts==1&&abs(P.p.gear.helix_deg)<1e-12&&P.installation.feed>0);
u=struct('ns',o.ns,'nz',o.nz,'order',o.order,'dt_s',o.dt_s,'substeps',o.substeps,...
 'forceStripWidth_mm',o.forceStripWidth_mm,'coefficients',coefficients,'saveFrames',false,'plot',false);
co=r30b_options(u);nt=ceil(o.endTime_s/o.dt_s);nc=ceil(nt/o.chunkIntervals);
if ~isfolder(outDir),mkdir(outDir);end
physicsOptions=rmfield(o,{'maxChunksPerCall','plot','progressEvery_s','checkpointEvery_s'});
meta=struct('kind','full_force_stream_v1','parameters',P,'options',physicsOptions,'coreOptions',co,...
 'intervalCount',nt,'chunkCount',nc,'preflightModelTimeRange_s',plan.modelTimeRange_s,'coreSignature',core_signature(coreDir),...
 'engineSource',fileread([mfilename('fullpath') '.m']),...
 'collectorSource',fileread(fullfile(fileparts(mfilename('fullpath')),'r30full_collect.m')));
meta.readerSource=fileread(fullfile(fileparts(mfilename('fullpath')),'r30full_read.m'));
mf=fullfile(outDir,'full_run_plan.mat');
if isfile(mf)
 D=load(mf,'meta','runToken');assert(isequaln(meta,D.meta),'R30FULL:Resume','代码/参数/离散设置有变化，请使用新的输出目录。');runToken=D.runToken;
else
 runToken=datestr(now,'yyyymmddTHHMMSSFFF');save(mf,'meta','runToken','-v7');
end
R=struct('completed',false,'fullGearComputed',false,'outputDir',outDir,'intervalCount',nt);
processed=0;clock=tic;progress=tic;checkpoint=tic;
fprintf('全程试算：%d齿、%d侧；t=0..%.9f s；%d区间；条带%.4g mm。\n',numel(o.toothIds),2*numel(o.toothIds),nt*o.dt_s,nt,o.forceStripWidth_mm);
for g=1:numel(o.toothIds)
 tooth=o.toothIds(g);leaf=fullfile(outDir,sprintf('tooth_%02d',tooth));chunks=fullfile(leaf,'chunks');
 if ~isfolder(chunks),mkdir(chunks);end
 cp=fullfile(leaf,'checkpoint.mat');
 if isfile(cp)
  D=load(cp,'S');S=D.S;clear D;assert(strcmp(S.runToken,runToken)&&S.tooth==tooth);
  if strcmp(S.phase,'done'),fprintf('复用已完成第%d齿。\n',tooth);continue;end
 else
  S=struct('runToken',runToken,'tooth',tooth,'phase','running','nextChunk',1,'Q',{{}},...
   'removedVolume_mm3',[0 0],'maxPartitionError_mm3',0,'maxSubstepError_mm3',0,...
   'initialStockHitCount',[0 0],'initialMinLevel',[Inf Inf],'geometryEvaluations',[0 0]);
  for f=1:2
   Q=r30b_init(P,co,[tooth f]);[ss,zz]=samples(Q);
   level=r30b_level(ss,zz,ones(size(ss))*Q.stock,0,Q.p,Q.A,tooth,f);
   S.initialStockHitCount(f)=nnz(level<-1e-10);S.initialMinLevel(f)=min(level(:));S.Q{f}=Q;
  end
  if o.requireClearStart,assert(~any(S.initialStockHitCount),'R30FULL:Start','第%d齿初始材料已经与砂轮相交，需重新检查行程起点。',tooth);end
  r30b_atomic(cp,'S',S);
 end
 fprintf('第%d齿：从块%d/%d继续；两侧材料状态分别继承。\n',tooth,S.nextChunk,nc);
 for j=S.nextChunk:nc
  first=(j-1)*o.chunkIntervals+1;last=min(nt,j*o.chunkIntervals);n=last-first+1;
  X=zeros(n,13,2);possible=false(1,2);
  for f=1:2,possible(f)=r30b_gate(P,[tooth f],(first-1)*o.dt_s,last*o.dt_s,co.normalDepthLimit_mm);end
  for f=1:2
   Q=S.Q{f};
   if ~possible(f),Q.inPass=false;S.Q{f}=Q;continue;end
   for k=first:last
    ta=(k-1)*o.dt_s;t=k*o.dt_s;
    if ~r30b_gate(P,[tooth f],ta,t,co.normalDepthLimit_mm),Q.inPass=false;continue;end
    before=r30b_heights(Q);subVolume=0;
    for sub=1:o.substeps
     ts=ta+(t-ta)*sub/o.substeps;
     if r30b_gate(P,[tooth f],ts,ts,co.normalDepthLimit_mm)
      if ~Q.inPass,Q.passStart=Q.h;Q.passId=Q.passId+1;Q.inPass=true;end
      [Q,v]=r30b_step(Q,ts,co);subVolume=subVolume+v;S.geometryEvaluations(f)=S.geometryEvaluations(f)+1;
     else,Q.inPass=false;end
    end
    after=r30b_heights(Q);dh=before-after;volume=sum(Q.weightsFine(:).*dh(:));area=sum(Q.weightsFine(dh>0));
    temp=Q;temp.h=Q.passStart;pass=r30b_heights(temp);
    try
     rows=r30b_force(r30b_strips(Q,before,after,pass,t,o.dt_s,co,false),co.coefficients);
    catch stripError
     failure=struct('tooth',tooth,'side',f,'intervalIndex',k,'t_s',t,'dt_s',o.dt_s,...
      'identifier',stripError.identifier,'message',stripError.message,'Q',Q,'before',before,'after',after,'passStart',pass,'coreOptions',co);
     r30b_atomic(fullfile(outDir,'last_force_error.mat'),'failure',failure);
     rethrow(stripError);
    end
    assigned=0;if ~isempty(rows),assigned=sum([rows.dV_mm3]);end
    pe=abs(assigned-volume);se=abs(subVolume-volume);assert(pe<1e-9&&se<1e-9,'R30FULL:Volume','条带/子步/材料体积不闭合。');
    S.maxPartitionError_mm3=max(S.maxPartitionError_mm3,pe);S.maxSubstepError_mm3=max(S.maxSubstepError_mm3,se);
    a=r30b_sum(rows,co.referenceOutward,t,o.dt_s,volume,area);
    X(k-first+1,:,f)=[volume area a.force_N a.normalVector_N a.moment_Nmm a.Fn_magnitude_sum_N a.nForceCells];
    if toc(progress)>o.progressEvery_s
     fprintf('  第%d齿第%d侧，t=%.6f/%.6f s；本次耗时%.1f s。\n',tooth,f,t,nt*o.dt_s,toc(clock));progress=tic;
    end
   end
   S.Q{f}=Q;
  end
  keep=find(any(any(X~=0,3),2));
  chunk=struct('runToken',runToken,'tooth',tooth,'chunkIndex',j,'firstIndex',first,'lastIndex',last,...
   'rowIndex',uint32(keep),'values',X(keep,:,:));
  r30b_atomic(fullfile(chunks,sprintf('chunk_%06d.mat',j)),'chunk',chunk);
  S.removedVolume_mm3=S.removedVolume_mm3+[sum(X(:,1,1)) sum(X(:,1,2))];S.nextChunk=j+1;processed=processed+1;
  if toc(checkpoint)>o.checkpointEvery_s,r30b_atomic(cp,'S',S);checkpoint=tic;end
  if processed>=o.maxChunksPerCall,r30b_atomic(cp,'S',S);R.status='PAUSED_AT_SAVED_CHUNK';return;end
 end
 finalState=struct('tooth',tooth,'t_s',nt*o.dt_s,'faces',{{}});faceStatistics=zeros(2,15);
 for f=1:2
  Q=S.Q{f};H=r30b_heights(Q);weights=Q.weightsFine;stateVolume=sum(weights(:).*(Q.stock-H(:)));
  closure=abs(stateVolume-S.removedVolume_mm3(f));assert(closure<max(1e-8,1e-10*abs(stateVolume)),'R30FULL:Closure','全程累计体积与最终高度不一致。');
  [ss,zz]=samples(Q);endLevel=r30b_level(ss,zz,ones(size(H))*Q.stock,nt*o.dt_s,Q.p,Q.A,tooth,f);
  tol=o.targetTolerance_mm;neg=zz<=-P.p.gear.faceWidth_mm/2+o.edgeBand_mm;pos=zz>=P.p.gear.faceWidth_mm/2-o.edgeBand_mm;
  faceStatistics(f,:)=[tooth f stateVolume closure max(H(:))*1000 max(0,-min(H(:)))*1000,...
   fraction(H<Q.stock-co.heightTolerance_mm,weights),fraction(abs(H)<=tol,weights),fraction(H>tol,weights),fraction(H<-tol,weights),...
   max_or_nan(H(neg))*1000 max_or_nan(H(pos))*1000 S.initialStockHitCount(f) nnz(endLevel<-1e-10) S.geometryEvaluations(f)];
  finalState.faces{f}=struct('height_mm',H,'support_s_mm',Q.support_s,'support_z_mm',Q.support_z,...
   'sample_s_mm',ss(:,1),'sample_z_mm',zz(1,:),'passId',Q.passId,'passStartHeight_mm',pass_height(Q));
 end
 r30b_atomic(fullfile(leaf,'final_state.mat'),'finalState',finalState);
 S.phase='done';diagnostics=rmfield(S,{'Q'});diagnostics.faceStatistics=faceStatistics;
 r30b_atomic(fullfile(leaf,'diagnostics.mat'),'diagnostics',diagnostics);
 r30b_atomic(cp,'S',S);fprintf('第%d齿完成；移除体积 %.6g / %.6g mm^3。\n',tooth,S.removedVolume_mm3);
end
R=r30full_collect(outDir);if o.plot,r30full_figures(outDir);end
end
function [ss,zz]=samples(Q)
[ss,zz]=ndgrid(reshape(Q.s',[],1),reshape(Q.z',[],1));
end
function H=pass_height(Q),Q.h=Q.passStart;H=r30b_heights(Q);end
function x=fraction(mask,w),x=sum(w(mask))/sum(w(:));end
function x=max_or_nan(v),if isempty(v),x=NaN;else,x=max(v);end;end
function s=core_signature(core)
files=dir(fullfile(core,'*.m'));s='';for j=1:numel(files),s=[s files(j).name char(10) fileread(fullfile(core,files(j).name))];end
end
