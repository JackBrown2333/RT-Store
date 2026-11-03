function C=run_r30_resolution_check(resultDir,coreDir,outDir,user)
% Short-window comparisons seeded from exact native replay states.
% Uses the current Fn(b,heq,lg) law. Full-pass ae_pass is NOT evaluated here.
% A smaller force interval is compared both raw and averaged to old intervals.
if nargin<4,user=struct();end
o=struct('windowDuration_s',.002,'targetIndices',[],'plot',true,'progressEvery_s',5);
names=fieldnames(user);for j=1:numel(names),assert(isfield(o,names{j}));o.(names{j})=user.(names{j});end
addpath(coreDir);D=load(fullfile(resultDir,'multitooth_pilot.mat'),'M');M=D.M;clear D;
dt=M.forceDt_s;nt=numel(M.time_s);nbase=min(nt,max(1,round(o.windowDuration_s/dt)));
assert(M.options.order==2&&M.coefficients.background_pressure==0);
if isfield(M,'preflight'),assert(M.preflight.candidateCoverageComplete);end
if ~isfolder(outDir),mkdir(outDir);end
jobDir=fullfile(outDir,'jobs');anchorDir=fullfile(outDir,'anchors');
if ~isfolder(jobDir),mkdir(jobDir);end;if ~isfolder(anchorDir),mkdir(anchorDir);end
signature=core_signature(coreDir);leaves=cell(numel(M.toothIds),1);leafOptions=leaves;common=0:nt;
for g=1:numel(M.toothIds)
 if ~M.calculatedHistoryForTooth(g),continue;end
 leaf=fullfile(resultDir,'animation_export','replay',sprintf('tooth_%02d',M.toothIds(g)));
 assert(isfile(fullfile(leaf,'r30b_result.mat')),'R30R:Replay','缺少本机原生重放：%s',leaf);
 D=load(fullfile(leaf,'r30b_result.mat'),'R30B');R=D.R30B;clear D;
 assert(isequaln(R.parameters,M.parameters)&&isequaln(R.forceCoefficients,M.coefficients));
 fields={'ns','nz','order','substeps','forceStripWidth_mm','dt_s'};
 for q=1:numel(fields),assert(isequal(R.options.(fields{q}),M.options.(fields{q})));end
 assert(isequal(R.faces,[double(M.toothIds(g)) 1;double(M.toothIds(g)) 2]));
 assert(numel(R.poses)==nt&&numel(R.frameFiles)==numel(R.frameIndices));
 delta=vertcat(R.poses.force_N)-reshape(M.toothForce_N(:,g,:),nt,3);assert(max0(abs(delta))<1e-8);
 relative=strrep(M.perToothResultFiles{g},'\',filesep);
 D=load(fullfile(fileparts(fullfile(resultDir,relative)),'r30b_checkpoint.mat'),'S');
 assert(strcmp(D.S.phase,'done')&&strcmp(D.S.meta.codeSignature,signature),'R30R:Code','核心代码与原计算不一致。');clear D;
 leaves{g}=struct('directory',leaf,'R',R);leafOptions{g}=R.options;
 common=intersect(common,[0 double(R.frameIndices(:)')]);
end
assert(~isempty(common));
if isempty(o.targetIndices)
 [~,k1]=max(M.F_resultant_N);[~,k2]=max(M.moment_Nmm(:,3));[~,k3]=min(M.moment_Nmm(:,3));targets=unique([k1 k2 k3]);
else,targets=unique(o.targetIndices(:)');end
assert(all(targets>=1&targets<=nt&mod(targets,1)==0));
windows=zeros(numel(targets),3);
for w=1:numel(targets)
 ideal=max(0,targets(w)-floor(nbase/2));anchor=max(common(common<=ideal));last=min(nt,anchor+nbase);
 assert(targets(w)<=last,'R30R:SparseReplay','代表点前的原生帧过稀，无法构造小窗口。');
 windows(w,:)=[targets(w) anchor last];
end
[~,keep]=unique(windows(:,2:3),'rows','stable');windows=windows(sort(keep),:);
ownSource=fileread([mfilename('fullpath') '.m']);
meta=struct('kind','r30_short_resolution_v1','parameters',M.parameters,'coefficients',M.coefficients,...
 'grid',[M.options.ns M.options.nz M.options.order],'time_s',M.time_s,'force_N',M.force_N,...
 'moment_Nmm',M.moment_Nmm,'windows',windows,'dt_s',dt,'stripWidth_mm',M.options.forceStripWidth_mm,...
 'coreSignature',signature,'leafOptions',{leafOptions},'source',ownSource);
mf=fullfile(outDir,'resolution_plan.mat');
if isfile(mf),D=load(mf,'meta');assert(isequaln(D.meta,meta),'R30R:Resume','输入、代码或窗口已改变，请使用新输出目录。');
else,save(mf,'meta','-v7');end
globalCases=cell(size(windows,1),3);maximumBaselineForce=0;maximumBaselineMoment=0;
heightDifference=zeros(size(windows,1),1);completedJobs=0;
for g=1:numel(leaves)
 if isempty(leaves{g}),continue;end
 R=leaves{g}.R;leaf=leaves{g}.directory;
 D=load(fullfile(leaf,R.initialStateFile),'initialState');initial=D.initialState;clear D;
 assert(abs(initial.t_s-M.window_s(1))<1e-10);
 H=cell(2,1);for f=1:2,H{f}=initial.faces{f}.height_mm;end
 snapshots=cell(size(windows,1),1);
 for w=1:size(windows,1),if windows(w,2)==0,snapshots{w}=H;end;end
 lastTime=initial.t_s;
 for j=1:numel(R.frameIndices)
  k=double(R.frameIndices(j));if k>max(windows(:,2)),break;end
  D=load(fullfile(leaf,strrep(R.frameFiles{j},'\',filesep)),'frame');frame=D.frame;clear D;
  assert(frame.index==k&&abs(frame.t_s-M.time_s(k))<1e-10);
  for f=1:2
   u=frame.faces{f}.stateUpdates;assert(abs(u.since_t_s-lastTime)<1e-10);
   ids=double(u.index(:));hh=H{f};v=u.height_mm(:);
   assert(numel(ids)==numel(v)&&numel(unique(ids))==numel(ids)&&all(v<=hh(ids)+1e-13));
   hh(ids)=v;H{f}=hh;
  end
  for w=1:size(windows,1),if windows(w,2)==k,snapshots{w}=H;end;end
  lastTime=frame.t_s;
 end
 for w=1:size(windows,1)
  assert(~isempty(snapshots{w}));anchor=struct('tooth',double(M.toothIds(g)),...
   'sourceIndex',windows(w,2),'height',{snapshots{w}});
  af=fullfile(anchorDir,sprintf('tooth_%02d_w%02d.mat',M.toothIds(g),w));
  if isfile(af),D=load(af,'anchor');assert(isequaln(D.anchor,anchor),'R30R:State','源材料状态改变，请使用新目录。');
  else,save(af,'anchor','-v7');end
  jf=fullfile(jobDir,sprintf('tooth_%02d_w%02d.mat',M.toothIds(g),w));
  if isfile(jf)
   D=load(jf,'job');job=D.job;clear D;fprintf('复用：第%d齿，第%d窗口。\n',M.toothIds(g),w);
  else
   fprintf('收敛对照：第%d齿，第%d/%d窗口；基准+条带减半、时间步减半。\n',M.toothIds(g),w,size(windows,1));
   job=struct('values',{{[],[],[]}},'baselineForceError_N',0,'baselineMomentError_Nmm',0,...
    'maximumFinalHeightDifference_mm',0,'elapsed_s',0);clock=tic;
   for f=1:2
    Q=r30b_init(R.parameters,R.options,R.faces(f,:));h0=anchor.height{f};
    assert(isequal(size(h0),size(Q.weightsFine)));
    base=local_case(Q,h0,R,M,windows(w,2:3),1,[M.options.forceStripWidth_mm M.options.forceStripWidth_mm/2],o.progressEvery_s);
    fine=local_case(Q,h0,R,M,windows(w,2:3),2,M.options.forceStripWidth_mm,o.progressEvery_s);
    expected=R.facePoses{f}(windows(w,2)+1:windows(w,3));
    eF=max0(abs(base.values(:,3:5,1)-vertcat(expected.force_N)));
    eN=max0(abs(base.values(:,6:8,1)-vertcat(expected.normalVector_N)));
    eM=max0(abs(base.values(:,9:11,1)-vertcat(expected.moment_Nmm)));
    eV=max0(abs(base.values(:,2,1)-[expected.dV_mm3]'));
    assert(max(eF,eN)<1e-7&&eM<1e-5&&eV<1e-9,'R30R:Baseline','基准不能重现原结果，已停止对照。');
    job.baselineForceError_N=max(job.baselineForceError_N,max(eF,eN));
    job.baselineMomentError_Nmm=max(job.baselineMomentError_Nmm,eM);
    job.maximumFinalHeightDifference_mm=max(job.maximumFinalHeightDifference_mm,max0(abs(fine.finalHeight-base.finalHeight)));
    job.values{1}=add_values(job.values{1},base.values(:,:,1));
    job.values{2}=add_values(job.values{2},fine.values(:,:,1));
    job.values{3}=add_values(job.values{3},base.values(:,:,2));
   end
   job.elapsed_s=toc(clock);r30b_atomic(jf,'job',job);
  end
  for c=1:3,globalCases{w,c}=add_values(globalCases{w,c},job.values{c});end
  maximumBaselineForce=max(maximumBaselineForce,job.baselineForceError_N);
  maximumBaselineMoment=max(maximumBaselineMoment,job.baselineMomentError_Nmm);
  heightDifference(w)=max(heightDifference(w),job.maximumFinalHeightDifference_mm);completedJobs=completedJobs+1;
 end
end
comparisons=cell(size(windows,1),1);
headers={'model_time_s','dV_mm3','Fx_N','Fy_N','Fz_N','Nx_N','Ny_N','Nz_N','Mx_Nmm','My_Nmm','Mz_Nmm','active_area_mm2'};
for w=1:size(windows,1)
 B=globalCases{w,1};F=globalCases{w,2};S=globalCases{w,3};assert(size(F,1)==2*size(B,1));
 ids=windows(w,2)+1:windows(w,3);
 assert(max0(abs(B(:,3:5)-M.force_N(ids,:)))<1e-7&&max0(abs(B(:,9:11)-M.moment_Nmm(ids,:)))<1e-5,...
  'R30R:GlobalBaseline','各齿基准叠加不能重现原合力/力矩。');
 FA=(F(1:2:end,:)+F(2:2:end,:))/2;FA(:,1)=B(:,1);FA(:,2)=F(1:2:end,2)+F(2:2:end,2);
 r=struct('sourceTarget',windows(w,1),'initialSourceFrame',windows(w,2),...
  'lastSourceInterval',windows(w,3),'modelTime_s',[M.window_s(1)+windows(w,2)*dt B(end,1)],...
  'timeHalf',compare(B,F,FA,dt,dt/2,M.parameters.installation),...
  'stripHalf',compare(B,S,S,dt,dt,M.parameters.installation),...
  'maximumTimeHalfFinalHeightDifference_mm',heightDifference(w));
 comparisons{w}=r;
 for c=1:3
  labels={'base','time_half','strip_half'};
  write_csv(fullfile(outDir,sprintf('w%02d_%s.csv',w,labels{c})),headers,globalCases{w,c});
 end
 if o.plot,plot_comparison(B,F,S,w,outDir,M.window_s(1)+windows(w,2)*dt);end
end
C=struct('kind','short_window_resolution_comparison_v1','status','COMPLETED_DIFFERENCES_REPORTED',...
 'outputDir',outDir,'windows',windows,'completedToothWindowJobs',completedJobs,...
 'maximumBaselineForceDifference_N',maximumBaselineForce,'maximumBaselineMomentDifference_Nmm',maximumBaselineMoment,...
 'comparisons',{comparisons},'materialHistoryBeforeLocalWindowsChanged',false,...
 'fullPassComputed',false,'physicalAccuracyValidated',false,...
 'passDepthScope','Local before-state is used as a proxy; this test does not verify full-pass ae_pass. The current force law does not use ae_pass.',...
 'scope','Representative short windows only. Not a proof of convergence over the full axial path or of the pre-window 50us history.');
save(fullfile(outDir,'resolution_curves.mat'),'C','globalCases','-v7');
write_text(fullfile(outDir,'resolution_metrics.json'),jsonencode(C));
fprintf('短窗口对照完成：%s\n',fullfile(outDir,'resolution_metrics.json'));
end

function out=local_case(Q,h0,R,M,range,divisor,widths,progressEvery)
o=R.options;Q.h=unpack(h0,Q);Q.passStart=Q.h;Q.inPass=false;
indices=range(1)+1:range(2);n=numel(indices)*divisor;values=zeros(n,12,numel(widths));
dt=M.forceDt_s/divisor;ix=0;clock=tic;
for k=indices
 ta=M.window_s(1)+(k-1)*M.forceDt_s;tb=M.window_s(1)+k*M.forceDt_s;
 for u=1:divisor
  ix=ix+1;left=ta+(u-1)*(tb-ta)/divisor;t=ta+u*(tb-ta)/divisor;
  if toc(clock)>progressEvery
   fprintf('  第%d齿第%d侧，时间细分%d：%d/%d区间。\n',Q.opts.tooth,Q.opts.side,divisor,ix,n);clock=tic;
  end
  values(ix,1,:)=t;
  face=[Q.opts.tooth Q.opts.side];
  if ~r30b_gate(R.parameters,face,left,t,o.normalDepthLimit_mm),Q.inPass=false;continue;end
  before=r30b_heights(Q);subVolume=0;
  for sub=1:o.substeps
   ts=left+(t-left)*sub/o.substeps;
   if r30b_gate(R.parameters,face,ts,ts,o.normalDepthLimit_mm)
    [Q,v]=r30b_step(Q,ts,o);subVolume=subVolume+v;
   end
  end
  after=r30b_heights(Q);dh=before-after;volume=sum(Q.weightsFine(:).*dh(:));area=sum(Q.weightsFine(dh>0));
  assert(abs(volume-subVolume)<1e-9);
  for c=1:numel(widths)
   oc=o;oc.forceStripWidth_mm=widths(c);
   rows=r30b_force(r30b_strips(Q,before,after,before,t,dt,oc,false),o.coefficients);
   if ~isempty(rows)
    assigned=sum([rows.dV_mm3]);f=o.coefficients;
    for j=1:numel(rows)
     rr=rows(j);heq=rr.dV_mm3/dt/rr.b_mm/rr.vc_mmps;
     expected=f.Kn*rr.b_mm*(heq/f.h_ref)^f.m*(rr.lg_mm/f.l_ref)^f.p;
     assert(abs(rr.Fn_N-expected)<1e-8,'R30R:Law','当前力公式不符合无ae_pass输入的假设。');
    end
   else,assigned=0;end
   assert(abs(assigned-volume)<1e-9);
   s=r30b_sum(rows,o.referenceOutward,t,dt,volume,area);
   values(ix,:,c)=[t volume s.force_N s.normalVector_N s.moment_Nmm area];
  end
 end
end
out=struct('values',values,'finalHeight',r30b_heights(Q));
assert(abs(sum(values(:,2,1))-sum(Q.weightsFine(:).*(h0(:)-out.finalHeight(:))))<1e-8);
end
function r=compare(base,raw,aligned,dtBase,dtRaw,A)
fb=vecnorm(base(:,3:5),2,2);fr=vecnorm(raw(:,3:5),2,2);
nb=vecnorm(base(:,6:8),2,2);nr=vecnorm(raw(:,6:8),2,2);
cb=base(:,11)/1000;cr=raw(:,11)/1000;
eb=-sum(pair_power(base,A))*dtBase;er=-sum(pair_power(raw,A))*dtRaw;
[~,ib]=max(fb);[~,ir]=max(fr);
r=struct('peakForceBase_N',max(fb),'peakForceTrial_N',max(fr),...
 'peakNormalBase_N',max(nb),'peakNormalTrial_N',max(nr),...
 'peakAbsTorqueBase_Nm',max(abs(cb)),'peakAbsTorqueTrial_Nm',max(abs(cr)),...
 'peakForceRelativeChange',rel(max(fr)-max(fb),max(fr)),...
 'peakNormalRelativeChange',rel(max(nr)-max(nb),max(nr)),...
 'peakAbsTorqueRelativeChange',rel(max(abs(cr))-max(abs(cb)),max(abs(cr))),...
 'maxAlignedForceVectorDifferenceRelative',rel(max(vecnorm(base(:,3:5)-aligned(:,3:5),2,2)),max(fr)),...
 'maxAlignedNormalVectorDifferenceRelative',rel(max(vecnorm(base(:,6:8)-aligned(:,6:8),2,2)),max(nr)),...
 'maxAlignedTorqueDifferenceRelative',rel(max(abs(base(:,11)-aligned(:,11)))/1000,max(abs(cr))),...
 'volumeBase_mm3',sum(base(:,2)),'volumeTrial_mm3',sum(raw(:,2)),...
 'volumeRelativeChange',rel(sum(raw(:,2))-sum(base(:,2)),sum(raw(:,2))),...
 'energyBase_J',eb,'energyTrial_J',er,'energyRelativeChange',rel(er-eb,er),...
 'peakForceTimeShift_us',(raw(ir,1)-base(ib,1))*1e6);
end
function p=pair_power(X,A)
F=X(:,3:5);mm=X(:,9:11);tw=-mm+cross(repmat([A.center 0 0],size(X,1),1),F,2);
p=(A.omegaG*mm(:,3)-A.feed*F(:,3)+tw*(A.omegaW*A.axisW(:)))/1000;
end
function x=rel(n,d),if abs(d)<1e-15,if abs(n)<1e-15,x=0;else,x=Inf;end;else,x=n/abs(d);end;end
function a=add_values(a,b)
if isempty(a),a=b;else,assert(max0(abs(a(:,1)-b(:,1)))<1e-10);a(:,2:end)=a(:,2:end)+b(:,2:end);end
end
function h=unpack(H,Q),h=permute(reshape(H,Q.n,Q.ns,Q.n,Q.nz),[2 4 1 3]);end
function x=max0(v),if isempty(v),x=0;else,x=max(v(:));end;end
function s=core_signature(core)
files=dir(fullfile(core,'*.m'));s='';for j=1:numel(files),s=[s files(j).name char(10) fileread(fullfile(core,files(j).name))];end
end
function write_csv(path,headers,data)
fid=fopen(path,'w');assert(fid>0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',strjoin(headers,','));
fprintf(fid,[repmat('%.17g,',1,size(data,2)-1) '%.17g\n'],data');
end
function write_text(path,value),fid=fopen(path,'w','n','UTF-8');assert(fid>0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',value);end
function plot_comparison(B,F,S,w,out,start)
h=figure('Visible','off','Color','w','Units','centimeters','Position',[2 2 16 6.8],'PaperPositionMode','auto');
CN='\fontname{SimSun}';EN='\fontname{Times New Roman}';
X={B,F,S};colors=[0 .145 1;1 .11 0;0 0 0];
for j=1:2
 ax=axes('Parent',h,'Position',[.09+(j-1)*.51 .30 .35 .64],'FontName','Times New Roman','FontSize',10,...
  'Box','on','TickDir','in','LineWidth',.8);hold(ax,'on');
 for c=1:3
  d=X{c};y=vecnorm(d(:,3:5),2,2);if j==2,y=d(:,11)/1000;end
  plot(ax,(d(:,1)-start)*1000,y,'Color',colors(c,:),'LineWidth',.65);
 end
 xlabel(ax,[CN '局部时间' EN '(ms)'],'Interpreter','tex');
 if j==1,ylabel(ax,[CN '合力' EN '(N)'],'Interpreter','tex');else,ylabel(ax,[CN '转矩' EN '(N\cdot m)'],'Interpreter','tex');end
 legend(ax,{[CN '基准'],[CN '时间步减半'],[CN '条带减半']},'Box','off','Location','northeast','FontSize',8.5,'Interpreter','tex');
 yl=ylim(ax);ylim(ax,[yl(1),yl(2)+.45*diff(yl)]);
 cap='(a)合力时程';if j==2,cap='(b)C轴载荷转矩';end
 annotation(h,'textbox',[.045+(j-1)*.51 .04 .44 .07],'String',cap,'FontName','SimSun','FontSize',10.5,'EdgeColor','none','HorizontalAlignment','center');
end
name=fullfile(out,sprintf('w%02d_resolution',w));savefig(h,[name '.fig']);
try,exportgraphics(h,[name '.pdf'],'ContentType','vector','BackgroundColor','white');exportgraphics(h,[name '.png'],'Resolution',400,'BackgroundColor','white');
catch,set(h,'PaperUnits','centimeters','PaperSize',[16 6.8],'PaperPosition',[0 0 16 6.8]);print(h,[name '.pdf'],'-dpdf','-painters');print(h,[name '.png'],'-dpng','-r400');end
close(h);
end
