function R30B=r30b_run(parameterFile,userOptions)
% New dual-flank solver: parameters -> t=0 history -> common-time forces.
% No old MAT results, interpolation history, translated cycles or fixed phase
% window are read. Shared analytic geometry is namespaced r30b_*.
if nargin<1||isempty(parameterFile),parameterFile=fullfile(fileparts(mfilename('fullpath')),'..','parameters.json');end
if nargin<2,userOptions=struct();end
o=r30b_options(userOptions);P=r30b_read_json(parameterFile);validate_parameters(P);
window=r30b_window(P,o);nt=round(diff(window)/o.dt_s);nf=size(o.faces,1);
assert(nt>=1,'R30B:Window','No output intervals');
if ~isfolder(o.outputDir),mkdir(o.outputDir);end
if o.saveFrames&&~isfolder(fullfile(o.outputDir,'frames')),mkdir(fullfile(o.outputDir,'frames'));end
files=dir(fullfile(fileparts(mfilename('fullpath')),'*.m'));signature='';
for k=1:numel(files),signature=[signature files(k).name char(10) fileread(fullfile(files(k).folder,files(k).name))];end %#ok<AGROW>
meta=struct('engineVersion','r30_dual_from_zero_v2','parameters',P,'options',o,'window_s',window,'codeSignature',signature);
cp=fullfile(o.outputDir,'r30b_checkpoint.mat');resultFile=fullfile(o.outputDir,'r30b_result.mat');
if o.resume&&isfile(cp)
 D=load(cp,'S');assert(isfield(D,'S')&&isfield(D.S,'meta')&&isequaln(D.S.meta,meta),...
  'R30B:CheckpointMismatch','Checkpoint code/parameters/options differ. Use a new output directory for the revised run.');S=D.S;
 if strcmp(S.phase,'done'),D=load(resultFile,'R30B');R30B=D.R30B;return;end
 fprintf('Resuming new-engine phase %s; no old material result was loaded.\n',S.phase);
else
 S=struct('meta',meta,'phase','history','historyNext',0,'windowNext',1,'Q',{{}},...
  'historyEvaluations',zeros(nf,1),'historyRejectedSamples',zeros(nf,1),...
  'facePoses',{cell(nt,nf)},'poses',{cell(nt,1)},'frameFiles',{{}},'frameIndices',[],...
  'maximumPartitionError_mm3',0,'maximumSubstepBalanceError_mm3',0,'initialFine',{{}},...
  'lastFrameFine',{{}},'lastFrameTime_s',window(1),'windowVolume_mm3',zeros(nf,1));
 for f=1:nf,S.Q{f}=r30b_init(P,o,o.faces(f,:));end
end
clock=tic;checkpointClock=tic;progressClock=tic;
fprintf('New calculation: %d faces, %dx%d cells x %d^2 points; window %.8f..%.8f s.\n',nf,o.ns,o.nz,o.order,window);
nh=ceil(window(1)/o.historyDt_s);block=max(1,floor(o.historyBlock_s/o.historyDt_s));
while strcmp(S.phase,'history')&&S.historyNext<=nh
 first=S.historyNext;last=min(nh,first+block-1);ta=min(first*o.historyDt_s,window(1));tb=min(last*o.historyDt_s,window(1));
 possible=false(nf,1);for f=1:nf,possible(f)=r30b_gate(P,o.faces(f,:),ta,tb,o.normalDepthLimit_mm);end
 if ~any(possible)
  for f=1:nf,S.Q{f}.inPass=false;S.historyRejectedSamples(f)=S.historyRejectedSamples(f)+last-first+1;end
  S.historyNext=last+1;
 else
  for h=first:last
   t=min(h*o.historyDt_s,window(1));
   for f=1:nf
    Q=S.Q{f};can=possible(f)&&r30b_gate(P,o.faces(f,:),t,t,o.normalDepthLimit_mm);
    if can
     if ~Q.inPass,Q.passStart=Q.h;Q.passId=Q.passId+1;Q.inPass=true;end
     [Q,v]=r30b_step(Q,t,o);Q.historyVolume_mm3=Q.historyVolume_mm3+v;S.historyEvaluations(f)=S.historyEvaluations(f)+1;
    else,Q.inPass=false;S.historyRejectedSamples(f)=S.historyRejectedSamples(f)+1;end
    S.Q{f}=Q;
   end
   S.historyNext=h+1;
  end
 end
 if toc(progressClock)>=o.progressEvery_s
  fprintf('History from zero: t=%.6f / %.6f s; elapsed %.1f s\n',min((S.historyNext-1)*o.historyDt_s,window(1)),window(1),toc(clock));progressClock=tic;
 end
 if toc(checkpointClock)>=o.checkpointEvery_s,r30b_atomic(cp,'S',S);checkpointClock=tic;end
end
if strcmp(S.phase,'history')
 initialState=struct('engineVersion',meta.engineVersion,'t_s',window(1),'faces',{{}});
 for f=1:nf
  Q=S.Q{f};H=r30b_heights(Q);S.initialFine{f}=H;S.lastFrameFine{f}=H;
  temp=Q;temp.h=Q.passStart;pass=r30b_heights(temp);
  initialState.faces{f}=struct('tooth_id',o.faces(f,1),'flank_side',o.faces(f,2),...
   'height_mm',H,'passStartHeight_mm',pass,'passId',Q.passId,'stock_nominal_mm',Q.stock,...
   'support_s_mm',Q.support_s,'support_z_mm',Q.support_z,...
   'sample_s_mm',reshape(Q.s',[],1),'sample_z_mm',reshape(Q.z',[],1));
 end
 r30b_atomic(fullfile(o.outputDir,'initial_state.mat'),'initialState',initialState);
 S.phase='window';r30b_atomic(cp,'S',S);checkpointClock=tic;
end
for k=S.windowNext:nt
 ta=window(1)+(k-1)*o.dt_s;t=window(1)+k*o.dt_s;capture=o.saveFrames&&(mod(k-1,o.frameStride)==0||k==nt);
 frame=struct('engineVersion',meta.engineVersion,'index',k,'t_start_s',ta,'t_s',t,'dt_s',o.dt_s,...
  'toolSampleTimes_s',ta+(1:o.substeps)*(t-ta)/o.substeps,'faces',{{}});
 allRows=struct([]);allVolume=0;allArea=0;
 for f=1:nf
  Q=S.Q{f};before=r30b_heights(Q);subVolume=0;subRecords=cell(o.substeps,1);
  for sub=1:o.substeps
   ts=ta+(t-ta)*sub/o.substeps;
   if r30b_gate(P,o.faces(f,:),ts,ts,o.normalDepthLimit_mm)
    if ~Q.inPass,Q.passStart=Q.h;Q.passId=Q.passId+1;Q.inPass=true;end
    if capture,[Q,v,~,change]=r30b_step(Q,ts,o);subRecords{sub}=change;
    else,[Q,v]=r30b_step(Q,ts,o);end
    subVolume=subVolume+v;
   else
    Q.inPass=false;
    if capture,subRecords{sub}=struct('t_s',ts,'index',zeros(0,1),'before_mm',zeros(0,1),'after_mm',zeros(0,1));end
   end
  end
  after=r30b_heights(Q);temp=Q;temp.h=Q.passStart;pass=r30b_heights(temp);dh=before-after;active=dh>0;
  volume=sum(Q.weightsFine(:).*dh(:));area=sum(Q.weightsFine(active));
  S.maximumSubstepBalanceError_mm3=max(S.maximumSubstepBalanceError_mm3,abs(volume-subVolume));
  rows=r30b_strips(Q,before,after,pass,t,o.dt_s,o,capture);rows=r30b_force(rows,o.coefficients);
  assigned=0;if ~isempty(rows),assigned=sum([rows.dV_mm3]);end
  difference=abs(assigned-volume);assert(difference<1e-9,'R30B:VolumeBalance','Force strips and material change disagree');
  S.maximumPartitionError_mm3=max(S.maximumPartitionError_mm3,difference);
  S.facePoses{k,f}=r30b_sum(rows,o.referenceOutward,t,o.dt_s,volume,area);S.windowVolume_mm3(f)=S.windowVolume_mm3(f)+volume;
  if ~isempty(rows)
   if isempty(allRows),allRows=rows;else,allRows=[allRows;rows(:)];end %#ok<AGROW>
  end
  allVolume=allVolume+volume;allArea=allArea+area;
  if capture
   changed=find(after<S.lastFrameFine{f});ids=find(active);[i,j]=ind2sub(size(after),ids);
   removal=struct('index',ids,'before_mm',before(ids),'after_mm',after(ids),'passStart_mm',pass(ids),...
    'rectangles_sz_mm',[Q.support_s(i),Q.support_s(i+1),Q.support_z(j),Q.support_z(j+1)],...
    'sample_sz_mm',[reshape_sample(Q.s,i),reshape_sample(Q.z,j)]);
   frame.faces{f}=struct('tooth_id',o.faces(f,1),'flank_side',o.faces(f,2),'passId',Q.passId,...
    'stateUpdates',struct('since_t_s',S.lastFrameTime_s,'index',changed,'height_mm',after(changed)),...
    'removal',removal,'substepUpdates',{subRecords},'strips',rows,'summary',S.facePoses{k,f});
   S.lastFrameFine{f}=after;
  end
  S.Q{f}=Q;
 end
 S.poses{k}=r30b_sum(allRows,o.referenceOutward,t,o.dt_s,allVolume,allArea);
 if capture
  frame.combined=S.poses{k};name=sprintf('frames/frame_%06d.mat',k);
  r30b_atomic(fullfile(o.outputDir,name),'frame',frame);S.frameFiles{end+1}=name;S.frameIndices(end+1)=k;S.lastFrameTime_s=t;
 end
 S.windowNext=k+1;
 if toc(progressClock)>=o.progressEvery_s
  fprintf('Common-time forces: %d/%d, t=%.7f s; strips=%d, Fn(pure)=%.6g N\n',k,nt,t,S.poses{k}.nForceCells,S.poses{k}.Fn_pure_resultant_N);progressClock=tic;
 end
 if toc(checkpointClock)>=o.checkpointEvery_s,r30b_atomic(cp,'S',S);checkpointClock=tic;end
end
finalState=struct('t_s',window(2),'faces',{{}});stateVolume=zeros(nf,1);historyClosure=zeros(nf,1);historyVolume=zeros(nf,1);
for f=1:nf
 Q=S.Q{f};H=r30b_heights(Q);finalState.faces{f}=struct('tooth_id',o.faces(f,1),'flank_side',o.faces(f,2),'height_mm',H);
 stateVolume(f)=sum(Q.weightsFine(:).*(S.initialFine{f}(:)-H(:)));
 beforeVolume=sum(Q.weightsFine(:).*(Q.stock-S.initialFine{f}(:)));historyVolume(f)=Q.historyVolume_mm3;
 historyClosure(f)=abs(beforeVolume-Q.historyVolume_mm3);
end
closure=abs(stateVolume-S.windowVolume_mm3);
assert(max([closure;historyClosure])<1e-8,'R30B:MaterialBalance','Material state and integrated volume disagree');
wholePassEnclosed=true;
for f=1:nf
 wholePassEnclosed=wholePassEnclosed&&~r30b_gate(P,o.faces(f,:),window(1),window(1),o.normalDepthLimit_mm)&&...
  ~r30b_gate(P,o.faces(f,:),window(2),window(2),o.normalDepthLimit_mm);
end
r30b_atomic(fullfile(o.outputDir,'final_state.mat'),'finalState',finalState);
R30B=struct('engineVersion',meta.engineVersion,'parameters',P,'options',o,'faces',o.faces,...
 'window_s',window,'timeEdges_s',window(1)+(0:nt)*o.dt_s,'poses',vertcat(S.poses{:}),...
 'facePoses',{cell(nf,1)},'frameFiles',{S.frameFiles},'frameIndices',S.frameIndices,...
 'initialStateFile','initial_state.mat','finalStateFile','final_state.mat','forceCoefficients',o.coefficients,...
 'fullGearComputed',false,'momentOrigin_mm',[0 0 0],'oldResultsUsed',false,'historyFromZero',true,...
 'wholePassEnclosed',wholePassEnclosed,...
 'historyVolume_mm3',historyVolume,'windowVolume_mm3',S.windowVolume_mm3,...
 'validation',struct('historyClosure_mm3',historyClosure,'windowClosure_mm3',closure,...
  'maximumPartitionError_mm3',S.maximumPartitionError_mm3,'maximumSubstepBalanceError_mm3',S.maximumSubstepBalanceError_mm3),...
 'historyEvaluations',S.historyEvaluations,'historyRejectedSamples',S.historyRejectedSamples,...
 'geometryRepresentation','positive-weight numerical support rectangles, not exact continuous intersection',...
 'forceTimeMeaning','interval volume-rate force with endpoint direction/geometry quadrature',...
 'coefficientMeaning','candidate macro law, Kn in N/mm, not an original Dietz fit',...
 'scope','Selected independent involute flanks above base circle; no root fillet or grain model');
for f=1:nf,R30B.facePoses{f}=vertcat(S.facePoses{:,f});end
r30b_atomic(resultFile,'R30B',R30B);r30b_outputs(R30B,o.outputDir);
S.phase='done';r30b_atomic(cp,'S',S);
fprintf('Finished from-zero calculation: %d faces, %d common intervals; result: %s\n',nf,nt,resultFile);
end

function sample=reshape_sample(q,index)
values=reshape(q',[],1);sample=values(index);
end

function validate_parameters(P)
assert(all(isfield(P,{'p','installation','cfg'})),'R30B:Parameters','Use geometry/kinematics parameters JSON, not a saved material result');
assert(P.p.gear.teeth==30&&P.p.wheel.starts==1&&abs(P.p.gear.helix_deg)<1e-12,...
 'R30B:Scope','This implementation supports the current single-start spur geometry');
A=P.installation;w=A.axisW(:);y=A.axisY(:);
assert(abs(norm(w)-1)<1e-10&&abs(norm(y)-1)<1e-10&&abs(dot(w,y))<1e-10&&abs(w(1))+abs(y(1))<1e-12,...
 'R30B:Installation','Wheel axes must be orthonormal and perpendicular to world X');
assert(A.omegaG~=0&&A.omegaW~=0&&P.cfg.stock_um>0&&P.cfg.stock_um<=500,...
 'R30B:Parameters','Invalid motion or initial normal stock');
rp=P.p.gear.module_mm*P.p.gear.teeth/2;rb=rp*cos(P.p.gear.pressure_deg*pi/180);
assert(abs(P.p.calc.gear.pitchRadius_mm-rp)<1e-8&&abs(P.p.calc.gear.baseDiameter_mm/2-rb)<1e-8,...
 'R30B:GearConsistency','Basic gear inputs and derived radii disagree; regenerate parameters before calculating');
assert(A.lead>0&&A.helicalP>0&&A.rRoot>0&&A.tipNew>A.rRoot&&A.axialLength>0&&...
 abs(A.lead-2*pi*A.helicalP)<1e-8,'R30B:WheelConsistency','Invalid or inconsistent ideal worm geometry');
assert(isfield(A,'profileType')&&strcmp(A.profileType,'involute_helicoid'),'R30B:Profile','Require the corrected ideal involute-helicoid profile');
end
