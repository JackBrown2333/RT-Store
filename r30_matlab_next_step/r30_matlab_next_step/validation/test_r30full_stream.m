function report=test_r30full_stream(sourceMat,coreDir,workDir)
% Numerical software fixture only. Geometry phase is shifted to a contact pose.
% requireClearStart=false is used ONLY in this synthetic regression test.
addpath(coreDir);D=load(sourceMat,'M');
M=struct('parameters',D.M.parameters,'options',D.M.options,'forceDt_s',D.M.forceDt_s,'coefficients',D.M.coefficients);clear D;
original=M.parameters;
P=original;A=P.installation;period=2*pi/abs(A.omegaW);shift=round(15.014/period)*period;
P.installation.phi0=A.phi0+A.omegaG*shift;P.installation.feedStart=A.feedStart+A.feed*shift;
M.parameters=P;M.options.ns=12;M.options.nz=24;M.options.order=2;
if ~isfolder(workDir),mkdir(workDir);end
save(fullfile(workDir,'fixture.mat'),'M','-v7');dt=M.forceDt_s;nt=16;
plan=struct('parameters',P,'modelTimeRange_s',[0 nt*dt]);save(fullfile(workDir,'plan.mat'),'plan','-v7');
o=struct('toothIds',[17 18],'chunkIntervals',5,'plot',false,'requireClearStart',false,'maxChunksPerCall',1);
out=fullfile(workDir,'resumed');R=run_r30_full_force(fullfile(workDir,'fixture.mat'),coreDir,out,fullfile(workDir,'plan.mat'),o);
assert(~R.completed&&strcmp(R.status,'PAUSED_AT_SAVED_CHUNK'));
o.maxChunksPerCall=Inf;R=run_r30_full_force(fullfile(workDir,'fixture.mat'),coreDir,out,fullfile(workDir,'plan.mat'),o);
assert(R.completed&&~R.fullGearComputed&&~R.machiningCoverageAccepted);
co=r30b_options(struct('ns',12,'nz',24,'order',2,'dt_s',dt,'substeps',M.options.substeps,...
 'forceStripWidth_mm',M.options.forceStripWidth_mm/2,'coefficients',M.coefficients,'plot',false));
reference=zeros(nt,16);stateError=0;
for tooth=[17 18]
 expected=zeros(nt,13,2);
 for f=1:2
  Q=r30b_init(P,co,[tooth f]);
  for k=1:nt
   before=r30b_heights(Q);ta=(k-1)*dt;t=k*dt;
   for sub=1:co.substeps
    ts=ta+(t-ta)*sub/co.substeps;
    if r30b_gate(P,[tooth f],ts,ts,co.normalDepthLimit_mm)
     if ~Q.inPass,Q.passStart=Q.h;Q.passId=Q.passId+1;Q.inPass=true;end
     [Q,~]=r30b_step(Q,ts,co);
    else,Q.inPass=false;end
   end
   after=r30b_heights(Q);temp=Q;temp.h=Q.passStart;pass=r30b_heights(temp);dh=before-after;
   volume=sum(Q.weightsFine(:).*dh(:));area=sum(Q.weightsFine(dh>0));
   rows=r30b_force(r30b_strips(Q,before,after,pass,t,dt,co,false),co.coefficients);a=r30b_sum(rows,co.referenceOutward,t,dt,volume,area);
   expected(k,:,f)=[volume area a.force_N a.normalVector_N a.moment_Nmm a.Fn_magnitude_sum_N a.nForceCells];
  end
  D=load(fullfile(out,sprintf('tooth_%02d',tooth),'final_state.mat'),'finalState');
  stateError=max(stateError,max(abs(reshape(D.finalState.faces{f}.height_mm-r30b_heights(Q),[],1))));
 end
 actual=r30full_read(out,1:nt,tooth);assert(max(abs(actual(:)-expected(:)))<1e-7);
 T=sum(expected,3);reference(:,1:13)=reference(:,1:13)+T;
 reference(:,14)=reference(:,14)+(T(:,1)>0);reference(:,15)=reference(:,15)+sum(expected(:,1,:)>0,3);
 reference(:,16)=reference(:,16)+(vecnorm(T(:,3:5),2,2)>.01);
end
actual=r30full_read(out,1:nt);err=max(abs(actual(:)-reference(:)));assert(err<1e-7&&stateError<1e-14&&R.rawPeakForce_N>0);
pick=[1 5 6 16];assert(isequal(r30full_read(out,pick),actual(pick,:)));
rejected=false;bad=o;bad.ns=13;
try,run_r30_full_force(fullfile(workDir,'fixture.mat'),coreDir,out,fullfile(workDir,'plan.mat'),bad);
catch e,rejected=strcmp(e.identifier,'R30FULL:Resume');end
assert(rejected);
% An entirely empty-contact run exercises zero-row chunks and endpoints.
M.parameters=original;plan.parameters=original;save(fullfile(workDir,'empty_fixture.mat'),'M','-v7');save(fullfile(workDir,'empty_plan.mat'),'plan','-v7');
o.toothIds=17;o.requireClearStart=true;emptyOut=fullfile(workDir,'empty');
E=run_r30_full_force(fullfile(workDir,'empty_fixture.mat'),coreDir,emptyOut,fullfile(workDir,'empty_plan.mat'),o);
empty=r30full_read(emptyOut,1:nt);assert(~any(empty(:))&&E.allEndpointStockSamplesClear&&E.rawPeakForce_N==0);
report=struct('syntheticOnly',true,'toothCount',2,'flankCount',4,'intervals',nt,'chunkIntervals',5,...
 'streamMatchesUnchunkedCore',true,'maxRawDifference',err,'maxFinalStateDifference_mm',stateError,...
 'resumeAndContinuousMaterialState',true,'parameterChangeRejected',true,'noncontiguousRawRead',true,...
 'emptySparseChunks',true,'all60FaceFullStrokeExecutedHere',false,'nativeMatlabGraphicsTested',false);
fid=fopen(fullfile(workDir,'software_test_report.json'),'w');fprintf(fid,'%s',jsonencode(report));fclose(fid);disp(report);
end
