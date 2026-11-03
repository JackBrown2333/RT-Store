function report=r30b_verify(resultFile,outputFile)
% Independent bookkeeping of summaries, saved footprints and sparse replay.
if nargin<2,outputFile='';end
D=load(resultFile,'R30B');R=D.R30B;root=fileparts(resultFile);
assert(strcmp(R.engineVersion,'r30_dual_from_zero_v2')&&R.historyFromZero&&~R.oldResultsUsed,'R30B:Version','Wrong engine or legacy-state reuse');
nf=size(R.faces,1);nt=numel(R.poses);maxVector=0;maxVolume=0;maxLaw=0;maxAnnotation=0;maxSubstepStateError=0;
for k=1:nt
 F=[0 0 0];N=F;M=F;v=0;fn=0;
 for f=1:nf
  s=R.facePoses{f}(k);assert(abs(s.t_s-R.poses(k).t_s)<1e-12&&abs(s.dt_s-R.options.dt_s)<1e-12);
  F=F+s.force_N;N=N+s.normalVector_N;M=M+s.moment_Nmm;v=v+s.dV_mm3;fn=fn+s.Fn_magnitude_sum_N;
 end
 a=R.poses(k);maxVector=max(maxVector,max(abs([F-a.force_N,N-a.normalVector_N,M-a.moment_Nmm,fn-a.Fn_magnitude_sum_N])));
 maxVolume=max(maxVolume,abs(v-a.dV_mm3));assert(abs(norm(N)-a.Fn_pure_resultant_N)<1e-9);
end
D=load(fullfile(root,R.initialStateFile),'initialState');initial=D.initialState;heights=cell(nf,1);weights=cell(nf,1);
for f=1:nf
 heights{f}=initial.faces{f}.height_mm;weights{f}=diff(initial.faces{f}.support_s_mm)*diff(initial.faces{f}.support_z_mm)';
end
for j=1:numel(R.frameFiles)
 D=load(fullfile(root,R.frameFiles{j}),'frame');frame=D.frame;k=frame.index;
 for f=1:nf
  x=frame.faces{f};u=x.stateUpdates;h=heights{f};
  assert(numel(unique(u.index))==numel(u.index)&&all(u.index>=1&u.index<=numel(h)),'R30B:Replay','Invalid material update indices');
  assert(all(u.height_mm<=h(u.index)+1e-13),'R30B:Replay','Replay restores removed material');h(u.index)=u.height_mm;heights{f}=h;
  cut=x.removal;assert(all(cut.before_mm>=cut.after_mm)&&all(cut.passStart_mm>=cut.before_mm-1e-12));
  rollback=h;microVolume=0;
  for sub=numel(x.substepUpdates):-1:1
   step=x.substepUpdates{sub};ids=step.index;
   if ~isempty(ids),maxSubstepStateError=max(maxSubstepStateError,max(abs(rollback(ids)-step.after_mm)));end
   assert(numel(unique(ids))==numel(ids)&&all(step.before_mm>=step.after_mm));
   microVolume=microVolume+sum(weights{f}(ids).*(step.before_mm-step.after_mm));rollback(ids)=step.before_mm;
  end
  if ~isempty(cut.index),maxSubstepStateError=max(maxSubstepStateError,max(abs(rollback(cut.index)-cut.before_mm)));end
  amount=sum(weights{f}(cut.index).*(cut.before_mm-cut.after_mm));
  maxVolume=max(maxVolume,abs(amount-microVolume));
  maxVolume=max(maxVolume,abs(amount-R.facePoses{f}(k).dV_mm3));
  rows=x.strips;local=0;
  for n=1:numel(rows)
   r=rows(n);c=R.forceCoefficients;heq=r.dV_mm3/(r.dt_s*r.b_mm*r.vc_mmps);
   force=c.Kn*r.b_mm*(heq/c.h_ref)^c.m*(r.lg_mm/c.l_ref)^c.p+c.background_pressure*r.contact_area_mm2;
   maxLaw=max(maxLaw,max(abs([heq-r.heq_mm,force-r.Fn_N,c.mu*force-r.Ft_N])));local=local+r.dV_mm3;
   g=r.geometry;maxAnnotation=max(maxAnnotation,max(abs([norm(diff(g.b_segment_world_mm))-r.b_mm,...
    norm(diff(g.lg_equivalent_segment_world_mm))-r.lg_mm,...
    norm(diff(g.ae_pass_segment_world_mm))-r.ae_pass_mm,...
    norm(diff(g.dh_step_segment_world_mm))-r.dh_step_mm])));
  end
  maxVolume=max(maxVolume,abs(local-amount));
 end
end
replayError=NaN;
if ~isempty(R.frameFiles)
 D=load(fullfile(root,R.finalStateFile),'finalState');replayError=0;
 for f=1:nf,replayError=max(replayError,max(abs(heights{f}(:)-D.finalState.faces{f}.height_mm(:))));end
 assert(replayError<1e-12,'R30B:Replay','Sparse updates do not reproduce final material state');
end
assert(maxVector<1e-8&&maxVolume<1e-9&&maxLaw<1e-7&&maxAnnotation<1e-9&&maxSubstepStateError<1e-12,...
 'R30B:Verification','Saved numerical relationships disagree');
report=struct('status','PASS','engineVersion',R.engineVersion,'faces',R.faces,'checkedIntervals',nt,...
 'checkedGeometryFrames',numel(R.frameFiles),'maxSummaryError',maxVector,'maxVolumeError_mm3',maxVolume,...
 'maxConstitutiveError',maxLaw,'maxDimensionAnchorError_mm',maxAnnotation,'materialReplayError_mm',replayError,...
 'maxSubstepStateError_mm',maxSubstepStateError,...
 'physicalCalibrationValidated',false,'productionConvergenceValidated',false);
if ~isempty(outputFile),fid=fopen(outputFile,'w','n','UTF-8');assert(fid>0);fprintf(fid,'%s',jsonencode(report));fclose(fid);end
disp(report);
end
