function report=r30b_test(parameterFile,outputFile)
% Actual independent-flank calculation from t=0 on a SMALL verification grid.
% This is a numerical implementation test, not production convergence.
if nargin<1,parameterFile=fullfile(fileparts(mfilename('fullpath')),'..','parameters.json');end
if nargin<2,outputFile='';end
P=r30b_read_json(parameterFile);o=r30b_options(struct('ns',10,'nz',20));
rootError=0;screenError=0;repeatedError=0;gateRejects=0;hitCount=zeros(2,1);
for side=1:2
 Q=r30b_init(P,o,[17 side]);assert(all(Q.w>0)&&abs(sum(Q.w)-1)<1e-14);
 assert(abs(sum(Q.weightsFine(:))-Q.smax*P.p.gear.faceWidth_mm)<1e-10);
 [ss,zz]=ndgrid(reshape(Q.s',[],1),reshape(Q.z',[],1));
 for t=linspace(14.94,15.08,71)
  initial=Q.stock+zeros(size(ss));hit=r30b_level(ss,zz,initial,t,Q.p,Q.A,17,side)<0;
  if any(hit(:))
   s=ss(hit);z=zz(hit);before=initial(hit);a=r30b_normal_update(s,z,before,t,Q.p,Q.A,17,side,Q.stock);
   lo=-ones(size(a));hi=before;assert(all(r30b_level(s,z,lo,t,Q.p,Q.A,17,side)>0));
   for k=1:48
    mid=(lo+hi)/2;inside=r30b_level(s,z,mid,t,Q.p,Q.A,17,side)<=0;
    hi(inside)=mid(inside);lo(~inside)=mid(~inside);
   end
   rootError=max(rootError,max(abs(a-hi)));again=r30b_normal_update(s,z,a,t,Q.p,Q.A,17,side,Q.stock);
   repeatedError=max(repeatedError,max(abs(again-a)));hitCount(side)=hitCount(side)+numel(a);
  end
 end
 for t=[14.98 15.005 15.03 15.05]
  [a,va]=r30b_step(Q,t,o,true);[b,vb]=r30b_step(Q,t,o,false);
  screenError=max(screenError,max(abs(a.h(:)-b.h(:))));assert(abs(va-vb)<1e-10);
 end
 for ta=linspace(0,15.1,181)
  tb=ta+.002;
  if ~r30b_gate(P,[17 side],ta,tb,o.normalDepthLimit_mm)
   gateRejects=gateRejects+1;
   for t=linspace(ta,tb,3)
    for h=[Q.stock 0 -o.normalDepthLimit_mm]
     v=r30b_level(ss,zz,h,t,Q.p,Q.A,17,side);assert(all(v(:)>=-1e-10),'R30B:GateTest','Enclosing gate rejected a sampled tool intersection');
    end
   end
  end
 end
end
assert(all(hitCount>0)&&rootError<1e-9&&screenError<1e-11&&repeatedError<1e-9&&gateRejects>0);
automatic=r30b_window(P,r30b_options());
for side=1:2
 assert(~r30b_gate(P,[17 side],automatic(1),automatic(1),1)&&~r30b_gate(P,[17 side],automatic(2),automatic(2),1));
end
folder=tempname;
small=struct('ns',12,'nz',24,'order',2,'window_s',[14.99 15.025],'dt_s',.001,...
 'substeps',1,'historyDt_s',.001,'forceStripWidth_mm',.3,'frameStride',7,...
 'plot',false,'resume',true,'outputDir',folder);
R=r30b_run(parameterFile,small);verification=r30b_verify(fullfile(folder,'r30b_result.mat'));
assert(isequal(R.faces,[17 1;17 2])&&R.historyFromZero&&~R.oldResultsUsed);
for side=1:2,assert(any([R.facePoses{side}.dV_mm3]>0)&&any([R.facePoses{side}.Fn_pure_resultant_N]>0));end
N1=vertcat(R.facePoses{1}.normalVector_N);N2=vertcat(R.facePoses{2}.normalVector_N);
assert(any(sum(N1.*N2,2)<0),'R30B:Directions','Independent flank load directions were not distinguished');
again=r30b_run(parameterFile,small);assert(isequaln(R.poses,again.poses),'R30B:Resume','Completed checkpoint changed the result');
bad=small;bad.coefficients=struct('m',.65);caught=false;
try,r30b_run(parameterFile,bad);catch e,caught=strcmp(e.identifier,'R30B:CheckpointMismatch');end
assert(caught,'R30B:CheckpointTest','Changed model must not reuse a checkpoint');
jsonFile=fullfile(folder,'animation_subset.json');r30b_export_json(fullfile(folder,'r30b_result.mat'),jsonFile,2);
J=jsondecode(fileread(jsonFile));D=load(fullfile(folder,'final_state.mat'),'finalState');
for f=1:2
 initialFace=get_item(J.initialState.faces,f);h=initialFace.height_mm;
 for k=1:numel(J.frames)
  frame=get_item(J.frames,k);ff=get_item(frame.faces,f);u=ff.stateUpdates;h(u.index)=u.height_mm;
 end
 assert(max(abs(h(:)-D.finalState.faces{f}.height_mm(:)))<1e-12,'R30B:JSONReplay','Exported subsampled material updates are incomplete');
end
% Manufactured cancellation checks time averaging BEFORE vector norms.
a=R;b=R;a.window_s=[0 .2];b.window_s=a.window_s;a.options.dt_s=.2;b.options.dt_s=.1;
a.poses=R.poses(1);a.poses.t_s=.2;a.poses.dt_s=.2;a.poses.dV_mm3=3;
a.poses.force_N=[0 0 0];a.poses.normalVector_N=[0 0 0];b.poses=repmat(a.poses,2,1);
for k=1:2
 b.poses(k).t_s=k*.1;b.poses(k).dt_s=.1;b.poses(k).dV_mm3=k;
 b.poses(k).force_N=[(-1)^k 0 0];b.poses(k).normalVector_N=[0 (-1)^k 0];
end
r30b_atomic(fullfile(folder,'compare_a.mat'),'R30B',a);r30b_atomic(fullfile(folder,'compare_b.mat'),'R30B',b);
comparison=r30b_compare_runs(fullfile(folder,'compare_a.mat'),fullfile(folder,'compare_b.mat'));
assert(comparison.maxForceVectorDifference_N==0&&comparison.maxPureNormalNormDifference_overFinePeak==0&&comparison.totalVolumeDifference_mm3==0);
report=struct('status','PASS','runtime',version,'nativeMatlabTested',exist('OCTAVE_VERSION','builtin')==0,...
 'fromZeroDualFlankRunTested',true,'coarseGrid',[12 24 2],'historyDt_s',.001,...
 'rootVsBisectionMaxError_mm',rootError,'spatialScreenMaxStateError_mm',screenError,...
 'repeatPoseMaxError_mm',repeatedError,'enclosingGateRejectedIntervalsTested',gateRejects,...
 'automaticWindow_s',automatic,'verification',verification,'checkpointGuardTested',true,...
 'subsampledJsonReplayTested',true,'commonBinVectorAverageTested',true,'productionConvergenceValidated',false,'experimentalValidation',false);
if ~isempty(outputFile),fid=fopen(outputFile,'w','n','UTF-8');assert(fid>0);fprintf(fid,'%s',jsonencode(report));fclose(fid);end
% Remove only this test's own private files, never a user-specified folder.
items=dir(fullfile(folder,'frames','*.mat'));for k=1:numel(items),delete(fullfile(folder,'frames',items(k).name));end
rmdir(fullfile(folder,'frames'));items=dir(folder);
for k=1:numel(items),if ~items(k).isdir,delete(fullfile(folder,items(k).name));end;end;rmdir(folder);
disp(report);
end

function x=get_item(array,k)
if iscell(array),x=array{k};else,x=array(k);end
end
