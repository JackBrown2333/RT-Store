function plan=run_r30_full_preflight(sourceMat,coreDir,outDir,user)
% Full-gear nominal-stock geometry scan, NOT a force or material-history run.
% Phase folding is used only for the identical, undeformed nominal stock.
% A future force solver must maintain 60 independent material histories.
if nargin<4,user=struct();end
o=struct('profileSamples',25,'axialSamples',41,'phaseSamples',12,...
 'axialStep_mm',.5,'clearanceLevel',.01,'tailMargin_mm',1,...
 'maxExtension_mm',20,'progressEvery_s',5);
names=fieldnames(user);for j=1:numel(names),assert(isfield(o,names{j}));o.(names{j})=user.(names{j});end
assert(o.profileSamples>=3&&o.axialSamples>=3&&o.phaseSamples>=4&&o.axialStep_mm>0);
addpath(coreDir);D=load(sourceMat,'M');M=D.M;P=M.parameters;A=P.installation;
assert(P.p.gear.teeth==30&&P.p.wheel.starts==1&&A.feed>0);
assert(abs(A.omegaW/abs(A.omegaG)-30)<1e-8,'R30P:Ratio','This stock-phase scan requires the current 30:1 synchronous motion.');
if ~isfolder(outDir),mkdir(outDir);end
rb=P.p.calc.gear.baseDiameter_mm/2;ra=P.p.calc.gear.tipDiameter_mm/2;
smax=(ra^2-rb^2)/(2*rb);s=smax*linspace(0,1,o.profileSamples).^2;
z=linspace(-P.p.gear.faceWidth_mm/2,P.p.gear.faceWidth_mm/2,o.axialSamples);
[ss,zz]=ndgrid(s,z);stock=ones(size(ss))*P.cfg.stock_um/1000;
wheelPeriod=2*pi/abs(A.omegaW);phaseTimes=(0:o.phaseSamples-1)*wheelPeriod/o.phaseSamples;
originalStart=A.feedStart;start=originalStart;finish=max(-originalStart,originalStart+P.p.gear.faceWidth_mm)+o.tailMargin_mm;
assert(finish>start);
clearFound=false;
for attempt=0:ceil(o.maxExtension_mm/o.axialStep_mm)
 [hits,mn]=scan_at(start);
 if hits==0&&(isnan(mn)||mn>=o.clearanceLevel),clearFound=true;break;end
 start=start-o.axialStep_mm;
end
assert(clearFound,'R30P:Start','Could not find a sampled clear start; inspect the path/geometry.');
clearFound=false;
for attempt=0:ceil(o.maxExtension_mm/o.axialStep_mm)
 [hits,mn]=scan_at(finish);
 if hits==0&&(isnan(mn)||mn>=o.clearanceLevel),clearFound=true;break;end
 finish=finish+o.axialStep_mm;
end
assert(clearFound,'R30P:End','Could not find a sampled clear end; inspect the path/geometry.');
nominalTravel=abs(P.p.process.axialEnd_mm-P.p.process.axialStart_mm);
nominalTime=nominalTravel/A.feed;nominalEnd=originalStart+nominalTravel;
positions=unique([start:o.axialStep_mm:finish,finish,originalStart,nominalEnd,0]);
positions=positions(positions>=start&positions<=finish);
scan=zeros(numel(positions),6);potential=false([size(ss) 2]);clock=tic;
for j=1:numel(positions)
 [hits,mn,cov,candidates]=scan_at(positions(j));potential=potential|cov;
 scan(j,:)=[positions(j),(positions(j)-originalStart)/A.feed,hits,mn,candidates,double(hits==0&&(isnan(mn)||mn>=o.clearanceLevel))];
 if toc(clock)>o.progressEvery_s,fprintf('全程几何预检：%d/%d轴向位置，进给坐标%.4f mm。\n',j,numel(positions),positions(j));clock=tic;end
end
[initialHits,initialMin]=actual_pose(0);
[nominalEndHits,nominalEndMin]=actual_pose(nominalTime);
timeRange=([start finish]-originalStart)/A.feed;
duration=diff(timeRange);nt=ceil(duration/M.forceDt_s);nFine=ceil(duration/(M.forceDt_s/2));
plan=struct('kind','full_gear_nominal_stock_preflight_v1','fullForceComputed',false,...
 'sourceMat',sourceMat,'options',o,'parameters',P,'sourceGrid',[M.options.ns M.options.nz M.options.order],...
 'allToothIds',1:30,'faceCount',60,'feed_mmps',A.feed,...
 'feedCoordinateMeaning','a=feedStart+feed*t_model; world gear point Z=z-a',...
 'sampledClearFeedCoordinates_mm',[start finish],'modelTimeRange_s',timeRange,...
 'elapsedDuration_s',duration,'nominalInputTravel_mm',nominalTravel,'nominalInputDuration_s',nominalTime,...
 'nominalDurationEndFeedCoordinate_mm',nominalEnd,...
 'actualOriginalStartStockHits',initialHits,'actualOriginalStartMinLevel',initialMin,...
 'actualNominalEndStockHits',nominalEndHits,'actualNominalEndMinLevel',nominalEndMin,...
 'forceIntervalsAtCurrentDt',nt,'forceIntervalsAtHalfDt',nFine,...
 'exampleBytesIf12DoublesPerFacePerInterval',double(nt)*60*12*8,...
 'nominalWheelDiameter_mm',P.p.calc.wheel.outerDiameter_mm,'implicitWheelDiameter_mm',2*A.tipNew,...
 'forceCoefficientsCalibrated',logical(M.coefficients.calibrated),...
 'potentialSampleCoverageFractionBySide',[mean(reshape(potential(:,:,1),[],1)),mean(reshape(potential(:,:,2),[],1))],...
 'scope','One monotone, constant-feed stroke; above-base involute flanks only.',...
 'qualification','Sampled nominal-stock geometry coverage only. Not actual remaining-stock coverage, continuous contact exclusion, force convergence, or physical validation.');
save(fullfile(outDir,'full_trial_plan.mat'),'plan','scan','s','z','potential','-v7');
write_text(fullfile(outDir,'full_preflight_report.json'),jsonencode(plan));
write_csv(fullfile(outDir,'full_preflight_axial_scan.csv'),...
 {'feed_coordinate_mm','model_time_s','stock_hit_samples_all_phases','min_level_NOT_distance','candidate_face_poses','sampled_clear'},scan);
v=zeros(numel(ss)*2,4);k=0;
for side=1:2
 cv=potential(:,:,side);ix=k+(1:numel(ss));v(ix,:)=[ss(:),zz(:),ones(numel(ss),1)*side,double(cv(:))];k=k+numel(ss);
end
write_csv(fullfile(outDir,'potential_stock_coverage.csv'),{'profile_s_mm','local_gear_z_mm','side','ever_sampled_inside_wheel'},v);
fprintf('预检建议模型时间：%.6f..%.6f s；单向行程时长约%.3f s。\n',timeRange,duration);
fprintf('注意：这只是几何采样预检；完整受力、最终余量和连续接触边界尚未计算。\n');

 function [hits,minLevel,coverage,candidateCount]=scan_at(feedPosition)
  hits=0;minLevel=Inf;coverage=false([size(ss) 2]);candidateCount=0;
  for ph=1:numel(phaseTimes)
   t=phaseTimes(ph);PP=P;PP.installation.feedStart=feedPosition-A.feed*t;AA=PP.installation;
   for tooth=1:30
    for side=1:2
     if ~r30b_gate(PP,[tooth side],t,t,1),continue;end
     candidateCount=candidateCount+1;level=r30b_level(ss,zz,stock,t,PP.p,AA,tooth,side);
     mask=level<0;hits=hits+nnz(mask);minLevel=min(minLevel,min(level(:)));coverage(:,:,side)=coverage(:,:,side)|mask;
    end
   end
  end
  if isinf(minLevel),minLevel=NaN;end
 end
 function [hits,minLevel]=actual_pose(t)
  hits=0;minLevel=Inf;
  for tooth=1:30
   for side=1:2
    level=r30b_level(ss,zz,stock,t,P.p,A,tooth,side);hits=hits+nnz(level<0);minLevel=min(minLevel,min(level(:)));
   end
  end
 end
end
function write_csv(path,headers,data)
fid=fopen(path,'w');assert(fid>0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',strjoin(headers,','));
fprintf(fid,[repmat('%.17g,',1,size(data,2)-1) '%.17g\n'],data');
end
function write_text(path,value)
fid=fopen(path,'w','n','UTF-8');assert(fid>0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',value);
end
