function R=r30full_collect(outDir)
% Assemble completed tooth chunks by vector sum; never repeat one tooth.
D=load(fullfile(outDir,'full_run_plan.mat'),'meta','runToken');m=D.meta;token=D.runToken;
o=m.options;P=m.parameters;nt=m.intervalCount;nc=m.chunkCount;dt=o.dt_s;
stats=zeros(2*numel(o.toothIds),15);maxPartition=0;maxSubstep=0;
for g=1:numel(o.toothIds)
 leaf=fullfile(outDir,sprintf('tooth_%02d',o.toothIds(g)));
 D=load(fullfile(leaf,'checkpoint.mat'),'S');assert(strcmp(D.S.phase,'done')&&strcmp(D.S.runToken,token));clear D;
 D=load(fullfile(leaf,'diagnostics.mat'),'diagnostics');a=D.diagnostics;
 stats(2*g-1:2*g,:)=a.faceStatistics;maxPartition=max(maxPartition,a.maxPartitionError_mm3);maxSubstep=max(maxSubstep,a.maxSubstepError_mm3);
end
folder=fullfile(outDir,'global_chunks');if ~isfolder(folder),mkdir(folder);end
bin=max(1,round(o.overviewBin_s/dt));nb=ceil(nt/bin);A=zeros(nb,12);count=zeros(nb,1);
A(:,[2 7])=Inf;A(:,[3 5 8])=-Inf;
peakF=-Inf;peakIndex=1;peakN=0;peakTorque=0;maxTeeth=0;maxFaces=0;
volume=0;energy=0;positiveEnergy=0;positiveIntervals=0;maxPower=-Inf;minPower=Inf;
for j=1:nc
 first=(j-1)*o.chunkIntervals+1;last=min(nt,j*o.chunkIntervals);n=last-first+1;
 path=fullfile(folder,sprintf('chunk_%06d.mat',j));
 if isfile(path)
  D=load(path,'chunk');c=D.chunk;assert(strcmp(c.runToken,token)&&c.chunkIndex==j&&c.firstIndex==first&&c.lastIndex==last);
  G=zeros(n,16);G(double(c.rowIndex),:)=c.values;
 else
  G=zeros(n,16);
  for tooth=o.toothIds
   D=load(fullfile(outDir,sprintf('tooth_%02d',tooth),'chunks',sprintf('chunk_%06d.mat',j)),'chunk');c=D.chunk;
   assert(strcmp(c.runToken,token)&&c.tooth==tooth&&c.chunkIndex==j&&c.firstIndex==first&&c.lastIndex==last);
   ix=double(c.rowIndex);T=sum(c.values,3);
   G(ix,1:13)=G(ix,1:13)+T;
   G(ix,14)=G(ix,14)+(T(:,1)>0);
   G(ix,15)=G(ix,15)+sum(c.values(:,1,:)>0,3);
   G(ix,16)=G(ix,16)+(vecnorm(T(:,3:5),2,2)>o.loadThreshold_N);
  end
  keep=find(any(G~=0,2));chunk=struct('runToken',token,'chunkIndex',j,'firstIndex',first,'lastIndex',last,...
   'rowIndex',uint32(keep),'values',G(keep,:));r30b_atomic(path,'chunk',chunk);
 end
 f=vecnorm(G(:,3:5),2,2);nn=vecnorm(G(:,6:8),2,2);tc=G(:,11)/1000;
 [pf,ii]=max(f);if pf>peakF,peakF=pf;peakIndex=first+ii-1;end
 peakN=max(peakN,max(nn));peakTorque=max(peakTorque,max(abs(tc)));maxTeeth=max(maxTeeth,max(G(:,14)));maxFaces=max(maxFaces,max(G(:,15)));
 power=pair_power(G,P.installation);volume=volume+sum(G(:,1));energy=energy-sum(power)*dt;
 positiveEnergy=positiveEnergy+sum(max(power,0))*dt;positiveIntervals=positiveIntervals+nnz(power>1e-9);
 maxPower=max(maxPower,max(power));minPower=min(minPower,min(power));
 bins=floor(((first:last)'-1)/bin)+1;
 for b=unique(bins)'
  q=bins==b;count(b)=count(b)+nnz(q);
  A(b,1)=A(b,1)+sum(f(q));A(b,2)=min(A(b,2),min(f(q)));A(b,3)=max(A(b,3),max(f(q)));
  A(b,4)=A(b,4)+sum(nn(q));A(b,5)=max(A(b,5),max(nn(q)));A(b,6)=A(b,6)+sum(tc(q));
  A(b,7)=min(A(b,7),min(tc(q)));A(b,8)=max(A(b,8),max(tc(q)));
  A(b,9)=max(A(b,9),max(G(q,14)));A(b,10)=max(A(b,10),max(G(q,15)));
  A(b,11)=A(b,11)+sum(G(q,1));A(b,12)=A(b,12)+sum(G(q,12));
 end
 if mod(j,50)==0||j==nc,fprintf('合成全齿载荷：块%d/%d。\n',j,nc);end
end
assert(abs(volume-sum(stats(:,3)))<max(1e-8,abs(volume)*1e-10),'R30FULL:GlobalVolume','合成体积与各齿材料体积不闭合。');
A(:,[1 4 6 12])=bsxfun(@rdivide,A(:,[1 4 6 12]),count);
first=(0:nb-1)'*bin+1;last=min(nt,(1:nb)'*bin);overview=[dt*(first+last)/2 A];
write_csv(fullfile(outDir,'full_force_overview.csv'),...
 {'model_time_s','mean_force_N','min_force_N','max_force_N','mean_normal_resultant_N','max_normal_resultant_N',...
 'mean_Tc_Nm','min_Tc_Nm','max_Tc_Nm','max_removing_teeth','max_removing_faces','removed_volume_mm3','mean_scalar_Fn_sum_N'},overview);
write_csv(fullfile(outDir,'final_face_statistics.csv'),...
 {'tooth','side','removed_mm3','volume_closure_mm3','max_remaining_um','max_overcut_um','touched_area_fraction',...
 'within_tolerance_area_fraction','remaining_above_tolerance_area_fraction','overcut_below_tolerance_area_fraction',...
 'negative_z_edge_max_remaining_um','positive_z_edge_max_remaining_um','initial_stock_hit_points','end_nominal_stock_hit_points','geometry_evaluations'},stats);
half=max(2,round(.005/dt));ix=max(1,peakIndex-half):min(nt,peakIndex+half);
peakWindow=struct('indices',ix,'modelTime_s',ix(:)*dt,'globalData',r30full_read(outDir,ix),...
 'toothIds',o.toothIds,'perToothForce_N',zeros(numel(ix),numel(o.toothIds),3));
for g=1:numel(o.toothIds)
 T=r30full_read(outDir,ix,o.toothIds(g));peakWindow.perToothForce_N(:,g,:)=reshape(sum(T(:,3:5,:),3),numel(ix),1,3);
end
save(fullfile(outDir,'peak_window.mat'),'peakWindow','-v7');
spectrum=struct('computed',false);
if nt>=16
 n=min(nt,round((4*pi/abs(P.installation.omegaG))/dt));n=floor(n/2)*2;
 center=-P.installation.feedStart/P.installation.feed;
 first=max(1,min(nt-n+1,round(center/dt-n/2)+1));ix=first:first+n-1;Y=r30full_read(outDir,ix);
 signals=[vecnorm(Y(:,3:5),2,2),Y(:,11)/1000];signals=bsxfun(@minus,signals,mean(signals,1));w=.5-.5*cos(2*pi*(0:n-1)'/(n-1));
 Z=fft(bsxfun(@times,signals,w));amp=abs(Z(1:n/2+1,:))/sum(w);amp(2:end-1,:)=2*amp(2:end-1,:);freq=(0:n/2)'/(n*dt);
 spectrum=struct('computed',true,'modelTimeRange_s',[(first-1)*dt (first+n-1)*dt],...
  'sampleCount',n,'sampleRate_Hz',1/dt,'window','symmetric Hann; mean removed; coherent-gain amplitude correction',...
  'frequency_Hz',freq,'forceNormAmplitude_N',amp(:,1),'torqueAmplitude_Nm',amp(:,2),...
  'meaning','Finite central-traverse segment, not a stationarity proof or structural vibration measurement.');
 take=freq<=min(1000,.5/dt);write_csv(fullfile(outDir,'central_segment_fft_0_1000Hz.csv'),...
  {'frequency_Hz','force_norm_amplitude_N','Tc_amplitude_Nm'},[freq(take),amp(take,:)]);
end
save(fullfile(outDir,'central_segment_fft.mat'),'spectrum','-v7');
tail=min(nt,max(1,round(.01/dt)));headData=r30full_read(outDir,1:tail);tailData=r30full_read(outDir,nt-tail+1:nt);
R=struct('completed',true,'dataComputed',true,'fullGearComputed',isequal(sort(o.toothIds),1:30),...
 'requestedTrajectoryComputed',true,'machiningCoverageAccepted',false,'physicalAccuracyValidated',false,...
 'outputDir',outDir,'modelTimeRange_s',[0 nt*dt],'intervalCount',nt,'toothIds',o.toothIds,...
 'forceDt_s',dt,'geometrySubstep_s',dt/o.substeps,'stripWidth_mm',o.forceStripWidth_mm,...
 'rawPeakForce_N',peakF,'rawPeakForceModelTime_s',peakIndex*dt,'rawPeakNormalResultant_N',peakN,...
 'rawPeakAbsTorque_Nm',peakTorque,'maximumRemovingTeeth',maxTeeth,'maximumRemovingFaces',maxFaces,...
 'totalRemovedVolume_mm3',volume,'totalDissipatedPairEnergy_J',energy,'positivePairPowerIntervals',positiveIntervals,...
 'positivePairPowerEnergy_J',positiveEnergy,'pairPowerRange_W',[minPower maxPower],...
 'maxFaceVolumeClosure_mm3',max(stats(:,4)),'maxStripPartitionError_mm3',maxPartition,'maxSubstepError_mm3',maxSubstep,...
 'maxRemainingStock_um',max(stats(:,5)),'maxOvercut_um',max(stats(:,6)),...
 'allEndpointStockSamplesClear',all(stats(:,13)==0&stats(:,14)==0),...
 'headTailIntervalCount',tail,'headTailRemovedVolume_mm3',[sum(headData(:,1)) sum(tailData(:,1))],...
 'targetTolerance_mm',o.targetTolerance_mm,'edgeBand_mm',o.edgeBand_mm,...
 'spectrumComputed',spectrum.computed,'overviewMeaning','Mean and extrema per bin; raw common-time loads are retained in global_chunks.',...
 'materialMeaning','Independent signed normal heights; thin-layer ds dz dh volume; no clipping to target.',...
 'scope','Current ideal geometry, one constant axial-feed stroke, above-base involute flanks, active removal only; force coefficients uncalibrated.');
save(fullfile(outDir,'full_run_result.mat'),'R','overview','stats','-v7');
fid=fopen(fullfile(outDir,'full_run_report.json'),'w','n','UTF-8');assert(fid>0);fprintf(fid,'%s',jsonencode(R));fclose(fid);
fprintf('全程载荷计算完成。请检查最终余量与过切，不能用completed字段判定磨削合格。\n');
end
function p=pair_power(G,A)
F=G(:,3:5);mm=G(:,9:11);tw=-mm+cross(repmat([A.center 0 0],size(G,1),1),F,2);
p=(A.omegaG*mm(:,3)-A.feed*F(:,3)+tw*(A.omegaW*A.axisW(:)))/1000;
end
function write_csv(path,header,X)
fid=fopen(path,'w');assert(fid>0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',strjoin(header,','));fprintf(fid,[repmat('%.17g,',1,size(X,2)-1) '%.17g\n'],X');
end
