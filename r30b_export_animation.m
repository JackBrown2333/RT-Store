function zipFiles = r30b_export_animation(resultDir,userOptions)
% R30B_EXPORT_ANIMATION 从现有R30B结果导出小型动画数据；不重算磨削。
% 放在run_r30_dual.m同级，在MATLAB执行：
%   zipFiles = r30b_export_animation;
% 也可指定结果目录：
%   zipFiles = r30b_export_animation('r30_dual_recomputed');
%
% 主要原则：
%  1. 所有源帧逐一回放；抽帧之间的材料变化累计后写入保留帧。
%  2. 仅显示表面选取部分原材料点，保留点高度使用原double值。
%  3. 接触轮廓取原活动支撑矩形的并集，不以凸包替代；保留孔洞和分岛。
%  4. 完整力曲线、保留帧的全部条带力和尺寸锚点均直接复制。
%  5. 不把局部时间加密结果拼接进原曲线；源步长写入元数据。
%  6. 不改原结果、检查点及frames；所有输出另存新目录。
%
% 可选参数（显示采样数不改变原计算网格）：
%   o=struct('frameStride',5,'profileSamples',96,'axialSamples',128,...
%     'framesPerChunk',24,'targetPartMiB',15,'keepCriticalFrames',true);
%   zipFiles=r30b_export_animation('r30_dual_recomputed',o);
%
% 依赖：MATLAB内置load/save/jsonencode/zip，无额外工具箱。

if nargin<1 || isempty(resultDir)
    resultDir=fullfile(fileparts(mfilename('fullpath')),'r30_dual_recomputed');
end
if nargin<2,userOptions=struct();end
o=export_options(userOptions);
[ok,a]=fileattrib(resultDir);assert(ok && a.directory,'R30B:ExportFolder','找不到结果目录：%s',resultDir);root=a.Name;
resultFile=fullfile(root,'r30b_result.mat');assert(isfile(resultFile),'R30B:ExportInput','缺少r30b_result.mat');
D=load(resultFile,'R30B');assert(isfield(D,'R30B'),'R30B:ExportInput','MAT文件没有R30B变量');R=D.R30B;
assert(strcmp(R.engineVersion,'r30_dual_from_zero_v2') && R.historyFromZero && ~R.oldResultsUsed,...
    'R30B:ExportVersion','要求新版从零历史结果，不接受旧版结果。');
nf=size(R.faces,1);nt=numel(R.poses);nsrc=numel(R.frameFiles);
assert(nsrc>0 && numel(R.frameIndices)==nsrc,'R30B:ExportFrames','结果中没有完整的帧目录索引。');
sourceIndices=double(R.frameIndices(:)');
assert(all(diff(sourceIndices)>0) && sourceIndices(1)>=1 && sourceIndices(end)==nt,...
    'R30B:ExportFrames','源帧序号有重复、顺序错误或缺少末帧。');
for j=1:nsrc
    assert(isfile(fullfile(root,R.frameFiles{j})),'R30B:ExportMissing',...
        '缺少源帧：%s。请使用原始完整目录，不能使用诊断抽帧包。',R.frameFiles{j});
end
for name={R.initialStateFile,R.finalStateFile}
    assert(isfile(fullfile(root,name{1})),'R30B:ExportMissing','缺少%s',name{1});
end
checks=check_summaries(R);
D=load(fullfile(root,R.initialStateFile),'initialState');initial=D.initialState;
assert(abs(initial.t_s-R.window_s(1))<1e-10,'R30B:ExportTime','初始状态时间不符。');
height=cell(nf,1);weights=height;selectedIds=height;previousDisplay=height;
displayFaces=cell(nf,1);supports=cell(nf,2);
for f=1:nf
    q=initial.faces{f};assert(q.tooth_id==R.faces(f,1) && q.flank_side==R.faces(f,2));
    h=q.height_mm;assert(all(isfinite(h(:))),'R30B:ExportState','材料高度存在非有限值。');
    [nr,nc]=size(h);se=q.support_s_mm(:);ze=q.support_z_mm(:);
    assert(numel(se)==nr+1 && numel(ze)==nc+1 && all(diff(se)>0) && all(diff(ze)>0));
    rows=unique(round(linspace(1,nr,min(nr,o.profileSamples))))';
    cols=unique(round(linspace(1,nc,min(nc,o.axialSamples))))';
    idx=bsxfun(@plus,rows,(cols'-1)*nr);idx=idx(:);
    height{f}=h;weights{f}=diff(se)*diff(ze)';selectedIds{f}=idx;
    previousDisplay{f}=reshape(h(idx),numel(rows),numel(cols));supports{f,1}=se;supports{f,2}=ze;
    displayFaces{f}=struct('tooth_id',q.tooth_id,'flank_side',q.flank_side,...
        'sourceSize',[nr nc],'sourceRowIndices',uint32(rows),'sourceColumnIndices',uint32(cols),...
        'displaySize',[numel(rows) numel(cols)],'sample_s_mm',q.sample_s_mm(rows),...
        'sample_z_mm',q.sample_z_mm(cols),'height_mm',previousDisplay{f},...
        'stock_nominal_mm',q.stock_nominal_mm);
end
clear initial D;
keep=false(nsrc,1);keep(unique([1:o.frameStride:nsrc nsrc]))=true;
if o.keepCriticalFrames
    wanted=critical_indices(R);
    for k=wanted(:)'
        [~,j]=min(abs(sourceIndices-k));keep(max(1,j-1):min(nsrc,j+1))=true;
    end
end
selectedSource=sourceIndices(keep);
if isempty(o.outputDir)
    base=fullfile(root,['animation_export_' datestr(now,'yyyymmdd_HHMMSS')]);outDir=base;c=0;
    while isfolder(outDir),c=c+1;outDir=[base '_' num2str(c)];end
else
    outDir=o.outputDir;assert(~isfolder(outDir),'R30B:ExportOutput','输出目录已存在，请指定新目录。');
end
mkdir(outDir);[ok,a]=fileattrib(outDir);assert(ok);outDir=a.Name;
header=struct('schemaVersion','r30_animation_compact_v1','bundleId',outDir_name(outDir),...
    'engineVersion',R.engineVersion,'faces',R.faces,'window_s',R.window_s,...
    'parameters',R.parameters,'sourceResolution',struct('ns',R.options.ns,'nz',R.options.nz,...
    'order',R.options.order,'forceOutputDt_s',R.options.dt_s,'geometrySubsteps',R.options.substeps,...
    'historyDt_s',R.options.historyDt_s,'forceStripWidth_mm',R.options.forceStripWidth_mm,...
    'savedFrameStride',R.options.frameStride),...
    'initialTime_s',R.window_s(1),'displayFaces',{displayFaces},'forceCoefficients',R.forceCoefficients,...
    'momentOrigin_mm',R.momentOrigin_mm,'referenceOutward',R.options.referenceOutward,...
    'historyFromZero',R.historyFromZero,'wholePassEnclosed',R.wholePassEnclosed,...
    'geometryRepresentation',R.geometryRepresentation,'forceTimeMeaning',R.forceTimeMeaning,...
    'coefficientMeaning',R.coefficientMeaning,'scope',R.scope,...
    'materialIndexing','MATLAB 1-based column-major; exported indices refer to DISPLAY grid',...
    'displayHeightMeaning','exact double heights at selected original material points; interpolation between display points is visualization only',...
    'materialHeightMeaning','signed normal height above target; no zero floor',...
    'lgMeaning','area divided by transverse width is an equivalent mean chord; actual chord segments are stored separately',...
    'boundaryMeaning','directed union boundary of original active positive-weight support rectangles in developed s,z; not exact continuous wheel intersection',...
    'boundaryWorldHeightMeaning','world overlay may interpolate the display surface; saved world dimension anchors and force locations are retained directly',...
    'forceDataChanged',false,'localRefinementSpliced',false,'fullGearComputed',false,...
    'timeConvergenceOfWholeWindowValidated',false,'exportOptions',o);
curves=struct('timeEdges_s',R.timeEdges_s,'faces',R.faces,'combined',R.poses,'perFace',{R.facePoses},...
    'forceOutputDt_s',R.options.dt_s,'forceTimeMeaning',R.forceTimeMeaning,'forceDataChanged',false);
save(fullfile(outDir,'animation_header.mat'),'header','-v7');
save(fullfile(outDir,'animation_curves.mat'),'curves','-v7');
write_text(fullfile(outDir,'README_zh.txt'),bundle_readme());
commonBytes=0;
for name={'animation_header.mat','animation_curves.mat','README_zh.txt'}
    info=dir(fullfile(outDir,name{1}));commonBytes=commonBytes+info.bytes;
end
chunkBudget=min(o.targetPartMiB*1024^2*.45,o.targetPartMiB*1024^2*.90-commonBytes-16384);
assert(chunkBudget>1024,'R30B:ExportPartSize','分包上限不足以容纳公共数据，请增加targetPartMiB。');
checks.sourceFramesRead=nsrc;checks.exportedFrames=nnz(keep);checks.sourceIntervals=nt;
checks.maxReplayVolumeDifference_mm3=0;checks.maxBoundaryAreaDifference_mm2=0;
checks.maxDisplayReplayError_mm=0;checks.maxFullReplayError_mm=0;
checks.forceDataChanged=false;checks.localRefinementSpliced=false;
lastSourceTime=R.window_s(1);lastIndex=0;lastDisplayTime=lastSourceTime;
batch={};checkpoint=previousDisplay;checkpointTime=lastDisplayTime;chunkNumber=1;chunkFiles={};
clock=tic;progress=tic;
for j=1:nsrc
    D=load(fullfile(root,R.frameFiles{j}),'frame');frame=D.frame;k=double(frame.index);
    assert(k==sourceIndices(j) && abs(frame.t_s-R.poses(k).t_s)<1e-10,...
        'R30B:ExportTime','源帧序号或时间不一致。');
    compactFaces=cell(nf,1);
    for f=1:nf
        x=frame.faces{f};assert(x.tooth_id==R.faces(f,1) && x.flank_side==R.faces(f,2));
        u=x.stateUpdates;ids=double(u.index(:));values=u.height_mm(:);h=height{f};
        assert(abs(u.since_t_s-lastSourceTime)<1e-10,'R30B:ExportReplay','材料增量的时间链断开。');
        assert(numel(ids)==numel(values) && numel(unique(ids))==numel(ids) && ...
            all(isfinite(ids)) && all(ids>=1 & ids<=numel(h) & mod(ids,1)==0) && all(isfinite(values)),...
            'R30B:ExportReplay','材料更新索引或高度无效。');
        assert(all(values<=h(ids)+1e-13),'R30B:ExportReplay','更新恢复了已去除的材料。');
        replayVolume=sum(weights{f}(ids).*(h(ids)-values));
        targetVolume=sum([R.facePoses{f}(lastIndex+1:k).dV_mm3]);
        checks.maxReplayVolumeDifference_mm3=max(checks.maxReplayVolumeDifference_mm3,abs(replayVolume-targetVolume));
        assert(abs(replayVolume-targetVolume)<1e-9,'R30B:ExportVolume','源帧回放体积与汇总不一致。');
        h(ids)=values;height{f}=h;
        if keep(j)
            sample=reshape(h(selectedIds{f}),size(previousDisplay{f}));old=previousDisplay{f};
            changed=find(sample~=old);assert(all(sample(changed)<=old(changed)+1e-13));
            updates=struct('since_t_s',lastDisplayTime,'index',uint32(changed),'height_mm',sample(changed));
            replay=old;replay(double(updates.index))=updates.height_mm;
            checks.maxDisplayReplayError_mm=max(checks.maxDisplayReplayError_mm,max(abs(replay(:)-sample(:))));
            [segments,footprintArea]=contact_outline(x.removal.index,supports{f,1},supports{f,2});
            error=abs(footprintArea-R.facePoses{f}(k).area_mm2);
            checks.maxBoundaryAreaDifference_mm2=max(checks.maxBoundaryAreaDifference_mm2,error);
            assert(error<1e-8,'R30B:ExportBoundary','导出轮廓面积与原活动支撑面积不一致。');
            compactFaces{f}=struct('tooth_id',x.tooth_id,'flank_side',x.flank_side,'passId',x.passId,...
                'stateUpdates',updates,'contactInterval_s',[frame.t_start_s frame.t_s],...
                'boundarySegments_sz_mm',segments,'contactArea_mm2',footprintArea,...
                'strips',compact_strips(x.strips),'summary',x.summary);
            previousDisplay{f}=sample;
        end
    end
    if keep(j)
        compact=struct('sourceIntervalIndex',k,'t_s',frame.t_s,'t_start_s',frame.t_start_s,...
            'sourceDt_s',frame.dt_s,'faces',{compactFaces},'combined',frame.combined);
        batch{end+1}=compact;lastDisplayTime=frame.t_s; %#ok<AGROW>
        if numel(batch)>=o.framesPerChunk || j==nsrc
            [files,chunkNumber]=write_chunks(outDir,batch,checkpoint,checkpointTime,chunkNumber,chunkBudget);
            chunkFiles=[chunkFiles files];batch={};checkpoint=previousDisplay;checkpointTime=lastDisplayTime; %#ok<AGROW>
        end
    end
    lastSourceTime=frame.t_s;lastIndex=k;
    if toc(progress)>5
        fprintf('动画导出：已读取%d/%d源帧，保留%d帧；用时%.1f s。\n',j,nsrc,sum(keep(1:j)),toc(clock));progress=tic;
    end
end
D=load(fullfile(root,R.finalStateFile),'finalState');
assert(abs(D.finalState.t_s-R.window_s(2))<1e-10,'R30B:ExportFinal','最终状态时间不符。');
for f=1:nf
    final=D.finalState.faces{f};assert(final.tooth_id==R.faces(f,1) && final.flank_side==R.faces(f,2));
    assert(isequal(size(height{f}),size(final.height_mm)));
    checks.maxFullReplayError_mm=max(checks.maxFullReplayError_mm,max(abs(height{f}(:)-final.height_mm(:))));
    sample=final.height_mm(selectedIds{f});
    checks.maxDisplayReplayError_mm=max(checks.maxDisplayReplayError_mm,max(abs(previousDisplay{f}(:)-sample(:))));
end
assert(checks.maxFullReplayError_mm<1e-12 && checks.maxDisplayReplayError_mm<1e-12,...
    'R30B:ExportFinal','完整材料回放或显示回放与最终状态不符。');
checks.status='PASS';checks.elapsed_s=toc(clock);checks.sourceForceOutputDt_s=R.options.dt_s;
checks.displaySamplingOnly=true;checks.wholeWindowTimeConvergenceValidated=false;
manifest=struct('schemaVersion','r30_animation_compact_v1','bundleId',header.bundleId,...
    'faces',R.faces,'window_s',R.window_s,'sourceIntervalCount',nt,'sourceSavedFrameCount',nsrc,...
    'selectedSourceIntervalIndices',selectedSource,'chunks',{chunkFiles},...
    'forceOutputDt_s',R.options.dt_s,'forceDataChanged',false,'parts',{{}});
write_json(fullfile(outDir,'export_report.json'),checks);
write_json(fullfile(outDir,'animation_manifest.json'),manifest);
shared={'animation_header.mat','animation_curves.mat','animation_manifest.json','export_report.json','README_zh.txt'};
groups=pack_groups(outDir,chunkFiles,shared,o.targetPartMiB);
names=cell(numel(groups),1);
for p=1:numel(groups),names{p}=sprintf('r30_animation_%02d.zip',p);end
manifest.parts=names;manifest.partChunks=groups;write_json(fullfile(outDir,'animation_manifest.json'),manifest);
zipFiles=cell(numel(groups),1);
for p=1:numel(groups)
    zipFiles{p}=fullfile(outDir,names{p});zip(zipFiles{p},[shared groups{p}],outDir);
    info=dir(zipFiles{p});assert(info.bytes<=o.targetPartMiB*1024^2,...
        'R30B:ExportPartSize','ZIP超出设定大小，请用更大的targetPartMiB导出到新目录。');
end
fprintf('\n导出核验：PASS。完整材料回放误差 %.3g mm；显示材料回放误差 %.3g mm。\n',...
    checks.maxFullReplayError_mm,checks.maxDisplayReplayError_mm);
fprintf('原力计算区间 %.6g us；此次仅导出，不改变原力数据。\n',R.options.dt_s*1e6);
fprintf('请上传下面所有ZIP，不需要上传原frames：\n');
for p=1:numel(zipFiles),info=dir(zipFiles{p});fprintf('%s  (%.2f MiB)\n',zipFiles{p},info.bytes/1024^2);end
end

function o=export_options(user)
o=struct('frameStride',5,'profileSamples',96,'axialSamples',128,'framesPerChunk',24,...
    'targetPartMiB',15,'keepCriticalFrames',true,'outputDir','');
names=fieldnames(user);for j=1:numel(names)
    assert(isfield(o,names{j}),'R30B:ExportOptions','未知选项：%s',names{j});o.(names{j})=user.(names{j});
end
n=[o.frameStride o.profileSamples o.axialSamples o.framesPerChunk];
assert(all(isfinite(n)) && all(n>=1) && all(mod(n,1)==0),'R30B:ExportOptions','采样数和步幅必须是正整数。');
assert(isscalar(o.targetPartMiB) && isfinite(o.targetPartMiB) && o.targetPartMiB>=.05,...
    'R30B:ExportOptions','targetPartMiB过小或无效。');
assert(isscalar(o.keepCriticalFrames) && ismember(o.keepCriticalFrames,[false true]));
end

function checks=check_summaries(R)
checks=struct('maxSummaryForceDifference_N',0,'maxSummaryMomentDifference_Nmm',0,'maxSummaryVolumeDifference_mm3',0);
for k=1:numel(R.poses)
    F=[0 0 0];N=F;M=F;V=0;
    for f=1:size(R.faces,1)
        q=R.facePoses{f}(k);assert(abs(q.t_s-R.poses(k).t_s)<1e-10);
        F=F+q.force_N;N=N+q.normalVector_N;M=M+q.moment_Nmm;V=V+q.dV_mm3;
    end
    q=R.poses(k);
    checks.maxSummaryForceDifference_N=max(checks.maxSummaryForceDifference_N,max(abs([F-q.force_N N-q.normalVector_N norm(N)-q.Fn_pure_resultant_N])));
    checks.maxSummaryMomentDifference_Nmm=max(checks.maxSummaryMomentDifference_Nmm,max(abs(M-q.moment_Nmm)));
    checks.maxSummaryVolumeDifference_mm3=max(checks.maxSummaryVolumeDifference_mm3,abs(V-q.dV_mm3));
end
assert(checks.maxSummaryForceDifference_N<1e-7 && checks.maxSummaryMomentDifference_Nmm<1e-5 &&...
    checks.maxSummaryVolumeDifference_mm3<1e-9,'R30B:ExportSummary','各侧汇总与合力不一致。');
end

function wanted=critical_indices(R)
nt=numel(R.poses);wanted=[1 nt];both=true(nt,1);
for f=1:size(R.faces,1)
    q=R.facePoses{f};v=[q.dV_mm3];fn=[q.Fn_pure_resultant_N];dh=[q.dh_step_peak_mm];
    active=find(v>0);both=both & (v(:)>0);if isempty(active),continue;end
    [~,p]=max(fn);[~,r]=max(v);[~,d]=max(dh);
    near=active(1):min(nt,active(1)+ceil(.003/R.options.dt_s));[~,e]=max(fn(near));
    wanted=[wanted active(1) active(end) p r d near(e)]; %#ok<AGROW>
end
overlap=find(both);if ~isempty(overlap),[~,q]=min([R.poses(overlap).Fn_pure_resultant_N]);wanted(end+1)=overlap(q);end
wanted=unique(wanted);
end

function [segments,area]=contact_outline(ids,se,ze)
% 四邻接支撑矩形的有向边界；只合并同向共线连续边，不填平孔洞。
nr=numel(se)-1;nc=numel(ze)-1;mask=false(nr,nc);ids=double(ids(:));
assert(numel(unique(ids))==numel(ids) && all(ids>=1 & ids<=numel(mask) & mod(ids,1)==0));
mask(ids)=true;pieces={};n=0;
D=diff([zeros(1,nc);double(mask);zeros(1,nc)],1,1);
for i=find(any(D~=0,2))'
    line=D(i,:);[a,b]=runs(line);if isempty(a),continue;end
    B=[repmat(se(i),numel(a),1) ze(a) repmat(se(i),numel(a),1) ze(b+1)];
    rev=line(a)>0;B(rev,:)=B(rev,[3 4 1 2]);n=n+1;pieces{n}=B; %#ok<AGROW>
end
D=diff([zeros(nr,1) double(mask) zeros(nr,1)],1,2);
for j=find(any(D~=0,1))
    line=D(:,j)';[a,b]=runs(line);if isempty(a),continue;end
    B=[se(a) repmat(ze(j),numel(a),1) se(b+1) repmat(ze(j),numel(a),1)];
    rev=line(a)<0;B(rev,:)=B(rev,[3 4 1 2]);n=n+1;pieces{n}=B; %#ok<AGROW>
end
if isempty(pieces),segments=zeros(0,4);area=0;
else
    segments=vertcat(pieces{:});area=sum(segments(:,1).*segments(:,4)-segments(:,3).*segments(:,2))/2;
end
end

function [a,b]=runs(line)
a=find(line~=0 & [true line(2:end)~=line(1:end-1)]);
b=find(line~=0 & [line(1:end-1)~=line(2:end) true]);
a=a(:);b=b(:);assert(numel(a)==numel(b));
end

function rows=compact_strips(rows)
% 原条带力、位置、真实截线与尺寸锚点保留；仅去除重建轮廓不需要的大多边形。
for k=1:numel(rows)
    assert(isfield(rows(k),'geometry'),'R30B:ExportGeometry','源帧缺少尺寸与截线几何。');
    g=rows(k).geometry;drop=intersect(fieldnames(g),{'clippedSupportPolygons','vertexCounts'});
    if ~isempty(drop),g=rmfield(g,drop);end;rows(k).geometry=g;
end
end

function [names,next]=write_chunks(outDir,frames,startState,startTime,index,budget)
% MAT已压缩后再判断大小；过大时二分，确保中间块有自己的材料起点。
name=sprintf('animation_chunk_%04d.mat',index);path=fullfile(outDir,name);
chunk=struct('schemaVersion','r30_animation_chunk_v1','chunkIndex',index,'startTime_s',startTime,...
    'startState',{startState},'frames',{frames},'firstSourceIntervalIndex',frames{1}.sourceIntervalIndex,...
    'lastSourceIntervalIndex',frames{end}.sourceIntervalIndex);
save(path,'chunk','-v7');info=dir(path);
if info.bytes>budget && numel(frames)>1
    delete(path);mid=floor(numel(frames)/2);state=startState;
    for j=1:mid
        for f=1:numel(state)
            u=frames{j}.faces{f}.stateUpdates;h=state{f};h(double(u.index))=u.height_mm;state{f}=h;
        end
    end
    [left,next]=write_chunks(outDir,frames(1:mid),startState,startTime,index,budget);
    [right,next]=write_chunks(outDir,frames(mid+1:end),state,frames{mid}.t_s,next,budget);
    names=[left right];
else,names={name};next=index+1;
end
end

function groups=pack_groups(outDir,chunks,shared,targetMiB)
target=targetMiB*1024^2;common=0;
for k=1:numel(shared),a=dir(fullfile(outDir,shared{k}));common=common+a.bytes;end
% 给最终清单增加的文件名以及ZIP头保留空间；实际ZIP写完后再次核对。
limit=target*.92-16384;assert(common<limit,'R30B:ExportPartSize',...
    '公共元数据超过分包上限，请增加targetPartMiB。');
groups={};current={};bytes=common;
for k=1:numel(chunks)
    a=dir(fullfile(outDir,chunks{k}));assert(common+a.bytes<limit,...
        'R30B:ExportPartSize','单帧块超过分包上限，请增加targetPartMiB。');
    if bytes+a.bytes>limit && ~isempty(current)
        groups{end+1}=current;current={};bytes=common; %#ok<AGROW>
    end
    current{end+1}=chunks{k};bytes=bytes+a.bytes; %#ok<AGROW>
end
if ~isempty(current),groups{end+1}=current;end
end

function name=outDir_name(path)
[~,name]=fileparts(path);
end

function write_json(path,value)
write_text(path,jsonencode(value));
end

function write_text(path,value)
fid=fopen(path,'w','n','UTF-8');assert(fid>0,'R30B:ExportOutput','无法写入%s',path);
guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',value);
end

function text=bundle_readme()
text=sprintf(['这是现有R30B结果的动画精简数据，不是新的磨削计算结果。\n'...
    '所有源帧逐一回放；抽帧之间的材料变化已累计进保留帧。\n'...
    '显示表面仅选取部分原材料点；这些点的高度保留原double值。\n'...
    '接触轮廓是原活动正权支撑矩形在s,z坐标中的并集边界，保留孔洞及分岛。\n'...
    '它不是精确连续砂轮求交线；世界坐标覆盖图可在显示表面上插值。\n'...
    '原条带力、尺寸世界坐标锚点、真实截线及完整汇总曲线均直接保留。\n'...
    '力数据未修改，步长见animation_header.mat的sourceResolution。\n'...
    '局部时间加密结果没有拼接进整段曲线。\n'...
    '每个ZIP都是独立压缩包，不是二进制分卷；请上传所有ZIP。\n'...
    '重复公共文件应一致；所有animation_chunk需要按chunkIndex顺序读取。\n'...
    '每个chunk包含自己的显示材料起点，可独立检查局部回放。\n'...
    '原始frames、initial_state、final_state继续保留在本机。\n']);
end
