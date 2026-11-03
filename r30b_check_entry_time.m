function report = r30b_check_entry_time(resultDir, window_s, factors)
% R30B_CHECK_ENTRY_TIME 读取原结果，在短入切区间做完整空间网格的时间加密。
% 不从t=0重算：从initial_state和前序frames恢复同一个局部起始材料状态。
% 不改动原结果、原frames、检查点和求解器；全部输出另存新目录。
%
% 用法（本文件放在run_r30_dual.m同级）：
%   addpath('matlab_recompute');
%   report = r30b_check_entry_time;
% 或显式指定：
%   report = r30b_check_entry_time('r30_dual_recomputed', ...
%              [14.9767 14.9790], [1 2 4 8]);
%
% 因子1复现原100us输出/50us子步；因子2、4、8分别加密输出和子步。
% 比较时先对同一个原100us区间内的力向量作时间平均，再取模。
% 本检查固定原来的空间网格和此前材料历史，不验证历史步长、空间
% 网格、条带宽度收敛，也不验证砂轮实测形状和力系数。
% 本文件在交付环境未运行MATLAB/Octave；相同局部算法已用Python
% 对上传帧做过独立复算。请以本机运行和baseline复现检查结果为准。

if nargin < 1 || isempty(resultDir)
    resultDir = fullfile(fileparts(mfilename('fullpath')), 'r30_dual_recomputed');
end
if nargin < 2, window_s = []; end
if nargin < 3 || isempty(factors), factors = [1 2 4 8]; end
assert(isfolder(resultDir), 'R30B:EntryFolder', '找不到原结果目录：%s', resultDir);
assert(exist('r30b_init','file') == 2, 'R30B:EntryPath', ...
    '先执行 addpath(''matlab_recompute'')。');
assert(isnumeric(factors) && isvector(factors) && factors(1) == 1 && ...
    all(isfinite(factors)) && all(factors >= 1) && all(mod(factors,1) == 0) && ...
    all(diff(factors) > 0), 'R30B:EntryFactors', 'factors必须从1开始、为严格递增整数。');
factors = factors(:)';
D = load(fullfile(resultDir,'r30b_result.mat'), 'R30B'); R = D.R30B;
assert(strcmp(R.engineVersion,'r30_dual_from_zero_v2') && R.historyFromZero && ...
    ~R.oldResultsUsed, 'R30B:EntryVersion', '只支持新版从零历史结果。');
assert(R.options.frameStride == 1 && ...
    isequal(R.frameIndices(:)',1:numel(R.poses)), 'R30B:EntryFrames', ...
    '此脚本要求原结果逐步保存完整frames；诊断抽帧包不能用于本检查。');
dt = R.options.dt_s; nf = size(R.faces,1); nt = numel(R.poses);
if isempty(window_s)
    f = find(R.faces(:,1)==17 & R.faces(:,2)==1,1);
    assert(~isempty(f),'R30B:EntryFace','未找到17齿1侧。');
    first = find([R.facePoses{f}.dV_mm3] > 0,1);
    assert(~isempty(first),'R30B:EntryNoCut','此结果没有17齿1侧的材料去除。');
    k0 = max(1,first-1); k1 = min(nt,first+21);
else
    assert(isnumeric(window_s) && numel(window_s)==2 && ...
        all(isfinite(window_s)) && window_s(2)>window_s(1), ...
        'R30B:EntryWindow','window_s应为递增的两个有限数。');
    q = (window_s(:)'-R.window_s(1))/dt;
    assert(max(abs(q-round(q)))<1e-7, 'R30B:EntryWindow', ...
        '局部窗口两端必须与原时间边界对齐。');
    k0 = round(q(1))+1; k1 = round(q(2));
end
assert(k0>=1 && k1<=nt && k1>=k0,'R30B:EntryWindow','局部窗口超出原结果。');
window_s = R.timeEdges_s([k0 k1+1]);

% 逐帧覆盖绝对高度；不把缺失帧当作零去除。此处不读取大检查点。
D = load(fullfile(resultDir,R.initialStateFile),'initialState');
initial = D.initialState; H = cell(nf,1);
for f=1:nf
    assert(initial.faces{f}.tooth_id==R.faces(f,1) && ...
        initial.faces{f}.flank_side==R.faces(f,2));
    H{f}=initial.faces{f}.height_mm;
end
lastTime = initial.t_s;
for k=1:k0-1
    frame = load_frame(resultDir,R,k);
    for f=1:nf
        x=frame.faces{f}; u=x.stateUpdates;
        assert(abs(u.since_t_s-lastTime)<1e-10,'R30B:EntryReplay','前序帧时间不连续。');
        h=H{f}; assert(all(u.height_mm<=h(u.index)+1e-12), ...
            'R30B:EntryReplay','前序材料更新使已去除材料恢复。');
        h(u.index)=u.height_mm; H{f}=h;
    end
    lastTime=frame.t_s;
end
assert(abs(lastTime-window_s(1))<1e-10,'R30B:EntryReplay','局部起点未正确恢复。');
Q0=cell(nf,1); o=R.options; o.saveFrames=false; o.plot=false; o.resume=false;
for f=1:nf
    Q=r30b_init(R.parameters,o,R.faces(f,:)); h=H{f};
    assert(isequal(size(h),[Q.ns*Q.n Q.nz*Q.n]),'R30B:EntryGrid','材料高度尺寸不符。');
    % r30b_heights的精确逆变换，保留每个Gauss材料点，不作插值。
    Q.h=permute(reshape(h,Q.n,Q.ns,Q.n,Q.nz),[2 4 1 3]);
    layoutDifference=r30b_heights(Q)-h;
    assert(max(abs(layoutDifference(:)))<1e-14,'R30B:EntryGrid','材料布局逆变换失败。');
    % passStart仅供条带函数的局部余量字段；本检查不导出ae_pass。
    % 力计算始终使用实际before-after体积，不依赖这个局部标记。
    Q.passStart=Q.h; Q0{f}=Q;
end
clear H D initial;

outDir=fullfile(resultDir,['entry_time_check_' datestr(now,'yyyymmdd_HHMMSS')]);
assert(~isfolder(outDir),'R30B:EntryOutput','输出目录已存在，请稍后重试。');
mkdir(outDir); runs=cell(numel(factors),1);
baseline=struct('maxNormalVectorDifference_N',0,'maxForceVectorDifference_N',0, ...
    'maxVolumeDifference_mm3',0,'maxSavedActiveHeightDifference_mm',0, ...
    'changedPointSetsAgree',true);
ncoarse=k1-k0+1;
for r=1:numel(factors)
    factor=factors(r); odt=dt/factor; Qs=Q0;
    records=cell(ncoarse*factor,1);
    fprintf('局部完整网格：输出步长 %.6g us，子步 %.6g us；无需重建前序历史。\n', ...
        odt*1e6,odt/o.substeps*1e6);
    for j=1:ncoarse*factor
        ta=window_s(1)+(j-1)*odt; tb=window_s(1)+j*odt;
        faceData=cell(nf,1); totalF=[0 0 0]; totalN=totalF; totalM=totalF; totalV=0;
        if factor==1, original=load_frame(resultDir,R,k0+j-1); end
        for f=1:nf
            Q=Qs{f}; before=r30b_heights(Q); microV=0;
            for sub=1:o.substeps
                ts=ta+(tb-ta)*sub/o.substeps;
                % 完整空间网格仍由原空间筛选处理；此局部检查不使用时间gate。
                [Q,v]=r30b_step(Q,ts,o); microV=microV+v;
            end
            after=r30b_heights(Q); dh=before-after;
            volume=sum(Q.weightsFine(:).*dh(:)); active=dh>0;
            assert(abs(volume-microV)<1e-9,'R30B:EntryVolume','子步与区间体积不一致。');
            localPass=r30b_heights(Q0{f});
            rows=r30b_strips(Q,before,after,localPass,tb,odt,o,false);
            rows=r30b_force(rows,R.forceCoefficients);
            s=r30b_sum(rows,o.referenceOutward,tb,odt,volume,sum(Q.weightsFine(active)));
            faceData{f}=selected_summary(s);
            totalF=totalF+s.force_N; totalN=totalN+s.normalVector_N;
            totalM=totalM+s.moment_Nmm; totalV=totalV+volume;
            if factor==1
                old=R.facePoses{f}(k0+j-1); cut=original.faces{f}.removal;
                baseline.maxNormalVectorDifference_N=max(baseline.maxNormalVectorDifference_N, ...
                    max(abs(s.normalVector_N-old.normalVector_N)));
                baseline.maxForceVectorDifference_N=max(baseline.maxForceVectorDifference_N, ...
                    max(abs(s.force_N-old.force_N)));
                baseline.maxVolumeDifference_mm3=max(baseline.maxVolumeDifference_mm3,abs(volume-old.dV_mm3));
                baseline.changedPointSetsAgree=baseline.changedPointSetsAgree && isequal(find(active),cut.index(:));
                if ~isempty(cut.index)
                    baseline.maxSavedActiveHeightDifference_mm=max(baseline.maxSavedActiveHeightDifference_mm, ...
                        max(abs(after(cut.index)-cut.after_mm)));
                end
            end
            Qs{f}=Q;
        end
        records{j}=struct('t_s',tb,'dt_s',odt,'faces',{faceData}, ...
            'combined',struct('normalVector_N',totalN,'force_N',totalF,'moment_Nmm',totalM, ...
            'dV_mm3',totalV,'Fn_pure_resultant_N',norm(totalN)));
    end
    if factor==1
        baseline.status='PASS';
        assert(baseline.changedPointSetsAgree && baseline.maxNormalVectorDifference_N<1e-6 && ...
            baseline.maxForceVectorDifference_N<1e-6 && baseline.maxVolumeDifference_mm3<1e-9 && ...
            baseline.maxSavedActiveHeightDifference_mm<1e-8, 'R30B:EntryBaseline', ...
            '未能复现原结果；请检查当前求解器是否与原计算版本一致。停止加密。');
    end
    % 先平均向量，再取模；不直接比较不同时间宽度的标量峰值。
    commonN=zeros(ncoarse,3); commonF=commonN; commonV=zeros(ncoarse,1);
    peak=0; peakTime=NaN; impulse=[0 0 0];
    for k=1:ncoarse
        pick=(k-1)*factor+(1:factor); nv=zeros(factor,3); fv=nv;
        for u=1:factor
            s=records{pick(u)}.combined; nv(u,:)=s.normalVector_N; fv(u,:)=s.force_N;
            commonV(k)=commonV(k)+s.dV_mm3; impulse=impulse+s.normalVector_N*odt;
            if s.Fn_pure_resultant_N>peak,peak=s.Fn_pure_resultant_N;peakTime=records{pick(u)}.t_s;end
        end
        commonN(k,:)=mean(nv,1); commonF(k,:)=mean(fv,1);
    end
    runs{r}=struct('factor',factor,'output_dt_us',odt*1e6, ...
        'substep_dt_us',odt/o.substeps*1e6,'peakPureNormal_N',peak, ...
        'peakIntervalEnd_s',peakTime,'normalImpulse_Ns',impulse, ...
        'commonBinNormalVector_N',commonN,'commonBinForceVector_N',commonF, ...
        'commonBinVolume_mm3',commonV,'records',{records});
end
fine=runs{end}; scale=max(sqrt(sum(fine.commonBinNormalVector_N.^2,2)));
for r=1:numel(runs)
    d=runs{r}.commonBinNormalVector_N-fine.commonBinNormalVector_N;
    runs{r}.maxCommonBinNormalVectorDifference_N=max(sqrt(sum(d.^2,2)));
    runs{r}.maxCommonBinNormalVectorDifference_overFinestPeak= ...
        runs{r}.maxCommonBinNormalVectorDifference_N/max(scale,realmin);
end
report=struct('status','finite_comparison_no_automatic_convergence_claim', ...
    'engineVersion',R.engineVersion,'faces',R.faces,'window_s',window_s, ...
    'fullSpatialGridUsed',true,'newPointsMayEnterBetweenOriginalSamples',true, ...
    'priorHistoryRecalculated',false,'spatialResolutionChanged',false, ...
    'forceCoefficientsCalibrated',false,'baseline',baseline,'runs',{runs}, ...
    'note','固定原局部起始材料状态；完整空间网格上的短窗口时间加密。未验证前序历史、空间网格、条带宽度或物理模型。');
save(fullfile(outDir,'entry_time_check.mat'),'report','-v7');
write_json(fullfile(outDir,'entry_time_check.json'),report);
write_csv(fullfile(outDir,'entry_time_check.csv'),runs,R.faces);
names={'entry_time_check.mat','entry_time_check.json','entry_time_check.csv'};
try
    fig=figure('Visible','off','Color','w'); guard=onCleanup(@() close(fig));
    subplot(2,1,1); hold on;
    for r=1:numel(runs)
        records=runs{r}.records;tt=cellfun(@(s)s.t_s,records);
        yy=cellfun(@(s)s.combined.Fn_pure_resultant_N,records);
        plot(tt,yy,'DisplayName',sprintf('dt %.4g us',runs{r}.output_dt_us));
    end
    ylabel('Pure normal resultant / N');legend('show');grid on;
    subplot(2,1,2);hold on;tt=R.timeEdges_s(k0+1:k1+1);
    for r=1:numel(runs)
        yy=sqrt(sum(runs{r}.commonBinNormalVector_N.^2,2));
        plot(tt,yy,'DisplayName',sprintf('dt %.4g us',runs{r}.output_dt_us));
    end
    xlabel('Time / s');ylabel('Common-bin mean vector norm / N');grid on;
    saveas(fig,fullfile(outDir,'entry_time_check.png'));names{end+1}='entry_time_check.png';
    clear guard;
catch err
    warning('R30B:EntryPlot','数值结果已保存，绘图未完成：%s',err.message);
end
zipFile=fullfile(outDir,'entry_time_check.zip');zip(zipFile,names,outDir);
fprintf('完成。baseline复现通过。只需上传：\n%s\n',zipFile);
end

function frame=load_frame(root,R,k)
path=fullfile(root,R.frameFiles{k});
assert(isfile(path),'R30B:EntryMissingFrame','缺少原始帧：%s',path);
D=load(path,'frame');frame=D.frame;
assert(frame.index==k && abs(frame.t_s-R.poses(k).t_s)<1e-10, ...
    'R30B:EntryFrameIndex','帧序号或时间不一致。');
end

function s=selected_summary(x)
% 不输出局部重置的ae_pass，避免误认为原整次扫掠余量。
s=struct('normalVector_N',x.normalVector_N,'force_N',x.force_N, ...
    'moment_Nmm',x.moment_Nmm,'dV_mm3',x.dV_mm3, ...
    'Fn_pure_resultant_N',x.Fn_pure_resultant_N);
end

function write_json(path,value)
fid=fopen(path,'w','n','UTF-8');assert(fid>0,'无法写入JSON');
guard=onCleanup(@() fclose(fid));fprintf(fid,'%s',jsonencode(value));
end

function write_csv(path,runs,faces)
fid=fopen(path,'w');assert(fid>0,'无法写入CSV');guard=onCleanup(@() fclose(fid));
fprintf(fid,'factor,tooth_id,flank_side,t_s,dt_s,Fn_pure_resultant_N,Nx_N,Ny_N,Nz_N,Fx_N,Fy_N,Fz_N,Vdot_mm3ps\n');
for r=1:numel(runs)
    records=runs{r}.records;
    for k=1:numel(records)
        q=records{k};
        for f=1:size(faces,1)+1
            if f<=size(faces,1),s=q.faces{f};face=faces(f,:);else,s=q.combined;face=[0 0];end
            fprintf(fid,'%d,%d,%d,%.15g,%.15g,%.15g,%.15g,%.15g,%.15g,%.15g,%.15g,%.15g,%.15g\n', ...
                runs{r}.factor,face(1),face(2),q.t_s,q.dt_s,s.Fn_pure_resultant_N, ...
                s.normalVector_N,s.force_N,s.dV_mm3/q.dt_s);
        end
    end
end
end
