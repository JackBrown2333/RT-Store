function zipFiles = r30b_pack_diagnostics(resultLocation, maxPartMB)
% R30B_PACK_DIAGNOSTICS 在本机核验大结果，只打包少量原始关键帧。
% 用法：把本文件放到原代码包根目录，在该目录执行：
%   zipFiles = r30b_pack_diagnostics;
% 也可指定结果文件夹或MAT文件，以及每个独立ZIP的目标大小（MiB）：
%   zipFiles = r30b_pack_diagnostics('r30_dual_recomputed',15);
%
% 本工具不运行材料/力求解器，不修改或删除原始结果，不读取检查点。
% 完整材料回放核验在本机执行。上传包只包含汇总、核验报告和关键帧，
% 不包含完整材料初末状态，不能作为全程动画或全程回放的数据包。

if nargin < 1 || isempty(resultLocation)
    resultLocation = fullfile(pwd,'r30_dual_recomputed');
end
if nargin < 2 || isempty(maxPartMB), maxPartMB = 15; end
assert(isscalar(maxPartMB) && isfinite(maxPartMB) && maxPartMB >= 0.125, ...
    'R30B:PackSize','每包目标大小至少为0.125 MiB。');
if isfolder(resultLocation)
    resultFile = fullfile(resultLocation,'r30b_result.mat');
else
    resultFile = resultLocation;
end
assert(isfile(resultFile),'R30B:PackInput','未找到结果，请指定r30b_result.mat或所在文件夹。');
[ok,attr] = fileattrib(resultFile); assert(ok); resultFile = attr.Name;
root = fileparts(resultFile);
D = load(resultFile,'R30B'); R = D.R30B;
assert(strcmp(R.engineVersion,'r30_dual_from_zero_v2') && ~isempty(R.frameFiles), ...
    'R30B:PackVersion','需要本次新版双侧计算且已保存几何帧。');

% 核验器来自已有求解代码；本工具须放在包根目录，不要放入matlab_recompute。
here = fileparts(mfilename('fullpath'));
solverDirs = {fullfile(here,'matlab_recompute'), ...
    fullfile(fileparts(root),'matlab_recompute')};
for k = 1:numel(solverDirs)
    if isfolder(solverDirs{k}), addpath(solverDirs{k}); end
end

% 每次新建输出目录，避免覆盖以前的诊断包或计算结果。
stem = fullfile(root,['diagnostics_' datestr(now,'yyyymmdd_HHMMSS')]);
out = stem; suffix = 1;
while isfolder(out), out = sprintf('%s_%d',stem,suffix); suffix = suffix+1; end
mkdir(out); mkdir(fullfile(out,'frames'));

fprintf('正在本机读取完整结果并核验；不重新计算磨削过程。\n');
fprintf('这一步需要顺序读取frames，1 GB数据可能需要数分钟。\n');
verification = struct('status','NOT_RUN');
try
    assert(exist('r30b_verify','file') == 2,'R30B:MissingVerifier', ...
        '未找到原包matlab_recompute/r30b_verify.m，请把本工具放在原代码包根目录。');
    verification = r30b_verify(resultFile);
catch err
    % 核验失败也保留诊断材料，便于定位错误，绝不把失败写成通过。
    verification = struct('status','FAIL_OR_INCOMPLETE', ...
        'identifier',err.identifier,'message',err.message);
    warning('R30B:PackVerification','本机核验未通过或未完成：%s。仍导出关键帧。',err.message);
end
write_json(fullfile(out,'verification_report.json'),verification);

% 自动选取：各侧首次入切、入切早期峰值及其邻帧、峰后、全局峰值、
% 退切前后、双侧载荷最接近的时刻，以及首末保存帧。
nt = numel(R.poses); dt = R.options.dt_s; desired = [];
normalLoads = zeros(nt,size(R.faces,1));
for f = 1:size(R.faces,1)
    S = R.facePoses{f}; normalLoads(:,f) = [S.Fn_pure_resultant_N]';
    active = find([S.dV_mm3] > 0);
    if isempty(active), continue; end
    first = active(1); last = active(end);
    early = first:min(last,first+max(1,round(0.005/dt)));
    [~,j] = max(normalLoads(early,f)); earlyPeak = early(j);
    [~,globalPeak] = max(normalLoads(:,f));
    desired = [desired, earlyPeak+(-2:2), first+[-1 0 1], ...
        earlyPeak+max(1,round(0.001/dt)), globalPeak, last+[-1 0 1]]; %#ok<AGROW>
end
if size(normalLoads,2) == 2
    mask = normalLoads(:,1) > max(normalLoads(:,1))*0.25 & ...
           normalLoads(:,2) > max(normalLoads(:,2))*0.25;
    candidates = find(mask);
    if ~isempty(candidates)
        ratio = abs(normalLoads(candidates,1)-normalLoads(candidates,2)) ./ ...
            sum(normalLoads(candidates,:),2);
        [~,j] = min(ratio); desired(end+1) = candidates(j);
    end
end
available = double(R.frameIndices(:));
desired = [desired,available(1),available(end)];
desired = unique(min(nt,max(1,round(desired))),'stable');
% 如果原计算曾抽帧，映射到实际存在的最近一帧，并在清单中写明。
chosen = zeros(size(desired));
for k = 1:numel(desired), [~,chosen(k)] = min(abs(available-desired(k))); end
chosen = unique(chosen,'stable');

names = cell(1,numel(chosen)); sizes = zeros(size(chosen)); selected = struct([]);
for k = 1:numel(chosen)
    j = chosen(k); index = available(j);
    rel = strrep(R.frameFiles{j},'\',filesep); rel = strrep(rel,'/',filesep);
    source = fullfile(root,rel);
    assert(isfile(source),'R30B:PackFrame','缺少原始帧：%s',source);
    names{k} = sprintf('frames/frame_%06d.mat',index);
    % 原样复制，保留全部材料before/after、子步、条带力和尺寸标注数据。
    copyfile(source,fullfile(out,names{k})); info = dir(source); sizes(k) = info.bytes;
    selected(k).index = index; selected(k).t_s = R.poses(index).t_s; %#ok<AGROW>
    selected(k).file = names{k}; selected(k).source_bytes = sizes(k);
end
copyfile(resultFile,fullfile(out,'r30b_result.mat'));
manifest = struct('schemaVersion','r30_diagnostic_selection_v1', ...
    'scope','Local full verification report plus lossless selected frames; not a full replay bundle', ...
    'engineVersion',R.engineVersion,'faces',R.faces,'window_s',R.window_s, ...
    'requestedIntervalIndices',desired,'selectedFrames',selected, ...
    'targetPartMiB',maxPartMB,'fullMaterialStatesIncluded',false, ...
    'fullHistoryRecalculated',false,'originalFilesModified',false);
write_json(fullfile(out,'selection.json'),manifest);
fid = fopen(fullfile(out,'README_zh.txt'),'w','n','UTF-8'); assert(fid>0);
fprintf(fid,['这是少量原始关键帧的诊断包，不是完整动画包。\n' ...
    'verification_report.json记录在原始完整目录执行的本机核验结果。\n' ...
    'r30b_result.mat保留全部汇总曲线；它引用的多数帧不在本包中。\n' ...
    '所选frames逐字节原样复制，未降精度、平滑或删减帧内数据。\n' ...
    'stateUpdates仍是原始保存间隔的增量，没有合并跳过帧的材料更新。\n' ...
    '不能用这些孤立帧从头回放完整材料表面，不能在本包上再次运行完整核验。\n' ...
    '可检查选中区间内的before/after、substepUpdates、条带力和几何。\n' ...
    '各ZIP是独立压缩包，不是需要拼接的分卷；全部上传即可。\n' ...
    '原始frames、initial_state和final_state请保留在本机供后续动画使用。\n']);
fclose(fid);

% 按源MAT压缩文件大小分组；留出ZIP目录开销，每包重复附上报告和汇总。
common = {'r30b_result.mat','verification_report.json','selection.json','README_zh.txt'};
commonBytes = 0;
for k = 1:numel(common), info = dir(fullfile(out,common{k})); commonBytes = commonBytes+info.bytes; end
budget = maxPartMB*1024^2-8*1024;
assert(commonBytes < budget,'R30B:PackMetadata','汇总文件超过目标大小，请增大第二个参数。');
groups = {}; current = {}; used = commonBytes;
for k = 1:numel(names)
    if used+sizes(k) > budget && ~isempty(current)
        groups{end+1} = current; current = {}; used = commonBytes; %#ok<AGROW>
    end
    if commonBytes+sizes(k) > budget
        warning('R30B:PackLargeFrame','单帧超过目标大小，将原样单独打包：%s。',names{k});
    end
    current{end+1} = names{k}; used = used+sizes(k); %#ok<AGROW>
end
if ~isempty(current), groups{end+1} = current; end
zipFiles = cell(numel(groups),1);
for k = 1:numel(groups)
    zipFiles{k} = fullfile(out,sprintf('r30_diagnostic_%02d.zip',k));
    zip(zipFiles{k},[common,groups{k}],out);
    info = dir(zipFiles{k});
    fprintf('上传：%s  （%.2f MiB）\n',zipFiles{k},info.bytes/1024^2);
end
fprintf('完成：选取%d/%d个已保存帧，生成%d个独立ZIP。\n',numel(chosen),numel(available),numel(zipFiles));
fprintf('核验状态：%s。只需上传上面列出的ZIP；原始frames保留在本机。\n',verification.status);
end

function write_json(filename,value)
% 显式UTF-8，保证中文路径和错误信息可在不同系统读取。
fid = fopen(filename,'w','n','UTF-8'); assert(fid>0);
cleanup = onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s',jsonencode(value));
end
