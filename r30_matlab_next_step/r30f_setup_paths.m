function C=r30f_setup_paths(task,allowPicker)
% Resolve only known package/project locations; never recursively scan a disk.
if nargin<2,allowPicker=true;end
assert(any(strcmp(task,{'plots','direction','check'})));
root=fileparts(mfilename('fullpath'));addpath(root,'-begin');rehash;
required={'run_r30_academic_figures.m','run_r30_contact_direction_audit.m','r30f_spectrum.m','r30f_thesis_figures.m'};
for j=1:numel(required),assert(isfile(fullfile(root,required{j})),'R30:Incomplete','工具包不完整，缺少%s。请完整解压ZIP。',required{j});end
base={root,fileparts(root),fileparts(fileparts(root)),pwd};snapshot=fullfile(root,'example_data','multitooth_pilot.mat');
originalResult='';
if isfile(snapshot)
 D=load(snapshot,'M');
 if isfield(D.M,'options')&&isfield(D.M.options,'outputDir'),originalResult=D.M.options.outputDir;base{end+1}=fileparts(originalResult);end
 if isfield(D.M,'options')&&isfield(D.M.options,'parameterFile'),base{end+1}=fileparts(D.M.options.parameterFile);end
end
project='';for j=1:numel(base),if ~isempty(base{j})&&isfolder(fullfile(base{j},'matlab_recompute')),project=base{j};break;end;end
if strcmp(task,'direction')&&isempty(project)
 assert(allowPicker,'R30:Project','未找到原工程目录。');
 project=uigetdir(root,'选择原工程根目录（里面应有matlab_recompute文件夹）');
 assert(~isequal(project,0),'R30:Cancelled','已取消目录选择。');
 assert(isfolder(fullfile(project,'matlab_recompute')),'R30:Project','所选目录内没有matlab_recompute，请选工程根目录。');
end
result='';candidates={};
if ~isempty(project),candidates{end+1}=fullfile(project,'r30_multitooth_pilot_15_v2_160x320');end
if ~isempty(originalResult),candidates{end+1}=originalResult;end
for j=1:numel(candidates),if isfile(fullfile(candidates{j},'multitooth_pilot.mat')),result=candidates{j};break;end;end
if isempty(result)&&strcmp(task,'direction')
 assert(allowPicker,'R30:Result','未找到完整多齿结果目录。');
 result=uigetdir(project,'选择完整多齿结果目录（里面应有multitooth_pilot.mat和animation_export）');
 assert(~isequal(result,0)&&isfile(fullfile(result,'multitooth_pilot.mat')),'R30:Result','所选目录缺少multitooth_pilot.mat。');
elseif isempty(result)
 result=fileparts(snapshot);assert(isfile(snapshot),'R30:Snapshot','找不到原结果或包内example_data/multitooth_pilot.mat。');
end
C=struct('packageDir',root,'projectDir',project,'resultDir',result,'coreDir','',...
 'figureFunction',which('run_r30_academic_figures'),'directionFunction',which('run_r30_contact_direction_audit'));
if ~isempty(project),C.coreDir=fullfile(project,'matlab_recompute');end
assert(strcmp(fileparts(C.figureFunction),root)&&strcmp(fileparts(C.directionFunction),root),'R30:Shadowed','MATLAB仍指向其他版本，请重新打开本包启动脚本。');
fprintf('工具目录：%s\n数据来源：%s\n',root,result);
end
