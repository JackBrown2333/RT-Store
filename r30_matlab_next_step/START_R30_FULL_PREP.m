% 将本启动脚本和两个run_*.m放入已有r30_matlab_next_step目录，打开本文件按F5。
% 本次运行：全齿行程几何预检 + 短窗口收敛对照。尚不计算全程力。
r30PrepLauncherDir=fileparts(mfilename('fullpath'));
assert(~isempty(r30PrepLauncherDir),'R30P:Launcher','请打开START_R30_FULL_PREP.m并按F5。');
addpath(r30PrepLauncherDir,'-begin');rehash;
assert(isfile(fullfile(r30PrepLauncherDir,'r30f_setup_paths.m')),'R30P:Install',...
 '请把这3个M文件放入原有r30_matlab_next_step工具目录，再运行本文件。');
clear run_r30_full_preflight run_r30_resolution_check;
R30_FULL_PREP_RESULT=struct('completed',false,'fullForceComputed',false);
r30PrepSetup=r30f_setup_paths('direction');
r30PrepOut=fullfile(r30PrepSetup.resultDir,'full_process_preparation_v1');
if ~isfolder(r30PrepOut),mkdir(r30PrepOut);end
r30PrepZip=fullfile(r30PrepOut,'r30_full_preparation_result.zip');
if isfile(r30PrepZip),delete(r30PrepZip);end
try
 r30PrepPlan=run_r30_full_preflight(fullfile(r30PrepSetup.resultDir,'multitooth_pilot.mat'),...
  r30PrepSetup.coreDir,fullfile(r30PrepOut,'preflight'));
 r30PrepCheck=run_r30_resolution_check(r30PrepSetup.resultDir,r30PrepSetup.coreDir,...
  fullfile(r30PrepOut,'resolution'));
 r30PrepSummary=struct('completed',true,'fullForceComputed',false,'outputDir',r30PrepOut,...
  'resultZipFile',r30PrepZip,'plan',r30PrepPlan,'resolution',r30PrepCheck);
 save(fullfile(r30PrepOut,'preparation_summary.mat'),'r30PrepSummary','-v7');
 r30PrepFiles={'preparation_summary.mat'};
 for r30PrepSub={'preflight','resolution'}
  r30PrepList=dir(fullfile(r30PrepOut,r30PrepSub{1}));
  for r30PrepJ=1:numel(r30PrepList)
   r30PrepItem=r30PrepList(r30PrepJ);if r30PrepItem.isdir,continue;end
   % 输入/代码签名计划和大材料检查点不上传；保留报告、曲线和图。
   if strcmp(r30PrepItem.name,'resolution_plan.mat'),continue;end
   r30PrepFiles{end+1}=fullfile(r30PrepSub{1},r30PrepItem.name);
  end
 end
 zip(r30PrepZip,r30PrepFiles,r30PrepOut);
 assert(isfile(r30PrepZip),'R30P:Zip','结果ZIP未生成。');
 R30_FULL_PREP_RESULT=r30PrepSummary;
catch r30PrepError
 fprintf(2,'本次准备未完成。已完成的齿/窗口保留检查点，可用同一入口继续。\n');
 fprintf(2,'本次尚未生成有效结果ZIP；也没有计算全程力。\n');
 rethrow(r30PrepError);
end
fprintf('\n本次准备完成，尚未计算全程力。请上传：\n%s\n',r30PrepZip);
fprintf('收敛报告列出差异，不会把“程序运行成功”当成“精度通过”。\n');
if ispc
 try,winopen(r30PrepOut);catch,fprintf('结果目录：%s\n',r30PrepOut);end
end
