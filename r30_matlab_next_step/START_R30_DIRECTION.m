% 打开本文件并按F5。自动定位原工程，缺少路径时弹出目录选择窗口。
r30LauncherDir=fileparts(mfilename('fullpath'));
assert(~isempty(r30LauncherDir),'R30:Launcher','请打开START_R30_DIRECTION.m并按F5，不要只复制其内容。');
addpath(r30LauncherDir,'-begin');rehash;
clear r30d_wheel_geometry run_r30_contact_direction_audit;
R30_DIRECTION_RESULT=struct('completed',false);
test_r30d_shapes(fullfile(r30LauncherDir,'example_data','multitooth_pilot.mat'));
r30Setup=r30f_setup_paths('direction');
fprintf('诊断入口：%s\n',which('run_r30_contact_direction_audit'));
try
 r30AuditResult=run_r30_contact_direction_audit(r30Setup.resultDir,struct('coreDir',r30Setup.coreDir));
catch r30AuditError
 fprintf(2,'本次诊断未完成，尚无本次可上传的结果ZIP；原工作区中的A不代表本次诊断结果。\n');
 rethrow(r30AuditError);
end
assert(isfile(r30AuditResult.resultZipFile),'R30:MissingZip','诊断返回但未找到结果ZIP。');
R30_DIRECTION_RESULT=r30AuditResult;R30_DIRECTION_RESULT.completed=true;
A=R30_DIRECTION_RESULT;
fprintf('已成功生成本次结果：\n%s\n',R30_DIRECTION_RESULT.resultZipFile);
if ispc
 try,winopen(R30_DIRECTION_RESULT.outputDir);catch,fprintf('结果目录：%s\n',R30_DIRECTION_RESULT.outputDir);end
end
