function fixedCore=r30full_prepare_corefix(originalCore,fixedCore)
% Copy complete, prebuilt functions; never splice MATLAB source at runtime.
root=fileparts(mfilename('fullpath'));patchDir=fullfile(root,'geometry_fix_v3');
geometry=fileread(fullfile(patchDir,'r30b_geometry.m'));
strips=fileread(fullfile(patchDir,'r30b_strips.m'));
reference=fileread(fullfile(patchDir,'original_r30b_strips_reference.txt'));
originalStrips=fileread(fullfile(originalCore,'r30b_strips.m'));
assert(strcmp(canonical(originalStrips),canonical(reference)),...
 'R30FULL:CoreCompatibility','原r30b_strips与已核验版本不同。原文件保持不变，请提供原r30b_strips.m核对。');
patchIdentity=struct('version','static_geometry_fix_v3','originalSignature',signature(originalCore),...
 'geometrySource',geometry,'stripsSource',strips);
manifestPath=fullfile(fixedCore,'geometry_fix_manifest.mat');
if isfolder(fixedCore)
 assert(isfile(manifestPath),'R30FULL:FixFolder','修复核心目录缺少版本清单，请使用新的专用目录。');
 D=load(manifestPath,'patchIdentitySaved','outputSignature');
 assert(isequaln(D.patchIdentitySaved,patchIdentity)&&strcmp(D.outputSignature,signature(fixedCore)),...
  'R30FULL:FixVersion','修复核心内容不一致，请使用新的专用目录。');
else
 parent=fileparts(fixedCore);if ~isfolder(parent),mkdir(parent);end
 stage=tempname(parent);mkdir(stage);
 files=dir(fullfile(originalCore,'*.m'));
 for k=1:numel(files),copyfile(fullfile(originalCore,files(k).name),fullfile(stage,files(k).name));end
 copyfile(fullfile(patchDir,'r30b_geometry.m'),fullfile(stage,'r30b_geometry.m'));
 copyfile(fullfile(patchDir,'r30b_strips.m'),fullfile(stage,'r30b_strips.m'));
 parse_core(stage);
 patchIdentitySaved=patchIdentity;outputSignature=signature(stage);
 save(fullfile(stage,'geometry_fix_manifest.mat'),'patchIdentitySaved','outputSignature','-v7');
 rmpath(stage);clear r30b_geometry r30b_strips r30b_kinematics;
 [ok,msg]=movefile(stage,fixedCore);assert(ok,'R30FULL:FixInstall','%s',msg);
end
parse_core(fixedCore);
fprintf('完整函数文件及MATLAB解析检查通过：%s\n',fixedCore);
end
function parse_core(folder)
addpath(folder,'-begin');rehash;clear r30b_geometry r30b_strips r30b_kinematics;
assert(strcmp(fileparts(which('r30b_geometry')),folder)&&strcmp(fileparts(which('r30b_strips')),folder),...
 'R30FULL:FixPath','MATLAB仍指向其他几何函数。');
% NARGIN/NARGOUT load and parse M-files without advancing material state.
assert(nargin('r30b_strips')==8&&nargout('r30b_strips')==1,...
 'R30FULL:StripAPI','Unexpected strip function interface.');
assert(nargin('r30b_geometry')==1&&nargout('r30b_geometry')==4,...
 'R30FULL:GeometryAPI','Unexpected polygon geometry interface.');
end
function s=canonical(s)
% Known reference has no percent character in string literals.
if ~isempty(s)&&double(s(1))==65279
 s=s(2:end);
elseif numel(s)>=3&&isequal(double(s(1:3)),[239 187 191])
 s=s(4:end);
end
s=regexprep(s,'%[^\r\n]*','');s=regexprep(s,'\s','');
end
function s=signature(folder)
f=dir(fullfile(folder,'*.m'));assert(~isempty(f));s='';
for j=1:numel(f),s=[s f(j).name char(10) fileread(fullfile(folder,f(j).name))];end
end
