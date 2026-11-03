function report=test_r30full_static_install(originalCore,workDir)
% Ensure source layout/line endings cannot affect the complete installed file.
root=fileparts(fileparts(mfilename('fullpath')));patch=fullfile(root,'geometry_fix_v3');
if ~isfolder(workDir),mkdir(workDir);end
source=fileread(fullfile(patch,'original_r30b_strips_reference.txt'));
variant=fullfile(workDir,'ascii_crlf_original');mkdir(variant);
files=dir(fullfile(originalCore,'*.m'));
for k=1:numel(files),copyfile(fullfile(originalCore,files(k).name),fullfile(variant,files(k).name));end
path=fullfile(variant,'r30b_strips.m');fid=fopen(path,'wb');assert(fid>0);
fwrite(fid,strrep(source,char(10),[char(13) char(10)]),'char');fclose(fid);
fixed=fullfile(workDir,'complete_static_core');r30full_prepare_corefix(variant,fixed);
assert(strcmp(fileread(fullfile(fixed,'r30b_strips.m')),fileread(fullfile(patch,'r30b_strips.m'))));
r30full_prepare_corefix(variant,fixed);assert(nargin('r30b_strips')==8&&nargout('r30b_geometry')==4);
fid=fopen(path,'wb');fwrite(fid,[239 187 191],'uint8');
fwrite(fid,strrep(source,char(10),[char(13) char(10)]),'char');fclose(fid);
r30full_prepare_corefix(variant,fullfile(workDir,'bom_static_core'));
% An unseen computation edit must be rejected before installation.
fid=fopen(path,'wb');fwrite(fid,strrep(source,'width=o.forceStripWidth_mm;','width=2*o.forceStripWidth_mm;'),'char');fclose(fid);
rejected=false;
try,r30full_prepare_corefix(variant,fullfile(workDir,'must_not_install'));
catch e,rejected=strcmp(e.identifier,'R30FULL:CoreCompatibility');end
assert(rejected&&~isfolder(fullfile(workDir,'must_not_install')));
report=struct('completeSourceCopiedVerbatim',true,'asciiCRLFInputAccepted',true,...
 'utf8BOMInputAccepted',true,...
 'installedFileParsedBeforeUse',true,'cachedInstallRechecked',true,'unseenCoreEditRejected',true,...
 'nativeMatlabExecutedHere',false);disp(report);
end
