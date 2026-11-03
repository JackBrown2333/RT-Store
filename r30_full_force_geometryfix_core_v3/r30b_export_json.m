function r30b_export_json(resultFile,filename,stride)
% Portable animation bundle. When dropping frames, merge material updates so
% the remaining frames still reproduce the saved material surface exactly.
if nargin<3,stride=1;end
assert(stride>=1&&mod(stride,1)==0);D=load(resultFile,'R30B');R=D.R30B;root=fileparts(resultFile);
assert(strcmp(R.engineVersion,'r30_dual_from_zero_v2')&&~isempty(R.frameFiles),'R30B:Frames','No new-engine geometry frames saved');
D=load(fullfile(root,R.initialStateFile),'initialState');initial=D.initialState;nf=size(R.faces,1);
height=cell(nf,1);previous=cell(nf,1);for f=1:nf,height{f}=initial.faces{f}.height_mm;previous{f}=height{f};end
fid=fopen(filename,'w','n','UTF-8');assert(fid>0,'R30B:Output','Cannot write export');cleanup=onCleanup(@()fclose(fid));
header=struct('schemaVersion','r30_dual_geometry_v2','parameters',R.parameters,'options',R.options,...
 'faces',R.faces,'window_s',R.window_s,'initialState',initial,'capabilities',...
 struct('independentMaterialHistory',true,'sparseMaterialReplay',true,'cuttingSupportRectangles',true,...
 'parameterDimensionAnchors',true,'localForceVectors',true,'exactContinuousIntersection',false),...
 'geometryRepresentation',R.geometryRepresentation,'forceTimeMeaning',R.forceTimeMeaning,...
 'coefficientMeaning',R.coefficientMeaning,'fullGearComputed',false,...
 'materialIndexing','MATLAB 1-based column-major: row=mod(index-1,nRows)+1; col=floor((index-1)/nRows)+1');
encoded=jsonencode(header);fprintf(fid,'%s,"frames":[',encoded(1:end-1));first=true;lastTime=R.window_s(1);
for j=1:numel(R.frameFiles)
 D=load(fullfile(root,R.frameFiles{j}),'frame');frame=D.frame;
 for f=1:nf,u=frame.faces{f}.stateUpdates;h=height{f};h(u.index)=u.height_mm;height{f}=h;end
 if mod(j-1,stride)~=0&&j~=numel(R.frameFiles),continue;end
 for f=1:nf
  ids=find(height{f}<previous{f});frame.faces{f}.stateUpdates=struct('since_t_s',lastTime,'index',ids,'height_mm',height{f}(ids));previous{f}=height{f};
 end
 if ~first,fprintf(fid,',');end;fprintf(fid,'%s',jsonencode(frame));first=false;lastTime=frame.t_s;
end
fprintf(fid,']}');
end
