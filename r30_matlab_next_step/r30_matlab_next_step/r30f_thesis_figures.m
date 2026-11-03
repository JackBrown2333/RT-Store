function names=r30f_thesis_figures(M,S,o,out)
% Reference: uploaded thesis, PDF p40 Fig3-7 and p54 Figs4-7..4-9.
% Change presentation only; data and spectral processing are unchanged.
o.chineseFontName=choose_cjk_font(o.chineseFontName);
C=[0 .145 1;1 .11 0;0 0 0;.05 .55 .15;.65 0 .65];
t=double(M.time_s(:))-double(M.window_s(1));limits=[0 numel(t)*double(M.forceDt_s)];
ids=double(M.toothIds(:)');n=numel(t);g=find(ids==o.singleTooth);pair=arrayfun(@(x)find(ids==x),o.pairTeeth);
F1=reshape(M.toothForce_N(:,pair(1),:),n,3);F2=reshape(M.toothForce_N(:,pair(2),:),n,3);Fp=F1+F2;
a1=vecnorm(F1,2,2);a2=vecnorm(F2,2,2);ap=vecnorm(Fp,2,2);upper=a1+a2;
eta=NaN(n,1);loaded=upper>1e-12;eta(loaded)=1-ap(loaded)./upper(loaded);
fMesh=abs(M.parameters.installation.omegaG)/(2*pi)*M.parameters.p.gear.teeth;
names={};
[h,ax]=paperfig(o,1,2);lines(ax(1),t,M.force_N,C,o);labels(ax(1),'时间','(s)','力','(N)',o);
leg(ax(1),{'F_x','F_y','F_z'},o);caption(h,ax(1),'(a)','全局力分量',o);
lines(ax(2),t,[vecnorm(M.force_N,2,2),vecnorm(M.normalVector_N,2,2)],C,o);labels(ax(2),'时间','(s)','合力','(N)',o);
leg(ax(2),{zh('全力合力',o),zh('纯法向合力',o)},o);caption(h,ax(2),'(b)','全力与纯法向合力',o);
names{end+1}=finish(h,ax,out,'F01_global_force',o,limits);

[h,ax]=paperfig(o,1,2);lines(ax(1),S.frequency_Hz(2:end),S.amplitude(2:end,1:3),C,o);
labels(ax(1),'频率','(Hz)','幅值','(N)',o);leg(ax(1),{'F_x','F_y','F_z'},o);meshline(ax(1),fMesh);caption(h,ax(1),'(a)','力分量幅值谱',o);
lines(ax(2),S.frequency_Hz(2:end),S.amplitude(2:end,4),C,o);labels(ax(2),'频率','(Hz)','幅值','(N\cdot m)',o);
meshline(ax(2),fMesh);caption(h,ax(2),'(b)','C轴转矩幅值谱',o);
names{end+1}=finish(h,ax,out,'F02_global_FFT',o,[0 min(o.fftMaxHz,1/(2*M.forceDt_s))]);

names{end+1}=single(M,g,t,limits,C,o,out);
if o.allActiveTeeth
 active=find(any(M.toothVdot_mm3ps>0,1));for k=active,if k~=g,names{end+1}=single(M,k,t,limits,C,o,out);end;end
end
[h,ax]=paperfig(o,2,2);
lines(ax(1),t,[a1 a2 ap],C,o);labels(ax(1),'时间','(s)','力','(N)',o);
leg(ax(1),{toothlabel(o.pairTeeth(1),o),toothlabel(o.pairTeeth(2),o),zh('矢量合力',o)},o);caption(h,ax(1),'(a)','单齿及双齿合力',o);
lines(ax(2),t,[upper ap],C,o);labels(ax(2),'时间','(s)','力','(N)',o);
leg(ax(2),{zh('模长之和',o),zh('矢量合力',o)},o);caption(h,ax(2),'(b)','标量上界与矢量合力',o);
lines(ax(3),t,Fp,C,o);labels(ax(3),'时间','(s)','力','(N)',o);leg(ax(3),{'F_x','F_y','F_z'},o);caption(h,ax(3),'(c)','双齿合力分量',o);
lines(ax(4),t,eta,C,o);labels(ax(4),'时间','(s)','相消指标','',o);ylim(ax(4),[0 1]);caption(h,ax(4),'(d)','方向相消指标',o);
names{end+1}=finish(h,ax,out,sprintf('F04_pair_%02d_%02d_superposition',o.pairTeeth),o,limits);

[h,ax]=paperfig(o,2,1);perN=sqrt(sum(M.toothNormalVector_N.^2,3));
imagesc(ax(1),t,1:numel(ids),perN');set(ax(1),'YDir','normal','YTick',1:numel(ids),'YTickLabel',arrayfun(@num2str,ids,'UniformOutput',false));
labels(ax(1),'时间','(s)','齿号','',o);colormap(h,parula(128));cb=colorbar(ax(1));cb.FontName=o.fontName;cb.FontSize=o.fontSize;
ylabel(cb,[zh('力',o) latin('(N)',o)],'Interpreter','tex');caption(h,ax(1),'(a)','各齿纯法向合力',o);
lines(ax(2),t,[sum(perN,2) vecnorm(M.normalVector_N,2,2)],C,o);labels(ax(2),'时间','(s)','力','(N)',o);
leg(ax(2),{zh('各齿模长之和',o),zh('多齿矢量合力',o)},o);caption(h,ax(2),'(b)','多齿法向载荷叠加',o);
names{end+1}=finish(h,ax,out,'F05_per_tooth_normal_contributions',o,limits);

[h,ax]=paperfig(o,1,2);stairs(ax(1),t,double(M.nRemovingTeeth),'Color',C(1,:),'LineWidth',o.lineWidth);hold(ax(1),'on');
stairs(ax(1),t,double(M.nRemovingFaces),'Color',C(2,:),'LineWidth',o.lineWidth);
labels(ax(1),'时间','(s)','数量','',o);leg(ax(1),{zh('去除齿数',o),zh('去除齿面数',o)},o);ylim(ax(1),[0 max(double(M.nRemovingFaces))+1]);caption(h,ax(1),'(a)','同时参与去除的齿与齿面',o);
lines(ax(2),t,M.moment_Nmm(:,3)/1000,C,o);labels(ax(2),'时间','(s)','转矩','(N\cdot m)',o);caption(h,ax(2),'(b)','C轴载荷转矩',o);
names{end+1}=finish(h,ax,out,'F06_removing_counts_and_torque',o,limits);

% Thesis Fig3-7 layout: signed torque time history next to its spectrum.
[h,ax]=paperfig(o,1,2);lines(ax(1),t,M.moment_Nmm(:,3)/1000,C,o);labels(ax(1),'时间','(s)','转矩','(N\cdot m)',o);
caption(h,ax(1),'(a)','时域转矩',o);lines(ax(2),S.frequency_Hz(2:end),S.amplitude(2:end,4),C,o);
labels(ax(2),'频率','(Hz)','幅值','(N\cdot m)',o);meshline(ax(2),fMesh);caption(h,ax(2),'(b)','频域转矩',o);
names{end+1}=finish(h,ax,out,'F07_torque_time_frequency',o,{limits,[0 min(o.fftMaxHz,1/(2*M.forceDt_s))]});
fprintf('论文风格：中文%s，数字%s；蓝红细线、内向刻度、图下分图说明。\n',o.chineseFontName,o.fontName);
end

function name=single(M,g,t,limits,C,o,out)
n=numel(t);tooth=double(M.toothIds(g));F=reshape(M.toothForce_N(:,g,:),n,3);N=reshape(M.toothNormalVector_N(:,g,:),n,3);
faces=find(M.faces(:,1)==tooth);assert(numel(faces)==2);fc=zeros(n,2);tc=zeros(n,2);
for j=1:2,fc(:,j)=sqrt(sum(M.faceForce_N(:,faces(j),:).^2,3));tc(:,j)=M.faceMoment_Nmm(:,faces(j),3)/1000;end
[h,ax]=paperfig(o,2,2);lines(ax(1),t,F,C,o);labels(ax(1),'时间','(s)','力','(N)',o);
leg(ax(1),{'F_x','F_y','F_z'},o);caption(h,ax(1),'(a)',sprintf('第%d齿力分量',tooth),o);
lines(ax(2),t,[vecnorm(F,2,2) vecnorm(N,2,2)],C,o);labels(ax(2),'时间','(s)','力','(N)',o);
leg(ax(2),{zh('全力合力',o),zh('纯法向合力',o)},o);caption(h,ax(2),'(b)','全力与纯法向合力',o);
lines(ax(3),t,fc,C,o);labels(ax(3),'时间','(s)','力','(N)',o);
leg(ax(3),{zh('齿面1',o),zh('齿面2',o)},o);caption(h,ax(3),'(c)','两齿面受力',o);
lines(ax(4),t,[tc sum(tc,2)],C,o);labels(ax(4),'时间','(s)','转矩','(N\cdot m)',o);
leg(ax(4),{zh('齿面1',o),zh('齿面2',o),zh('合计',o)},o);caption(h,ax(4),'(d)','两齿面转矩叠加',o);
name=finish(h,ax,out,sprintf('F03_tooth_%02d_force',tooth),o,limits);
end
function [h,ax]=paperfig(o,nr,nc)
vis='off';if o.visible,vis='on';end
height=6.8;if nr==2,height=12.8;end
h=figure('Visible',vis,'Color','w','Units','centimeters','Position',[2 2 16 height],...
 'PaperPositionMode','auto','InvertHardcopy','off','Renderer','painters');
ax=gobjects(nr*nc,1);
for j=1:nr*nc
 row=ceil(j/nc);col=mod(j-1,nc)+1;
 if nc==1,x=.13;width=.73;else,x=.09+(col-1)*.51;width=.35;end
 if nr==1,y=.30;hh=.64;else,y=.655-(row-1)*.45;hh=.30;end
 ax(j)=axes('Parent',h,'Position',[x y width hh],'FontName',o.fontName,'FontSize',o.fontSize,...
  'LineWidth',.8,'Box','on','TickDir','in','TickLength',[.02 .02],'XGrid','off','YGrid','off');
end
end
function lines(ax,x,Y,C,o)
hold(ax,'on');for j=1:size(Y,2),plot(ax,x,Y(:,j),'Color',C(mod(j-1,size(C,1))+1,:),'LineWidth',o.lineWidth);end
end
function labels(ax,xcn,xunit,ycn,yunit,o)
xlabel(ax,[zh(xcn,o) latin(xunit,o)],'Interpreter','tex','FontSize',o.fontSize,'FontWeight','normal');
ylabel(ax,[zh(ycn,o) latin(yunit,o)],'Interpreter','tex','FontSize',o.fontSize,'FontWeight','normal');
end
function text=zh(s,o),text=['\fontname{' o.chineseFontName '}' s];end
function text=latin(s,o),text=['\fontname{' o.fontName '}' s];end
function text=toothlabel(k,o),text=[zh('第',o) latin(num2str(k),o) zh('齿',o)];end
function leg(ax,strings,o)
lg=legend(ax,strings,'Location','northeast','Box','off','FontName',o.fontName,...
 'FontSize',o.fontSize-1,'Interpreter','tex','Color','none');
components=numel(strings)==3&&all(cellfun(@(s)~isempty(regexp(s,'^F_[xyz]$','once')),strings));
if components,set(lg,'NumColumns',3);pad=.20;elseif numel(strings)>2,pad=.55;else,pad=.35;end
setappdata(ax,'R30LegendTopPadding',pad);
end
function caption(h,ax,letter,text,o)
pos=get(ax,'Position');y=.04;if pos(4)<.4,y=pos(2)-.13;end
annotation(h,'textbox',[pos(1)-.045 y pos(3)+.09 .07],...
 'String',[latin(letter,o) zh(text,o)],'Interpreter','tex','FontName',o.chineseFontName,...
 'FontSize',o.fontSize,'FontWeight','normal','EdgeColor','none','HorizontalAlignment','center','VerticalAlignment','middle');
end
function meshline(ax,f)
yl=ylim(ax);plot(ax,[f f],yl,':','Color',[.55 .55 .55],'LineWidth',.6,'HandleVisibility','off');
end
function name=finish(h,ax,out,name,o,limits)
for j=1:numel(ax)
 if iscell(limits),xlim(ax(j),limits{j});else,xlim(ax(j),limits);end
 if isappdata(ax(j),'R30LegendTopPadding'),yl=ylim(ax(j));ylim(ax(j),[yl(1),yl(2)+getappdata(ax(j),'R30LegendTopPadding')*diff(yl)]);end
end
savefig(h,fullfile(out,[name '.fig']));
try
 exportgraphics(h,fullfile(out,[name '.pdf']),'ContentType','vector','BackgroundColor','white');
 exportgraphics(h,fullfile(out,[name '.png']),'Resolution',o.pngDpi,'BackgroundColor','white');
catch
 pos=get(h,'Position');set(h,'PaperUnits','centimeters','PaperSize',pos(3:4),'PaperPosition',[0 0 pos(3:4)]);
 print(h,fullfile(out,[name '.pdf']),'-dpdf','-painters');print(h,fullfile(out,[name '.png']),'-dpng',sprintf('-r%d',o.pngDpi));
end
if ~o.visible,close(h);end
end
function name=choose_cjk_font(name)
try
 available=listfonts;
 if any(strcmpi(available,name)),return;end
 candidates={'SimSun','宋体','Songti SC','Noto Serif CJK SC','Noto Sans CJK SC'};
 for j=1:numel(candidates),if any(strcmpi(available,candidates{j})),name=candidates{j};return;end;end
 warning('R30F:ChineseFont','未检测到宋体或常用中文字体，请安装宋体或设置chineseFontName。');
catch
 % Font discovery is optional; MATLAB can resolve the requested font itself.
end
end
