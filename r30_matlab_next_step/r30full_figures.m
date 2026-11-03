function r30full_figures(outDir)
% Paper-style MATLAB figures; raw extrema are retained in overview bins.
D=load(fullfile(outDir,'full_run_result.mat'),'R','overview','stats');R=D.R;O=D.overview;S=D.stats;
figDir=fullfile(outDir,'figures');if ~isfolder(figDir),mkdir(figDir);end
CN='\fontname{SimSun}';EN='\fontname{Times New Roman}';blue=[0 .145 1];red=[1 .11 0];
h=newfig();
for j=1:2
 ax=panel(h,j);x=O(:,1);
 if j==1,y=O(:,2);lo=O(:,3);hi=O(:,4);yl=[CN '合力' EN '(N)'];cap='(a)全程合力及区间极值';
 else,y=O(:,7);lo=O(:,8);hi=O(:,9);yl=[CN '转矩' EN '(N\cdot m)'];cap='(b)全程C轴载荷转矩';end
 fill(ax,[x;flipud(x)],[lo;flipud(hi)],[.83 .87 1],'EdgeColor','none');plot(ax,x,y,'Color',blue,'LineWidth',.7);
 xlabel(ax,[CN '模型时间' EN '(s)'],'Interpreter','tex');ylabel(ax,yl,'Interpreter','tex');caption(h,j,cap);
 legend(ax,{[CN '区间最小至最大值'],[CN '区间均值']},'Interpreter','tex','Box','off','FontSize',8,'Location','best');
end
savepaper(h,figDir,'full_force_overview');
D=load(fullfile(outDir,'peak_window.mat'),'peakWindow');W=D.peakWindow;F=W.globalData(:,3:5);t=(W.modelTime_s-W.modelTime_s(1))*1000;
h=newfig();ax=panel(h,1);plot(ax,t,vecnorm(F,2,2),'Color',blue,'LineWidth',.7);plot(ax,t,vecnorm(W.globalData(:,6:8),2,2),'k','LineWidth',.7);
xlabel(ax,[CN '局部时间' EN '(ms)'],'Interpreter','tex');ylabel(ax,[CN '力' EN '(N)'],'Interpreter','tex');
legend(ax,{[CN '全齿合力'],[CN '纯法向向量合力']},'Interpreter','tex','Box','off','FontSize',8,'Location','best');caption(h,1,'(a)全程合力峰值附近');
loads=zeros(numel(W.toothIds),1);for g=1:numel(loads),loads(g)=max(vecnorm(reshape(W.perToothForce_N(:,g,:),numel(t),3),2,2));end
[~,rank]=sort(loads,'descend');a=rank(1);b=rank(min(2,numel(rank)));Fa=reshape(W.perToothForce_N(:,a,:),numel(t),3);Fb=reshape(W.perToothForce_N(:,b,:),numel(t),3);
ax=panel(h,2);plot(ax,t,vecnorm(Fa,2,2),'Color',blue,'LineWidth',.7);
if numel(rank)>=2
 plot(ax,t,vecnorm(Fb,2,2),'Color',red,'LineWidth',.7);plot(ax,t,vecnorm(Fa+Fb,2,2),'k','LineWidth',.7);
 legend(ax,{sprintf('%s第%d齿',CN,W.toothIds(a)),sprintf('%s第%d齿',CN,W.toothIds(b)),[CN '两齿矢量叠加']},'Interpreter','tex','Box','off','FontSize',8,'Location','best');
end
xlabel(ax,[CN '局部时间' EN '(ms)'],'Interpreter','tex');ylabel(ax,[CN '合力' EN '(N)'],'Interpreter','tex');caption(h,2,'(b)主要受力齿与叠加');savepaper(h,figDir,'peak_multitooth_forces');
h=newfig();ax=panel(h,1);stairs(ax,O(:,1),O(:,10),'Color',blue,'LineWidth',.7);stairs(ax,O(:,1),O(:,11),'Color',red,'LineWidth',.7);
xlabel(ax,[CN '模型时间' EN '(s)'],'Interpreter','tex');ylabel(ax,[CN '数量'],'Interpreter','tex');
legend(ax,{[CN '去除齿数'],[CN '去除齿面数']},'Interpreter','tex','Box','off','FontSize',8,'Location','best');caption(h,1,'(a)每个显示区间的最大参与数');
ax=panel(h,2);plot(ax,S(1:2:end,1),S(1:2:end,5),'Color',blue,'LineWidth',.7);plot(ax,S(2:2:end,1),S(2:2:end,5),'Color',red,'LineWidth',.7);
xlabel(ax,[CN '齿号'],'Interpreter','tex');ylabel(ax,[CN '最大剩余余量' EN '(\mum)'],'Interpreter','tex');
legend(ax,{[CN '侧1'],[CN '侧2']},'Interpreter','tex','Box','off','FontSize',8,'Location','best');caption(h,2,'(b)各齿两侧最终最大余量');savepaper(h,figDir,'contact_and_remaining_stock');
D=load(fullfile(outDir,'central_segment_fft.mat'),'spectrum');SP=D.spectrum;
if SP.computed
 h=newfig();take=SP.frequency_Hz<=1000;
 for j=1:2
  ax=panel(h,j);y=SP.forceNormAmplitude_N;yl=[CN '合力幅值' EN '(N)'];cap='(a)合力模的单边幅值谱';
  if j==2,y=SP.torqueAmplitude_Nm;yl=[CN '转矩幅值' EN '(N\cdot m)'];cap='(b)C轴转矩的单边幅值谱';end
  plot(ax,SP.frequency_Hz(take),y(take),'Color',blue,'LineWidth',.65);
  xlabel(ax,[CN '频率' EN '(Hz)'],'Interpreter','tex');ylabel(ax,yl,'Interpreter','tex');caption(h,j,cap);
 end
 savepaper(h,figDir,'central_segment_fft');
end
% Signed height maps for every calculated tooth; zero is the target flank.
cmap=[linspace(0,1,128)' linspace(.2,1,128)' ones(128,1);ones(128,1) linspace(1,0,128)' linspace(1,0,128)'];
for tooth=R.toothIds
 D=load(fullfile(outDir,sprintf('tooth_%02d',tooth),'final_state.mat'),'finalState');V=D.finalState;
 limit=max([1;abs(V.faces{1}.height_mm(:))*1000;abs(V.faces{2}.height_mm(:))*1000]);h=newfig();
 for j=1:2
  ax=panel(h,j);v=V.faces{j};
  % The profile mesh is quadratic, not uniform: draw actual weighted supports.
  [Z,SS]=meshgrid(v.support_z_mm,v.support_s_mm);C=zeros(size(Z));
  C(1:end-1,1:end-1)=v.height_mm*1000;C(end,:)=C(end-1,:);C(:,end)=C(:,end-1);
  surf(ax,Z,SS,zeros(size(Z)),C,'FaceColor','flat','EdgeColor','none');view(ax,2);axis(ax,'tight');
  set(ax,'YDir','normal');caxis(ax,[-limit limit]);colormap(ax,cmap);
  xlabel(ax,[CN '齿宽坐标' EN '(mm)'],'Interpreter','tex');ylabel(ax,[CN '齿廓弧长' EN '(mm)'],'Interpreter','tex');
  cb=colorbar(ax);cb.FontName='Times New Roman';cb.FontSize=8;cb.Label.String=[CN '余量' EN '(\mum)'];cb.Label.Interpreter='tex';
  caption(h,j,sprintf('(%c)第%d齿侧%d',96+j,tooth,j));
 end
 savepaper(h,figDir,sprintf('final_tooth_%02d',tooth));
end
fprintf('MATLAB图已保存：%s\n',figDir);
end
function h=newfig()
h=figure('Visible','off','Color','w','Units','centimeters','Position',[2 2 16 6.8],'PaperPositionMode','auto');
end
function ax=panel(h,j)
ax=axes('Parent',h,'Position',[.09+(j-1)*.51 .30 .35 .64],'FontName','Times New Roman','FontSize',10,'Box','on','TickDir','in','LineWidth',.8);hold(ax,'on');
end
function caption(h,j,s)
annotation(h,'textbox',[.045+(j-1)*.51 .035 .44 .08],'String',s,'FontName','SimSun','FontSize',10.5,'EdgeColor','none','HorizontalAlignment','center');
end
function savepaper(h,folder,name)
path=fullfile(folder,name);savefig(h,[path '.fig']);
try,exportgraphics(h,[path '.pdf'],'ContentType','vector','BackgroundColor','white');exportgraphics(h,[path '.png'],'Resolution',400,'BackgroundColor','white');
catch,set(h,'PaperUnits','centimeters','PaperSize',[16 6.8],'PaperPosition',[0 0 16 6.8]);print(h,[path '.pdf'],'-dpdf','-painters');print(h,[path '.png'],'-dpng','-r400');end
close(h);
end
