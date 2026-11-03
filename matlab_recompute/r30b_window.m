function interval=r30b_window(P,o)
if ~isempty(o.window_s)
 interval=o.window_s(:)';assert(numel(interval)==2&&all(isfinite(interval))&&interval(1)>=0&&interval(2)>interval(1),...
  'R30B:Window','window_s must be [start end] in seconds');
else
 period=2*pi/abs(P.installation.omegaG);a=max(0,o.windowCenter_s-period/2);b=o.windowCenter_s+period/2;
 edges=linspace(a,b,ceil((b-a)/o.windowSearchStep_s)+1);possible=false(numel(edges)-1,1);
 for k=1:numel(possible)
  for f=1:size(o.faces,1)
   possible(k)=possible(k)||r30b_gate(P,o.faces(f,:),edges(k),edges(k+1),o.normalDepthLimit_mm);
  end
 end
 active=find(possible);assert(~isempty(active),'R30B:NoWindow','No enclosing-cylinder engagement found near the requested center');
 assert(active(1)>1&&active(end)<numel(possible),'R30B:WindowBoundary',...
  'Search interval cuts an enclosing engagement; change windowCenter_s or provide window_s');
 % Snap outwards to the history grid, with clear intervals on both sides.
 interval=[floor(edges(active(1)-1)/o.historyDt_s)*o.historyDt_s,...
  ceil(edges(active(end)+2)/o.historyDt_s)*o.historyDt_s];
end
steps=(interval(2)-interval(1))/o.dt_s;
if isempty(o.window_s),interval(2)=interval(1)+ceil(steps-1e-8)*o.dt_s;
else,assert(abs(steps-round(steps))<1e-7,'R30B:Window','Window must contain whole output intervals');end
end
