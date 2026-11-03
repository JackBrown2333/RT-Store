function half=r30b_worm_width(r,A)
% 理想渐开螺旋砂轮轴向半齿厚。当前参数不是实测修整型面。
rb=A.rPitch/hypot(1,A.rPitch*A.slope/A.helicalP);
assert(all(r(:)>=rb),'砂轮渐开螺旋面不能延伸到基圆柱以内');
u=sqrt(max((r/rb).^2-1,0));up=sqrt((A.rPitch/rb)^2-1);
half=A.lead/4-A.helicalP*((u-atan(u))-(up-atan(up)));
end
