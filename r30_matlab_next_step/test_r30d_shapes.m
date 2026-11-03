function report=test_r30d_shapes(sourceM)
% Regression cases for MATLAB empty logical selections and one-point inputs.
% No toolboxes, material-history replay or graphics required.
if nargin<1,sourceM=fullfile(fileparts(mfilename('fullpath')),'example_data','multitooth_pilot.mat');end
D=load(sourceM,'M');A=D.M.parameters.installation;clear D;t0=15.014937746;
half=A.lead/4; % At the pitch radius the involute-width correction is zero.
r=[A.rPitch;A.tipNew+1;A.rPitch;A.rRoot-1;A.rPitch;A.rRoot-1];
ax=[0;0;A.axialLength/2+1;0;0;-A.axialLength/2-1];
phase=[half+.05;0;0;0;-half-.05;0];
phi=(ax-phase)/A.helicalP+A.omegaW*t0;
P=repmat([A.center 0 0],6,1)+[-r.*cos(phi),zeros(6,2)]+...
 (r.*sin(phi))*A.axisY(:)'+ax*A.axisW(:)';
expected=[1;2;3;4;1;3];G=r30d_wheel_geometry(P,t0,A);
assert(isequal(G.branch,expected),'R30D:TestBranch','Controlled geometry did not cover all branches.');
names={'level','branch','branchGap','gradient','outwardNormal','partialTimeDerivative','rigidTransportResidual'};
calls=0;maxDifference=0;
for orient=1:2
 B=A;if orient==1,B.axisW=A.axisW(:)';B.axisY=A.axisY(:)';else,B.axisW=A.axisW(:);B.axisY=A.axisY(:);end
 for j=1:6
  J=r30d_wheel_geometry(P(j,:),t0,B);calls=calls+1;
  for f=1:numel(names),jv=J.(names{f});gv=G.(names{f});maxDifference=max(maxDifference,max(abs(jv(:)-reshape(gv(j,:),[],1))));end
  assert(isequal(size(J.gradient),[1 3])&&isequal(size(J.level),[1 1]));
 end
 for i=1:6
  for j=1:6
   ids=[i j];
   for shape=1:3
    tt=t0;if shape==2,tt=[t0;t0];elseif shape==3,tt=[t0 t0];end
    J=r30d_wheel_geometry(P(ids,:),tt,B);calls=calls+1;
    for f=1:numel(names),jv=J.(names{f});gv=G.(names{f});v=jv-gv(ids,:);maxDifference=max(maxDifference,max(abs(v(:))));end
    assert(isequal(size(J.gradient),[2 3])&&isequal(size(J.level),[2 1]));
   end
  end
 end
 for shape=1:2
  tt=t0;if shape==2,tt=zeros(0,1);end
  E=r30d_wheel_geometry(zeros(0,3),tt,B);calls=calls+1;
  for f=1:numel(names)
   width=1;if any(strcmp(names{f},{'gradient','outwardNormal'})),width=3;end
   assert(isequal(size(E.(names{f})),[0 width]));
  end
 end
end
assert(maxDifference<1e-9);
epsX=1e-5;fd=zeros(size(P));
for j=1:3
 pp=P;pm=P;pp(:,j)=pp(:,j)+epsX;pm(:,j)=pm(:,j)-epsX;
 gp=r30d_wheel_geometry(pp,t0,A);gm=r30d_wheel_geometry(pm,t0,A);fd(:,j)=(gp.level-gm.level)/(2*epsX);
end
gradientError=max(abs(reshape(fd-G.gradient,[],1)));assert(gradientError<1e-5);
epsT=1e-7;gp=r30d_wheel_geometry(P,t0+epsT,A);gm=r30d_wheel_geometry(P,t0-epsT,A);
timeError=max(abs((gp.level-gm.level)/(2*epsT)-G.partialTimeDerivative));assert(timeError<1e-3);
transportError=max(abs(G.rigidTransportResidual));assert(transportError<1e-8);
assert(max(abs(sqrt(sum(G.outwardNormal.^2,2))-1))<1e-12);
report=struct('kind','empty_singleton_branch_regression_v3','all_four_geometry_branches',true,...
 'single_and_empty_inputs',true,'row_and_column_axes_and_times',true,'shape_test_calls',calls,...
 'batch_vs_single_max_difference',maxDifference,'gradient_finite_difference_error',gradientError,...
 'time_derivative_finite_difference_error',timeError,'rigid_transport_error',transportError);
disp(report);fprintf('R30方向诊断空分支/单点回归检查：PASS。\n');
end
