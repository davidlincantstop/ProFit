%this program is about prior knowledge fitting: a direct implementation of 
%profit algorithm introduced in "Profit: two-dimensional prior-knowledge
%fitting of J-resolved spectra by RF Schulte - 2006
%The algorithm executes three iterations of a decided fitting procedure.
%each iteration involves a first part of non-linear least-squares fitting
%with constraints where the concentrations parameters are substituted which
%directly reduced the number of fitting parameters by the number of
%metabolite concentrations to be fitted.The results will be used to
%update the initial values of the parameters
%A second part of the fitting uses a linear least-square optimisation
%method,for the fitting of the concentration parameters . the results will
%go on to update the basis data.
%The freedom of the fitting will be increased after each iteration.
%The first two iterations only include the three predominant
%singlets of NAA ,creatine and choline while the difference being the
%second iteration will include the lineshape distortion function.
%The last iteration will be fitting with all the metabolites.

%% initialization
clear all;
close all;
clc;
%%

% first iteration fit 
% x0=[delta1,phi0,theta1,theta2,delta2,sigmae,sigmag]

% non-linear parameters initialisation
gs = GlobalSearch('TolX',0.1,'MaxTime',5000,'Display','iter','NumStageOnePoints',1000,'NumTrialPoints',3000,'PlotFcns',@gsplotbestf,'StartPointsToRun','bounds');
x0=[zeros(1,9)];
lb(1:9,1)=[-1;-pi/90;-1;-1;-1;0;0;0;0];
ub(1:9,1)=[1;pi/90;1;1;1;1;1;1;1];
options=optimset('Algorithm','interior-point');
% xstart = randn(3,1);
% nonliear fitting
problem = createOptimProblem('fmincon', ...
    'x0',x0,'objective',@fcost, ...
    'lb',lb,'ub',ub,'options',options);





[x , ~ , ~, output] = run(gs,problem);

% linear least-squares fitting
lbc(1:3,1)=0.8;ubc(1:3,1)=1.2;
c0=ones(3,1);
[~,Bt_r,sr]=fcost(x);
c = lsqlin(Bt_r,sr,[],[],[],[],lbc,ubc,c0);



%%
%second interation
% non-linear fitting
gs = GlobalSearch('TolX',0.1,'MaxTime',5000,'Display','iter','NumStageOnePoints',1000,'NumTrialPoints',3000,'PlotFcns',@gsplotbestf,'StartPointsToRun','bounds');

lb(1:11,1)=[-1;-pi/90;-pi/90;-pi/90;-1;-1;-1;0;0;0;0];
ub(1:11,1)=[1;pi/90;pi/90;pi/90;1;1;1;1;1;1;1];


% lb(1:7,1)=-inf;lb(8:11,1)=0;
% ub(1:11,1)=inf;
x=[x(1:2) 0 0 x(3:9)];
problem = createOptimProblem('fmincon', ...
    'x0',x,'objective',@fcost, ...
    'lb',lb,'ub',ub,'options',options);

[x , ~ , ~, output] = run(gs,problem);


% x = fmincon(@fcost,x,[],[],[],[],lb',ub',[],options);
%linear fitting
[~,Bt_r,sr]=fcost(x);
c = lsqlin(Bt_r,sr,[],[],[],[],lbc,ubc,c);

%%
% third iteration fitting
% non-linear fitting
gs = GlobalSearch('TolX',0.1,'MaxTime',10000,'Display','iter','NumStageOnePoints',1000,'NumTrialPoints',3000,'PlotFcns',@gsplotbestf,'StartPointsToRun','bounds');

lb(1:31,1)=[-1;-pi/90;-pi/90;-pi/90;-ones(13,1);zeros(14,1)];
ub(1:31,1)=[1;pi/90;pi/90;pi/90;ones(27,1)];

x=[x(1:7) zeros(1,10) x(8:10) zeros(1,10) x(11)];
% x=[ zeros(1,21)  ones(1,12) ];
% options=optimset('Algorithm','interior-point');
problem = createOptimProblem('fmincon', ...
    'x0',x,'objective',@fcost, ...
    'lb',lb,'ub',ub,'options',options);
[x fval eflag output] = run(gs,problem);


% x = fmincon(@fcost,x,[],[],[],[],lb',ub',[],options);
%linear fitting
[~,Bt_r,sr]=fcost(x);
%c=[NAAS creatine303 tcholine   Creatine Gln Glu  NAAG PCr 
%   PCho Tau mI];
lbc(1:13,1)=0.8;ubc(1:13,1)=1.2;
c = lsqlin(Bt_r,sr,[],[],[],[],lbc,ubc,[c ;ones(10,1)]);

