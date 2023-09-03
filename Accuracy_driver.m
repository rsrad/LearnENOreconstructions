% To check scheme accuracy test for smooth initial condition  
% Define some parameters
close all
clear all
global choice_of_r non_osc_method
format short e
param.test=1; %  from initial condition mat file 
modelchoice =[1,2,3,4,5];
param.model = modelchoice(4);
%global choice_of_r
choice_of_r =1;
hidden_layer_n =7;

% 1= Ms = Model trained on smooth data set
% 2= Md = Model trained on discontinuos data set
% 3 = Mu = union model
% 4 = ML = flux limiter
% 5 = Mh = hybrid model 
param.rk_method='SSP3'; % ODE solver
param.bc='Periodic';% Boundary condition
param.gc=4; % Number of ghost cell on each side

param.animation=0;
ntable=5;
non_osc_method= 1;  %1: classical, 2: length based eno, 3: non_osc_method
if non_osc_method ==1
    method= 'eno3';
    data_dir = '../ENO3data/Details/';
    prename ='NNENO-';
    param.gc=4;% no. of ghost cell
    param.reconstruction= @ NN_eno3; %NN_eno; % options: NN_eno, enol3, weno5, wenol5,
    figno=5;
elseif non_osc_method ==2
    method= 'enol3';
    data_dir = '../ENOL3data/Details/';
    prename ='NNENO-';
    param.gc=4;% no. of ghost cell
    param.reconstruction= @ NN_eno3L; %NN_eno; % options: NN_eno, enol3, weno5, wenol5,
    figno=6;
else
    method= 'weno3c';
    data_dir = '../WENO3data/Details/';
    prename ='NNWENO-';
    param.gc=2;% no. of ghost cell
    param.reconstruction= @ NN_weno;%NN_weno; %NN_eno; % options: NN_eno, enol3, weno5, wenol5,
    figno=6;
end
smooth_weight_file = data_dir + "trained_weights_"+"Hn_"+num2str(hidden_layer_n) + "smoothdata_"+choice_of_r+".mat";
smooth_bias_file = data_dir +"trained_biases_"+"Hn_"+num2str(hidden_layer_n) + "smoothdata_" + choice_of_r+ ".mat";

nonsmooth_weight_file = data_dir +"trained_weights_"+"Hn_"+num2str(hidden_layer_n) + "nonsmoothdata_"+choice_of_r+".mat";
nonsmooth_bias_file = data_dir +"trained_biases_"+"Hn_"+num2str(hidden_layer_n) + "nonsmoothdata_" + choice_of_r+ ".mat";

uniondata_weight_file = data_dir +"trained_weights_"+"Hn_"+num2str(hidden_layer_n) + "union_data_"+choice_of_r+".mat";
uniondata_bias_file = data_dir +"trained_biases_"+"Hn_"+num2str(hidden_layer_n) + "union_data_" + choice_of_r+ ".mat";

% smooth weights and biases
 load(smooth_weight_file,'WEIGHT_smoothdata')
 load(smooth_bias_file,'BIAS_smoothdata')
 param.weight=WEIGHT_smoothdata;
 param.bias=BIAS_smoothdata;

% non_smooth weights and biases
 load(nonsmooth_weight_file,'WEIGHT_nonsmoothdata')
 load(nonsmooth_bias_file,'BIAS_nonsmoothdata')
 param.weight1=WEIGHT_nonsmoothdata; % weight1 and bias1 for non smooth
 param.bias1=BIAS_nonsmoothdata; 

% union(smooth+discontnuous) weights and biases
 load(uniondata_weight_file,'WEIGHT_union_data')
 load(uniondata_bias_file,'BIAS_union_data')
 param.weight2= WEIGHT_union_data; % weight1 and bias1 for smooth
 param.bias2= BIAS_union_data;
 
[L1_error,Linf_error,L1_order,Linf_order]=deal(zeros(1,ntable)); % pre-allocate memory
disp('% accuracy table')
for k=1:ntable
    param.nx= 20*2^(k-1); % Number of grids
    [t,u,param]=scalar1D_solver(param); % Main solver
%     Exact=sin(2*pi*(param.x-t));
    P_exact=param.exact_sol(param.x,t);
% P_exact =Exact;
    u_exact=P_exact;
    L1_error(k)=(1/param.nx)*vecnorm(u_exact-u,1);
    Linf_error(k)=vecnorm(u_exact-u,Inf);
    if k==1
        L1_order(k)=0.0;
        Linf_order(k)=0.0;
        disp('%% Exporting latex table format of accuracy test')
        disp('\begin{table}[htb!]')
        disp('\centering')
        disp('\begin{tabular}{|c|c|c|c|c|}')
        disp('\hline ')
        fprintf(' N  & $L^\\infty$ error  &  Order &    $L^1$ error    & Order \\\\ \n')
        disp('\hline ')
        fprintf('%% N & Linf err & Order &  L1 err  & Order \\\\ \n')
        fprintf(' %d & %1.2E &  ...  & %1.2E & ...   \\\\ \n', param.nx, Linf_error(k), L1_error(k))
    else
        L1_order(k)=log2(L1_error(k-1)/L1_error(k));
        Linf_order(k)=log2(Linf_error(k-1)/Linf_error(k));
        fprintf(' %d & %1.2E &  %.2f & %1.2E & %.2f  \\\\ \n',param.nx,Linf_error(k), Linf_order(k), L1_error(k), L1_order(k))
    end
end
disp('\hline')
disp('\end{tabular} ')
disp('\caption{}')
disp('\label{tab:}')
fprintf('%% The solution of test %d at T=%.2f for CFL=%.2f. \n',param.test,param.tf,param.cfl)
disp(['% It was exucuted through ',cd,' on ','date'])
disp('\end{table}')