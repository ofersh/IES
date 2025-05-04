clear all; close all;
global Hess c0 ConditionNum N LB UB loggerToFiles problemName fid mainindex
NRUNS=50; MaxITER=1000000;
problemName='IQP-Sphere'; fopt = 0;
LB=-1000; UB=1000; sigma = (UB-LB)/6;
DIM = [10,30,80]; %N = 80; 
loggerToFiles=1;

writeFiles=0; %epochHistory = MaxITER/1000;
%Hess = eye(N); ConditionNum = 1;
xi1=[7,-7]';
%mode='1p1DG'; mutFunc = 'geometric_mutation';
%mode='1p1TN'; mutFunc = 'normal_mutation';
%mode='1p1DU'; mutFunc = 'uniform_mutation';
mode='1p1SB'; mutFunc = 'binomial_mutation';
%
for N=DIM,
    c0 = repmat(xi1,N/2,1);
    Hess = eye(N); ConditionNum = 1;
    for j=1:N, DE{j}=[LB:UB]; end
    if loggerToFiles
        % filename = ['./output/',mode,'_',problemName,num2str(ConditionNum),'_DIM',num2str(N),'.csv'];
        filename = ['./output/',mode,'_',problemName,'_DIM',num2str(N),'.csv'];
        fid=fopen(filename,'wt');
        fprintf(fid,'EvalCounter, FuncValue, FuncID, AlgoID, DIM, RunID \n');
    end
    for mainindex=1:NRUNS
        [bX{mainindex},bF{mainindex},bT{mainindex},HF{mainindex},HS{mainindex}] = OnePlusOneIntegerES('problem', 'IntegerEllipsoid', 'N', N, 'sigma0', sigma, 'mutation', mutFunc,'mode',mode,'numIter', MaxITER, 'fstop', fopt, 'gap', 1e-3, 'DE', DE, 'sigmaMin', 1e-3);
        disp(['Best attained: ',num2str(bF{mainindex}),'; function evaluations used: ',num2str(bT{mainindex}),'.']);
    end
    disp('---===---');
    if loggerToFiles, fclose(fid); end
    if writeFiles
        X = cell2mat(bX); F = cell2mat(bF); T = cell2mat(bT); %hF = cell2mat(HF); hS = cell2mat(HS);
        basemat = ['bestG_',problemName,num2str(N),'_K',num2str(Kcnstr),'.mat'];
        save(basemat,'X', 'F', 'T', 'HF', 'HS', 'fopt');
    end
end
%-----------------------
