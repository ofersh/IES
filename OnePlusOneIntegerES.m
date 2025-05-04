function [xBest,fBest,tBest,history_f,history_sigma] = OnePlusOneIntegerES(varargin)
%(1+1) Integer Evolution Strategy 
%    numIter                 number of iterations
%    N                       dimension of a problem (number of variables)
%    problem                 problem name
%    DE                      constraints on parameters
%Default parameters:
SIGMA_MIN = 1.0;
numIter = 1000000;
N = 80;
fstop = 0;
problem = 'IntegerConstrainedEllipsoid';
gap = 1e-6; %epochHistory = 100;
for i=1:N, DE{i}=[0:1]; end
sigma = 1/6;
% Resolve parameters via varargin:
i=1;
while i<=length(varargin),
    argok = 1;
    switch varargin{i},
        case 'numIter', i=i+1; numIter = varargin{i};
        case 'N',       i=i+1; N = varargin{i};
        case 'fstop',    i=i+1; fstop = varargin{i};
        case 'gap',     i=i+1; gap = varargin{i};
        case 'problem', i=i+1; problem = varargin{i};
        case 'DE',      i=i+1; DE = varargin{i};
        case 'mutation', i=i+1; mutFunc = varargin{i};
        case 'mode', i=i+1; mode= varargin{i};
        case 'sigma0', i=i+1; sigma = varargin{i};
        case 'sigmaMin', i=i+1; SIGMA_MIN = varargin{i};
        case 'seed', i=i+1; xBest=varargin{i};
        otherwise argok=0;
    end
    if ~argok,
        disp(['(mies) Ignoring invalid argument #' num2str(i+1)]);
    end
    i = i+1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global loggerToFiles verbose problemName fid mainindex;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
epoch = 50; %numIter/10;
k_sigma = 0.827;
history_sigma = NaN(1,numIter); %[];
history_f = NaN(1,numIter);
xBest=init1p1IES(DE);
fBest = feval(problem,xBest);
tBest = NaN;
history_f(1) = fBest;
history_sigma(1) = sigma;
t = 2; osuccess = 0;
while (t <= numIter) && (fBest > fstop+gap)
    %
    X = feval(mutFunc,sigma,xBest,DE);
    F = feval(problem,X);
    if F < fBest
        fBest = F;
        xBest = X;
        osuccess = osuccess+1;
        tBest = t;
    end
    if mod(t,epoch) == 0
        disp(['XBest at iteration ' num2str(t),': ',num2str(fBest), ' [',num2str(X(1:N,1)'),']']);
        ps = osuccess/epoch;
        if (ps < 0.2)
            sigma = sigma * k_sigma;
        elseif (ps > 0.2)
            sigma = sigma / k_sigma;
        end
        sigma = max(sigma,SIGMA_MIN);
        osuccess = 0;
        if loggerToFiles, fprintf(fid,'%d, %f, %s, %s, %d, %d \n', t, fBest, problemName, mode, N, mainindex ); end
    end
    history_sigma(t) = sigma;
    history_f(t) = fBest;
    t=t+1;
end
disp(['XBest at iteration ' num2str(t),': ',num2str(fBest), ' [',num2str(X(1:N,1)'),']']);
if loggerToFiles, fprintf(fid,'%d, %f, %s, %s, %d, %d \n', t, fBest, problemName, mode, N, mainindex ); end
end

%--------------------------------------------------------------------------
function X = init1p1IES(DE)
N = length(DE);             %Number of variables
X = zeros(N ,1);            %Integers only

for j=1:N
    X(j) = round(rand()*(DE{j}(end)-DE{j}(1)) + DE{j}(1));
end
end
%------------------------   Mutation operators   ------------------------------
function Xnew = normal_mutation(sigma,X,~)
Xnew = X + round(sigma*randn(size(X)));
end
%
function Xnew = geometric_mutation(sigma,X,~)
N = length(X);
pGeo = 1 - (sigma/N) / (1 + sqrt(1 + (sigma/N)^2));
Geo_Step = floor(log(1-rand(N,1))/log(1-pGeo)) - floor(log(1-rand(N,1))/log(1-pGeo));
Xnew = X + Geo_Step;
end
%
function Xnew = uniform_mutation(sigma,X,~)
S = sigma/length(X);
MAX_RANGE = round(0.5*(2*S - 1 + sqrt(1+S^2)));
Xnew = X + randi([-MAX_RANGE MAX_RANGE],size(X));
end
%
function Xnew = binomial_mutation(S,X,~)
N = ceil((pi / 2) * (S / length(X))^2);
binomialZ = binornd(N, 0.5, size(X));
Xnew = X + (2*binomialZ - N);
end
%------------------------   Test functions   ------------------------------
function F = IntegerEllipsoid(X)
global Hess c0 ConditionNum
%disp(X)
xc1 = (X-c0)';
F = (xc1*Hess*xc1') / ConditionNum;
end
%
function F = IntegerConstrainedEllipsoid(X)
global Hess1 Hess2 c1 c2 Kcnstr penalty ConditionNum
%N = (size(X,1))/2;
%n = N/2;
xc1 = (X-c1)'; xc2 = (X-c2)';
F = (xc1*Hess1*xc1')/ConditionNum;
G = (xc2*Hess2*xc2')/ConditionNum;
if G > Kcnstr
    F = F + penalty*(G - Kcnstr)^2;
end
end
function F = MixedVarsConstrainedEllipsoid(X)
global Hess1 Hess2 c1 c2 c3 c4 Kcnstr penalty ConditionNum
N = size(X,1);
n = N/2;
F = nan; G = nan;
sepF = length(Hess1)==n;
sepG = length(Hess2)==n;
if sepF %Raw objective function
    xc1 = (X(1:n)-c1)'; zc1 = (X(n+1:N)-c1)';
    F = (xc1*Hess1*xc1' + zc1*Hess1*zc1')/ConditionNum;
else
    xc1 = (X(1:N)-c3)';
    F = (xc1*Hess1*xc1')/ConditionNum;
end
if sepG %Parabolic constraint
    xc2 = (X(1:n)-c2)'; zc2 = (X(n+1:N)-c2)';
    G = (xc2*Hess2*xc2' + zc2*Hess2*zc2')/ConditionNum;
else
    xc2 = (X(1:N)-c4)';
    G = (xc2*Hess2*xc2')/ConditionNum;
end
if G > Kcnstr
    F = F + penalty*(G - Kcnstr)^2;
end
end