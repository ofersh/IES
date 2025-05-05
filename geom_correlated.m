function z = geom_correlated(alpha_matrix)
n = length(alpha_matrix);
S = diag(alpha_matrix);
sfix=sqrt(pi/2); %via a continuous approximation to the expected L1 norm
% Generate uncorrelated double-geometric variates
su = NaN(n,1); %S .* randn(n,1);
for j=1:n
    p = 1.0 - ( (S(j)/sfix) / (1.0+sqrt(1+((S(j)/sfix)^2))) );
    su(j) = geornd(p);%-geornd(p);
end
% Apply rotations to introduce correlations
nq = n*(n-1)/2;
for k=1:n-1,
    n1 = n-k;
    n2 = n;
    for i=1:k,
        d1 = su(n1);
        d2 = su(n2);
        su(n1) = (d1*cos(alpha_matrix(n1,n2))- d2*sin(alpha_matrix(n1,n2)));
        su(n2) = (d1*sin(alpha_matrix(n1,n2))+ d2*cos(alpha_matrix(n1,n2)));
        n2 = n2-1;
        nq = nq-1;
    end
end
z = su; %continuous - prior to rounding
end