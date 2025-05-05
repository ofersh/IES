function z = corr_sigma(alpha_matrix)
n = length(alpha_matrix);
S = diag(alpha_matrix);

% Generate uncorrelated standard normal variates
su = S .* randn(n,1);
%alpha = triu(alpha_matrix,1);

% Apply rotations to introduce correlations
nq = n*(n-1)/2;
for k=1:n-1,
    n1 = n-k;
    n2 = n;
    for i=1:k,
        %disp(['nq = ',num2str(nq), '; n1 = ', num2str(n1), '; n2 = ', num2str(n2)]);
        d1 = su(n1);
        d2 = su(n2);
        su(n1) = d1*cos(alpha_matrix(n1,n2))- d2*sin(alpha_matrix(n1,n2));
        su(n2) = d1*sin(alpha_matrix(n1,n2))+ d2*cos(alpha_matrix(n1,n2));
        n2 = n2-1;
        nq = nq-1;
    end
end
z = su;
end