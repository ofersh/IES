function A_psd = enforceSymmetryAndPSD(A)
    % Step 1: Enforce symmetry
    A_sym = (A + A') / 2;  % Make the matrix symmetric
    
    % Step 2: Perform Eigenvalue Decomposition
    [Q, D] = eig(A_sym);  % Q contains eigenvectors, D is diagonal matrix of eigenvalues
    
    % Step 3: Set negative eigenvalues to zero
    %D(D <= 0) = 1e-3;  % Replace negative eigenvalues with 1e-3
    for i=1:length(D), if D(i,i)<=0, D(i,i)=1e-3; end, end
    % Step 4: Reconstruct the matrix with non-negative eigenvalues
    A_psd = Q * D * Q';
    
    % Ensure symmetry again to mitigate numerical errors
    A_psd = (A_psd + A_psd') / 2;  % Double check symmetry
    
    % Step 5: Output the positive semi-definite matrix
end