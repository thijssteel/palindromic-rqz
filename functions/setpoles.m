function [A,Q, info] = setpoles(A, xi)
%SETPOLES Set half of the poles of the matrix A to xi

    n = size(A,1);
    Q = eye(n);
    assert(size(xi,1) < n/2)
    nb_swaps = 0;
    
    for xi_index=size(xi,1):-1:1
        
        %type 1 update
        T = xi(xi_index, 2)*A(n-1:n, 1:1+1) - xi(xi_index, 1)*A(1:1+1,n-1:n)';    
        [Qt, ~] = ct_zero(T(2,1),T(1,1));
        Qt = flip(flip(Qt, 2),1);

        A(n-1:n,:) = Qt*A(n-1:n,:);
        A(:,n-1:n) = A(:,n-1:n)*Qt';
        Q(:,n-1:n) = Q(:,n-1:n)*Qt';
        
        nb_swaps = nb_swaps + 1;
        
        for k=2:xi_index
            [~,~,Q1,Z1] = swap_single( A(n-k:n-k+1,k-1:k), A(k-1:k,n-k:n-k+1)' );
            Z2 = Q1';

            A(n-k:n,k-1:k) = A(n-k:n,k-1:k)*Z1;
            A(k-1:k,n-k:n) = Z1'*A(k-1:k,n-k:n);
            A(n-k:n-k+1,k-1:n) = Z2'*A(n-k:n-k+1,k-1:n);
            A(k-1:n,n-k:n-k+1) = A(k-1:n,n-k:n-k+1)*Z2;
            Q(:,k-1:k) = Q(:,k-1:k)*Z1;
            Q(:,n-k:n-k+1) = Q(:,n-k:n-k+1)*Z2;
            A(n-k,k-1) = 0;
            A(k-1,n-k) = 0;
        
            nb_swaps = nb_swaps + 1;
        end
        
        
    end
    
    info = {};
    info.nb_swaps = nb_swaps;
end