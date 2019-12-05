function [A,Q, info] = type2e_pal(A)
%TYPE2E Performs a type2e move on the matrix A

info = {};
info.no_convergence = false;
info.refinement_steps = 0;

version = 2;
refine = 1;
normA = norm(A, 'fro');

Q = eye(3);

if(version == 1)
    z1 = [A(1,3), A(3,1)'];
    H = z1(2)*A - z1(1)*A';

    [Z1, ~] = ct_zero(H(2,3), H(2,2));
    Z1(1,1) = conj(Z1(1,1)); Z1(2,2) = conj(Z1(2,2));

    H(:,2:3) = H(:,2:3)*Z1;
    Q(:,2:3) = Q(:,2:3)*Z1;


    [Z2, ~] = ct_zero(H(3,2), H(3,1));
    Z2(1,1) = conj(Z2(1,1)); Z2(2,2) = conj(Z2(2,2));

    H(:,1:2) = H(:,1:2)*Z2;
    Q(:,1:2) = Q(:,1:2)*Z2;

    A = Q'*A*Q;
    [Z3, ~] = ct_zero(A(1,3), A(1,2));
    Z3(1,1) = conj(Z3(1,1)); Z3(2,2) = conj(Z3(2,2));
    A(:,2:3) = A(:,2:3)*Z3;
    A(2:3,:) = Z3'*A(2:3,:);
    Q(:,2:3) = Q(:,2:3)*Z3;
elseif(version == 2)
    
    x21 = - ( A(3,1)'*A(2,3) - A(1,3)*A(3,2)' ) / ( A(3,1)'*A(2,2) - A(1,3)*A(2,2)' );
    x32 = - ( A(2,2)'*A(3,2) - A(2,2)*A(2,3)' ) / ( A(2,2)'*A(3,1) - A(2,2)*A(1,3)' );
    t = [A(3,3);A(3,3)']+x21*[A(3,2);A(2,3)'];
    x31 = - ( A(3,1)'*t(1) - A(1,3)*t(2) ) / ( A(3,1)'*A(3,1) - A(1,3)*A(1,3)' );
    
    X = [
    1,0,0;
    x21,1,0;
    x31,x32,1
    ];

    [Q, ~] = qr(flip(X,1));
    A = Q'*A*Q;
    
else
    
end

if(refine == 1)
    tol = 10*eps;
    iter = 0;
    maxit = 10;
    while (abs(A(1,1)) > normA*tol || abs(A(2,1)) > normA*tol || abs(A(1,2)) > normA*tol )
        if(iter >= maxit)
%           warning('swap did not converge')
            info.no_convergence = true;
            break;
        end
        
        x31 = ( -A(1,1)'*A(3,1) + A(1,1)*A(1,3)' ) / ( A(3,1)'*A(3,1) - A(1,3)*A(1,3)' );
        y13 = ( -A(3,1)'*A(1,1) + A(1,3)*A(1,1)' ) / ( A(3,1)'*A(3,1) - A(1,3)*A(1,3)' );
        
        t = [-A(1,2)-y13*A(3,2);-A(2,1)'-y13*A(2,3)'];
        x32 = ( t(2)*A(2,2) - t(1)*A(2,2)' ) / ( A(3,1)'*A(2,2) - A(1,3)*A(2,2)' );
        
        t = [-A(2,1)-x31*A(2,3);-A(1,2)'-x31*A(3,2)'];
        x21 = ( t(2)*A(3,1) - t(1)*A(1,3)' ) / ( A(2,2)'*A(3,1) - A(2,2)*A(1,3)' );
        
        X = [
            1,0,0;
            x21,1,0;
            x31,x32,1
        ];

        [Q2, ~] = qr(X);
        A = Q2'*(A*Q2);
        Q = Q*Q2;

        iter = iter + 1;

    end
        
    info.refinement_steps = iter;
end



end

