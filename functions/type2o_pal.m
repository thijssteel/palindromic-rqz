function [A, Q, info] = type2o_pal(A)
%TYPE2O Performs a type2o move on the matrix A
    info = {};
    info.no_convergence = false;
    info.refinement_steps = 0;
    refine = 1;
    normA = norm(A,'fro');
    
    z = [A(1,2), A(2,1)'];
    H = z(2)*A(2,:) - z(1)*A(:,2)';

    [Q, ~] = ct_zero(H(2), H(1));
    Q(1,1) = conj(Q(1,1)); Q(2,2) = conj(Q(2,2));


    A = Q'*(A*Q);
    
    if(refine == 1)    
        tol = 10*eps;
        iter = 0;
        maxit = 10;
        while (abs(A(1,1)) > normA*tol)
            if(iter >= maxit)
                info.no_convergence = true;
                break;
            end
            
            H = z(2)*A - z(1)*A';

            [Q2, ~] = ct_zero(H(1,2), H(1,1));
            Q2(1,1) = conj(Q2(1,1)); Q2(2,2) = conj(Q2(2,2));

            A = Q2'*(A*Q2);
            Q = Q*Q2;

            iter = iter + 1;

        end
        
        info.refinement_steps = iter;
    end
    

end

