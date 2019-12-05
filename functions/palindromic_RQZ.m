function [A,Q,e,info] = palindromic_RQZ(A, maxit)
%PALINDROMIC_RQZ Solves the palindromic eigenvalue problem
%  Input:  - A:      nxn matrix, assumed to be in anti-hessenberg form
%                    without checking
%          - maxit:  maximum number of iterations (I recommend about 2.7 iterations per eigenvalue as a minimum)
%  Output: - A:      A = Q'*A_input*Q, in anti-triangular form

tol = 3*eps;
n = size(A,1);

indexes = [];
res = [];
e = [];

istart = 1;
iend = n;

Q = eye(n);

it_count = 0;
nb_swaps = 0;
nb_mid_swaps = 0;
nb_refinement = 0;
last_deflation = 0;
while (iend - istart > 3 && it_count < maxit)
    
    %check convergence
    if(abs(A(istart,iend-1)) < tol*(abs(A(istart,iend)) + abs(A(istart+1,iend-1))) && abs(A(iend-1,istart)) < tol*( abs(A(iend-1,istart+1)) + abs(A(iend,istart)) ) )
        A(istart,iend-1) = 0;
        A(iend-1,istart) = 0;
        istart = istart+1;
        iend = iend-1;
        last_deflation = it_count;
        e = [e; A(istart,iend)./A(iend,istart)'; A(iend,istart)./A(istart,iend)'];
    elseif( abs(A(istart+1,iend-2)) < tol*(abs(A(istart+1,iend-1)) + abs(A(istart+2,iend-2))) && abs(A(iend-2,istart+1)) < tol*( abs(A(iend-2,istart+2)) + abs(A(iend-1,istart+1)) ) )
        %TODO: solve 2x2 system before deflating
        istart = istart+2;
        iend = iend-2;
        last_deflation = it_count;
        e = [e; eig(A(iend-1:iend,istart:istart+1), A(istart:istart+1, iend-1:iend)'); ...
            eig(A(istart:istart+1, iend-1:iend), A(iend-1:iend,istart:istart+1)' )];
    end
    
    if(iend - istart <= 3)
        break;
    end

    %store some info
%     indexes = [indexes;istart,iend]; %#ok<AGROW>
%     res = [res, zeros(n-1,1)]; %#ok<AGROW>
%     for j=1:n-1
%         res(j,it_count + 1) = abs(A(n-j,j))/(abs(A(n-j+1,j)) + abs(A(n-j,j+1))) + abs(A(j,n-j))/(abs(A(j+1,n-j)) + abs(A(j,n-j+1)));
%     end

    
    %determine shifts as eigenvalues of top right 2x2 block
    s = eig(A(istart:istart+1,iend-1:iend), A(iend-1:iend, istart:istart+1)');
    %take eigenvalue closest to rightmost entry
    [~,idx] = min(abs(s - A(istart,iend)/A(iend,istart)'));
    z = [s(idx), 1];
    if(abs(z(1)) > 10^15)
        z = [sign(z(1)),0];
    end
    if(abs(z(1)) < 10^-15)
        z = [0,sign(z(1))];
    end
    if(it_count - last_deflation >= 6 && mod(it_count - last_deflation, 6) == 0 )
        if(iend - istart < 15)
            %small matrix with slow convergence, check if eigenvalues are
            %one
            e2 = eig(A(istart:iend,istart:iend),A(istart:iend,istart:iend)');
            if(norm(abs(e2)-1) < 10^-13)
                break;
            end
        end
        %slow convergence detected, try a random shift
        z = 0.001*[rand(1) + 1i*rand(1), 1];
    end
    
    %calculate type 1 update
    T = z(2)*A(iend-1:iend, istart:istart+1) - z(1)*A(istart:istart+1,iend-1:iend)';    
    [Qt, ~] = ct_zero(T(2,1),T(1,1));
    Qt = flip(flip(Qt, 2),1);

    %apply type 1 update
    A(iend-1:iend,n-iend+1:n) = Qt*A(iend-1:iend,n-iend+1:n);
    A(n-iend+1:n,iend-1:iend) = A(n-iend+1:n,iend-1:iend)*Qt';
    Q(:,iend-1:iend) = Q(:,iend-1:iend)*Qt';
    
    nb_swaps = nb_swaps + 1;
    
    accuracy = abs(z(1) - A(iend-1,istart)/A(istart,iend-1)');
    if(accuracy > 10^-10)
%         warning('inaccurate swap, %d', accuracy)
    end
    
    
    mid_swap_done = false;
    swap_failed = false;
    for k=istart+1:iend-2
        if(swap_failed)
            break;
        end
        if(k < n/2 || k > (n/2 + 1) )
            % swap shifts k and k-1
            
            %check for deflation
            if(abs(A(n-k,k)) < tol*(abs(A(n-k+1,k)) + abs(A(n-k,k+1))) ...
                    && abs(A(k,n-k)) < tol*(abs(A(k+1,n-k)) + abs(A(k,n-k+1))) )
                %shift k is undefined, just skip this iteration, we will
                %reintroduce the shift next iteration
                continue;
            end
            if(abs(A(n-k+1,k-1)) < tol*(abs(A(n-k+1,k)) + abs(A(n-k+2,k-1))) ...
                && abs(A(k-1,n-k+1)) < tol*(abs(A(k,n-k+1)) + abs(A(k-1,n-k+2))) )
                %shift k-1 is undefined, we can now use a type 1
                %transformation to reintroduce the shift
                
                istart2 = k;
                iend2 = n-k+1;
                
                T = z(2)*A(iend2-1:iend2, istart2:istart2+1) - z(1)*A(istart2:istart2+1,iend2-1:iend2)';    
                [Qt, ~] = ct_zero(T(2,1),T(1,1));
                Qt = flip(flip(Qt, 2),1);
                A(iend2-1:iend2,n-iend2+1:n) = Qt*A(iend2-1:iend2,n-iend2+1:n);
                A(n-iend2+1:n,iend2-1:iend2) = A(n-iend2+1:n,iend2-1:iend2)*Qt';
                Q(:,iend2-1:iend2) = Q(:,iend2-1:iend2)*Qt';
                
                continue;
            end
            
            %the actual swapping
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
            
            
            accuracy = abs(z(1) - A(n-k,k)/A(k,n-k)');
            if(accuracy > 10^-10)
%                 warning('inaccurate swap, %d', accuracy)
            end
            
            nb_swaps = nb_swaps + 1;
            
            %check for deflation
            if(abs(A(n-k+1,k-1)) < tol*(abs(A(n-k+1,k)) + abs(A(n-k+2,k-1))) ...
                && abs(A(k-1,n-k+1)) < tol*(abs(A(k,n-k+1)) + abs(A(k-1,n-k+2))) )
                %shift k-1 is undefined, we can now use a type 1
                %transformation to reintroduce the shift
                
                istart2 = k;
                iend2 = n-k+1;
                
                T = z(2)*A(iend2-1:iend2, istart2:istart2+1) - z(1)*A(istart2:istart2+1,iend2-1:iend2)';    
                [Qt, ~] = ct_zero(T(2,1),T(1,1));
                Qt = flip(flip(Qt, 2),1);
                A(iend2-1:iend2,n-iend2+1:n) = Qt*A(iend2-1:iend2,n-iend2+1:n);
                A(n-iend2+1:n,iend2-1:iend2) = A(n-iend2+1:n,iend2-1:iend2)*Qt';
                Q(:,iend2-1:iend2) = Q(:,iend2-1:iend2)*Qt';
                
                continue;
            end
            
            
        else
            if(~mid_swap_done && ~swap_failed)
                %special swap in the middle
                if(mod(n,2) == 1)
                    [~,Qt,info] = type2o_pal(A(k-1:k,k-1:k));

                    if(info.no_convergence)
                        swap_failed = true;
                        continue;
                    end

                    A(k-1:n,k-1:k) = A(k-1:n,k-1:k)*Qt;
                    A(k-1:k,k-1:n) = Qt'*A(k-1:k,k-1:n);
                    Q(:,k-1:k) = Q(:,k-1:k)*Qt;
                    A(k-1,k-1) = 0;
                else
                    [At,Qt,info] = type2e_pal(A(k-1:k+1,k-1:k+1));

                    if(info.no_convergence)
                        swap_failed = true;
                        continue;
                    end

                    A(k-1:n,k-1:k+1) = A(k-1:n,k-1:k+1)*Qt;
                    A(k-1:k+1,k-1:n) = Qt'*A(k-1:k+1,k-1:n);
                    Q(:,k-1:k+1) = Q(:,k-1:k+1)*Qt;
                    A(k-1,k-1) = 0;
                    A(k-1,k) = 0;
                    A(k,k-1) = 0;
                end
                if(info.refinement_steps > 0)
                    nb_refinement = nb_refinement + 1;
                end
                nb_mid_swaps = nb_mid_swaps + 1;
                nb_swaps = nb_swaps + 1;
                mid_swap_done = true;
            end
        end

    end
    it_count = it_count + 1;

end

info = {};
info.indexes = indexes;
info.nb_swaps = nb_swaps;
info.nb_iterations = it_count;
info.istart = istart;
info.iend = iend;
info.residuals = res;
info.nb_mid_swaps = nb_mid_swaps;
info.nb_refinement = nb_refinement;



end

