function [A,B,Q,Z] = swap_single(A,B)
% swap the poles of (A,B) via (AA, BB) = Q*(A,B)*Z

    A = flip(A,1);
    B = flip(B,1);

    if( abs(B(2,2)*A(1,1)) > abs(B(1,1)*A(2,2)) )

        H = B(2,2)*A - A(2,2)*B;

        [Z, ~] = ct_zero(H(1,2), H(1,1));
        Z(1,1) = conj(Z(1,1)); Z(2,2) = conj(Z(2,2));

        temp = B*Z(:,1);

        [Q, ~] = ct_zero(temp(1,1), temp(2,1));

        A = Q*(A*Z);
        B = (Q*B)*Z;

    else

        H = B(2,2)*A - A(2,2)*B;

        [Z, ~] = ct_zero(H(1,2), H(1,1));
        Z(1,1) = conj(Z(1,1)); Z(2,2) = conj(Z(2,2));

        temp = A*Z(:,1);

        [Q, ~] = ct_zero(temp(1,1), temp(2,1));

        A = Q*(A*Z);
        B = (Q*B)*Z;

    end
    
    A = flip(A,1);
    B = flip(B,1);
    Q = flip(flip(Q,1),2);

end

