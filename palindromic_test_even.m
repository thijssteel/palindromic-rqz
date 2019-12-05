%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  Tests on random matrices of even sizes %%%%%%%%
%%%%% results are saved in 'results/pal_test_3.csv %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear;
n_samples = 1;
n_array = [100;200;400;800;1600];
n_experiments = length(n_array);

backward_errors_pal = zeros(n_experiments,1);
backward_errors_palz = zeros(n_experiments,1);

it_count_pal = zeros(n_experiments,1);
it_count_palz = zeros(n_experiments,1);

swap_count_pal = zeros(n_experiments,1);
swap_count_palz = zeros(n_experiments,1);
swap_count_pal_with_reduction = zeros(n_experiments,1);

mid_swap_count_pal = zeros(n_experiments,1);
mid_swap_count_palz = zeros(n_experiments,1);

refinement_count_pal = zeros(n_experiments,1);
refinement_count_palz = zeros(n_experiments,1);

for n_index=1:n_experiments
    n = n_array(n_index);
    maxit = 20*n;
    
    for k=1:n_samples

        A = 2*rand(n,n) + 1i*rand(n,n);
        for i=1:n
            for j=1:n-(i-1) - 2
                A(j,i) = 0;
            end
        end

        % test pal
        [AA, Q, e2, infoPal1] = palindromic_RQZ(A, maxit);

        backward_errors_palz(n_index) = max(backward_errors_palz(n_index), norm (Q'*(A*Q) - AA)/norm(A, 'fro'));
        it_count_palz(n_index) = it_count_palz(n_index) + infoPal1.nb_iterations;
        swap_count_palz(n_index) = swap_count_palz(n_index) + infoPal1.nb_swaps;
        mid_swap_count_palz(n_index) = mid_swap_count_palz(n_index) + infoPal1.nb_mid_swaps;
        refinement_count_palz(n_index) = refinement_count_palz(n_index) + infoPal1.nb_refinement;

        Xi = [ones( (n-2)/2,1 ), zeros((n-2)/2,1 )];
        [AAA, QQ, infoSetPoles] = setpoles(A,Xi);

        [AAAA, QQQ, e2, infoPal2] = palindromic_RQZ(AAA, maxit);

        backward_errors_pal(n_index) = max(norm (QQQ'*(QQ'*A*QQ)*QQQ - AAAA)/norm(A, 'fro'), backward_errors_pal(n_index));
        it_count_pal(n_index) = it_count_pal(n_index) + infoPal2.nb_iterations;
        swap_count_pal(n_index) = swap_count_pal(n_index) + infoPal2.nb_swaps;
        swap_count_pal_with_reduction(n_index) = swap_count_pal_with_reduction(n_index) + infoPal2.nb_swaps + infoSetPoles.nb_swaps;
        mid_swap_count_pal(n_index) = mid_swap_count_pal(n_index) + infoPal2.nb_mid_swaps;
        refinement_count_pal(n_index) = refinement_count_pal(n_index) + infoPal2.nb_refinement;
    
    end
    it_count_palz(n_index) = it_count_palz(n_index)./n_samples;
    it_count_pal(n_index) = it_count_pal(n_index)./n_samples;
    swap_count_palz(n_index) = swap_count_palz(n_index)./n_samples;
    swap_count_pal(n_index) = swap_count_pal(n_index)./n_samples;
    swap_count_pal_with_reduction(n_index) = swap_count_pal_with_reduction(n_index)./n_samples;
    mid_swap_count_palz(n_index) = mid_swap_count_palz(n_index)./n_samples;
    refinement_count_palz(n_index) = refinement_count_palz(n_index)./n_samples;
    mid_swap_count_pal(n_index) = mid_swap_count_pal(n_index)./n_samples;
    refinement_count_pal(n_index) = refinement_count_pal(n_index)./n_samples;
    
    table(n_array, backward_errors_pal, it_count_pal, swap_count_pal,...
        swap_count_pal_with_reduction, mid_swap_count_pal, refinement_count_pal,...
        backward_errors_palz, it_count_palz, swap_count_palz, mid_swap_count_pal,...
        refinement_count_pal)
end


writetable(table(n_array, backward_errors_pal, it_count_pal, swap_count_pal,...
        swap_count_pal_with_reduction, mid_swap_count_pal, refinement_count_pal,...
        backward_errors_palz, it_count_palz, swap_count_palz, mid_swap_count_pal,...
        refinement_count_pal), 'results/pal_test_3.csv')