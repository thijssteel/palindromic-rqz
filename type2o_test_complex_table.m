%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Tests of type IIo moves on complex matrices %%%%%%%%%%%%%%%%%%%%
%%%%% results are saved in results/type2o_refinement_complex.csv' %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
buckets = [-15, -12;
            -12,-9;
            -9, 0;
            0, 3
            ];
n_buckets = size(buckets, 1);
n_experiments = 100000; %the number of experiments per bucket


errors = zeros(n_buckets, 1);
average_refinement_steps = zeros(n_buckets, 1);
max_refinement_steps = zeros(n_buckets, 1);
min_refinement_steps = zeros(n_buckets, 1);

for i=1:n_buckets
    
    for j=1:n_experiments
        g = rand(1);
        g = buckets(i,1) + g*(buckets(i,2) - buckets(i,1));
        g = 10^g;
        
        a = randlog() + 1i*randlog();
        c = randlog() + 1i*randlog();

        A = [0,a;(a)*(1+g),c];

        [AA, Q, info] = type2o_pal(A);
        
        errors(i) = max(errors(i), abs(AA(1,1))/norm(A, 'fro'));
        average_refinement_steps(i) = average_refinement_steps(i) + info.refinement_steps;
        max_refinement_steps(i) = max(max_refinement_steps(i), info.refinement_steps);
        min_refinement_steps(i) = min(min_refinement_steps(i), info.refinement_steps);
    end
    
end
average_refinement_steps = average_refinement_steps./n_experiments;

t = table(buckets, errors, max_refinement_steps, average_refinement_steps, min_refinement_steps)

writetable(t, 'results/type2o_table_complex.csv')