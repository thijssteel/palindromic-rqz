%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Tests of type IIe moves on complex matrices %%%%%%%%%%%%%%%%%%%%
%%%%% results are saved in results/type2e_refinement_complex.csv' %%%%
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
        b = randlog() + 1i*randlog();
        c = randlog() + 1i*randlog();
        d = randlog() + 1i*randlog();
        e = randlog() + 1i*randlog();
        f = randlog() + 1i*randlog();


        A = [0,0,b;
            0,a,d;
            b*(1+g),e,f];

        [AA, Q,info] = type2e_pal(A);
        
        errors(i) = max(errors(i), max( [abs(AA(1,1)),abs(AA(2,1)),abs(AA(1,2)) ])/norm(A, 'fro'));
        average_refinement_steps(i) = average_refinement_steps(i) + info.refinement_steps;
        max_refinement_steps(i) = max(max_refinement_steps(i), info.refinement_steps);
        min_refinement_steps(i) = min(min_refinement_steps(i), info.refinement_steps);
    end
    
end
average_refinement_steps = average_refinement_steps./n_experiments;

t = table(buckets, errors, max_refinement_steps, average_refinement_steps, min_refinement_steps)

writetable(t, 'results/type2e_table_complex.csv')