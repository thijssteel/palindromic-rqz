%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Tests of type IIo moves on complex matrices %%%%%%%%%%%%%%%%%%%%
%%%%% results are saved in results/type2o_refinement_complex.csv' %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_gaps = 100; %the number of values for g that are tested
n_experiments = 10000; %the number of experiments per value of g
g_array = logspace(-15,3,n_gaps)';


errors = zeros(n_gaps, 1);
no_convergence = zeros(n_gaps, 1);
average_refinement_steps = zeros(n_gaps, 1);
max_refinement_steps = zeros(n_gaps, 1);
min_refinement_steps = zeros(n_gaps, 1);

for i=1:n_gaps
    g = g_array(i);
    
    for j=1:n_experiments
        a = randlog() + 1i*randlog();
        c = randlog() + 1i*randlog();

        A = [0,a;(a)*(1+g),c];

        [AA, Q, info] = type2o_pal(A);
        
        errors(i) = max(errors(i), abs(AA(1,1))/norm(A, 'fro'));
        no_convergence(i) = no_convergence(i) + info.no_convergence;
        average_refinement_steps(i) = average_refinement_steps(i) + info.refinement_steps;
        max_refinement_steps(i) = max(max_refinement_steps(i), info.refinement_steps);
        min_refinement_steps(i) = min(min_refinement_steps(i), info.refinement_steps);
    end
    
end
no_convergence = no_convergence./n_experiments;
average_refinement_steps = average_refinement_steps./n_experiments;

figure(1);clf;
loglog(g_array, errors);
xlabel('g')
ylabel('maximum error')

figure(2);clf;
semilogx(g_array, average_refinement_steps);
hold on;
semilogx(g_array, max_refinement_steps);
semilogx(g_array, min_refinement_steps);
xlabel('g')
ylabel('refinement steps')
legend('average', 'max', 'min')

writetable(table(g_array, max_refinement_steps, average_refinement_steps, min_refinement_steps), 'results/type2o_refinement_complex.csv')
writetable(table(g_array, errors), 'results/type2o_errors_complex.csv')
