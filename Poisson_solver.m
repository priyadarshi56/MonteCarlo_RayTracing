%----- Poisson solver -----------------------------------------------------
if strcmpi(include_pores,'Yes')  && strcmpi(include_electrolyte,'Yes')
    tic
    Poisson2D_SC_Electrolyte;
    poisson_time=toc;
    % plot_Poisson_results;
    fprintf('\n2D Poisson solver time= %g seconds\n\n',poisson_time);
    avg_Ec_poisson=mean(mean(Ec_matrix));
    avg_n_poisson=mean(mean(n_matrix));
elseif strcmpi(include_ordered_poly,'Yes')
    tic
    Poisson2D_SC_SC_Poly;
    poisson_time=toc;
    % plot_Poisson_results_poly;
    fprintf('\n2D Poisson solver time= %g seconds\n \n',poisson_time);
    avg_Ec_poisson=mean(mean(Ec_matrix));
    avg_n_poisson=mean(mean(n_matrix));
elseif strcmpi(include_gb_doping,'Yes')
    tic
    Poisson2D_SC_SC_GB;
    poisson_time=toc;
    % plot_Poisson_results_poly;
    fprintf('\n2D Poisson solver time= %g seconds\n\n',poisson_time);
    avg_Ec_poisson=mean(mean(Ec_matrix));
    avg_n_poisson=mean(mean(n_matrix));
else
    poisson_time=0;
end