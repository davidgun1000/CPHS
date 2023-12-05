% xlswrite('PG_real_leverage_N100_T1000.xlsx',[Post.B1(2001:end,:),Post.mu(2001:end,:),Post.tau(2001:end,:),Post.phi(2001:end,:),Post.rho(2001:end,:),Post.ctraj1(2001:end,:),Post.ctraj2(2001:end,:),Post.ctraj10(2001:end,:),...
%     Post.ctraj11(2001:end,:),Post.ctraj12(2001:end,:),Post.ctrajfactor(2001:end,:)]);
% xlswrite('PG_real_leverage_N250_T1000.xlsx',[Post.B1(2001:end,:),Post.mu(2001:end,:),Post.tau(2001:end,:),Post.phi(2001:end,:),Post.rho(2001:end,:),Post.ctraj1(2001:end,:),Post.ctraj2(2001:end,:),Post.ctraj10(2001:end,:),...
%     Post.ctraj11(2001:end,:),Post.ctraj12(2001:end,:),Post.ctrajfactor(2001:end,:)]);
% xlswrite('PG_real_leverage_N500_T1000.xlsx',[Post.B1(2001:end,:),Post.mu(2001:end,:),Post.tau(2001:end,:),Post.phi(2001:end,:),Post.rho(2001:end,:),Post.ctraj1(2001:end,:),Post.ctraj2(2001:end,:),Post.ctraj10(2001:end,:),...
%     Post.ctraj11(2001:end,:),Post.ctraj12(2001:end,:),Post.ctrajfactor(2001:end,:)]);
% xlswrite('PG_real_leverage_N1000_T1000.xlsx',[Post.B1(2001:end,:),Post.mu(2001:end,:),Post.tau(2001:end,:),Post.phi(2001:end,:),Post.rho(2001:end,:),Post.ctraj1(2001:end,:),Post.ctraj2(2001:end,:),Post.ctraj10(2001:end,:),...
%     Post.ctraj11(2001:end,:),Post.ctraj12(2001:end,:),Post.ctrajfactor(2001:end,:)]);
% 
% xlswrite('mix_real_leverage_N100_T1000.xlsx',[Post.B1(2001:end,:),Post.mu(2001:end,:),Post.tau(2001:end,:),Post.phi(2001:end,:),Post.rho(2001:end,:),Post.ctraj1(2001:end,:),Post.ctraj2(2001:end,:),Post.ctraj10(2001:end,:),...
%     Post.ctraj11(2001:end,:),Post.ctraj12(2001:end,:),Post.ctrajfactor1(2001:end,:)]);
% xlswrite('mix_real_leverage_N250_T1000.xlsx',[Post.B1(2001:end,:),Post.mu(2001:end,:),Post.tau(2001:end,:),Post.phi(2001:end,:),Post.rho(2001:end,:),Post.ctraj1(2001:end,:),Post.ctraj2(2001:end,:),Post.ctraj10(2001:end,:),...
%     Post.ctraj11(2001:end,:),Post.ctraj12(2001:end,:),Post.ctrajfactor1(2001:end,:)]);
% xlswrite('mix_real_leverage_N500_T1000.xlsx',[Post.B1(2001:end,:),Post.mu(2001:end,:),Post.tau(2001:end,:),Post.phi(2001:end,:),Post.rho(2001:end,:),Post.ctraj1(2001:end,:),Post.ctraj2(2001:end,:),Post.ctraj10(2001:end,:),...
%     Post.ctraj11(2001:end,:),Post.ctraj12(2001:end,:),Post.ctrajfactor1(2001:end,:)]);
% xlswrite('mix_real_leverage_N1000_T1000.xlsx',[Post.B1(2001:end,:),Post.mu(2001:end,:),Post.tau(2001:end,:),Post.phi(2001:end,:),Post.rho(2001:end,:),Post.ctraj1(2001:end,:),Post.ctraj2(2001:end,:),Post.ctraj10(2001:end,:),...
%     Post.ctraj11(2001:end,:),Post.ctraj12(2001:end,:),Post.ctrajfactor1(2001:end,:)]);
% 
% xlswrite('corrmix_real_leverage_N100_T1000.xlsx',[Post.B1(2001:end,:),Post.mu(2001:end,:),Post.tau(2001:end,:),Post.phi(2001:end,:),Post.rho(2001:end,:),Post.ctraj1(2001:end,:),Post.ctraj2(2001:end,:),Post.ctraj10(2001:end,:),...
%     Post.ctraj11(2001:end,:),Post.ctraj12(2001:end,:),Post.ctrajfactor(2001:end,:)]);

xlswrite('standardmix_real_leverage_N1000_T3000.xlsx',[Post.phi(1001:end,1),Post.mu(1001:end,1),Post.tau(1001:end,1),Post.rho(1001:end,1),Post.ctraj(1001:end,:)]);
xlswrite('standardmix_real_leverage_N2000_T3000.xlsx',[Post.phi(1001:end,1),Post.mu(1001:end,1),Post.tau(1001:end,1),Post.rho(1001:end,1),Post.ctraj(1001:end,:)]);
xlswrite('standardmix_real_leverage_N5000_T3000.xlsx',[Post.phi(1001:end,1),Post.mu(1001:end,1),Post.tau(1001:end,1),Post.rho(1001:end,1),Post.ctraj(1001:end,:)]);

xlswrite('PG_univ_real_leverage_N500_T3000.xlsx',[Post.phi(5001:end,1),Post.mu(5001:end,1),Post.tau(5001:end,1),Post.rho(5001:end,1),Post.ctraj(5001:end,:)]);
xlswrite('PG_univ_real_leverage_N1000_T3000.xlsx',[Post.phi(5001:end,1),Post.mu(5001:end,1),Post.tau(5001:end,1),Post.rho(5001:end,1),Post.ctraj(5001:end,:)]);

xlswrite('corrPMMH_N20_T3000_leverage.xlsx',[Post.phi(5001:end,1),Post.mu(5001:end,1),Post.tau(5001:end,1),Post.rho(5001:end,1),Post.ctraj(5001:end,:)]);
xlswrite('corrPMMH_N50_T3000_leverage.xlsx',[Post.phi(5001:end,1),Post.mu(5001:end,1),Post.tau(5001:end,1),Post.rho(5001:end,1),Post.ctraj(5001:end,:)]);
xlswrite('corrPMMH_N100_T3000_leverage.xlsx',[Post.phi(5001:end,1),Post.mu(5001:end,1),Post.tau(5001:end,1),Post.rho(5001:end,1),Post.ctraj(5001:end,:)]);

csvwrite('corrMix_N100_T2000_leverage.csv',[Post.phi(5001:end,1),Post.mu(5001:end,1),Post.tau(5001:end,1),Post.rho(5001:end,1),Post.ctraj(5001:end,:)]);
csvwrite('corrMix_N100_T4000_leverage.csv',[Post.phi(5001:end,1),Post.mu(5001:end,1),Post.tau(5001:end,1),Post.rho(5001:end,1),Post.ctraj(5001:end,:)]);
csvwrite('corrMix_N100_T6000_leverage.csv',[Post.phi(5001:end,1),Post.mu(5001:end,1),Post.tau(5001:end,1),Post.rho(5001:end,1),Post.ctraj(5001:end,:)]);
csvwrite('corrMix_N100_T8000_leverage.csv',[Post.phi(5001:end,1),Post.mu(5001:end,1),Post.tau(5001:end,1),Post.rho(5001:end,1),Post.ctraj(5001:end,:)]);
csvwrite('corrMix_N100_T10000_leverage.csv',[Post.phi(5001:end,1),Post.mu(5001:end,1),Post.tau(5001:end,1),Post.rho(5001:end,1),Post.ctraj(5001:end,:)]);
csvwrite('corrMix_N100_T12000_leverage.csv',[Post.phi(5001:end,1),Post.mu(5001:end,1),Post.tau(5001:end,1),Post.rho(5001:end,1),Post.ctraj(5001:end,:)]);
csvwrite('corrMix_N100_T14000_leverage.csv',[Post.phi(5001:end,1),Post.mu(5001:end,1),Post.tau(5001:end,1),Post.rho(5001:end,1),Post.ctraj(5001:end,:)]);
csvwrite('corrMix_N100_T16000_leverage.csv',[Post.phi(5001:end,1),Post.mu(5001:end,1),Post.tau(5001:end,1),Post.rho(5001:end,1),Post.ctraj(5001:end,:)]);
csvwrite('corrMix_N100_T18000_leverage.csv',[Post.phi(5001:end,1),Post.mu(5001:end,1),Post.tau(5001:end,1),Post.rho(5001:end,1),Post.ctraj(5001:end,:)]);
csvwrite('corrMix_N100_T20000_leverage.csv',[Post.phi(5001:end,1),Post.mu(5001:end,1),Post.tau(5001:end,1),Post.rho(5001:end,1),Post.ctraj(5001:end,:)]);

csvwrite('corrPMMH_N100_T2000_leverage.csv',[Post.phi(5001:end,1),Post.mu(5001:end,1),Post.tau(5001:end,1),Post.rho(5001:end,1),Post.ctraj(5001:end,:)]);
csvwrite('corrPMMH_N100_T4000_leverage.csv',[Post.phi(5001:end,1),Post.mu(5001:end,1),Post.tau(5001:end,1),Post.rho(5001:end,1),Post.ctraj(5001:end,:)]);
csvwrite('corrPMMH_N100_T6000_leverage.csv',[Post.phi(5001:end,1),Post.mu(5001:end,1),Post.tau(5001:end,1),Post.rho(5001:end,1),Post.ctraj(5001:end,:)]);
csvwrite('corrPMMH_N100_T8000_leverage.csv',[Post.phi(5001:end,1),Post.mu(5001:end,1),Post.tau(5001:end,1),Post.rho(5001:end,1),Post.ctraj(5001:end,:)]);
csvwrite('corrPMMH_N100_T10000_leverage.csv',[Post.phi(5001:end,1),Post.mu(5001:end,1),Post.tau(5001:end,1),Post.rho(5001:end,1),Post.ctraj(5001:end,:)]);
csvwrite('corrPMMH_N100_T12000_leverage.csv',[Post.phi(5001:end,1),Post.mu(5001:end,1),Post.tau(5001:end,1),Post.rho(5001:end,1),Post.ctraj(5001:end,:)]);
csvwrite('corrPMMH_N100_T14000_leverage.csv',[Post.phi(5001:end,1),Post.mu(5001:end,1),Post.tau(5001:end,1),Post.rho(5001:end,1),Post.ctraj(5001:end,:)]);
csvwrite('corrPMMH_N100_T16000_leverage.csv',[Post.phi(5001:end,1),Post.mu(5001:end,1),Post.tau(5001:end,1),Post.rho(5001:end,1),Post.ctraj(5001:end,:)]);
csvwrite('corrPMMH_N100_T18000_leverage.csv',[Post.phi(5001:end,1),Post.mu(5001:end,1),Post.tau(5001:end,1),Post.rho(5001:end,1),Post.ctraj(5001:end,:)]);
csvwrite('corrPMMH_N100_T20000_leverage.csv',[Post.phi(5001:end,1),Post.mu(5001:end,1),Post.tau(5001:end,1),Post.rho(5001:end,1),Post.ctraj(5001:end,:)]);



% %dim2
% csvwrite('corrMix_LGSS_N100_T500_dim2.csv',[Post.theta(5001:end,1),Post.ctraj1{1,1}(5001:end,:),Post.ctraj1{2,1}(5001:end,:)]);
% csvwrite('corrPMMH_LGSS_N250_T500_dim2.csv',[Post.theta(5001:end,1),Post.ctraj1{1,1}(5001:end,:),Post.ctraj1{2,1}(5001:end,:)]);
% 
% %dim3
% csvwrite('corrMix_LGSS_N100_T500_dim3.csv',[Post.theta(5001:end,1),Post.ctraj1{1,1}(5001:end,:),Post.ctraj1{2,1}(5001:end,:),Post.ctraj1{3,1}(5001:end,:)]);
% csvwrite('corrPMMH_LGSS_N250_T500_dim3.csv',[Post.theta(5001:end,1),Post.ctraj1{1,1}(5001:end,:),Post.ctraj1{2,1}(5001:end,:),Post.ctraj1{3,1}(5001:end,:)]);
% 
% %dim4
% csvwrite('corrMix_LGSS_N250_T500_dim4.csv',[Post.theta(5001:end,1),Post.ctraj1{1,1}(5001:end,:),Post.ctraj1{2,1}(5001:end,:),Post.ctraj1{3,1}(5001:end,:),Post.ctraj1{4,1}(5001:end,:)]);
% csvwrite('corrPMMH_LGSS_N500_T500_dim4.csv',[Post.theta(5001:end,1),Post.ctraj1{1,1}(5001:end,:),Post.ctraj1{2,1}(5001:end,:),Post.ctraj1{3,1}(5001:end,:),Post.ctraj1{4,1}(5001:end,:)]);
% 
% %dim2
% csvwrite('corrMix_LGSS_N100_T1000_dim2.csv',[Post.theta(5001:end,1),Post.ctraj1(5001:end,:),Post.ctraj2(5001:end,:)]);
% csvwrite('corrPMMH_LGSS_N250_T1000_dim2.csv',[Post.theta(5001:end,1),Post.ctraj1{1,1}(5001:end,:),Post.ctraj1{2,1}(5001:end,:)]);
% 
% %dim3
% csvwrite('corrMix_LGSS_N250_T1000_dim3.csv',[Post.theta(5001:end,1),Post.ctraj1{1,1}(5001:end,:),Post.ctraj1{2,1}(5001:end,:),Post.ctraj1{3,1}(5001:end,:)]);
% csvwrite('corrPMMH_LGSS_N500_T1000_dim3.csv',[Post.theta(5001:end,1),Post.ctraj1{1,1}(5001:end,:),Post.ctraj1{2,1}(5001:end,:),Post.ctraj1{3,1}(5001:end,:)]);
% 
% %dim4
% csvwrite('corrMix_LGSS_N250_T1000_dim4.csv',[Post.theta(5001:end,1),Post.ctraj1{1,1}(5001:end,:),Post.ctraj1{2,1}(5001:end,:),Post.ctraj1{3,1}(5001:end,:),Post.ctraj1{4,1}(5001:end,:)]);
% csvwrite('corrPMMH_LGSS_N500_T1000_dim4.csv',[Post.theta(5001:end,1),Post.ctraj1{1,1}(5001:end,:),Post.ctraj1{2,1}(5001:end,:),Post.ctraj1{3,1}(5001:end,:),Post.ctraj1{4,1}(5001:end,:)]);
% 
% 
% csvwrite('corrMix_PHS_GARCH.csv',[Post.theta(5001:end,1),Post.ctraj1{1,1}(5001:end,:),Post.ctraj1{2,1}(5001:end,:)]);

csvwrite('CorrMixGARCH.csv',[Post.B1(5001:end,:),Post.B2(5001:end,2:end),Post.B3(5001:end,3:end),Post.B4(5001:end,4:end),Post.a(5001:end,:),Post.mu(5001:end,:),Post.tau(5001:end,:),...
    Post.phi(5001:end,:),Post.tau_factor(5001:end,:)]);
csvwrite('PGGARCH500.csv',[Post.B1(5001:end,:),Post.B2(5001:end,2:end),Post.B3(5001:end,3:end),Post.B4(5001:end,4:end),Post.a(5001:end,:),Post.mu(5001:end,:),Post.tau(5001:end,:),...
    Post.phi(5001:end,:),Post.tau_factor(5001:end,:)]);
csvwrite('PGGARCH1000.csv',[Post.B1(5001:end,:),Post.B2(5001:end,2:end),Post.B3(5001:end,3:end),Post.B4(5001:end,4:end),Post.a(5001:end,:),Post.mu(5001:end,:),Post.tau(5001:end,:),...
    Post.phi(5001:end,:),Post.tau_factor(5001:end,:)]);
csvwrite('STANDARDGARCH500.csv',[Post.B1(5001:end,:),Post.B2(5001:end,2:end),Post.B3(5001:end,3:end),Post.B4(5001:end,4:end),Post.a(5001:end,:),Post.mu(5001:end,:),Post.tau(5001:end,:),...
    Post.phi(5001:end,:),Post.tau_factor(5001:end,:)]);
csvwrite('STANDARDGARCH1000.csv',[Post.B1(5001:end,:),Post.B2(5001:end,2:end),Post.B3(5001:end,3:end),Post.B4(5001:end,4:end),Post.a(5001:end,:),Post.mu(5001:end,:),Post.tau(5001:end,:),...
    Post.phi(5001:end,:),Post.tau_factor(5001:end,:)]);

mean_B1=mean(x(1:26))
max_B1=max(x(1:26))

mean_B2=mean(x(27:51))
max_B2=max(x(27:51))

mean_B3=mean(x(52:75))
max_B3=max(x(52:75))

mean_B4=mean(x(76:98))
max_B4=max(x(76:98))

mean_a=mean(x(99:124))
max_a=max(x(99:124))

mean_mu=mean(x(125:150))
max_mu=max(x(125:150))

mean_tau=mean([x(151:176);x(181:184)])
max_tau=max([x(151:176);x(181:184)])

mean_phi=mean(x(177:180))
max_phi=max(x(177:180))

mean_param=mean(x(1:184))
max_param=max(x(1:184))

% mean_tau_factor=mean(x(181:184))
% max_tau_factor=max(x(181:184))



% xlswrite('corr_real_leverage_N20.xlsx',[Post.B1(2001:end,:),Post.B2(2001:end,:),Post.B3(2001:end,:),Post.B4(2001:end,:),Post.mu(2001:end,:),Post.tau(2001:end,:),Post.phi(2001:end,:),Post.rho(2001:end,:),Post.ctraj1(2001:end,:),Post.ctraj2(2001:end,:),Post.ctraj10(2001:end,:),...
%     Post.ctraj11(2001:end,:),Post.ctraj12(2001:end,:)]);
% xlswrite('corr_real_leverage_N50.xlsx',[Post.B1(2001:end,:),Post.B2(2001:end,:),Post.B3(2001:end,:),Post.B4(2001:end,:),Post.mu(2001:end,:),Post.tau(2001:end,:),Post.phi(2001:end,:),Post.rho(2001:end,:),Post.ctraj1(2001:end,:),Post.ctraj2(2001:end,:),Post.ctraj10(2001:end,:),...
%     Post.ctraj11(2001:end,:),Post.ctraj12(2001:end,:)]);
% xlswrite('corr_real_leverage_N100.xlsx',[Post.B1(2001:end,:),Post.B2(2001:end,:),Post.B3(2001:end,:),Post.B4(2001:end,:),Post.mu(2001:end,:),Post.tau(2001:end,:),Post.phi(2001:end,:),Post.rho(2001:end,:),Post.ctraj1(2001:end,:),Post.ctraj2(2001:end,:),Post.ctraj10(2001:end,:),...
%     Post.ctraj11(2001:end,:),Post.ctraj12(2001:end,:)]);
% xlswrite('corr_real_leverage_N200.xlsx',[Post.B1(2001:end,:),Post.B2(2001:end,:),Post.B3(2001:end,:),Post.B4(2001:end,:),Post.mu(2001:end,:),Post.tau(2001:end,:),Post.phi(2001:end,:),Post.rho(2001:end,:),Post.ctraj1(2001:end,:),Post.ctraj2(2001:end,:),Post.ctraj10(2001:end,:),...
%     Post.ctraj11(2001:end,:),Post.ctraj12(2001:end,:)]);
% 
% xlswrite('PG_real_leverage_N20.xlsx',[Post.B1(2001:end,:),Post.B2(2001:end,:),Post.B3(2001:end,:),Post.B4(2001:end,:),Post.mu(2001:end,:),Post.tau(2001:end,:),Post.phi(2001:end,:),Post.rho(2001:end,:),Post.ctraj1(2001:end,:),Post.ctraj2(2001:end,:),Post.ctraj10(2001:end,:),...
%     Post.ctraj11(2001:end,:),Post.ctraj12(2001:end,:)]);
% xlswrite('PG_real_leverage_N50.xlsx',[Post.B1(2001:end,:),Post.B2(2001:end,:),Post.B3(2001:end,:),Post.B4(2001:end,:),Post.mu(2001:end,:),Post.tau(2001:end,:),Post.phi(2001:end,:),Post.rho(2001:end,:),Post.ctraj1(2001:end,:),Post.ctraj2(2001:end,:),Post.ctraj10(2001:end,:),...
%     Post.ctraj11(2001:end,:),Post.ctraj12(2001:end,:)]);
% xlswrite('PG_real_leverage_N100.xlsx',[Post.B1(2001:end,:),Post.B2(2001:end,:),Post.B3(2001:end,:),Post.B4(2001:end,:),Post.mu(2001:end,:),Post.tau(2001:end,:),Post.phi(2001:end,:),Post.rho(2001:end,:),Post.ctraj1(2001:end,:),Post.ctraj2(2001:end,:),Post.ctraj10(2001:end,:),...
%     Post.ctraj11(2001:end,:),Post.ctraj12(2001:end,:)]);
% xlswrite('PG_real_leverage_N200.xlsx',[Post.B1(2001:end,:),Post.B2(2001:end,:),Post.B3(2001:end,:),Post.B4(2001:end,:),Post.mu(2001:end,:),Post.tau(2001:end,:),Post.phi(2001:end,:),Post.rho(2001:end,:),Post.ctraj1(2001:end,:),Post.ctraj2(2001:end,:),Post.ctraj10(2001:end,:),...
%     Post.ctraj11(2001:end,:),Post.ctraj12(2001:end,:)]);
% 
% xlswrite('corr_real_leverage_N20_factor.xlsx',[Post.ctrajfactor1(2001:end,:),Post.ctrajfactor2(2001:end,:),Post.ctrajfactor3(2001:end,:),Post.ctrajfactor4(2001:end,:)]);
% xlswrite('corr_real_leverage_N50_factor.xlsx',[Post.ctrajfactor1(2001:end,:),Post.ctrajfactor2(2001:end,:),Post.ctrajfactor3(2001:end,:),Post.ctrajfactor4(2001:end,:)]);
% xlswrite('corr_real_leverage_N100_factor.xlsx',[Post.ctrajfactor1(2001:end,:),Post.ctrajfactor2(2001:end,:),Post.ctrajfactor3(2001:end,:),Post.ctrajfactor4(2001:end,:)]);
% xlswrite('corr_real_leverage_N200_factor.xlsx',[Post.ctrajfactor1(2001:end,:),Post.ctrajfactor2(2001:end,:),Post.ctrajfactor3(2001:end,:),Post.ctrajfactor4(2001:end,:)]);
% 
% xlswrite('PG_real_leverage_N20_factor.xlsx',[Post.ctrajfactor1(2001:end,:),Post.ctrajfactor2(2001:end,:),Post.ctrajfactor3(2001:end,:),Post.ctrajfactor4(2001:end,:)]);
% xlswrite('PG_real_leverage_N50_factor.xlsx',[Post.ctrajfactor1(2001:end,:),Post.ctrajfactor2(2001:end,:),Post.ctrajfactor3(2001:end,:),Post.ctrajfactor4(2001:end,:)]);
% xlswrite('PG_real_leverage_N100_factor.xlsx',[Post.ctrajfactor1(2001:end,:),Post.ctrajfactor2(2001:end,:),Post.ctrajfactor3(2001:end,:),Post.ctrajfactor4(2001:end,:)]);
% xlswrite('PG_real_leverage_N200_factor.xlsx',[Post.ctrajfactor1(2001:end,:),Post.ctrajfactor2(2001:end,:),Post.ctrajfactor3(2001:end,:),Post.ctrajfactor4(2001:end,:)]);
% 


% xlswrite('mix_real_N100.xlsx',[Post.B1(2001:end,:),Post.mu(2001:end,:),Post.tau(2001:end,:),Post.phi(2001:end,:),Post.ctraj1(2001:end,:),Post.ctraj2(2001:end,:),Post.ctraj10(2001:end,:),...
%     Post.ctraj11(2001:end,:),Post.ctraj12(2001:end,:)]);
% xlswrite('mix_real_N250.xlsx',[Post.B1(2001:end,:),Post.mu(2001:end,:),Post.tau(2001:end,:),Post.phi(2001:end,:),Post.ctraj1(2001:end,:),Post.ctraj2(2001:end,:),Post.ctraj10(2001:end,:),...
%     Post.ctraj11(2001:end,:),Post.ctraj12(2001:end,:)]);
% xlswrite('mix_real_N500.xlsx',[Post.B1(2001:end,:),Post.mu(2001:end,:),Post.tau(2001:end,:),Post.phi(2001:end,:),Post.ctraj1(2001:end,:),Post.ctraj2(2001:end,:),Post.ctraj10(2001:end,:),...
%     Post.ctraj11(2001:end,:),Post.ctraj12(2001:end,:)]);
% xlswrite('mix_real_N1000.xlsx',[Post.B1(2001:end,:),Post.mu(2001:end,:),Post.tau(2001:end,:),Post.phi(2001:end,:),Post.ctraj1(2001:end,:),Post.ctraj2(2001:end,:),Post.ctraj10(2001:end,:),...
%     Post.ctraj11(2001:end,:),Post.ctraj12(2001:end,:)]);
% xlswrite('PG_EIS_real_N100.xlsx',[Post.B1(2001:end,:),Post.mu(2001:end,:),Post.tau(2001:end,:),Post.phi(2001:end,:),Post.ctraj1(2001:end,:),Post.ctraj2(2001:end,:),Post.ctraj10(2001:end,:),...
%     Post.ctraj11(2001:end,:),Post.ctraj12(2001:end,:)]);
% xlswrite('PG_real_N100_leverage.xlsx',[Post.B1(2001:end,:),Post.mu(2001:end,:),Post.tau(2001:end,:),Post.phi(2001:end,:),Post.rho(2001:end,:),Post.ctraj1(2001:end,:),Post.ctraj2(2001:end,:),Post.ctraj10(2001:end,:),...
%     Post.ctraj11(2001:end,:),Post.ctraj12(2001:end,:)]);
% xlswrite('corr_real_N100_leverage.xlsx',[Post.B1(2001:end,:),Post.mu(2001:end,:),Post.tau(2001:end,:),Post.phi(2001:end,:),Post.rho(2001:end,:),Post.ctraj1(2001:end,:),Post.ctraj2(2001:end,:),Post.ctraj10(2001:end,:),...
%     Post.ctraj11(2001:end,:),Post.ctraj12(2001:end,:)]);

