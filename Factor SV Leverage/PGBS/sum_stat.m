%N20
mean_B1=mean(IACT_res(1:26))
max_B1=max(IACT_res(1:26))

mean_B2=mean(IACT_res(28:52))
max_B2=max(IACT_res(28:52))

mean_B3=mean(IACT_res(55:78))
max_B3=max(IACT_res(55:78))

mean_B4=mean(IACT_res(82:104))
max_B4=max(IACT_res(82:104))

mean_mu=mean(IACT_res(105:130))
max_mu=max(IACT_res(105:130))

mean_phi=mean(IACT_res(131:160))
max_phi=max(IACT_res(131:160))

mean_tau=mean(IACT_res(161:190))
max_tau=max(IACT_res(161:190))

mean_rho=mean(IACT_res(191:216))
max_rho=max(IACT_res(191:216))

mean_h1=mean(IACT_res(217:316))
max_h1=max(IACT_res(217:316))

mean_h2=mean(IACT_res(317:416))
max_h2=max(IACT_res(317:416))

mean_h10=mean(IACT_res(417:516))
max_h10=max(IACT_res(417:516))

mean_h11=mean(IACT_res(517:616))
max_h11=max(IACT_res(517:616))

mean_h12=mean(IACT_res(617:716))
max_h12=max(IACT_res(617:716))

mean_lambda1=mean(IACT_res(717:816))
max_lambda1=max(IACT_res(717:816))

mean_lambda2=mean(IACT_res(817:916))
max_lambda2=max(IACT_res(817:916))

mean_lambda3=mean(IACT_res(917:1016))
max_lambda3=max(IACT_res(917:1016))

mean_lambda4=mean(IACT_res(1017:1116))
max_lambda4=max(IACT_res(1017:1116))


mean_param=mean([IACT_res(1:26);IACT_res(28:52);IACT_res(55:78);IACT_res(82:104);IACT_res(105:130);IACT_res(131:160);IACT_res(161:190);
    IACT_res(191:216)])
max_param=max([IACT_res(1:26);IACT_res(28:52);IACT_res(55:78);IACT_res(82:104);IACT_res(105:130);IACT_res(131:160);IACT_res(161:190);
    IACT_res(191:216)])


% mean_B1=mean(IACT_res(1:26))
% max_B1=max(IACT_res(1:26))
% 
% mean_mu=mean(IACT_res(27:52))
% max_mu=max(IACT_res(27:52))
% 
% mean_tau=mean(IACT_res(53:79))
% max_tau=max(IACT_res(53:79))
% 
% mean_phi=mean(IACT_res(80:106))
% max_phi=max(IACT_res(80:106))
% 
% mean_rho=mean(IACT_res(107:132))
% max_rho=max(IACT_res(107:132))
% 
% mean_h1=mean(IACT_res(133:232))
% max_h1=max(IACT_res(133:232))
% 
% mean_h2=mean(IACT_res(233:332))
% max_h2=max(IACT_res(233:332))
% 
% mean_h10=mean(IACT_res(333:432))
% max_h10=max(IACT_res(333:432))
% 
% mean_h11=mean(IACT_res(433:532))
% max_h11=max(IACT_res(433:532))
% 
% mean_h12=mean(IACT_res(533:632))
% max_h12=max(IACT_res(533:632))
% 
% mean_lambda1=mean(IACT_res(633:732))
% max_lambda1=max(IACT_res(633:732))
% 
% mean_param=mean([IACT_res(1:26);IACT_res(27:52);IACT_res(53:79);IACT_res(80:106);IACT_res(107:132)])
% max_param=max([IACT_res(1:26);IACT_res(27:52);IACT_res(53:79);IACT_res(80:106);IACT_res(107:132)])

% max_param=max([IACT_PG20(1:26);IACT_PG20(28:52);IACT_PG20(55:78);IACT_PG20(82:104);IACT_PG20(105:130);IACT_PG20(131:160);IACT_PG20(161:190);
%     IACT_PG20(191:216)])


% %N20
% min_B1=min(IACT_PG20(1:26))
% mean_B1=mean(IACT_PG20(1:26))
% max_B1=max(IACT_PG20(1:26))
% 
% min_B2=min(IACT_PG20(28:52))
% mean_B2=mean(IACT_PG20(28:52))
% max_B2=max(IACT_PG20(28:52))
% 
% min_B3=min(IACT_PG20(55:78))
% mean_B3=mean(IACT_PG20(55:78))
% max_B3=max(IACT_PG20(55:78))
% 
% min_B4=min(IACT_PG20(82:104))
% mean_B4=mean(IACT_PG20(82:104))
% max_B4=max(IACT_PG20(82:104))
% 
% min_mu=min(IACT_PG20(105:130))
% mean_mu=mean(IACT_PG20(105:130))
% max_mu=max(IACT_PG20(105:130))
% 
% min_phi=min(IACT_PG20(131:160))
% mean_phi=mean(IACT_PG20(131:160))
% max_phi=max(IACT_PG20(131:160))
% 
% min_tau=min(IACT_PG20(161:190))
% mean_tau=mean(IACT_PG20(161:190))
% max_tau=max(IACT_PG20(161:190))
% 
% min_rho=min(IACT_PG20(191:216))
% mean_rho=mean(IACT_PG20(191:216))
% max_rho=max(IACT_PG20(191:216))
% 
% min_h1=min(IACT_PG20(217:316))
% mean_h1=mean(IACT_PG20(217:316))
% max_h1=max(IACT_PG20(217:316))
% 
% min_h2=min(IACT_PG20(317:416))
% mean_h2=mean(IACT_PG20(317:416))
% max_h2=max(IACT_PG20(317:416))
% 
% min_h10=min(IACT_PG20(417:516))
% mean_h10=mean(IACT_PG20(417:516))
% max_h10=max(IACT_PG20(417:516))
% 
% min_h11=min(IACT_PG20(517:616))
% mean_h11=mean(IACT_PG20(517:616))
% max_h11=max(IACT_PG20(517:616))
% 
% min_h12=min(IACT_PG20(617:716))
% mean_h12=mean(IACT_PG20(617:716))
% max_h12=max(IACT_PG20(617:716))
% 
% mean_param=mean([IACT_PG20(1:26);IACT_PG20(28:52);IACT_PG20(55:78);IACT_PG20(82:104);IACT_PG20(105:130);IACT_PG20(131:160);IACT_PG20(161:190);
%     IACT_PG20(191:216)])
% max_param=max([IACT_PG20(1:26);IACT_PG20(28:52);IACT_PG20(55:78);IACT_PG20(82:104);IACT_PG20(105:130);IACT_PG20(131:160);IACT_PG20(161:190);
%     IACT_PG20(191:216)])
% 
% %N50
% min_B1=min(IACT_PG50(1:26))
% mean_B1=mean(IACT_PG50(1:26))
% max_B1=max(IACT_PG50(1:26))
% 
% min_B2=min(IACT_PG50(28:52))
% mean_B2=mean(IACT_PG50(28:52))
% max_B2=max(IACT_PG50(28:52))
% 
% min_B3=min(IACT_PG50(55:78))
% mean_B3=mean(IACT_PG50(55:78))
% max_B3=max(IACT_PG50(55:78))
% 
% min_B4=min(IACT_PG50(82:104))
% mean_B4=mean(IACT_PG50(82:104))
% max_B4=max(IACT_PG50(82:104))
% 
% min_mu=min(IACT_PG50(105:130))
% mean_mu=mean(IACT_PG50(105:130))
% max_mu=max(IACT_PG50(105:130))
% 
% min_phi=min(IACT_PG50(131:160))
% mean_phi=mean(IACT_PG50(131:160))
% max_phi=max(IACT_PG50(131:160))
% 
% min_tau=min(IACT_PG50(161:190))
% mean_tau=mean(IACT_PG50(161:190))
% max_tau=max(IACT_PG50(161:190))
% 
% min_rho=min(IACT_PG50(191:216))
% mean_rho=mean(IACT_PG50(191:216))
% max_rho=max(IACT_PG50(191:216))
% 
% min_h1=min(IACT_PG50(217:316))
% mean_h1=mean(IACT_PG50(217:316))
% max_h1=max(IACT_PG50(217:316))
% 
% min_h2=min(IACT_PG50(317:416))
% mean_h2=mean(IACT_PG50(317:416))
% max_h2=max(IACT_PG50(317:416))
% 
% min_h10=min(IACT_PG50(417:516))
% mean_h10=mean(IACT_PG50(417:516))
% max_h10=max(IACT_PG50(417:516))
% 
% min_h11=min(IACT_PG50(517:616))
% mean_h11=mean(IACT_PG50(517:616))
% max_h11=max(IACT_PG50(517:616))
% 
% min_h12=min(IACT_PG50(617:716))
% mean_h12=mean(IACT_PG50(617:716))
% max_h12=max(IACT_PG50(617:716))
% 
% mean_param=mean([IACT_PG50(1:26);IACT_PG50(28:52);IACT_PG50(55:78);IACT_PG50(82:104);IACT_PG50(105:130);IACT_PG50(131:160);IACT_PG50(161:190);
%     IACT_PG50(191:216)])
% max_param=max([IACT_PG50(1:26);IACT_PG50(28:52);IACT_PG50(55:78);IACT_PG50(82:104);IACT_PG50(105:130);IACT_PG50(131:160);IACT_PG50(161:190);
%     IACT_PG50(191:216)])
% 
% %N100
% min_B1=min(IACT_PG100(1:26))
% mean_B1=mean(IACT_PG100(1:26))
% max_B1=max(IACT_PG100(1:26))
% 
% min_B2=min(IACT_PG100(28:52))
% mean_B2=mean(IACT_PG100(28:52))
% max_B2=max(IACT_PG100(28:52))
% 
% min_B3=min(IACT_PG100(55:78))
% mean_B3=mean(IACT_PG100(55:78))
% max_B3=max(IACT_PG100(55:78))
% 
% min_B4=min(IACT_PG100(82:104))
% mean_B4=mean(IACT_PG100(82:104))
% max_B4=max(IACT_PG100(82:104))
% 
% min_mu=min(IACT_PG100(105:130))
% mean_mu=mean(IACT_PG100(105:130))
% max_mu=max(IACT_PG100(105:130))
% 
% min_phi=min(IACT_PG100(131:160))
% mean_phi=mean(IACT_PG100(131:160))
% max_phi=max(IACT_PG100(131:160))
% 
% min_tau=min(IACT_PG100(161:190))
% mean_tau=mean(IACT_PG100(161:190))
% max_tau=max(IACT_PG100(161:190))
% 
% min_rho=min(IACT_PG100(191:216))
% mean_rho=mean(IACT_PG100(191:216))
% max_rho=max(IACT_PG100(191:216))
% 
% min_h1=min(IACT_PG100(217:316))
% mean_h1=mean(IACT_PG100(217:316))
% max_h1=max(IACT_PG100(217:316))
% 
% min_h2=min(IACT_PG100(317:416))
% mean_h2=mean(IACT_PG100(317:416))
% max_h2=max(IACT_PG100(317:416))
% 
% min_h10=min(IACT_PG100(417:516))
% mean_h10=mean(IACT_PG100(417:516))
% max_h10=max(IACT_PG100(417:516))
% 
% min_h11=min(IACT_PG100(517:616))
% mean_h11=mean(IACT_PG100(517:616))
% max_h11=max(IACT_PG100(517:616))
% 
% min_h12=min(IACT_PG100(617:716))
% mean_h12=mean(IACT_PG100(617:716))
% max_h12=max(IACT_PG100(617:716))
% 
% mean_param=mean([IACT_PG100(1:26);IACT_PG100(28:52);IACT_PG100(55:78);IACT_PG100(82:104);IACT_PG100(105:130);IACT_PG100(131:160);IACT_PG100(161:190);
%     IACT_PG100(191:216)])
% max_param=max([IACT_PG100(1:26);IACT_PG100(28:52);IACT_PG100(55:78);IACT_PG100(82:104);IACT_PG100(105:130);IACT_PG100(131:160);IACT_PG100(161:190);
%     IACT_PG100(191:216)])
% 
% %N200
% min_B1=min(IACT_PG200(1:26))
% mean_B1=mean(IACT_PG200(1:26))
% max_B1=max(IACT_PG200(1:26))
% 
% min_B2=min(IACT_PG200(28:52))
% mean_B2=mean(IACT_PG200(28:52))
% max_B2=max(IACT_PG200(28:52))
% 
% min_B3=min(IACT_PG200(55:78))
% mean_B3=mean(IACT_PG200(55:78))
% max_B3=max(IACT_PG200(55:78))
% 
% min_B4=min(IACT_PG200(82:104))
% mean_B4=mean(IACT_PG200(82:104))
% max_B4=max(IACT_PG200(82:104))
% 
% min_mu=min(IACT_PG200(105:130))
% mean_mu=mean(IACT_PG200(105:130))
% max_mu=max(IACT_PG200(105:130))
% 
% min_phi=min(IACT_PG200(131:160))
% mean_phi=mean(IACT_PG200(131:160))
% max_phi=max(IACT_PG200(131:160))
% 
% min_tau=min(IACT_PG200(161:190))
% mean_tau=mean(IACT_PG200(161:190))
% max_tau=max(IACT_PG200(161:190))
% 
% min_rho=min(IACT_PG200(191:216))
% mean_rho=mean(IACT_PG200(191:216))
% max_rho=max(IACT_PG200(191:216))
% 
% min_h1=min(IACT_PG200(217:316))
% mean_h1=mean(IACT_PG200(217:316))
% max_h1=max(IACT_PG200(217:316))
% 
% min_h2=min(IACT_PG200(317:416))
% mean_h2=mean(IACT_PG200(317:416))
% max_h2=max(IACT_PG200(317:416))
% 
% min_h10=min(IACT_PG200(417:516))
% mean_h10=mean(IACT_PG200(417:516))
% max_h10=max(IACT_PG200(417:516))
% 
% min_h11=min(IACT_PG200(517:616))
% mean_h11=mean(IACT_PG200(517:616))
% max_h11=max(IACT_PG200(517:616))
% 
% min_h12=min(IACT_PG200(617:716))
% mean_h12=mean(IACT_PG200(617:716))
% max_h12=max(IACT_PG200(617:716))
% 
% mean_param=mean([IACT_PG200(1:26);IACT_PG200(28:52);IACT_PG200(55:78);IACT_PG200(82:104);IACT_PG200(105:130);IACT_PG200(131:160);IACT_PG200(161:190);
%     IACT_PG200(191:216)])
% max_param=max([IACT_PG200(1:26);IACT_PG200(28:52);IACT_PG200(55:78);IACT_PG200(82:104);IACT_PG200(105:130);IACT_PG200(131:160);IACT_PG200(161:190);
%     IACT_PG200(191:216)])
% 
% 
% %---------------------------------
% %corrPMMH+PG
% %N20
% min_B1=min(IACT_corr20(1:26))
% mean_B1=mean(IACT_corr20(1:26))
% max_B1=max(IACT_corr20(1:26))
% 
% min_B2=min(IACT_corr20(28:52))
% mean_B2=mean(IACT_corr20(28:52))
% max_B2=max(IACT_corr20(28:52))
% 
% min_B3=min(IACT_corr20(55:78))
% mean_B3=mean(IACT_corr20(55:78))
% max_B3=max(IACT_corr20(55:78))
% 
% min_B4=min(IACT_corr20(82:104))
% mean_B4=mean(IACT_corr20(82:104))
% max_B4=max(IACT_corr20(82:104))
% 
% min_mu=min(IACT_corr20(105:130))
% mean_mu=mean(IACT_corr20(105:130))
% max_mu=max(IACT_corr20(105:130))
% 
% min_phi=min(IACT_corr20(131:160))
% mean_phi=mean(IACT_corr20(131:160))
% max_phi=max(IACT_corr20(131:160))
% 
% min_tau=min(IACT_corr20(161:190))
% mean_tau=mean(IACT_corr20(161:190))
% max_tau=max(IACT_corr20(161:190))
% 
% min_rho=min(IACT_corr20(191:216))
% mean_rho=mean(IACT_corr20(191:216))
% max_rho=max(IACT_corr20(191:216))
% 
% min_h1=min(IACT_corr20(217:316))
% mean_h1=mean(IACT_corr20(217:316))
% max_h1=max(IACT_corr20(217:316))
% 
% min_h2=min(IACT_corr20(317:416))
% mean_h2=mean(IACT_corr20(317:416))
% max_h2=max(IACT_corr20(317:416))
% 
% min_h10=min(IACT_corr20(417:516))
% mean_h10=mean(IACT_corr20(417:516))
% max_h10=max(IACT_corr20(417:516))
% 
% min_h11=min(IACT_corr20(517:616))
% mean_h11=mean(IACT_corr20(517:616))
% max_h11=max(IACT_corr20(517:616))
% 
% min_h12=min(IACT_corr20(617:716))
% mean_h12=mean(IACT_corr20(617:716))
% max_h12=max(IACT_corr20(617:716))
% 
% mean_param=mean([IACT_corr20(1:26);IACT_corr20(28:52);IACT_corr20(55:78);IACT_corr20(82:104);IACT_corr20(105:130);IACT_corr20(131:160);IACT_corr20(161:190);
%     IACT_corr20(191:216)])
% max_param=max([IACT_corr20(1:26);IACT_corr20(28:52);IACT_corr20(55:78);IACT_corr20(82:104);IACT_corr20(105:130);IACT_corr20(131:160);IACT_corr20(161:190);
%     IACT_corr20(191:216)])
% 
% %N50
% min_B1=min(IACT_corr50(1:26))
% mean_B1=mean(IACT_corr50(1:26))
% max_B1=max(IACT_corr50(1:26))
% 
% min_B2=min(IACT_corr50(28:52))
% mean_B2=mean(IACT_corr50(28:52))
% max_B2=max(IACT_corr50(28:52))
% 
% min_B3=min(IACT_corr50(55:78))
% mean_B3=mean(IACT_corr50(55:78))
% max_B3=max(IACT_corr50(55:78))
% 
% min_B4=min(IACT_corr50(82:104))
% mean_B4=mean(IACT_corr50(82:104))
% max_B4=max(IACT_corr50(82:104))
% 
% min_mu=min(IACT_corr50(105:130))
% mean_mu=mean(IACT_corr50(105:130))
% max_mu=max(IACT_corr50(105:130))
% 
% min_phi=min(IACT_corr50(131:160))
% mean_phi=mean(IACT_corr50(131:160))
% max_phi=max(IACT_corr50(131:160))
% 
% min_tau=min(IACT_corr50(161:190))
% mean_tau=mean(IACT_corr50(161:190))
% max_tau=max(IACT_corr50(161:190))
% 
% min_rho=min(IACT_corr50(191:216))
% mean_rho=mean(IACT_corr50(191:216))
% max_rho=max(IACT_corr50(191:216))
% 
% min_h1=min(IACT_corr50(217:316))
% mean_h1=mean(IACT_corr50(217:316))
% max_h1=max(IACT_corr50(217:316))
% 
% min_h2=min(IACT_corr50(317:416))
% mean_h2=mean(IACT_corr50(317:416))
% max_h2=max(IACT_corr50(317:416))
% 
% min_h10=min(IACT_corr50(417:516))
% mean_h10=mean(IACT_corr50(417:516))
% max_h10=max(IACT_corr50(417:516))
% 
% min_h11=min(IACT_corr50(517:616))
% mean_h11=mean(IACT_corr50(517:616))
% max_h11=max(IACT_corr50(517:616))
% 
% min_h12=min(IACT_corr50(617:716))
% mean_h12=mean(IACT_corr50(617:716))
% max_h12=max(IACT_corr50(617:716))
% 
% mean_param=mean([IACT_corr50(1:26);IACT_corr50(28:52);IACT_corr50(55:78);IACT_corr50(82:104);IACT_corr50(105:130);IACT_corr50(131:160);IACT_corr50(161:190);
%     IACT_corr50(191:216)])
% max_param=max([IACT_corr50(1:26);IACT_corr50(28:52);IACT_corr50(55:78);IACT_corr50(82:104);IACT_corr50(105:130);IACT_corr50(131:160);IACT_corr50(161:190);
%     IACT_corr50(191:216)])
% 
% %N100
% min_B1=min(IACT_corr100(1:26))
% mean_B1=mean(IACT_corr100(1:26))
% max_B1=max(IACT_corr100(1:26))
% 
% min_B2=min(IACT_corr100(28:52))
% mean_B2=mean(IACT_corr100(28:52))
% max_B2=max(IACT_corr100(28:52))
% 
% min_B3=min(IACT_corr100(55:78))
% mean_B3=mean(IACT_corr100(55:78))
% max_B3=max(IACT_corr100(55:78))
% 
% min_B4=min(IACT_corr100(82:104))
% mean_B4=mean(IACT_corr100(82:104))
% max_B4=max(IACT_corr100(82:104))
% 
% min_mu=min(IACT_corr100(105:130))
% mean_mu=mean(IACT_corr100(105:130))
% max_mu=max(IACT_corr100(105:130))
% 
% min_phi=min(IACT_corr100(131:160))
% mean_phi=mean(IACT_corr100(131:160))
% max_phi=max(IACT_corr100(131:160))
% 
% min_tau=min(IACT_corr100(161:190))
% mean_tau=mean(IACT_corr100(161:190))
% max_tau=max(IACT_corr100(161:190))
% 
% min_rho=min(IACT_corr100(191:216))
% mean_rho=mean(IACT_corr100(191:216))
% max_rho=max(IACT_corr100(191:216))
% 
% min_h1=min(IACT_corr100(217:316))
% mean_h1=mean(IACT_corr100(217:316))
% max_h1=max(IACT_corr100(217:316))
% 
% min_h2=min(IACT_corr100(317:416))
% mean_h2=mean(IACT_corr100(317:416))
% max_h2=max(IACT_corr100(317:416))
% 
% min_h10=min(IACT_corr100(417:516))
% mean_h10=mean(IACT_corr100(417:516))
% max_h10=max(IACT_corr100(417:516))
% 
% min_h11=min(IACT_corr100(517:616))
% mean_h11=mean(IACT_corr100(517:616))
% max_h11=max(IACT_corr100(517:616))
% 
% min_h12=min(IACT_corr100(617:716))
% mean_h12=mean(IACT_corr100(617:716))
% max_h12=max(IACT_corr100(617:716))
% 
% mean_param=mean([IACT_corr100(1:26);IACT_corr100(28:52);IACT_corr100(55:78);IACT_corr100(82:104);IACT_corr100(105:130);IACT_corr100(131:160);IACT_corr100(161:190);
%     IACT_corr100(191:216)])
% max_param=max([IACT_corr100(1:26);IACT_corr100(28:52);IACT_corr100(55:78);IACT_corr100(82:104);IACT_corr100(105:130);IACT_corr100(131:160);IACT_corr100(161:190);
%     IACT_corr100(191:216)])
% 
% %N200
% min_B1=min(IACT_corr200(1:26))
% mean_B1=mean(IACT_corr200(1:26))
% max_B1=max(IACT_corr200(1:26))
% 
% min_B2=min(IACT_corr200(28:52))
% mean_B2=mean(IACT_corr200(28:52))
% max_B2=max(IACT_corr200(28:52))
% 
% min_B3=min(IACT_corr200(55:78))
% mean_B3=mean(IACT_corr200(55:78))
% max_B3=max(IACT_corr200(55:78))
% 
% min_B4=min(IACT_corr200(82:104))
% mean_B4=mean(IACT_corr200(82:104))
% max_B4=max(IACT_corr200(82:104))
% 
% min_mu=min(IACT_corr200(105:130))
% mean_mu=mean(IACT_corr200(105:130))
% max_mu=max(IACT_corr200(105:130))
% 
% min_phi=min(IACT_corr200(131:160))
% mean_phi=mean(IACT_corr200(131:160))
% max_phi=max(IACT_corr200(131:160))
% 
% min_tau=min(IACT_corr200(161:190))
% mean_tau=mean(IACT_corr200(161:190))
% max_tau=max(IACT_corr200(161:190))
% 
% min_rho=min(IACT_corr200(191:216))
% mean_rho=mean(IACT_corr200(191:216))
% max_rho=max(IACT_corr200(191:216))
% 
% min_h1=min(IACT_corr200(217:316))
% mean_h1=mean(IACT_corr200(217:316))
% max_h1=max(IACT_corr200(217:316))
% 
% min_h2=min(IACT_corr200(317:416))
% mean_h2=mean(IACT_corr200(317:416))
% max_h2=max(IACT_corr200(317:416))
% 
% min_h10=min(IACT_corr200(417:516))
% mean_h10=mean(IACT_corr200(417:516))
% max_h10=max(IACT_corr200(417:516))
% 
% min_h11=min(IACT_corr200(517:616))
% mean_h11=mean(IACT_corr200(517:616))
% max_h11=max(IACT_corr200(517:616))
% 
% min_h12=min(IACT_corr200(617:716))
% mean_h12=mean(IACT_corr200(617:716))
% max_h12=max(IACT_corr200(617:716))
% 
% mean_param=mean([IACT_corr200(1:26);IACT_corr200(28:52);IACT_corr200(55:78);IACT_corr200(82:104);IACT_corr200(105:130);IACT_corr200(131:160);IACT_corr200(161:190);
%     IACT_corr200(191:216)])
% max_param=max([IACT_corr200(1:26);IACT_corr200(28:52);IACT_corr200(55:78);IACT_corr200(82:104);IACT_corr200(105:130);IACT_corr200(131:160);IACT_corr200(161:190);
%     IACT_corr200(191:216)])