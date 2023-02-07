
load('fluorophore_spec2.mat')

figure
plot(fluop.collagen_ex(:,1),fluop.collagen_ex(:,2),'linewidth',5)		
hold all
plot(fluop.elastin_ex(:,1),fluop.elastin_ex(:,2),'linewidth',5)	
plot(fluop.flavin_ex(:,1),fluop.flavin_ex(:,2),'linewidth',5)	
plot(fluop.lipopigm_ex(:,1),fluop.lipopigm_ex(:,2),'linewidth',5)	
plot(fluop.nadh_ex(:,1),fluop.nadh_ex(:,2),'linewidth',5)	
plot(fluop.pirodox_ex(:,1),fluop.pirodox_ex(:,2),'linewidth',5)	
plot(fluop.porphyr_ex(:,1),fluop.porphyr_ex(:,2),'linewidth',5)	
plot(fluop.tryptop_ex(:,1),fluop.tryptop_ex(:,2),'linewidth',5)	
legend('Collagen','Elastin','Flavin','Lipopigments','NADH','Piridoxine','Porphyrins','Tryptophan','FontSize',14)
xlabel('\lambda [nm]','FontSize',17)
ylabel('Fluorescence "Absorption" [arb. units]','FontSize',17)
title('Fluorophore Excitation Spectra','FontSize',20);
ylim([0 1])
set(gca,'FontSize',20)

figure
plot(fluop.collagen_em(:,1),fluop.collagen_em(:,2),'linewidth',5)		
hold all
plot(fluop.elastin_em(:,1),fluop.elastin_em(:,2),'linewidth',5)	
plot(fluop.flavin_em(:,1),fluop.flavin_em(:,2),'linewidth',5)	
plot(fluop.lipopigm_em(:,1),fluop.lipopigm_em(:,2),'linewidth',5)	
plot(fluop.nadh_em(:,1),fluop.nadh_em(:,2),'linewidth',5)	
plot(fluop.pirodox_em(:,1),fluop.pirodox_em(:,2),'linewidth',5)	
plot(fluop.porphyr_em(:,1),fluop.porphyr_em(:,2),'linewidth',5)	
plot(fluop.tryptop_em(:,1),fluop.tryptop_em(:,2),'linewidth',5)	
legend('Collagen','Elastin','Flavin','Lipopigments','NADH','Piridoxine','Porphyrins','Tryptophan','FontSize',14)
xlabel('\lambda [nm]','FontSize',17)
ylabel('Fluorescence "Emission" [arb. units]','FontSize',17)
title('Fluorophore Emission Spectra','FontSize',20);
ylim([0 1])
set(gca,'FontSize',20)
