%% CHROMOPHORE PREPROCESS

% Interpolation, fitting and stitching used to preprocess original chromophore spectra (chromophores4.mat) to allow for
% extension from 430 nm to 340 nm (chromophores.mat) using omlc.org data
% and MetHb.csv file.

% ORIGINAL DATA (chromophores4.mat)
% Rami Nachabé, Benno H. W. Hendriks, Marjolein van der Voort, Adrien E. Desjardins, and Henricus J. C. M. Sterenborg, "Estimation of biological chromophores using diffuse optical spectroscopy: benefit of extending the UV-VIS wavelength range to include 1000 to 1600 nm," Biomed. Opt. Express 1, 1432-1442 (2010)

% WATER: https://omlc.org/spectra/water/data/buiteveld94.txt
% H. Buiteveld and J. M. H. Hakvoort and M. Donze, "The optical properties
% of pure water," in SPIE Proceedings on Ocean Optics XII, edited by J. S.
% Jaffe, 2258, 174--183, (1994).

% Hb and HbO2: https://omlc.org/spectra/hemoglobin/summary.html
% W. B. Gratzer, Med. Res. Council Labs, Holly Hill, London
% N. Kollias, Wellman Laboratories, Harvard Medical School, Boston

% BILIRUBIN (main Bile chromophore): https://omlc.org/spectra/PhotochemCAD/data/119-abs.txt
% H. Du, R. A. Fuh, J. Li, A. Corkan, J. S. Lindsey, "PhotochemCAD: A computer-aided design
% and research tool in photochemistry," Photochem. Photobiol., 68, 141-142, 1998. 
% J. M. Dixon, M. Taniguchi and J. S. Lindsey "PhotochemCAD 2. A refined program with accompanying 
% spectral databases for photochemical calculations", Photochem. Photobiol., 81, 212-213, 2005.
% The spectral absorption measurements were made by J. Li on 12-11-1997 using a Cary 3.
% The reported molar extinction coefficient is from
% Agati, G. and F. Fusi (1990) New trends in photobiology recent advances in bilirubin photophysics.  J.  Photochem. Photobiol. 7, 1-14.

% MetHb.csv:
% Generating S-Nitrosothiols from Hemoglobin
% Roche, Camille J. et al.
% Journal of Biological Chemistry, Volume 288, Issue 31, 22408 - 22425

%%

% chrom_new.methb = table2array(readtable('MetHb.csv'));
% chrom_new.bile(chrom_new.bile(:,2)<0,:) = [];  % delete rows with negative values (omlc data)!

%% Lipid left stitch (exponential)

% wl_left = (340:1:floor(chrom.wl)-1)';

%cftool(chrom.wl(1:21),chrom.lipid(1:21))
% [a,b] = deal(19.24,-0.01267);
% lipid_left = a*exp(b*wl_left);
% 
% chrom_new.lipid = cat(1,lipid_left,chrom.lipid(1:271));
% chrom_new.wl = (340:1:700)';

%% Fit exponentials to Bilirubin and Bile from 430 nm to 550 nm, subtract scattering component and fit difference

% lambda_resampled_bile = 340:1:550;
% ind_chrom = size(lambda_resampled_bile);
% for icount = 1:length(lambda_resampled_bile)
%      [~,ind_chrom(icount)] = min(abs(lambda_resampled_bile(icount) - chrom_new.bile(:,1)));
% end

%cftool(chrom.wl(1:121),chrom.bile(1:121))
% [a,b] = deal(6.971e4,-0.01673);
% bile_left = a*exp(b*lambda_resampled_bile)';
% temp_bile = chrom.bile(1:121) - bile_left(91:end);

%cftool(chrom_new.bile(762:1242,1),chrom_new.bile(762:1242,2))
% [a,b] = deal(9.041e8,-0.02197);
% bile_new_left = a*exp(b*chrom_new.bile(ind_chrom,1));
% temp_bile_new = chrom_new.bile(ind_chrom(91:end),2) - bile_new_left(91:end);

% fun_min_bile = @(c_bile) (temp_bile - c_bile*temp_bile_new);
% options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt','Display','iter','TolFun',1e-10,'StepTolerance',1e-10);
% [opt_ans_bile,resnorm_bile,residual_bile,exitflag_water,output_bile] = lsqnonlin(fun_min_bile,1e-5,2e-4,1e-3,options);  % x0,lb,ub

% abs_bile = chrom_new.bile(ind_chrom(1:90),2) - bile_new_left(1:90);
% stitched_left = cat(1,bile_left(1:90)+opt_ans_bile*abs_bile,chrom.bile(1:271));
% chrom_new.bile = [];
% chrom_new.bile = stitched_left;

%% Interpolate 

% fieldnames = {'hb','hbo2'};
% lambda_resampled_hb_hbo2 = 340:1:700;
% for i = 1:length(fieldnames)
%     chrom_new.(fieldnames{i}) = interp1(chrom_new.(fieldnames{i})(46:226,1),chrom_new.(fieldnames{i})(46:226,2),lambda_resampled_hb_hbo2)';
%     chrom_new.(fieldnames{i})(:,2) = lambda_resampled_hb_hbo2';     
%     chrom_new.(fieldnames{i})(:,[1 2]) = chrom_new.(fieldnames{i})(:,[2 1]);  % swap columns
% end

% lambda_resampled_water = 340:1:700;
% chrom_new.water = interp1(chrom_new.water(21:201,1),chrom_new.water(21:201,2),lambda_resampled_water)';
% chrom_new.water(:,2) = lambda_resampled_water';
% chrom_new.water(:,[1 2]) = chrom_new.water(:,[2 1]);  % swap columns

% lambda_resampled_methb = 374:1:698;
% chrom_new.methb = interp1(chrom_new.methb(:,1),chrom_new.methb(:,2),lambda_resampled_methb)';
% chrom_new.methb(:,2) = lambda_resampled_methb';
% chrom_new.methb(:,[1 2]) = chrom_new.methb(:,[2 1]);  % swap columns

% Set remaining rows
% chrom_new.methb(326,2) = chrom_new.methb(end,2);
% chrom_new.methb(327,2) = chrom_new.methb(end,2);

% chrom_new.methb(326,1) = 699;
% chrom_new.methb(327,1) = 700;

%% Water, Hb, HbO2 fits to omlc.org data, MetHb to .csv file data

% lambda_resampled = 430:1:700;  % maximal fitting range (our data starts at 430 nm and water omlc data goes to 700 nm, methb goes to 698.8540)

% ind_chrom = size(lambda_resampled);
% for icount = 1:length(lambda_resampled)
%    [~,ind_chrom(icount)] = min(abs(lambda_resampled(icount) - chrom_new.water(:,1)));
%    [~,ind_chrom(icount)] = min(abs(lambda_resampled(icount) - chrom_new.hb(:,1)));
%    [~,ind_chrom(icount)] = min(abs(lambda_resampled(icount) - chrom_new.hbo2(:,1)));
%    [~,ind_chrom(icount)] = min(abs(lambda_resampled(icount) - chrom_new.methb(:,1)));
% end

% fun_min_water = @(c_water) (c_water*chrom_new.water(ind_chrom,2) - chrom.water(1:271));
% options = optimoptions('lsqnonlin','Display','iter','TolFun',1e-10,'StepTolerance',1e-10);
% [opt_ans_water,resnorm_water,residual_water,exitflag_water,output_water] = lsqnonlin(fun_min_water,1e-1,1,1e1,options);
 
% fun_min_hb = @(c_hb) (c_hb*chrom_new.hb(ind_chrom,2) - chrom.hb(1:271));
% options = optimoptions('lsqnonlin','Display','iter','TolFun',1e-10,'StepTolerance',1e-10);
% [opt_ans_hb,resnorm_hb,residual_hb,exitflag_hb,output_hb] = lsqnonlin(fun_min_hb,1e-5,1e-3,1,options);
  
% fun_min_hbo2 = @(c_hbo2) (c_hbo2*chrom_new.hbo2(ind_chrom,2) - chrom.hbo2(1:271));
% options = optimoptions('lsqnonlin','Display','iter','TolFun',1e-10,'StepTolerance',1e-10);
% [opt_ans_hbo2,resnorm_hbo2,residual_hbo2,exitflag_hbo2,output_hbo2] = lsqnonlin(fun_min_hbo2,1e-5,1e-3,1,options);

% fun_min_methb = @(c_methb) (c_methb*chrom_new.methb(ind_chrom,2) - chrom.methb(1:271));
% options = optimoptions('lsqnonlin','Display','iter','TolFun',1e-10,'StepTolerance',1e-10);
% [opt_ans_methb,resnorm_methb,residual_methb,exitflag_methb,output_methb] = lsqnonlin(fun_min_methb,1,1e2,1e3,options);

% figure
% semilogy(chrom.wl(1:271), chrom.bile(1:271),'-*');
% hold all
% semilogy(chrom_new.bilirubin(ind_chrom,1),opt_ans_bile*chrom_new.bilirubin(ind_chrom,2));

% semilogy(chrom.wl(1:271), chrom.water(1:271),'-*');
% hold all
% semilogy(chrom_new.water(ind_chrom,1),opt_ans_water*chrom_new.water(ind_chrom,2));

% semilogy(chrom.wl(1:271), chrom.hb(1:271),'-*');
% hold all
% semilogy(chrom_new.hb(ind_chrom,1),opt_ans_hb*chrom_new.hb(ind_chrom,2));

% semilogy(chrom.wl(1:271), chrom.hbo2(1:271),'-*');
% hold all
% semilogy(chrom_new.hbo2(ind_chrom,1),opt_ans_hbo2*chrom_new.hbo2(ind_chrom,2));

% semilogy(chrom.wl(1:271), chrom.methb(1:271),'-*');
% hold all
% semilogy(chrom_new.methb(ind_chrom,1),opt_ans_methb*chrom_new.methb(ind_chrom,2));

%% Finalise new data with fitting parameter and delete wavelength vector [340,700] nm

% chrom_new.hb(:,2) = opt_ans_hb*chrom_new.hb(:,2);
% chrom_new.hb(:,1) = [];

% chrom_new.hbo2(:,2) = opt_ans_hb*chrom_new.hbo2(:,2);
% chrom_new.hbo2(:,1) = [];

% chrom_new.water(:,2) = opt_ans_water*chrom_new.water(:,2);
% chrom_new.water(:,1) = [];

% chrom_new.methb(:,2) = opt_ans_methb*chrom_new.methb(:,2);

% Fit gaussian to left side of methb and stitch

%cftool(chrom_new.methb(1:38,1),chrom_new.methb(1:38,2))
% methb_left = (340:1:chrom_new.methb(1,1)-1)';
% [a1,b1,c1,a2,b2,c2] = deal(279.5,407.7,7.645,1.188e17,2502,368.4);
% methb_left(:,2) = (a1*exp(-((methb_left(:,1)-b1)/c1).^2) + a2*exp(-((methb_left(:,1)-b2)/c2).^2))';
% chrom_new.methb = cat(1,methb_left,chrom_new.methb); 
% chrom_new.methb(:,1) = [];

% save('chromophores.mat','chrom_new')

%% Plotting chromophore spectra

figure
semilogy(chrom_new.wl,chrom_new.hbo2,'linewidth',2);
hold on
% semilogy(chrom.wl,chrom.hbo2,'-*','linewidth',2);
% hold on
semilogy(chrom_new.wl,chrom_new.hb,'linewidth',2);
hold on
% semilogy(chrom.wl,chrom.hb,'-*','linewidth',2);
% hold on
% semilogy(chrom_new.wl,chrom_new.water,'linewidth',2);
% hold on
% semilogy(chrom.wl,chrom.water,'-*','linewidth',2);
% hold on
% semilogy(chrom_new.wl,chrom_new.lipid,'linewidth',2);
% hold on
% semilogy(chrom.wl,chrom.lipid,'-*','linewidth',2);
% hold on
% semilogy(chrom_new.wl,chrom_new.bile,'linewidth',2);
% hold on
% semilogy(chrom.wl,chrom.bile,'-*','linewidth',2);
% hold on
% semilogy(chrom_new.wl,chrom_new.methb,'linewidth',2);
% hold on
% semilogy(chrom.wl,chrom.methb,'-*','linewidth',2);
% hold on
xlabel('\lambda [nm]','FontSize',17)
ylabel(' \mu_a [cm^{-1}]','FontSize',17)
title('Chromophore Spectra','FontSize',20);
legend('Hbo2','Hb','FontSize',14);
set(gca,'FontSize',20)
xlim([340 700])

clear