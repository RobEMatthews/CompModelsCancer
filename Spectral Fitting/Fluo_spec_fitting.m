%% Preprocess/Read in data

% Create structure containing all patient spectra
% fieldnames = {'fl_fat_45','fl_fat_46','fl_fat_47','fl_fat_48','fl_fat_49','fl_fat_50',...
%     'fl_lmc_45','fl_lmc_46','fl_lmc_47','fl_lmc_48','fl_lmc_49','fl_lmc_50',...
%     'fl_lml_45','fl_lml_46','fl_lml_47','fl_lml_48','fl_lml_49','fl_lml_50','fl_msw_47',...
%     'fl_tmb_45','fl_tmb_46','fl_tmb_49','fl_tum_45','fl_tum_46','fl_tum_47','fl_tum_48','fl_tum_50'};
% fieldvars = {fl_fat_45,fl_fat_46,fl_fat_47,fl_fat_48,fl_fat_49,fl_fat_50,...
%     fl_lmc_45,fl_lmc_46,fl_lmc_47,fl_lmc_48,fl_lmc_49,fl_lmc_50,...
%     fl_lml_45,fl_lml_46,fl_lml_47,fl_lml_48,fl_lml_49,fl_lml_50,fl_msw_47,...
%     fl_tmb_45,fl_tmb_46,fl_tmb_49,fl_tum_45,fl_tum_46,fl_tum_47,fl_tum_48,fl_tum_50};
% for i = 1:length(fieldnames)
%     patient_spectra.(fieldnames{i}) = fieldvars{i};
% end

% save('patient_spectra.mat','patient_spectra')

% Read in LUT database
% LUT = readmatrix('LUT.txt');
 
% calibration_dataset.mua_ex_lin = LUT(:,1);
% calibration_dataset.mus_ex_lin = LUT(:,2);
% calibration_dataset.mua_em_lin = LUT(:,3);
% calibration_dataset.mus_em_lin = LUT(:,4);
% calibration_dataset.PHD_int = LUT(:,5);
% calibration_dataset.refrac_water = refrac_water;

% clear LUT

calibration_dataset.calib_fact = 1.5e-5;  % theoretical optimal: 1.4 x 10^-4
% save('calibration_dataset.mat','calibration_dataset')

%% Load data

% close all
% clear

% load('wavelen.mat')
% load('chromophores.mat')
% load('fluorophores.mat')
% load('calibration_dataset.mat')
% load('patient_spectra.mat')

g = 0.9;  % anisotropy factor
n_air = 1;  % refractive index of air assumed constant

calib_fact = calibration_dataset.calib_fact;
mua_em_lin = calibration_dataset.mua_em_lin;
mus_em_lin = calibration_dataset.mus_em_lin;
mua_ex_lin = calibration_dataset.mua_ex_lin;
mus_ex_lin = calibration_dataset.mus_ex_lin;
PHD_int = calibration_dataset.PHD_int;

PHD_int = permute(reshape(PHD_int,10,20,10,20),[4,3,2,1]);  % correctly format integrated value into 4D array

% Interpolate refractive index of water data
refrac_water = calibration_dataset.refrac_water;
lambda_resampled_refrac_water = 200:1:7000;
    temp = interp1(refrac_water(:,1),refrac_water(:,2),lambda_resampled_refrac_water)';
refrac_water = [];
refrac_water(:,1) = lambda_resampled_refrac_water;
refrac_water(:,2) = temp;

global diff_vector

global Fluo_generated

%% Resample experimental wavelength

lambda_exp = wavelen;
lambda_resampled = 500:2:580;  % choose wavelength vector to optimize the fitting time

ind_exp = size(lambda_resampled);
for icount = 1:length(lambda_resampled)
    [~,ind_exp(icount)] = min(abs(lambda_resampled(icount) - lambda_exp)); 
    % take indices of experimental wavelength vector/array at resampled wavelengths
end

lamda_exp_sampled = lambda_exp(ind_exp);

%% Resample water refractive index and calculate log of critical angle

ind_refrac = size(lambda_resampled);
for icount = 1:length(lambda_resampled)
    [~,ind_refrac(icount)] = min(abs(lambda_resampled(icount) - refrac_water(:,1))); 
    % take indices of refractive water index vector/array at resampled wavelengths
end

refrac_water_sampled = refrac_water(ind_refrac,2);
crit_angle = log10(asind(n_air./refrac_water_sampled));

%% Resample chromophore data

lambda_resampled_chrom = 340:2:700;
ind_chrom = size(lambda_resampled_chrom);
for icount = 1:length(lambda_resampled_chrom)
    [~,ind_chrom(icount)] = min(abs(lambda_resampled_chrom(icount) - chrom.wl));
    % take indices of the chromophore concentration vectors/arrays at resampled wavelengths
end

lambda_chrom = chrom.wl(ind_chrom);  % resampled wavelengths
mua_hbo2 = chrom.hbo2(ind_chrom);  % resampled mua of HbO2
mua_hb = chrom.hb(ind_chrom);
mua_water = chrom.water(ind_chrom);
mua_lipid = chrom.lipid(ind_chrom);
mua_bile = chrom.bile(ind_chrom);
mua_methb = chrom.methb(ind_chrom);

%% Resample fluorophore data

ind_fluop = size(lambda_resampled);
for icount = 1:length(lambda_resampled)
    %[~,ind_fluop(icount)] = min(abs(lambda_resampled(icount) - fluop.wl_ex));
    [~,ind_fluop(icount)] = min(abs(lambda_resampled(icount) - fluop.wl_em));
    % take indices of the fluorophore concentration vectors/arrays at resampled wavelengths
end

lambda_fluop = fluop.wl_em(ind_fluop);  % resampled emission wavelengths
collagen_em = fluop.collagen_em(ind_fluop);
elastin_em = fluop.elastin_em(ind_fluop);
flavin_em = fluop.flavin_em(ind_fluop);
lipopigm_em = fluop.lipopigm_em(ind_fluop);
nadh_em = fluop.nadh_em(ind_fluop);
pirodox_em = fluop.pirodox_em(ind_fluop);
porphyr_em = fluop.porphyr_em(ind_fluop);
tryptop_em = fluop.tryptop_em(ind_fluop);

% Muliply by excitation fluorophore concentration at 340 nm to get
% concentration response
F_coll_resp = collagen_em * fluop.collagen_ex(281);
F_elas_resp = elastin_em * fluop.elastin_ex(281);
F_flav_resp = flavin_em * fluop.flavin_ex(281);
F_lipop_resp = lipopigm_em * fluop.lipopigm_ex(281);
F_nadh_resp = nadh_em * fluop.nadh_ex(281);
F_pirod_resp = pirodox_em * fluop.pirodox_ex(281);
F_porp_resp = porphyr_em * fluop.porphyr_ex(281);
F_trypt_resp = tryptop_em * fluop.tryptop_ex(281);

%% Plotting chromophore spectra

% figure
% semilogy(lambda_chrom,mua_hbo2,'-*');
% hold on
% semilogy(lambda_chrom,mua_hb,'-*');
% hold on
% semilogy(lambda_chrom,mua_water,'-*');
% hold on
% plot(lambda_chrom,mua_lipid,'-*');
% hold on
% semilogy(lambda_chrom,mua_bile,'-*');
% hold on
% semilogy(lambda_chrom,mua_methb,'-*');
% xlabel('\lambda [nm]')
% ylabel(' \mu_a [cm^{-1}]')
% title('Resampled chromophore spectra');
% legend('Hbo2','Hb','Water','Lipid','Bile','MetHb');

%% Plotting fluorophore response spectra

% figure
% semilogy(lambda_fluop,F_coll_resp,'-*');
% hold on
% semilogy(lambda_fluop,F_elas_resp,'-*');
% hold on
% semilogy(lambda_fluop,F_flav_resp,'-*');
% hold on
% semilogy(lambda_fluop,F_lipop_resp,'-*');
% hold on
% semilogy(lambda_fluop,F_nadh_resp,'-*');
% hold on
% semilogy(lambda_fluop,F_pirod_resp,'-*');
% hold on
% semilogy(lambda_fluop,F_porp_resp,'-*');
% hold on
% semilogy(lambda_fluop,F_trypt_resp,'-*');
% xlabel('\lambda [nm]')
% ylabel('F_{Response} [-]')
% title('Resampled fluorophore spectra');
% legend('Collagen','Elastin','Flavin','Lipopigm','Nadh','Pirodox','Porphyr','Tryptop');

%% Choose fluorescence spectra

Fitting_spectra = reshape(patient_spectra.fl_fat_45(4,3,:),[1,1003]);

for icount = 1:size(Fitting_spectra,1)

diff_vector = 0;  % difference vector to be minimized

F_exp = Fitting_spectra;  % experimental spectra to be fitted

%% Plot experimental measurements

% figure
% plot(lambda_exp,F_exp)
% xlabel('\lambda [nm]')
% ylabel(' F [-]')
% title('Raw experimental spectrum');

F_exp_sampled = F_exp(ind_exp);  % resampled fluorescence spectrum (experimental 
% fluorescence spectrum with lower wavelength resolution)
% 
% %% Check interpolated fluorescence
% 
% figure
% plot(lamda_exp_sampled,F_exp_sampled,'-*')
% xlabel('\lambda [nm]')
% ylabel(' F [-]')
% title('Resampled experimental spectrum');

%% Define mus_dash, mua_em, mus_em and F_resp parametrically 

musdash_tissue = @(a_dash,fray,bMie) (log10(a_dash.*(fray.*((lambda_chrom)./(500)).^(-4) + (1-fray).*((lambda_chrom)./(500)).^(-bMie))*(1-g)));
mua_em_tissue = @(R,f_blood,sato2,percentage_conc_water_lipid,f_lipid,f_bile,f_methb) (log10(((1-exp(-2.*R.*(sato2.*mua_hbo2+(1-sato2).*mua_hb)))./(2.*R.*(sato2.*mua_hbo2+(1-sato2).*mua_hb))).*f_blood.*(sato2.*mua_hbo2+(1-sato2).*mua_hb)+percentage_conc_water_lipid.*(f_lipid.*mua_lipid + (1-f_lipid).*mua_water)+f_bile.*mua_bile+f_methb.*mua_methb));
F_resp = @(f_coll,f_elas,f_flav,f_lipop,f_nadh,f_pirod,f_porp,f_trypt) (log10(F_coll_resp.*f_coll + F_elas_resp.*f_elas + F_flav_resp.*f_flav + F_lipop_resp.*f_lipop + F_nadh_resp.*f_nadh + F_pirod_resp.*f_pirod + F_porp_resp.*f_porp + F_trypt_resp.*f_trypt));

fun_min = @(argu) (obj_fun_tissue_log_improved_fluo(crit_angle,g,mua_ex_lin,mus_ex_lin,mua_em_lin,mus_em_lin,PHD_int,F_exp_sampled,calib_fact,mua_em_tissue,musdash_tissue,F_resp,lambda_fluop,argu(1),argu(2),argu(3),argu(4),argu(5),argu(6),argu(7),argu(8),argu(9),argu(10),argu(11),argu(12),argu(13),argu(14),argu(15),argu(16),argu(17),argu(18)));  % to make the extra parameters for optimisation implicit

% argu = R,f_blood,sato2,percentage_conc_water_lipid,f_lipid,f_bile,f_methb,a_dash,fray,bMie,f_coll,f_elas,f_flav,f_lipop,f_nadh,f_pirod,f_porp,f_trypt

%% Optimization 

options = optimoptions('lsqnonlin','Display','iter','TolFun',1e-10,'StepTolerance',1e-10,'Algorithm','levenberg-marquardt');
x0=[25e-6,0.05,0.64,2,1,0.01,0.01,40,0.56,0.1,0.01,2,0.1,0.09,0.2,0.1,0,0];  % initial fitting parameters
lb=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];  % lower boundary of fitting parameters
ub=[0.01,1,1,10,5,3,0.1,100,1,4,10,10,10,10,10,10,0,0];  % upper boundary of fitting parameters
[opt_ans,resnorm,residual,exitflag,output] = lsqnonlin(fun_min,x0,lb,ub,options);

average_vessel_radius=opt_ans(1);
f_blood=opt_ans(2);
sato2=opt_ans(3);
percentage_conc_water_lipid=opt_ans(4);
f_lipid=opt_ans(5);
f_bile=opt_ans(6);
f_methb=opt_ans(7);
a_dash=opt_ans(8);
fray=opt_ans(9);
bMie=opt_ans(10);
f_coll=opt_ans(11);
f_elas=opt_ans(12);
f_flav=opt_ans(13);
f_lipop=opt_ans(14);
f_nadh=opt_ans(15);
f_pirod=opt_ans(16);
f_porp=opt_ans(17);
f_trypt=opt_ans(18);





% % 
% % 
% % 
% % f_hbo2=sato2*f_blood;
% % f_hb=(1-sato2)*f_blood;
% % f_water=1-f_lipid;
% % 
% % blood_packaging_factor=((1-exp(-2.*average_vessel_radius.*(sato2.*mua_hbo2+(1-sato2).*mua_hb)))./(2.*average_vessel_radius.*(sato2.*mua_hbo2+(1-sato2).*mua_hb)));
% % average_blood_packaging_factor=mean((1-exp(-2.*average_vessel_radius.*(sato2.*mua_hbo2+(1-sato2).*mua_hb)))./(2.*average_vessel_radius.*(sato2.*mua_hbo2+(1-sato2).*mua_hb)));
% % f_hbo2_homogenous=sato2*f_blood*average_blood_packaging_factor;
% % f_hb_homogenous=(1-sato2)*f_blood*average_blood_packaging_factor;
% % f_lipid_homogenous=percentage_conc_water_lipid*f_lipid;
% % f_water_homogenous=percentage_conc_water_lipid*(1-f_lipid);
% % 
% % sprintf('Inverse MC parameters are : \n f_hbo2=%f\n f_hb=%f\n f_water=%f\n f_lipid=%f\n f_bile=%f\n f_methb=%f\n a_dash=%f\n fray=%f\n bMie=%f\n',f_hbo2,f_hb,f_water,f_lipid,f_bile,f_methb,a_dash,fray,bMie)
% % 
% % 
% % Plotting tissue optical properties
% % 
% % figure
% % plot(lambda_chrom,musdash_tissue(a_dash,fray,bMie),'-*');
% % xlabel('\lambda [nm]')
% % ylabel('\mu_s^{red} [cm^{-1}]')
% % title('Inverse Monte Carlo \mu_s^{red}')
% % 
% % figure
% % plot(lambda_chrom,mua_tissue(f_hbo2,f_hb,f_water,f_lipid,f_bile),'-*');
% % xlabel('\lambda [nm]')
% % ylabel('\mu_a [cm^{-1}]')
% % title('Inverse Monte Carlo \mu_a')
% % 
% % 
% % tissue_musdash_log10 = musdash_tissue(a_dash,fray,bMie);
% % tissue_mua_log10 = mua_tissue(average_vessel_radius,f_blood,sato2,percentage_conc_water_lipid,f_lipid,f_bile,f_methb);
% % 
% % tissue_musdash=10.^tissue_musdash_log10;
% % tissue_mua=10.^tissue_mua_log10;

end