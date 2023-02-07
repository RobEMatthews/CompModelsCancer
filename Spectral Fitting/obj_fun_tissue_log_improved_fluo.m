function [diff_vector] = obj_fun_tissue_log_improved_fluo(crit_angle,g,mua_ex_lin,mus_ex_lin,mua_em_lin,mus_em_lin,PHD_int,F_exp_sampled,calib_fact,mua_em_tissue,musdash_tissue,F_resp,lambda_fluop,R,f_blood,sato2,percentage_conc_water_lipid,f_lipid,f_bile,f_methb,a_dash,fray,bMie,f_coll,f_elas,f_flav,f_lipop,f_nadh,f_pirod,f_porp,f_trypt)

global diff_vector
global Fluo_generated

% Key functions
mus_em_val = musdash_tissue(a_dash,fray,bMie);
mua_em_val = mua_em_tissue(R,f_blood,sato2,percentage_conc_water_lipid,f_lipid,f_bile,f_methb);
F_resp_val = F_resp(f_coll,f_elas,f_flav,f_lipop,f_nadh,f_pirod,f_porp,f_trypt);

% Remain within corners of LUT
mua_em_val(mua_em_val>max(max(mua_em_lin))) = max(max(mua_em_lin));
mus_em_val(mus_em_val>max(max(mus_em_lin))) = max(max(mus_em_lin));

mua_em_val(mua_em_val<min(min(mua_em_lin))) = min(min(mua_em_lin));
mus_em_val(mus_em_val<min(min(mus_em_lin))) = min(min(mus_em_lin));

% Interpolate using just emission optical property profiles (excitation is
% value at 340 nm)

for i=81:121
 Fluo_generated(i-80) = interpn(log10(unique(mua_ex_lin)),log10(unique(mus_ex_lin)),log10(unique(mua_em_lin)),log10(unique(mus_em_lin)),log10(PHD_int),mua_em_val(1),mus_em_val(1)-log10(1-g),mua_em_val(i),mus_em_val(i)-log10(1-g),'linear');  % key interpolation
end

Fluo_convolv = (F_resp_val.*crit_angle).*Fluo_generated';  % convolve with fluorescence response and critical angle at tissue interface
diff_vector = (calib_fact.*F_exp_sampled - (10.^Fluo_convolv));
%diff_vector = 100.*abs(calib_fact.*F_exp_sampled - (10.^Fluo_convolv)) ./ (calib_fact.*F_exp_sampled);

subplot(2,1,1)
plot(lambda_fluop,(10.^Fluo_convolv),'-r','linewidth',5);
hold on
plot(lambda_fluop,calib_fact.*F_exp_sampled,'-*');
xlabel('\lambda  [nm]')
ylabel('fluorescence')
legend('F inverse MC','F experimental');
hold off

subplot(2,1,2)
plot(lambda_fluop,100.*abs(calib_fact.*F_exp_sampled-(10.^Fluo_convolv'))./(calib_fact.*F_exp_sampled),'-r');
xlabel('\lambda  [nm]')
ylabel('Error in % ')

drawnow;
norm_val=sum(abs(diff_vector).^2);

end