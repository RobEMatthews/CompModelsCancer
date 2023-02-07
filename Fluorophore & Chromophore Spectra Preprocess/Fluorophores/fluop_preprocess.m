%% FLUOROPHORE PREPROCESS

% Interpolation, Gaussian fitting, stitching and data smoothening
% (Savitzky-Golay filtering) used to preprocess original fluorophore
% spectra (fluorophores_spec2.mat) to allow for extension over full
% excitation and emission wavelength ranges.

% ORIGINAL DATA (fluorophores_spec2.mat)
% Wagnieres, G.A., Star, W.M. and Wilson, B.C. (1998), In Vivo Fluorescence Spectroscopy and Imaging for Oncological Applications. Photochemistry and Photobiology, 68: 603-632. https://doi.org/10.1111/j.1751-1097.1998.tb02521.x

%% Normalize fluorophores in excitation and emission by maximum of one fluorophore

load('fluorophore_spec2.mat')

fluop.collagen_ex(:,2) = fluop.collagen_ex(:,2)/max(fluop.lipopigm_ex(:,2));
fluop.elastin_ex(:,2) = fluop.elastin_ex(:,2)/max(fluop.lipopigm_ex(:,2));
fluop.flavin_ex(:,2) = fluop.flavin_ex(:,2)/max(fluop.lipopigm_ex(:,2));
fluop.lipopigm_ex(:,2) = fluop.lipopigm_ex(:,2)/max(fluop.lipopigm_ex(:,2));
fluop.nadh_ex(:,2) = fluop.nadh_ex(:,2)/max(fluop.lipopigm_ex(:,2));
fluop.pirodox_ex(:,2) = fluop.pirodox_ex(:,2)/max(fluop.lipopigm_ex(:,2));
fluop.porphyr_ex(:,2) = fluop.porphyr_ex(:,2)/max(fluop.lipopigm_ex(:,2));
fluop.tryptop_ex(:,2) = fluop.tryptop_ex(:,2)/max(fluop.lipopigm_ex(:,2));

fluop.collagen_em(:,2) = fluop.collagen_em(:,2)/max(fluop.pirodox_em(:,2));
fluop.elastin_em(:,2) = fluop.elastin_em(:,2)/max(fluop.pirodox_em(:,2));
fluop.flavin_em(:,2) = fluop.flavin_em(:,2)/max(fluop.pirodox_em(:,2));
fluop.lipopigm_em(:,2) = fluop.lipopigm_em(:,2)/max(fluop.pirodox_em(:,2));
fluop.nadh_em(:,2) = fluop.nadh_em(:,2)/max(fluop.pirodox_em(:,2));
fluop.pirodox_em(:,2) = fluop.pirodox_em(:,2)/max(fluop.pirodox_em(:,2));
fluop.porphyr_em(:,2) = fluop.porphyr_em(:,2)/max(fluop.pirodox_em(:,2));
fluop.tryptop_em(:,2) = fluop.tryptop_em(:,2)/max(fluop.pirodox_em(:,2));

%% Average duplicated values

fieldnames_ex = {'collagen_ex','elastin_ex','flavin_ex','lipopigm_ex','nadh_ex','pirodox_ex','porphyr_ex','tryptop_ex'};
fieldnames_em = {'collagen_em','elastin_em','flavin_em','lipopigm_em','nadh_em','pirodox_em','porphyr_em','tryptop_em'};
for i = 1:length(fieldnames_ex)
    [wl_ex_vals,~,index] = unique(fluop.(fieldnames_ex{i})(:,1));
    fluo_ex_vals = accumarray(index,fluop.(fieldnames_ex{i})(:,2),[],@mean);
    fluop.(fieldnames_ex{i}) = [wl_ex_vals,fluo_ex_vals]; 
    
    [wl_em_vals,~,index] = unique(fluop.(fieldnames_em{i})(:,1));
    fluo_em_vals = accumarray(index,fluop.(fieldnames_em{i})(:,2),[],@mean);
    fluop.(fieldnames_em{i}) = [wl_em_vals,fluo_em_vals]; 
end

%% Create linspaces for each fluorophore (excitation and emission), fit Gaussians (or superpositions) and stitch

% EXCITATION (WAVELENGTHS 200:0.5:500)

%cftool(fluop.collagen_ex(1:43,1),fluop.collagen_ex(1:43,2))
collagen_ex_left = (200:0.5:floor(fluop.collagen_ex(1,1)))';
[a1,b1,c1,a2,b2,c2,a3,b3,c3,a4,b4,c4] = deal(0.4257,334.5,23.21,-0.03298,329.7,10.8,-0.01438,320.1,6.821,0.122,303.3,169.7);
collagen_ex_left(:,2) = (a1*exp(-((collagen_ex_left(:,1)-b1)/c1).^2) + a2*exp(-((collagen_ex_left(:,1)-b2)/c2).^2) + a3*exp(-((collagen_ex_left(:,1)-b3)/c3).^2) + a4*exp(-((collagen_ex_left(:,1)-b4)/c4).^2))';

%cftool(fluop.collagen_ex(43:end,1),fluop.collagen_ex(43:end,2))
collagen_ex_right = (ceil(fluop.collagen_ex(end,1)):0.5:500)';
[a1,b1,c1] = deal(0.5179,335.4,27.85);
collagen_ex_right(:,2) = (a1*exp(-((collagen_ex_right(:,1)-b1)/c1).^2))';

fluop.collagen_ex = cat(1,collagen_ex_left,fluop.collagen_ex,collagen_ex_right);
fluop.collagen_ex = (abs(fluop.collagen_ex) + fluop.collagen_ex)/2;

%cftool(fluop.elastin_ex(1:33,1),fluop.elastin_ex(1:33,2))
elastin_ex_left = (200:0.5:floor(fluop.elastin_ex(1,1)))';
[a1,b1,c1,a2,b2,c2,a3,b3,c3,a4,b4,c4,a5,b5,c5] = deal(0.02881,331.2,8.336,0.08425,342.2,17.96,-0.0008165,326.8,2.105,0.117,344.4,20.92,2.463,586.7,168.1);
elastin_ex_left(:,2) = (a1*exp(-((elastin_ex_left(:,1)-b1)/c1).^2) + a2*exp(-((elastin_ex_left(:,1)-b2)/c2).^2) + a3*exp(-((elastin_ex_left(:,1)-b3)/c3).^2) + a4*exp(-((elastin_ex_left(:,1)-b4)/c4).^2) + a5*exp(-((elastin_ex_left(:,1)-b5)/c5).^2))';

%cftool(fluop.elastin_ex(34:end,1),fluop.elastin_ex(34:end,2))
elastin_ex_right = (ceil(fluop.elastin_ex(end,1)):0.5:500)';
[a1,b1,c1,a2,b2,c2] = deal(0.0779,353.6,25.88,0.4586,338.3,57.87);
elastin_ex_right(:,2) = (a1*exp(-((elastin_ex_right(:,1)-b1)/c1).^2) + a2*exp(-((elastin_ex_right(:,1)-b2)/c2).^2))';

fluop.elastin_ex = cat(1,elastin_ex_left,fluop.elastin_ex,elastin_ex_right);
fluop.elastin_ex = (abs(fluop.elastin_ex) + fluop.elastin_ex)/2;

%cftool(fluop.flavin_ex(1:9,1),fluop.flavin_ex(1:9,2))
flavin_ex_left = (200:0.5:floor(fluop.flavin_ex(1,1)))';
[a1,b1,c1,a2,b2,c2] = deal(0.4751,215.2,21.34,0.04942,202.3,7.012);
flavin_ex_left(:,2) = (a1*exp(-((flavin_ex_left(:,1)-b1)/c1).^2) + a2*exp(-((flavin_ex_left(:,1)-b2)/c2).^2))';

%cftool(fluop.flavin_ex(135:end,1),fluop.flavin_ex(135:end,2))
flavin_ex_right = (ceil(fluop.flavin_ex(end,1)):0.5:500)';
[a1,b1,c1,a2,b2,c2,a3,b3,c3] = deal(0.2693,449.6,34.15,-0.07825,514.6,14.68,0.06951,477.1,21.49);
flavin_ex_right(:,2) = (a1*exp(-((flavin_ex_right(:,1)-b1)/c1).^2) + a2*exp(-((flavin_ex_right(:,1)-b2)/c2).^2) + a3*exp(-((flavin_ex_right(:,1)-b3)/c3).^2))';

fluop.flavin_ex = cat(1,flavin_ex_left,fluop.flavin_ex,flavin_ex_right);
fluop.flavin_ex = (abs(fluop.flavin_ex) + fluop.flavin_ex)/2;

%cftool(fluop.lipopigm_ex(1:71,1),fluop.lipopigm_ex(1:71,2))
lipopigm_ex_left = (200:0.5:floor(fluop.lipopigm_ex(1,1)))';
[a1,b1,c1,a2,b2,c2,a3,b3,c3] = deal(1.003,343,23.55,-0.004689,336.2,0.2509,0.251,321,10.44);
lipopigm_ex_left(:,2) = (a1*exp(-((lipopigm_ex_left(:,1)-b1)/c1).^2) + a2*exp(-((lipopigm_ex_left(:,1)-b2)/c2).^2) + a3*exp(-((lipopigm_ex_left(:,1)-b3)/c3).^2))';

%cftool(fluop.lipopigm_ex(132:end,1),fluop.lipopigm_ex(132:end,2))
lipopigm_ex_right = (ceil(fluop.lipopigm_ex(end,1)):0.5:500)';
[a1,b1,c1,a2,b2,c2,a3,b3,c3] = deal(0.09807,423.8,7.767,0.01972,434.6,5.083,0.3667,437.8,29.1);
lipopigm_ex_right(:,2) = (a1*exp(-((lipopigm_ex_right(:,1)-b1)/c1).^2) + a2*exp(-((lipopigm_ex_right(:,1)-b2)/c2).^2) + a3*exp(-((lipopigm_ex_right(:,1)-b3)/c3).^2))';

fluop.lipopigm_ex = cat(1,lipopigm_ex_left,fluop.lipopigm_ex,lipopigm_ex_right);
fluop.lipopigm_ex = (abs(fluop.lipopigm_ex) + fluop.lipopigm_ex)/2;

%cftool(fluop.nadh_ex(1:19,1),fluop.nadh_ex(1:19,2))
nadh_ex_left = (200:0.5:floor(fluop.nadh_ex(1,1)))';
[a1,b1,c1,a2,b2,c2] = deal(0.6878,261.8,11.96,0.1119,255,3.626);
nadh_ex_left(:,2) = (a1*exp(-((nadh_ex_left(:,1)-b1)/c1).^2) + a2*exp(-((nadh_ex_left(:,1)-b2)/c2).^2))';

%cftool(fluop.nadh_ex(83:end,1),fluop.nadh_ex(83:end,2))
nadh_ex_right = (ceil(fluop.nadh_ex(end,1)):0.5:500)';
[a1,b1,c1,a2,b2,c2] = deal(0.1092,344,12.19,0.2716,360.1,22.88);
nadh_ex_right(:,2) = (a1*exp(-((nadh_ex_right(:,1)-b1)/c1).^2) + a2*exp(-((nadh_ex_right(:,1)-b2)/c2).^2))';

fluop.nadh_ex = cat(1,nadh_ex_left,fluop.nadh_ex,nadh_ex_right);
fluop.nadh_ex = (abs(fluop.nadh_ex) + fluop.nadh_ex)/2;

%cftool(fluop.pirodox_ex(1:50,1),fluop.pirodox_ex(1:50,2))
pirodox_ex_left = (200:0.5:floor(fluop.pirodox_ex(1,1)))';
[a1,b1,c1,a2,b2,c2,a3,b3,c3] = deal(-0.01706,298.7,5.16,0.7691,329.8,43.17,0.103,299.2,12.13);
pirodox_ex_left(:,2) = (a1*exp(-((pirodox_ex_left(:,1)-b1)/c1).^2) + a2*exp(-((pirodox_ex_left(:,1)-b2)/c2).^2) + a3*exp(-((pirodox_ex_left(:,1)-b3)/c3).^2))';

%cftool(fluop.pirodox_ex(50:end,1),fluop.pirodox_ex(50:end,2))
pirodox_ex_right = (ceil(fluop.pirodox_ex(end,1)):0.5:500)';
[a1,b1,c1] = deal(0.6672,312.9,22.78);
pirodox_ex_right(:,2) = (a1*exp(-((pirodox_ex_right(:,1)-b1)/c1).^2))';

fluop.pirodox_ex = cat(1,pirodox_ex_left,fluop.pirodox_ex,pirodox_ex_right);
fluop.pirodox_ex = (abs(fluop.pirodox_ex) + fluop.pirodox_ex)/2;

%cftool(fluop.porphyr_ex(1:62,1),fluop.porphyr_ex(1:62,2))
porphyr_ex_left = (200:0.5:floor(fluop.porphyr_ex(1,1)))';
[a1,b1,c1,a2,b2,c2] = deal(-0.01793,378,8.811,0.9562,421.2,67.27);
porphyr_ex_left(:,2) = (a1*exp(-((porphyr_ex_left(:,1)-b1)/c1).^2) + a2*exp(-((porphyr_ex_left(:,1)-b2)/c2).^2))';

%cftool(fluop.porphyr_ex(116:end,1),fluop.porphyr_ex(116:end,2))
porphyr_ex_right = (ceil(fluop.porphyr_ex(end,1)):0.5:500)';
[a1,b1,c1] = deal(0.2104,532.5,47.06);
porphyr_ex_right(:,2) = (a1*exp(-((porphyr_ex_right(:,1)-b1)/c1).^2))';

fluop.porphyr_ex = cat(1,porphyr_ex_left,fluop.porphyr_ex,porphyr_ex_right);
fluop.porphyr_ex = (abs(fluop.porphyr_ex) + fluop.porphyr_ex)/2;

%cftool(fluop.tryptop_ex(1:20,1),fluop.tryptop_ex(1:20,2))
tryptop_ex_left = (200:0.5:floor(fluop.tryptop_ex(1,1)))';
[a1,b1,c1,a2,b2,c2,a3,b3,c3] = deal(0.7141,215.8,20.61,0.0395,208.7,4.903,0.02541,204,2.034);
tryptop_ex_left(:,2) = (a1*exp(-((tryptop_ex_left(:,1)-b1)/c1).^2) + a2*exp(-((tryptop_ex_left(:,1)-b2)/c2).^2) + a3*exp(-((tryptop_ex_left(:,1)-b3)/c3).^2))';

%cftool(fluop.tryptop_ex(73:end,1),fluop.tryptop_ex(73:end,2))
tryptop_ex_right = (ceil(fluop.tryptop_ex(end,1)):0.5:500)';
[a1,b1,c1,a2,b2,c2,a3,b3,c3] = deal(0.303,272.6,9.205,0.1173,282.6,5.998,0.1481,290.1,10.05);
tryptop_ex_right(:,2) = (a1*exp(-((tryptop_ex_right(:,1)-b1)/c1).^2) + a2*exp(-((tryptop_ex_right(:,1)-b2)/c2).^2) + a3*exp(-((tryptop_ex_right(:,1)-b3)/c3).^2))';

fluop.tryptop_ex = cat(1,tryptop_ex_left,fluop.tryptop_ex,tryptop_ex_right);
fluop.tryptop_ex = (abs(fluop.tryptop_ex) + fluop.tryptop_ex)/2;

% EMISSION (WAVELENGTHS: 300:0.5:700)

%cftool(fluop.collagen_em(1:30,1),fluop.collagen_em(1:30,2))
collagen_em_left = (300:0.5:floor(fluop.collagen_em(1,1)))';
[a1,b1,c1,a2,b2,c2] = deal(0.3525,402.7,31.52,0.2143,366.1,34.16);
collagen_em_left(:,2) = (a1*exp(-((collagen_em_left(:,1)-b1)/c1).^2) + a2*exp(-((collagen_em_left(:,1)-b2)/c2).^2))';

%cftool(fluop.collagen_em(30:end,1),fluop.collagen_em(30:end,2))
collagen_em_right = (ceil(fluop.collagen_em(end,1)):0.5:700)';
[a1,b1,c1,a2,b2,c2,a3,b3,c3] = deal(0.2158,380.5,27.95,0.328,413.7,42.88,0.07416,465.4,32.05);
collagen_em_right(:,2) = (a1*exp(-((collagen_em_right(:,1)-b1)/c1).^2) + a2*exp(-((collagen_em_right(:,1)-b2)/c2).^2) + a3*exp(-((collagen_em_right(:,1)-b3)/c3).^2))';

fluop.collagen_em = cat(1,collagen_em_left,fluop.collagen_em,collagen_em_right);
fluop.collagen_em = (abs(fluop.collagen_em) + fluop.collagen_em)/2;

%cftool(fluop.elastin_em(1:33,1),fluop.elastin_em(1:33,2))
elastin_em_left = (300:0.5:floor(fluop.elastin_em(1,1)))';
[a1,b1,c1] = deal(0.4849,439.1,77.06);
elastin_em_left(:,2) = (a1*exp(-((elastin_em_left(:,1)-b1)/c1).^2))';

%cftool(fluop.elastin_em(33:end,1),fluop.elastin_em(33:end,2))
elastin_em_right = (ceil(fluop.elastin_em(end,1)):0.5:700)';
[a1,b1,c1,a2,b2,c2] = deal(0.2531,392.1,51.47,0.2672,444.7,69.41);
elastin_em_right(:,2) = (a1*exp(-((elastin_em_right(:,1)-b1)/c1).^2) + a2*exp(-((elastin_em_right(:,1)-b2)/c2).^2))';

fluop.elastin_em = cat(1,elastin_em_left,fluop.elastin_em,elastin_em_right);
fluop.elastin_em = (abs(fluop.elastin_em) + fluop.elastin_em)/2;

%cftool(fluop.flavin_em(1:31,1),fluop.flavin_em(1:31,2))
flavin_em_left = (300:0.5:floor(fluop.flavin_em(1,1)))';
[a1,b1,c1,a2,b2,c2,a3,b3,c3,a4,b4,c4] = deal(0.3851,553.8,21.62,0.004659,541.2,2.667,0.1369,534.4,13.67,0.1288,520.1,14.7);
flavin_em_left(:,2) = (a1*exp(-((flavin_em_left(:,1)-b1)/c1).^2) + a2*exp(-((flavin_em_left(:,1)-b2)/c2).^2) + a3*exp(-((flavin_em_left(:,1)-b3)/c3).^2) + a4*exp(-((flavin_em_left(:,1)-b4)/c4).^2))';

%cftool(fluop.flavin_em(31:end,1),fluop.flavin_em(31:end,2))
flavin_em_right = (ceil(fluop.flavin_em(end,1)):0.5:700)';
[a1,b1,c1,a2,b2,c2,a3,b3,c3] = deal(0.05049,621.7,48.5,0.1613,555.5,25.43,0.4023,474.6,110.6);
flavin_em_right(:,2) = (a1*exp(-((flavin_em_right(:,1)-b1)/c1).^2) + a2*exp(-((flavin_em_right(:,1)-b2)/c2).^2) + a3*exp(-((flavin_em_right(:,1)-b3)/c3).^2))';

fluop.flavin_em = cat(1,flavin_em_left,fluop.flavin_em,flavin_em_right);
fluop.flavin_em = (abs(fluop.flavin_em) + fluop.flavin_em)/2;

%cftool(fluop.lipopigm_em(1:69,1),fluop.lipopigm_em(1:69,2))
lipopigm_em_left = (300:0.5:floor(fluop.lipopigm_em(1,1)))';
[a1,b1,c1,a2,b2,c2,a3,b3,c3,a4,b4,c4] = deal(0.8554,566.4,26.84,0.2314,536.9,18.44,0.003125,528.1,3.506,0.5785,516.7,36.61);
lipopigm_em_left(:,2) = (a1*exp(-((lipopigm_em_left(:,1)-b1)/c1).^2) + a2*exp(-((lipopigm_em_left(:,1)-b2)/c2).^2) + a3*exp(-((lipopigm_em_left(:,1)-b3)/c3).^2) + a4*exp(-((lipopigm_em_left(:,1)-b4)/c4).^2))';

%cftool(fluop.lipopigm_em(69:end,1),fluop.lipopigm_em(69:end,2))
lipopigm_em_right = (ceil(fluop.lipopigm_em(end,1)):0.5:700)';
[a1,b1,c1,a2,b2,c2] = deal(0.02375,555.4,8.216,0.9859,564.2,64.1);
lipopigm_em_right(:,2) = (a1*exp(-((lipopigm_em_right(:,1)-b1)/c1).^2) + a2*exp(-((lipopigm_em_right(:,1)-b2)/c2).^2))';

fluop.lipopigm_em = cat(1,lipopigm_em_left,fluop.lipopigm_em,lipopigm_em_right);
fluop.lipopigm_em = (abs(fluop.lipopigm_em) + fluop.lipopigm_em)/2;

%cftool(fluop.nadh_em(1:34,1),fluop.nadh_em(1:34,2))
nadh_em_left = (300:0.5:floor(fluop.nadh_em(1,1)))';
[a1,b1,c1,a2,b2,c2] = deal(0.2676,481.5,26.71,0.253,454.5,35.16);
nadh_em_left(:,2) = (a1*exp(-((nadh_em_left(:,1)-b1)/c1).^2) + a2*exp(-((nadh_em_left(:,1)-b2)/c2).^2))';

%cftool(fluop.nadh_em(34:end,1),fluop.nadh_em(34:end,2))
nadh_em_right = (ceil(fluop.nadh_em(end,1)):0.5:700)';
[a1,b1,c1] = deal(0.4493,466,42.5);
nadh_em_right(:,2) = (a1*exp(-((nadh_em_right(:,1)-b1)/c1).^2))';

fluop.nadh_em = cat(1,nadh_em_left,fluop.nadh_em,nadh_em_right);
fluop.nadh_em = (abs(fluop.nadh_em) + fluop.nadh_em)/2;

%cftool(fluop.pirodox_em(1:58,1),fluop.pirodox_em(1:58,2))
pirodox_em_left = (300:0.5:floor(fluop.pirodox_em(1,1)))';
[a1,b1,c1,a2,b2,c2,a3,b3,c3] = deal(0.6,404.8,17.11,0.01391,364,5.93,0.6969,385.2,21.29);
pirodox_em_left(:,2) = (a1*exp(-((pirodox_em_left(:,1)-b1)/c1).^2) + a2*exp(-((pirodox_em_left(:,1)-b2)/c2).^2) + a3*exp(-((pirodox_em_left(:,1)-b3)/c3).^2))';

%cftool(fluop.pirodox_em(58:end,1),fluop.pirodox_em(58:end,2))
pirodox_em_right = (ceil(fluop.pirodox_em(end,1)):0.5:700)';
[a1,b1,c1,a2,b2,c2] = deal(0.7642,394.3,34.13,0.273,414.9,51.6);
pirodox_em_right(:,2) = (a1*exp(-((pirodox_em_right(:,1)-b1)/c1).^2) + a2*exp(-((pirodox_em_right(:,1)-b2)/c2).^2))';

fluop.pirodox_em = cat(1,pirodox_em_left,fluop.pirodox_em,pirodox_em_right);
fluop.pirodox_em = (abs(fluop.pirodox_em) + fluop.pirodox_em)/2;

%cftool(fluop.porphyr_em(1:41,1),fluop.porphyr_em(1:41,2))
porphyr_em_left = (300:0.5:floor(fluop.porphyr_em(1,1)))';
[a1,b1,c1,a2,b2,c2,a3,b3,c3,a4,b4,c4] = deal(1.126,638.4,18.09,0.05018,621.9,1.307,0.06446,617.9,2.08,0.2398,624.1,4.08);
porphyr_em_left(:,2) = (a1*exp(-((porphyr_em_left(:,1)-b1)/c1).^2) + a2*exp(-((porphyr_em_left(:,1)-b2)/c2).^2) + a3*exp(-((porphyr_em_left(:,1)-b3)/c3).^2) + a4*exp(-((porphyr_em_left(:,1)-b4)/c4).^2))';

%cftool(fluop.porphyr_em(87:end,1),fluop.porphyr_em(87:end,2))
fluop.porphyr_em(end,:) = [];  % remove last row outside of 700 nm range
porphyr_em_right = (ceil(fluop.porphyr_em(end,1)):0.5:700)';
[a1,b1,c1,a2,b2,c2,a3,b3,c3,a4,b4,c4] = deal(0,704.8,0.1026,0.3361,684.8,6.821,0.03017,689.9,1.706,0.2308,699.4,22.05);
porphyr_em_right(:,2) = (a1*exp(-((porphyr_em_right(:,1)-b1)/c1).^2) + a2*exp(-((porphyr_em_right(:,1)-b2)/c2).^2) + a3*exp(-((porphyr_em_right(:,1)-b3)/c3).^2) + a4*exp(-((porphyr_em_right(:,1)-b4)/c4).^2))';

fluop.porphyr_em = cat(1,porphyr_em_left,fluop.porphyr_em,porphyr_em_right);
fluop.porphyr_em = (abs(fluop.porphyr_em) + fluop.porphyr_em)/2;

%cftool(fluop.tryptop_em(1:45,1),fluop.tryptop_em(1:45,2))
tryptop_em_left = (300:0.5:floor(fluop.tryptop_em(1,1)))';
[a1,b1,c1,a2,b2,c2] = deal(0.8642,342.4,29.33,0.2056,313.2,17.74);
tryptop_em_left(:,2) = (a1*exp(-((tryptop_em_left(:,1)-b1)/c1).^2) + a2*exp(-((tryptop_em_left(:,1)-b2)/c2).^2))';

%cftool(fluop.tryptop_em(45:end,1),fluop.tryptop_em(45:end,2))
tryptop_em_right = (ceil(fluop.tryptop_em(end,1)):0.5:700)';
[a1,b1,c1,a2,b2,c2] = deal(0.2213,342.3,20.52,0.6675,338.3,52.23);
tryptop_em_right(:,2) = (a1*exp(-((tryptop_em_right(:,1)-b1)/c1).^2) + a2*exp(-((tryptop_em_right(:,1)-b2)/c2).^2))';

fluop.tryptop_em = cat(1,tryptop_em_left,fluop.tryptop_em,tryptop_em_right);
fluop.tryptop_em = (abs(fluop.tryptop_em) + fluop.tryptop_em)/2;

save('fluorophores.mat','fluop')
clearvars -except fieldnames_ex fieldnames_em

load('fluorophores.mat')

%% Dealing with dips: Savitzky-Golay filtering and removal of rows

fluop.elastin_ex(206,:) = [];
fluop.lipopigm_ex(194,:) = [];
fluop.pirodox_ex(106,:) = [];
fluop.pirodox_ex(208,:) = [];  % array size decreased by one due to previous line!
fluop.porphyr_ex(204,:) = [];
fluop.porphyr_ex(327,:) = [];  % array size decreased by one due to previous line!
fluop.tryptop_ex(104,:) = [];

fluop.collagen_em(130,:) = [];
fluop.collagen_em(127:131,2) = sgolayfilt(fluop.collagen_em(127:131,2),2,5);
fluop.elastin_em(60,:) = [];
fluop.elastin_em(130,:) = [];  % array size decreased by one due to previous line!
fluop.lipopigm_em(302,:) = [];
fluop.lipopigm_em(436,:) = [];  % array size decreased by one due to previous line!
fluop.nadh_em(210,:) = [];
fluop.pirodox_em(104,:) = [];
fluop.porphyr_em(600,:) = [];
fluop.porphyr_em(702,:) = [];  % array size decreased by one due to previous line!
fluop.tryptop_em(4,:) = [];

%% Interpolate 

lambda_resampled_fluo_ex = 200:0.5:500;
lambda_resampled_fluo_em = 300:0.5:700;
for i = 1:length(fieldnames_ex)
    fluop.(fieldnames_ex{i}) = interp1(fluop.(fieldnames_ex{i})(:,1),fluop.(fieldnames_ex{i})(:,2),lambda_resampled_fluo_ex)';
    fluop.(fieldnames_em{i}) = interp1(fluop.(fieldnames_em{i})(:,1),fluop.(fieldnames_em{i})(:,2),lambda_resampled_fluo_em)';
end
fluop.wl_ex = lambda_resampled_fluo_ex';
fluop.wl_em = lambda_resampled_fluo_em';

save('fluorophores.mat','fluop')

%% Plot spectra 

figure
plot(lambda_resampled_fluo_ex,fluop.collagen_ex,'linewidth',5)		
hold all
plot(lambda_resampled_fluo_ex,fluop.elastin_ex,'linewidth',5)	
plot(lambda_resampled_fluo_ex,fluop.flavin_ex,'linewidth',5)	
plot(lambda_resampled_fluo_ex,fluop.lipopigm_ex,'linewidth',5)	
plot(lambda_resampled_fluo_ex,fluop.nadh_ex,'linewidth',5)	
plot(lambda_resampled_fluo_ex,fluop.pirodox_ex,'linewidth',5)	
plot(lambda_resampled_fluo_ex,fluop.porphyr_ex,'linewidth',5)	
plot(lambda_resampled_fluo_ex,fluop.tryptop_ex,'linewidth',5)	
legend('Collagen','Elastin','Flavin','Lipopigments','NADH','Piridoxine','Porphyrins','Tryptophan','FontSize',14)
xlabel('\lambda [nm]','FontSize',17)
ylabel('Fluorescence "Absorption" [arb. units]','FontSize',17)
title('Fluorophore Excitation Spectra','FontSize',20);
ylim([0 1])
set(gca,'FontSize',20)

figure
plot(lambda_resampled_fluo_em,fluop.collagen_em,'linewidth',5)		
hold all
plot(lambda_resampled_fluo_em,fluop.elastin_em,'linewidth',5)	
plot(lambda_resampled_fluo_em,fluop.flavin_em,'linewidth',5)	
plot(lambda_resampled_fluo_em,fluop.lipopigm_em,'linewidth',5)	
plot(lambda_resampled_fluo_em,fluop.nadh_em,'linewidth',5)	
plot(lambda_resampled_fluo_em,fluop.pirodox_em,'linewidth',5)	
plot(lambda_resampled_fluo_em,fluop.porphyr_em,'linewidth',5)	
plot(lambda_resampled_fluo_em,fluop.tryptop_em,'linewidth',5)	
legend('Collagen','Elastin','Flavin','Lipopigments','NADH','Piridoxine','Porphyrins','Tryptophan','FontSize',14)
xlabel('\lambda [nm]','FontSize',17)
ylabel('Fluorescence "Emission" [arb. units]','FontSize',17)
title('Fluorophore Emission Spectra','FontSize',20);
ylim([0 1])
set(gca,'FontSize',20)

clear