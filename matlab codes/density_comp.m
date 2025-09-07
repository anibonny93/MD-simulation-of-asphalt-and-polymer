close all
clear all
atomistic_data = readlammpslog('/mnt/Shared_Data/anic/LONI_output/New_files/Kspace_pppm_long/2025/Density_measure/AP_oxi_d30/512/AP_oxi_d30_atomistic.log');
atomistic_mass = ReadMassFromLAMMPSData('/mnt/Shared_Data/anic/LONI_output/New_files/Kspace_pppm_long/2025/Density_measure/AP_oxi_d30/512/AP_oxi_d30_atomistic_replicated.data');
atomistic_mass = atomistic_mass*256;
atomistic_density = atomistic_mass*1.661./atomistic_data.Volume;


opts=delimitedTextImportOptions('Delimiter',' ','ConsecutiveDelimitersRule','join');
opts.VariableTypes='double'

cg_data = readmatrix('/mnt/Shared_Data/anic/LONI_output/New_files/Kspace_pppm_long/2025/Density_measure/AAA1_all_oxidized/variable_temp/AAA1_cg_all_oxidized.log',opts);
vol=cg_data(1:end,3);
vol=vol(~isnan(vol))
vol=vol(vol>1e4);
vol=vol(vol<3e6);

%cg_density = atomistic_mass*1.661./vol;
cg_density =263052.153*1.661./vol;
avg_atomistic_density = sum(atomistic_density)/length(atomistic_density);
avg_cg_density = sum(cg_density(900:5800))/length(cg_density(900:5800));


figure
plot(1:length(cg_density(900:6388)),cg_density(900:6388));
hold on
plot(1:length(atomistic_density),atomistic_density);
set(gca,'XDir','reverse');
xlim([2000 5500]);
ylim([1.1,1.3]);

xlabel('1000 femtoseconds');
ylabel('Density (g/cm^3)');
legend(sprintf('Average coarse-grained density %.4f',avg_cg_density), sprintf('Average atomistic density %.4f',avg_atomistic_density),'Location','best');
title('Asphaltene phenol density comparison');