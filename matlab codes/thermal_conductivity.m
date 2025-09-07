tke = 62155.65/2;
tke_j = (tke*4184)/6.022e23;
box_dim = 62.51e-10;
cross_sec = box_dim^2;
time_fs = 2e6;
time = time_fs*1e-15;
heat_flux = tke_j/(time*cross_sec);
delta_T = 62.022;
delta_Z = box_dim/2;

therm_cond = heat_flux/(delta_T/delta_Z)