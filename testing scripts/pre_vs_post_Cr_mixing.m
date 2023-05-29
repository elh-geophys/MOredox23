% Thinking about Cr oxidation
% EL
% April 2, 2023

Mm = 1;

melt = 0.3*Mm;
solid = 0.7*Mm;

Fe_melt_pre = 0.1200;
Fe_solid_post = 0.03;

Fe_melt_post = Fe_melt_pre - 0.35/8.05;

mix_postCr = (Fe_melt_post*melt + Fe_solid_post*solid) / Mm;

Fe_solid_pre = Fe_solid_post + 0.35/8.05;

mix_preCr = (Fe_melt_pre*melt + Fe_solid_pre*solid) / Mm - 0.35/8.05;
