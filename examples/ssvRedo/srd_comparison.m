gridsize = [26, 54, 108];




lts_fv1_v       = [5.5853279394542333E-003, 1.1149525213654070E-003, 2.3549655340437793E-004];
srd1_fv1_v      = [6.0794582719505622E-003, 1.3042574491654634E-003, 2.6755711249439187E-004];
lts_fv1q_v      = [2.4844256695815614E-003, 5.8801730618267332E-004, 1.4221501276227199E-004]; 
srd1q_fv1q_v    = [2.7913230728114189E-003, 5.6644065885654710E-004, 1.3791945588817643E-004]; %, 3.2539184498343274E-005

lts_fv1_b       = [7.8707896557172308E-002, 2.6110416966706587E-002, 9.3126781246204837E-003];
srd1_fv1_b      = [8.2552342711314800E-002, 2.8236528741063492E-002, 1.0141414806058742E-002];
lts_fv1q_b      = [3.8628512141920035E-002, 1.2501956972788070E-002, 4.8114614519942081E-003]; 
srd1q_fv1q_b    = [3.8809709936814089E-002, 1.1594032084576614E-002, 4.4020764964905748E-003]; 

lw = 3;
ms = 25;
loglog(gridsize, srd1_fv1_b, '.--b', 'markersize', ms, 'linewidth', lw); hold on;
loglog(gridsize, srd1q_fv1q_b, '.--r', 'markersize', ms, 'linewidth', lw); hold on;
loglog(gridsize, lts_fv1_b, '.--c', 'markersize', ms, 'linewidth', lw); hold on;
loglog(gridsize, lts_fv1q_b, '.--k', 'markersize', ms, 'linewidth', lw); hold on;

loglog(gridsize, srd1_fv1_v, '.-b', 'markersize', ms, 'linewidth', lw); hold on;
loglog(gridsize, srd1q_fv1q_v, '.-r', 'markersize', ms, 'linewidth', lw); hold on;
loglog(gridsize, lts_fv1_v, '.-c', 'markersize', ms, 'linewidth', lw); hold on;
loglog(gridsize, lts_fv1q_v, '.-k', 'markersize', ms, 'linewidth', lw); hold on;

grid on;
set(gca, "linewidth", 2);
xlabel("Grid size (N)", "fontsize",12);
ylabel("L1 error", "fontsize",12);
legend("first order gradients + SRD (L1 bdry)", "second order gradients + SRD (L1 bdry)", "local time stepping with first order gradients (L1 bdry)", "local time stepping with second order gradients (L1 bdry)", ...
       "first order gradients + SRD (L1 volume)", "second order gradients + SRD (L1 volume)", "local time stepping with first order gradients (L1 volume)", "local time stepping with second order gradients (L1 volume)","fontsize", 12);