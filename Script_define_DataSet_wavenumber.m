clear all
close all
%%
load('w_hat.mat')

w_hat = w_hat / max( abs( w_hat(:) ) );
w_hat = fftshift(w_hat);

data.description = cat(2,...
'This record contains a simulated wavefield in wavenumber space of one filtered frequency,\n which was evaluated on the surface of a 2 mm thick steel plate.\n',...
'The simulation contains a rectangular defect at a depth of 1mm,\n the corner points of which are noted in the data set.\n',...
'The simulation has an equidistant grid in both spatial directions.\n',...
'The first step size corresponds to the row index,\n while the second step size corresponds to the column index.\n',...
'In addition, the approximate wave numbers for different plate thicknesses are noted.\n');

data.simulation.spectrum.description = 'complex spectrum of the wavefield in the shifted wavenumber space by filtering one frequency';
data.simulation.spectrum.value       = w_hat;
data.simulation.spectrum.unit        = 'normalized';

data.simulation.frequency.description = 'the frequency, which was filtered out';
data.simulation.frequency.value       = 150e3;
data.simulation.frequency.unit        = 'Hz';

data.simulation.N_POINT.description = 'number of points in the equidistant grid';
data.simulation.N_POINT.value       = [300, 300];
data.simulation.N_POINT.unit        = '';

data.simulation.dx.description = 'step size of the equidistant grid';
data.simulation.dx.value       = [1e-3, 1e-3];
data.simulation.dx.unit        = 'm';

data.simulation.thickness.description = 'thickness of the plate';
data.simulation.thickness.value       = 2e-3;
data.simulation.thickness.unit        = 'm';

data.simulation.material.description = 'steel';

data.defect.position.description = 'corner positions of a defect';
data.defect.position.value(1,:)  = [0.08, 0.12, 0.12, 0.08];
data.defect.position.value(2,:)  = [0.13, 0.13, 0.17, 0.17];
data.defect.position.value(3,:)  = [1e-3, 1e-3, 1e-3, 1e-3];
data.defect.position.unit        = 'm';

data.relation.description      = 'relationship between wavenumbers and plate thickness for a given frequency';
data.relation.frequency.value  = 150e3;
data.relation.frequency.unit   = 'Hz';
data.relation.wavenumber.value = [600, 650, 700, 750, 800];
data.relation.wavenumber.unit  = '1/m';
data.relation.thickness.value  = [2, 1.75, 1.5, 1.25, 1];
data.relation.thickness.unit   = 'mm';

save('Data.mat','-struct','data')
