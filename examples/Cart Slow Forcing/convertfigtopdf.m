%Processing figures in Matlab
fig = openfig('Parameter_Conservative_Plots.fig');
origUnits = fig.Units;
fig.Units = fig.PaperUnits;
fig.PaperSize = fig.Position(3:4);
fig.Units = origUnits;
exportgraphics(fig, 'Parameter_Conservative_Plots.pdf');

