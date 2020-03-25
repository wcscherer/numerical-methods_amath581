xspan = [-1 1];
npoints = 8;

[A1 span] = Dx_1D(xspan,npoints,0);

A1 = full(A1);

[A2 span] = Dx_1D(xspan,npoints,1);

A2 = full(A2);
