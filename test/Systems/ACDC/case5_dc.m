% based on NESTA v0.6.0
% tests dc line with costs
% tests generator and dc line voltage setpoint warnings

function mpc = nesta_case5_pjm
mpc.version = '2';
mpc.baseMVA = 100.0;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	 1	   0.0	   0.00	 0.0	 0.0	 1	    1.06355	    2.87619	 230.0	 1	    1.10000	    0.90000;
	2	 1	 300.0	  98.61	 0.0	 0.0	 1	    1.08009	   -0.79708	 230.0	 1	    1.10000	    0.90000;
	3	 1	 300.0	  98.61	 0.0	 0.0	 1	    1.10000	   -0.64925	 230.0	 1	    1.10000	    0.90000;
	4	 3	 400.0	 131.47	 0.0	 0.0	 1	    1.06414	    0.00000	 230.0	 1	    1.10000	    0.90000;
	5	 1	   0.0	   0.00	 0.0	 0.0	 1	    1.05304	    3.70099	 230.0	 1	    1.10000	    0.90000;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
mpc.gen = [
	1	 40.0	 30.0	 30.0	 -30.0	 1.06355	 100.0	 1	 40.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0;
	1	 170.0	 127.5	 127.5	 -127.5	 1.07762	 100.0	 1	 170.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0;
	3	 333.6866	 390.0	 390.0	 -390.0	 1.1	 100.0	 1	 520.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0;
	4	 0.0	 -70.8186	 150.0	 -150.0	 1.06414	 100.0	 1	 200.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0;
	5	 463.5555	 -184.9224	 450.0	 -450.0	 1.06907	 100.0	 1	 600.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0;
];

%% generator cost data
%	2	startup	shutdown	n	c(n-1)	...	c0
mpc.gencost = [
	2	 0.0	 0.0	 3	   0.000000	  14.000000	   0.000000;
	2	 0.0	 0.0	 3	   0.000000	  15.000000	   0.000000;
	2	 0.0	 0.0	 3	   0.000000	  30.000000	   0.000000;
	2	 0.0	 0.0	 3	   0.000000	  40.000000	   0.000000;
	2	 0.0	 0.0	 3	   0.000000	  10.000000	   0.000000;
];

%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [
	1	 2	 0.00281	 0.0281	 0.00712	 400.0	 400.0	 400.0	 0.0	 0.0	 1	 -30.0	 30.0;
	1	 4	 0.00304	 0.0304	 0.00658	 426	 426	 426	 0.0	 0.0	 1	 -30.0	 30.0;
	1	 5	 0.00064	 0.0064	 0.03126	 426	 426	 426	 0.0	 0.0	 1	 -30.0	 30.0;
	2	 3	 0.00108	 0.0108	 0.01852	 426	 426	 426	 0.0	 0.0	 1	 -30.0	 30.0;
	3	 4	 0.00297	 0.0297	 0.00674	 426	 426	 426	 0.0	 0.0	 1	 -30.0	 30.0;
	4	 5	 0.00297	 0.0297	 0.00674	 240.0	 240.0	 240.0	 0.0	 0.0	 1	 -30.0	 30.0;
];

%% dcline data
%	fbus	tbus	status	Pf	Pt	Qf	Qt	Vf	Vt	Pmin	Pmax	QminF	QmaxF	QminT	QmaxT	loss0	loss1
mpc.dcline = [
	3	5	1	10	8.9	99.9934	-10.4049	1.1	1.05555	10	100 	-100	100	-100 100	1	0.01;
];

%% dcline cost data
%	2	startup	shutdown	n	c(n-1)	...	c0
mpc.dclinecost = [
	2	 0.0	 0.0	 3	   0.000000	  40.000000	   0.000000;
];

% adds current ratings to branch matrix
%column_names%	c_rating_a
mpc.branch_currents = [
400;426;426;426;426;240;
];
