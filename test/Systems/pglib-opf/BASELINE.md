# Baseline Results
This section presents the results of running 
[PowerModels.jl](https://github.com/lanl-ansi/PowerModels.jl) 
on the OPF benchmarks from the
[PGLib Archive](https://github.com/power-grid-lib/pglib-opf). 
providing baseline results for these test cases. All models were solved using 
[IPOPT](https://link.springer.com/article/10.1007/s10107-004-0559-y).

Note that the displayed solve times are only approximate.
These times do not include Julia's JIT time, around 2-5 seconds, and
use the [HSL](http://www.hsl.rl.ac.uk/ipopt/) ma27 solver in IPOPT.
The default linear solver will increase the runtime by 2-6x.

## Software Versions
**[PowerModels.jl](https://github.com/lanl-ansi/PowerModels.jl):** v0.19.1

**[Ipopt.jl](https://github.com/JuliaOpt/Ipopt.jl):** v0.9.1

**[PGLib OPF](https://github.com/power-grid-lib/pglib-opf):** v21.07

**Hardware:** Dual Intel 2.10GHz CPUs, 128GB RAM


## Typical Operating Conditions (TYP)
| **Case Name** | **Nodes** | **Edges** | **DC (\$/h)** | **AC (\$/h)** | **QC Gap (%)** | **SOC Gap (%)** | **DC Time (sec.)** | **AC Time (sec.)** | **QC Time (sec.)** | **SOC Time (sec.)** |
| ------------- | --------- | --------- | ------------- | ------------- | -------------- | --------------- | ------------------ | ------------------ | ------------------ | ------------------- |
| pglib_opf_case3_lmbd | 3 | 3 | 5.6959e+03 | 5.8126e+03 | 1.22 | 1.32 | <1 | <1 | <1 | <1 |
| pglib_opf_case5_pjm | 5 | 6 | 1.7480e+04 | 1.7552e+04 | 14.55 | 14.55 | <1 | <1 | <1 | <1 |
| pglib_opf_case14_ieee | 14 | 20 | 2.0515e+03 | 2.1781e+03 | 0.11 | 0.11 | <1 | <1 | <1 | <1 |
| pglib_opf_case24_ieee_rts | 24 | 38 | 6.1001e+04 | 6.3352e+04 | 0.02 | 0.02 | <1 | <1 | <1 | <1 |
| pglib_opf_case30_as | 30 | 41 | 7.6760e+02 | 8.0313e+02 | 0.06 | 0.06 | <1 | <1 | <1 | <1 |
| pglib_opf_case30_ieee | 30 | 41 | 7.4728e+03 | 8.2085e+03 | 18.81 | 18.84 | <1 | <1 | <1 | <1 |
| pglib_opf_case39_epri | 39 | 46 | 1.3689e+05 | 1.3842e+05 | 0.55 | 0.56 | <1 | <1 | <1 | <1 |
| pglib_opf_case57_ieee | 57 | 80 | 3.4773e+04 | 3.7589e+04 | 0.16 | 0.16 | <1 | <1 | <1 | <1 |
| pglib_opf_case60_c | 60 | 88 | 9.0700e+04 | 9.2694e+04 | 0.06 | 0.07 | <1 | <1 | <1 | <1 |
| pglib_opf_case73_ieee_rts | 73 | 120 | 1.8300e+05 | 1.8976e+05 | 0.04 | 0.04 | <1 | <1 | <1 | <1 |
| pglib_opf_case89_pegase | 89 | 210 | 1.0504e+05 | 1.0729e+05 | 0.75 | 0.75 | <1 | <1 | <1 | <1 |
| pglib_opf_case118_ieee | 118 | 186 | 9.3101e+04 | 9.7214e+04 | 0.79 | 0.91 | <1 | <1 | <1 | <1 |
| pglib_opf_case162_ieee_dtc | 162 | 284 | 1.0146e+05 | 1.0808e+05 | 5.84 | 5.95 | <1 | <1 | <1 | <1 |
| pglib_opf_case179_goc | 179 | 263 | 7.5188e+05 | 7.5427e+05 | 0.16 | 0.16 | <1 | <1 | <1 | <1 |
| pglib_opf_case200_activ | 200 | 245 | 2.7480e+04 | 2.7558e+04 | 0.01 | 0.01 | <1 | <1 | <1 | <1 |
| pglib_opf_case240_pserc | 240 | 448 | 3.2714e+06 | 3.3297e+06 | 2.73 | 2.78 | <1 | 4 | 6 | 2 |
| pglib_opf_case300_ieee | 300 | 411 | 5.1785e+05 | 5.6522e+05 | 2.58 | 2.63 | <1 | <1 | 2 | <1 |
| pglib_opf_case500_goc | 500 | 733 | 4.4055e+05 | 4.5495e+05 | 0.25 | 0.25 | <1 | <1 | 2 | <1 |
| pglib_opf_case588_sdet | 588 | 686 | 3.1013e+05 | 3.1314e+05 | 1.91 | 2.14 | <1 | <1 | 2 | <1 |
| pglib_opf_case793_goc | 793 | 913 | 2.5831e+05 | 2.6020e+05 | 1.32 | 1.33 | <1 | 2 | 3 | <1 |
| pglib_opf_case1354_pegase | 1354 | 1991 | 1.2182e+06 | 1.2588e+06 | 1.56 | 1.57 | <1 | 5 | 8 | 4 |
| pglib_opf_case1888_rte | 1888 | 2531 | 1.3529e+06 | 1.4025e+06 | 2.05 | 2.05 | <1 | 8 | 10 | 48 |
| pglib_opf_case1951_rte | 1951 | 2596 | 2.0316e+06 | 2.0856e+06 | 0.13 | 0.14 | <1 | 17 | 12 | 8 |
| pglib_opf_case2000_goc | 2000 | 3639 | 9.4304e+05 | 9.7343e+05 | 0.31 | 0.31 | <1 | 6 | 12 | 6 |
| pglib_opf_case2312_goc | 2312 | 3013 | 4.4033e+05 | 4.4133e+05 | 1.90 | 1.90 | <1 | 6 | 13 | 5 |
| pglib_opf_case2383wp_k | 2383 | 2896 | 1.8041e+06 | 1.8682e+06 | 0.97 | 1.04 | <1 | 7 | 13 | 7 |
| pglib_opf_case2736sp_k | 2736 | 3504 | 1.2760e+06 | 1.3080e+06 | 0.30 | 0.31 | <1 | 5 | 14 | 6 |
| pglib_opf_case2737sop_k | 2737 | 3506 | 7.6401e+05 | 7.7773e+05 | 0.26 | 0.27 | <1 | 4 | 11 | 5 |
| pglib_opf_case2742_goc | 2742 | 4673 | 2.5970e+05 | 2.7571e+05 | 1.33 | 1.35 | <1 | 28 | 17 | 7 |
| pglib_opf_case2746wop_k | 2746 | 3514 | 1.1782e+06 | 1.2083e+06 | 0.36 | 0.37 | <1 | 5 | 12 | 5 |
| pglib_opf_case2746wp_k | 2746 | 3514 | 1.5814e+06 | 1.6317e+06 | 0.32 | 0.33 | <1 | 6 | 12 | 6 |
| pglib_opf_case2848_rte | 2848 | 3776 | 1.2677e+06 | 1.2866e+06 | 0.12 | 0.13 | <1 | 16 | 17 | 10 |
| pglib_opf_case2853_sdet | 2853 | 3921 | 2.0370e+06 | 2.0524e+06 | 0.87 | 0.91 | <1 | 8 | 16 | 8 |
| pglib_opf_case2868_rte | 2868 | 3808 | 1.9667e+06 | 2.0096e+06 | 0.10 | 0.10 | <1 | 14 | 18 | 11 |
| pglib_opf_case2869_pegase | 2869 | 4582 | 2.3864e+06 | 2.4628e+06 | 1.01 | 1.01 | <1 | 12 | 25 | 10 |
| pglib_opf_case3012wp_k | 3012 | 3572 | 2.5090e+06 | 2.6008e+06 | 0.98 | 1.03 | <1 | 9 | 18 | 15 |
| pglib_opf_case3022_goc | 3022 | 4135 | 5.9922e+05 | 6.0138e+05 | 2.76 | 2.77 | <1 | 9 | 17 | 6 |
| pglib_opf_case3120sp_k | 3120 | 3693 | 2.0880e+06 | 2.1480e+06 | 0.55 | 0.56 | <1 | 9 | 16 | 7 |
| pglib_opf_case3375wp_k | 3374 | 4161 | 7.3170e+06 | 7.4382e+06 | 0.54 | 0.55 | <1 | 10 | 61 | 9 |
| pglib_opf_case3970_goc | 3970 | 6641 | 9.3422e+05 | 9.6099e+05 | 0.34 | 0.34 | 2 | 19 | 28 | 16 |
| pglib_opf_case4020_goc | 4020 | 6988 | 7.9506e+05 | 8.2225e+05 | 1.22 | 1.23 | 2 | 20 | 43 | 17 |
| pglib_opf_case4601_goc | 4601 | 7199 | 7.9381e+05 | 8.2624e+05 | 0.53 | 0.54 | 2 | 22 | 36 | 18 |
| pglib_opf_case4619_goc | 4619 | 8150 | 4.5744e+05 | 4.7670e+05 | 0.88 | 0.91 | 2 | 19 | 41 | 30 |
| pglib_opf_case4661_sdet | 4661 | 5997 | 2.2163e+06 | 2.2513e+06 | 1.89 | 1.99 | <1 | 15 | 34 | 15 |
| pglib_opf_case4837_goc | 4837 | 7765 | 8.5040e+05 | 8.7226e+05 | 0.47 | 0.47 | 2 | 19 | 36 | 14 |
| pglib_opf_case4917_goc | 4917 | 6726 | 1.3837e+06 | 1.3878e+06 | 2.50 | 2.50 | <1 | 17 | 35 | 12 |
| pglib_opf_case6468_rte | 6468 | 9000 | 1.9828e+06 | 2.0697e+06 | 1.12 | 1.13 | 2 | 61 | 106 | 40 |
| pglib_opf_case6470_rte | 6470 | 9005 | 2.1361e+06 | 2.2376e+06 | 1.75 | 1.76 | 2 | 34 | 57 | 33 |
| pglib_opf_case6495_rte | 6495 | 9019 | 2.5618e+06 | 3.0678e+06 | 15.09 | 15.11 | 2 | 66 | 80 | 34 |
| pglib_opf_case6515_rte | 6515 | 9037 | 2.5593e+06 | 2.8255e+06 | 6.39 | 6.40 | 2 | 51 | 57 | 31 |
| pglib_opf_case8387_pegase | 8387 | 14561 | 2.5028e+06 | 2.7714e+06 | 54.61 | 64.11 | 3 | 45 | 140 | 39 |
| pglib_opf_case9241_pegase | 9241 | 16049 | 6.0287e+06 | 6.2431e+06 | 1.71 | 2.54 | 3 | 46 | 135 | 49 |
| pglib_opf_case9591_goc | 9591 | 15915 | 1.0309e+06 | 1.0617e+06 | 0.61 | 0.62 | 3 | 58 | 147 | 58 |
| pglib_opf_case10000_goc | 10000 | 13193 | 1.3461e+06 | 1.3540e+06 | 1.51 | 1.54 | 3 | 46 | 125 | 38 |
| pglib_opf_case10480_goc | 10480 | 18559 | 2.2158e+06 | 2.3146e+06 | 1.15 | 1.23 | 4 | 63 | 143 | 55 |
| pglib_opf_case13659_pegase | 13659 | 20467 | 8.7699e+06 | 8.9480e+06 | 0.98 | 1.39 | 3 | 59 | 157 | 81 |
| pglib_opf_case19402_goc | 19402 | 34704 | 1.8978e+06 | 1.9778e+06 | 1.18 | 1.19 | 9 | 132 | 400 | 144 |
| pglib_opf_case24464_goc | 24464 | 37816 | 2.5128e+06 | 2.6295e+06 | 0.95 | 0.97 | 7 | 117 | 401 | 129 |
| pglib_opf_case30000_goc | 30000 | 35393 | 1.0921e+06 | 1.1423e+06 | 2.88 | 2.89 | 12 | 265 | 615 | 169 |


## Congested Operating Conditions (API)
| **Case Name** | **Nodes** | **Edges** | **DC (\$/h)** | **AC (\$/h)** | **QC Gap (%)** | **SOC Gap (%)** | **DC Time (sec.)** | **AC Time (sec.)** | **QC Time (sec.)** | **SOC Time (sec.)** |
| ------------- | --------- | --------- | ------------- | ------------- | -------------- | --------------- | ------------------ | ------------------ | ------------------ | ------------------- |
| pglib_opf_case3_lmbd__api | 3 | 3 | 1.0444e+04 | 1.1236e+04 | 5.60 | 9.27 | <1 | <1 | <1 | <1 |
| pglib_opf_case5_pjm__api | 5 | 6 | 7.5433e+04 | 7.6377e+04 | 4.09 | 4.09 | <1 | <1 | <1 | <1 |
| pglib_opf_case14_ieee__api | 14 | 20 | 4.7976e+03 | 5.9994e+03 | 5.13 | 5.13 | <1 | <1 | <1 | <1 |
| pglib_opf_case24_ieee_rts__api | 24 | 38 | 1.2336e+05 | 1.3494e+05 | 13.01 | 17.88 | <1 | <1 | <1 | <1 |
| pglib_opf_case30_as__api | 30 | 41 | 3.0921e+03 | 4.9962e+03 | 44.61 | 44.61 | <1 | <1 | <1 | <1 |
| pglib_opf_case30_ieee__api | 30 | 41 | 1.6147e+04 | 1.8044e+04 | 5.46 | 5.46 | <1 | <1 | <1 | <1 |
| pglib_opf_case39_epri__api | 39 | 46 | 2.4462e+05 | 2.4967e+05 | 1.70 | 1.72 | <1 | <1 | <1 | <1 |
| pglib_opf_case57_ieee__api | 57 | 80 | 4.7421e+04 | 4.9290e+04 | 0.08 | 0.08 | <1 | <1 | <1 | <1 |
| pglib_opf_case60_c__api | 60 | 88 | 1.7618e+05 | 1.8524e+05 | 2.28 | 2.30 | <1 | <1 | <1 | <1 |
| pglib_opf_case73_ieee_rts__api | 73 | 120 | 3.6236e+05 | 4.2263e+05 | 11.05 | 12.87 | <1 | <1 | <1 | <1 |
| pglib_opf_case89_pegase__api | 89 | 210 | 1.1729e+05 | 1.3017e+05 | 23.07 | 23.11 | <1 | <1 | <1 | <1 |
| pglib_opf_case118_ieee__api | 118 | 186 | 2.2448e+05 | 2.4224e+05 | 29.74 | 29.97 | <1 | <1 | <1 | <1 |
| pglib_opf_case162_ieee_dtc__api | 162 | 284 | 1.1166e+05 | 1.2099e+05 | 4.33 | 4.36 | <1 | <1 | <1 | <1 |
| pglib_opf_case179_goc__api | 179 | 263 | 1.8900e+06 | 1.9320e+06 | 5.93 | 9.88 | <1 | <1 | 2 | <1 |
| pglib_opf_case200_activ__api | 200 | 245 | 3.5008e+04 | 3.5701e+04 | 0.03 | 0.03 | <1 | <1 | <1 | <1 |
| pglib_opf_case240_pserc__api | 240 | 448 | 4.5760e+06 | 4.6406e+06 | 0.64 | 0.67 | <1 | 4 | 6 | 3 |
| pglib_opf_case300_ieee__api | 300 | 411 | 6.6114e+05 | 6.8499e+05 | 0.83 | 0.85 | <1 | <1 | 2 | <1 |
| pglib_opf_case500_goc__api | 500 | 733 | 6.5084e+05 | 6.9241e+05 | 3.44 | 3.44 | <1 | 2 | 2 | <1 |
| pglib_opf_case588_sdet__api | 588 | 686 | 3.8834e+05 | 3.9476e+05 | 1.39 | 1.61 | <1 | <1 | 3 | 2 |
| pglib_opf_case793_goc__api | 793 | 913 | 3.1086e+05 | 3.1885e+05 | 13.21 | 13.43 | <1 | 2 | 4 | <1 |
| pglib_opf_case1354_pegase__api | 1354 | 1991 | 1.4558e+06 | 1.4983e+06 | 0.54 | 0.55 | <1 | 4 | 9 | 5 |
| pglib_opf_case1888_rte__api | 1888 | 2531 | 1.9019e+06 | 1.9539e+06 | 0.22 | 0.23 | <1 | 7 | 19 | 9 |
| pglib_opf_case1951_rte__api | 1951 | 2596 | 2.3361e+06 | 2.4108e+06 | 0.51 | 0.57 | <1 | 8 | 11 | 7 |
| pglib_opf_case2000_goc__api | 2000 | 3639 | 1.3934e+06 | 1.4686e+06 | 2.07 | 2.07 | <1 | 7 | 10 | 5 |
| pglib_opf_case2312_goc__api | 2312 | 3013 | 5.3368e+05 | 5.7152e+05 | 13.08 | 13.09 | <1 | 10 | 11 | 5 |
| pglib_opf_case2383wp_k__api | 2383 | 2896 | 2.7913e+05 | 2.7913e+05 | 0.01 | 0.01 | <1 | 2 | 4 | 2 |
| pglib_opf_case2736sp_k__api | 2736 | 3504 | 6.1114e+05 | 6.5394e+05 | 10.83 | 10.84 | <1 | 7 | 12 | 5 |
| pglib_opf_case2737sop_k__api | 2737 | 3506 | 3.4557e+05 | 3.6717e+05 | 5.89 | 5.89 | <1 | 6 | 5 | 3 |
| pglib_opf_case2742_goc__api | 2742 | 4673 | 5.2504e+05 | 6.4219e+05 | 24.40 | 24.45 | <1 | 22 | 20 | 7 |
| pglib_opf_case2746wop_k__api | 2746 | 3514 | 5.1166e+05 | 5.1166e+05 | 0.01 | 0.01 | <1 | 2 | 4 | 2 |
| pglib_opf_case2746wp_k__api | 2746 | 3514 | 5.8183e+05 | 5.8183e+05 | 0.01 | 0.00 | <1 | 3 | 6 | 3 |
| pglib_opf_case2848_rte__api | 2848 | 3776 | 1.4624e+06 | 1.4970e+06 | 0.25 | 0.25 | <1 | 26 | 23 | 9 |
| pglib_opf_case2853_sdet__api | 2853 | 3921 | 2.4300e+06 | 2.4578e+06 | 1.92 | 1.96 | <1 | 10 | 16 | 8 |
| pglib_opf_case2868_rte__api | 2868 | 3808 | 2.2344e+06 | 2.2946e+06 | 0.17 | 0.18 | <1 | 22 | 19 | 9 |
| pglib_opf_case2869_pegase__api | 2869 | 4582 | 2.8484e+06 | 2.9296e+06 | 0.99 | 1.00 | <1 | 12 | 26 | 10 |
| pglib_opf_case3012wp_k__api | 3012 | 3572 | 7.2887e+05 | 7.2887e+05 | 0.00 | 0.00 | <1 | 4 | 7 | 3 |
| pglib_opf_case3022_goc__api | 3022 | 4135 | 6.3985e+05 | 6.5189e+05 | 9.80 | 9.82 | <1 | 9 | 18 | 6 |
| pglib_opf_case3120sp_k__api | 3120 | 3693 | 8.8346e+05 | 9.3692e+05 | 23.90 | 23.96 | 2 | 12 | 19 | 5 |
| pglib_opf_case3375wp_k__api | 3374 | 4161 | 5.7582e+06 | 5.8478e+06 | 9.35 | -- | 2 | 11 | 27 | 681 |
| pglib_opf_case3970_goc__api | 3970 | 6641 | 1.0784e+06 | 1.4557e+06 | 34.53 | 34.58 | 2 | 36 | 28 | 12 |
| pglib_opf_case4020_goc__api | 4020 | 6988 | 1.0853e+06 | 1.2979e+06 | 17.56 | 17.58 | 2 | 37 | 35 | 12 |
| pglib_opf_case4601_goc__api | 4601 | 7199 | 7.1860e+05 | 7.9253e+05 | 20.39 | 20.51 | 2 | 28 | 32 | 12 |
| pglib_opf_case4619_goc__api | 4619 | 8150 | 9.5405e+05 | 1.0299e+06 | 8.24 | 8.30 | 3 | 27 | 39 | 20 |
| pglib_opf_case4661_sdet__api | 4661 | 5997 | 2.6410e+06 | 2.6953e+06 | 2.54 | 2.64 | <1 | 16 | 32 | 123 |
| pglib_opf_case4837_goc__api | 4837 | 7765 | 1.0693e+06 | 1.1578e+06 | 8.10 | 8.15 | 2 | 35 | 44 | 16 |
| pglib_opf_case4917_goc__api | 4917 | 6726 | 1.5417e+06 | 1.5479e+06 | 8.52 | 8.54 | 2 | 17 | 34 | 12 |
| pglib_opf_case6468_rte__api | 6468 | 9000 | 2.2215e+06 | 2.3135e+06 | 0.80 | 0.82 | 2 | 79 | 74 | 302 |
| pglib_opf_case6470_rte__api | 6470 | 9005 | 2.4983e+06 | 2.6065e+06 | 1.19 | 1.20 | 2 | 50 | 49 | 98 |
| pglib_opf_case6495_rte__api | 6495 | 9019 | 2.8138e+06 | 2.9750e+06 | 1.99 | 2.04 | 2 | 51 | 80 | 30 |
| pglib_opf_case6515_rte__api | 6515 | 9037 | 2.8952e+06 | 3.0617e+06 | 1.87 | 1.89 | 2 | 58 | 62 | 29 |
| pglib_opf_case8387_pegase__api | 8387 | 14561 | 5.3453e+06 | 5.5571e+06 | 8.57 | 10.33 | 4 | 45 | 109 | 40 |
| pglib_opf_case9241_pegase__api | 9241 | 16049 | 6.7974e+06 | 7.0112e+06 | 1.84 | 2.67 | 3 | 64 | 240 | 1530 |
| pglib_opf_case9591_goc__api | 9591 | 15915 | 1.3274e+06 | 1.4259e+06 | 13.79 | 13.84 | 6 | 136 | 143 | 35 |
| pglib_opf_case10000_goc__api | 10000 | 13193 | 2.2110e+06 | 2.3728e+06 | 7.60 | 7.65 | 4 | 51 | 96 | 36 |
| pglib_opf_case10480_goc__api | 10480 | 18559 | 2.6217e+06 | 2.7627e+06 | 4.72 | 4.78 | 4 | 76 | 117 | 46 |
| pglib_opf_case13659_pegase__api | 13659 | 20467 | 9.0372e+06 | 9.2842e+06 | 1.19 | 1.85 | 3 | 54 | 169 | 106 |
| pglib_opf_case19402_goc__api | 19402 | 34704 | 2.2861e+06 | 2.3987e+06 | 4.71 | 4.75 | 11 | 187 | 310 | 124 |
| pglib_opf_case24464_goc__api | 24464 | 37816 | 2.3433e+06 | 2.4723e+06 | 3.93 | 4.00 | 9 | 134 | 344 | 187 |
| pglib_opf_case30000_goc__api | 30000 | 35393 | 1.3219e+06 | 1.3530e+06 | 24.74 | 24.79 | 19 | 522 | 741 | 160 |


## Small Angle Difference Conditions (SAD)
| **Case Name** | **Nodes** | **Edges** | **DC (\$/h)** | **AC (\$/h)** | **QC Gap (%)** | **SOC Gap (%)** | **DC Time (sec.)** | **AC Time (sec.)** | **QC Time (sec.)** | **SOC Time (sec.)** |
| ------------- | --------- | --------- | ------------- | ------------- | -------------- | --------------- | ------------------ | ------------------ | ------------------ | ------------------- |
| pglib_opf_case3_lmbd__sad | 3 | 3 | 5.8560e+03 | 5.9593e+03 | 1.42 | 3.75 | <1 | <1 | <1 | <1 |
| pglib_opf_case5_pjm__sad | 5 | 6 | inf. | 2.6109e+04 | 0.99 | 3.62 | <1 | <1 | <1 | <1 |
| pglib_opf_case14_ieee__sad | 14 | 20 | inf. | 2.7768e+03 | 21.48 | 21.53 | <1 | <1 | <1 | <1 |
| pglib_opf_case24_ieee_rts__sad | 24 | 38 | 7.8122e+04 | 7.6918e+04 | 2.93 | 9.55 | <1 | <1 | <1 | <1 |
| pglib_opf_case30_as__sad | 30 | 41 | inf. | 8.9735e+02 | 2.31 | 7.88 | <1 | <1 | <1 | <1 |
| pglib_opf_case30_ieee__sad | 30 | 41 | inf. | 8.2085e+03 | 5.94 | 9.70 | <1 | <1 | <1 | <1 |
| pglib_opf_case39_epri__sad | 39 | 46 | 1.5067e+05 | 1.4834e+05 | 0.21 | 0.67 | <1 | <1 | <1 | <1 |
| pglib_opf_case57_ieee__sad | 57 | 80 | inf. | 3.8663e+04 | 0.35 | 0.71 | <1 | <1 | <1 | <1 |
| pglib_opf_case60_c__sad | 60 | 88 | inf. | 1.1350e+05 | 2.28 | 4.37 | <1 | <1 | <1 | <1 |
| pglib_opf_case73_ieee_rts__sad | 73 | 120 | 2.3268e+05 | 2.2760e+05 | 2.54 | 6.73 | <1 | <1 | <1 | <1 |
| pglib_opf_case89_pegase__sad | 89 | 210 | inf. | 1.0729e+05 | 0.71 | 0.73 | <1 | <1 | <1 | <1 |
| pglib_opf_case118_ieee__sad | 118 | 186 | inf. | 1.0516e+05 | 6.79 | 8.17 | <1 | <1 | <1 | <1 |
| pglib_opf_case162_ieee_dtc__sad | 162 | 284 | 1.0629e+05 | 1.0869e+05 | 6.25 | 6.48 | <1 | <1 | <1 | <1 |
| pglib_opf_case179_goc__sad | 179 | 263 | inf. | 7.6253e+05 | 1.01 | 1.12 | <1 | <1 | <1 | <1 |
| pglib_opf_case200_activ__sad | 200 | 245 | inf. | 2.7558e+04 | 0.01 | 0.01 | <1 | <1 | <1 | <1 |
| pglib_opf_case240_pserc__sad | 240 | 448 | inf. | 3.4054e+06 | 4.37 | 4.93 | <1 | 4 | 5 | 2 |
| pglib_opf_case300_ieee__sad | 300 | 411 | 5.2729e+05 | 5.6570e+05 | 2.43 | 2.61 | <1 | <1 | 2 | <1 |
| pglib_opf_case500_goc__sad | 500 | 733 | inf. | 4.8740e+05 | 6.45 | 6.67 | <1 | 2 | 3 | <1 |
| pglib_opf_case588_sdet__sad | 588 | 686 | inf. | 3.2936e+05 | 5.98 | 6.67 | <1 | <1 | 2 | <1 |
| pglib_opf_case793_goc__sad | 793 | 913 | inf. | 2.8580e+05 | 6.26 | 7.97 | <1 | 2 | 3 | <1 |
| pglib_opf_case1354_pegase__sad | 1354 | 1991 | inf. | 1.2588e+06 | 1.53 | 1.57 | <1 | 5 | 8 | 4 |
| pglib_opf_case1888_rte__sad | 1888 | 2531 | 1.3532e+06 | 1.4139e+06 | 2.81 | 2.82 | <1 | 9 | 11 | 23 |
| pglib_opf_case1951_rte__sad | 1951 | 2596 | inf. | 2.0924e+06 | 0.42 | 0.46 | 2 | 16 | 12 | 7 |
| pglib_opf_case2000_goc__sad | 2000 | 3639 | inf. | 9.9288e+05 | 1.35 | 1.52 | 2 | 7 | 14 | 6 |
| pglib_opf_case2312_goc__sad | 2312 | 3013 | inf. | 4.6235e+05 | 3.40 | 3.93 | 3 | 7 | 12 | 5 |
| pglib_opf_case2383wp_k__sad | 2383 | 2896 | inf. | 1.9112e+06 | 1.93 | 2.86 | 3 | 8 | 12 | 7 |
| pglib_opf_case2736sp_k__sad | 2736 | 3504 | inf. | 1.3266e+06 | 1.32 | 1.58 | 2 | 7 | 14 | 6 |
| pglib_opf_case2737sop_k__sad | 2737 | 3506 | inf. | 7.9094e+05 | 1.70 | 1.88 | 2 | 6 | 12 | 5 |
| pglib_opf_case2742_goc__sad | 2742 | 4673 | 2.5970e+05 | 2.7571e+05 | 1.34 | 1.36 | <1 | 27 | 19 | 7 |
| pglib_opf_case2746wop_k__sad | 2746 | 3514 | inf. | 1.2337e+06 | 1.97 | 2.32 | 2 | 6 | 11 | 5 |
| pglib_opf_case2746wp_k__sad | 2746 | 3514 | inf. | 1.6669e+06 | 1.65 | 2.19 | 2 | 6 | 12 | 7 |
| pglib_opf_case2848_rte__sad | 2848 | 3776 | inf. | 1.2890e+06 | 0.23 | 0.25 | 2 | 16 | 17 | 8 |
| pglib_opf_case2853_sdet__sad | 2853 | 3921 | inf. | 2.0692e+06 | 1.65 | 1.70 | 3 | 9 | 14 | 8 |
| pglib_opf_case2868_rte__sad | 2868 | 3808 | inf. | 2.0213e+06 | 0.56 | 0.59 | 2 | 16 | 18 | 9 |
| pglib_opf_case2869_pegase__sad | 2869 | 4582 | inf. | 2.4687e+06 | 1.02 | 1.12 | 4 | 11 | 25 | 10 |
| pglib_opf_case3012wp_k__sad | 3012 | 3572 | inf. | 2.6195e+06 | 1.35 | 1.56 | 3 | 10 | 18 | 17 |
| pglib_opf_case3022_goc__sad | 3022 | 4135 | 5.9922e+05 | 6.0143e+05 | 2.76 | 2.77 | <1 | 10 | 17 | 6 |
| pglib_opf_case3120sp_k__sad | 3120 | 3693 | inf. | 2.1749e+06 | 1.38 | 1.52 | 3 | 11 | 22 | 9 |
| pglib_opf_case3375wp_k__sad | 3374 | 4161 | 7.3196e+06 | 7.4382e+06 | 0.50 | 0.55 | <1 | 10 | 82 | 9 |
| pglib_opf_case3970_goc__sad | 3970 | 6641 | inf. | 9.6555e+05 | 0.63 | 0.76 | 5 | 21 | 30 | 16 |
| pglib_opf_case4020_goc__sad | 4020 | 6988 | inf. | 8.8969e+05 | 8.45 | 8.66 | 5 | 26 | 39 | 17 |
| pglib_opf_case4601_goc__sad | 4601 | 7199 | 1.1956e+06 | 8.7818e+05 | 6.25 | 6.37 | 3 | 28 | 36 | 18 |
| pglib_opf_case4619_goc__sad | 4619 | 8150 | inf. | 4.8435e+05 | 1.95 | 2.00 | 7 | 20 | 40 | 21 |
| pglib_opf_case4661_sdet__sad | 4661 | 5997 | inf. | 2.2610e+06 | 1.79 | 1.96 | 2 | 15 | 33 | 16 |
| pglib_opf_case4837_goc__sad | 4837 | 7765 | inf. | 8.7712e+05 | 0.73 | 0.85 | 8 | 20 | 35 | 14 |
| pglib_opf_case4917_goc__sad | 4917 | 6726 | 1.3843e+06 | 1.3890e+06 | 2.45 | 2.52 | <1 | 17 | 32 | 12 |
| pglib_opf_case6468_rte__sad | 6468 | 9000 | 1.9828e+06 | 2.0697e+06 | 1.10 | 1.12 | 2 | 60 | 90 | 43 |
| pglib_opf_case6470_rte__sad | 6470 | 9005 | 2.1391e+06 | 2.2416e+06 | 1.87 | 1.91 | 2 | 35 | 55 | 30 |
| pglib_opf_case6495_rte__sad | 6495 | 9019 | 2.5618e+06 | 3.0678e+06 | 14.99 | 15.11 | 2 | 65 | 71 | 32 |
| pglib_opf_case6515_rte__sad | 6515 | 9037 | 2.5595e+06 | 2.8699e+06 | 7.82 | 7.85 | 2 | 57 | 60 | 30 |
| pglib_opf_case8387_pegase__sad | 8387 | 14561 | 2.5977e+06 | 2.8039e+06 | 47.49 | 59.72 | 4 | 46 | 117 | 38 |
| pglib_opf_case9241_pegase__sad | 9241 | 16049 | inf. | 6.3185e+06 | 2.40 | 2.47 | 9 | 51 | 118 | 48 |
| pglib_opf_case9591_goc__sad | 9591 | 15915 | inf. | 1.1674e+06 | 9.35 | 9.56 | 13 | 79 | 142 | 57 |
| pglib_opf_case10000_goc__sad | 10000 | 13193 | inf. | 1.4902e+06 | 5.35 | 6.77 | 19 | 42 | 89 | 36 |
| pglib_opf_case10480_goc__sad | 10480 | 18559 | inf. | 2.3147e+06 | 1.14 | 1.23 | 26 | 61 | 132 | 54 |
| pglib_opf_case13659_pegase__sad | 13659 | 20467 | inf. | 9.0422e+06 | 1.63 | 1.68 | 19 | 61 | 161 | 60 |
| pglib_opf_case19402_goc__sad | 19402 | 34704 | 1.9102e+06 | 1.9838e+06 | 1.44 | 1.49 | 11 | 149 | 372 | 139 |
| pglib_opf_case24464_goc__sad | 24464 | 37816 | 2.5549e+06 | 2.6540e+06 | 1.77 | 1.84 | 12 | 148 | 381 | 129 |
| pglib_opf_case30000_goc__sad | 30000 | 35393 | 1.2532e+06 | 1.2686e+06 | 9.76 | 11.05 | 10 | 408 | 509 | 169 |


