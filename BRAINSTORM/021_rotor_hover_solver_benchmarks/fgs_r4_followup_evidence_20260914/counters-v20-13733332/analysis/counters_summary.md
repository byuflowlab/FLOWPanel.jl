CLK_TCK = 100

## Task A: perf counters (one warmed prepared solve, FIFO-gated, user-space only)

| arm | cycles | instructions | IPC | cache-refs | cache-misses | miss ratio | task-clock (ms) | multiplexed events |
|---|---|---|---|---|---|---|---|---|
| j4-b1 | 119647448660 | 429941766403 | 3.593 | 11162174910 | 643295318 | 5.763% | 37562.31 | none (all 100%) |
| j64-b1 | 139893307008 | 456662547317 | 3.264 | 11172378918 | 667422025 | 5.974% | 44828.49 | none (all 100%) |

## Task B: stage thread activity (per arm, ordered by sequence)

### j4-b1

| sequence | stage | span_seconds | n_threads | total_cpu_ticks | busy_cpu_seconds | avg_active_threads | n_excluded(-1) | pct_incomplete |
|---|---|---|---|---|---|---|---|---|
| 2 | initialization | 0.998285 | 6 | 101 | 1.0100 | 1.012 | 0 | 0.00% |
| 3 | fmm | 0.245936 | 6 | 96 | 0.9600 | 3.903 | 0 | 0.00% |
| 4 | influence_mapping | 0.000486 | 6 | 0 | 0.0000 | 0.000 | 0 | 0.00% |
| 5 | residual | 0.028727 | 6 | 3 | 0.0300 | 1.044 | 0 | 0.00% |
| 6 | nearfield_update | 0.371375 | 6 | 37 | 0.3700 | 0.996 | 0 | 0.00% |
| 7 | fmm | 0.244003 | 6 | 96 | 0.9600 | 3.934 | 0 | 0.00% |
| 8 | influence_mapping | 0.000480 | 6 | 0 | 0.0000 | 0.000 | 0 | 0.00% |
| 9 | residual | 0.005310 | 6 | 1 | 0.0100 | 1.883 | 0 | 0.00% |
| 10 | nearfield_update | 0.347560 | 6 | 34 | 0.3400 | 0.978 | 0 | 0.00% |
| 11 | fmm | 0.244173 | 6 | 97 | 0.9700 | 3.973 | 0 | 0.00% |
| 12 | influence_mapping | 0.000478 | 6 | 0 | 0.0000 | 0.000 | 0 | 0.00% |
| 13 | residual | 0.005364 | 6 | 0 | 0.0000 | 0.000 | 0 | 0.00% |
| 14 | nearfield_update | 0.347034 | 6 | 35 | 0.3500 | 1.009 | 0 | 0.00% |
| 15 | fmm | 0.243709 | 6 | 95 | 0.9500 | 3.898 | 0 | 0.00% |
| 16 | influence_mapping | 0.000468 | 6 | 0 | 0.0000 | 0.000 | 0 | 0.00% |
| 17 | residual | 0.005301 | 6 | 1 | 0.0100 | 1.887 | 0 | 0.00% |
| 18 | nearfield_update | 0.348595 | 6 | 34 | 0.3400 | 0.975 | 0 | 0.00% |
| 19 | fmm | 0.243515 | 6 | 97 | 0.9700 | 3.983 | 0 | 0.00% |
| 20 | influence_mapping | 0.000466 | 6 | 0 | 0.0000 | 0.000 | 0 | 0.00% |
| 21 | residual | 0.005318 | 6 | 0 | 0.0000 | 0.000 | 0 | 0.00% |
| 22 | nearfield_update | 0.346978 | 6 | 35 | 0.3500 | 1.009 | 0 | 0.00% |
| 23 | fmm | 0.243868 | 6 | 96 | 0.9600 | 3.937 | 0 | 0.00% |
| 24 | influence_mapping | 0.000481 | 6 | 0 | 0.0000 | 0.000 | 0 | 0.00% |
| 25 | residual | 0.005343 | 6 | 2 | 0.0200 | 3.743 | 0 | 0.00% |
| 26 | nearfield_update | 0.346909 | 6 | 34 | 0.3400 | 0.980 | 0 | 0.00% |
| 27 | fmm | 0.290269 | 6 | 101 | 1.0100 | 3.480 | 0 | 0.00% |
| 28 | influence_mapping | 0.000477 | 6 | 0 | 0.0000 | 0.000 | 0 | 0.00% |
| 29 | residual | 0.005336 | 6 | 0 | 0.0000 | 0.000 | 0 | 0.00% |
| 30 | nearfield_update | 0.348025 | 6 | 35 | 0.3500 | 1.006 | 0 | 0.00% |
| 31 | fmm | 0.248022 | 6 | 98 | 0.9800 | 3.951 | 0 | 0.00% |
| 32 | influence_mapping | 0.000485 | 6 | 0 | 0.0000 | 0.000 | 0 | 0.00% |
| 33 | residual | 0.005371 | 6 | 0 | 0.0000 | 0.000 | 0 | 0.00% |
| 34 | nearfield_update | 0.350029 | 6 | 35 | 0.3500 | 1.000 | 0 | 0.00% |
| 35 | fmm | 0.245073 | 6 | 96 | 0.9600 | 3.917 | 0 | 0.00% |
| 36 | influence_mapping | 0.000477 | 6 | 0 | 0.0000 | 0.000 | 0 | 0.00% |
| 37 | residual | 0.005353 | 6 | 1 | 0.0100 | 1.868 | 0 | 0.00% |
| 38 | nearfield_update | 0.347011 | 6 | 35 | 0.3500 | 1.009 | 0 | 0.00% |
| 39 | fmm | 0.243219 | 6 | 95 | 0.9500 | 3.906 | 0 | 0.00% |
| 40 | influence_mapping | 0.000475 | 6 | 0 | 0.0000 | 0.000 | 0 | 0.00% |
| 41 | residual | 0.005359 | 6 | 0 | 0.0000 | 0.000 | 0 | 0.00% |
| 42 | nearfield_update | 0.346985 | 6 | 35 | 0.3500 | 1.009 | 0 | 0.00% |
| 43 | fmm | 0.244455 | 6 | 96 | 0.9600 | 3.927 | 0 | 0.00% |
| 44 | influence_mapping | 0.000477 | 6 | 0 | 0.0000 | 0.000 | 0 | 0.00% |
| 45 | residual | 0.005326 | 6 | 1 | 0.0100 | 1.878 | 0 | 0.00% |
| 46 | nearfield_update | 0.347115 | 6 | 34 | 0.3400 | 0.980 | 0 | 0.00% |
| 47 | fmm | 0.243964 | 6 | 97 | 0.9700 | 3.976 | 0 | 0.00% |
| 48 | influence_mapping | 0.000492 | 6 | 0 | 0.0000 | 0.000 | 0 | 0.00% |
| 49 | residual | 0.005338 | 6 | 0 | 0.0000 | 0.000 | 0 | 0.00% |
| 50 | nearfield_update | 0.347784 | 6 | 35 | 0.3500 | 1.006 | 0 | 0.00% |
| 51 | fmm | 0.244004 | 6 | 96 | 0.9600 | 3.934 | 0 | 0.00% |
| 52 | influence_mapping | 0.000488 | 6 | 0 | 0.0000 | 0.000 | 0 | 0.00% |
| 53 | residual | 0.005436 | 6 | 1 | 0.0100 | 1.840 | 0 | 0.00% |
| 54 | nearfield_update | 0.353020 | 6 | 35 | 0.3500 | 0.991 | 0 | 0.00% |
| 55 | fmm | 0.247219 | 6 | 96 | 0.9600 | 3.883 | 0 | 0.00% |
| 56 | influence_mapping | 0.000649 | 6 | 0 | 0.0000 | 0.000 | 0 | 0.00% |
| 57 | residual | 0.005382 | 6 | 0 | 0.0000 | 0.000 | 0 | 0.00% |
| 58 | nearfield_update | 0.348602 | 6 | 35 | 0.3500 | 1.004 | 0 | 0.00% |
| 59 | fmm | 0.244871 | 6 | 96 | 0.9600 | 3.920 | 0 | 0.00% |
| 60 | influence_mapping | 0.000540 | 6 | 0 | 0.0000 | 0.000 | 0 | 0.00% |
| 61 | residual | 0.005373 | 6 | 1 | 0.0100 | 1.861 | 0 | 0.00% |
| 62 | nearfield_update | 0.347070 | 6 | 34 | 0.3400 | 0.980 | 0 | 0.00% |
| 63 | fmm | 0.247267 | 6 | 97 | 0.9700 | 3.923 | 0 | 0.00% |
| 64 | influence_mapping | 0.000541 | 6 | 0 | 0.0000 | 0.000 | 0 | 0.00% |
| 65 | residual | 0.005357 | 6 | 0 | 0.0000 | 0.000 | 0 | 0.00% |
| 66 | nearfield_update | 0.347185 | 6 | 35 | 0.3500 | 1.008 | 0 | 0.00% |
| 67 | fmm | 0.244149 | 6 | 96 | 0.9600 | 3.932 | 0 | 0.00% |
| 68 | influence_mapping | 0.000536 | 6 | 0 | 0.0000 | 0.000 | 0 | 0.00% |
| 69 | residual | 0.005422 | 6 | 1 | 0.0100 | 1.844 | 0 | 0.00% |
| 70 | nearfield_update | 0.349398 | 6 | 35 | 0.3500 | 1.002 | 0 | 0.00% |
| 71 | fmm | 0.244675 | 6 | 96 | 0.9600 | 3.924 | 0 | 0.00% |
| 72 | influence_mapping | 0.000535 | 6 | 0 | 0.0000 | 0.000 | 0 | 0.00% |
| 73 | residual | 0.005380 | 6 | 1 | 0.0100 | 1.859 | 0 | 0.00% |
| 74 | nearfield_update | 0.348608 | 6 | 34 | 0.3400 | 0.975 | 0 | 0.00% |
| 75 | fmm | 0.244357 | 6 | 97 | 0.9700 | 3.970 | 0 | 0.00% |
| 76 | influence_mapping | 0.000532 | 6 | 0 | 0.0000 | 0.000 | 0 | 0.00% |
| 77 | residual | 0.005368 | 6 | 0 | 0.0000 | 0.000 | 0 | 0.00% |
| 78 | nearfield_update | 0.373163 | 6 | 37 | 0.3700 | 0.992 | 0 | 0.00% |
| 79 | fmm | 0.245280 | 6 | 95 | 0.9500 | 3.873 | 0 | 0.00% |
| 80 | influence_mapping | 0.000537 | 6 | 1 | 0.0100 | 18.637 | 0 | 0.00% |
| 81 | residual | 0.012390 | 6 | 1 | 0.0100 | 0.807 | 0 | 0.00% |
| 82 | nearfield_update | 0.347505 | 6 | 35 | 0.3500 | 1.007 | 0 | 0.00% |
| 83 | fmm | 0.246256 | 6 | 96 | 0.9600 | 3.898 | 0 | 0.00% |
| 84 | influence_mapping | 0.000515 | 6 | 0 | 0.0000 | 0.000 | 0 | 0.00% |
| 85 | residual | 0.005390 | 6 | 0 | 0.0000 | 0.000 | 0 | 0.00% |
| 86 | nearfield_update | 0.347872 | 6 | 35 | 0.3500 | 1.006 | 0 | 0.00% |
| 87 | fmm | 0.244582 | 6 | 96 | 0.9600 | 3.925 | 0 | 0.00% |
| 88 | influence_mapping | 0.000500 | 6 | 0 | 0.0000 | 0.000 | 0 | 0.00% |
| 89 | residual | 0.005422 | 6 | 1 | 0.0100 | 1.844 | 0 | 0.00% |
| 90 | nearfield_update | 0.347931 | 6 | 35 | 0.3500 | 1.006 | 0 | 0.00% |
| 91 | fmm | 0.244759 | 6 | 97 | 0.9700 | 3.963 | 0 | 0.00% |
| 92 | influence_mapping | 0.000487 | 6 | 0 | 0.0000 | 0.000 | 0 | 0.00% |
| 93 | residual | 0.005338 | 6 | 0 | 0.0000 | 0.000 | 0 | 0.00% |
| 94 | nearfield_update | 0.347042 | 6 | 35 | 0.3500 | 1.009 | 0 | 0.00% |
| 95 | fmm | 0.244785 | 6 | 96 | 0.9600 | 3.922 | 0 | 0.00% |
| 96 | influence_mapping | 0.000493 | 6 | 0 | 0.0000 | 0.000 | 0 | 0.00% |
| 97 | residual | 0.005354 | 6 | 1 | 0.0100 | 1.868 | 0 | 0.00% |
| 98 | nearfield_update | 0.347055 | 6 | 34 | 0.3400 | 0.980 | 0 | 0.00% |
| 99 | fmm | 0.243997 | 6 | 97 | 0.9700 | 3.975 | 0 | 0.00% |
| 100 | influence_mapping | 0.000505 | 6 | 0 | 0.0000 | 0.000 | 0 | 0.00% |
| 101 | residual | 0.005424 | 6 | 0 | 0.0000 | 0.000 | 0 | 0.00% |
| 102 | nearfield_update | 0.364530 | 6 | 37 | 0.3700 | 1.015 | 0 | 0.00% |
| 103 | fmm | 0.243826 | 6 | 95 | 0.9500 | 3.896 | 0 | 0.00% |
| 104 | influence_mapping | 0.000522 | 6 | 0 | 0.0000 | 0.000 | 0 | 0.00% |
| 105 | residual | 0.020670 | 6 | 2 | 0.0200 | 0.968 | 0 | 0.00% |
| 106 | nearfield_update | 0.347071 | 6 | 35 | 0.3500 | 1.008 | 0 | 0.00% |
| 107 | fmm | 0.245372 | 6 | 96 | 0.9600 | 3.912 | 0 | 0.00% |
| 108 | influence_mapping | 0.000512 | 6 | 0 | 0.0000 | 0.000 | 0 | 0.00% |
| 109 | residual | 0.005344 | 6 | 0 | 0.0000 | 0.000 | 0 | 0.00% |
| 110 | nearfield_update | 0.347922 | 6 | 35 | 0.3500 | 1.006 | 0 | 0.00% |
| 111 | fmm | 0.244406 | 6 | 96 | 0.9600 | 3.928 | 0 | 0.00% |
| 112 | influence_mapping | 0.000503 | 6 | 0 | 0.0000 | 0.000 | 0 | 0.00% |
| 113 | residual | 0.005448 | 6 | 1 | 0.0100 | 1.836 | 0 | 0.00% |
| 114 | strength_copy | 0.000664 | 6 | 0 | 0.0000 | 0.000 | 0 | 0.00% |

### j64-b1

| sequence | stage | span_seconds | n_threads | total_cpu_ticks | busy_cpu_seconds | avg_active_threads | n_excluded(-1) | pct_incomplete |
|---|---|---|---|---|---|---|---|---|
| 2 | initialization | 0.335353 | 96 | 42 | 0.4200 | 1.252 | 0 | 0.00% |
| 3 | fmm | 0.037532 | 96 | 118 | 1.1800 | 31.440 | 0 | 0.00% |
| 4 | influence_mapping | 0.000611 | 96 | 0 | 0.0000 | 0.000 | 0 | 0.00% |
| 5 | residual | 0.029198 | 96 | 3 | 0.0300 | 1.027 | 0 | 0.00% |
| 6 | nearfield_update | 0.382724 | 96 | 38 | 0.3800 | 0.993 | 0 | 0.00% |
| 7 | fmm | 0.038728 | 96 | 120 | 1.2000 | 30.985 | 0 | 0.00% |
| 8 | influence_mapping | 0.000665 | 96 | 0 | 0.0000 | 0.000 | 0 | 0.00% |
| 9 | residual | 0.004765 | 96 | 1 | 0.0100 | 2.099 | 0 | 0.00% |
| 10 | nearfield_update | 0.355774 | 96 | 35 | 0.3500 | 0.984 | 0 | 0.00% |
| 11 | fmm | 0.038927 | 96 | 121 | 1.2100 | 31.083 | 0 | 0.00% |
| 12 | influence_mapping | 0.000735 | 96 | 1 | 0.0100 | 13.607 | 0 | 0.00% |
| 13 | residual | 0.005438 | 96 | 0 | 0.0000 | 0.000 | 0 | 0.00% |
| 14 | nearfield_update | 0.353707 | 96 | 36 | 0.3600 | 1.018 | 0 | 0.00% |
| 15 | fmm | 0.039566 | 96 | 118 | 1.1800 | 29.824 | 0 | 0.00% |
| 16 | influence_mapping | 0.000630 | 96 | 0 | 0.0000 | 0.000 | 0 | 0.00% |
| 17 | residual | 0.004726 | 96 | 0 | 0.0000 | 0.000 | 0 | 0.00% |
| 18 | nearfield_update | 0.353281 | 96 | 36 | 0.3600 | 1.019 | 0 | 0.00% |
| 19 | fmm | 0.039708 | 96 | 121 | 1.2100 | 30.472 | 0 | 0.00% |
| 20 | influence_mapping | 0.000726 | 96 | 0 | 0.0000 | 0.000 | 0 | 0.00% |
| 21 | residual | 0.005392 | 96 | 1 | 0.0100 | 1.855 | 0 | 0.00% |
| 22 | nearfield_update | 0.354319 | 96 | 35 | 0.3500 | 0.988 | 0 | 0.00% |
| 23 | fmm | 0.039830 | 96 | 127 | 1.2700 | 31.886 | 0 | 0.00% |
| 24 | influence_mapping | 0.000707 | 96 | 0 | 0.0000 | 0.000 | 0 | 0.00% |
| 25 | residual | 0.005451 | 96 | 1 | 0.0100 | 1.835 | 0 | 0.00% |
| 26 | nearfield_update | 0.354565 | 96 | 35 | 0.3500 | 0.987 | 0 | 0.00% |
| 27 | fmm | 0.039655 | 96 | 124 | 1.2400 | 31.270 | 0 | 0.00% |
| 28 | influence_mapping | 0.000746 | 96 | 0 | 0.0000 | 0.000 | 0 | 0.00% |
| 29 | residual | 0.005354 | 96 | 1 | 0.0100 | 1.868 | 0 | 0.00% |
| 30 | nearfield_update | 0.354184 | 96 | 36 | 0.3600 | 1.016 | 0 | 0.00% |
| 31 | fmm | 0.040272 | 96 | 119 | 1.1900 | 29.549 | 0 | 0.00% |
| 32 | influence_mapping | 0.000666 | 96 | 1 | 0.0100 | 15.025 | 0 | 0.00% |
| 33 | residual | 0.004787 | 96 | 0 | 0.0000 | 0.000 | 0 | 0.00% |
| 34 | nearfield_update | 0.354239 | 96 | 36 | 0.3600 | 1.016 | 0 | 0.00% |
| 35 | fmm | 0.040641 | 96 | 122 | 1.2200 | 30.019 | 0 | 0.00% |
| 36 | influence_mapping | 0.000675 | 96 | 0 | 0.0000 | 0.000 | 0 | 0.00% |
| 37 | residual | 0.004788 | 96 | 1 | 0.0100 | 2.089 | 0 | 0.00% |
| 38 | nearfield_update | 0.353773 | 96 | 36 | 0.3600 | 1.018 | 0 | 0.00% |
| 39 | fmm | 0.040628 | 96 | 122 | 1.2200 | 30.028 | 0 | 0.00% |
| 40 | influence_mapping | 0.000747 | 96 | 0 | 0.0000 | 0.000 | 0 | 0.00% |
| 41 | residual | 0.005379 | 96 | 0 | 0.0000 | 0.000 | 0 | 0.00% |
| 42 | nearfield_update | 0.353800 | 96 | 36 | 0.3600 | 1.018 | 0 | 0.00% |
| 43 | fmm | 0.047758 | 96 | 123 | 1.2300 | 25.755 | 0 | 0.00% |
| 44 | influence_mapping | 0.000742 | 96 | 1 | 0.0100 | 13.484 | 0 | 0.00% |
| 45 | residual | 0.005485 | 96 | 1 | 0.0100 | 1.823 | 0 | 0.00% |
| 46 | nearfield_update | 0.354294 | 96 | 35 | 0.3500 | 0.988 | 0 | 0.00% |
| 47 | fmm | 0.039838 | 96 | 119 | 1.1900 | 29.871 | 0 | 0.00% |
| 48 | influence_mapping | 0.000665 | 96 | 0 | 0.0000 | 0.000 | 0 | 0.00% |
| 49 | residual | 0.004786 | 96 | 1 | 0.0100 | 2.089 | 0 | 0.00% |
| 50 | nearfield_update | 0.353933 | 96 | 35 | 0.3500 | 0.989 | 0 | 0.00% |
| 51 | fmm | 0.040325 | 96 | 123 | 1.2300 | 30.502 | 0 | 0.00% |
| 52 | influence_mapping | 0.000672 | 96 | 0 | 0.0000 | 0.000 | 0 | 0.00% |
| 53 | residual | 0.004806 | 96 | 0 | 0.0000 | 0.000 | 0 | 0.00% |
| 54 | nearfield_update | 0.353483 | 96 | 36 | 0.3600 | 1.018 | 0 | 0.00% |
| 55 | fmm | 0.039977 | 96 | 125 | 1.2500 | 31.268 | 0 | 0.00% |
| 56 | influence_mapping | 0.000690 | 96 | 0 | 0.0000 | 0.000 | 0 | 0.00% |
| 57 | residual | 0.004823 | 96 | 0 | 0.0000 | 0.000 | 0 | 0.00% |
| 58 | nearfield_update | 0.355189 | 96 | 36 | 0.3600 | 1.014 | 0 | 0.00% |
| 59 | fmm | 0.037393 | 96 | 121 | 1.2100 | 32.359 | 0 | 0.00% |
| 60 | influence_mapping | 0.000666 | 96 | 1 | 0.0100 | 15.022 | 0 | 0.00% |
| 61 | residual | 0.004897 | 96 | 1 | 0.0100 | 2.042 | 0 | 0.00% |
| 62 | nearfield_update | 0.355772 | 96 | 36 | 0.3600 | 1.012 | 0 | 0.00% |
| 63 | fmm | 0.037803 | 96 | 119 | 1.1900 | 31.479 | 0 | 0.00% |
| 64 | influence_mapping | 0.000691 | 96 | 0 | 0.0000 | 0.000 | 0 | 0.00% |
| 65 | residual | 0.005370 | 96 | 1 | 0.0100 | 1.862 | 0 | 0.00% |
| 66 | nearfield_update | 0.351660 | 96 | 35 | 0.3500 | 0.995 | 0 | 0.00% |
| 67 | fmm | 0.038490 | 96 | 127 | 1.2700 | 32.995 | 0 | 0.00% |
| 68 | influence_mapping | 0.000630 | 96 | 0 | 0.0000 | 0.000 | 0 | 0.00% |
| 69 | residual | 0.004812 | 96 | 1 | 0.0100 | 2.078 | 0 | 0.00% |
| 70 | nearfield_update | 0.353415 | 96 | 35 | 0.3500 | 0.990 | 0 | 0.00% |
| 71 | fmm | 0.037228 | 96 | 121 | 1.2100 | 32.502 | 0 | 0.00% |
| 72 | influence_mapping | 0.000622 | 96 | 0 | 0.0000 | 0.000 | 0 | 0.00% |
| 73 | residual | 0.004840 | 96 | 1 | 0.0100 | 2.066 | 0 | 0.00% |
| 74 | nearfield_update | 0.354545 | 96 | 35 | 0.3500 | 0.987 | 0 | 0.00% |
| 75 | fmm | 0.038347 | 96 | 116 | 1.1600 | 30.250 | 0 | 0.00% |
| 76 | influence_mapping | 0.000622 | 96 | 0 | 0.0000 | 0.000 | 0 | 0.00% |
| 77 | residual | 0.004843 | 96 | 1 | 0.0100 | 2.065 | 0 | 0.00% |
| 78 | nearfield_update | 0.353412 | 96 | 36 | 0.3600 | 1.019 | 0 | 0.00% |
| 79 | fmm | 0.038402 | 96 | 120 | 1.2000 | 31.248 | 0 | 0.00% |
| 80 | influence_mapping | 0.000603 | 96 | 0 | 0.0000 | 0.000 | 0 | 0.00% |
| 81 | residual | 0.004798 | 96 | 1 | 0.0100 | 2.084 | 0 | 0.00% |
| 82 | nearfield_update | 0.353378 | 96 | 35 | 0.3500 | 0.990 | 0 | 0.00% |
| 83 | fmm | 0.037429 | 96 | 121 | 1.2100 | 32.328 | 0 | 0.00% |
| 84 | influence_mapping | 0.000595 | 96 | 0 | 0.0000 | 0.000 | 0 | 0.00% |
| 85 | residual | 0.004829 | 96 | 1 | 0.0100 | 2.071 | 0 | 0.00% |
| 86 | nearfield_update | 0.353895 | 96 | 35 | 0.3500 | 0.989 | 0 | 0.00% |
| 87 | fmm | 0.038435 | 96 | 123 | 1.2300 | 32.002 | 0 | 0.00% |
| 88 | influence_mapping | 0.000578 | 96 | 1 | 0.0100 | 17.315 | 0 | 0.00% |
| 89 | residual | 0.004783 | 96 | 0 | 0.0000 | 0.000 | 0 | 0.00% |
| 90 | nearfield_update | 0.354298 | 96 | 37 | 0.3700 | 1.044 | 0 | 0.00% |
| 91 | fmm | 0.038631 | 96 | 122 | 1.2200 | 31.581 | 0 | 0.00% |
| 92 | influence_mapping | 0.000601 | 96 | 0 | 0.0000 | 0.000 | 0 | 0.00% |
| 93 | residual | 0.004835 | 96 | 0 | 0.0000 | 0.000 | 0 | 0.00% |
| 94 | nearfield_update | 0.354194 | 96 | 35 | 0.3500 | 0.988 | 0 | 0.00% |
| 95 | fmm | 0.038186 | 96 | 123 | 1.2300 | 32.211 | 0 | 0.00% |
| 96 | influence_mapping | 0.000597 | 96 | 0 | 0.0000 | 0.000 | 0 | 0.00% |
| 97 | residual | 0.004810 | 96 | 0 | 0.0000 | 0.000 | 0 | 0.00% |
| 98 | nearfield_update | 0.354277 | 96 | 36 | 0.3600 | 1.016 | 0 | 0.00% |
| 99 | fmm | 0.038417 | 96 | 122 | 1.2200 | 31.757 | 0 | 0.00% |
| 100 | influence_mapping | 0.000601 | 96 | 0 | 0.0000 | 0.000 | 0 | 0.00% |
| 101 | residual | 0.004793 | 96 | 1 | 0.0100 | 2.086 | 0 | 0.00% |
| 102 | nearfield_update | 0.353736 | 96 | 36 | 0.3600 | 1.018 | 0 | 0.00% |
| 103 | fmm | 0.038842 | 96 | 120 | 1.2000 | 30.894 | 0 | 0.00% |
| 104 | influence_mapping | 0.000633 | 96 | 0 | 0.0000 | 0.000 | 0 | 0.00% |
| 105 | residual | 0.004725 | 96 | 2 | 0.0200 | 4.233 | 0 | 0.00% |
| 106 | nearfield_update | 0.353148 | 96 | 35 | 0.3500 | 0.991 | 0 | 0.00% |
| 107 | fmm | 0.039873 | 96 | 128 | 1.2800 | 32.102 | 0 | 0.00% |
| 108 | influence_mapping | 0.000627 | 96 | 0 | 0.0000 | 0.000 | 0 | 0.00% |
| 109 | residual | 0.004814 | 96 | 0 | 0.0000 | 0.000 | 0 | 0.00% |
| 110 | nearfield_update | 0.354267 | 96 | 36 | 0.3600 | 1.016 | 0 | 0.00% |
| 111 | fmm | 0.038594 | 96 | 119 | 1.1900 | 30.833 | 0 | 0.00% |
| 112 | influence_mapping | 0.000636 | 96 | 0 | 0.0000 | 0.000 | 0 | 0.00% |
| 113 | residual | 0.004760 | 96 | 0 | 0.0000 | 0.000 | 0 | 0.00% |
| 114 | strength_copy | 0.000422 | 96 | 1 | 0.0100 | 23.680 | 0 | 0.00% |

## Task C: cross-checks

| arm | baseline diagnostic_s | activity diagnostic_s | counters diagnostic_s | sum(span_seconds) | span_sum - activity_diag |
|---|---|---|---|---|---|
| j4-b1 | 17.57 | 18.35 | 17.4 | 17.57 | -0.7855 |
| j64-b1 | 11.92 | 12.37 | 11.12 | 11.2 | -1.169 |

## Anomalies

- j4-b1: stage 'influence_mapping' (seq 4) has zero total cpu_ticks
- j4-b1: stage 'influence_mapping' (seq 8) has zero total cpu_ticks
- j4-b1: stage 'influence_mapping' (seq 12) has zero total cpu_ticks
- j4-b1: stage 'residual' (seq 13) has zero total cpu_ticks
- j4-b1: stage 'influence_mapping' (seq 16) has zero total cpu_ticks
- j4-b1: stage 'influence_mapping' (seq 20) has zero total cpu_ticks
- j4-b1: stage 'residual' (seq 21) has zero total cpu_ticks
- j4-b1: stage 'influence_mapping' (seq 24) has zero total cpu_ticks
- j4-b1: stage 'influence_mapping' (seq 28) has zero total cpu_ticks
- j4-b1: stage 'residual' (seq 29) has zero total cpu_ticks
- j4-b1: stage 'influence_mapping' (seq 32) has zero total cpu_ticks
- j4-b1: stage 'residual' (seq 33) has zero total cpu_ticks
- j4-b1: stage 'influence_mapping' (seq 36) has zero total cpu_ticks
- j4-b1: stage 'influence_mapping' (seq 40) has zero total cpu_ticks
- j4-b1: stage 'residual' (seq 41) has zero total cpu_ticks
- j4-b1: stage 'influence_mapping' (seq 44) has zero total cpu_ticks
- j4-b1: stage 'influence_mapping' (seq 48) has zero total cpu_ticks
- j4-b1: stage 'residual' (seq 49) has zero total cpu_ticks
- j4-b1: stage 'influence_mapping' (seq 52) has zero total cpu_ticks
- j4-b1: stage 'influence_mapping' (seq 56) has zero total cpu_ticks
- j4-b1: stage 'residual' (seq 57) has zero total cpu_ticks
- j4-b1: stage 'influence_mapping' (seq 60) has zero total cpu_ticks
- j4-b1: stage 'influence_mapping' (seq 64) has zero total cpu_ticks
- j4-b1: stage 'residual' (seq 65) has zero total cpu_ticks
- j4-b1: stage 'influence_mapping' (seq 68) has zero total cpu_ticks
- j4-b1: stage 'influence_mapping' (seq 72) has zero total cpu_ticks
- j4-b1: stage 'influence_mapping' (seq 76) has zero total cpu_ticks
- j4-b1: stage 'residual' (seq 77) has zero total cpu_ticks
- j4-b1: stage 'influence_mapping' (seq 84) has zero total cpu_ticks
- j4-b1: stage 'residual' (seq 85) has zero total cpu_ticks
- j4-b1: stage 'influence_mapping' (seq 88) has zero total cpu_ticks
- j4-b1: stage 'influence_mapping' (seq 92) has zero total cpu_ticks
- j4-b1: stage 'residual' (seq 93) has zero total cpu_ticks
- j4-b1: stage 'influence_mapping' (seq 96) has zero total cpu_ticks
- j4-b1: stage 'influence_mapping' (seq 100) has zero total cpu_ticks
- j4-b1: stage 'residual' (seq 101) has zero total cpu_ticks
- j4-b1: stage 'influence_mapping' (seq 104) has zero total cpu_ticks
- j4-b1: stage 'influence_mapping' (seq 108) has zero total cpu_ticks
- j4-b1: stage 'residual' (seq 109) has zero total cpu_ticks
- j4-b1: stage 'influence_mapping' (seq 112) has zero total cpu_ticks
- j4-b1: stage 'strength_copy' (seq 114) has zero total cpu_ticks
- j64-b1: stage 'influence_mapping' (seq 4) has zero total cpu_ticks
- j64-b1: stage 'influence_mapping' (seq 8) has zero total cpu_ticks
- j64-b1: stage 'residual' (seq 13) has zero total cpu_ticks
- j64-b1: stage 'influence_mapping' (seq 16) has zero total cpu_ticks
- j64-b1: stage 'residual' (seq 17) has zero total cpu_ticks
- j64-b1: stage 'influence_mapping' (seq 20) has zero total cpu_ticks
- j64-b1: stage 'influence_mapping' (seq 24) has zero total cpu_ticks
- j64-b1: stage 'influence_mapping' (seq 28) has zero total cpu_ticks
- j64-b1: stage 'residual' (seq 33) has zero total cpu_ticks
- j64-b1: stage 'influence_mapping' (seq 36) has zero total cpu_ticks
- j64-b1: stage 'influence_mapping' (seq 40) has zero total cpu_ticks
- j64-b1: stage 'residual' (seq 41) has zero total cpu_ticks
- j64-b1: stage 'influence_mapping' (seq 48) has zero total cpu_ticks
- j64-b1: stage 'influence_mapping' (seq 52) has zero total cpu_ticks
- j64-b1: stage 'residual' (seq 53) has zero total cpu_ticks
- j64-b1: stage 'influence_mapping' (seq 56) has zero total cpu_ticks
- j64-b1: stage 'residual' (seq 57) has zero total cpu_ticks
- j64-b1: stage 'influence_mapping' (seq 64) has zero total cpu_ticks
- j64-b1: stage 'influence_mapping' (seq 68) has zero total cpu_ticks
- j64-b1: stage 'influence_mapping' (seq 72) has zero total cpu_ticks
- j64-b1: stage 'influence_mapping' (seq 76) has zero total cpu_ticks
- j64-b1: stage 'influence_mapping' (seq 80) has zero total cpu_ticks
- j64-b1: stage 'influence_mapping' (seq 84) has zero total cpu_ticks
- j64-b1: stage 'residual' (seq 89) has zero total cpu_ticks
- j64-b1: stage 'influence_mapping' (seq 92) has zero total cpu_ticks
- j64-b1: stage 'residual' (seq 93) has zero total cpu_ticks
- j64-b1: stage 'influence_mapping' (seq 96) has zero total cpu_ticks
- j64-b1: stage 'residual' (seq 97) has zero total cpu_ticks
- j64-b1: stage 'influence_mapping' (seq 100) has zero total cpu_ticks
- j64-b1: stage 'influence_mapping' (seq 104) has zero total cpu_ticks
- j64-b1: stage 'influence_mapping' (seq 108) has zero total cpu_ticks
- j64-b1: stage 'residual' (seq 109) has zero total cpu_ticks
- j64-b1: stage 'influence_mapping' (seq 112) has zero total cpu_ticks
- j64-b1: stage 'residual' (seq 113) has zero total cpu_ticks
