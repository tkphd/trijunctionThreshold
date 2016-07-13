# Vertex resolution in triangular grain collapse

## run1: TK_R3_p139
* created ideal/vertex/run1/vertex.dat
* width = 14.5
* dt = 0.01
* t = 150,000&times;dt
* 6 cores: runtime 9829 seconds
* Grain 0 collapsed at t = 110,141&times;dt, v = 71&times;dV

## run2: TK_R3_p141
* created ideal/vertex/run2/vertex.dat
* width = 2
* dt = 0.01
* t = 1,000,000&times;dt
* 6 cores: runtime 54059 seconds
* Grain 0 collapsed at t = 0&times;dt, v = 0&times;dV
* Note: &delta;=2 *does not* establish a diffuse interface! The bottom edge vanished.

## run3: TK_R3_p141
* created ideal/vertex/run3/vertex.dat
* width = 3
* dt = 0.01
* t = 1,100,000&times;dt
* 6 cores: runtime 45277 seconds
* Grain 0 collapsed at t = 1,068,077&times;dt, v = 4.90118&times;dV

## run4: TK_R3_p141
* created ideal/vertex/run4/vertex.dat
* width = 4.0
* dt = 0.01
* t = 550,000&times;dt
* 6 cores: runtime 31107 seconds
* Grain 0 collapsed at t = 521,720&times;dt, v = 6&times;dV

## run6: TK_R3_p141
* created ideal/vertex/run6/vertex.dat
* width = 6.0
* dt = 0.01
* t = 300,000&times;dt
* 6 cores: runtime 17487 seconds
* Grain 0 collapsed at t = 286,363&times;dt, v = 13&times;dV

## run8: TK_R3_p141
* created ideal/vertex/run8/vertex.dat
* width = 8.0
* dt = 0.01
* t = 250,000&times;dt
* 6 cores: runtime 14778 seconds
* Grain 0 collapsed at t = 206,493&times;dt, v = 24&times;dV

## run10: TK_R3_p141
* created ideal/vertex/run10/vertex.dat
* width = 10.0
* dt = 0.01
* t = 200,000&times;dt
* 6 cores: runtime 12292 seconds
* Grain 0 collapsed at t = 162,713&times;dt, v = 39&times;dV

## run12: TK_R3_p141
* created ideal/vertex/run12/vertex.dat
* width = 12.0
* dt = 0.01
* t = 200,000&times;dt
* 6 cores: runtime 12665 seconds
* Grain 0 collapsed at t = 134,319&times;dt, v = 54&times;dV

## run14: TK_R3_p141
* created ideal/vertex/run14/vertex.dat
* width = 14.0
* dt = 0.01
* t = 200,000&times;dt
* 6 cores: runtime 13082 seconds
* Grain 0 collapsed at t = 114,270&times;dt, v = 76&times;dV

## run16: TK_R3_p141
* created ideal/vertex/run16/vertex.dat
* width = 16.0
* dt = 0.01
* t = 200,000&times;dt
* 6 cores: runtime 13392 seconds
* Grain 0 collapsed at t = 99,335&times;dt, v = 97&times;dV



| &delta; (&Delta;x) | Steps (&Delta;t) | Area (&Delta;x²) | Separation (&Delta;x) | Time/100k (sec) | Ideal Time (sec) |
| -----------------: | ---------------: | ---------------: | --------------------: | --------------: | ---------------: |
|  2.0               |      --          | --               | --                    | 5,406           |     --           |
|  3.0               | 1,068,077        | 4.90118          |  3.364341             | 4,116           | 43,962           |
|  4.0               |   521,720        | 6.41772          |  3.849817             | 5,655           | 29,503           |
|  6.0               |   286,363        | 13.8114          |  5.64766              | 5,873           | 16,818           |
|  8.0               |   206,493        | 24.5994          |  7.537233             | 5,911           | 12,205           |
| 10.0               |   162,713        | 39.9948          |  9.610621             | 6,146           | 10,000           |
| 12.0               |   134,319        | 54.1383          | 11.18155              | 6,332           |  8,505           |
| 14.0               |   114,270        | 76.4811          | 13.29005              | 6,541           |  7,474           |
| 14.5               |   110,141        | 71.3438          | 12.83594              | 6,552           |  7,216           |
| 16.0               |    99,335        | 97.773           | 15.02655              | 6,696           |  6,651           |
