# Vertex resolution in triangular grain collapse

## run1: TK_R3_p139
* created drag/vertex/run1/vertex.dat
* Annulus: x = 0.6422, s = 0.0175
* m<sub>tj</sub> = 0.01
* width = 14.5
* dt = 0.01
* t = 140,000&times;dt
* 6 cores: runtime 9600 seconds
* Grain 0 collapsed at t = 134,986&times;dt, V = 64&times;dV

## run10: TK_R3_p141
* created drag/vertex/run10/vertex.dat
* Annulus: x = 0.6422, s = 0.0175
* m<sub>tj</sub> = 0.01
* width = 10.0
* dt = 0.01
* t = 200,000&times;dt
* 6 cores: runtime 13829 seconds
* Grain 0 collapsed at t = 191,383&times;dt, V = 31&times;dV

## run12: TK_R3_p141
* created drag/vertex/run12/vertex.dat
* Annulus: x = 0.6422, s = 0.0175
* m<sub>tj</sub> = 0.01
* width = 12.0
* dt = 0.01
* t = 200,000&times;dt
* 6 cores: runtime 14191 seconds
* Grain 0 collapsed at t = 159,939&times;dt, V = 49&times;dV

## run14: TK_R3_p141
* created drag/vertex/run14/vertex.dat
* Annulus: x = 0.6422, s = 0.0175
* m<sub>tj</sub> = 0.01
* width = 14.0
* dt = 0.01
* t = 200,000&times;dt
* 6 cores: runtime 14575 seconds
* Grain 0 collapsed at t = 139,415&times;dt, V = 65&times;dV

## run16: TK_R3_p141
* created drag/vertex/run16/vertex.dat
* Annulus: x = 0.6422, s = 0.0175
* m<sub>tj</sub> = 0.01
* width = 16.0
* dt = 0.01
* t = 200,000&times;dt
* 6 cores: runtime 14895 seconds
* Grain 0 collapsed at t = 125,618&times;dt, V = 87&times;dV


## run4: TK_R3_p141
* created drag/vertex/run4/vertex.dat
* Annulus: x = 0.6422, s = 0.0175
* m<sub>tj</sub> = 0.01
* width = 4.0
* dt = 0.01
* t = 550,000&times;dt
* 6 cores: runtime 33997 seconds
* Grain 0 collapsed at t = 515,790&times;dt, V = 5&times;dV

## run6: TK_R3_p141
* created drag/vertex/run6/vertex.dat
* Annulus: x = 0.6422, s = 0.0175
* m<sub>tj</sub> = 0.01
* width = 6.0
* dt = 0.01
* t = 350,000&times;dt
* 6 cores: runtime 21969 seconds
* Grain 0 collapsed at t = 312,306&times;dt, V = 13&times;dV

## run8: TK_R3_p141
* created drag/vertex/run8/vertex.dat
* Annulus: x = 0.6422, s = 0.0175
* m<sub>tj</sub> = 0.01
* width = 8.0
* dt = 0.01
* t = 250,000&times;dt
* 6 cores: runtime 15365 seconds
* Grain 0 collapsed at t = 233,302&times;dt, V = 20&times;dV

## run2: TK_R3_p141
* created drag/vertex/run2/vertex.dat
* Annulus: x = 0.6422, s = 0.0175
* m<sub>tj</sub> = 0.01
* width = 2.0
* dt = 0.01
* t = 1,000,000&times;dt
* 6 cores: runtime 45973 seconds
* Grain 0 collapsed at t = &times;dt, V = &times;dV
* Note: &delta;=2 *does not* establish a diffuse interface! The bottom edge vanished.

## run3: TK_R3_p141
* created drag/vertex/run3/vertex.dat
* Annulus: x = 0.6422, s = 0.0175
* m<sub>tj</sub> = 0.01
* width = 3.0
* dt = 0.01
* t = 1,100,000&times;dt
* 6 cores: runtime 48789 seconds
* Grain 0 collapsed at t = 1,084,350&times;dt, V = 4.46304&times;dV


## run125: TK_R3_p144
* created drag/vertex/run124/vertex.dat
* Annulus: x = 0.6422, s = 0.0175
* m<sub>tj</sub> = 0.01
* width = 12.5
* dt = 0.01
* t = 160,000&times;dt
* 6 cores: runtime 11494 seconds
* Grain 0 collapsed at t = 156,030&times;dt, V = 54.4&times;dV


# mmsp2graph results
| &delta; (&Delta;x) | Steps (&Delta;t) | Area (&Delta;x²) | Separation (&Delta;x) | Time/100k (sec) | Ideal Time (sec) |
| -----------------: | ---------------: | ---------------: | --------------------: | --------------: | ---------------: |
|  2.0               |      --          | --               | --                    | 4,597           | --               |
|  3.0               | 1,084,350        |  4.46304         |  3.210444             | 4,435           | 48,091           |
|  4.0               |   515,790        |  5.77896         |  3.653209             | 6,181           | 31,881           |
|  6.0               |   312,306        | 13.2454          |  5.530727             | 6,484           | 20,249           |
|  8.0               |   233,302        | 20.9526          |  6.956145             | 6,146           | 14,338           |
| 10.0               |   191,383        | 31.8567          |  8.57729              | 6,914           | 13,232           |
| 12.0               |   159,939        | 49.7543          | 10.71926              | 7,095           | 11,347           |
| 12.5               |   156,030        | 54.3842          | 11.20691              | 7,184           | 11,209           |
| 14.0               |   139,415        | 65.4545          | 12.29474              | 7,287           | 10,159           |
| 14.5               |   134,986        | 64.5379          | 12.20835              | 6,400           |  8,639           |
| 16.0               |   125,618        | 87.4255          | 14.20917              | 7,447           |  9,354           |



# mmsp2topo results
| &delta; (&Delta;x) | Steps (&Delta;t) | Area (&Delta;x²) | Separation (&Delta;x) |
| -----------------: | ---------------: | ---------------: | --------------------: |
|  2.0               |   --             | --               | --                    |
|  3.0               |   --             | --               | --                    |
|  4.0               |  515,657         |  9.29988         |  4.634345             |
|  6.0               |  312,321         | 12.0003          |  5.264362             |
|  8.0               |  233,339         | 18.1691          |  6.477634             |
| 10.0               |  191,459         | 24.4907          |  7.520562             |
| 12.0               |  160,088         | 34.8029          |  8.965147             |
| 14.0               |  139,558         | 44.4354          | 10.13011              |
| 16.0               |  125,834         | 57.0602          | 11.47932              |

