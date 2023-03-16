# Sample Results

```matlab
K = 20;
scp_iters = 10;
wvc = 50;
wtr = 0.1;

Kstr = round(K/3);
c_d = 0.01;
rmax = 10;
vmax = 5;
umax = 8;
```

Sub-problem type: SOCP with trust-region penalty formulated as SOCs.

Potential objectives:
 - Maximize time available to defer decision.
 - Minimize time of maneuver for one of the targets.
 - Minimize cumulative fuel consumption.

## 2D -- 2 targets

```matlab
cost_factor = 0.3;
cost_bound = [20,50];
g = [0;0];
rK = [2   5;
      7   0];        
vK = 0.1*[1  0;
          0  1];
```

## 2D -- 3 targets

```matlab
cost_factor = 0.2;
cost_bound = [20,70,50]
g = [0;0];
rK = [2   5   0;
      7   0   10];
vK = 0.1*[1  0  0;
          0  1  0];
```

## 3D -- 3 targets 

```matlab
cost_factor = 0.4;
cost_bound = [50,30,40];
g = [0;0;0.5];
rK = [2   5   0;
      2   0   10;
      8   1   3];              
vK = 0.1*[1  0  0;
          0  1  0;
          0  0  1];
```
