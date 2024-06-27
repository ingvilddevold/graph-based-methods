This repository holds the MATLAB code developed during my master's project on hybrid, graph-based methods for reduced-order and data-driven reservoir modeling.

## Installation ##
1. Install [MRST](http://www.mrst.no). 
2. Clone this repository, e.g., into the modules folder of MRST.
3. Create a file named `startup_user.m` within the MRST folder, at the same level as `startup.m`, and add the line
```
mrstPath('register', 'graph-based-methods', 'path/to/repo/graph-based-methods')
```

### Use ###

After following the steps above, the module can be used like any other MRST module. That is, first start MRST by running the `startup.m` file. Then the module is activated by the command 

~~~~~
mrstModule add graph-based-methods 
~~~~~
