# OGPCH

Matlab code for genetic marker selection using linear programming and Python code for evaluating selected genes.

The version of Matlab is R2020a, the required toolkit: cplex12.10

There are two main functions for finding markers using MATLAB. 
First, OGPCH.pre_OGPCH_HIE(dataset,cons1=5000, cons2=5000), data preprocessing and saving the variables required for subsequent solving.
- dataset: the name of gene dataset
- cons1: the number of constraints(considering high-level cell labels, default 5000)
- cons2: the number of constraints(considering low-level cell labels, default 5000)

second, OGPCH.OGPCH_HIE(dataset,mu,nu), Solving linear programming models。
- dataset: the name of gene dataset
- mu: the coefficient of the error considering low-level labels in the objective function
- nu: the coefficient of the error considering high-level labels in the objective function







The version of Python is 3.10.2, the required module: numpy, pandas, scipy, sklearn, matplotlib, seaborn
