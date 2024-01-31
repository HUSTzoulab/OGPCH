# OGPCH

**MATLAB code for genetic marker selection using linear programming and Python code for evaluating selected genes.**

The version of MATLAB is R2020a, the required toolkit: cplex12.10  
The version of Python is 3.10.2, the required module: numpy, pandas, scipy, sklearn, matplotlib, seaborn

There are four main functions for finding markers using MATLAB.  
- First, OGPCH.pre_OGPCH_FLAT(dataset,cons=5000), data preprocessing and saving the variables required for subsequent solving.
  - dataset: the name of gene dataset
  - cons: the number of constraints(default 5000)

- Second, OGPCH.OGPCH_FLAT(dataset,mu), Solving linear programming models.
  - dataset: the name of gene dataset
  - mu: the coefficient of the error in the objective function

- Third, OGPCH.pre_OGPCH_HIE(dataset,cons1=5000,cons2=5000), data preprocessing and saving the variables required for subsequent solving.
  - dataset: the name of gene dataset
  - cons1: the number of constraints(considering high-level cell labels, default 5000)
  - cons2: the number of constraints(considering low-level cell labels, default 5000)

- Fourth, OGPCH.OGPCH_HIE(dataset,mu,nu), Solving linear programming models.
  - dataset: the name of gene dataset
  - mu: the coefficient of the error considering low-level labels in the objective function
  - nu: the coefficient of the error considering high-level labels in the objective function


## PBMC3K example

Data from https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz.  


```matlab
# finding markers
pre_OGPCH_FLAT('PBMC3K',5000);
OGPCH_FLAT('PBMC3K',0.0074);
```

```python
get_accuracy('PBMC3K','FLAT','knn','OGPCH')
accuracy_line_chart('PBMC3K','knn')
```
![png](https://github.com/HUSTzoulab/OGPCH/blob/main/pictures/PBMC3K/PBMC3K_knn_Accuracy_Compare.png)

```python
kappa(('PBMC3K','FLAT','knn','OGPCH')
kappa_line_chart('PBMC3K','knn')
```
![png](https://github.com/HUSTzoulab/OGPCH/blob/main/pictures/PBMC3K/PBMC3K_knn_Kappa_Compare.png)

To prove the robustness of random selection constraints, randomly selecte 100 constraints and calculate the kappa coefficients for classification.

```python
tuning_kappa_sam('PBMC3K', 'FLAT', 'knn')
robust_voilin('PBMC3K')
```
![png](https://github.com/HUSTzoulab/OGPCH/blob/main/pictures/PBMC3K/PBMC3K_violin.png)

#TSNE plot

```python
Compare_25_plot('PBMC3K')
```
![png](https://github.com/HUSTzoulab/OGPCH/blob/main/pictures/PBMC3K/PBMC3K_compare_25.png)




## zeisel example

Data from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60631.  
This is a hierarchical dataset.


```matlab
# finding markers
pre_OGPCH_HIE('zeisel',5000);
OGPCH_HIE('zeisel',0.0074);
```

```python
get_accuracy('zeisel','HIE','knn','OGPCH')
accuracy_line_chart('zeisel','knn')
```
![png](https://github.com/HUSTzoulab/OGPCH/blob/main/pictures/zeisel/zeisel_knn_Accuracy_Compare.png)

```python
kappa(('zeisel','HIE','knn','OGPCH')
kappa_line_chart('zeisel','knn')
```
![png](https://github.com/HUSTzoulab/OGPCH/blob/main/pictures/zeisel/zeisel_knn_Kappa_Compare.png)

To prove the robustness of random selection constraints, randomly selecte 100 constraints and calculate the kappa coefficients for classification.
```python
tuning_kappa_sam('zeisel', 'HIE', 'knn')
robust_voilin('zeisel')
```
![png](https://github.com/HUSTzoulab/OGPCH/blob/main/pictures/zeisel/zeisel_violin.png)

Restore the original hierarchical structure using 45 genes with a 100% accuracy rate.
```python
tree_plot('zeisel')
```
![png](https://github.com/HUSTzoulab/OGPCH/blob/main/pictures/zeisel/zeisel_tree.png)


#TSNE plot
```python
Compare_25_plot('zeisel')
```
![png](https://github.com/HUSTzoulab/OGPCH/blob/main/pictures/zeisel/zeisel_compare_25.png)









