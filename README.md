# OGPCH

**MATLAB code for genetic marker selection using linear programming and Python code for evaluating selected genes.**

The version of MATLAB is R2020a, the required toolkit: cplex12.10  
The version of Python is 3.10.2, the required module: numpy1.23.5, pandas1.5.3, scipy1.10.0, sklearn0.0.post5, matplotlib3.7.0, seaborn0.12.2

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

<div align=center><img src="https://z4a.net/images/2024/01/31/PBMC3K_knn_Accuracy_Compare.png" width = "400" height = "270"  /></div>

```python
kappa(('PBMC3K','FLAT','knn','OGPCH')
kappa_line_chart('PBMC3K','knn')
```
<div align=center><img src="https://z4a.net/images/2024/01/31/PBMC3K_knn_Kappa_Compare.png" width = "400" height = "270"  /></div>

To prove the robustness of random selection constraints, randomly selecte 100 constraints and calculate the kappa coefficients for classification.

```python
tuning_kappa_sam('PBMC3K', 'FLAT', 'knn')
robust_voilin('PBMC3K')
```
<div align=center><img src="https://z4a.net/images/2024/01/31/PBMC3K_violin.png" width = "200" height = "300"  /></div>



```python
#TSNE plot
Compare_25_plot('PBMC3K')
```
![png](https://z4a.net/images/2024/01/31/PBMC3K_compare_25.png)




## zeisel example

Data from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60631.  
This is a hierarchical dataset.


```matlab
# finding markers
pre_OGPCH_HIE('zeisel',5000,5000);
OGPCH_HIE('zeisel',0.000006,0.000006);
```

```python
get_accuracy('zeisel','HIE','knn','OGPCH')
accuracy_line_chart('zeisel','knn')
```
<div align=center><img src="https://z4a.net/images/2024/01/31/zeisel_knn_Accuracy_Compare.png" width = "400" height = "270"  /></div>

```python
kappa('zeisel','HIE','knn','OGPCH')
kappa_line_chart('zeisel','knn')
```
<div align=center><img src="https://z4a.net/images/2024/01/31/zeisel_knn_Kappa_Compare.png" width = "400" height = "270"  /></div>

To prove the robustness of random selection constraints, randomly selecte 100 constraints and calculate the kappa coefficients for classification.
```python
tuning_kappa_sam('zeisel', 'HIE', 'knn')
robust_voilin('zeisel')
```
<div align=center><img src="https://z4a.net/images/2024/01/31/zeisel_violin.png" width = "200" height = "300"  /></div>

Restore the original hierarchical structure using 45 genes with a 100% accuracy rate.
```python
tree_plot('zeisel')
```
<div align=center><img src="https://z4a.net/images/2024/01/31/zeisel_tree.png" width = "600" height = "700"  /></div>



```python
#TSNE plot
Compare_25_plot('zeisel')
```
![png](https://z4a.net/images/2024/01/31/zeisel_compare_25.png)




## Mouse_E9.5 example

Data from https://tome.gs.washington.edu/.    
This is a huge hierarchical dataset with over one hundred thousand cells. OGPCH can find its markers very quickly.


```matlab
# finding markers
pre_OGPCH_HIE('Mouse_E9_5',5000,5000);
OGPCH_HIE('Mouse_E9_5',0.0003,0.00001);
```

```python
get_accuracy('Mouse_E9_5','HIE','knn','OGPCH')
accuracy_line_chart('Mouse_E9_5','knn')
```
<div align=center><img src="https://z4a.net/images/2024/01/31/Mouse_E9_5_knn_Accuracy_Compare.png" width = "400" height = "270"  /></div>

```python
kappa('Mouse_E9_5','HIE','knn','OGPCH')
kappa_line_chart('Mouse_E9_5','knn')
```
<div align=center><img src="https://z4a.net/images/2024/01/31/Mouse_E9_5_knn_Kappa_Compare.png" width = "400" height = "270"  /></div>

To prove the robustness of random selection constraints, randomly selecte 100 constraints and calculate the kappa coefficients for classification.
```python
tuning_kappa_sam('Mouse_E9_5', 'HIE', 'knn')
robust_voilin('Mouse_E9_5')
```
<div align=center><img src="https://z4a.net/images/2024/01/31/Mouse_E9_5_violin.png" width = "200" height = "300"  /></div>

Restore the original hierarchical structure using 100 genes with a good accuracy rate.
```python
tree_plot('Mouse_E9_5')
```
<div align=center><img src="https://z4a.net/images/2024/01/31/Mouse_E9_5_tree.png" width = "600" height = "700"  /></div>



In order to show the classification effect more clearly the original samples are sampled, due to the extreme imbalance of the samples.
```python
#TSNE plot
Compare_25_plot('Mouse_E9_5')
```
![png](https://z4a.net/images/2024/01/31/Mouse_E9_5_compare_25_sam.png)










