# Statmodel

Modification of the algorithm for Statistical Modeling and Analysis of Experiments without ANOVA [1] applied to 
gene regulatory network (GRN) inference. 
 
Tested on MATLAB R2016a with the Parallel Computing Toolbox. 

Required Input 
--------------

* genes = cell array (# of genes,1) of genes names.
* regulators = cell array (# of regulators, 1)  of regulators names.
* expressiondata = numeric matrix (# of genes, # of conditions) of the gene expression data.

Output
------

* Net = Inferred network write as 3-column cell array 
	* Column 1 = Regulator 
	* Column 2 = Target Gene 
	* Column 3 = score of the interaction (the higher the score, the more reliable the interaction)
		 
Bibliography		 
------------

[1] Hernandez, Hugo. "Statistical Modeling and Analysis of Experiments without ANOVA." ForsChem Research Reports 5 (2018). DOI: 10.13140/RG.2.2.21499.00803

License
-------

This project is licensed under the GNU General Public License. For the exact terms please see the [LICENSE file]
