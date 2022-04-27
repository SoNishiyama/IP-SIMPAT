# IP-SIMPAT
IP-SIMPAT (***I***nterger ***P***rogramming for ***S***electing ***I***nformative ***M***arkers in ***PAT***ernity inference)
is designed for optimizing marker set for paternity inference in half-sib family.
This repo includes a MATLAB implementation and analytical codes used in [Nishiyama et al. 2021](https://github.com/SoNishiyama/IP-MARS), 

This implementation has three prerequisites:  
1. the maternal individual is known 
2. the potential paternal individuals have been identified (i.e., a closed population) 
3. the potential parents have been genotyped  

Under these assumptions, IP-SIMPAT selects small number of markers that could distinguish the paternal parent crossed wtih a particular maternal parent from an offspring genotype.  

### Overview
`IP-SIMPAT-main/` containts all the MATLAB implementation for the proposed optimization method. `util/` contains several scripts 
that would be used for input preparation and genotype simulation. `test-data/` contains test data that is adoption of valuable apple genotype data by [Muranty et al. 2020](https://doi.org/10.1186/s12870-019-2171-6) 
with some process presented in [Nishiyama et al. 2021](https://github.com/SoNishiyama/IP-MARS), and could be used for optimization test.

## Citation
XXX XXX

## Step-by-step tutorial
Here we are going to find an optimized marker set to identify a paternal parent of offspring whose maternal parent is an apple cultivar 'Fuji'. The test data `Fuji-small` includes genotype data 
of 200 biallelic SNPs for 29 candidates paternal individuals plus 'Fuji'. For the detail, please refer to [Nishiyama et al. 2021](https://github.com/SoNishiyama/IP-MARS).

### input preparation (optional)
Here we process `Fuji-small.vcf`, containing initial genotype data for genome-wide markers, for subsequent marker selection step using a custom python script.
```
$ python SIMPAT-input-prep-vcf.py Fuji-small.vcf Fuji-small 318
```
`318` is an individual code representing target maternal parent, 'Fuji', used in Muranty et al (2020). This needs to be changed according to the representation of your data.
The above command would produce two outputs: `Fuji-small.num.txt`, including genotype information, and `Fuji-small.pos.txt`, including marker position.
  
IP-SIMPAT requires additional input, reprenting chromosome length. An example of this input would be found at `test-data/apple-contig-length`.

### optimization
All the optimization codes would be found at `IP-SIMPAT-main/`.  
We proposed two implementations, namely `greedy+ns` and `intlinprog_quad` in the manuscript, and here provide MATLAB live script for each. 

|methods|description|function|live script|
|:--|:--|:--|:--|
|intlinprog_quad|MATLAB optimization toolbox function intlinprog, slow but maybe more accurate, takes very long time for ≥ 10,000 initial markers |pat_intlinprog_quad.m|intlinprog_quad_example.mlx|
|greedy+ns|implementation of proposed greedy algorithm and subsequent neighborhood search, fast and applicable for large-scale data |pat_greedy.m, pat_nsearch.m|greedy_ns_example.mlx|

```
parameters:

h, heterozygosity weight          A weight controls assingment power of homozygous (AA or BB) - heterozygous AB pairs; h homozygous - heterozygous pairs have equal discrimination power as the homozygous AA – homozygous BB pairs. The oprimal h depends on genetic relatedness of parent population but h = 8-16 provides good optimization in most cases.  
w, adjacency weight control       A test weight that controls negative power on closely located marker pairs. q_kl = (|B_k-B_l |/D_k )^(-1/w). Large w reduces the negative power and thus physically close pairs can be more likely selected. My test yielded that w = 4, which add 10 times weight for a 3 kb distance pair located in a 30 Mb chromosome than a linkage equibrium pair, can be a good option to produce reasonable optimization, and the smaller produces conservative solution. 
constthre                         constraint threshold, how many distinct loci (homoz ref vs homoz alt, plus info from heterogyzous loci depending on h) is needed for discriminating a pair of individual. Setting "1" worked in every optimization in our test. The larger produces conservative solution, while more markers are required.
flip_frac, flip fraction v        fraction v of markers to be flipped in neighborhood search
```

The test function `intlinprog_single` without applying the adjacency weight **Q**, proposed in [Nishiyama et al. 2021](https://github.com/SoNishiyama/IP-MARS),
was also included in `IP-SIMPAT-main/`.

### input preparation for Cervus
Next we perform parentage analysis using a likelihood-based software [Cervus](http://www.fieldgenetics.com/pages/aboutCervus_Overview.jsp). A script for formatting into Cervus is provided in `util/`.

```
$ python cervus_input_prep.py Fuji-small.num.txt Fuji-small.ind optimized_set_h4 fuji.h4
```

### production of simulated offspring
To see the effectiveness of the optimized marker set, it might be interesting to perform simulated crossing between the target maternal parent and each of potential paternal parent, and validate if the set successfully produce the parentage assignment.
Here a script to perform simulated crossing is provided in `util/`

```
$ python mating_fixed_mother.py Fuji-small.num.txt 318 Fuji-small.pat.ind optimized_set_h4 5 fuji.h4.sim
```
The above command will produce simulated genotype of five individuals for each of crosses between Fuji (318) and each of 29 candidate paternal individuals.
Three outputs would be produced: `.cervus.gen.txt` and `.offspring.txt` can be formatted and directly applied for Cervus. `.truecross.txt` contains true combination of simulated crosses, and this would be used for the subsequent test.

### run Cervus
Here we use Cervus, one of the most common software for parentage assignment. Please refer its [documentation](http://www.fieldgenetics.com/pages/aboutCervus_Overview.jsp) upon running Cervus.

### Accuracy test
A script that automatically tests the inference accuracy was provided in `util/`, and we would use this for evaluating optimized marker set.
```
$ python cervus_res_validation.py input.list.txt output.tsv
```

## references
- Muranty, H, Denancé, C, Feugey, L et al. (2020) Using whole-genome SNP data to reconstruct a large multi-generation pedigree in apple germplasm. BMC Plant Biol 20: 2. https://doi.org/10.1186/s12870-019-2171-6
- Kalinowski, ST, Taper, ML & Marshall, TC (2007) Revising how the computer program CERVUS accommodates genotyping error increases success in paternity assignment. Molecular Ecology 16: 1099-1106. http://dx.doi.org/10.1111/j.1365-294x.2007.03089.x

