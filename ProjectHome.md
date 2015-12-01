## Downloads ##
  * **Mac** users please download [SimRare for Mac OS](http://code.google.com/p/simrare/downloads/detail?name=SimRareMac.dmg&can=2&q=) and double-click the downloaded dmg file to open it up. A Finder window will appear. Locate SimRare application in the new Finder window. Double-clicking that will launch SimRare. Drag and drop it into Application directory.
```
Please Note that it only works on Mac OS 10.6+
```


  * **Windows** users please download SimRare for Windows [32-bit](http://code.google.com/p/simrare/downloads/detail?name=SimRareWin32.exe&can=2&q=) or [64-bit](http://code.google.com/p/simrare/downloads/detail?name=SimRareWin64.exe&can=2&q=) (Windows7 tested) and double-click the downloaded exe file to launch the application.

  * **Linux** (Ubuntu & Linux-mint tested and preferred) users please download SimRare for Linux [32-bit](http://code.google.com/p/simrare/downloads/detail?name=SimRareLinux_32bit.tar.gz&can=2&q=) or [64-bit](http://code.google.com/p/simrare/downloads/detail?name=SimRareLinux_64bit.tar.gz&can=2&q=) and open up a new terminal. 'cd' to the directory where SimRareLinux.tar.gz is downloaded. Run commands
```
>  tar -xvzf SimRareLinux.tar.gz
>  cd SimRareLinux
>  ./SimRare
```
> to launch the application. For systems that may miss libssl0.9.8 (or libssl1.0.0) package, if first time launching SimRare returns error about libssl, please install libssl0.9.8 and libssl0.9.8-dbg (or libssl1.0.0 and libssl1.0.0-dbg). In Debian-like box (eg. Ubuntu) Linux system, please run the following command:
> (x.x.x = 0.9.8 or 1.0.0)
```
> sudo apt-get install libsslx.x.x libsslx.x.x-dbg
```
> otherwise, please run
```
> sudo yum install libsslx.x.x libsslx.x.x-dbg
```

  * **Variant data** - Examples of pre-simulated variant data, which include _maf_, _sel_ and _pos_ files, using a variety of demographic models are available to download for the ease of quickly launching all SimRare modules. Refer to [#2.1\_Generation\_of\_variant\_data](#2.1_Generation_of_variant_data.md) for details about variant data simulation and format of output files.
    * [Boyko African model with gene length 1800 and number of replicates 2000](http://www.bioinformatics.org/simped/simrare/Boyko2008African1800Rep2000.zip)
    * [Boyko European model with gene length 1800 and number of replicates 2000](http://www.bioinformatics.org/simped/simrare/Boyko2008European1800Rep2000.zip)
    * [Kryukov European model with gene length 1800 and number of replicates 500](http://www.bioinformatics.org/simped/simrare/Kryukov2009European1800Rep500.zip)

---


## Citation ##
Please cite
Biao Li, Gao Wang and Suzanne M. Leal (2012) [SimRare: a program to generate and analyze sequence-based data for association studies of quantitative and qualitative traits](http://bioinformatics.oxfordjournals.org/content/28/20/2703.long)

## Introduction ##
Currently it is difficult to compare rare variant association methods because there is no standard to generate data and often the comparisons are biased. As a unified simulation platform with user-friendly interface for rare variant studies, SimRare provides an unbiased and easy manner to evaluate association methods, including novel methods. It consists of three modules, variant data simulator, genotype/phenotype generator and association method evaluator. SimRare generates variant data for gene regions using forward-time simulation which incorporates realistic population demographic and evolutionary scenarios. For phenotype data it is capable of generating both case-control and quantitative traits. The phenotypic effects of variants can be detrimental, protective or non-causal. SimRare has a graphical user interface which allows for easy entry of genetic and phenotypic parameters. Simulated data can be written into external files in a standard format. For novel association method implemented in R it can be imported into SimRare, which has been equipped built in functions to evaluate performance of new method and visually compare it with currently available ones in an unbiased manner.

---


  * [#1\_Software\_overview](#1_Software_overview.md)
    * [#1.1\_Features](#1.1_Features.md)
    * [#1.2\_Applications](#1.2_Applications.md)
  * [#2\_Quick\_start\_guide](#2_Quick_start_guide.md)
    * [#2.1\_Generation\_of\_variant\_data](#2.1_Generation_of_variant_data.md)
      * [#2.1.1\_Basics](#2.1.1_Basics.md)
      * [#2.1.2\_Demography](#2.1.2_Demography.md)
      * [#2.1.3\_Genetic\_forces](#2.1.3_Genetic_forces.md)
      * [#2.1.4\_Save\_&\_Run](#2.1.4_Save_&_Run.md)
      * [#2.1.5\_Output\_files](#2.1.5_Output_files.md)
    * [#2.2\_Generation\_of\_genotype\_&\_phenotype\_data](#2.2_Generation_of_genotype_&_phenotype_data.md)
      * [#2.2.1\_Initialization](#2.2.1_Initialization.md)
      * [#2.2.2\_Model\_parameters](#2.2.2_Model_parameters.md)
      * [#2.2.3\_Mimic\_genotyping](#2.2.3_Mimic_genotyping.md)
      * [#2.2.4\_Save\_&\_Run](#2.2.4_Save_&_Run.md)
      * [#2.2.5\_Output\_files](#2.2.5_Output_files.md)
    * [#2.3\_Association\_tests](#2.3_Association_tests.md)
      * [#2.3.1\_Configuration](#2.3.1_Configuration.md)
      * [#2.3.2\_Existing\_methods](#2.3.2_Existing_methods.md)
      * [#2.3.3\_Novel\_method](#2.3.3_Novel_method.md)
      * [#2.3.4\_Run](#2.3.4_Run.md)
      * [#2.3.5\_Output\_graphs\_&\_files](#2.3.5_Output_graphs_&_files.md)
  * [#3\_Downloads](#3_Downloads.md)
  * [#4\_Cookbook\_example](#4_Cookbook_example.md)
    * [#4.1\_Case-control\_sample\_analysis](#4.1_Case-control_sample_analysis.md)
      * [#4.1.1\_Create\_variant\_pool](#4.1.1_Create_variant_pool.md)
      * [#4.1.2\_Choose\_study\_design](#4.1.2_Choose_study_design.md)
      * [#4.1.3\_Evaluate\_association\_methods](#4.1.3_Evaluate_association_methods.md)

---




# 1\_Software\_overview #
SimRare is stand-alone executable software with user-friendly graphical interface implemented in Python/C++ for rare variant association studies. It is designed as a unified simulation framework to provide an unbiased and easy manner to evaluate association methods, including novel methods, under a broad range of choice of biological contexts.

## 1.1\_Features ##
  * **Superior**  -- allow realistic set-up of evolutionary scenarios for simulating sequence-based variant data on gene regions; assign phenotype on generated genotype with a wide range of choice of underlying phenotypic/disease models

  * **Flexible** -- novel rare variant association method written in R can be imported into SimRare for evaluation of type I error, power and comparison with existing methods

  * **Multifunctional** -- integrate genotype generator, phenotype generator and association method evaluator; three modules into one program and each module is individually functional

  * **Standard** -- generated data output in standard _ped_ (Linkage) format; simulation result output in tables or figures comparing different association methods; provide a standard unified platform with everything needed to impartially evaluate novel association method

  * **Intuitive** -- interactive and graphical user-friendly interface; no dependency software installation required; no learning curve


## 1.2\_Applications ##
  * **Benchmark** -- fairly and easily evaluate rare variant association methods to assess their power and robustness

  * **Study design** -- aid parameter analysis of underlying disease model to test a wide variety of scenarios or to achieve adequate power

  * **Method development** -- boost statistical geneticists in design of novel association methods without having to develop special software to generate data in order to evaluate their methods


# 2\_Quick\_start\_guide #
SimRare consists of three modules with tasks of generating sequencing variant data pool, simulating genotype/phenotype data and performing rare variant association analysis, respectively. The program is coupled with graphical interfaces with self-explanatory comment along with each input slot. It should be relatively straightforward to start using. Any wrong input value will be erased immediately and the corresponded input slot will be highlighted by yellow until it receives an allowed one. For more details in how to start running the program and/or specify user inputs for any of the three modules, please see as followed 2.1, 2.2 & 2.3.

## 2.1\_Generation\_of\_variant\_data ##
This module repeatedly evolves a population of gene-based sequences using forward-time simulation with mutation, natural selection and demography. At the end of each repeated simulation, variant sites information such as minor allele frequencies, selection coefficients, position information will be saved. Overall, running this module will create a variant haplotype data pool for a gene region based on any user-defined evolutionary scenario. This simulation module is featured by its capability of incorporating multi-stage population expansion/bottleneck model, handling multi-locus selection model of fitness and random or locus-specific selection coefficient to novel variant, and simplicity of creating gene-based variant data pool based on which individuals' genotype of multiple sample replicates can be easily assigned without performing time-consuming evolutionary simulations round after round again. Eventually running this module, it creates reusable variant data pool based on the user-defined evolutionary scenario.

By launching SimRare program the GUI of this module is activated after triggering **'Need generate rare variant data pool ? --> Click Here'**,  which is shown on the SimRare main window.

A snapshot of the interface is like:

![http://i1050.photobucket.com/albums/s416/libiaospe/varGen3.png](http://i1050.photobucket.com/albums/s416/libiaospe/varGen3.png)

This module is implemented in _SimuPOP_ ([Peng and Kimmel 2005](http://simupop.sourceforge.net)), a general purpose forward-time individual-based population genetics simulation framework and adapts features from _srv_ ([Peng and Liu 2011](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3164177/)), a program implemented in _SimuPOP_ to simulate sequences of the human genome with rare variants.

### 2.1.1\_Basics ###
  * Users may choose to load previously saved configuration file or use parameter defaults to avoid inputting parameters manually. For loading a configuration file, check **'Use Saved Input Parameters?'** and click **'Load File'** to choose the relevant _par_ (short term of parameter) file. For using default input parameters, simply check 'Use Default Input Parameters?'.

```
  Note that: first time user always need to input parameters or choose to use defaults. The par file can be 
  generated by clicking 'Save Input Parameters' and kept for future use. 
```

The manual input of parameters begins with inputting basic information, such as gene length, name of output files and number of replicates.

  * **'Gene Length Mode'** -- check **'Fixed'** to specify a fixed number as the gene length of each replicate; check **'Random'** to choose a random number between a specified range for each replicate as its gene length

```
  Note that: 'Gene Length Mode' need to be selected first in order to activate other input slots.
```

  * **'Gene length'** -- if **'Fixed'** is checked, input the fixed gene length as a positive integer, e.g. **2500**; otherwise, input the range of gene length, e.g. **2500 | 5000**

  * **'Output Files Name (prefix)'** -- a preferred file name as the prefix of output files, e.g. **simVarPool**

  * **'Number of Replicates**' -- input a positive integer, which is equivalent to the size of the data pool, e.g. **100**

### 2.1.2\_Demography ###
For a N-stage demographic model, effective population sizes at the beginning of each stage and at the end of the last stage need to be specified along with specification of number of generations at each stage.

  * **'Effective Population Sizes'** -- need to be a sequence of (N+1) positives integers separated by either space (' ') or comma (,) or semicolon (;), where the first N numbers specify the population sizes at the beginning of each stage and the last number indicates the population size at the end of the last stage. For any two adjacent numbers, _i_ & _j_, if _i = j_; the population size remains constant through the corresponded stage, such as a burn-in stage, etc,; if _i > j_: the population size is reduced to _j_ instantly by entering that stage to mimic bottleneck effect; if _i < j_: the population expands exponentially from _i_ to _j_ through that stage. E.g. **5000, 5000, 800, 30000**, which simulates a three-stage demographic model where a population beginning with 5000 individuals, first undergoes a burn-in stage with constant population size 5000, then goes through a bottleneck of 800 individuals, and after that expands exponentially to a size of 30000.

```
  Note that: if 'Number of Generations Per Stage' has already been specified as a sequence of M numbers, 
  'Effective Population Sizes' need to be a sequence of (M+1) numbers.
```

  * **'Number of Generations Per Stage'** -- this slot requires an input of a series of N positive integers separated by space or comma or semicolon, where the _ith_ value represents the number of generations at the _ith_ stage. E.g.
**300, 100, 200** specifies a three-stage demographic model where the population will be evolved 300, 100, 200 generations through the 1st, 2nd and 3rd stage, respectively.

```
  Note that: if 'Effective Population Sizes' has already been specified as a sequence of (N+1) numbers,
  'Number of Generation Per Stage' need to be a sequence of N numbers.
```

### 2.1.3\_Genetic\_forces ###
Here, users opt to determine mutation model, mutation rate, locus-specific selection model, selection coefficient distribution model and recombination.

  * **'Mutation Model'** -- check **'infinite\_sites'** for selecting infinite-sites mutation model, where mutation can only occur at wild-type loci; check **'finite\_sites'** for finite-sites mutation model, where mutation can occur at any site.

  * **'Mutation Rate'** -- mutation rate per base pair, a float number between (0, 1), e.g. 1.8e-08

  * **'Revert Fixed Sites'** -- check to allow the program to revert fixed sites to wildtype alleles during the evolutionary simulation, otherwise fixed sites will remain as variant sites

  * **'Selection Coefficient Distribution'** -- choose one of most commonly used models as the distribution to assign selection coefficient to any new mutant in order to determine its site-specific fitness of genotypes _AA, Aa/aA, aa_, where _A_ is wildtype allele. By choosing **Constant** ([Williamson\_2005](http://www.pnas.org/content/102/22/7882.full)) it gives all mutant sites a constant selection coefficient as 0.01; whereas any other incorporated choice uses the estimated selection parameters from the corresponded demographic model, [Boyko\_2008\_European](http://www.plosgenetics.org/article/info%3Adoi%2F10.1371%2Fjournal.pgen.1000083), [Boyko\_2008\_African](http://www.plosgenetics.org/article/info%3Adoi%2F10.1371%2Fjournal.pgen.1000083), [Eyre-Walker\_2006](http://www.genetics.org/content/173/2/891.full), [Kyrukov\_2009](http://www.pnas.org/content/106/10/3871.full).

```
  Note that: In single-locus selection model, the fitness of genotypes AA, Aa/aA and aa are given by 
  1, 1-hs and 1-s, respectively, where s is the selection coefficient and h is the dominance coefficient 
  (default 0.5 for additivity). Also, s > 0 results in site being negatively selected, s = 0 in neutral and 
  s < 0 in positively selected. If you would like to define your own selection model, please input your own
  model parameters into 'Customized Selection Coefficient'
```

  * **'Multi-locus Selection Model'** -- choose one from **additive, multiplicative and exponential** to determine the fitness of an individual over all its variant sites

  * **'Customized Selection Coefficient'** -- one can specify any selection coefficient distribution model other than choosing one of those provided by **Selection Coefficient Distribution**. check **'Need Customize Selection Coefficient?'** to enable customized input of selection coefficient distribution model. There are three modeling schemes available (_constant, gamma distributed and mixed-gamma distributed_). For _constant_, it requires a list of two values, **s, h**, constant selection coefficient and dominance coefficient, e.g. **-0.001, 0** specifies a recessive model with fixed positive selection. For the need of _gamma distribution_, a list of three values is required as input, **k, d, h**, where k, d are shape and scale parameters of gamma distribution and h is the dominance coefficient; e.g. inputting **0.206, 0.292, 0.5** specifies a gamma distribution with parameters 0.206, 0.292 for selection coefficient and h=0.5, which mimics [Boyko\_2008\_European](http://www.plosgenetics.org/article/info%3Adoi%2F10.1371%2Fjournal.pgen.1000083) model and has the same effect as choosing [Boyko\_2008\_European](http://www.plosgenetics.org/article/info%3Adoi%2F10.1371%2Fjournal.pgen.1000083) from **Selection Coefficient Distribution**. For specifying _mixed-gamma_ distribution model, an acceptable input is a list of either 5 or 7 parameters which are **p, s, k, d, h** or **p, s, k, d, h, l, u**, where p is the probability of having the selection coefficient equal to s; k, d are gamma distribution coefficients; h is the dominance coefficient and l, u (optional) are lower and upper boundaries of the selection coefficient.

  * **'Recombination Rate'** -- check **'Need Recombination ?'** to enable the input slot for recombination rate per base pair.
```
Note that: if recombination rate times gene length is greater than 0.5, 
a rate of 0.5 divided by gene length will be used. 
```

### 2.1.4\_Save_&_Run ###

  * **'Screen Output Mode'** -- specify screen output mode during running time; choose one from **quiet** (no screen output) , **minimum** (minimum output of simulation progress and time spent for each replicate) or **regular** (regular screen output of statistics, simulation progress and time spent, etc.)

  * **'Detailed Screen Output Interval per Stage'** -- enabled only if **regular** is selected as **Screen Output Mode**; input a single number or a list of n numbers, where n equals to number of demographic stages, as intervals of number of generations at which statistic information will be calculated and output. If left unspecified, screen output from the beginning to the end of every generation of each stage will be shown.

```
Note that:  'output Interval of generations per stage' need to be a single number 
or a list of numbers. 'Numbers of generations per stage' need to be specified first
```

  * **'Save Genotype for Replicates'** -- optional output of genotype information. This option is disabled by default because this format is not efficient in storing rare variants. If enabled and specified by a positive number n, genotype in standard _ped_ format for the first n replicates will  be output to n files. E.g. **1**: genotype of the first replicate will be saved to file; **3**: genotype of the first three replicates will be output to three files. In particular, if specified as **-1**, genotype information of ALL simulation replicates will be saved.

  * **'Save Statistics for Replicates'** -- this optional parameter (disabled by default) may output statistics to _stat_ files. If enabled it should be specified in the same manner as **Save Genotype for Replicates** requires. The statistics are outputted into 7 columns in the order of: 1. generation number, 2. population size (a list), 3. number of segregation sites, 4. average number of segregation sites per individual, 5. average allele frequency x 100, 6. average fitness value, 7. minimal fitness value of the parental population.

  * **'Set Path to Output Files'** -- click **'Path'** to select a desired directory to which output files will be saved. The default path is current directory where SimRare is launched.

  * **'Save Input Parameters'** -- check to save user inputs into an external _par_ file

  * **'Run'** -- click to start generating variant pool.

### 2.1.5\_Output\_files ###

The standard output of this module will be four files with prefix all the same as **Output Files Name (prefix)** specifies and suffix _maf_, _sel_, _pos_ and _len_, respectively. All of them are data files with each one's line number equal to the **Number of Replicates** (assume numReps=n). For _i = 1,2,..., n_; the _ith_ line in _len_ file has a single number which tells the gene/sequence length of the _ith_ replicate; the _ith_ line in _pos_ file contains _n(i)_ values of position information of mutant sites where _n(i)_ is equal to the number of total mutant sites for the _ith_ replicate; the _ith_ line in _maf_ and _sel_ files also has _n(i)_ values which save minor allele frequencies and site-specific selection coefficients on corresponded _n(i)_ mutant sites for the _ith_ replicate.

Additionally, while running any replicate during simulation users have the option of saving genotype information of all individuals into a file in standard _ped_ format with file name as 'fileName\_rep\_i.ped', where 'fileName' is replaced by **Output Files Name (prefix)**. See 2.1.4 for more details.

The _maf_, _sel_ and _pos_ output files are essentially needed to move further to the next module for simulating genotype&phenotype. The output _ped_ files have nothing to do with further analysis and are just to be viewed to check if the program performs correctly in obtaining rare variant sites information.


## 2.2\_Generation\_of\_genotype_&_phenotype\_data ##

Based on the generated genotype from the simulated data pool of sequence-based variants, this module can generate phenotype data for both case-control and quantitative traits. Users have a wide range of choice among available underlying study design models with distinguished parameter settings, such as case-control odds ratio and prevalence ([Risch 1990](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1684987/?tool=pubmed)), population attributable risk ([Madsen and Browning 2009](http://www.plosgenetics.org/article/info:doi/10.1371/journal.pgen.1000384)), quantitative traits ([Kryukov et al. 2009](http://www.pnas.org/content/106/10/3871.full)), QTL extreme traits sampling ([Huang and Lin 2007](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1821103/)) and Mendelian traits ([Ott 1999](http://onlinelibrary.wiley.com/doi/10.1111/j.1469-1809.2000.ahg641_0089_2.x/abstract)) simulation models.

Launching SimRare application pops up the interface for this module immediately.

A snapshot of the interface is like:

![http://i1050.photobucket.com/albums/s416/libiaospe/phenoGen2.png](http://i1050.photobucket.com/albums/s416/libiaospe/phenoGen2.png)

### 2.2.1\_Initialization ###
In order to generate variant sites based genotype and phenotype data one need to first import the variant data pool which is represented by _maf_ (minor allele frequencies), _sel_ (site-specific selection coefficients) and _pos_ (variants positions information) files generated and output by the previous module. Click **load --- file** to load them.
```
Note that: having these files and being able to load them is the prerequisite for the current module to be properly 
functioning. Please click 'Need generate rare variant data pool -->Click Here' if you have never run that before. 
Alternatively, variant sites information estimated from real data or generated by other programs can also be adapted 
to use but have to be firstly converted into 'maf', 'sel' and 'pos' file formats as required, see '2.1.5_Output_files' 
for format details.
```

  * **'Use Saved Parameters'** -- click to open a saved _par_ (input parameters) file and import previously specified inputs to the user interface. See **2.2.4\_Save_&_Run** for how to save a _par_ file that stores user's inputs.

  * **'Choose a study design model'** -- select a study design from the drop-down list on the right of **How do you want to establish genotype-phenotype associations?** Based on the generated genotype and selected study design model phenotype data will be simulated. Currently there are six available study designs to choose from: _Case-control prevalence model of the population; Case-control prevalence model of a sample; Case-control population attributable risk model of a sample; Quantitative-traits (QT) simulation; QTL extreme traits sampling; Mendelian traits simulation_. Once the model has been selected required inputs for model-dependent parameters will be enabled accordingly shown in section **'Simulation Parameters Settings'**

  * **'Use Default Input Parameters'** -- check to use default inputs. The checkbox remains disabled until a study design model is selected.

### 2.2.2\_Model\_parameters ###

  * **List of Parameters**

```
Note that: Here below it lists information on all input slots that are shown on the module interface within 
section 'Simulation Parameters Settings'. But, only those that are enabled are required input parameters 
for the selected study design model. Please refer accordingly. 
```

  * **'Proportion of Detrimental RVs'** -- proportion of functional detrimental/deleterious variants : float number between [0, 1]

  * **'Proportion of Protective RVs'** -- proportion of functional protective/beneficial variants : float number between [0, 1]

```
    Note that: each variant site is determined as synonymous, detrimental or protective by its corresponding 
    selection coefficient
```

  * **'# Cases'** -- number of cases : integer > 0

  * **'# Controls'** -- number of controls : integer > 0

  * **'# Unphenotyped'** -- number of unphenotyped cohort controls : integer > 0

  * **'# Samples'** -- number of total samples : integer > 0

  * **'Fixed Effect Model/Variable Effect Model'** -- choose one : if **Fixed Effect Model** is selected, odds ratio per causal variant is constant; if **Variable Effect Model** is selected, odds ratio per casual variant is determined by its minor allele frequency

  * **'Prevalence'** -- disease prevalence which will be used as baseline penetrance of a gene : float number between (0, .5]

  * **'Odds Ratio for Common Mutations'** -- Odds ratio for common variants : float number (1.0 for neutral, > 1.0 for deleterious and < 1.0 for protective)

  * **'Odds Ratio for Protective Mutations'** -- if **Fixed Effect Model** has been selected, specify a constant (float number between (0, 1]) as odds ratios for all protective variants; if **Variable Effect Model** has been selected, specify a minimum and a maximum (both float numbers between (0, 1]) as the allowed range of odds ratio for protective mutants.

  * **'Odds Ratio for Detrimental Mutations'** -- similar as **Odds Ratio for Protective Mutations**, specify with float numbers >= 1 as odd ratios for detrimental variants.

  * **'Mode of inheritance'** -- choose one mode of inheritance under which the phenotype data is simulated

  * **'Attributable Risk for Detrimental Mutations'** -- total population attributable risk for deleterious variants : float number between [0, 1)

  * **'Attributable Risk for Protective Mutations'** -- total population attributable risk for protective variants : float number between [0, 1)

  * **'Constant Parameters?'** -- check to set locus-specific population attributable risk inversely proportional to its minor allele frequency rather than uniformly distributed

  * **'QT Coefficient for Common Variants'** -- mean value shift for common variants : float number >= 0.0

  * **'QT Coefficient for Causal Variants'**  -- if **Fixed Effect Model** has been selected, specify a constant (float number >= 0.0) as the mean value shift for all variants; if **Variable Effect Model** has been selected, specify a minimum and a maximum (both float numbers >= 0.0) as the allowed range of mean value shift per variant

  * **'QTL Cutoffs'** -- specify a lower and an upper percentile cutoff for quantitative traits in extreme QT sampling : both float numbers between (0, 1)

  * **'Mark Case-Control?'** -- if checked, convert extreme quantitative traits into binary traits by re-coding extreme quantitative traits using the binary coding

  * **'Percentage of Causal RVs'** -- percentage of rare variants being causal in Mendelian traits simulation : float number between (0, 1]

  * **'Proportion of Heterogeneous Cases'** -- proportion of cases that do not carry the disease allele at the gene region : float number between (0, 1]

  * **'Allelic Heterogeneity?'** -- if left unchecked there is no allelic heterogeneity for Mendelian traits, otherwise it will fix the causal variants of the Mendelian trait to the one that has the (proportion of heterogeneous cases)x100%-th smallest minor allele frequency

  * **Study design model specificity**
    * **Case-control prevalence model of the population** -- need to specify proportions of functional RVs (detrimental and protective), # samples, prevalence, effects model, odds ratios (common, detrimental & protective variants) and mode of inheritance.
    * **Case-control prevalence model of a sample** -- need to specify proportions of functional RVs, # cases, # controls, # unphenotyped, effects model, prevalence, odds ratios and mode of inheritance.
    * **Case-control population attributable risk model of a sample** -- need to specify proportions of functional RVs, # cases, # controls, # unphenotyped, mode of inheritance, attributable risks (detrimental and protective mutants) and if locus-specific risk is constant.
    * **Quantitative-traits (QT) simulation** -- need to specify proportions of functional RVs, # samples, effects model, QT coefficients (common and causal variants)
    * **QTL extreme traits sampling** -- need to specify proportions of functional RVs, # unphenotyped, # samples (or # cases and # controls), effects model, QT coefficients, QTL cutoffs, if quantitative traits will be marked as case-control.

```
    Note that: under this study design model, there are two sampling schemes to sample either a number of 
    cases & controls or a number of cohort samples. If # cases and # controls are specified the simulator 
    generates exact numbers of cases and controls as specified that have QT extreme traits, otherwise if 
    # samples is specified, it draws out all individuals that have extreme traits as cases or controls from 
    the cohort samples (those having QT trait values > upper percentile cutoff as cases, < lower percentile
    cutoff as controls), where the number of cases or controls is unknown before the simulation completes
    but definitely smaller than the specified number of cohort samples.
```

  * **Mendelian traits simulation** -- need to specify # cases, # controls, mode of inheritance, % causal RVs, if allelic heterogeneity is envoked and proportion of heterogeneous cases.


### 2.2.3\_Mimic\_genotyping ###
To mimic the actual genotyping process that involves missing sites the estimated proportions of different types of missingness can be specified within section **'Mimic Genotyping Process'**.

Input proportions of **missing detrimental variants, missing protective variants, missing non-causal mutants and missing synonymous mutants**, respectively. Missing genotypes will be encoded as wildtypes by default.

  * Check **'Mark Missing RVs'** to recode missing data from wildtype genotype (**0**) to **-9** to indicate the missingness

### 2.2.4\_Save_&_Run ###

Section **'Save Simulation Results'** and **'Show Data'** are used to determine how to save and show generated genotype/phenotype data.

  * Click **'Set Path for Output File'** to change directory where result and parameter files will be saved. Input a file name (prefix) to specify or change output file name (prefix)

  * To save current interface inputs check the box **'Save Input Parameters?'** and input the desired _par_ (parameter) file name then click **'Save'**. A saved _par_ file can be imported by **'Use Saved Input Parameters'**.

```
  Note that: to launch the module of performing association tests it requires to first import a configuration file, 
  which is the saved _par_ file created here. It is critical to make sure that all related inputs have been correctly 
  specified then click 'Save Input Parameters', name the parameter file and then click 'Save'.
```

  * **'Is Syno Trimmed?'** -- check to have synonymous variant sites removed otherwise they will be kept

  * **'Is CV Trimmed?'** -- check to have common variant sites removed

  * **'Is Ped Written?'** -- check to save the simulated genotype/phenotype dataset to file in _ped_ format.

  * **'Print Genotypes?'** -- check to print out the generated genotype on screen

  * **'Print Phenotypes?'** -- check to print out the generated phenotype on screen

  * **'Print Minor Allele Frequencies?'** -- check to print out variants' MAFs on screen


### 2.2.5\_Output\_files ###

Running this module creates a _ped_ file that stores simulated genotypes and phenotypes if **'Is Ped Written?'** is checked. The optional output of _ped_ file is for viewing purpose only and besides that has nothing to do here with further analysis such as evaluation of association tests, which requires simulating genotype/phenotype dataset repeatedly under the specified modeling scenario and applying association test on each replicate. Therefore, here the _ped_ file can be regarded as data of one replicate.

The configuration/_par_ file that stores information about what model is selected and model related parameters can be written by clicking **Save Input Parameters**. The generated _par_ file containing configurations is the prerequisite for using **module of association tests** (see [2.3\_Association\_tests](http://code.google.com/p/simrare/#2.3_Association_tests) -> [2.3.1\_Configuration](http://code.google.com/p/simrare/#2.3.1_Configuration)).


## 2.3\_Association\_tests ##

To analyze associations between generated genotypes and phenotypes given the underlying demographical and disease models, here it provides a wide range of most commonly used complex trait rare variant association methods to choose from and to compare each other's performance, such as 1.) Combined Multivariate and Collapsing Method (CMC) ([Li and Leal 2008](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2842185/?tool=pubmed)); 2.) Weighted Sum Statistic (WSS) ([Madsen and Browning 2009](http://www.plosgenetics.org/article/info:doi/10.1371/journal.pgen.1000384)); 3.) Variable Threshold (VT) ([Price et al. 2010](http://www.cell.com/AJHG/abstract/S0002-9297(10)00207-7)); etc. It is also capable of incorporating any novel association method that is developed in R program with existing methods into the comparative evaluation on type I and type II errors. Users need to specify testing thresholds, such as significance level and number of replicates to be generated. A large number of existing tests are permutation based because the significance of their association methods cannot be evaluated analytically. To apply those tests users also need to specify number of permutations and should expect a large load of computational burden.

To activate this module, click **'Need apply existing/novel association methods? --> Click Here'** shown on the upper right corner of SimRare main interface.
A snapshot is like:

![http://i1050.photobucket.com/albums/s416/libiaospe/assocTest2.png](http://i1050.photobucket.com/albums/s416/libiaospe/assocTest2.png)
### 2.3.1\_Configuration ###

Click **'Load configuration File'** to choose the desired parameter (_.par_) file that has been created before by SimRare's main interface (module of genotype & phenotype simulation, see [2.2.5\_Output\_files](http://code.google.com/p/simrare/#2.2.5_Output_files)).
```
Note that: loading configuration file successfully and correctly is the first and crucial step.  
```

  * **'Rare variant frequency bound'** -- set lower and upper bounds of observed sample minor allele frequency : float numbers between [0, 1), loci having observed MAF < lower bound or > upper bound will not be analyzed

  * **'Specify number of replicates for power calculation'** -- integer > 0

  * **'Specify Significance level'** -- at which power will be evaluated : float number (0, 0.5]

  * **'Specify number of permutations'** -- only applicable to permutation based methods : integer > 0

### 2.3.2\_Existing\_methods ###

Select existing association methods that are going to be evaluated.

  * **'CMC'** -- ([Li and Leal 2008](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2842185/?tool=pubmed)), Combined Multivariate and Collapsing, choose **'CMC-one'** for one-sided test, choose **'CMC-QT'** for test of quantitative traits

  * **'MZ'** -- ([Morris and Zeggini 2010](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2962811/)), searching for RV accumulations within the same functional unit, choose **'MZ-one'** for one-sided test, choose **'MZ-QT'** for test of quantitative traits

  * **'WSS'** -- ([Madsen and Browning 2009](http://www.plosgenetics.org/article/info:doi/10.1371/journal.pgen.1000384)), Weighted Sum Statistic, choose **'WSS-one'** for one-sided test

  * **'VT'** -- ([Price et al. 2010](http://www.cell.com/AJHG/abstract/S0002-9297(10)00207-7)), Variable Threshold, choose **'VT-one'** for one-sided test

  * **'KBAC'** -- ([Liu and Leal 2010](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2962811/)), Kernel-Based Adaptive Cluster, choose **'KBAC-one'** for one-sided test

```
  Note that: 'CMC-QT', 'MZ-one', 'MZ', 'KBAC-one', 'KBAC', 'VT-one', 'VT' are permutation based tests which cost
  a lot of computational time
```

### 2.3.3\_Novel\_method ###

Novel association method implemented in R program can be imported and evaluated with other existing methods. Check **'Load R file of a new method'** and then click **'Load File'** to choose a _r_ file.

In order to successfully apply loaded R script on the simulated genotype/phenotype data, above all, the user need to make sure that R is in system PATH,
```
open a terminal (or command prompt in Windows) and try if command 'Rscript' works.
If it is not recognized as a command, set R to system PATH first. 
```
and then follow the steps below:
  1. Include a few lines of code on top of the R script to read **pedFile**, which denotes the full path to the file that saves the simulated data, and to retrieve phenotypes and genotypes from it. For example, add these lines to the very beginning of your R script
```
  # read ped file as a data matrix
  simData <- as.matrix(read.table(pedFile, header=FALSE, sep=" "))
  # first 5 columns are irrelevant to our need, phenotypes are stored in the 6th col 
  # genotypes are stored from the 7th to the last column 
  pheno = simData[, 6]
  geno = simData[, -(1:6)]
```
```
  Note that: "pedFile" is a fixed name not a name variable. Do **NOT** alter it to anything else. 
```
  1. Establish associations between phenotypes and genotypes using the novel testing method and return the calculated p-value. Note that genotypes are unphased and allelic status of each locus is encoded by 0, 1, 2 mode, where 0 for wildtype homozygote, 1 for heterozygote and 2 for variant homozygote
  1. At the bottom of the R script, attribute the calculated p-value to a R variable that **MUST** be named by **pValue** in order for the p-value to be passed into plot and output functions. E.g.
```
  pValue <- R_function_to_calculate_p_value(geno, pheno, otherArgs)
```

### 2.3.4\_Run ###

Click **'Run'** to simultaneously begin generating genotype/phenotype replicates based on the configuration and evaluating type I and type II errors for selected existing methods and/or the novel method. If running error is encountered double-check if every step is strictly followed as [2.3.3\_Novel\_method](http://code.google.com/p/simrare/#2.3.3_Novel_method) requires while creating the R script.

### 2.3.5\_Output\_graphs_&_files ###

Click **'Set Path for Output File'** to set path for output files, otherwise the current working directory will be used.

Once the simulation finishes it outputs two graphical files to show comparison of association methods, one draws QQ plots for type I errors and the other is a bar graph to display power comparison. It also stores type I errors and powers into two additional text files, respectively.

```
Note that: in order to draw plots R must be in the system PATH
```

# 3\_Downloads #

[click here to download instructions](http://code.google.com/p/simrare/#Downloads)

[click here to download files](http://code.google.com/p/simrare/downloads/list)


# 4\_Cookbook\_example #

## 4.1\_Case-control\_sample\_analysis ##
This example shows how to use SimRare to evaluate association methods on simulated case-control samples. It generates the variant/haplotype pool using the demographic model based on [Kryukov et al. 2009](http://www.pnas.org/content/106/10/3871.full) and each individual's phenotype based on its genotype and disease-related parameters given by the case-control prevalence model of a sample.

### 4.1.1\_Create\_variant\_pool ###
We will create a variant pool with 50 replicates (50 independent realizations of gene sequence evolution given the same demographic parameters as input).
  * Launch SimRare and activate the module by clicking **Need generate rare variant data pool ? --> Click Here** on the main interface.
  * -> Select **Use Default Input Parameters?** on the upper right corner of the window (Default inputs are based on Kryukov's demographic model)
  * -> Change **Number of Replicates** into **50**
  * -> Click **Path** button and select a folder where the output files will be saved
  * -> Click **Run**
A snapshot looks like:

![http://i1050.photobucket.com/albums/s416/libiaospe/fig4_1.jpg](http://i1050.photobucket.com/albums/s416/libiaospe/fig4_1.jpg)

After it finishes running you will find files, such as simVarPool.maf, simVarPool.sel and simVarPool.pos, inside the selected **Path** folder. Those 3 files are the created variant/haplotype pool for 50 replicates. For meaning of each parameter shown on the interface please refer to [#2.1\_Generation\_of\_variant\_data](#2.1_Generation_of_variant_data.md)

### 4.1.2\_Choose\_study\_design ###
In this cookbook example, we focus on the case-control prevalence model of a sample and establish genotype-phenotype simulation for both under the null (H0) and under the alternative (H1) hypothesis
```
Note: under H0, it assumes that there is no relationship between rare variants and the complex trait. 
```

  * On SimRare main module (genotype & phenotype simulation) choose **Case-control prevalence model of a sample** from the drop-down list on the right side of **How do you want to establish genotype-phenotype associations?**
  * -> Select  **Use Default Input Parameters?** below
  * -> Click **Load .maf file** in the **Initialization** area and select the simVarPool.maf file generated in 4.1.1 (SimRare will automatically load **.sel and**.pos files given that they are saved under the same folder by the same prefix)
  * -> In the **Simulation Parameter Settings** area change both **# Cases** and **# Controls** to be **1000**
  * -> Change **Odds Ratio for Detrimental Mutations** to be **1.6**
  * -> Change **Odds Ratio for Protective Mutations** to be **0.9**
  * -> In the **Save Simulation Results** area click **Set Path for Output File** and choose a folder to save output files (can be the same folder selected in 4.1.1)
  * -> Select **Save Input Parameters?**, specify a file name, e.g. **MySimu\_H1.par** to save current user inputs, hit return key and click **Save** button next to it

A snapshot will be looking like:
![http://i1050.photobucket.com/albums/s416/libiaospe/fig4_2_new.png](http://i1050.photobucket.com/albums/s416/libiaospe/fig4_2_new.png)


```
  Note: Here, (1) we do not need to run or save any generated dataset into external files since our ultimate goal 
           is to evaluate and compare association methods in 4.1.3. Thus, it is sufficient to just save those inputs 
           as model parameters into an external *.par (configuration) file for later use but without the need to 
           click *Run*. For more details please refer to section 2.2.5
           (2) given those numbers above specified for odds ratios it will generate dataset under the alternative 
           hypothesis, therefore, we name the *.par file as *_H1.par 
```
  * -> Now repeat those steps above but put **1.0** to both **Odds Ratio for Detrimental Mutations** and **Odds Ratio for Protective Mutations**, specify a different parameter file name, such as **MySimu\_H0.par**, hit return key and click **Save** button
  * -> in the case of odds ratio for any type of variant being equal to 1 it will generate data under H0

A snapshot of the interface here will be looking like:
![http://i1050.photobucket.com/albums/s416/libiaospe/fig4_3.png](http://i1050.photobucket.com/albums/s416/libiaospe/fig4_3.png)

### 4.1.3\_Evaluate\_association\_methods ###
Association tests will be performed using parameters settings under H0 and H1, respectively. You may choose to evaluate any one or more test methods. Here we will select CMC-one, WS-one, KBAC-one and MZ-one.
  * Activate the module of association tests
  * -> Click **Load Configuration file** and select the **?H0.par** file generated by 4.1.2
  * -> Select existing association methods **CMC-one, MZ-one, KBAC-one and WS-one**
  * -> By default, **number of replicates for power calculation** is **2000**, you may change it to be smaller, e.g. 1000, if you do not want to generate too many replicates
  * -> Click **Set path to output file** to specify any folder where output figures will be saved (it can be the same folder selected for in 4.1.1 and 4.1.2)
  * -> Specify an output file name (prefix), such as **test\_H0**
  * -> Click **Run**

A snapshot of the interface here will be looking like:

![http://i1050.photobucket.com/albums/s416/libiaospe/fig4_4.png](http://i1050.photobucket.com/albums/s416/libiaospe/fig4_4.png)

```
  Note: running this step requires relatively long time particularly if number of replicates and/or number of
           permutations is large, please be patient...
           Also, you may choose to save datasets of generated genotype/phenotype for all replicates into 
           external files (*_rep1.ped, *_rep2.ped, ..., *._repx.ped) by clicking "Output geno/pheno info for all
           replicates" at the lower left corner of the interface. But this is NOT recommended if your purpose 
           is just to evaluate association methods because it will slow down the running speed while creating 
           a large number of files and outputting to your local folder.
```

  * -> Now repeat those steps above but load the other .par file (**?H1.par**) generated by 4.1.2
  * -> Specify a different file name, such as **test\_H1** and click **Run**

A snapshot of the interface here will be looking like:

![http://i1050.photobucket.com/albums/s416/libiaospe/fig4_5.png](http://i1050.photobucket.com/albums/s416/libiaospe/fig4_5.png)

Find plots of power comparison for both cases (under **H0** and **H1**) from the output file folder.
There you can open test\_H0\_powerPlot.pdf,
![http://i1050.photobucket.com/albums/s416/libiaospe/test_H0_powerPlot.png](http://i1050.photobucket.com/albums/s416/libiaospe/test_H0_powerPlot.png)

and test\_H1\_powerPlot.pdf,
![http://i1050.photobucket.com/albums/s416/libiaospe/test_H1_powerPlot.png](http://i1050.photobucket.com/albums/s416/libiaospe/test_H1_powerPlot.png)

It is reasonable to see that if data is generated under H0, the test power for any association method is around the significance level, which is 0.05 in this example. Meanwhile, if data is generated under H1 based on the odds ratio specified above, all methods have approximately equal power to detect associations, while CMC and WS are slightly better than the others in this particular scenario of study design.

Also, you may check p-values by viewing QQ plot for all association methods under H0. For example, CMC-one,
you can open test\_H0\_QQPlot\_CMC-one.pdf
![http://i1050.photobucket.com/albums/s416/libiaospe/test_H0_QQPlot_CMC-one.png](http://i1050.photobucket.com/albums/s416/libiaospe/test_H0_QQPlot_CMC-one.png)




