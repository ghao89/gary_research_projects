# gary_research_projects

Folder "DLSMT" contains the code, data and results for the book chapter "Multiple testing for detecting differentially methylation in sequencing data".

* __Requires "fdrDiscreteNull" package__
* __Workflow__
    1. Use "simulation_setup.R" to simulate data for discrete multiple testing. The result will be stored as "dat.Rdata". There are a few parameters that you could tune based on your interest.
        1. 'm'. The total number of tests. Default is 5000.
        2. 'n'. The number of samples in both groups. Default is 2.
        3. 'pi_0'. The proportion of true null hypotheses. Default is 0.6.
        4. 'mu_pois'. The mean of the Poisson distribution used to generate the number of trials. Default is 20.
    2. Use the demo files to try different discrete FDR control procedures included in the chapter separately.
        1. "fdrDiscreteNull_demo.R" contains BH procedure, BHH method, aBH method, aBHH method and Habiger's method.
        2. "gilbert_demo.R" contains Gilbert's method.
        3. "mcf_demo.R" contains MCF method.
        4. "storey_demo.R" contains Storey's q-value method.
        5. "liang_demo.R" contains Liang's method.
    3. Use the "replicate_book_chapter_results.Rmd" to replicate the plots in the chapter.