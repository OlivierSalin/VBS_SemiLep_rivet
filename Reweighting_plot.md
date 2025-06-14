Can have access to the different plot using YODA

Pipeline: YODA --> Root --> Plot 

What we have is the different plot for each weight variation:
In rivet_job.py:
    We can have access to the weight names, indices and we can configure the job option in such a way that the job variation is only on the selected weighted (Reweighted operators)
    Now we would like to compare plots


Remaining questions:
 - Handling of weights --> use YODA or Root inside rivet
 - How to access the weight name inside of Rivet to directly use root Ttree
 - How cross section are handled in rivet/ how to normalise the distribution

