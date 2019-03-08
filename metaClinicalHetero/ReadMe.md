# ReadMe

This is R code to run a meta-analysis method with adjusted maximum likelihood method for the between-study variance (Published: ***). 

-Tittle: 
Guarantee clinical heterogeneity in random-effect meta-analysis under small number of studies: an adjusted maximum likelihood approach for solving a between-study boundary estimate problem

-Authors:
Daisuke Yoneoka (DY) and Masayuki Henmi (MH)
DY are mainly responsible for writing the code.

-Emails:
DY: blue.sky.sea.dy@gmail.com


-Session info:
>> sessionInfo()
>R version 3.3.2 (2016-10-31)
>Platform: x86_64-apple-darwin13.4.0 (64-bit)
>Running under: macOS Sierra 10.12.3
>
>locale:
>[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
>
>attached base packages:
>[1] stats     graphics  grDevices utils     datasets  methods   base    


-How to run:
The main function is re.reml() in case_study.R, and thus please modify the code if one wants to implement our method in their own dataset.
In addition, we will publish R package called “metaClinicalHetero” from CRAN.

For replication of the application section, please run case_study.R.
For replication of the simulation section, please run simulation.R.

 




