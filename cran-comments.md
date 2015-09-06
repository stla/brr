## Release summary

This is my first submission.

## Test environments

* ubuntu 14.04, R 3.2.2
* windows 7 64bit, R 3.0.2
* online win-builder.r-project.org

## R CMD check results

There were no ERRORs or WARNINGs. 

There were 3 NOTEs.

* checking CRAN incoming feasibility ... NOTE
Maintainer: ‘Stéphane Laurent <laurent_step@yahoo.fr>’
New submission

* checking top-level files ... NOTE
File README.md cannot be checked without ‘pandoc’ being installed.
  
* checking R code for possible problems ... NOTE

brr_generic: no visible binding for global variable ‘.’
coef.brr: no visible binding for global variable ‘.’
confint.brr: no visible binding for global variable ‘.’
plot.brr: no visible binding for global variable ‘.’

plot.brr: no visible binding for global variable ‘S’
plot.brr: no visible binding for global variable ‘y’

prior: no visible binding for global variable ‘a’
prior: no visible binding for global variable ‘b’
prior: no visible binding for global variable ‘d’

  The dot ‘.’ is related to the pipe operator %>% of the magrittr package.
  The messages about the plot.brr() anf prior() functions are due to the use of assign() in the code. I think it is OK.





