*********************
******************
*
* Bloque3a. Estimacion del AIDS
*
******************
*********************

set more off
use "$LOutput\00Bases\1Enig\0DBaseEstAIDS",clear

rename lnexpfd lnexp

set more off

* Referencia barranquilla
drop dx132 
* Refencia estrato 3
drop dx93

*logaritmo ingreso
replace x10=ln(x10)


global HHControls x1 x2 x3 x4 x5 x6 x7 x8 x10 dx9* dx13*
global PAids lnp1 lnp2 lnp3 lnp4 


********
* regresiones log-log
********

/*
forvalues i=1/4 {
reg lne`i' $PAids $HHControls
}
*/

********
* Ecuacion de seleccion
********

probit di1 $PAids $HHControls,vce(robust)
outreg2 using "$LOutput\30Resu\1Enig\AIDSSelEq.tex",replace  e(N r2_p chi2) ctitle(Coefficient) dec(3) addtext(Controles de estrato, Si, Controles ciudad, Si)  bracket label
* AMEs – Average Marginal Effects
margins, dydx(_all) vce(unconditional) post
outreg2 using "$LOutput\30Resu\1Enig\AIDSSelEq.tex", append  e( N )  ctitle(Margins) dec(3) addtext(Controles de estrato, Si, Controles ciudad, Si) bracket label


set more off
forvalues i=2/4 {
probit di`i' $PAids $HHControls ,vce(robust)
outreg2 using "$LOutput\30Resu\1Enig\AIDSSelEq.tex",append  e(N r2_p chi2) ctitle(Coefficient) dec(3) addtext(Controles de estrato, Si, Controles ciudad, Si)  bracket label
* AMEs – Average Marginal Effects
margins, dydx(_all) vce(unconditional) post
outreg2 using "$LOutput\30Resu\1Enig\AIDSSelEq.tex", append  e( N )  ctitle(Margins) dec(3) addtext(Controles de estrato, Si, Controles ciudad, Si) bracket label
}
set more on

*genera variables auxiliares

set more off
forvalues i=1/4 {
probit di`i' $PAids $HHControls
predict f_d`i', xb
label var f_d`i' "Valor ajustado para el bien `i'"
generate pdf`i'=normalden(f_d`i')
generate cdf`i'=normal(f_d`i')
generate mills`i'=(di`i'*(pdf`i'/cdf`i'))+ ((1-di`i')*(pdf`i'/1-cdf`i'))
}
set more on


********
* Estima QUAIDS censurado
********


drop x1 x2
rename x6 x1
rename x10 x2


*keep if x12==3


cap program drop nlsurquaidsNNP
program nlsurquaidsNNP

version 13

syntax varlist(min=16 max=16) if, at(name)

tokenize `varlist'
args w1 w2 w3 lnp1 lnp2 lnp3 lnp4 lnm x1 x2 pdf1 pdf2 pdf3 cdf1 cdf2 cdf3

// With four goods, there are 15 parameters that can be 
// estimated, after eliminating one of the goods and 
// imposing adding up, symmetry, and homogeneity
// constraints, in the QUAIDS model
// Here, we extract those parameters from the `at'
// vector, and impose constraints as we go along

tempname a1 a2 a3 a4
scalar `a1' = `at'[1,1] 
scalar `a2' = `at'[1,2] 
scalar `a3' = `at'[1,3] 
scalar `a4' = 1 - `a1' - `a2' - `a3'

tempname b1 b2 b3 b4
scalar `b1' = `at'[1,4]
scalar `b2' = `at'[1,5]
scalar `b3' = `at'[1,6]
scalar `b4' = -`b1' - `b2' - `b3'

tempname g11 g12 g13 g14
tempname g21 g22 g23 g24
tempname g31 g32 g33 g34
tempname g41 g42 g43 g44
scalar `g11' = `at'[1,7]
scalar `g12' = `at'[1,8]
scalar `g13' = `at'[1,9]
scalar `g14' = -`g11' - `g12' - `g13'

scalar `g21' = `g12'
scalar `g22' = `at'[1,10]
scalar `g23' = `at'[1,11]
scalar `g24' = -`g21' - `g22' - `g23'

scalar `g31' = `g13'
scalar `g32' = `g23'
scalar `g33' = `at'[1,12]
scalar `g34' = -`g31' - `g32' - `g33'

scalar `g41' = `g14'
scalar `g42' = `g24'
scalar `g43' = `g34'
scalar `g44' = -`g41' - `g42' - `g43'

tempname l1 l2 l3 l4
scalar `l1' = `at'[1,13]
scalar `l2' = `at'[1,14]
scalar `l3' = `at'[1,15]
scalar `l4' = -`l1' - `l2' - `l3'

// constant and household demographics
tempname r11 r12
tempname r21 r22
tempname r31 r32

scalar `r11' = `at'[1,16]
scalar `r12' = `at'[1,17]
scalar `r21' = `at'[1,18]
scalar `r22' = `at'[1,19]
scalar `r31' = `at'[1,20]
scalar `r32' = `at'[1,21]

// pdf
tempname d1 d2 d3
scalar `d1' = `at'[1,22]
scalar `d2' = `at'[1,23]
scalar `d3' = `at'[1,24]

// calculate the expenditure shares. 
quietly {
// First get the price index
// I set a_0 = 6
tempvar lnpindex
gen double `lnpindex' = 10 + `a1'*`lnp1' + `a2'*`lnp2' ///
+ `a3'*`lnp3' + `a4'*`lnp4'
forvalues i = 1/4 {
forvalues j = 1/4 {
replace `lnpindex' = `lnpindex' + ///
0.5*`g`i'`j''*`lnp`i''*`lnp`j''
}
}
// The b(p) term in the QUAIDS model:
tempvar bofp
gen double `bofp' = 0
forvalues i = 1/4 {
replace `bofp' = `bofp' + `lnp`i''*`b`i''
}
replace `bofp' = exp(`bofp')
// Finally, the expenditure shares for 3 of the 4
// goods (the fourth is dropped to avoid singularity)
replace `w1' = (`a1' + `g11'*`lnp1' + `g12'*`lnp2' + ///
`g13'*`lnp3' + `g14'*`lnp4' + ///
`b1'*(`lnm' - `lnpindex') + ///
`l1'/`bofp'*(`lnm' - `lnpindex')^2 + ///
`r11'*`x1' +`r12'*`x2')*`cdf1'+ ///
`d1'*`pdf1'

replace `w2' = (`a2' + `g21'*`lnp1' + `g22'*`lnp2' + ///
`g23'*`lnp3' + `g24'*`lnp4' + ///
`b2'*(`lnm' - `lnpindex') + ///
`l2'/`bofp'*(`lnm' - `lnpindex')^2 + ///
`r21'*`x1' +`r22'*`x2')*`cdf2' + ///
`d2'*`pdf2'

replace `w3' = (`a3' + `g31'*`lnp1' + `g32'*`lnp2' + ///
`g33'*`lnp3' + `g34'*`lnp4' + ///
`b3'*(`lnm' - `lnpindex') + ///
`l3'/`bofp'*(`lnm' - `lnpindex')^2 + ///
`r31'*`x1' +`r32'*`x2')*`cdf3' + ///
`d3'*`pdf3' 
}

end

nlsur quaidsNNP @ w1 w2 w3 lnp1 lnp2 lnp3 lnp4 lnexp ///
x1 x2 pdf1 pdf2 pdf3 cdf1 cdf2 cdf3, ifgnls nequations(3) ///
param(a1 a2 a3 ///
g11 g21 g31 ///
g22 g32 ///
g33 ///
b1 b2 b3 ///
l1 l2 l3 ///
r1 r2 ///
m11 m12 m13 ///
m21 m22 m23 ///
d1 d2 d3 ) nolog 

est store quaidsNNP2


*set trace on
*set tracedepth 4

*Share, price and total expenditure means 
quietly {
foreach x of varlist w* lnp* lnexp {
sum `x'
scalar `x'mean=r(mean)
}
* Price indexes
glo asum "_b[a1]*lnp1mean"
forv i=2(1)4 {
glo asum "${asum} + _b[a`i']*lnp`i'mean"
}
glo gsum ""
forv i=1(1)4 {
forv j=1(1)4 {
glo gsum "${gsum} + 0.5*_b[g`i'`j']*lnp`i'mean*lnp`j'mean"
}
}
glo ap "6 + ${asum} ${gsum}"
glo bp "_b[b1]*lnp1mean"
forv i=2(1)4 {
glo bp "${bp} + _b[b`i']*lnp`i'mean"
}
glo bp "(exp(${bp}))"
* Mus
forv i=1(1)4 {
glo mu`i' "_b[b`i'] + 2*_b[l`i']/${bp}*(lnexpmean-(${ap}))"
}
forv j=1(1)4 {
glo gsum2`j' ""
forv k=1(1)4 {
glo gsum2`j' "${gsum2`j'} + _b[g`j'`k']*lnp`k'mean"
}
}
}
nlcom (a1:_b[/a1])(a2:_b[/a2])(a3:_b[/a3])(a4:1-_b[/a1]-_b[/a2]-_b[/a3]) ///
(b1:_b[/b1])(b2:_b[/b2])(b3:_b[/b3])(b4:-_b[/b1]-_b[/b2]-_b[/b3]) ///
(g11:_b[/g11])(g12:_b[/g21])(g13:_b[/g31]) ///
(g21:_b[/g21])(g22:_b[/g22])(g23:_b[/g32]) ///
(g31:_b[/g31])(g32:_b[/g32])(g33:_b[/g33]) ///
(g14:-_b[/g11]-_b[/g21]-_b[/g31]) ///
(g24:-_b[/g21]-_b[/g22]-_b[/g32]) ///
(g34:-_b[/g31]-_b[/g32]-_b[/g33]) ///
(g41:-_b[/g11]-_b[/g21]-_b[/g31]) ///
(g42:-_b[/g21]-_b[/g22]-_b[/g32]) ///
(g43:-_b[/g31]-_b[/g32]-_b[/g33]) ///
(g44:-(-_b[/g11]-_b[/g21]-_b[/g31])-(-_b[/g21]-_b[/g22]-_b[/g32])-(-_b[/g31]-_b[/g32]-_b[/g33])) ///
(l1:_b[/l1])(l2:_b[/l2])(l3:_b[/l3])(l4:-_b[/l1]-_b[/l2]-_b[/l3]), post noheader

est store quaidsNNP2
outreg2 using "$LOutput\30Resu\1Enig\AIDSParameters.tex",replace dec(3) bracket 

forv i=1(1)4 {
forv j=1(1)4 {
glo delta=cond(`i'==`j',1,0)
glo mu`i'`j' "_b[g`i'`j'] - ${mu`i'}*(_b[a`j'] ${gsum2`j'})-_b[l`i']*_b[b`j']/${bp}*(lnexpmean - (${ap}))^2"
cap nlcom (elasexp`i': ${mu`i'}/w`i'mean + 1) (mu`i'`j': ${mu`i'`j'}), post noheader
if _rc {
est restore quaidsNNP2
qui nlcom (elasexp`i': ${mu`i'}/w`i'mean + 1) (mu`i'`j'f: (1e+2)*(${mu`i'`j'})), post noheader
qui nlcom (elasexp`i': _b[elasexp`i']) (mu`i'`j':_b[mu`i'`j'f]/(1e+2)), post noheader
}
* Uncompensated price elasticity
nlcom (elasexp`i': _b[elasexp`i']) (elu`i'`j':_b[mu`i'`j']/w`i'mean - ${delta}) , post noheader
* Compensated price elasticity
nlcom (elc`i'`j': _b[elu`i'`j'] + _b[elasexp`i']*w`j'mean), noheader
qui est restore quaidsNNP2
}
}

outreg2 using "$LOutput\30Resu\1Enig\Elasticity.tex", coefastr se

set more on
