*! Version 2021.02.23
********************************************************************************
* Author: Carlos Eduardo Veras Neves								           *
*                                                                              *
* Estimation of Technical Efficiency with SFA and DEA for Brazilian            *
* Federal Toll Road Concessions 
*                                                                              *
********************************************************************************
* configuração inicial
set trace off

clear all

set more off

capture log close

// muda diretório de trabalho
cd "...\Stata"

// carrega a base de dados original
//use "dados.2012_2016_adj.dta"
import delimited "data\dados_2012_2016_adj.csv"


// create output folder with current time and date included in the name
local D = c(current_date)
local T = c(current_time)
local T = subinstr("`T'",":","_",.)

*local fname = "SFA3.`D'`T'"
local fname = "SFA"

// cria diretório para armazenar arquivos de saída (FAZER CÓDIGO PARA VERIFICAR EXISTÊNCIA DO DIRETÓRIO)
mata : st_numscalar("OK", direxists("`fname'"))

if scalar(OK) == 0{
	mkdir "`fname'"  
}


// entra no diretorio dos arquivos de simulação
cd "`fname'"

// log file
log using "log`fname'.`D'`T'.smcl", replace

/************************************************************************
(3) ANÁLISE SFA, TOMANDO COMO VARIÁVEIS:
- INPUTS: REC; CUST;
- OUTPUTS: AVG 

OUTPUT ORIENTED
*************************************************************************/
* set the directory to the location of .ado files 
* .ado files from:
*
* A Practitioner's Guide to Stochastic Frontier Analysis Using Stata         
* Subal C. Kumbhakar, Hung-Jen Wang and Alan P. Horncastle                     
*                                                                               
//adopath ++ "C:\ado\personal\sfbook\sfbook_ado"

// apaga qualquer unidade que tenha valor de receita ou custo nulo 
//drop if rec  <=0 
//drop if cust <=0

/*
subsititui os id's originais por id's para dados em painel (um id para cada concessionária)
*/
// lista de concessionárias na base de dados 
levelsof conc, local(concs)
// contador local para gerar id
local count 1
// loop na lista de concessionárias
foreach v in `concs'{
	
	// substitui o id original, para os casos em que há match entre o nome da
	// concessionária na base e na lista de nomes 
	quietly replace id = `count' if conc == "`v'"
	local ++count // atualiza o contador de id
}

// gera variáveis dummy para marcar as etapas
quietly gen dEtapa2 = (etp==2)
quietly gen dEtapa3 = (etp==3)

// Como a análise na parte (1) indicou não haver variação relevante nas variáveis ao longo do tempo, optou-se por calcular as médias do período das variáveis de input e output
quietly bysort id: egen custm = mean(cust_km_adj)
quietly bysort id: egen recm = mean(rec_km_adj)
quietly bysort id: egen avgm = mean(avg)
quietly bysort id: egen tarm = mean(tar_adj)
// identifica duplicatas
quietly by id: gen dup = cond(_N==1, 0, _n)
// apaga duplicatas
drop if dup > 1

// gera índices para as DMUs
quietly gen dmu = id

// apaga variáveis dispensáveis para a análise DEA
//quietly drop id ano anoconc cust cust_km rec rec_km avg dup

// salva dados para simulação DEA
save data_sfa, replace

/*
gera logaritmo das variáveis para a função de produção Cobb-Douglas

*/
gen l_rec = log(recm)
gen l_cust = log(custm)
gen l_avg = log(avgm)

// define inputs e outputs da análise SFA
global inp l_cust l_rec
global out l_avg

* variáveis explicativas
global reglist npp tarm dEtapa2 dEtapa3 

********************************************************************************************
* Testa se o modelo especificado é válido antes de proceder com as estimativas
* Foi utilizado o teste do terceiro momento do resíduo OLS (M3T) de Coelli (1995)
* para checar se os resíduos são assimétricos para a esquerda, o que confirmaria
* a possibilidade de ser utilizada a especificação para a função de produção. 
* Para uma função de produção de fronteira estocástica, o erro composto(v - u), u > 0 e v
* distribuído simetricamente em torno de zero, deve ter os resíduos da estimatica OLS
* enviesadas para a esquerda. 
********************************************************************************************
//reg $out $inp
reg $out l_cust l_rec npp
predict e, residual
sum e, detail
sktest e, noadj

qnorm e
kdensity e

* notar que a distribuição é assimétrica para a esquerda (Skewness < 0)
quietly sum e
local e_mean = r(mean)
local N = r(N)
egen double m2 = mean((e - `e_mean')^2)
egen double m3 = mean((e - `e_mean')^3)
gen double M3T = m3/sqrt(6*((m2)^3)/`N')
//di "Como ", M3T[1], " é negativo, então os resíduos são assimétricos para a esquerda (rejeição de H0: não há assimetria). "
di  M3T[1]

*********************************************************************************************
// (1) OLS - Fronteira determinística
reg $out $inp 
eststo OLS

esttab OLS using "sfaolsest.tex" , beta tex replace

predict e_ols, residual /* salva o resíduo da estimativa OLS na variável e */
quietly summarize e_ols /* obtém as estatísticas sem imprimir nada na tela */
gen double eff_ols = exp(e_ols)
summarize eff_ols

*********************************************************************************************
// (2) COLS - Fronteira determinística (OLS Corrigido)
reg $out $inp 
eststo COLS

esttab COLS using "sfacolsest.tex" , beta tex replace


predict e_cols, residual /* salva o resíduo da estimativa OLS na variável e */
quietly summarize e_cols /* obtém as estatísticas sem imprimir nada na tela */
gen double u_star = -(e_cols - r(max)) /* obtém ineficiência a partir dos resíduos corrigidos */ 
gen double eff_cols = exp(-u_star) /* calcula eficiência técnica */
summarize eff_cols
list conc u_star eff_cols

*********************************************************************************************
* (3) CMAD - Corrected Mean absolute Deviation 
qreg $out $inp /* regressão quantílica - utiliza a mediana */
eststo CMAD

esttab CMAD using "sfacmadest.tex" , beta tex replace

estout CMAD, style(tex)

predict e_cmad, residual /* salva o resíduo da estimativa OLS na variável e */
quietly summarize e_cmad /* obtém as estatísticas sem imprimir nada na tela */
gen double eta_star_q = -(e_cmad - r(max)) /* obtém ineficiência a partir dos resíduos corrigidos */ 
gen double eff_cmad = exp(-eta_star_q) /* calcula eficiência técnica */


*********************************************************************************************
* (4) Cobb-Douglas Half-Normal Model
frontier $out $inp 
eststo sfa_hn

esttab sfa_hn using "sfa_hnest.tex" , beta tex replace

predict e_hn, u /* salva o resíduo da estimativa OLS na variável e */
quietly summarize e_hn /* obtém as estatísticas sem imprimir nada na tela */
gen double eta_hn= -(e_hn - r(max)) /* obtém ineficiência a partir dos resíduos corrigidos */ 
gen double eff_hn = exp(-eta_hn) /* calcula eficiência técnica */



*********************************************************************************************
* (5) Cobb-Douglas Exponential Model
frontier $out $inp, distribution(exponential) 
eststo sfa_exp

esttab sfa_exp using "sfa_expest.tex" , beta tex replace

predict e_exp, u /* salva o resíduo da estimativa OLS na variável e */
quietly summarize e_exp /* obtém as estatísticas sem imprimir nada na tela */
gen double eta_exp= -(e_exp - r(max)) /* obtém ineficiência a partir dos resíduos corrigidos */ 
gen double eff_exp = exp(-eta_exp) /* calcula eficiência técnica */

/* não converge!!!!
*********************************************************************************************
* (6) Cobb-Douglas Trucated Normal
frontier $out $inp, distribution(tnormal) cm($reglist)
eststo sfa_exp

predict e_exp, u /* salva o resíduo da estimativa OLS na variável e */
quietly summarize e_exp /* obtém as estatísticas sem imprimir nada na tela */
gen double eta_exp= -(e_exp - r(max)) /* obtém ineficiência a partir dos resíduos corrigidos */ 
gen double eff_exp = exp(-eta_exp) /* calcula eficiência técnica */

*/

*********************************************************************************************
* (7) Cobb-Douglas Half-Normal Model - com variáveis explicativas (uhet)
*
* Estimador da matrix de variância-covariância (vce) do tipo "clustered sandwich estimator",
* assumindo que as observações entre os grupos são inpendentes, porém dentro dos grupos
* pode não haver independência. 
frontier $out $inp, uhet($reglist) vce(cluster etp)  
eststo sfa_hn_ex

esttab sfa_hn_ex using "sfa_hn_exest.tex" , beta tex replace

predict e_hnx, u /* salva o resíduo da estimativa OLS na variável e */
quietly summarize e_hnx /* obtém as estatísticas sem imprimir nada na tela */
gen double eta_hnx= -(e_hnx - r(max)) /* obtém ineficiência a partir dos resíduos corrigidos */ 
gen double eff_hnx = exp(-eta_hnx) /* calcula eficiência técnica */

*********************************************************************************************
* (7) Cobb-Douglas Half-Normal Model - com variáveis explicativas (uhet)
*
* Estimador da matrix de variância-covariância (vce) do tipo "clustered sandwich estimator",
* assumindo que as observações entre os grupos são inpendentes, porém dentro dos grupos
* pode não haver independência.  (MESMA ESPECIFICAÇÃO DO MODELO EXP_EX)

frontier $out $inp, uhet(npp tarm) vce(cluster etp)  
eststo sfa_hn_ex2

esttab sfa_hn_ex2 using "sfa_hn_ex2est.tex" , beta tex replace

predict e_hnx2, u /* salva o resíduo da estimativa OLS na variável e */
quietly summarize e_hnx2 /* obtém as estatísticas sem imprimir nada na tela */
gen double eta_hnx2= -(e_hnx2 - r(max)) /* obtém ineficiência a partir dos resíduos corrigidos */ 
gen double eff_hnx2 = exp(-eta_hnx2) /* calcula eficiência técnica */

*
*********************************************************************************************
* (8) Cobb-Douglas Exponential Model- com variáveis explicativas (uhet)
*
*  Como o modelo aceita somente uma variável explicativa, foi escolhida a extensão, na
* medida que estão sendo utilizadas as receitas e custos sem ser por km
frontier $out $inp, distribution(exponential) uhet(npp tarm ) vce(cluster etp) 
eststo sfa_exp_ex

esttab sfa_exp_ex using "sfa_exp_exest.tex" , beta tex replace


predict e_expx, u /* salva o resíduo da estimativa OLS na variável e */
quietly summarize e_expx /* obtém as estatísticas sem imprimir nada na tela */
gen double eta_expx= -(e_expx - r(max)) /* obtém ineficiência a partir dos resíduos corrigidos */ 
gen double eff_expx = exp(-eta_expx) /* calcula eficiência técnica */


// saída de dados:
summarize eff_ols eff_cols eff_cmad eff_hn eff_exp eff_hnx eff_hnx2 eff_expx 

esttab using "summary_eff_sfa.tex" , stat(mean sd count) tex replace

list dmu conc eff_ols eff_cols eff_cmad eff_hn eff_hnx eff_hnx2 eff_exp  eff_expx 
estimates table OLS COLS CMAD sfa_hn sfa_hn_ex sfa_hn_ex2 sfa_exp sfa_exp_ex, b  star (.01 .05 .1)   //stats(N r2 r2_o r2_b r2_w F chi2)drop(_cons)

esttab using "summary_eff_sfa_all.tex" , tex replace

save results_sfa, replace

savesome dmu conc etp eff_ols eff_cols eff_cmad eff_hn eff_hnx eff_hnx2 eff_exp eff_expx using eff_sfa, replace

est clear  // clear the est locals
summarize eff_ols eff_cols eff_cmad eff_hn eff_hnx eff_hnx2 eff_exp eff_expx


est clear  // clear the est locals
estpost tabstat eff_ols eff_cols eff_cmad eff_hn eff_hnx eff_hnx2 eff_exp eff_expx, c(stat) stat(mean sd min max n)
esttab using "summary_eff_sfa.tex" , ///
 cells("mean(fmt(%13.2fc)) sd(fmt(%13.2fc)) min(fmt(%13.2fc)) max(fmt(%13.2fc)) count(fmt(%13.0fc))") nonumber ///
  nomtitle nonote noobs label collabels("Mean" "SD" "Min" "Max" "N") tex replace

est clear  // clear the est locals
by etp, sort: summarize eff_ols eff_cols eff_cmad eff_hn eff_hnx eff_hnx2 eff_exp eff_expx

esttab using "summary_eff_sfa_byetp.tex" , tex replace

/*est clear  // clear the est locals
estpost tabstat eff_ols eff_cols eff_cmad eff_hn eff_hnx eff_hnx2 eff_exp eff_expx, c(stat) stat(mean sd min max n)
esttab using "summary_eff_sfa.tex" , ///
 cells("mean(fmt(%13.2fc)) sd(fmt(%13.2fc)) min(fmt(%13.2fc)) max(fmt(%13.2fc)) count(fmt(%13.0fc))") nonumber ///
  nomtitle nonote noobs label collabels("Mean" "SD" "Min" "Max" "N") tex replace*/




export delimited conc  tar etp npp ext dEtapa2 dEtapa3 custm recm avgm tarm dmu eff_ols eff_cols eff_cmad eff_hn eff_hnx eff_hnx2 eff_exp eff_expx using "sfa_results.csv",replace

* close log file
log close
exit
>>>>> end here

