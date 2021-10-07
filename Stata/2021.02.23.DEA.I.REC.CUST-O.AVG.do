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
//use "dados.2012_2016.dta"
. import delimited "data\dados_2012_2016_adj.csv"

// create output folder with current time and date included in the name
local D = c(current_date)
local T = c(current_time)
local T = subinstr("`T'",":","_",.)

*local fname = "DEA3.`D'`T'"

local fname = "DEA"

// cria diretório para armazenar arquivos de saída
mata : st_numscalar("OK", direxists("`fname'"))

if scalar(OK) == 0{
	mkdir "`fname'"  
}  

// entra no diretorio dos arquivos de simulação
cd "`fname'"

// log file
log using "log`fname'.`D'`T'.smcl", replace

** /////////////////////////////////////////////////////////////////////////////
/************************************************************************
(2) ANÁLISE DEA, TOMANDO COMO VARIÁVEIS:
- INPUTS: REC; CUST;
- OUTPUTS:  AVG;


*************************************** Comentário anterior ************************************************************
* A Escolha das variáveis leva em conta a análise realizada na parte (1), em que se 
verifica a relação positiva entre Custos Totais por km e a avaliação geral do trecho, quando se avalia o efeito de cada variável por etapa no modelo RE. Por outro lado, o aumento da receita não implica em melhora na avaliação geral. Isto posto, o modelo parte da premissa de que, por meio da abordagem DEA orientada a input, é possível identificar aquelas unidades em há melhor equilíbrio relativo entre receitas, custos e avaliação geral. 
A ideia seria classificar como ineficiente, aquela unidade em que as receitas estão em um nível muito elevado em relação aos custos totais e a avaliação geral do trecho.

Vale a observação de que Profeta (PROFETA, 2014: REGULAÇÃO E EFICIÊNCIA DOS MODELOS DE CONCESSÕES DE RODOVIAS NO BRASIL, TESE DE DOUTORADO - UFV/MG) entendeu que a eficiência estaria relacionada simplesmente à minimização dos custos totais. Assim elegeu como INPUT custos totais, e como OUTPUS Receitas Totais e Avaliação Geral. Contudo, aparentemente, Profeta não se atentou ao fato de que, de acordo com os relatórios financeiros das concessionárias, os custos totais incluem investimentos, que naturalmente implicam em melhora da condição da rodovia (e por consequência na sua avaliação. Assim, pode ser que não faça sentido a eficiência estar baseada na minimização desses custos.
Os resultados da parte (1) parecem indicar também esse fato e, por isso, optou-se por uma configuração diferente de INPUTs e OUTPUTs.
*************************************** Comentário atualizado **********************************************************

OUTPUT ORIENTED (para compatibilizar com SFA)


*************************************************************************/

// apaga qualquer unidade que tenha valor de receita ou custo nulo 
drop if rec  <=0 
drop if cust <=0

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
quietly bysort id: gen custm = log(cust_km_adj)
quietly bysort id: gen recm = log(rec_km_adj)
quietly bysort id: gen avgm = log(avg)
quietly bysort id: gen tarm = tar_adj
// identifica duplicatas
quietly by id: gen dup = cond(_N==1, 0, _n)
// apaga duplicatas
drop if dup > 1

// gera índices para as DMUs
quietly gen dmu = id

// apaga variáveis dispensáveis para a análise DEA
//quietly drop id ano anoconc cust cust_km rec rec_km avg dup

// salva dados para simulação DEA
save data_dea, replace

// define inputs e outputs da análise DEA (colocar em log???????)
global inp recm custm
global out avgm 


/***********************************************************************
ANÁLISE DEA CONSIDERANDO AS ORIENTAÇÕES INPUT E OUTPUT
************************************************************************/
foreach direction in i o{

	* DEA - variable returns to scale
	quietly: dea $inp = $out, rts(vrs) ort(`direction') stage(2) saving(dea_results)
	mat deascores = r(dearslt)
	mat deascores = deascores[1...,"theta"]
	sort dmu
	cap drop dea1
	svmat deascores, names(dea)
	rename dea1 deascore_`direction'
	gen efficient_`direction' = deascore_`direction' == 1

}

* Escores de eficiência após primeiro estágio
tabstat deascore_i deascore_o efficient_i efficient_o, columns(statistics) statistics(mean sd median min max) format(%7.0g)


di " ##### ANÁLISE DEA - INPUTS: REC;CUST - OUTPUTS: AVG ##### "

list dmu conc deascore_i deascore_o
/************************************************************************
SECOND-STAGE REGRESSION ANALYSIS - TOBIT

*************************************************************************/
* variáveis explicativas
global reglist npp tar dEtapa2 dEtapa3

tabstat $reglist, columns(statistics) statistics(mean sd median min max) format(%7.0g)

* (Naive) Censored Regression Analysis - eficiência orientada a input
/*
quietly: stepwise, pr(.20): tobit deascore_i $reglist, ul(1) /* backward selection */
eststo tob_inp
*/

* (Naive) Censored Regression Analysis - eficiência orientada a output
//quietly: stepwise, pr(.20): tobit deascore_o $reglist, ul(1) /* backward selection */
//eststo tob_out

/************************************************************************
SECOND-STAGE REGRESSION ANALYSIS - Procedimento Simar & Wilson (2007)

*************************************************************************/
* Simar & Wilson (2007) Procedure - orientado a input
// simarwilson deascore_i $reglist, unit dots reps(1999) 
// eststo sw_inp
/*logscore makes simarwilson use the natural logarithm of the efficiency score as the left-hand-side
        variable in the truncated regressions.  With logscore specified, truncation is at 0 rather than
        at 1 and is always one sided.  If externally estimated scores are used, do not take the
        logarithm beforehand, but let the original score enter the estimation procedure as depvar.*/


* Simar & Wilson (2007) Procedure - orientado a output
//simarwilson deascore_o $reglist, logscore nounit dots reps(5000)
//eststo sw_out

/*estimates table tob_inp sw_inp tob_out sw_out, b stats(N r2 r2_o r2_b r2_w F chi2) star (.05 .1 .2)*/

//estimates table tob_out sw_out, b stats(N r2 r2_o r2_b r2_w F chi2) star (.05 .1 .2)
 

********************************************************************************
* adotando outro código, o algoritmo fica bem mais simples:
********************************************************************************

teradialbc $out = $inp, rts(vrs) base(o) heterogeneous tename(score) tebc(tobscore) biasboot(te1bias) varboot(te1var) biassqvar(te1biassqvar) telower(te1lower) teupper(te1upper) reps(2000)

* (Naive) Censored Regression Analysis - eficiência orientada a output
tobit tobscore $reglist,  ll(1) /* backward selection */
eststo tob_out1

* Algoritmo #2 de Simar-Wilson (2007)
//simarwilson scorebc $reglist, logscore nounit reps(5000) dots
//eststo sw_out1

//predict ohat

//Simar and Wilson (2007) two-stage efficiency analysis with internal bias correction, that is,
// algorithm 2
simarwilson ( $out = $inp ) $reglist, logscore nounit algorithm(2) reps(2000) bcreps(1500) dots rts(vrs) base(output) saveall(BBSTR) tebc(swscore)
eststo sw_out2

predict ohat_2


* valores de eficiência entre 0 e 1
gen deascore_bc = 1./tobscore
replace deascore_bc = 1 if deascore_bc > 1

//gen deascore_sw = 1-exp(ohat) //1./ohat
//replace deascore_sw = 1 if deascore_sw > 1

gen deascore_sw = 1./swscore
//replace deascore_sw = 1 if deascore_sw_2 > 1


tabstat deascore_o deascore_bc deascore_sw, columns(statistics) statistics(mean sd median min max) format(%7.0g)

//est clear  // clear the est locals
estpost tabstat deascore_o deascore_bc deascore_sw, c(stat) stat(mean sd min max n)
esttab using "summary_eff_dea.tex" , ///
 cells("mean(fmt(%13.2fc)) sd(fmt(%13.2fc)) min(fmt(%13.2fc)) max(fmt(%13.2fc)) count(fmt(%13.0fc))") nonumber ///
  nomtitle nonote noobs label collabels("Mean" "SD" "Min" "Max" "N") tex replace


//estimates table tob_out1 sw_out2, b stats(N r2 r2_o r2_b r2_w F chi2) star (.05 .1 .2)

//est clear  // clear the est locals
estimates table tob_out1 sw_out2, b  star (.01 .05 .1)   //stats(N r2 r2_o r2_b r2_w F chi2)drop(_cons)
esttab using "regression_dea2s.tex" , tex replace



/*
#############################################################################
Conclusão da análise:

- Em nenhum dos modelos se pode concluir que a 5% de significância pelo menos 
um dos coeficientes não é igual a zero. 
- Para o nível de significância de 5%, em princípio, não é necessário realizar
correções por meio das variáveis de ambiente.

*/



//by etp, sort: summarize npp tar

by etp, sort: summarize deascore_o deascore_bc deascore_sw


savesome dmu conc etp deascore_o deascore_bc deascore_sw using eff_dea, replace

export delimited conc tar etp npp ext dEtapa2 dEtapa3 custm recm avgm tarm dmu deascore_o deascore_bc deascore_sw  using "dea_results.csv", replace

* close log file
log close
exit
>>>>> end here	

