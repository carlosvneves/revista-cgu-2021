** versão de fevereiro de 2021, com valores deflacionados pelo IPCA (jan base 2016)

// limpa dados da seção do Stata
clear all

// muda diretório de trabalho
cd "...\data"

// carrega a base de dados
import delimited "dados_2012_2016_adj.csv"

//exclui variáveis desnecessárias
drop v1 x

** /////////////////////////////////////////////////////////////////////////////
/************************************************************************
(1) AVALIAR SE AS VARIÁVEIS ESCOLHIDAS SÃO RELEVANTES PARA A ANÁLISE DE 
DE DESEMPENHO DAS CONCESSIONÁRIAS E TAMBÉM SE HÁ UMA VARIAÇÃO RELEVANTE
DO DESEMPENHO EM TERMOS DE NOTA DA PESQUISA CNT NO TEMPO.
*************************************************************************/

/*
subsititui os id's originais por id's para dados em painel (um id para cada concessionária)
*/
// lista de concessionárias na base de dados 
levelsof conc, local(concs)
// contador local para gerar id
/*local count 1
// loop na lista de concessionárias
foreach v in `concs'{
	
	// substitui o id original, para os casos em que há match entre o nome da
	// concessionária na base e na lista de nomes 
	replace id = `count' if conc == "`v'"
	local ++count // atualiza o contador de id
}*/


// configura painel
xtset id ano

preserve
// Gráficos para explorar os dados

// gráfico de variação within de avg
xtdata, fe clear
graph twoway scatter avg ano || lfit avg ano
graph save Graph "Graph_within_var_fev21.gph", replace
restore

// gráfico de variação between de avg
preserve
xtdata, be clear
graph twoway scatter avg ano || lfit avg ano
graph save Graph "Graph_between_var_fev21.gph", replace
restore

// Agrupados por etapa
//Custo/km no tempo
graph twoway scatter cust_km_adj ano || lfit cust_km_adj ano, by(etp)
graph save Graph "Graph_cust_km_byetp_fev21.gph", replace
//Receita/km no tempo
graph twoway scatter rec_km_adj ano || lfit rec_km_adj ano, by(etp)
graph save Graph "Graph_rec_km_byetp_fev21.gph", replace
//Avaliação Geral no tempo
graph twoway scatter avg ano || lfit avg ano, by(etp)
graph save Graph "Graph_avg_byetp_fev21.gph", replace

// Por Concessionária
//Custo/km no tempo
graph twoway scatter cust_km_adj ano || lfit cust_km_adj ano, by(conc)
graph save Graph "Graph_cust_km_byconc_fev21.gph", replace
//Receita/km no tempo
graph twoway scatter rec_km_adj ano || lfit rec_km_adj ano, by(conc)
graph save Graph "Graph_rec_km_byconc_fev21.gph", replace
//Avaliação Geral no tempo
graph twoway scatter avg ano || lfit avg ano, by(conc)
graph save Graph "Graph_avg_km_byconc_fev21.gph", replace

// gera variáveis dummy para marcar as etapas
gen dEtapa2 = (etp==2)
gen dEtapa3 = (etp==3)

global x cust_km_adj rec_km_adj 
// modelo POLS com Erros-Padrão Robustos Clusterizados
reg avg $x dEtapa2 dEtapa3, vce (cluster id)
eststo	 POLS_rob

// modelo com estimador between
xtreg avg $x dEtapa2 dEtapa3, be
eststo BE

// modelo com estimador efeitos fixos
xtreg avg $x, fe
scalar Fstat = e(F)
di "Estatística F:", Fstat // Teste de Chow: indica que a H0 de que os interceptos são todos iguais é rejeitada (POLS não mais adequado que FE)
eststo FE

// modelo com estimador efeitos fixos com Erros-Padrão Robustos Clusterizados
xtreg avg $x, fe vce (cluster id)
eststo FE_rob

// modelo com estimador efeitos aleatórios
//xtreg avg $x dEtapa2 dEtapa3, re
xtreg avg $x, re
xttest0 // Teste Breush- Pagan:indica que a H0 de que o modelo POLS é mais adequado que o RE é rejeitada
eststo RE

// modelo com estimador efeitos aleatórios com Erros-Padrão Robustos Clusterizados
quietly xtreg avg $x dEtapa2 dEtapa3, re vce (cluster id)
eststo RE_rob

//estimates table POLS_rob BE FE FE_rob RE RE_rob, b stats(N r2 r2_o r2_b r2_w F chi2) star (.05 .10 .15)
estimates table POLS_rob BE FE FE_rob , b stats(N r2 r2_o r2_b r2_w F chi2) star (.05 .10 .15)

// Obtendo a decomposição de variância para cada variável
xtsum id ano avg cust_km_adj rec_km_adj tar npp ext  

// verificando se modelo mais adequado é o de FE ou de RE
hausman FE RE, sigmamore

*** TESTE ROBUSTO DE HAUSMAN
quietly xtreg avg $x, re
sort id ano
by id: gen T= _N
gen theta = 1-sqrt(e(sigma_e)^2/(e(sigma_e)^2+T*e(sigma_u)^2))
foreach var of varlist avg cust_km rec_km tar{

by id: egen mean`var' = mean(`var')
gen `var'_re = `var' - theta * mean`var'
gen `var'_fe = `var' - mean`var'

}

qui reg avg_re cust_km_re rec_km_re tar_re cust_km_fe rec_km_fe tar_fe, vce(cluster id)
test cust_km_fe rec_km_fe tar_fe 


// usando RE que é mais consistente segundo o teste de Hausman:
preserve
statsby, by(etp) clear: xtreg avg $x, re
list, clean
restore

// salva resultados
save output_phase1, replace

// limpa espaço de trabalho
clear all

/* Conclusões: 
- os gráficos demonstram que no tempo, a avaliação é praticamente constante, sendo mais relevante a diferença entre concessionárias;
- é possível realizar a análise correlacionando as variáveis avg, cust_km, rec_km (existe
relação relevante entre elas);
- resultado xtsum indica que a variação entre concessionárias é mais relevante que no tempo (a avaliação de eficiência considerando a média dos dados do período parece adequada);
- Teste de Chow: indica que a H0 de que os interceptos são todos iguais é rejeitada (POLS não mais adequado que FE);
- Teste Breush- Pagan:indica que a H0 de que o modelo POLS é mais adequado que o RE é rejeitada;
- Teste de Hausman: não é possível rejeitar a H0 de que modelo RE oferece estimativas dos parâmetros mais consistentes (0.05 de significância).
- Coeficientes para cada etapa de concessão:

. statsby, by(etp) clear: xtreg avg $x, re
(running xtreg on estimation sample)

      command:  xtreg avg cust_km rec_km tar, re
           by:  etp

Statsby groups
----+--- 1 ---+--- 2 ---+--- 3 ---+--- 4 ---+--- 5 
...

. list, clean

       etp   _b_cus~m   _b_rec_km      _b_tar    _b_cons  
  1.     1   .0001197   -.0000441    -.009915   3.857822  
  2.     2   .0004137   -1.06e-06    .2419112   2.806615  
  3.     3   .0002248   -.0000439   -.0806115    4.04382 

  - O aumento da receita para todas as etapas influencia negativamente a avaliação da rodovia;
  - O aumento dos custos totais (incluindo investimentos) influencia positivamente a avaliação da rodovia;
  - O aumento da tarifa influencia de modo distinto 1a e 3a etapas, e 2a etapa. Para esta, o aumento da tarifa influencia positivamente a avaliação da rodovia.

**/



