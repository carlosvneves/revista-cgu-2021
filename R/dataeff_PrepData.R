#######################################################
# Corrige receitas e custos das concessionárias de
# acordo com a inflação - índice : IPCA
#######################################################
setwd("D:/Google Drive/_DOUTORADO/__TESE/__CODEandDATA/__TESE_code-data/data")

library(sidrar)
library(tidyverse)
library(tsibble)

# Acessa IPCA número índice para deflacionar séries diretamente da base SIDRA-IBGE
ipca <- get_sidra(api='/t/1737/p/all/v/2266/N1/all')%>%
  mutate(date = as.yearmon(`Mês (Código)`, format='%Y%m')) %>%
  dplyr::select(date, Valor) %>%
  as_tibble()

ipca <- ts(ipca$Valor, start=c(1979,12), frequency = 12)
ipca <- window(ipca, start=2012, end=2017)
ipca <- as_tsibble(ipca)

# filtra dados até dezembro de 2016
ipca <- ipca %>% filter_index("2012-01" ~ "2016-12")

# calcula os fatores de correção 
ipca <- ipca %>% mutate(fator = tail(ipca$value,1)/ipca$value) 

# seleciona os valores dos índices para os anos da base de dados
ipca_yr <- ipca %>% filter_index("2012-01","2013-01","2014-01","2015-01","2016-01" )
ipca_yr <- tibble(ipca_yr)

# Lê a base de dados
dados_2012_2016 <- read.csv("dados_2012_2016.csv")

# transforma base de dados em tsibble para manipulação
data_tbl <- tibble(dados_2012_2016)

# acrescenta índices de correção na tabela
data_tbl <- data_tbl %>% mutate(index = yearmonth(paste0(ano,"-01")))

# constrói a tabela com os fatores de inflação
data_tbl <- merge(data_tbl, ipca_yr,all=TRUE)

# calcula os valores ajustados tendo como base janeiro de 2016
data_tbl <- data_tbl %>% mutate(coef_adj = 1 + (fator - tail(fator,1))/tail(fator,1)) %>% 
                      mutate(cust_adj = cust * coef_adj) %>%
                      mutate(rec_adj = rec * coef_adj) %>% 
                      mutate(cust_km_adj = cust_adj / ext) %>%
                      mutate(rec_km_adj = rec_adj / ext)  %>%
                      mutate(tar_adj = tar * coef_adj)

# ordena pelo id da DMU de modo ascendente
data_tbl <- data_tbl %>% arrange(id)  


# escreve dados em arquivo para a análise de eficiência
write.csv(data_tbl, "dados_2012_2016_adj.csv")

########################################
# TODO: executar a análise por ano e a análise com os valores médios
# e comparar as conclusões
########################################
library(plm)

datapanel <- data_tbl

# gera dummies para o ano e para a etapa
datapanel <- datapanel %>% mutate(d2012 = if_else(ano == 2012, 1, 0)) %>%
                           mutate(d2013 = if_else(ano == 2013, 1, 0)) %>%  
                           mutate(d2014 = if_else(ano == 2014, 1, 0)) %>%
                           mutate(d2015 = if_else(ano == 2015, 1, 0)) %>%
                           mutate(d2016 = if_else(ano == 2016, 1, 0)) %>%
                           mutate(dEtapa2 = if_else(etp == 2, 1, 0))  %>%
                           mutate(dEtapa3 = if_else(etp == 3, 1, 0))   


# gera o id por dmu
datapanel$dmu <- datapanel %>% group_by(conc) %>% group_indices(conc)

# gera o dataframe para análise de dados em painel
datapanel <- pdata.frame(datapanel, index =("dmu"))

# dados para testes com as receitas e custos não parametrizados por km
datatest <- datapanel[,c("dmu","ano","rec_adj","cust_adj",
                         "avg","ext","tar","mat","d2012","d2013","d2014","d2015","d2016")]

modelo1 = plm(log(avg) ~ rec_adj + cust_adj + ext + tar + d2012 + d2013 + d2014 + d2015 + d2016, data = datatest, model = "pooling")

summary(modelo1)

plmtest(modelo1, effect='twoways' , type='ghm' )

#Estamos, por suposto, rejeitando a hipótese nula de ausência de efeitos

modelo2 = plm(log(avg) ~  rec_adj + cust_adj + ext + tar, data = datatest, model = "within")

summary(modelo2)

pFtest(modelo2, modelo1)

modelo3 = plm(log(avg) ~ rec_adj + cust_adj + ext + tar - 1, 
              data = datatest, effect = "twoways", model = "random")

summary(modelo3)

phtest(modelo2, modelo3)



# Os modelos sugerem que a variação entre indivíduos é relevante, enquanto que a variação no tempo é irrelevante



#####################################################################################
# dados para testes com as receitas e custos parametrizados por km
datatest <- datapanel[,c("dmu","ano","rec_km_adj","cust_km_adj","avg","d2012","d2013","d2014","d2015","d2016")]

modelo1 = plm(log(avg) ~ log(rec_km_adj) + log(cust_km_adj) + d2012 + d2013 + d2014 + d2015 + d2016, data = datatest, model = "pooling")

summary(modelo1)

plmtest(modelo1, effect='twoways' , type='ghm' )

#Estamos, por suposto, rejeitando a hipótese nula de ausência de efeitos

modelo2 = plm(log(avg) ~  log(rec_km_adj) + log(cust_km_adj), data = datatest, model = "within")

summary(modelo2)

pFtest(modelo2, modelo1)

modelo3 = plm(log(avg) ~ log(rec_km_adj) + log(cust_km_adj), data = datatest,  effect = "twoways", model = "random")

summary(modelo3)

phtest(modelo2, modelo3)



modelo4 <- plm(avg ~ rec_km_adj + cust_km_adj, 
                      data = datatest,
                      index = c("dmu","ano"),
                      model = "within", 
                      effect = "twoways")

lmtest::coeftest(modelo4, vcov = vcovHC, type = "HC1")


#######################################################################
# Como a análise dos dados em formato de painel não demonstram
# que há variação considerável dos dados. Logo, trabalharemos com valores
# médios e corte em seção transversal

data_med <- datapanel %>% group_by(conc) %>% summarize(avg = mean(avg),
                                           rec = mean(rec_km_adj), 
                                           cust = mean(cust_km_adj), 
                                           tar = mean(tar_adj), 
                                           ext = mean(ext),
                                           npp = mean(npp),
                                           etp = mean(etp))


data_med <- data_med %>% mutate(dEtp2 = if_else(etp == 2,1,0)) %>%
                         mutate(dEtp3 = if_else(etp == 3,1,0))  

# escreve dados em arquivo para a análise de eficiência
write.csv(data_med, "data_eff.csv")

library(Benchmarking)

y <- log(data_med$avg)
rec <- log(data_med$rec)
cust <- log(data_med$cust)
x <- cbind(rec,cust)
o <- Benchmarking::sfa(x,y)


e <- dea(x,y, RTS="vrs", ORIENTATION="out",)
1/eff(e)


dea.plot.frontier(x[,1],y,txt=TRUE)

sb <- dea.add(x,y,RTS="vrs")


sol_MM <- stoned(x, y)
sol_PSL <- stoned(x, y, METHOD="PSL")

plot(x[,1],y)
points(x[,1],sol_MM$front, col="red")
points(x[,1],sol_PSL$front, col="blue", pch=16, cex=.6)


plot(x[,2],y)
points(x[,2],sol_MM$front, col="red")
points(x[,2],sol_PSL$front, col="blue", pch=16, cex=.6)

y <- as.matrix(y)

ec <- dea(x,y, RTS="crs", ORIENTATION = "out")
Ec <- 1./eff(ec)
ev <- dea(x,y, RTS="vrs", ORIENTATION = "out")
Ev <- 1./eff(ev)
# The test statistic; equation (6.1)
S <- sum(Ec)/sum(Ev)

# To calculate CRS and VRS efficiencies in the same bootstrap replicas
# we reset the random number generator before each call of the
# function dea.boot.

# To get the an initial value for the random number generating process
# we save its state (seed)
save.seed <- sample.int(1e9,1)

# The bootstrap and calculate CRS and VRS under the assumption that
# the true technology is CRS (the null hypothesis) and such that the
# results correponds to the case where CRS and VRS are calculated for
# the same reference set of firms; to make this happen we set the
# random number generator to the same state before the calls.
set.seed(save.seed)

nrep <- 10000

bc <- dea.boot(x,y, nrep,RTS="crs", ORIENTATION = "out")
set.seed(save.seed)
bv <- dea.boot(x,y, nrep,RTS="vrs",ORIENTATION = "out", XREF=x,YREF=y, EREF=ec$eff)

# Calculate the statistic for each bootstrap replica
bs <- colSums(bc$boot)/colSums(bv$boot)
# The critical value for the test (default size \code{alpha} of test is 5%)
critValue(bs, alpha=.1)
S
# Accept the hypothesis at 10% level?
critValue(bs, alpha=.1) <= S

# The probability of observering a smaller value of S when the
# hypothesis is true; the p--value.
typeIerror(S, bs)
# Accept the hypothesis at size level 10%?
typeIerror(S, bs) >= .10



sol_MM <- stoned(x, y, RTS = "vrs", MULT=1, METHOD="MM")
sol_PSL <- stoned(x, y, RTS = "vrs", MULT=1, METHOD="PSL")
