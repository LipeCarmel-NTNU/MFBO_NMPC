# Parecer técnico — *Runtime-Aware Multi-Fidelity Bayesian Optimisation for NMPC Tuning*

## 1. Escopo do parecer

Este documento apresenta uma avaliação técnica, metodológica e editorial do manuscrito:

> *Runtime-Aware Multi-Fidelity Bayesian Optimisation for NMPC Tuning*

Arquivo avaliado: `C:\Users\patri\Downloads\main (1).pdf`

Data da avaliação: 28 de julho de 2026.

O manuscrito foi lido integralmente e inspecionado visualmente em suas 14 páginas. A avaliação foi conduzida sob quatro perspectivas:

1. contribuição, novidade e posicionamento científico;
2. formulação de NMPC e validade das afirmações de controle;
3. otimização Bayesiana multiobjetivo, multifidelidade e modelagem de custo;
4. evidência experimental, reprodutibilidade, narrativa e prontidão editorial.

Também foram realizadas duas leituras técnicas independentes: uma com foco em NMPC/controle de processos e outra com foco em MFBO/MOBO. Os pareceres independentes convergiram nos principais bloqueios relatados abaixo.

### Convenção epistemológica

- **OBSERVADO:** diretamente presente no manuscrito.
- **REPORTADO:** declarado no manuscrito, mas não verificado em dados ou código.
- **INFERIDO:** interpretação técnica plausível, ainda não demonstrada causalmente.
- **VERIFY:** ponto que depende de inspeção dos artefatos primários.
- **AÇÃO MÍNIMA:** correção necessária para sustentar uma versão aplicada do artigo.
- **AÇÃO FORTE:** análise ou reformulação capaz de elevar a contribuição metodológica.

## 2. Veredito executivo

Eu gosto do artigo. A ideia é relevante, a aplicação é interessante e existe uma contribuição publicável. O elemento mais forte é a integração de:

- tuning biobjetivo de NMPC;
- horizonte de simulação como fidelidade contínua;
- extrapolação de custos acumulados;
- modelagem explícita do runtime;
- aquisição Pareto sensível ao custo;
- reavaliação dos candidatos em fidelidade completa;
- seleção de um controlador para implantação a partir da frente de Pareto.

No estado atual, contudo, a recomendação é **major revision forte**.

Se o trabalho for posicionado como uma nova metodologia MFBO ou como uma demonstração comprovada de eficiência, robustez e validade externa, a recomendação mais provável seria **reject and resubmit**. Se for reposicionado como um workflow aplicado e exploratório de MOBO com avaliações fechadas truncadas, existe um artigo potencialmente bom após as correções.

| Aspecto | Estado atual | Potencial após revisão |
|---|---:|---:|
| Qualidade da ideia | 4,5/5 | 4,5/5 |
| Novidade como integração aplicada | 3,5/5 | 4/5 |
| Novidade algorítmica | 2,5/5 | 3,5/5 |
| Rigor experimental | 2/5 | 4/5 |
| Validação de controle | 2,5/5 | 4/5 |
| Reprodutibilidade | 2/5 | 4/5 |
| Escrita e apresentação | 4/5 | 4,5/5 |

## 3. Contribuição científica reconstruída

### 3.1 Questão de engenharia

Como ajustar parâmetros, pesos e horizontes de um NMPC quando cada avaliação exige uma simulação fechada custosa contendo a solução repetida de problemas de programação não linear?

### 3.2 Gargalo

O tuning combina:

- espaço misto discreto–contínuo;
- objetivos conflitantes de tracking e atividade de controle;
- custo de avaliação dependente do próprio controlador;
- possibilidade de avaliações longas, ruidosas e numericamente difíceis.

### 3.3 Método proposto

O manuscrito combina:

- qLogNEHVI para otimização multiobjetivo;
- horizonte truncado parametrizado por uma fidelidade contínua \(z\);
- extrapoladores polinomiais para estimar custos de horizonte completo;
- GP adicional para o runtime;
- penalização de custo e incentivo à fidelidade na aquisição;
- refino posterior em \(z=1\).

### 3.4 Resultado principal reportado

O workflow recupera uma frente de Pareto entre tracking e incrementos quadráticos de controle. Alguns candidatos igualam aproximadamente o custo de tracking de um controlador de referência com custo de incrementos 8–22 vezes menor no cenário de tuning. Em um schedule retido de 30 h, a vantagem qualitativa permanece, mas sua magnitude cai para aproximadamente 2,01 vezes e 1,12 vezes nas duas políticas avaliadas.

### 3.5 Limite de generalidade

As evidências atuais pertencem a:

- um único modelo nominal de biorreator CSTR;
- uma realização principal de ruído;
- duas parametrizações de controlador;
- um run de BO por parametrização;
- um schedule retido, mas ainda dentro do mesmo simulador nominal.

Esses limites precisam aparecer explicitamente no abstract, na discussão e nas conclusões.

## 4. Pontos fortes

### 4.1 Problema relevante e bem motivado

**OBSERVADO:** o manuscrito explica com clareza por que o tuning de NMPC é um problema black-box caro e por que o custo varia entre candidatos.

O problema é genuinamente adequado a PSE, otimização e controle de processos.

### 4.2 Formulação biobjetivo

**OBSERVADO:** tracking e atividade do atuador permanecem separados em vez de serem combinados por uma scalarização fixa.

Essa escolha é forte porque:

- expõe diretamente o trade-off;
- evita uma preferência escondida em pesos externos;
- permite uma interface de seleção para o engenheiro;
- oferece uma narrativa de decisão mais transparente.

### 4.3 Validação além do cenário de tuning

**OBSERVADO:** o artigo não se limita a reavaliar o mesmo ponto operacional. Ele usa um schedule multi-setpoint de 30 h e duas políticas de referência para o substrato.

Isso é melhor do que uma validação puramente in-sample, mesmo que ainda não constitua validade externa.

### 4.4 Refino em fidelidade completa

**REPORTADO:** candidatos são reavaliados em \(z=1\), e a frente refinada permanece próxima da frente estimada.

O refino é uma salvaguarda importante e deve ser preservado. O problema atual está na falta de rastreabilidade sobre quantos candidatos foram refinados e quanto isso custou.

### 4.5 Interpretação física

**OBSERVADO:** o artigo tenta relacionar os pesos \(R_u\), \(R_{\Delta u}\), \(Q\), \(N_c\) e \(N_p\) ao comportamento de tracking e suavização.

Essa interpretação dá legibilidade à frente de Pareto. Ela deve, contudo, ser apresentada como hipótese mecanística enquanto não for confirmada por ablações controladas.

### 4.6 Qualidade de escrita e figuras

O título é claro, a organização é boa, as equações são legíveis e as figuras têm aparência profissional. A narrativa central é compreensível e o manuscrito já possui maturidade editorial.

## 5. Achados críticos

## Critical 1 — A eficiência do MFBO não foi demonstrada

**Localização:** Abstract; Seções 3.3–3.4; Seção 4.1; Discussão; Conclusões.

**Problema:** há somente um run por parametrização, com 20 avaliações de inicialização e 101 avaliações adaptativas. Não existem baselines ou repetições.

Não foram comparados:

- qLogNEHVI somente em \(z=1\);
- MOBO truncado sem penalização de runtime;
- BO cost-aware sem multifidelidade;
- Sobol ou random search;
- outro método MFBO;
- múltiplas seeds do próprio método.

**Consequência:** a queda de runtime após a inicialização mostra que a aquisição passou a selecionar avaliações mais baratas. Ela não demonstra que o método produziu uma frente melhor por unidade de custo, nem que superou uma alternativa.

As afirmações *searched efficiently*, *made tuning tractable* e equivalentes são atualmente **REPORTADAS**, não demonstradas.

**AÇÃO MÍNIMA:** executar, com o mesmo orçamento de wall-clock ou CPU-hours:

1. Sobol;
2. qLogNEHVI em \(z=1\);
3. fidelidade variável sem runtime penalty;
4. método completo.

Usar pelo menos 5 seeds; 10 seria preferível. Reportar:

- hypervolume versus tempo;
- distância à frente validada;
- número de candidatos Pareto corretamente identificados;
- custo total;
- mediana e intervalos entre seeds.

**AÇÃO FORTE:** desenho de ablação 2×2:

| Configuração | Fidelidade/extrapolação | Penalização de runtime |
|---|---:|---:|
| A | Não | Não |
| B | Sim | Não |
| C | Não | Sim |
| D | Sim | Sim |

Adicionar uma aquisição MF ou cost-aware contemporânea como baseline.

## Critical 2 — A aquisição multifidelidade é heurística e não explicitamente orientada ao alvo \(z=1\)

**Localização:** Seções 3.1–3.4; Eqs. (17), (18), (20) e (21).

**OBSERVADO:** \(z\) aparece junto aos parâmetros do controlador no vetor de entrada do GP. O qLogNEHVI opera sobre objetivos já extrapolados para horizonte completo, e os termos \(a_z(z)\) e \(E[t(x,z)]\) recompensam fidelidade e penalizam custo.

**Problema conceitual:** em MFBO target-fidelity, a consulta \((x,z)\) deveria ser escolhida pelo valor da informação que fornece sobre \(f(x,1)\) ou sobre a frente Pareto em \(z=1\). A Eq. (20) não calcula explicitamente esse valor.

O algoritmo pode explorar um erro sistemático da extrapolação em certo \(z\), em vez de selecionar fidelidade pelo ganho informativo sobre o alvo final.

**Consequência:** chamar a abordagem de multifidelidade não é categoricamente errado, mas a contribuição é melhor descrita como uma variante heurística de MFBO/MOBO.

**AÇÃO MÍNIMA:** separar claramente:

- \(x\): parâmetros do controlador;
- \(z\): variável de consulta/fidelidade;
- \(f(x,z)\): observação truncada;
- \(f(x,1)\): objetivo de implantação.

Reposicionar o método como *runtime-aware MOBO with truncated closed-loop evaluations* e declarar a natureza heurística da Eq. (20).

**AÇÃO FORTE:** usar uma aquisição target-fidelity, por exemplo:

- MF-KG;
- MF-MES;
- expected improvement/hypervolume sobre o posterior em \(z=1\), normalizado pelo custo.

## Critical 3 — O Pareto precisa ser recalculado e auditado

**Localização:** algoritmo da Seção 3.4; Seções 4.2, 4.4 e 4.5; Figuras 4, 8, 10 e 11.

**OBSERVADO:** a Seção 3.4 exclui os pontos de inicialização do conjunto não dominado.

**Problema:** um ponto Sobol inicial pode dominar um ponto proposto pela aquisição. O Pareto empírico deve incluir todas as avaliações válidas.

Também aparecem contagens incompatíveis:

- 11 candidatos na frente combinada, sendo 6 + 5;
- 21 controladores em outra análise, sendo 11 + 10;
- dez controladores na validação;
- candidatos inicialmente não Pareto promovidos após o refino.

**Consequência:** não é possível reconstruir quais controladores sustentam cada figura e conclusão.

**AÇÃO MÍNIMA:** recomputar a dominância incluindo todos os pontos DOE e BO. Fornecer um ledger:

```text
242 avaliados
  -> Pareto estimado incluindo DOE
  -> candidatos enviados a z=1
  -> Pareto validado em z=1
  -> candidatos testados no schedule de 30 h
  -> controlador selecionado
```

Criar uma tabela suplementar por ID com:

- parâmetros;
- seed;
- \(z\) original;
- origem DOE/BO;
- objetivos extrapolados;
- objetivos em \(z=1\);
- status de factibilidade;
- pertencimento a cada frente.

## Critical 4 — A comparação 12D versus 8D não isola dimensionalidade

**Localização:** Seção 3.1; Seções 4.2–4.3; Discussão.

**OBSERVADO:** do Caso 1 para o Caso 2, o estudo:

- fixa \(q_1\);
- remove \(R_u\);
- altera os mecanismos de suavização;
- reduz a dimensão;
- muda a classe de controladores.

Além disso, \(R_u=0\) não pertence ao domínio do Caso 1, cujo limite inferior é \(10^{-3}\).

**Consequência:** as diferenças observadas podem resultar de:

- remoção de \(R_u\);
- fixação de \(q_1\);
- redução dimensional;
- seed;
- interações entre esses fatores.

Não é possível atribuir causalmente a especialização da frente à dimensionalidade.

**AÇÃO MÍNIMA:** chamar os casos de **duas parametrizações do controlador**, não de teste direto de dimensionalidade.

**AÇÃO FORTE:** fazer ablação:

1. espaço completo;
2. somente \(q_1\) fixo;
3. somente \(R_u=0\);
4. ambos;
5. múltiplas seeds para cada configuração.

## Critical 5 — Afirmações teóricas de unicidade, estabilidade e viabilidade são excessivas

**Localização:** Seções 2.1.1–2.1.2; Eqs. (1)–(6) e (10)–(11); validação com referência móvel.

**OBSERVADO:** o texto sugere que penalidades positivas em \(R_u\) ou \(R_{\Delta u}\) preservam unicidade e que a penalidade terminal baseada em LQR assegura estabilidade e viabilidade recursiva.

**Problema:** o NMPC apresentado é não linear e não convexo, resolvido localmente por SQP. Não há:

- conjunto terminal;
- lei terminal aplicada em região terminal;
- prova de invariância;
- prova de viabilidade recursiva;
- demonstração de unicidade global;
- tratamento formal da referência móvel.

**Consequência:** essas afirmações são potencialmente incorretas.

**AÇÃO MÍNIMA:** substituir por linguagem como:

- *promotes local regularisation and numerical conditioning*;
- *terminal penalty motivated by a local LQR design*;
- *formal stability and recursive-feasibility guarantees are outside the scope*.

Remover *preserves uniqueness* e *ensuring stability and recursive feasibility*.

**AÇÃO FORTE:** introduzir conjunto terminal, lei local, condições de invariância e uma prova compatível com input blocking e mudança de setpoints.

## Critical 6 — O extrapolador é o elo metodológico mais frágil

**Localização:** Seção 3.2; Eq. (17); Figura 1; Seção 4.4; Figuras 9–10.

**OBSERVADO:** dois polinômios de quinta ordem são ajustados offline usando apenas três trajetórias completas e aplicados a todos os controladores.

A construção pressupõe aproximadamente:

\[
J_{\mathrm{part}}(x,z) \approx \phi(z)J_{\mathrm{full}}(x).
\]

**Problema:** \(\phi(z)\) é tratada como universal por objetivo, embora controladores lentos, oscilatórios, saturados ou constraint-active possam apresentar frações acumuladas muito diferentes.

O \(R^2\) elevado é evidência insuficiente porque:

- curvas acumuladas são fortemente correlacionadas no tempo;
- somente três controladores são usados;
- não está clara a independência entre ajuste e validação;
- clipping em \([0.01,1]\) pode amplificar erros em baixa fidelidade;
- não há incerteza do extrapolador propagada ao GP.

**AÇÃO MÍNIMA:** informar:

- IDs e critérios de seleção das três trajetórias;
- coeficientes dos polinômios;
- custo offline do ajuste;
- independência em relação ao BO;
- erro por faixa de \(z\);
- validação controller-held-out.

Adicionar:

- Spearman e Kendall;
- inversões de ranking;
- precision/recall do Pareto;
- erro de hypervolume;
- erro relativo por fidelidade.

**AÇÃO FORTE:** aprender discrepancy dependente de \((x,z)\) por co-kriging, GP autoregressivo, modelo multi-task ou extrapolador hierárquico com incerteza.

## Critical 7 — Robustez e validade externa estão superalegadas

**Localização:** Abstract; Seção 3.2; Seção 3.7; Seção 4.5; Discussão; Conclusões.

**OBSERVADO:** a mesma realização de ruído é reutilizada entre candidatos. O schedule de validação usa o mesmo modelo nominal.

**Interpretação:** common random numbers é uma boa estratégia para comparações menos ruidosas. Porém, isso produz uma função aproximadamente determinística condicionada a uma única trajetória de ruído.

O estudo não inclui:

- seeds independentes de ruído;
- incerteza paramétrica;
- process noise;
- delay;
- bias de sensores;
- dinâmica de atuadores;
- plant–model mismatch;
- perturbações fora do modelo nominal.

**Consequência:** o schedule de 30 h demonstra generalização para um cenário nominal retido, não validade externa. Duas parametrizações com um run cada também não demonstram robustez à redução do espaço.

**AÇÃO MÍNIMA:** substituir:

- *external validity* por *held-out nominal-scenario transfer*;
- *robust* por *retained performance in the tested schedule*;
- *robust to how the space is reduced* por uma descrição factual das duas parametrizações.

**AÇÃO FORTE:** Monte Carlo dos candidatos refinados com:

- múltiplas noise seeds;
- perturbações em \(\mu_{\max}\), \(K_s\), \(k_d\), \(Y_{XS}\) e \(S_{\mathrm{in}}\);
- condições iniciais;
- atraso e bias;
- taxa de falhas e violações;
- média, dispersão, quantis e CVaR dos KPIs.

## Critical 8 — Runtime de tuning offline é confundido com viabilidade online

**Localização:** Abstract; Seções 3.3–3.5; Seções 4.1, 4.3 e 4.5; Discussão.

**OBSERVADO:** o surrogate de custo usa wall-clock da avaliação fechada completa em HPC, com pool MATLAB de 31 processos.

Esse custo inclui:

- número de instantes simulados;
- horizonte do NMPC;
- custo de cada NLP;
- paralelismo;
- overhead;
- possível contenção de hardware.

**Problema:** isso mede custo de tuning offline, não latência online do controlador.

**AÇÃO MÍNIMA:** usar consistentemente *offline closed-loop evaluation cost*.

**AÇÃO FORTE:** medir para cada controlador:

- tempo por solução NMPC;
- mediana;
- p95/p99;
- máximo;
- deadline misses com \(T_s=1\) min;
- iterações do SQP;
- falhas e recuperações;
- CPU-hours;
- hardware utilizado.

## Critical 9 — O custo do refino em \(z=1\) está ambíguo

**Localização:** Seções 4.1 e 4.4.

**VERIFY:** o manuscrito apresenta totais de 14,15 h e 21,96 h para os dois runs, mas também descreve 242 arquivos processados e reavaliações em fidelidade completa.

Se todos os 242 candidatos foram reexecutados em \(z=1\), esse custo pode ser substancial e parece não estar incluído nos totais principais. Se não foram reexecutados, deve ser explicado como foram obtidos seus endpoints full-horizon.

**AÇÃO MÍNIMA:** criar ledger de custo:

| Etapa | Número de avaliações | Wall-clock | CPU-hours | Necessária no workflow? |
|---|---:|---:|---:|---|
| Ajuste offline do extrapolador |  |  |  |  |
| DOE |  |  |  |  |
| Otimização adaptativa |  |  |  |  |
| Refino da frente em \(z=1\) |  |  |  |  |
| Auditoria de todos os pontos |  |  |  |  |
| Schedule de 30 h |  |  |  |  |

Separar custo necessário para implantação do custo adicional realizado apenas para validar o artigo.

## 6. Achados maiores

## Major 1 — O benchmark não é plenamente equivalente

**Localização:** Seção 4.3; Tabela 1; Conclusões.

O benchmark possui:

- \(N_p=61\);
- \(N_c=6\);
- \(P=0\);
- pesos e estrutura distintos;
- origem em trabalho anterior.

O controlador selecionado usa penalidade terminal e pertence a uma parametrização diferente.

**Consequência:** o resultado mostra que esse controlador específico é dominado no simulador e nos objetivos atuais. Não demonstra superioridade de MFBO sobre tuning convencional.

**AÇÃO:** incluir, sob o mesmo orçamento e arquitetura:

- single-fidelity BO;
- Sobol/random;
- tuning manual refeito;
- eventualmente CMA-ES ou NSGA-II;
- benchmark com o mesmo tratamento de \(P\).

Substituir *controller to be replaced* por formulação específica ao cenário testado.

## Major 2 — \(J_{\Delta u}\) não é total variation convencional

**Localização:** Eqs. (9) e (16); figuras; abstract; discussão.

\[
J_{\Delta u}=\sum_k\|u_k-u_{k-1}\|_2^2
\]

é soma quadrática dos incrementos, não total variation \(L_1\).

**Consequência:** dizer “8–22 vezes menos control variation” pode ser interpretado como redução linear do movimento físico do atuador, o que não decorre do custo quadrático.

**AÇÃO:** usar consistentemente:

- *sum of squared input increments*;
- *quadratic move activity*;
- *move energy*.

Complementar com:

- \(\sum|\Delta u|\);
- RMS de \(\Delta u\);
- pico de \(\Delta u\);
- tempo em saturação;
- número de mudanças relevantes.

## Major 3 — A vantagem no schedule retido é muito menor

**Localização:** Seções 4.3 e 4.5; Abstract; Conclusões.

No tuning, a redução reportada é 8–22 vezes. No schedule de 30 h, o selecionado apresenta aproximadamente 2,01 vezes e 1,12 vezes menos custo de incrementos que o benchmark.

**Consequência:** a vantagem qualitativa persiste, mas a magnitude contrai fortemente. O abstract e a conclusão atualmente aproximam demais os dois resultados.

**AÇÃO:** relatar explicitamente:

- 8–22× no cenário de tuning;
- 2,01× e 1,12× nas duas políticas do schedule;
- posição do controlador selecionado entre os demais Pareto;
- tracking e constraints por estado.

## Major 4 — Regra de seleção do controlador não está operacionalizada

**Localização:** Seção 3.7 e Seção 4.3.

A regra *lowest \(J_{\mathrm{track}}\) among the top three Pareto picks relative to the benchmark* é ambígua.

**AÇÃO:** definir antes da validação:

- como os “top three” são obtidos;
- distância ou tolerância em relação ao benchmark;
- preferência sobre \(J_{\Delta u}\);
- critério de knee;
- desempate;
- tratamento de incerteza.

Uma regra escolhida depois de ver o schedule retido geraria risco de seleção pós-hoc.

## Major 5 — Escala da Eq. (15) é descrita de forma inconsistente

Se o residual é multiplicado por \([10,1,1]\) antes da norma quadrática, o erro de volume recebe peso efetivo 100 no custo, não 10.

**AÇÃO:** distinguir:

- fator 10 aplicado ao residual;
- fator 100 aplicado à contribuição quadrática.

Revisar todas as interpretações state-wise.

## Major 6 — Redundância e não identificabilidade dos pesos

Uma multiplicação comum de todos os pesos positivos do objetivo idealmente não altera a política ótima. Isso cria uma direção redundante no espaço e dificulta a interpretação absoluta de \(Q\), \(R_u\) e \(R_{\Delta u}\).

**AÇÃO:** considerar:

- fixar um peso de referência;
- impor média geométrica unitária;
- otimizar razões;
- analisar condição numérica e falhas do SQP;
- separar invariância teórica de efeitos causados por tolerâncias numéricas.

## Major 7 — Mecanismos causais dos pesos e horizontes são apenas hipóteses

As associações observadas entre \(R_u\), \(R_{\Delta u}\), \(N_p\), \(N_c\) e as caudas da frente são interessantes. Porém, os casos diferem simultaneamente em vários fatores.

**AÇÃO:** testar os mesmos candidatos com intervenções:

- ligar/desligar \(R_u\);
- fixar \(N_p\) e \(N_c\);
- perturbar pesos em ±0,5 década;
- comparar predição linear e não linear;
- aplicar análise de sensibilidade ou functional ANOVA.

Até então, usar *suggests*, *is consistent with* ou *we hypothesise*, não linguagem causal.

## Major 8 — Segurança, factibilidade e falhas numéricas não são reportadas

O manuscrito fornece limites de estados e entradas, mas não apresenta:

- violações;
- slacks;
- infeasibilidade;
- terminações do SQP;
- saturações;
- avaliações falhas;
- política de penalização na BO;
- convergência numérica.

**AÇÃO:** tabela para benchmark e candidatos finais contendo:

- IAE/ISE por estado;
- overshoot;
- settling;
- violações máximas e acumuladas;
- saturação de entradas;
- falhas de NLP;
- p95 de solve time;
- offset;
- status terminal.

## Major 9 — A referência virtual negativa para substrato precisa de justificativa

**OBSERVADO:** \(S_{sp}(X)\) pode ficar negativo, embora \(S\geq0\).

Isso pode ser uma referência virtual deliberada, mas deve ser explicado como tal. Também torna ainda mais frágil qualquer afirmação de garantia formal associada ao alvo terminal.

**AÇÃO:** justificar a referência negativa e sua segurança, ou usar clamp em \([0,3]\). Fazer sensibilidade das duas escolhas.

## Major 10 — Método GP/BO insuficientemente reproduzível

Faltam detalhes sobre:

- kernel e ARD;
- likelihood e noise floor;
- dependência/independência entre objetivos;
- priors;
- MC samples;
- raw samples e restarts;
- batch \(q\);
- tolerâncias;
- versões de BoTorch, PyTorch, MATLAB e solver;
- tratamento de duplicatas após rounding;
- avaliações inviáveis;
- seeds;
- warm start e reinicialização.

**AÇÃO:** adicionar pseudocódigo reprodutível e tabela completa de hiperparâmetros.

## Major 11 — Possível ambiguidade em `qLogNEHVI`

**VERIFY:** em BoTorch, o objeto `qLogNEHVI` já retorna uma aquisição no domínio log. A equação do manuscrito escreve algo equivalente a \(\log(\mathrm{qLogNEHVI})\).

Isso pode ser apenas uma imprecisão de notação ou um segundo log na implementação.

**AÇÃO:** auditar o código. Se o objeto já retorna log-NEHVI, escrever:

\[
\log \alpha = \mathrm{qLogNEHVI}_{\text{output}}+\ldots
\]

e não \(\log(\mathrm{qLogNEHVI})\).

## Major 12 — Variáveis inteiras por relaxação e arredondamento

O acquisition optimizer opera na relaxação contínua e depois arredonda horizontes. Isso pode:

- degradar o máximo da aquisição;
- gerar duplicatas;
- avaliar gradientes em pontos sem significado físico;
- enviesar a seleção perto de fronteiras inteiras.

**AÇÃO:** enumerar combinações admissíveis de horizontes ou usar aquisição mixed-variable. No mínimo, reportar duplicate rejection e reotimização após snapping.

## Major 13 — Referência de hypervolume dinâmica

A referência é atualizada com os dados correntes. Isso pode ser aceitável para a aquisição, mas impede comparar diretamente valores de hypervolume ao longo de iterações, métodos e seeds.

**AÇÃO:** usar referência fixa pós-hoc para todos os gráficos e testes comparativos.

## 7. Achados menores e itens de verificação

### Minor 1

A expressão “within-run frontiers are not true Pareto sets” é imprecisa. Elas são conjuntos não dominados dentro de cada amostra; apenas podem deixar de ser globalmente não dominadas após a união.

### Minor 2

Definir precisamente o arredondamento em \(N(z)\), o caso \(z=0.01\), o número mínimo de passos e qualquer efeito off-by-one.

### Minor 3

Há referência cruzada ao modelo de ruído como introduzido na Seção 3.1, embora ele apareça na Seção 3.2.

### Minor 4

“Two independent BO runs” pode sugerir replicações independentes. Na prática, são dois casos únicos com parametrizações distintas e noise realization compartilhada.

### Minor 5

O critério de knee é visual. Se “pronounced knee” for central, definir um método de curvatura, distância ou knee detection.

### Minor 6

O abstract possui aproximadamente 257 palavras e está denso. Deve ser encurtado e ter as alegações de eficiência, robustez e tratabilidade moderadas.

### Minor 7

Padronizar:

- \(J_{TV}\);
- \(J_{\Delta u}\);
- *control variation*;
- *input movement*;
- *move energy*.

### Minor 8

As captions são informativas, mas algumas estão longas. Alguns rótulos e detalhes das figuras ficam pequenos em página de duas colunas.

### Minor 9

O manuscrito ainda contém:

- “Author Name”;
- afiliação incompleta;
- Acknowledgements vazio.

### Minor 10

Faltam:

- Data Availability;
- Code Availability;
- Funding;
- Conflict of Interest;
- CRediT;
- eventual declaração de uso de IA conforme a política do periódico.

### Verify 1

Explicar como \(P\) é tratado quando \(A_{cl}\) não for Schur ou quando a equação de Lyapunov não produzir uma matriz adequada.

### Verify 2

Esclarecer a calibração do GP de runtime e separar ruído intrínseco do candidato de variação de carga da infraestrutura HPC.

### Verify 3

Auditar referências de 2025–2026, preprints, DOIs e status de publicação antes da submissão. Esta avaliação não constitui auditoria bibliográfica completa.

## 8. Novidade e posicionamento

### 8.1 Minha avaliação

A novidade é **boa como integração de engenharia** e **moderada a baixa como avanço algorítmico fundamental**.

O núcleo potencialmente publicável é:

> Um workflow de otimização Bayesiana multiobjetivo, sensível a runtime, que usa simulações fechadas truncadas e refino em horizonte completo para tuning de NMPC.

O artigo não demonstra ainda que a Eq. (20) é uma nova aquisição superior, nem que a estratégia multifidelidade causa a economia observada.

### 8.2 Literatura que precisa ser discutida

A revisão deve posicionar o trabalho em relação a:

- trace-aware/multi-fidelity knowledge gradient;
- early stopping e extrapolação de curvas;
- cost-aware BO contemporâneo;
- MFBO target-fidelity;
- controller tuning com digital twins;
- MFBO de simulação para experimento;
- cost-aware multiobjective BO.

Fontes relevantes para posicionamento:

- Wu et al., *Practical Multi-fidelity Bayesian Optimization for Hyperparameter Tuning*:  
  https://proceedings.mlr.press/v115/wu20a/wu20a.pdf
- Nobar et al., *Guided Multi-Fidelity Bayesian Optimization for Data-Driven Controller Tuning With Digital Twins*:  
  https://www.research-collection.ethz.ch/entities/publication/90245e4c-a04a-4558-8984-c24645f0416c
- Zhao et al., *Efficient Learning of Vehicle Controller Parameters via Multi-Fidelity Bayesian Optimization*:  
  https://arxiv.org/abs/2506.08719
- Xie et al., *Cost-aware Bayesian Optimization via the Pandora's Box Gittins Index*:  
  https://proceedings.neurips.cc/paper_files/paper/2024/hash/d14c355d5e88cff437a6303d2d716252-Abstract-Conference.html

Essas referências são sugestões de posicionamento; sua inclusão deve ser confirmada pelos autores após leitura integral e avaliação de suporte às afirmações específicas.

## 9. Avaliação editorial

### 9.1 Título

O título é forte e informativo. Se o método permanecer heurístico, considerar:

> *Runtime-Aware Multi-Objective Bayesian Optimisation with Truncated Closed-Loop Evaluations for NMPC Tuning*

O título atual pode ser mantido se a semântica multifidelidade for cuidadosamente formalizada.

### 9.2 Abstract

O abstract apresenta método, resultados e conclusão, mas promete mais do que a evidência sustenta.

Até que existam baselines e repetições, remover ou moderar:

- *efficiently*;
- *robust*;
- *operationally tractable*;
- *external validity*;
- qualquer sugestão de superioridade geral.

### 9.3 Narrativa

A pergunta secundária 12D versus 8D ocupa espaço demais e compete com a contribuição principal. Como está confundida com a mudança de estrutura do controlador, há duas opções:

1. rebaixá-la a análise exploratória secundária;
2. redesenhar como ablação fatorial e elevá-la a contribuição real.

### 9.4 Journal fit

**Computers & Chemical Engineering** apresenta excelente encaixe temático, incluindo otimização, computação, dinâmica e controle:

https://www.sciencedirect.com/journal/computers-and-chemical-engineering

Entretanto, o periódico provavelmente exigirá evidência comparativa mais forte se o artigo for apresentado como contribuição metodológica.

**Digital Chemical Engineering** pode ser uma opção mais segura para posicionamento como workflow aplicado/case study:

https://www.sciencedirect.com/journal/digital-chemical-engineering

Uma revista de controle mais teórico provavelmente exigiria provas, baselines e análises de estabilidade substancialmente mais fortes.

## 10. Plano de revisão recomendado

### Prioridade 0 — Integridade do resultado

1. Recalcular o Pareto incluindo DOE.
2. Reconciliar 11, 21, 10, 202 e 242 candidatos/avaliações.
3. Criar ledger de candidatos e de custo.
4. Confirmar quais pontos foram reavaliados em \(z=1\).
5. Propagar as correções para figuras, abstract e conclusão.

### Prioridade 1 — Correções científicas

1. Remover alegações incorretas de unicidade, estabilidade e viabilidade.
2. Separar \(x\) e \(z\) conceitualmente.
3. Reposicionar a aquisição como heurística ou implementar target-fidelity.
4. Chamar os casos de duas parametrizações.
5. Corrigir terminologia de \(J_{\Delta u}\).

### Prioridade 2 — Evidência experimental

1. Rodar Sobol, single-fidelity MOBO e ablações.
2. Usar orçamento comum em wall-clock/CPU-hours.
3. Executar múltiplas seeds.
4. Reportar hypervolume com referência fixa.
5. Validar ranking e classificação Pareto do extrapolador.

### Prioridade 3 — Controle e robustez

1. Monte Carlo com noise seeds independentes.
2. Perturbação de parâmetros e mismatch.
3. Constraints, saturação, falhas e settling.
4. Tempo online por NLP e deadline misses.
5. Métricas por estado e por atuador.

### Prioridade 4 — Reprodutibilidade

1. Versões de software e solver.
2. Kernels, priors, tolerâncias e seeds.
3. Política de falhas e duplicatas.
4. Dados brutos e processados.
5. Checkpoints da otimização.
6. Scripts para reproduzir figuras e tabelas.

### Prioridade 5 — Reescrita final

1. Atualizar literatura.
2. Reescrever abstract.
3. Reorganizar discussão.
4. Moderar generalizações.
5. Completar declarações editoriais.

## 11. Matriz resumida de decisão

| ID | Severidade | Tema | Estado | Ação necessária |
|---|---|---|---|---|
| C1 | Critical | Eficiência | Não demonstrada | Baselines, seeds e orçamento comum |
| C2 | Critical | Semântica MFBO | Heurística | Reposicionar ou usar target-fidelity |
| C3 | Critical | Pareto | Inconsistente | Incluir DOE e criar ledger |
| C4 | Critical | 12D vs 8D | Confundido | Chamar de parametrizações ou fazer ablação |
| C5 | Critical | Garantias de controle | Excessivas/incorretas | Remover ou provar |
| C6 | Critical | Extrapolador | Evidência frágil | Holdout, ranking e incerteza |
| C7 | Critical | Robustez | Superalegada | Moderar ou executar Monte Carlo |
| C8 | Critical | Runtime | Offline vs online | Separar e medir solve time |
| C9 | Critical/Verify | Custo de refino | Ambíguo | Ledger de custo end-to-end |
| M1 | Major | Benchmark | Não equivalente | Baseline comparável |
| M2 | Major | Métrica de controle | Nomenclatura inadequada | Renomear e complementar |
| M3 | Major | Schedule de 30 h | Efeito reduzido | Reportar magnitude corretamente |
| M4 | Major | Seleção | Regra ambígua | Pré-especificar decisão |
| M5 | Major | Escala | Fator 10 vs 100 | Corrigir interpretação |
| M6 | Major | Pesos | Não identificáveis | Normalizar ou otimizar razões |
| M7 | Major | Mecanismos | Inferenciais | Ablar e moderar causalidade |
| M8 | Major | Segurança | Evidência ausente | Reportar constraints e falhas |
| M9 | Major | Referência negativa | Precisa justificativa | Justificar ou limitar |
| M10 | Major | Reprodutibilidade BO | Incompleta | Especificar implementação |
| M11 | Major/Verify | qLogNEHVI | Possível duplo log | Auditar código e notação |
| M12 | Major | Inteiros | Relaxação/rounding | Mixed-variable ou enumeração |
| M13 | Major | Hypervolume | Referência dinâmica | Referência fixa para comparação |

## 12. Recomendação final

O manuscrito apresenta uma integração promissora e potencialmente útil de otimização Bayesiana multiobjetivo, truncamento de horizonte e modelagem de custo para tuning de NMPC. A evidência atual mostra que o workflow encontrou controladores interessantes no estudo nominal, mas ainda não demonstra que a componente multifidelidade/runtime-aware é responsável pelo ganho, nem robustez ou garantias formais de estabilidade.

A recomendação é **major revision**, com:

- correção das alegações teóricas;
- reconciliação completa do pipeline Pareto/full-fidelity;
- contabilização end-to-end do custo;
- ablações replicadas contra baselines sob o mesmo orçamento;
- validação mais forte do extrapolador;
- separação entre custo de tuning e viabilidade online;
- moderação das afirmações de robustez e validade externa.

Após essas revisões, o trabalho pode se tornar uma contribuição sólida e interessante para otimização e controle em Process Systems Engineering.

