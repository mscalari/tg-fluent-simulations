# Curvas de desempenho de conversores de torque na literatura

Arquivo de dados: `refs/conversor_literatura.csv` (separador `;`, ponto decimal, 115 linhas).
Uso: comparação de **tendência adimensional** (TR x SR, eta x SR, K/K_estol x SR). Os conversores são
diferentes do conversor do TG, então não se comparam valores absolutos.

Definições do TG: SR = w_turbina/w_bomba, TR = -T_turbina/T_bomba, eta = SR*TR, K = n_bomba/sqrt(T_bomba).

Todas as curvas foram **lidas manualmente de figuras** (renderização do PDF a 300-450 dpi, calibração dos eixos
pelos tiques e leitura do centro dos marcadores). Nenhum dos artigos obtidos traz tabela numérica das curvas de
desempenho. Nenhum valor foi tirado de memória.

## Convenção da coluna K

- Os artigos de Zhang, Guo e Ran não dão K = n/sqrt(T) e sim o **fator de capacidade adimensional CC**
  (CC = T_bomba / (rho * N_bomba^2 * D^5); Ran escreve com g e unidade 1e-6 min^2 rev^-2 m^-1). Como CC é proporcional a
  T/n^2, vale **K ~ 1/sqrt(CC)**, logo **K/K_estol = sqrt(CC_estol / CC)**. Esse valor está no campo `obs`
  como `K/K_estol(derivado)`. Para o estol, usar SR = 0 (Guo, Ran) ou SR = 0.01 (Zhang).
- Nadeem: K_cap = n_bomba[rpm]/sqrt(T_bomba[N.m]) (definição da tese, eq. 3.17), igual ao K do TG.
- `eta`: onde está vazio, a figura não traz a curva ou ela é inconsistente com SR*TR; use SR*TR. Em Guo, Ran e Nadeem
  o eta do CSV é o produto SR*TR (derivado, indicado em `obs`).

## Fontes obtidas

### 1. zhang2023 (10.3390/machines11040489, Machines 11(4):489) - 30 linhas
- Título: Simulation and Validation of Cavitating Flow in a Torque Converter with Scale-Resolving Methods
  (Zhang, Yan, Liu, Guo, Wei; Beijing Institute of Technology).
- Acesso: PDF aberto (mdpi-res.com), texto e figuras lidos.
- Extraído: **Figura 8 (p. 13 de 26)**: K (= TT/TP, ou seja, o TR), eta (%) e CC versus SR, em 10 pontos
  (SR = 0.01, 0.1, ..., 0.9). Origens: `exp` (teste), `cfd_rans` (SST) e `cfd_les` (WALE). SBES, DDES e SAS-SST estão
  na figura mas não foram digitalizados (curvas sobrepostas, sem ganho para o objetivo).
- Condições (Tabela 2 p. 7, Tabela 3 e texto p. 8-9): N_bomba = 2000 rpm, N_estator = 0, N_turbina = 20-1800 rpm;
  óleo rho = 835.2 kg/m3, mu = 0.0146 Pa.s; pressão de referência 0.4 MPa; pás 29/25/22 (bomba/turbina/estator).
  **O diâmetro não é dado.** O artigo-irmão de Guo (item 2) usa o mesmo óleo, rotação e pressão com o conversor de
  315 mm; isso sugere o mesmo conversor, mas é inferência, não está escrito.
- CFD: ANSYS CFX, cavitação (Zwart), SST (RANS), WALE (LES) e híbridos SBES/DDES/SAS.
- Leitura: incerteza estimada TR +-0.015 (exp) e +-0.02 (CFD), eta +-0.7 ponto percentual, CC +-0.05 (exp) e +-0.07 (CFD).
- Limitações: o último ponto é SR = 0.9 com TR = 0.865 (< 1), então não há ponto de acoplamento nem o platô de
  TR ~ 1. O SR mínimo é 0.01. A curva de eta dos modelos CFD não é consistente com SR*TR dos próprios modelos
  (por exemplo SST em SR = 0.4), por isso o eta CFD ficou em branco. Os dados experimentais são consistentes
  (eta = SR*TR confere).
- TR de estol (exp, SR = 0.01): **2.155**. eta máximo (exp): ~0.85 em SR = 0.7. TR cruza 1 em SR ~ 0.83
  (interpolação entre 0.8 e 0.9, sem ponto de acoplamento verdadeiro).

### 2. guo2021 (10.3390/app11083643, Appl. Sci. 11(8):3643; no pedido aparece como "li2021") - 27 linhas
- Título: Evaluation and Validation of Viscous Oil Cavitation Model Used in Torque Converter
  (Guo, Liu, Yan, Ke, Wei, Li). **Os autores não incluem Li como primeiro autor; o DOI do pedido corresponde a este artigo.**
- Acesso: PDF aberto (mdpi-res.com).
- Extraído: **Figura 7 (p. 9 de 23)**: K (TR), eta e CC versus SR, SR = 0 a 0.8 (9 pontos). Origens: `exp` (K_Test, CC_Test),
  `cfd_rans` sem cavitação (K_No Cavitation) e `cfd_rans` com cavitação (K_Cavitation, CC_Cavitation). O CC "sem cavitação"
  só é legível em SR <= 0.3 (curvas sobrepostas) e não foi gravado.
- Geometria e condições (Seção 2.2, 3.1, Tabela 2): conversor base YD315, diâmetro do torus **315 mm**; N_bomba = 2000 rpm,
  N_turbina = 0-1600 rpm, estator travado; óleo rho = 835.2 kg/m3, mu = 0.0146 Pa.s; pressão interna 0.4 MPa; óleo a 90 C;
  bancada com motor/gerador de 300 kW, sensores com exatidão de 0.25%.
- CFD: ANSYS CFX 19.2, RANS SST, **permanente com frozen rotor** (análogo ao MRF do TG), com e sem cavitação (Zwart).
- Leitura: TR +-0.01, CC +-0.05; eta = SR*TR (a curva eta_Test coincide visualmente).
- Limitações: SR máximo 0.8; sem acoplamento. A leitura de CC_Test em SR = 0.3-0.4 tem a incerteza maior
  (marcadores sobrepostos à linha tracejada). Checagem cruzada com o texto: erro de CC no estol 6.7% (lido: 7.2%)
  e queda de CC de 21.9% por cavitação (lido: 22%). Isso valida a calibração dos eixos.
- TR de estol: exp **2.175**, CFD sem cavitação 2.26, CFD com cavitação 2.175. TR cruza 1 em SR ~ 0.86 (extrapolação,
  fora da faixa medida). eta no último ponto (SR = 0.8): 0.86.

### 3. ran2021 (10.3390/app11062798, Appl. Sci. 11(6):2798) - 27 linhas
- Título: 3D Cavitation Shedding Dynamics: Cavitation Flow-Fluid Vortex Formation Interaction in a Hydrodynamic Torque Converter
  (Ran, Ma, Liu; Jilin University).
- Acesso: PDF aberto (mdpi-res.com).
- Extraído: **Figura 6a e 6b (p. 9 de 24)**: TR e CC versus SR, SR = 0 a 0.8 (9 pontos). Origens: `exp` (Test), `cfd_rans`
  com cavitação (Cavi) e `cfd_rans` sem cavitação (NOC).
- Condições (Tabela 2 p. 8, Tabela 3): N_bomba = 2000 rpm, N_turbina 0-1600 rpm, estator parado; óleo
  rho = 860 kg/m3, mu = 0.0258 Pa.s; pressão de carga 0.8 MPa, contrapressão 0.4 MPa; pás 29/24/22. **Diâmetro do torus não informado.**
  Os dados de teste e a bancada são remetidos a Liu et al. (2015), Math. Probl. Eng. (ref. [26]).
- CFD: ANSYS FLUENT, SST k-omega, transiente, passo 1e-4 s, com modelo Zwart-Gerber-Belamri.
- Leitura: TR +-0.01, CC +-0.05; eta = SR*TR (não há curva de eta no artigo).
- Limitações: SR máximo 0.8; sem acoplamento; a curva NOC fica coberta pela Cavi em SR >= 0.3, então esses
  pontos NOC são cópia da leitura Cavi. Texto confirma erro de CC no estol de 13.23% (NOC) e < 5% (Cavi).
- TR de estol: exp **2.40**, CFD Cavi 2.48, NOC 2.57. TR cruza 1 em SR ~ 0.81 (extrapolação linear, fora da faixa medida).

### 4. nadeem2020 (tese de mestrado, Michigan Technological University, 2020) - 22 linhas
- Título: Design of an Anti-Jerk Controller for Both Locked and Slipping Torque Converter Conditions in a Vehicle (Syed Ahmad Nadeem).
  PDF aberto em https://sites.ualberta.ca/~mahdi/Docs/A_Nadeem_MScThesis.pdf. **Não estava na lista do pedido**; foi achada
  na busca pelos coeficientes de Kotwicki.
- Extraído: **Figuras 3.6 e 3.7 (p. 76)**: TR e K_cap versus SR em 22 pontos (SR = 0 ... 0.9495). Os SR de cada ponto vêm
  dos rótulos do eixo (0, 0.1001, 0.2001, ..., 0.9495) e os valores foram lidos nos marcadores.
- Condições: conversor automotivo de um projeto patrocinado; os dados "experimentais" foram fornecidos pelo patrocinador
  (Seção 3.4.1.2) e replotados na tese. **Fonte secundária; diâmetro, óleo e rotação não informados.**
  Marcada como `exp` no CSV com essa ressalva.
- Leitura: TR +-0.01, K_cap +-0.5 (rpm/sqrt(N.m)).
- Limitações: TR lido ~0.99 para SR >= 0.93 (provável TR = 1.00 real; erro de leitura). Cobre o regime de conversor e a
  transição para o acoplamento, que é o intervalo que falta nos artigos acima.
- TR de estol: **1.85**; K_cap de estol: 164.9 rpm/sqrt(N.m) (mínimo 148.7 em SR = 0.5, depois sobe e dispara depois de SR = 0.93).
  SR de acoplamento (TR -> 1): **~0.92** (TR = 1.0 em SR = 0.9196-0.9297); a partir daí a tese considera o conversor travado.

## Modelo de Kotwicki

- Referência original: Kotwicki, A. J., "Dynamic Models for Torque Converter Equipped Vehicles", SAE Technical Paper 820393, 1982,
  DOI 10.4271/820393. **Fonte original NÃO obtida** (paga na SAE MOBILUS; só a página de resumo foi acessível, sem coeficientes).
- Forma do modelo citada em Nadeem (2020), Seção 3.4.1.3, eqs. 3.20-3.22, p. 78-79 (refere-se a Kotwicki, ref. [36] da tese).
  Modo conversor (estator preso pela roda livre): T_bomba = a1*w_b^2 + a2*w_b*w_t + a3*w_t^2 e
  T_turbina = b1*w_b^2 + b2*w_b*w_t + b3*w_t^2. Modo acoplamento: T_bomba = T_turbina = c1*w_b^2 + c2*w_b*w_t + c3*w_t^2.
- **Coeficientes disponíveis** (Nadeem 2020, Tabelas 3.1 e 3.2, p. 80-81, unidades N.m/(rad/s)^2; modo conversor, SR 0-0.95):

  | Coeficiente | Bomba (impeller) | Turbina |
  |---|---|---|
  | 1 | 0.77e-3 | 5.04e-3 |
  | 2 | 12.44e-3 | 6.36e-3 |
  | 3 | -11.52e-3 | -9.94e-3 |

  **Esses valores NÃO são de Kotwicki.** Foram ajustados por Nadeem (estimador de parâmetros do MATLAB) a um conversor
  automotivo específico, a partir de sua tabela SR-TR-K_cap. Os coeficientes do modo acoplamento não foram identificados na tese.
  Consequência para a tendência: dividindo as formas, TR(SR) = (b1+b2*SR+b3*SR^2)/(a1+a2*SR+a3*SR^2) e
  K(SR) = (60/(2*pi))/sqrt(a1+a2*SR+a3*SR^2). Conferido contra a própria curva da tese: bom ajuste para SR >= 0.3
  (por exemplo SR = 0.5: modelo TR = 1.395 e K = 149.0 contra lido 1.389 e 148.7), mas **ruim perto do estol** (modelo TR(0) = 6.5
  contra 1.85 medido). Por isso o CSV traz linhas `modelo` só para SR = 0.3-0.95 (fonte `kotwicki_forma_nadeem2020`).
  A ressalva vale também para as faixas de SR > 0.9: o modelo dá TR = 0.96 em SR = 0.95.
- Outras fontes tentadas sem sucesso para os coeficientes originais de Kotwicki: artigos de identificação online e patentes
  da Ford citados na busca usam a forma quadrática, mas com coeficientes ajustados a conversores específicos; não foram incluídos.

## Fontes não obtidas

| Fonte | DOI | Motivo |
|---|---|---|
| liu2015 (Liu, Liu, Ma. RANS, DES e LES de escoamentos internos de conversores, EACFM 9(1):1-12) | 10.1080/19942060.2015.1004814 | Taylor & Francis bloqueia download automatizado (HTTP 403 e verificação anti-bot do Cloudflare). Não foi contornado. Os dados não foram vistos. |
| liu2018 (Liu et al. LES em conversor hidrodinâmico, EACFM) | 10.1080/19942060.2018.1489896 | idem |
| liu2019 (Desenvolvimento e validação de procedimento de otimização CFD para cascata de conversor, EACFM) | 10.1080/19942060.2018.1562383 | idem |
| jung2011 (Jung, Kang, Hur. Progress in CFD 11:261-268) | 10.1504/PCFD.2011.041027 | Paga (Inderscience, HTTP 403). Só o resumo foi lido em repositório institucional: compara frozen rotor, sliding mesh e mixing plane com dados experimentais e estuda a variação da viscosidade com a temperatura. **Sem números de TR, eta ou K.** Nada foi gravado no CSV. |
| Kotwicki 1982 (SAE 820393) | 10.4271/820393 | Paga; só o resumo está acessível. Ver seção do modelo. |

Os três artigos da família Liu (Jilin University) são acesso aberto (CC-BY); basta baixar manualmente o PDF no navegador e
passar o caminho para digitalizar as figuras. As informações de apoio indicam o mesmo grupo (C. Liu, W. Ma) de Ran 2021,
mas isso não foi verificado.

## Observações para o gráfico de tendência

- As curvas experimentais de Zhang, Guo e Ran descem em TR de 2.15-2.40 (SR ~ 0) para ~1.0 em SR 0.8-0.9, com
  eta máximo próximo de SR 0.7-0.8; o fator de capacidade (CC) sobe ligeiramente em SR baixo (efeito da cavitação no estator)
  e cai depois de SR ~ 0.3. A curva de Nadeem tem TR de estol menor (1.85) e K_cap com mínimo em SR ~ 0.5 (equivalente ao máximo de CC).
- A cavitação faz o CC (e o K_cap = n/sqrt(T)) mudar de sentido em baixo SR; um CFD sem cavitação, como o MRF do TG,
  deve reproduzir o CC_NOC (decrescente de SR = 0 para 0.8 em Ran e Guo), não o teste.

## Revisão independente (2026-11-21, agente Sonnet, por amostragem; Nadeem conferido ponto a ponto no vetor do PDF)

- Todos os pontos conferidos estão dentro da incerteza declarada; séries e rótulos corretos; K/K_estol = √(CC_estol/CC) correto em Zhang, Guo e Ran; coeficientes de Nadeem (Tab. 3.1/3.2) transcritos corretamente.
- **Corrigido:** K/K_estol de Nadeem estava gravado como √(K_estol/K). Como K_cap = n/√T já é o análogo de 1/√CC, o correto é K/K_estol = K_cap/164,9 (agora no campo `obs`).
- **Corrigido:** TR < 1 em SR 0,93–0,95 em Nadeem não é erro de leitura: o dado vetorial tem TR = 0,994–0,996 (patamar real).
- **Corrigido:** definição de CC: Zhang e Ran usam T_P/(ρ·g·N²·D⁵); Guo usa T_P·D⁻⁵·ω_P⁻² (sem ρ). Não afeta K/K_estol.
- **Acrescentado:** CC sem cavitação de Guo em SR 0–0,3 (leitura do revisor: 6,63 / 6,56 / 6,40 / 6,21, ±0,05), com K/K_estol.
- O eixo X de Nadeem é categórico: os 22 pontos são equidistantes na figura, mas os SR são irregulares. Usar os SR do CSV, nunca a posição.
- Ran: a legenda "Figure 6" está na p. 10 do PDF; as imagens, na p. 9.
