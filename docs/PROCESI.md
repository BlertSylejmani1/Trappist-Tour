# Procesi i punës: nga fillimi (64k) deri në 5.815M

Ky dokument përshkruan rrugëtimin e projektit me radhë: çfarë provuam, çfarë funksionoi, çfarë jo,
dhe pse. Kodi i fazave të hershme dhe i eksperimenteve nuk është në repo; këtu ruhen vetëm vlerat
dhe përfundimet. Në repo është vetëm versioni final.

## Përmbledhje

| Faza | Ideja kryesore | Hypervolume |
|---|---|---:|
| 0 | Algoritme evolucionare (NSGA-II, R-NSGA-II, CMA-ES) mbi gjithë vektorin 34-D | 2.33M |
| 1 | Analiza fizike + beam search këmbë pas këmbe | 5.792M |
| 2 | Heuristikë për të ardhmen (si A\*) + rioptimizimi i fundit të turës | 5.812M |
| 2b | Bashkimi i të gjitha run-eve, me verifikim zyrtar dhe zgjedhje optimale të 100 pikave | **5.815M** |
| 3 | Eksperimente për vendin e parë (rezultat negativ) | – |

Për krahasim, renditja e garës zyrtare: 5.870M (vendi i 1-rë), ~5.73M (i 2-ti), ~5.60M (i 3-ti).

## Faza 0 – Pikënisja: algoritme evolucionare mbi gjithë vektorin (2.33M)
Qasja e parë e trajtoi problemin si optimizim të përgjithshëm shumë-objektivësh me variabla të
përziera (7 numra të plotë për radhën e planetëve + 27 të vazhdueshëm), me vlerësuesin zyrtar të
ESA-s të pandryshuar. Algoritmet u shtuan njëri pas tjetrit, sa herë që i mëparshmi ngecte:

| Eksperimenti | Ndryshimi në metodë | Hypervolume | Interpretimi |
|---|---|---:|---|
| Shembulli i ESA-s | një zgjidhje e vlefshme e dhënë nga ESA | 64,155 | vetëm pikë referimi |
| 2 | NSGA-II me kryqëzim/mutacion që ruan permutacionin + kërkim lokal Nelder–Mead | 428,506 | përmirësimi i parë i dukshëm |
| 3 | + rafinim adaptiv Gaussian (i frymëzuar nga CMA-ES) | 440,572 | përmirësim i vogël |
| 4 | ndërtues i parë këmbë pas këmbe | 0 | i refuzuar: asnjë zgjidhje brenda kutisë së referimit |
| 5 | R-NSGA-II (NSGA-II me pika referimi) mbi radhë të fiksuar planetësh | 1,048,562 | hap i madh: kërkimi përqendrohet te zona që shpërblen score-i |
| 6 | R-NSGA-II më i madh + CMA-ES | 1,105,301 | +5.4% |
| 7 | "zbulim i permutacioneve": testohen të gjitha radhët që nisin me planetin b, pastaj R-NSGA-II për familjet më të mira | 1,449,988 | +31%; katër familje të vlefshme rrugësh |
| 8 | buxhet më i madh, shumë familje rrugësh + CMA-ES | 2,285,569 | +58% |
| 9 | kërkim lokal "triplet" (3 variabla njëkohësisht mbi dritare 3-planetëshe) | **2,331,795** | +2%; kërkimi u përqendrua te një familje e vetme |
| 10 | 8 radhë planetësh nga literatura, secila e zgjidhur veçmas | 2,098,938 | nuk e kaloi eksperimentin 9 |

Zgjidhja më e mirë e kësaj faze: **2,331,795**, me ~437 m/s karburant te turat më të lira (~1855 ditë).

Idetë kryesore të kësaj faze erdhën nga literatura: R-NSGA-II dhe ndarja sipas permutacioneve
(P. Heine, Bielefeld, 2022), optimizimi me radhë të fiksuar planetësh dhe kërkimi lokal "triplet"
(Shimane & Iiyama, AAS 23-431, 2023), dhe kombinimi algoritëm evolucionar → CMA-ES (Honda Research
Institute, SpOC 2022).

**Mësimi:** algoritmet evolucionare mbi gjithë vektorin 34-D përparojnë ngadalë. Problemi është shumë
"kaotik": një ndryshim i vogël në fillim të turës e ndryshon krejtësisht fundin, dhe shumica e
vlerësimeve shkojnë te trajektore me faza të gabuara (anija mbërrin kur planeti s'është aty).
Përmirësimet e mëdha erdhën vetëm kur kërkimi u bë më i strukturuar (radhë e fiksuar, pika
referimi), dhe kjo tregoi drejtimin për fazën tjetër: të shfrytëzohet struktura e problemit
**këmbë pas këmbe**.

## Faza 1 – Analiza fizike dhe beam search (5.792M)
**Vlerësues i shpejtë.** Rishkruam funksionin zyrtar me numba: rreth 4 herë më i shpejtë, me gabim
~1e-13 krahasuar me pykep. Pa këtë, miliona vlerësimet e kërkimit nuk do të ishin të mundshme.

**Analiza fizike** ndryshoi drejtimin e punës:
- Hyrja nga 10 AU me 10 km/s zgjat natyrshëm ~1500 ditë pa karburant; shkurtimi kushton ~6–9 m/s/ditë.
- Anija mbërrin me ~118 km/s. Fly-by-i i parë më i mirë (planeti b) e kap në orbitë me periudhë
  **≥ ~134 ditë**, pavarësisht drejtimit të hyrjes. Prandaj këmba e parë pas hyrjes zgjat ~130–150 ditë.
- Pra problemi ndahet natyrshëm **këmbë pas këmbe**: çdo fly-by e ul energjinë dhe e shkurton orbitën.

**Algoritmi:** beam search ku çdo këmbë optimizohet si problem 4-D (Nelder–Mead), në dritare kohore.

| Hapi | Hypervolume |
|---|---:|
| Beam i vogël (B = 40) | 5.630M |
| B = 300, të gjithë planetët hyrës | 5.733M |
| + fillime "fizike" (kohët e kalimit të orbitës) | 5.735M |
| + dritare kohore të imta, planeti b i pari | 5.783M |
| Bashkimi i 7 run-eve + polish lokal | 5.792M |

**Zbulim:** fronti nis gjithmonë me b → c. Kjo është shpjegimi fizik i kapjes, jo rastësi.

**Kthime në rënie:** një run 2-orësh shtesë (4 run në 4 bërthama) dha vetëm **+749**
(5,791,586 → 5,792,335). Përfundim: më shumë kohë llogaritjeje nuk mjafton; duhet algoritëm më i mirë.

## Faza 2 – Heuristika si A\* (5.812M → 5.815M)
**Problemi i fazës 1:** beam-i i krahasonte gjendjet vetëm me atë që kishte ndodhur (ΔV dhe koha deri
tani), jo me atë që do të ndodhte. Një gjendje me orbitë të "kapur" mirë do të ketë këmbë të shkurtra
më pas, por beam-i nuk e dinte.

**Ideja:** pas çdo planeti llogaritet P_min, periudha më e shkurtër që mund të japë fly-by-i. Mbi
turat tona, koha e mbetur e turës parashikohet nga P_min me korrelacion **0.88–0.97**. Beam-i
rendit sipas `ΔV + λ·(t + h)`, si funksioni `f = g + h` te A\*.

**Krahasim i drejtë** (të njëjtat cilësime, B = 150, ~5 min):

| | Hypervolume |
|---|---:|
| Pa heuristikë | 5.660M |
| Me heuristikë | **5.772M (+112k)** |

**Run i gjatë:** 20 orë në 4 bërthama, 54 run (beam + rioptimizimi i fundit të turës):
5,795,295 → 5,811,765. Rekordet e reja erdhën shpesh në 8 orët e para, pastaj rrallë: konvergjencë
(grafiku `konvergjenca.png`).

**Gabim i gjetur dhe i rregulluar:** disa tura në kufi me distancën minimale nga ylli kalonin
kontrollin e shpejtë, por jo atë zyrtar (dallime numerike ~1e-13). Kur pastrohej arkivi, ato
"fshihnin" tura të vlefshme. Pas rregullimit (verifikim zyrtar i frontit para pastrimit) dhe
bashkimit të të gjitha run-eve (1.27 milionë tura): **5,814,770**.

**Zgjedhja e 100 pikave:** krahas metodës grykëse, implementuam zgjedhjen **optimale** me programim
dinamik (Bringmann et al., GECCO 2014). Në rastin tonë dhanë të njëjtin rezultat.

## Faza 3 – Eksperimente për vendin e parë (rezultat negativ)
Lexuam prezantimin e ekipit fitues (Space Coders, workshop i ESA-s, 2022) dhe provuam idetë e tyre:

| Ideja | Rezultati | Shpjegimi |
|---|---|---|
| Optimizim lokal me data mbërritjeje (ILS) | 0 | zgjidhjet tona janë në minimume të mprehta: 40 s zhvendosje = +60 m/s |
| Shkurtimi i hyrjes me ΔV (T0 − 60 ditë) | +1k | anija mbërrin më shpejt dhe kapet më keq; këmbët e tjera zgjaten |
| Beam me sekuencë fikse (0-1-5-4-3-2-6) | +2.3k | sekuencat ishin tashmë të eksploruara mirë |
| Kërkim vetëm te "gjuri" (ΔV ≤ 20 m/s) | më keq (1810–1822 ditë kundrejt 1802) | kufizimi i fortë pret rrugët e mira |
| Data e lirë e mbërritjes te planeti i parë | P_min 133.9 → 133.3 ditë | kufiri i kapjes është fizik |
| Polish lokal i gjithë vektorit (SLSQP, CMA-ES) | +2k | ndryshimet e vogla "shpërthejnë" në fund të turës |

"Breaking maneuver"-in (ΔV menjëherë pas fly-by-it të parë), idenë kryesore të fituesve, algoritmi
ynë e kishte gjetur vetë: e gjithë pjesa e frontit nën ~1790 ditë e përdor.

## Përfundimi
- **Rezultati final: 5,814,770** = 99.05% e rekordit të njohur (5,870,436), mbi vendin e 2-të.
- Metoda konvergjon (kthime në rënie të qarta në run-et e gjata).
- Turat pa karburant janë te ~1802–1806 ditë; kufiri fizik i kapjes (~133 ditë) e vendos minimumin
  rreth ~1800 ditë. Vlera ~1773 ditë e fituesve vjen nga një model i thjeshtuar.
- Kufiri teorik global i problemit nuk dihet (problem jokonveks, pa vërtetim optimaliteti).

## Referenca
- ESA Advanced Concepts Team, *SpOC 2022: Trappist Tour* – problemi dhe kodi zyrtar i vlerësimit.
- S. Goulet, V. Debout (CS Group, "Space Coders"), prezantimi i metodës fituese, workshop SpOC, ESA, 2022.
- P. Heine, *Solving the SpOC GECCO Trappist Tour challenge*, Bielefeld, 2022.
- Y. Shimane, K. Iiyama, *Methods for Dual-Objective High-Energy Tour Design*, AAS 23-431, 2023.
- K. Bringmann, T. Friedrich, P. Klitzke, *Two-Dimensional Subset Selection for Hypervolume and
  Epsilon-Indicator*, GECCO 2014.
- D. Izzo, *Revisiting Lambert's problem*, Celestial Mechanics and Dynamical Astronomy, 2015.
- K. Deb et al., *A Fast and Elitist Multiobjective Genetic Algorithm: NSGA-II*, IEEE TEC, 2002.
