# Procesi i punës: nga 2.4M në 5.815M

Ky dokument përshkruan rrugëtimin e projektit me radhë: çfarë provuam, çfarë funksionoi, çfarë jo,
dhe pse. Kodi i fazave të hershme dhe i eksperimenteve nuk është në repo; këtu ruhen vetëm vlerat
dhe përfundimet. Në repo është vetëm versioni final.

## Përmbledhje

| Faza | Ideja kryesore | Hypervolume |
|---|---|---:|
| 0 | Algoritme gjenetike (NSGA-II / R-NSGA-II) mbi gjithë vektorin 34-D | ~2.4M |
| 1 | Analiza fizike + beam search këmbë pas këmbe | 5.792M |
| 2 | Heuristikë për të ardhmen (si A\*) + rioptimizimi i fundit të turës | 5.812M |
| 2b | Bashkimi i të gjitha run-eve, me verifikim zyrtar dhe zgjedhje optimale të 100 pikave | **5.815M** |
| 3 | Eksperimente për vendin e parë (rezultat negativ) | – |

Për krahasim, renditja e garës zyrtare: 5.870M (vendi i 1-rë), ~5.73M (i 2-ti), ~5.60M (i 3-ti).

## Faza 0 – Pikënisja: algoritme gjenetike (~2.4M)
- NSGA-II dhe R-NSGA-II mbi të 34 variablat njëkohësisht.
- Rezultati: zgjidhjet më të mira kishin ~440 m/s karburant te ~1850 ditë. Score ~2.4M.
- **Mësimi:** problemi është shumë "kaotik" për kërkim të verbër. Një ndryshim i vogël në fillim
  të turës e ndryshon krejtësisht fundin, ndaj kërkimi gjenetik humbet në hapësirën 34-D.

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
