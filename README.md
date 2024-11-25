# Magistrska naloga: Detekcija kratkih petidov
To je zbirka skript, ki sem jih uporabila v okviru magistrske naloge na programu Biokemija, UL FKKT. Naloga se osredotoča na detekcijo kratkih petidov in vključuje analizo podatkov različnih programov za masno spektrometrijo.

**Mentor**  
izr. prof. dr. Tomaž Curk

**Avtorica**  
Tanja Gošnjak

## Vsebina repozitorija
Repozitorij vsebuje naslednje datoteke:  
- short_peptide_parser.ipynb - Skripta za pridobitev podatkov o številu identificiranih peptidov po uvedenem ključu (zaporedje, protein(i), vrsta mutacije).  
- results.Rmd - Skripta za vizualizacijo podatkov pridobljenih s short_peptide_parser.ipynb.  
- pepXML_parser.ipynb - Skripta za pridobitev podatkov o modificiranih peptidih iz pepXML datoteke.  
- parser_plots.ipynb - Skripta za primerjavo identificiranih peptidov in vizualizacijo rezultatov z Vennovim diagramom.  
- utils.py - Različne pomožne funkcije za obdelavo podatkov.  
- README.md - Opis projekta.  

## Namen magistrske naloge
Namen magistrske naloge je bil raziskati in primerjati različne prosto dostopne bioinformatske programe za analizo masnih spektrov peptidov, s posebnim poudarkom na identifikaciji peptidov in njihovih modifikacij ter mutacij. Želeli smo oceniti zmogljivosti, uporabnost in jasnost dokumentacije teh programov, da bi raziskovalcem omogočili izbiro najprimernejšega orodja za njihove potrebe. S tem namenom smo uporabili štiri prosto dostopne programe: MSFragger (1), MetaMorpheus (2), AlphaPept (3) in MaxQuant (4), in izvedli primerjalno analizo njihovih rezultatov na podlagi podatkov iz podatkovne zbirke PRIDE. V okviru naloge smo se osredotočili na število identificiranih peptidov in občutljivost na post-translacijske modifikacije ter točkovne mutacije. Zanimala nas je primerljivost navedenih rezultatov med programi in z objavljenimi rezultati.

(1) *Kong AT, Leprevost F V., Avtonomov DM, Mellacheruvu D, Nesvizhskii AI. MSFragger: Ultrafast and comprehensive peptide identification in mass spectrometry-based proteomics. Nat Methods. 2017 Apr 27;14(5):513–20.*  
(2) *Solntsev SK, Shortreed MR, Frey BL, Smith LM. Enhanced Global Post-translational Modification Discovery with MetaMorpheus. J Proteome Res. 2018 May 4;17(5):1844–51.*  
(3) *Strauss MT, Bludau I, Zeng WF, Voytik E, Ammar C, Schessner JP, et al. AlphaPept: a modern and open framework for MS-based proteomics. Nat Commun. 2024 Dec 1;15(1).*  
(4) *Cox J, Mann M. MaxQuant enables high peptide identification rates, individualized p.p.b.-range mass accuracies and proteome-wide protein quantification. Nat Biotechnol. 2008 Dec;26(12):1367–72.*
