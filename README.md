Travail d'étude et de recherche - M1 de Mathématiques
=====================================================
***Calcul rapide des polynômes cyclotomiques***

**Auteur :** *Jean-Philippe MERX* / 
**Enseignant :** *Pierre-Vincent Koseleff*



# Rappel du thème du TER
Dans [AM10](documents/cyclotomic-highperf.pdf), Arnold et Monagan proposent une méthode pour calculer efficacement et rapidement les polynômes cyclotomiques Φn.
Il s’agit de comprendre les méthodes utilisées et mises en œuvre, ainsi que d’expliquer, dans le cadre du programme du Master de Mathématiques, divers points restés obscurs ou imprécis : propriétés des polynômes cyclotomiques, irréductibilité, taille des coéfficients ; multiplications et divisions rapides, transformation de Fourier rapide.
Les ouvrages [GG13, AECF, SM] pourront utilement être consultés. Le travail peut mener à des calculs explicites et un premier contact avec des logiciels de calcul formel (Sage, Maple, etc).

## Références
- [[AM10]](http://wayback.cecm.sfu.ca/~ada26/cyclotomic/PDFs/highperf.pdf) A. Arnold, M. Monagan, A high-performance algorithm for calculating cyclotomic polynomials. In Proceedings of the 4th International Workshop on Parallel and Symbolic Computation, pp 112–120, 2010.
- [[AECF]](https://hal.archives-ouvertes.fr/AECF/) A. Bostan etal., Algorithmes Efficaces en Calcul Formel, 2017
- [GG13] J. von zur. Gathen, J. Gerhard, Modern Computer Algebra. Cambridge University Press, troisième édition, 2013.
- [[SM]](https://doc.sagemath.org/pdf/en/tutorial/sage_tutorial.pdf) SageMath Reference Manual, 2023


# Structure du [repo](https://github.com/jpmjpmjpm/cyclotom-fast) :
- Répertoire [documents](documents) :
    - [Article à analyser](documents/cyclotomic-highperf.pdf) : document à analyser constituant le thème du TER.
    - [Mémoire TER M1](documents/jpmerx-polynoms-cyclotomics.pdf) : document PDF du mémoire du TER.
    - [SPS4_64.c](documents/SPS4_64.c) : programme C de calcul récursif de l'article initial.
- Répertoire [cyclotomfast](cyclotomfast) :
    - [sps4.py](cyclotomfast/sps4.py) : programmes de calculs des polynômes cyclotomiques.
    - [test_sps4.py](cyclotomfast/test_sps4.py) : programmes de tests.

# Miscellaneous
### [SPS4_64.c](documents/SPS4_64.c) program
Comments on some variables:
- Dg: the degree of the computed product polynomial 