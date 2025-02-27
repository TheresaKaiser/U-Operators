# $U$-Operators
This repository contains the calculations explained in Section 5 of the Article " $U$-Operators Acting on Harmonic Cocycles for $\mathrm{GL}_3$ and Their Slopes" (Gebhard Böckle, Peter Mathias Gräf, Theresa Kaiser). 
It studies the action of the two natural $U$-operators acting on $\Gamma$-invariant spaces of harmonic cocycles for $\mathrm{GL}_3$ for certain congruence subgroups $\Gamma$, in a positive characteristic setting. The cocycle spaces we consider are conjecturally isomorphic to spaces of Drinfeld cusp forms of rank $3$ and level $\Gamma$ via an analogue of Teitelbaum's residue map. We give explicit descriptions of the spaces of harmonic cocycles as subspaces of the vector space of coefficients, and of the resulting $U$- and Hecke operators acting on these. We then explain how to implement these formulas in a computer algebra system. Using the resulting data of slopes (and characteristic polynomials) for the Hecke actions, we observe several patterns and interesting phenomena present in our slope tables. This appears to be the first such study in a $\mathrm{GL}_3$ setting.

The preprint can be found on arXiv.org and theoretical explanations as well as pseudocode are given there. The calculations were realized in the computer algebra system Magma (https://magma.maths.usyd.edu.au/magma/).

## Program files
The calculations are carried out in `main.mg`. **Caution!** The code is not optimized and takes several days to run. We therefore recommend to start the program with the command `nohup magma.mg &`. 
The program creates `.txt` files as output that are meant to be processed further before being read by humans. 

For the output files containing slopes, the first improvement we suggest is to typeset fractions with the LaTeX command `\frac`. In order to do this, open the file in an editor that is capable of search- and replace-actions using Regular Expressions, such as Geany or Notepad++. The syntax for the replace field might vary slightly.
- Search for `([0-9]+)/([0-9]+)` and replace it with `\\frac{\1}{\2}`.\* 

For the output files containing factorized characteristic polynomials, carry out the following. The actions where RegEx Search should be activated are marked with a star: 
- Search for `<` and replace it with `(`. 
- Search for `, ([0-9]+)>` and replace it with `\)^\{\1\}`.\*
- Search for `(X)` and replace it with `X`.
- Search for `,\n` and replace it with ` `. \* Ensure multiline-matching is activated.
- Search for `[\n` and replace it with `$`. \* Ensure multiline-matching is activated.
- Search for `\n]` and replace it with `$`. \* Ensure multiline-matching is activated.

After these steps are completed, the file can be processed with the utility program `toTable.mg`. This in turn creates a new file whose content can be directly copied into a LaTeX table. In `toTable.mg`, the user can choose which columns should appear in the final table.

## Result files
We will upload files containing longer result tables than in the preprint shortly.
