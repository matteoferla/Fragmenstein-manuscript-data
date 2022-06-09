# Fragmenstein-figures
Analyses and figures for the Fragmenstein manuscript (to be revealed later).

> :construction: This repository is still under construction.

This folder contains code to make the figures.
It is not intended for greater dissemination.


## Placement of examples

### NUDT7

In the deposition group G_1002045 are several molecules (see [NUDT7 notes](NUDT7/README.md)).

The JACS paper calls `NUOOA000181` `NUDT7-REV-1`, while the ChemMedChem paper calls it `NU181`.
This was merged with `PCM-0102716` to make `NUDT7-COV-1`/`NU443` in the JACS paper.

![NUDT7](NUDT7/NUDT7_NU443.png)
_Comparison of crystallised and placed followup from NUDT7 targets._

The two hits NU181 (in lavender, LHS) and PCM0102716 (puce, LHS) were merged in <paper> yielding the merger NU443.
The crystal structure of NU443 (sky blue, RHS) overlayed with the placement predicted by Fragmenstein (green).
PCM0102716 and NU443 are covalent with Cys73 via a acryloyl warhead._
The RMSD between the placed model and the crystal structure of NU443 is 0.28 Å, 
while the combined RMSD values of the model and the structure against the inspiration hits are 0.65 and 0.61 Å. 
In the placement process the hydroxyl of NU181 was automatically discarded from the mapping 
as it would otherwise impede the mapping of the second hit (PCM0102716) due to the greater proximity of 
the NU181 hydroxyl to the oxygen of the acryloyl warhead of PCM0102716 
rather than to the carbon bonded the benzene ring in PCM0102716.

But there appear to be two chiralities: `NU443` is 2S-pyrrolidine and `NU442` is 2R-pyrrolidine.
The latter, were it reacted would be `'*CC(=O)N1CCC[C@@H]1c1cccc(CC(=O)Nc2cccc(Cl)c2)c1'`.

The latter scores 1.1 kcal/mol higher (–9.2 kcal/mol) than the former (–10.3 kcal/mol)
and the SG atom in NU443 is only .3 Å away from where it is in PCM-0102716,
while in NU442 it is .9 Å away.
but this is far from a smoking gun to explain why it did not react.
The RMSD against the hits went from 0.6 Å to 0.9 Å.

## Color names
The colours at present are the divergent colo(u)r of Fragmenstein as found in `fragmenstein.branding.divergent_colors`.
The primary colour is <span style="color:#AED882;">#AED882, feijoa</span>, but I'll call it green.
As they are divergent colours, they are not quite the same as crayola colours.
I am calling <span style="color:#48DEE9;">#48DEE9, sky blue</span>,
<span style="color:#E9B9FF;">#E9B9FF, mauve</span>,
<span style="color:#FFB7AA;">#FFB7AA, lavender</span>.
Mauve is not a common colour (sounds fancy), but it's far from mushroom, orange or coral.
