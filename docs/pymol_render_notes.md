# PyMOL render notes for the structure panels

Structure panels are ray-traced renders, not scripted plots. This is the minimal
command reference used to make them. Input structures (AlphaFold3 predictions and
RCSB crystal structures) and the saved PyMOL sessions (.pse) are archived with
the raw data (see RAW_DATA.md). Panels: Fig. 1b (cpWHalo vs HaloTag7-JF669, PDB
8SW8), Fig. S4b (Ecad-WF), Fig. S4e (aCat-WF), Fig. S6b (mouse vinculin, worm DEB-1).

## pLDDT coloring (AlphaFold confidence is stored in the B-factor column)
```pml
set_color af_very_low, [1.0, 0.49, 0.27]
set_color af_low,      [1.0, 0.86, 0.30]
set_color af_conf,     [0.40, 0.76, 1.0]
set_color af_high,     [0.00, 0.33, 0.80]
color af_very_low, (af and b < 50)
color af_low,      (af and b > 49.9) and (af and b < 70)
color af_conf,     (af and b > 69.9) and (af and b < 90)
color af_high,     (af and b > 89.9)
```

## Overlay AlphaFold on crystal, style crystal as transparent reference
```pml
cealign xtal, af
hide everything, xtal and solvent
show cartoon, xtal
color gray70, xtal
set cartoon_transparency, 0.5, xtal
```

## Ligand in the 8SW8 crystal (Fig. 1b)
The dye ligand residue name in 8SW8 is XSR.
```pml
color magenta, xtal and resn XSR
show sticks, xtal and resn XSR
set stick_radius, 0.2, xtal and resn XSR
```

## Domain colors and residue ranges (Okabe-Ito palette)
```pml
# alpha-catenin (aCat-WF, chain A; same numbering as human 4IGG chain B)
set_color DD,    [0x00, 0x9E, 0x73]
set_color VBS,   [0xCC, 0x79, 0xA7]
set_color Mfrag, [0x56, 0xB4, 0xE9]
set_color FABD,  [0xE6, 0x9F, 0x00]
set_color WF2,   [0x33, 0x22, 0x88]
color DD,    aCat-WF and chain A and resi 82-259
color VBS,   aCat-WF and chain A and resi 277-393
color Mfrag, aCat-WF and chain A and resi 399-630
color WF2,   aCat-WF and chain A and resi 655-978
color FABD,  aCat-WF and chain A and resi 992-1166

# E-cadherin (Ecad-WF, chain A)
set_color EC,  [0x56, 0xB4, 0xE9]
set_color JMD, [0x00, 0x9E, 0x73]
set_color CBD, [0xE6, 0x9F, 0x00]
color EC,     Ecad-WF and chain A and resi 157-697
color black,  Ecad-WF and chain A and resi 709-733
color JMD,    Ecad-WF and chain A and resi 758-775
color WF2,    Ecad-WF and chain A and resi 781-1104
color CBD,    Ecad-WF and chain A and resi 1156-1185

# mouse vinculin (mmVCL, chain A)
set_color Head,   [0x00, 0x9E, 0x73]
set_color Linker, [0x93, 0x95, 0x98]
set_color Tail,   [0xE6, 0x9F, 0x00]
set_color Pro,    [0x56, 0xB4, 0xE9]
color Head,   mmVCL and chain A and resi 1-835
color Linker, mmVCL and chain A and resi 836-895
color Tail,   mmVCL and chain A and resi 896-1066
color Pro,    mmVCL and chain A and resi 840-847
color Pro,    mmVCL and chain A and resi 860-878

# worm vinculin (DEB-1d, chain A)
color Head,   DEB-1d and chain A and resi 1-810
color Linker, DEB-1d and chain A and resi 811-910
color Tail,   DEB-1d and chain A and resi 911-1090
color Pro,    DEB-1d and chain A and resi 831-837
color Pro,    DEB-1d and chain A and resi 866-887
```

## Render
```pml
bg_color white
set ray_shadows, off
set antialias, 2
ray 4000, 3000
png panel.png, dpi=300
```

Notes: avoid `>=`/`<=` in selections (use numeric offsets as above); always `ray`
before exporting; save the session with `save panel.pse` so a render is reproducible.
