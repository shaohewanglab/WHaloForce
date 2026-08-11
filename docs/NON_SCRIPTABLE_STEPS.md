# Non-scriptable and instrument steps

Some steps in the analysis are done in instrument or GUI software and are not
reproduced by code in this repository. They are documented here. The scripted
pipeline begins from the exports these steps produce.

## FLIM acquisition (Leica LAS X)
Time-correlated single-photon-counting FLIM (TCSPC-FLIM) was acquired on a Leica
STELLARIS built on a DMi8 inverted stand, controlled with Leica Application Suite
X (LAS X v4.5.0.25531), using a 40X 1.25 NA glycerol-immersion objective (Leica
11506422) and a fast-response HyD X2 or HyD X4 detector. Samples were excited
with an 80-MHz pulsed white-light laser; unless otherwise indicated, excitation
was at 561 nm and emission collected from 565 to 650 nm. Acquisition settings
were held constant within each experiment and matched between groups, and were
chosen to reach at least 100 detected photons per pixel at the structures of
interest (50 for exceptionally dim samples).

## FLIM decay fitting (Leica LAS X)
Fluorescence decays were fit in LAS X (v4.8.1.29271 -- a newer version than the
acquisition software above) by reconvolution to a two- or three-exponential
model, as appropriate for the fluorophore and experiment. The fit returns the
component lifetimes (tau_i), amplitudes (alpha_i), and the amplitude- and
intensity-weighted mean lifetimes. These per-image fit results are exported as
spreadsheets, and the raw photon data are exported as PTU files. The scripted
pipeline in this repository starts from those exports.

Intensity-weighted mean lifetime is the quantity plotted throughout, computed as
tau_int = sum_i (a_i * tau_i^2) / sum_i (a_i * tau_i).

Amplitude-weighted mean lifetime, which is proportional to the fluorescence
quantum yield, is also plotted for in vitro characterization (Fig. 1d), computed
as tau_amp = sum_i (a_i * tau_i) / sum_i (a_i).

The instrument response function is recovered from a Rhodamine B reference (2 uM
in ethanol) acquired at 37 C or room temperature with detector and spectral
settings matched to the corresponding experiment. The reference is fitted in LAS
X against the instrument IRF to obtain its mono-exponential lifetime; the
dataset-specific fitted values used here range from 2.088 to 2.196 ns. Which
value goes with which dataset is listed in src/pipeline/README.md.

For region-of-interest measurements, the raw photon decays of all qualifying
pixels in an ROI are summed and fit once by non-negative least squares (NNLS,
scipy.optimize) against an instrument response function recovered from a
Rhodamine B reference. This pooled-decay estimate is the principal ROI value.
That part *is* scripted: see src/pipeline/ (flim_pipeline.py, then flim_fa_fit.py
or flim_roi_fit.py) and the runnable check in examples/flim_roi_fit/.

## Segmentation (Cellpose-SAM 4.1.1)
Cells, focal adhesions, and junctions were segmented from intensity images using
Cellpose-SAM 4.1.1, either the base model or models fine-tuned in the Cellpose GUI.
The fine-tuning itself is a GUI step: hand-corrected crops from these datasets
were used as training data, and the resulting weights are deposited with the raw
data.

Inference itself is scripted (src/pipeline/cellpose_inference.py) and reproduces
well: re-running the published weights on the same image on CPU instead of the
original CUDA GPU recovered all 131 ROIs with median per-ROI IoU 1.000. The label
maps are deposited alongside the raw data only because the model weights are
~1.2 GB each and are impractical to commit here.

What does matter for reproduction is using the weights that go with each image.
The 3T3 focal-adhesion datasets applied a different fine-tune per treatment (DMSO
vs blebbistatin), and substituting one for the other yields plausible but
different adhesions, and hence different lifetimes. The model used for every image
is recorded in the model column of each deposited
cellpose_labels/<dataset>/inference_log.csv, and cellpose_inference.py refuses to
run when its per-condition --model patterns are ambiguous for an image.
examples/segmentation/ ships one mask with its intensity image and validates it;
see that folder's README.

## Structure renders (PyMOL and AlphaFold3)
Structure panels (Fig. 1b, S4b, S4e, S6b) are ray-traced renders made in PyMOL
from AlphaFold3 predictions and RCSB PDB crystal structures. See
docs/pymol_render_notes.md for the minimal command reference and the input
structure files. These are images, not quantitative plots; no source-data CSV
applies.

## Sequence alignment (Fig. S6a)
The vinculin linker alignment was made in Jalview using the T-Coffee method,
colored by the ClustalX scheme. NCBI accession numbers are listed in the figure
legend. This is a GUI step; no code is provided.

## Super-resolution imaging and processing (Fig. 5e)
The single-z super-resolution immunofluorescence images were acquired on a
separate instrument: a home-built line-rescan confocal system, to be fully
described in a forthcoming manuscript. A 60X 1.2 NA water-immersion objective
(Nikon CFI Plan Apo VC 60XC WI) was used with 405/488/561/637 nm laser lines
(Coherent OBIS CellX), four bandpass emission filters (Semrock FF02-447/60,
FF03-525/50, FF01-600/52, FF01-680/42) and an sCMOS camera (Hamamatsu
ORCA-Fusion BT). Pixel size 53.5 nm, axial step 200 nm, 40 ms per channel per
frame; a 32-um stack took 22.3 s.

Images were deconvolved with 30 iterations of Richardson-Lucy in MATLAB. That
deconvolution code is a separate published tool, not part of this repository:
https://github.com/eguomin/regDeconProject

This panel is a representative image, not a quantitative plot; no source-data
CSV applies.

## Single-molecule optical-tweezers force spectroscopy (Fig. 1f-j, S2c, S2d)
Force and confocal-fluorescence data were acquired on a LUMICKS C-Trap (optical
tweezers combined with confocal microscopy) using a five-channel microfluidic
flow cell, and exported as HDF5 (.h5).

Single molecules of cpWHalo-554 or cpHalo-554 were tethered between two 2.14 um
streptavidin beads (Spherotech SVP-20-5) through two dsDNA handles. Each handle
is a 7.5 kb PCR product carrying a dual-biotin at one end and an ALFA-tag
nanobody at the other (coupled via an azide-functionalized primer to DBCO-ALFA-Nb,
NanoTag N1505-DBCO), giving 15 kb of DNA in the assembled tether. Tethers whose
extension matched that 15 kb expectation under the extensible worm-like chain
model were selected for correlated force and fluorescence measurement. Imaging
used 575 nm excitation at ~5 uW, 100 nm pixels and 0.1-0.5 ms per pixel, in an
oxygen-scavenging buffer (4 mM Trolox, 2.5 mM protocatechuic acid, 0.2% v/v rPCO
in 50 mM HEPES / 250 mM NaCl, pH 7.2).

Kymograph intensity traces were extracted from the .h5 files upstream of this
repository, using scripts based on the LUMICKS Pylake API. This repository plots
the extracted force and intensity traces (the committed per-trace CSVs); the
dwell-time analysis of the quenched state is in
src/fig1/QU_transition_analysis.py. Full flow-cell passivation, force-calibration
and handle-conjugation protocols are in the manuscript methods.
