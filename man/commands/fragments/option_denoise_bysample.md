`--denoise_bysample`
: Treat the input file as several samples, identified by the `sample=`
  (or `barcodelabel=`) annotation of each read header, or, when neither
  is present, by the start of the header up to its first character
  that is not a letter, a digit or an underscore (as
  `--otutabout` does). Each sample is dereplicated and denoised on its
  own. The error model is learnt from the samples as DADA2's
  `learnErrors()` learns it: at each self-consistency round, every
  sample is denoised on its own, the observed transitions of all
  samples are added up, and the next model is fitted to the sum (see
  `--denoise_nbases`). With `--denoise_errin`, the given model is
  applied to each sample. Denoised sequences are written sample by
  sample, in order of first appearance, under the header of their
  first read, which carries the sample annotation. With several
  threads, samples are denoised in parallel (see `--threads`).
