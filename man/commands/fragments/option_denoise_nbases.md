`--denoise_nbases` *integer*
: With `--denoise_bysample`, learn the error model from whole samples,
  taken in order of first appearance, until their reads add up to at
  least *integer* bases. The default is 100000000, as DADA2's
  `learnErrors()` (`nbases = 1e8`). Zero means all samples. All samples
  are denoised, whichever are used to learn. Ignored without
  `--denoise_bysample`, or when `--denoise_errin` is used.
