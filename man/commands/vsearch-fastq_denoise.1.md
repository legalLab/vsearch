% vsearch-fastq_denoise(1) version 2.32.0 | vsearch manual
% Torbjørn Rognes, Tomás Flouri, and Frédéric Mahé
#(./fragments/date.md)

# NAME

vsearch \-\-fastq_denoise --- correct sequencing errors in amplicon reads with an error model learnt from the reads


# SYNOPSIS

| **vsearch** **\-\-fastq_denoise** _fastqfile_ (**\-\-fastqout** | **\-\-fastaout** | **\-\-denoise_errout**) _filename_ \[_options_]


# DESCRIPTION

The vsearch command `--fastq_denoise` corrects substitution errors
(and, optionally, indel errors) in amplicon reads, using the divisive
partitioning algorithm of DADA2 (Callahan et al. 2016). Once corrected,
each distinct sequence is an inferred biological sequence, and
clustering by similarity is no longer needed to absorb sequencing
noise.

The reads are dereplicated, and the mean quality score at each position
is recorded for each unique sequence. All unique sequences start in a
single partition, whose center is the most abundant one. Each sequence
is aligned with the center, and the probability *lambda* that a read of
the center is observed as that sequence is computed as the product,
over the positions of the sequence, of the error rates given by the
error model for the facing nucleotides and the quality score at that
position. If the partition holds *n* reads, the number of reads of the
sequence expected from errors alone follows a Poisson distribution of
mean *n* x *lambda*. The sequence with the smallest abundance p-value
(the probability of observing at least as many reads, given that at
least one was observed) becomes the center of a new partition if that
p-value, multiplied by the number of unique sequences, is below
`--denoise_omega_a`. All sequences are then compared with the new
center, and each joins the partition for which its expected number of
reads is the highest. The process is repeated until no sequence is
significant. Each read is finally corrected to the center of its
partition. Singletons are never significant, as their abundance
carries no evidence.

The error model is a table of error rates for each of the sixteen
transitions (true nucleotide to observed nucleotide) and each quality
score. It is learnt from the input reads: starting from a model in
which all rates are one and a single partition (so that every
difference with the most abundant sequence is counted as an error),
partitioning and estimation of the rates from the reads of each
partition, compared with its center (including, as in DADA2, the reads
that are not corrected), are alternated until the model repeats itself, or for at most
`--denoise_maxconsist` rounds. Rates are smoothed over quality scores
by a weighted local regression (loess, span 0.75, degree 2) of the
log10 of the observed rates, and bounded to [1e-7, 0.25]. The model can
be saved with `--denoise_errout`, and a model can be given with
`--denoise_errin`, in which case nothing is learnt. Since qualities
are averaged within unique sequences, the rates are not estimates of
the per-base error rates of the instrument, and should not be used as
such.

By default the model knows substitutions only, as in DADA2, and
sequences that differ from the center of their partition by
insertions or deletions only are merged into it. With
`--denoise_indels model`, indels are errors with learnt rates too, and
such sequences are retained when their abundance is not explained by
indel errors. This is recommended for markers in which indels are
diagnostic, such as 12S and 16S rRNA.

Input reads should have been quality-filtered and truncated (see
[`vsearch-fastx_filter(1)`](./vsearch-fastx_filter.1.md)), and primers
should have been removed. Reads with ambiguous nucleotides are not
corrected. Paired-end reads should be denoised before merging, forward
and reverse reads separately, as merging recomputes the quality
scores on which the error model depends; pairs can then be
re-synchronized and merged (see EXAMPLES). The input file is read
twice, and cannot be a pipe.

By default all the reads of the input file are pooled: they are
denoised together, and the error model is learnt from all of them.
With `--denoise_bysample`, the input file holds several samples,
identified by the `sample=` annotation of the read headers: each
sample is denoised on its own, and the error model is learnt from the
samples as DADA2's `learnErrors()` learns it, adding up at each round
the transitions observed in each sample denoised on its own. This is
the recommended way to denoise a sequencing run: the model draws on all
samples, and its cost grows with the sum, over samples, of their
unique sequences multiplied by their partitions, rather than with the
much larger product for the pooled run. The samples are denoised side
by side, one per thread, the largest first; with fewer samples than
threads, each sample spreads its comparisons over an equal share of
the threads. The results do not depend on the number of threads. A
self-consistency round lasts at least as long as its largest sample
takes on one thread, and threads that share a CPU core (simultaneous
multithreading) slow each other down: if a run is slower with all
available threads, try `--threads` set to the number of physical
cores.

Computation time grows with the number of unique sequences multiplied
by the number of partitions. Memory usage is about 2 kB per unique
sequence, plus 8 bytes per read, plus 8 bytes per read of the most
abundant unique sequence.


# OPTIONS

## mandatory options

At least one of `--fastqout`, `--fastaout` and `--denoise_errout` must
be specified.

#(./fragments/option_denoise_errout.md)

#(./fragments/option_fastaout_denoise.md)

#(./fragments/option_fastqout_denoise.md)


## core options

#(./fragments/option_denoise_bysample.md)

#(./fragments/option_denoise_errin.md)

#(./fragments/option_denoise_indels.md)

#(./fragments/option_denoise_maxconsist.md)

#(./fragments/option_denoise_nbases.md)

#(./fragments/option_denoise_omega_a.md)

#(./fragments/option_denoise_omega_c.md)

#(./fragments/option_fastq_ascii.md)

#(./fragments/option_fastq_qmax.md)

#(./fragments/option_fastq_qmin.md)

#(./fragments/option_fastqout_discarded_denoise.md)

#(./fragments/option_sizeout_denoise.md)


## secondary options

#(./fragments/option_bzip2_decompress.md)

#(./fragments/option_fasta_width.md)

#(./fragments/option_gzip_decompress.md)

#(./fragments/option_label_suffix.md)

#(./fragments/option_lengthout.md)

#(./fragments/option_log.md)

#(./fragments/option_no_progress.md)

#(./fragments/option_quiet.md)

#(./fragments/option_relabel.md)

#(./fragments/option_relabel_keep.md)

#(./fragments/option_relabel_md5.md)

#(./fragments/option_relabel_self.md)

#(./fragments/option_relabel_sha1.md)

#(./fragments/option_threads.md)

#(./fragments/option_xee.md)

#(./fragments/option_xlength.md)

#(./fragments/option_xsize.md)


# EXAMPLES

Correct the reads of a fastq file, treating indels as potential
variants, and write the corrected reads, the denoised sequences with
their abundances, and the error model:

```sh
vsearch \
    --fastq_denoise R1.filtered.fastq \
    --denoise_indels model \
    --fastqout R1.denoised.fastq \
    --fastaout R1.sequences.fasta \
    --sizeout \
    --denoise_errout R1.errors.tsv
```

Denoise forward and reverse reads separately, restore the pairing
(reads that were not corrected are missing from one file or the
other), then merge, allowing no mismatch in the overlap since both
reads are now assumed to be error-free:

```sh
vsearch --fastq_denoise R1.filtered.fastq --fastqout R1.denoised.fastq
vsearch --fastq_denoise R2.filtered.fastq --fastqout R2.denoised.fastq

vsearch \
    --fastx_syncpairs R1.denoised.fastq \
    --reverse R2.denoised.fastq \
    --fastqout R1.synced.fastq \
    --fastqout_rev R2.synced.fastq

vsearch \
    --fastq_mergepairs R1.synced.fastq \
    --reverse R2.synced.fastq \
    --fastq_maxdiffs 0 \
    --fastqout merged.fastq
```

Denoise all the samples of a run in one call: label the reads of each
sample, pool them in one file, then learn one error model from all
samples (as DADA2's `learnErrors()` does) and denoise each sample on
its own with it:

```sh
for sample in sample1 sample2 sample3 ; do
    vsearch \
        --fastq_filter ${sample}.R1.filtered.fastq \
        --sample ${sample} \
        --fastqout - ;
done > run.R1.fastq

vsearch \
    --fastq_denoise run.R1.fastq \
    --denoise_bysample \
    --fastaout run.R1.sequences.fasta \
    --sizeout \
    --denoise_errout run.R1.errors.tsv
```

Learn the error model on a subsample of a sequencing run, then apply it
to a sample of that run:

```sh
vsearch \
    --fastx_subsample run.R1.fastq \
    --sample_size 500000 \
    --fastqout subsample.R1.fastq

vsearch \
    --fastq_denoise subsample.R1.fastq \
    --denoise_errout run.R1.errors.tsv

vsearch \
    --fastq_denoise sample1.R1.fastq \
    --denoise_errin run.R1.errors.tsv \
    --fastqout sample1.R1.denoised.fastq
```


# SEE ALSO

[`vsearch-cluster_unoise(1)`](./vsearch-cluster_unoise.1.md),
[`vsearch-fastq_mergepairs(1)`](./vsearch-fastq_mergepairs.1.md),
[`vsearch-fastx_filter(1)`](./vsearch-fastx_filter.1.md),
[`vsearch-fastx_subsample(1)`](./vsearch-fastx_subsample.1.md),
[`vsearch-fastx_syncpairs(1)`](./vsearch-fastx_syncpairs.1.md),
[`vsearch-fastx_uniques(1)`](./vsearch-fastx_uniques.1.md),
[`vsearch-fastq(5)`](../formats/vsearch-fastq.5.md)


# CITATION

Callahan BJ, McMurdie PJ, Rosen MJ, Han AW, Johnson AJA, Holmes SP
(2016) DADA2: High-resolution sample inference from Illumina amplicon
data. *Nature Methods*, 13, 581-583. doi:
[10.1038/nmeth.3869](https://doi.org/10.1038/nmeth.3869)

Rosen MJ, Callahan BJ, Fisher DS, Holmes SP (2012) Denoising PCR-amplified
metagenome data. *BMC Bioinformatics*, 13, 283. doi:
[10.1186/1471-2105-13-283](https://doi.org/10.1186/1471-2105-13-283)


#(./fragments/footer.md)
