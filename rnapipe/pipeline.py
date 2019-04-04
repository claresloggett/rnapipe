'''
Build the pipeline workflow by plumbing the stages together.
'''

from ruffus import Pipeline, suffix, formatter, add_inputs, output_from
from rnapipe.stages import PipelineStages
from rnapipe.pipeline_base.utils import safe_make_dir
from rnapipe.experiment import Experiment
from os import path
from rnapipe.utils import *
from rnapipe.constants import *
import logging
import re
import sys

logging.basicConfig(format='[%(levelname)s] %(message)s', level=logging.DEBUG)

def make_pipeline(state):
    '''Build the pipeline by constructing stages and connecting them together'''
    # Build an empty pipeline
    pipeline = Pipeline(name="rnapipe")

    # Get the details of the experiment (samples, config, inputs, ...)
    experiment = Experiment(state)

    # Get reference file locations
    reference_genome = state.config.get_options("reference_genome")
    gene_ref = state.config.get_options("gene_ref")

    # Print out samples
    sample_text = [s.info() for s in experiment.sample_list]
    logging.info("Analysis samples:\n{}".format("\n".join(sample_text)))

    # Stages are dependent on the state. Experiment object is also passed so
    # we can access metadata later.
    stages = PipelineStages(state, experiment=experiment)

    # Make directories
    output_dir = get_output_paths(
            results_dir=state.config.get_options("results_dir"),
            default_paths=OUTPUT_PATHS)
    make_output_dirs(output_dir)
    logging.debug(output_dir)

    # The original FASTQ files
    # This is a dummy stage. It is useful because it makes a node in the
    # pipeline graph, and gives the pipeline an obvious starting point.
    pipeline.originate(
        task_func=stages.do_nothing,
        name="original_fastqs",
        output=experiment.R1_files)

    # Create reference index for alignment
    if not experiment.index_provided:
        pipeline.originate(
            task_func=stages.do_nothing,
            name="reference_genome",
            output=reference_genome)

        if experiment.alignment_method == "star":
            # Create reference index for STAR
            pipeline.transform(
                task_func=stages.create_star_index,
                name="create_star_index",
                input=output_from("reference_genome"),
                filter=formatter(".*"),
                add_inputs=add_inputs(gene_ref),
                output=path_list_join(output_dir["star_index"],
                           ["SA", "Genome", "genomeParameters.txt"]),
                extras=[output_dir["star_index"]])

        elif experiment.alignment_method == "hisat2":
            # Create reference index for HISAT2
            hisat_basename = path.join(output_dir["hisat_index"], "genome")
            pipeline.transform(
                task_func=stages.create_hisat_index,
                name="create_hisat_index",
                input=output_from("reference_genome"),
                filter=formatter(".*"),
                add_inputs=add_inputs(gene_ref),
                output=path_list_join(output_dir["hisat_index"],
                           ["genome.1.ht2", "genome.2.ht2"]),
                extras=[hisat_basename])
    else:
        # Don't create index if index is supplied
        if experiment.alignment_method == "star":
            output_dir["star_index"] = state.config.get_options("star_index")
            pipeline.originate(
                task_func=stages.do_nothing,
                name="create_star_index",
                output=path_list_join(output_dir["star_index"],
                           ["SA", "Genome", "genomeParameters.txt"]))
        elif experiment.alignment_method == "hisat2":
            hisat_basename = state.config.get_options("hisat_index")
            output_dir["hisat_index"] = path.dirname(hisat_basename)
            prefix = path.basename(hisat_basename)
            pipeline.originate(
                task_func=stages.do_nothing,
                name="create_hisat_index",
                output=path_list_join(output_dir["hisat_index"],
                           ["{prefix}.1.ht2".format(prefix=prefix),
                            "{prefix}.2.ht2".format(prefix=prefix)]))

    # Pre-trim FastQC
    if experiment.paired_end:
        pipeline.transform(
            task_func=stages.fastqc,
            name="fastqc",
            input=output_from("original_fastqs"),
            filter=formatter(".+/(?P<sample>[a-zA-Z0-9-_]+)_R1.fastq.gz"),
            add_inputs=add_inputs("{path[0]}/{sample[0]}_R2.fastq.gz"),
            output=path_list_join(output_dir["fastqc"],
                       ["{sample[0]}_R1_fastqc.zip",
                        "{sample[0]}_R2_fastqc.zip"]),
            extras=[output_dir["fastqc"]])
    else:
        pipeline.transform(
            task_func=stages.fastqc,
            name="fastqc",
            input=output_from("original_fastqs"),
            filter=suffix(".fastq.gz"),
            output="_fastqc.zip",
            output_dir=output_dir["fastqc"],
            extras=[output_dir["fastqc"]])

    # Trimmomatic
    if experiment.trim_reads and experiment.paired_end:
        pipeline.transform(
            task_func=stages.trim_reads,
            name="trim_reads",
            input=output_from("original_fastqs"),
            # Get R1 file and the corresponding R2 file
            filter=formatter(".+/(?P<sample>[a-zA-Z0-9-_]+)_R1.fastq.gz"),
            add_inputs=add_inputs("{path[0]}/{sample[0]}_R2.fastq.gz"),
            output=path_list_join(output_dir["seq"],
                       ["{sample[0]}_R1.trimmed.fastq.gz",
                        "{sample[0]}_R2.trimmed.fastq.gz"]),
            extras=path_list_join(output_dir["seq"],
                       ["{sample[0]}_R1.unpaired.fastq.gz",
                        "{sample[0]}_R2.unpaired.fastq.gz"]))
    elif experiment.trim_reads:
        pipeline.transform(
            task_func=stages.trim_reads,
            name="trim_reads",
            input=output_from("original_fastqs"),
            filter=formatter(".+/(?P<sample>[a-zA-Z0-9-_]+)_R1.fastq.gz"),
            output=path.join(output_dir["seq"],
                             "{sample[0]}_R1.trimmed.fastq.gz"))

    # Post-trim FastQC
    if experiment.paired_end and experiment.trim_reads:
        pipeline.transform(
            task_func=stages.fastqc,
            name="post_trim_fastqc",
            input=output_from("trim_reads"),
            filter=formatter(".+/(?P<sample>[a-zA-Z0-9-_]+)_R1.trimmed.fastq.gz"),
            output=path_list_join(output_dir["post_trim_fastqc"],
                       ["{sample[0]}_R1.trimmed_fastqc.gz",
                        "{sample[0]}_R2.trimmed_fastqc.gz"]),
            extras=["results/qc/post_trim_fastqc/"])
    elif experiment.trim_reads:
        pipeline.transform(
            task_func=stages.fastqc,
            name="post_trim_fastqc",
            input=output_from("trim_reads"),
            filter=suffix(".trimmed.fastq.gz"),
            output=".trimmed_fastqc.gz",
            output_dir=output_dir["post_trim_fastqc"],
            extras=[output_dir["post_trim_fastqc"]])

    # If there are technical replicates, each is mapped independently.
    # This is so each technical replicate maintains a separate read group.
    if experiment.alignment_method == "star":
        align_task_name = "star_align"
        if experiment.trim_reads:
            (pipeline.transform(
                task_func=stages.star_align,
                name=align_task_name,
                input=output_from("trim_reads"),
                filter=formatter(".+/(?P<sample>[a-zA-Z0-9-_]+)" \
                                 "_R[12](.trimmed)?.fastq.gz"),
                output="%s/{sample[0]}/{sample[0]}.star.Aligned.out.bam" \
                        % output_dir["alignments"],
                extras=[output_dir["star_index"], "{sample[0]}"])
            ).follows("create_star_index")
        else:
            (pipeline.transform(
                task_func=stages.star_align,
                name=align_task_name,
                input=output_from("original_fastqs"),
                filter=formatter(".+/(?P<sample>[a-zA-Z0-9-_]+)" \
                                 "_R[12](.trimmed)?.fastq.gz"),
                output="%s/{sample[0]}/{sample[0]}.star.Aligned.out.bam" \
                        % output_dir["alignments"],
                extras=[output_dir["star_index"], "{sample[0]}"])
            ).follows("create_star_index")

    if experiment.alignment_method == "hisat2":
        align_task_name = "hisat_align"
        if experiment.trim_reads:
            (pipeline.transform(
                task_func=stages.hisat_align,
                name="hisat_align",
                input=output_from("trim_reads"),
                filter=formatter(".+/(?P<sample>[a-zA-Z0-9-_]+)" \
                                 "_R[12](.trimmed)?.fastq.gz"),
                output="%s/{sample[0]}/{sample[0]}.hisat2.bam" \
                        % output_dir["alignments"],
                extras=[hisat_basename, "{sample[0]}"])
            ).follows("create_hisat_index")
        else:
            (pipeline.transform(
                task_func=stages.hisat_align,
                name="hisat_align",
                input=output_from("original_fastqs"),
                filter=formatter(".+/(?P<sample>[a-zA-Z0-9-_]+)" \
                                 "_R[12](.trimmed)?.fastq.gz"),
                output="%s/{sample[0]}/{sample[0]}.hisat2.bam" \
                        % output_dir["alignments"],
                extras=[hisat_basename, "{sample[0]}"])
            ).follows("create_hisat_index")

    # Sort BAM by coordinates
    pipeline.transform(
        task_func=stages.sort_bam_by_coordinate,
        name="sort_bam_by_coordinate",
        input=output_from(align_task_name),
        filter=formatter(".+/(?P<sample>[a-zA-Z0-9-_]+)\.(?P<method>(star|hisat2))\..*bam"),
        output=["{path[0]}/{sample[0]}.{method[0]}.sorted.bam",
                "{path[0]}/{sample[0]}.{method[0]}.sorted.bam.bai"])

    # Merge files with the same sample name
    if experiment.multiple_technical_replicates:
        pipeline.collate(
            task_func=stages.merge_bams,
            name="merge_bams",
            input=output_from("sort_bam_by_coordinate"),
            filter=formatter(".+/(SM_)?(?P<sm>[a-zA-Z0-9-]+)[^.]*\.(?P<method>(star|hisat2)).sorted.bam"),
            output=path_list_join(output_dir["alignments"],
                ["{sm[0]}.{method[0]}.bam", "{sm[0]}.{method[0]}.bam.bai"]))
    else:
        pipeline.transform(
            task_func=stages.create_symlinks,
            name="merge_bams",
            input=output_from("sort_bam_by_coordinate"),
            filter=formatter(".+/(SM_)?(?P<sm>[a-zA-Z0-9-]+)[^.]*\.(?P<method>(star|hisat2)).sorted.bam"),
            output=path_list_join(output_dir["alignments"],
                ["{sm[0]}.{method[0]}.bam", "{sm[0]}.{method[0]}.bam.bai"]))

    # Sort BAM by name for counting features
    pipeline.transform(
        task_func=stages.sort_bam_by_name,
        name="sort_bam_by_name",
        input=output_from("merge_bams"),
        filter=suffix(".bam"),
        output=".nameSorted.bam")

    # Count features with HTSeq-count
    pipeline.transform(
        task_func=stages.htseq_count,
        name="htseq_count",
        input=output_from("sort_bam_by_name"),
        filter=suffix(".nameSorted.bam"),
        output_dir=output_dir["counts"],
        output=".htseq.txt")

    # Count features with featureCounts
    pipeline.transform(
        task_func=stages.featurecounts,
        name="featurecounts",
        input=output_from("sort_bam_by_name"),
        filter=suffix(".nameSorted.bam"),
        output_dir=output_dir["counts"],
        output=".featureCounts.txt")

    # TODO: add multiqc step

#     # Stringtie assembly
#     pipeline.transform(
#         task_func=stages.stringtie_assembly,
#         name="stringtie_assembly",
#         input=output_from("merge_bams"),
#         filter=suffix(".bam"),
#         output_dir=output_dir["stringtie_assembly"],
#         output=".gtf")

    # Stringtie estimates
    pipeline.transform(
        task_func=stages.stringtie_estimates,
        name="stringtie_estimates",
        input=output_from("merge_bams"),
        filter=formatter(".+/(?P<sm>[a-zA-Z0-9-]+)\.(?P<method>(star|hisat2)).bam"),
        output=path_list_join(output_dir["stringtie_estimates"],
                              ["{sm[0]}/{sm[0]}.gtf", "{sm[0]}/e_data.ctab"]))

    # Stringtie counts
    pipeline.collate(
        task_func=stages.stringtie_prepDE,
        name="stringtie_prepDE",
        input=output_from("stringtie_estimates"),
        filter=formatter(".+\.gtf"),
        output=path_list_join(output_dir["stringtie_estimates"],
                              ["gene_count_matrix.csv", "transcript_count_matrix.csv"]))
    return pipeline
