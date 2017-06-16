# SV_Pipeline
WGS Structural Variation Pipeline
Software Version:
python/3.6.1
extract_sv_reads/1.1.2
lumpy-sv/v0.2.13
svtyper/v0.1.4

Computational requirements:
hoffman2 4-core 8gb-ram

Run Time:
2-3hr/sample Read Extraction + Lumpy
~10hr/sample svtyper

Example run:

lumpy_submission_shell.sh (edit paths to scratch / list of bams to call (line seperated), all bams should be indexed, calls SGE job array, call lumpy_variant.sh per sample)

sort.py, merge.py part of lumpy package

bnd_remover.py -i input_list -o output_suffix -d directory_path

sv_typer_submission_shell.sh (edit paths to scratch / list of bams to call (line seperated), all bams should be indexed, calls SGE job array, call svtyper_shell.sh per sample)
