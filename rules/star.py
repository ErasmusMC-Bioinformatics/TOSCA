import os
import tempfile
from snakemake.shell import shell

sorted = snakemake.params.get("sorted", "")
twopass = snakemake.params.get("twopass", "")
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

if "--outSAMtype BAM SortedByCoordinate" in sorted:
    stdout = "BAM_SortedByCoordinate"

with tempfile.TemporaryDirectory() as tmpdir:
    shell(
        "STAR "
        " --runThreadN {snakemake.threads}"
        " --genomeDir {snakemake.input.idx}"
        " --readFilesIn {snakemake.input.fq1}"
        " --readFilesCommand gunzip -c"
        " {twopass}"
        " {sorted}"
        " --outTmpDir {tmpdir}/STARtmp"
        " --outFileNamePrefix {tmpdir}/"
        " --outStd {stdout}"
        " > {snakemake.output.aln}"
        " {log}")
    shell("cat {tmpdir}/SJ.out.tab > {snakemake.output.sj:q}")
    shell("cat {tmpdir}/Log.out > {snakemake.output.log:q}")