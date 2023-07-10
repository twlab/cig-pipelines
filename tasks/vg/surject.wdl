version development

import "../../structs/runenv.wdl"

task run_surject {
  input {
     File gam
     String sample
     String library
     File gbz
     RunEnv runenv
  }

  String output_bam = sub(basename(gam), ".gam$", ".bam")
  # options:
  # -x, --xg-name FILE       use this graph or xg index (required)
  # -t, --threads N          number of threads to use
  # -p, --into-path NAME     surject into this path or its subpaths (many allowed, default: reference, then non-alt generic)
  # -F, --into-paths FILE    surject into path names listed in HTSlib sequence dictionary or path list FILE
  # -i, --interleaved        GAM is interleaved paired-ended, so when outputting HTS formats, pair reads
  # -M, --multimap           include secondary alignments to all overlapping paths instead of just primary
  # -G, --gaf-input          input file is GAF instead of GAM
  # -m, --gamp-input         input file is GAMP instead of GAM
  # -c, --cram-output        write CRAM to stdout
  # -b, --bam-output         write BAM to stdout
  # -s, --sam-output         write SAM to stdout
  # -l, --subpath-local      let the multipath mapping surjection produce local (rather than global) alignments
  # -P, --prune-low-cplx     prune short and low complexity anchors during realignment
  # -a, --max-anchors N      use no more than N anchors per target path (default: 200)
  # -S, --spliced            interpret long deletions against paths as spliced alignments
  # -A, --qual-adj           adjust scoring for base qualities, if they are available
  # -N, --sample NAME        set this sample name for all reads
  # -R, --read-group NAME    set this read group for all reads
  # -f, --max-frag-len N     reads with fragment lengths greater than N will not be marked properly paired in SAM/BAM/CRAM
  # -L, --list-all-paths     annotate SAM records with a list of all attempted re-alignments to paths in SS tag
  # -C, --compression N      level for compression [0-9]
  # -V, --no-validate        skip checking whether alignments plausibly are against the provided graph
  # -w, --watchdog-timeout N warn when reads take more than the given number of seconds to surject
  command <<<
    vg surject -t ~{runenv.cpu - 1} -b -N ~{sample} -R ~{library} -x ~{gbz} ~{gam} > ~{output_bam}
  >>>

  output {
    File bam = glob("*.bam")[0]
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: "${runenv.memory} GB"
  }
}
