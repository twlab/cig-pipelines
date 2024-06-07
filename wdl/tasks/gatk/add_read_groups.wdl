version development

import "../../structs/runenv.wdl"

task run_add_read_groups {
  input {
    File bam
    String sample
    String pl = "ILLUMINA"
    String pu = "FC.L.SB"
    RunEnv runenv
  }

  String output_bam = basename(bam)
  # -I <String>      Input file (BAM or SAM or a GA4GH url).  Required. 
  # -O <File>        Output file (SAM, BAM or CRAM).  Required. 
  # -LB <String>     Read-Group library  Required. 
  # -PL <String>     Read-Group platform (e.g. ILLUMINA, SOLID)  Required. 
  # -PU <String>     Read-Group platform unit (eg. run barcode)  Required. 
  # -SM <String>     Read-Group sample name  Required. 
  int javamem = runenv.memory - 1
  command <<<
    /gatk/gatk --java-options -Xmx~{javamem}g AddOrReplaceReadGroups \
      -I ~{bam} \
      -SM ~{sample} \
      -LB ~{sample + ".lib"} \
      -PL ~{pl} \
      -PU ~{pu} \
      -O ~{output_bam}
  >>>

  output {
    File output_bam = "~{output_bam}"
  }

  runtime {
    docker: runenv.docker
    cpu : runenv.cpu
    memory : "~{runenv.memory} GB"
    #disks:
  }
}
