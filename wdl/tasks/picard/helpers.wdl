version development

import "../../structs/runenv.wdl"

task determine_sort_order {
  input {
    File bam
    RunEnv runenv
  }

  command <<<
    java -Xmx~{runenv.memory - 1}g -jar /usr/picard/picard.jar ViewSam \
       --INPUT ~{bam} \
       --HEADER_ONLY true | head -1 |tr "\t" "\n" | grep SO | sed 's/SO://' > out
  >>>

  # unknown (default) / unsorted / queryname / coordinate
  #String so = read_string("out")
  output {
    String sort_order = read_string("out")
    Boolean is_name_sorted = if read_string("out") == "queryname" then true else false
    Boolean is_coordinate_sorted = if read_string("out") == "coordiante" then true else false
    #Boolean is_name_sorted = if "~{so}" == "queryname" then true else false
    #Boolean is_coordinate_sorted = if "~{so}" == "coordiante" then true else false
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: runenv.memory + " GB"
    #disks: runenv.disks
  }
}

task determine_read_groups {
  input {
    File bam
    RunEnv runenv
  }

  command <<<
    java -Xmx~{runenv.memory - 1}g -jar /usr/picard/picard.jar ViewSam \
       --INPUT ~{bam} \
       --HEADER_ONLY true | grep '@RG' > out
  >>>

  Array[String] rgs = read_lines("out")
  output {
    String sort_order = rgs
    Boolean has_read_groups = if length(rgs) > 0 then true else false
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: runenv.memory + " GB"
    #disks: runenv.disks
  }
}
