version development

import "../../structs/runenv.wdl"

task run_stats {
  input {
     File gam
     RunEnv runenv
  }

  String output_stats = "~{basename(gam)}.stats"
  # -z, --size             size of graph
  # -N, --node-count       number of nodes in graph
  # -E, --edge-count       number of edges in graph
  # -l, --length           length of sequences in graph
  # -L, --self-loops       number of self-loops
  # -s, --subgraphs        describe subgraphs of graph
  # -H, --heads            list the head nodes of the graph
  # -T, --tails            list the tail nodes of the graph
  # -e, --nondeterm        list the nondeterministic edge sets
  # -c, --components       print the strongly connected components of the graph
  # -A, --is-acyclic       print if the graph is acyclic or not
  # -n, --node ID          consider node with the given id
  # -d, --to-head          show distance to head for each provided node
  # -t, --to-tail          show distance to head for each provided node
  # -a, --alignments FILE  compute stats for reads aligned to the graph
  # -r, --node-id-range    X:Y where X and Y are the smallest and largest node id in the graph, respectively
  # -o, --overlap PATH    for each overlapping path mapping in the graph write a table:
  #                           PATH, other_path, rank1, rank2
  #                       multiple allowed; limit comparison to those provided
  # -O, --overlap-all     print overlap table for the cartesian product of paths
  # -R, --snarls          print statistics for each snarl
  # -F, --format          graph format from {VG-Protobuf, PackedGraph, HashGraph, XG}. Can't detect Protobuf if graph read from stdin
  # -D, --degree-dist     print degree distribution of the graph.
  # -b, --dist-snarls FILE print the sizes and depths of the snarls in a given distance index.
  # -p, --threads N       number of threads to use [all available]
  # -v, --verbose         output longer reports
  command <<<
    vg stats -p ~{runenv.cpu - 1} -a ~{gam} > ~{output_stats}
  >>>

  output {
    File stats = output_stats
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: "${runenv.memory} GB"
  }
}
