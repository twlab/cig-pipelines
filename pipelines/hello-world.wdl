version development

workflow hello_world {
  meta {
    author: "Eddie Belter"
    version: "1.0"
    description: "Hello World, an examp[le workflow"
  }

  input {
    String name
  }

  call run_hello_world as hw { input:
    name=name,
  }

  output {
    File out_fn = hw.out_fn
  }
}

task run_hello_world {
  input {
    String name
  }

  command <<<
    echo Hello world, ~{name}! > hw.out
  >>>

  runtime {
    docker: "ebelter/linux-tk:latest"
    cpu: 1
    memory: "2 GB"
    disks : "local-disk 10 SSD"
  }

  output {
    File out_fn = glob("hw.out")[0]
  }
}
